
inputdir = config["input"]
outputdir = config["output"]
program_folder = config["program_folder"]

if inputdir[-1] == "/":
    inputdir = inputdir[:-1]
if outputdir[-1] == "/":
    outputdir = outputdir[:-1]


input_files =  glob_wildcards(f"{inputdir}/{{prefix, [^/]+}}.fam")

def expand_chromosomes(number_of_chrs):
    chromosomes_ = [str(i) for i in range(1, int(number_of_chrs)+1)]
    chromosomes = [str(item).zfill(2) for item in chromosomes_]
    return(chromosomes)


def get_chr_nr(sp):
    if sp == "cow":
        return(29)
    if sp == "pig":
        return(18)
    if sp == "sheep":
        return(26)
    if sp == "goat":
        return(29)
    if sp == "dog":
        return(38)
    if sp == "cat":
        return(18)
    if sp == "alpaca":
        return(36)


rule all:
    input:
        expand("{out}/locations_{prefix}.csv",prefix=input_files.prefix,  out=config["output"]), 
        expand("{out}/{prefix}_homozygosity.csv", prefix=input_files.prefix, out=config["output"]),
        expand("{out}/{prefix}_plot.pdf", prefix=input_files.prefix, out=config["output"]),
        expand("{out}/plinkfiles/{prefix}_founder.hom", prefix=input_files.prefix, out=config["output"]),
        expand("{out}/{prefix}_all_common.csv", prefix=input_files.prefix, out=config["output"]),
        expand("{out}/{prefix}_common_without_paternal_homozygosity.csv", prefix=input_files.prefix, out=config["output"])


rule filter:
    message:
        """Filtering the input files with maf and mind, checking for mendelian errors with me"""
    input:
        snp = f"{inputdir}/{{prefix}}.bed", 
        lgen = f"{inputdir}/{{prefix}}.bim", 
        sample = f"{inputdir}/{{prefix}}.fam"
    output:
        bam = outputdir + "/plinkfiles/{prefix}_filtered.bed",
        bim = outputdir + "/plinkfiles/{prefix}_filtered.bim",
        fam = outputdir + "/plinkfiles/{prefix}_filtered.fam"
    params:
        name = lambda wc :  f"{inputdir}/{wc.prefix}",
        output =  outputdir + "/plinkfiles/{prefix}_filtered",
        chrs = get_chr_nr(config["species"]),
        maf = config["filter"]["maf"],
        mind = config["filter"]["mind"],
        geno = config["filter"]["geno"],
        mendel = config["filter"]["me"]
    resources:
        mem="900M",
        time="00:02:10",
        cpus=1
    shell:
        """
        module load PLINK 
        plink --maf {params.maf} --mind {params.mind} --me {params.mendel} --geno {params.geno} --bfile {params.name} --out {params.output} --make-bed --chr-set {params.chrs}
        """
        
    
rule recode:
    message:
        """Recode bim, bam and fam to ped"""
    input:
        input_ = rules.filter.output.bim
    output:
        output = ("{outputdir}/plinkfiles/{prefix}_recoded.ped"),
        map_ = (  "{outputdir}/plinkfiles/{prefix}_recoded.map")
    params:
        input_ =  "{outputdir}/plinkfiles/{prefix}_filtered",
        output =  "{outputdir}/plinkfiles/{prefix}_recoded",
        chrs = get_chr_nr(config["species"])
    resources:
        mem="900M",
        time="00:05:05",
        cpus=1
    shell:
        """
        module load PLINK 
        plink --recode 12 --bfile {params.input_} --out {params.output} --chr-set {params.chrs} --tab 
        """


rule rename_genes:
    message:
        """ Renaming genes to standardised format.  map  -> csv """
    input:
        input_ = "{outputdir}/plinkfiles/{prefix}_filtered.bim"
    output:
        filtered_csv = "{outputdir}/plinkfiles/{prefix}_filtered.csv"
    resources:
        mem="4G",
        time="00:06:06",
        cpus=1
    shell:
        """
        python3 {program_folder}/code/rename_genes.py \
        --input {input.input_} \
        --output {output.filtered_csv}
        """


rule paternal_haplotypes:
    message:
        """Finding paternal haplotypes in halfsibs"""
    input:
        ped = rules.recode.output.output,
        gene_map = rules.rename_genes.output.filtered_csv
    output:
        paternal_haplotypes = "{outputdir}/haplotypes/{prefix}_{chr}_output.csv"
    params:
        folder = "{outputdir}/",
        chromosomes = lambda wc: wc.get('chr')
    resources:
        mem="20G",
        time="00:10:00",
        cpus=3
    shell:
        """
        python3 {program_folder}/code/halfsib_paternal_haplotypes.py \
        --ped {input.ped} \
        --output {output.paternal_haplotypes} \
        --markers {input.gene_map} \
        --chr {params.chromosomes}
        """


rule linked_haplotypes:
    message:
        """Finding linked locations..."""
    input:
        haplotype = rules.paternal_haplotypes.output.paternal_haplotypes,
        map_ = rules.rename_genes.output
    output:
        linked = ("{outputdir}/{prefix}_locations_{chr}.csv")
    params:
        min_length = config['variants']['min_var_length'],
        n_fraction_max = config["variants"]["n_fraction_max"],
        fuse_adjacent = config["variants"]["fuse_adjacent"],
        fuse_adjacent_nr = config["variants"]["fuse_adjacent_nr"],
        min_markers = config["variants"]["min_markers"],
        chrs = lambda wc: wc.get('chr')
    resources:
        mem="10G",
        time="00:05:06",
        cpus=2
    shell:
        """
        python3 {program_folder}/code/common_locations.py \
        --map {input.map_} \
        --hapl {input.haplotype} \
        --min_markers {params.min_markers} \
        --length {params.min_length} \
        --n_fraction_max {params.n_fraction_max} \
        --fuse_adjacent {params.fuse_adjacent} \
        --fuse_adjacent_nr {params.fuse_adjacent_nr} \
        --output {output.linked} \
        --chr {params.chrs}
        """

rule merge_linked:
    message:
        """Merging linked files for all chrs"""
    input:
        linked = expand("{{outputdir}}/{{prefix}}_locations_{chr}.csv", chr=expand_chromosomes(get_chr_nr(config["species"])))
    output:
        linked = "{outputdir}/locations_{prefix}.csv"
    resources:
        mem="900M",
        time="00:05:04",
        cpus=1
    shell:
        """
        echo "CHR;BP1;BP2\n" >> {output}
        cat {input} >> {output}
        """

rule filter_founder:
    message:
        """Writing file to filter for founder"""
    input:
        rules.filter.output.fam
    output:
        "{outputdir}/{prefix}_IDlist.txt"
    resources:
        mem="900M",
        time="00:05:04",
        cpus=1
    shell:
        """
        python3 {program_folder}/code/founder.py \
        --input {input} \
        --output {output}
        """
    

rule founder:
    message:
        """Filtering non-founder data based on provided ID list"""
    input:
        input_ =  "{outputdir}/plinkfiles/{prefix}_recoded.map",
        ids_to_remove = rules.filter_founder.output
    output:
        output = "{outputdir}/plinkfiles/{prefix}_founder.bed"
    params:
        in_ =  "{outputdir}/plinkfiles/{prefix}_recoded",
        out = "{outputdir}/plinkfiles/{prefix}_founder",
        chrs = get_chr_nr(config["species"])
    resources:
        mem="900M",
        time="00:05:04",
        cpus=1
    shell:
        """
        module load PLINK 
        plink --file {params.in_} --remove {input.ids_to_remove} --out {params.out} --make-bed --chr-set {params.chrs} 
        """
        
rule homozygosity:
    message:
        """Find homozygous locations based on founder data"""
    input:
        input = rules.founder.output.output
    output:
        output = "{outputdir}/plinkfiles/{prefix}_founder.hom"
    params:
        input_ = "{outputdir}/plinkfiles/{prefix}_founder",
        output = "{outputdir}/plinkfiles/{prefix}_founder",
        chrs = get_chr_nr(config["species"]),
        density = config["homozygosity_params"]["density"],
        kb = config["homozygosity_params"]["kb"],
        snp = config["homozygosity_params"]["snp"],
        window_missing = config["homozygosity_params"]["window_missing"],
        window_snp = config["homozygosity_params"]["window_snp"],
        window_het = config["homozygosity_params"]["window_het"]
    resources:
        mem="900M",
        time="00:05:03",
        cpus=1
    shell:
        """
        module load PLINK 
        plink --bfile {params.input_} --homozyg --chr-set {params.chrs} --homozyg-density {params.density} --homozyg-kb {params.kb} --homozyg-snp {params.snp} --homozyg-window-missing {params.window_missing} --homozyg-window-snp {params.window_snp} --homozyg-window-het {params.window_het} --out {params.output}
        """

rule reformat_homozygosity:
    message:
        """Reformatting homozygous areas to a csv file for annotation"""
    input:
        input_ = rules.homozygosity.output
    output:
        output =  "{outputdir}/{prefix}_homozygosity.csv"
    resources:
        mem="500M",
        time="00:05:05",
        cpus=1
    shell:
        """
        python3 {program_folder}/code/homozygosity.py \
        --input {input.input_} \
        --output {output.output}
        """

rule output_table:
    message:
        """Finding variants which are not homozygous in the paternal genome"""
    input:
        homozygous = rules.reformat_homozygosity.output,
        linked = rules.merge_linked.output
    output:
        table =  "{outputdir}/{prefix}_common_without_paternal_homozygosity.csv",
        all_common = "{outputdir}/{prefix}_all_common.csv"
    resources:
        mem="500M",
        time="00:05:04",
        cpus=1
    shell:
        """
        python3 {program_folder}/code/write_tables.py \
        --homozygosity {input.homozygous} \
        --linked {input.linked} \
        --with_homozyg {output.all_common} \
        --common {output.table}
        """

rule plot:
    message:
        """constructing chromosome map plot with homozygosity and linked haplotypes"""
    input:
        linked = rules.merge_linked.output,
        chr_map = f"{program_folder}/chr_maps/" + config['species'] + "_chr_map.csv",
        homozyg = rules.reformat_homozygosity.output,
        linked_singular = expand("{outputdir}/{prefix}_locations_{chr}.csv", outputdir = config["output"],prefix=input_files.prefix,  chr=expand_chromosomes(get_chr_nr(config["species"]))),
    output:
        plot = "{outputdir}/{prefix}_plot.pdf"
    resources:
        mem="500M",
        time="00:05:05",
        cpus=1
    params:
        chrs = get_chr_nr(config["species"])
    shell:
        """
        module load matplotlib
        python3 {program_folder}/code/plot.py \
        --chr_file {input.chr_map} \
        --chr_nr {params.chrs} \
        --linked {input.linked} \
        --homozyg {input.homozyg} \
        --plot {output.plot}
        rm -f {input.linked_singular}
        """

