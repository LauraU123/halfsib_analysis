
#configfile: "config/configfile.yaml"

def expand_chromosomes(number_of_chrs):
    chromosomes_ = [str(i) for i in range(1, int(number_of_chrs)+1)]
    chromosomes = [str(item).zfill(2) for item in chromosomes_]
    return(chromosomes)


def get_chr_nr(species):
    if species == "cow" or "Bos taurus":
        return(29)
    elif species == "pig" or "Sus scrofa":
        return(31)
    elif species == "sheep" or "Ovis aries":
        return(26)
    elif species == "goat" or "Capra hircus":
        return(29)
    elif species == "dog" or "Canis familiaris":
        return(38)
    elif species == "cat" or "Felis cattus":
        return(18)
    elif species == "alpaca" or "Vicugna pacos":
        return(36)


rule all:
    input:
        expand("{eg}/locations.csv",  eg=config["output"]), 
        expand("{eg}/homozygosity.csv", eg=config["output"]),
        expand("{eg}/plot.pdf", eg=config["output"]),
        expand("{eg}/founder.hom", eg=config["output"]),
        expand("{eg}/all_common.csv", eg=config["output"]),
        expand( "{eg}/common_without_paternal_homozygosity.csv", eg=config["output"])

input_files = glob_wildcards("data/{prefix}.fam")

rule filter:
    message:
        """Filtering the input files with maf and mind, checking for mendelian errors with me"""
    input:
        snp = expand("{{example}}/{prefix}.bed", prefix=input_files.prefix),
        lgen = expand("{{example}}/{prefix}.bim", prefix=input_files.prefix),
        sample = expand("{{example}}/{prefix}.fam", prefix=input_files.prefix)
    output:
        bam = "{example}/filtered.bed",
        bim = "{example}/filtered.bim",
        fam = "{example}/filtered.fam"
    params:
        name =  "{example}/{prefix}",
        output =  "{example}/filtered",
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
        plink --maf{params.maf} --mind {params.mind} --me {params.mendel} --geno {params.geno} --bfile {params.name} --out {params.output} --make-bed --chr-set {params.chrs}
        """
        
    
rule recode:
    message:
        """Recode bim, bam and fam to ped"""
    input:
        input_ = rules.filter.output.bim
    output:
        output = ("{example}/recoded.ped"),
        map_ = (  "{example}/recoded.map")
    params:
        input_ =  "{example}/filtered",
        output =  "{example}/recoded",
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
        input_ = "{example}/filtered.bim"
    output:
        filtered_csv = "{example}/filtered.csv"
    resources:
        mem="4G",
        time="00:06:06",
        cpus=1
    shell:
        """
        python3 code/rename_genes.py \
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
        paternal_haplotypes = "{example}/{chr}_output.csv"
    params:
        folder = "{example}/",
        chromosomes = lambda wc: wc.get('chr')
    resources:
        mem="20G",
        time="00:10:00",
        cpus=3
    shell:
        """
        python3 code/halfsib_paternal_haplotypes.py \
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
        linked = ("{example}/locations_{chr}.csv")
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
        python3 code/common_locations.py \
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
        linked = expand("{{example}}/locations_{chr}.csv", chr=expand_chromosomes(get_chr_nr(config["species"])))
    output:
        linked = "{example}/locations.csv"
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
        "{example}/IDlist.txt"
    resources:
        mem="900M",
        time="00:05:04",
        cpus=1
    shell:
        """
        python3 code/founder.py \
        --input {input} \
        --output {output}
        """
    

rule founder:
    message:
        """Filtering non-founder data based on provided ID list"""
    input:
        input_ =  "{example}/recoded.map",
        ids_to_remove = rules.filter_founder.output
    output:
        output = "{example}/founder.bed"
    params:
        in_ =  "{example}/recoded",
        out = "{example}/founder",
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
        output = "{example}/founder.hom"
    params:
        input_ = "{example}/founder",
        output = "{example}/founder",
        chrs = get_chr_nr(config["species"]),
        density = config["homozygosity_params"]["density"],
        kb = config["homozygosity_params"]["kb"],
        snp = config["homozygosity_params"]["snp"],
        window_missing = config["homozygosity_params"]["window_missing"],
        window_snp = config["homozygosity_params"]["window_snp"]
    resources:
        mem="900M",
        time="00:05:03",
        cpus=1
    shell:
        """
        module load PLINK 
        plink --bfile {params.input_} --homozyg --chr-set {params.chrs} --homozyg-density {params.density} --homozyg-kb {params.kb} --homozyg-snp {params.snp} --homozyg-window-missing {params.window_missing} --homozyg-window-snp {params.window_snp} --out {params.output}
        """

rule reformat_homozygosity:
    message:
        """Reformatting homozygous areas to a csv file for annotation"""
    input:
        input_ = rules.homozygosity.output
    output:
        output =  "{example}/homozygosity.csv"
    resources:
        mem="500M",
        time="00:05:05",
        cpus=1
    shell:
        """
        python3 code/homozygosity.py \
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
        table =  "{example}/common_without_paternal_homozygosity.csv",
        all_common = "{example}/all_common.csv"
    resources:
        mem="500M",
        time="00:05:04",
        cpus=1
    shell:
        """
        python3 code/write_tables.py \
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
        chr_map_cattle = lambda wc: f"{wc.get('species')}_chr_map.csv",
        homozyg = rules.reformat_homozygosity.output,
    output:
        plot = "{example}/plot.pdf"
    resources:
        mem="500M",
        time="00:05:05",
        cpus=1
    params:
        chrs = get_chr_nr(config["species"])
    shell:
        """
        module load matplotlib
        python3 code/plot.py \
        --chr_file {input.chr_map_cattle} \
        --chr_nr {params.chrs} \
        --linked {input.linked} \
        --homozyg {input.homozyg} \
        --plot {output.plot}
        """

rule cleanup:
    input:
        a= expand("{eg}/locations.csv",  eg=config["output"]), 
        b= expand("{eg}/homozygosity.csv", eg=config["output"]),
        c = expand("{eg}/plot.pdf", eg=config["output"]),
        e = expand("{eg}/all_common.csv", eg=config["output"]),
        f = expand( "{eg}/common_without_paternal_homozygosity.csv", eg=config["output"]),
        linked = expand("{eg}/locations_{chr}.csv", eg = config["output"], chr=expand_chromosomes(get_chr_nr(config["species"]))),
        haplotypes = expand("{eg}/{chr}_output.csv", eg = config["output"], chr=expand_chromosomes(get_chr_nr(config["species"])))
    resources:
        mem="500M",
        time="00:05:05",
        cpus=1
    shell:
        "rm -f {input.linked} {input.haplotypes}"     