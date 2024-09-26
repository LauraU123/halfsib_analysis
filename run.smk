
configfile: "configfile.yaml"
outputdir = "results/"
inputdir = "data/"

def expand_chromosomes(number_of_chrs):
    chromosomes_ = [str(i) for i in range(1, int(number_of_chrs)+1)]
    chromosomes = [str(item).zfill(2) for item in chromosomes_]
    return(chromosomes)


rule all:
    input:
        expand(outputdir + "{eg}/locations.csv",  eg=config["example"]), 
        expand(outputdir + "{eg}/homozygosity.csv", eg=config["example"]),
        expand(outputdir  + "{eg}/plot.pdf", eg=config["example"]),
        expand(outputdir  + "{eg}/founder.hom", eg=config["example"]),
        expand(outputdir +  "{eg}/only_homozygous.csv", eg=config["example"])

rule filter:
    message:
        """Filtering the input files with maf and mind"""
    input:
        snp = inputdir + "{example}/plink.bed",
        lgen = inputdir + "{example}/plink.bim",
        sample = inputdir + "{example}/plink.fam"
    output:
        bam = outputdir + "{example}/filtered.bed",
        bim = outputdir + "{example}/filtered.bim",
        fam = outputdir + "{example}/filtered.fam"
    params:
        name = inputdir + "{example}/plink",
        output = outputdir + "{example}/filtered",
        chrs = config["chrs"]
    shell:
        "module load snakemake"
        "module load matplotlib"
        "module load PLINK"
        "plink --maf 0.01 --mind 0.1 --bfile {params.name} --out {params.output} --make-bed --chr-set {params.chrs}"
        


rule mendel:
    message:
        """Filtering mendelian errors with plink"""
    input:
        bam = rules.filter.output.bam,
        bim = rules.filter.output.bim,
        fam = rules.filter.output.fam
    output:
        mendel = outputdir + "{example}/filtered_mendelian.bim"
    params:
        input_ = outputdir + "{example}/filtered",
        output = outputdir + "{example}/filtered_mendelian",
        chrs = config["chrs"]
    shell:
        """
        plink --me 0.05 0.1 --bfile {params.input_}  --chr-set {params.chrs} --out {params.output} --make-bed
        """
    
rule recode:
    message:
        """Recode bim, bam and fam to ped"""
    input:
        input_ = rules.mendel.output.mendel
    output:
        output = outputdir + "{example}/recoded.ped",
        map_ = outputdir +  "{example}/recoded.map"
    params:
        input_ = outputdir +  "{example}/filtered_mendelian",
        output = outputdir + "{example}/recoded",
        chrs = config["chrs"]
    shell:
        """
        plink --recode 12 --bfile {params.input_} --out {params.output} --chr-set {params.chrs} --tab 
        """


rule rename_genes:
    message:
        """
        Rename genes to standardised format. filtered map file to csv file
        """
    input:
        input_ = outputdir +  "{example}/filtered_mendelian.bim"
    output:
        filtered_csv = outputdir + "{example}/filtered.csv"
    shell:
        """
        python3 code/rename_genes.py \
        --input {input.input_} \
        --output {output.filtered_csv}
        """


rule halfsib:
    message:
        """Find common halfsibs from input files"""
    input:
        ped = rules.recode.output.output,
        gene_map = rules.rename_genes.output.filtered_csv
    output:
        common_sequences = outputdir + "{example}/{chr}_output.csv"
    params:
        folder = outputdir + "{example}/",
        chromosomes = config["chrs"]
    shell:
        """
        python3 code/halfsib_v2.py \
        --ped {input.ped} \
        --folder {params.folder} \
        --markers {input.gene_map} \
        --chr {params.chromosomes}
        """


rule linked:
    message:
        """Finding linked locations..."""
    input:
        input_ = rules.halfsib.output.common_sequences,
        map_ = rules.rename_genes.output
    output:
        linked = outputdir +  "locations_{chr}.csv"
    params:
        min_length = config['variants']['min_var_length'],
        folder = outputdir + "{example}/",
        n_fraction_max = config["variants"]["n_fraction_max"],
        fuse_adjacent = config["variants"]["fuse_adjacent"],
        fuse_adjacent_nr = config["variants"]["fuse_adjacent_nr"],
        min_markers = config["variants"]["min_markers"]
    shell:
        """
        python3 code/locations_v2.py \
        --map {input.map_} \
        --locations {output.variants} \
        --folder {params.folder} \
        --min_markers {params.min_markers} \
        --length {params.min_length} \
        --n_fraction_max {params.n_fraction_max} \
        --fuse_adjacent {params.fuse_adjacent} \
        --fuse_adjacent_nr {params.fuse_adjacent_nr} \
        --output {output}
        """

rule merge_linked:
    message:
        """Merging linked files for all chrs"""
    input:
        linked = expand("locations_{chr}.csv", chr=expand_chromosomes(config["chrs"]))
    output:
        linked = outputdir +  "{example}/locations.csv"
    shell:
        """
        cat {input} > {output}
        """
    

rule founder:
    input:
        input_ = outputdir +  "{example}/recoded.map",
        ids_to_remove = inputdir + "{example}/IDlist.txt"
    output:
        output = outputdir + "{example}/founder.bed"
    params:
        in_ = outputdir +  "{example}/recoded",
        out = outputdir +  "{example}/founder",
        chrs = config["chrs"]
    shell:
        """
        plink --file {params.in_} --remove {input.ids_to_remove} --out {params.out} --make-bed --chr-set {params.chrs} 
        """
        
rule homozygosity:
    message:
        """Find the homozygosity locations"""
    input:
        input = rules.founder.output.output
    output:
        output = outputdir +  "{example}/founder.hom"
    params:
        input_ = outputdir +  "{example}/founder",
        output = outputdir +  "{example}/founder",
        chrs = config["chrs"],
        density = config["homozygosity_params"]["density"],
        kb = config["homozygosity_params"]["kb"],
        snp = config["homozygosity_params"]["snp"],
        window_missing = config["homozygosity_params"]["window_missing"],
        window_snp = config["homozygosity_params"]["window_snp"]
    shell:
        """
        plink --bfile {params.input_} --homozyg --chr-set {params.chrs} --homozyg-density {params.density} --homozyg-kb {params.kb} --homozyg-snp {params.snp} --homozyg-window-missing {params.window_missing} --homozyg-window-snp {params.window_snp} --out {params.output}
        """

rule reformat_homozygosity:
    message:
        """Reformatting homozygous areas to a csv file for annotation"""
    input:
        input_ = rules.homozygosity.output
    output:
        output = outputdir +  "{example}/homozygosity.csv"
    shell:
        """
        python3 code/homozygosity.py \
        --input {input.input_} \
        --output {output.output}
        """

rule only_variants:
    message:
        """Finding variants which are not homozygous in the paternal genome"""
    input:
        homozygous = rules.reformat_homozygosity.output,
        variants = rules.merge_linked.output
    output:
        only_variants = outputdir +  "{example}/only_homozygous.csv"
    shell:
        """
        python3 code/write_tables.py \
        --homozygosity {input.homozygous} \
        --variants {input.variants} \
        --output_csv {output.only_variants}
        """

rule plot:
    message:
        """constructing chromosome map plot with homozygosity and common sequences"""
    input:
        variants = rules.merge_linked.output,
        chr_map_cattle = inputdir + "chr_map.csv",
        homozyg = rules.reformat_homozygosity.output,
    output:
        plot = outputdir + "{example}/plot.pdf"
    shell:
        """
        python3 code/plot.py \
        --chr {input.chr_map_cattle} \
        --variants {input.variants} \
        --homozyg {input.homozyg} \
        --plot {output.plot}
        """