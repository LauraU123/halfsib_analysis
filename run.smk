#this is the snakefile for running all the input examples
#config: "configfile.yaml"
#outputdir: "results/"

chromosomes = ["01","02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29"]


rule all:
    input:
        locations = "results/example1/locations.csv", 
        homozygosity = "results/example1/homozygosity.csv"

rule filter:
    message:
        """Filtering the input files with maf and mind"""
    input:
        snp = "data/" + "{example}/plink.bed",
        lgen = "data/" + "{example}/plink.bim",
        sample = "data/" + "{example}/plink.fam"
    output:
        bam = "results/{example}/filtered.bed",
        bim = "results/{example}/filtered.bim",
        fam = "results/{example}/filtered.fam"
    params:
        name = "data/{example}/plink",
        output = "results/{example}/filtered"
    shell:
        """
        plink --maf 0.01 --mind 0.1 --bfile {params.name} --out {params.output} --make-bed --chr-set 29
        """

rule mendel:
    message:
        """Filtering mendelian errors with plink"""
    input:
        bam = rules.filter.output.bam,
        bim = rules.filter.output.bim,
        fam = rules.filter.output.fam
    output:
        mendel = "results/{example}/filtered_mendelian.bim"
    params:
        input_ = "results/{example}/filtered",
        output = "results/{example}/filtered_mendelian"
    shell:
        """
        plink --me 0.05 0.1 --bfile {params.input_}  --chr-set 29 --out {params.output} --make-bed
        """
    

rule recode:
    message:
        """Recode bim, bam and fam to ped"""
    input:
        input_ = rules.mendel.output.mendel
    output:
        output = "results/{example}/recoded.ped"
    params:
        input_ = "results/{example}/filtered_mendelian",
        output = "results/{example}/recoded"
    shell:
        """
        plink --recode 12 --bfile {params.input_} --out {params.output} --chr-set 29 --tab 
        """

rule rename_genes:
    message:
        """
        Rename genes to standardised format. filtered map file to csv file
        """
    input:
        input_ = "results/{example}/filtered_mendelian.bim"
    output:
        filtered_csv = "results/{example}/filtered.csv"
    shell:
        """
        python3 code/rename_genes.py \
        --input {input.input_} \
        --output {output.output}
        """

rule halfsib:
    message:
        """Find common halfsibs from input files"""
    input:
        ped = rules.recode.output,
        gene_map = rules.rename_genes.output.filtered_csv
    output:
        common_sequences = expand("results/{{example}}/{chr}_output.csv", chr=chromosomes)
    shell:
        """
        python3 code/halfsib_v2.py \
        --ped {input.ped} \
        --markers {input.gene_map} \
        """

rule locations:
    message:
        """Find the common locations from the files"""
    input:
        input_ = rules.halfsib.output.common_sequences,
        map_ = rules.rename_genes.output
    output:
        locations = "results/{example}/locations.csv",
        output = "results/{example}/chromosomes.csv",
        comparison = "results/{example}/compare.csv"
    params:
        min_length = 900000
    shell:
        """
        python3 code/locations_v2.py \
        --map {input.map_} \
        --markers {input.input_} \
        --locations {output.locations} \
        --length {params.min_length} \
        --top {output.output} \
        --comparison {output.comparison}
        """

rule founder:
    input:
        input_ = "results/{example}/recoded.map",
        ids_to_remove = "data/{example}/IDlist.txt"
    output:
        output =    "results/{example}/founder.bim"
    params:
        in_ = "results/{example}/recoded",
        out = "results/{example}/founder"
    shell:
        """
        plink --file {params.in_} --remove {input.ids_to_remove} --out {params.out} --make-bed --chr-set 29 
        """
        
rule homozygosity:
    message:
        """Find the homozygosity locations"""
    input:
        rules.founder.output
    output:
        output = "results/{example}/founder.hom"
    params:
        input_ = "results/{example}/founder",
        output = "results/{example}/founder"
    shell:
        """
        plink --bfile {params.input_} --homozyg --chr-set 29 --homozyg-density 1000 --homozyg-kb 100 --homozyg-snp 50 --homozyg-window-missing 3 --homozyg-window-snp 50 --out {params.output}
        """

rule reformat_homozygosity:
    message:
        """Reformat homozygous areas to a csv file for annotation"""
    input:
        input_ = rules.homozygosity.output
    output:
        output = "results/{example}/homozygosity.csv"
    shell:
        """
        python3 code/homozygosity.py \
        --input {input.input_} \
        --output {output.output}
        """

rule plot:
    message:
        """constructing chromosome map plot with homozygosity and common sequences"""
    input:
        homozygosity = rules.reformat_homozygosity.output,
        chr_map_cattle = "data/chr_map.csv",
        common = rules.locations.output.locations
    output:
        plot = "results/{example}/plot.png"
    shell:
        """
        python3 code/plot.py \
        --homozygosity {input.homozygosity} \
        --chr_map {input.chr_map_cattle} \
        --common {input.common} \
        --output {output.plot}        
        """
