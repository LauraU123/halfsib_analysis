#this is the snakefile for running all the input examples
#config: "configfile.yaml"
#outputdir: "results/"

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
        bam = rules.filtered.output.bam,
        bim = rules.filtered.output.bim,
        fam = rules.filtered.output.fam
    output:
        mendel = outputdir + "filtered.mendel"
    params:
        input_ = "results/{example}/filtered",
        output = "results/{example}/filtered_mendelian"
    shell:
        """
        plink --mendel --bfile {params.input_}  --chr-set 29 --out {params.output}
        """
        
rule filter_mendel:
    message:
        """Filtering mendelian errors with plink"""
    input:
        bam = rules.filtered.output.bam,
        bim = rules.filtered.output.bim,
        fam = rules.filtered.output.fam
        mendel = rules.mendel.output
    output:
        a = outputdir + "filtered_mendel.bam",
        b = outputdir + "filtered_mendel.bim",
        c = outputdir + "filtered_mendel.fam"
    params:
        input_ = "results/{example}/filtered",
        output = "results/{example}/filtered_mendelian"
    shell:
        """
        plink --mendel --bfile {params.input_} --make-bed --chr-set 29 --out {params.output}
        """

rule recode:
    message:
        """Recode bim, bam and fam to ped"""
    input:
        input_ = rules.filter.output
    output:
        output = "results/{example}/recoded.ped"
    params:
        input_ = "results/{example}/filtered",
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
        input_ = "results/{example}/filtered.bed"
    output:
        output = "results/{example}/filtered.csv"
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
        gene_map = rules.rename_genes.output.output
    output:
        common_sequences = "results/{example}/{config.chromosomes}_output.csv"
    shell:
        """
        python3 code/halfsib_v2.py \
        --ped {input.ped} \
        --markers {input.gene_map} \
        --output {output.common_sequences} 
        """

rule locations:
    message:
        """Find the common locations from the files"""
    input:
        input_ = rules.halfsib.output.common_sequences
    output:
        output = "results/{example}/locations.csv"
    shell:
        """
        python3 code/locations_v2.py \
        --input {input.input_} \
        --output {output.output}
        """

rule founder:
    input:
        input_ = "results/{example}/filtered.bim",
        ids_to_remove = "data/{example}/IDlist.txt"
    output:
        output =    "results/{example}/founder.bim"
    params:
        in_ = "results/{example}/filtered",
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
        input_ = "results/{example}/filtered",
        output = "results/{example}/founder"
    shell:
        """
        plink --file {params.input_} --homozyg --keep data/example1/founder.txt --chr-set 29 --homozyg-density 1000 --homozyg-kb 100 --homozyg-snp 50 --homozyg-window-missing 3 --homozyg-window-snp 50
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