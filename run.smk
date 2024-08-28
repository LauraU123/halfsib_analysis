#this is the snakefile for running all the input examples
config: "configfile.yaml"
outputdir: "results/"
rule all:
    input:
        locations = "results/example1/locations.csv", 
        homozygosity = "results/example1/homozygosity.csv"

rule inputs:
    message:
        """This rule generates plink input files from raw data"""
    input:
        snp = "results/{example}/SNP_Map.txt",
        sample = "results/{example}/Sample_Map.txt",
        lgen = "results/{example}/FinalReport.txt"
    output:
        snp = "results/" + "{example}/plink.map",
        lgen = "results/" + "{example}/plink.bim",
        sample = "results/" + "{example}/plink.fam",
    params:
        output = "results/{example}/plink"
    shell:
        """
        python3 code/SNP_generator.py \
            --SNP_in {input.snp} \
            --SNP_out {input.snp} \
            --experiment_in {input.lgen} \
            --experiment {output.lgen} \
            --sample_in {input.sample} \
            --sample_out {output.sample}
        """

rule filter:
    message:
        """Filtering the input files with maf and mind"""
    input:
        rules.inputs.output.snp,
        rules.inputs.output.lgen,
        rules.inputs.output.sample
    output:
        bam = "results/{example}/filtered.bam",
        bim = "results/{example}/filtered.bim",
        fam = "results/{example}/filtered.fam",
    params:
        name = "results/{example}/plink",
        output = "results/{example}/filtered"
    shell:
        """plink --maf 0.01 --mind 0.1 --lfile {params.name} --make-bed --chr-set 29 --out {params.output}"""

rule recode:
    message:
        """Recode bim, bam and fam to ped"""
    input:
        input_ = rules.filter.output
    output:
        output = "results/{example}/recoded.ped",
    shell:
        """
        plink --recode 12 --bfile {input.input_} --tab --out {output.output} --chr-set 29
        """

rule rename_genes:
    message:
        """Rename genes to standardised format. filtered map file to csv file"""
    input:
        input_ = "results/{example}/filtered.map"
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


rule homozygosity:
    message:
        """Find the homozygosity locations"""
    input:
        input_ = "results/{example}/filtered"
    output:
        output = "results/{example}/plink.hom"
    params:
        input_ = "results/{example}/filtered",
        output = "results/{example}/plink"
    shell:
        """plink --file {params.input_} --homozyg {config.founder} --chr-set 29 --homozyg-density 1000 --homozyg-kb 100 --homozyg-snp 50 --homozyg-window-missing 3 --homozyg-window-snp 50"""

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
