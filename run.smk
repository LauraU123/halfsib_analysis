#this is the snakefile for running all the input examples
configfile: "configfile.yaml"
outputdir = "results/"
inputdir = "data/"

chromosomes = ["01","02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29"]

rule all:
    input:
        locations = outputdir + "example1/locations.csv", 
        homozygosity = outputdir +  "example1/homozygosity.csv",
        plot = outputdir + "example1/plot.pdf",
        founder = outputdir +  "example1/founder.hom"

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
        output = outputdir + "{example}/filtered"
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
        mendel = outputdir + "{example}/filtered_mendelian.bim"
    params:
        input_ = outputdir + "{example}/filtered",
        output = outputdir + "{example}/filtered_mendelian"
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
        output = outputdir + "{example}/recoded.ped",
        map_ = outputdir +  "{example}/recoded.map"
    params:
        input_ = outputdir +  "{example}/filtered_mendelian",
        output = outputdir + "{example}/recoded"
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
        common_sequences = expand(outputdir + "{{example}}/{chr}_output.csv", chr=chromosomes)
    shell:
        """
        python3 code/halfsib_v2.py \
        --ped {input.ped} \
        --markers {input.gene_map} 
        """

rule locations:
    message:
        """Find the common locations from the files"""
    input:
        input_ = rules.halfsib.output.common_sequences,
        map_ = rules.rename_genes.output
    output:
        locations = outputdir +  "{example}/locations.csv"
    params:
        min_length = 1200000
    shell:
        """
        python3 code/locations_v2.py \
        --map {input.map_} \
        --locations {output.locations} \
        --length {params.min_length} 
        """

rule founder:
    input:
        input_ = outputdir +  "{example}/recoded.map",
        ids_to_remove = inputdir + "{example}/IDlist.txt"
    output:
        output = outputdir + "{example}/founder.bed"
    params:
        in_ = outputdir +  "{example}/recoded",
        out = outputdir +  "{example}/founder"
    shell:
        """
        plink --file {params.in_} --remove {input.ids_to_remove} --out {params.out} --make-bed --chr-set 29 
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
        output = outputdir +  "{example}/founder"
    shell:
        """
        plink --bfile {params.input_} --homozyg --chr-set 29 --homozyg-density 50 --homozyg-kb 1000 --homozyg-snp 20 --homozyg-window-missing 5 --homozyg-window-snp 50 --out {params.output}
        """

rule reformat_homozygosity:
    message:
        """Reformat homozygous areas to a csv file for annotation"""
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
        variants = rules.locations.output.locations
    output:
        only_variants = outputdir +  "{example}/only_homozygous.csv"
    shell:
        """
        python3 code/only_variants.py \
        --homozygosity {input.homozygous} \
        --variants {input.variants} \
        --output_csv {output.only_variants}
        """

rule plot:
    message:
        """constructing chromosome map plot with homozygosity and common sequences"""
    input:
        only_variants = rules.only_variants.output,
        chr_map_cattle = inputdir + "chr_map.csv",
    output:
        plot = outputdir + "{example}/plot.pdf"
    shell:
        """
        python3 code/plot.py \
        --chr {input.chr_map_cattle} \
        --locations {input.only_variants} \
        --plot {output.plot}
        """