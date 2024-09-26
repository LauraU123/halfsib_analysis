
configfile: "config/configfile.yaml"
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
        expand(outputdir +  "{eg}/table.csv", eg=config["example"])

rule filter:
    message:
        """Filtering the input files with maf and mind, checking for mendelian errors with me"""
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
        """
        module load PLINK 
        plink --maf 0.01 --mind 0.1 --me 0.05 0.1 --bfile {params.name} --out {params.output} --make-bed --chr-set {params.chrs}
        """
        
    
rule recode:
    message:
        """Recode bim, bam and fam to ped"""
    input:
        input_ = rules.filter.output.bim
    output:
        output = outputdir + "{example}/recoded.ped",
        map_ = outputdir +  "{example}/recoded.map"
    params:
        input_ = outputdir +  "{example}/filtered",
        output = outputdir + "{example}/recoded",
        chrs = config["chrs"]
    shell:
        """
        module load PLINK 
        plink --recode 12 --bfile {params.input_} --out {params.output} --chr-set {params.chrs} --tab 
        """


rule rename_genes:
    message:
        """ Renaming genes to standardised format.  map  -> csv """
    input:
        input_ = outputdir +  "{example}/filtered.bim"
    output:
        filtered_csv = outputdir + "{example}/filtered.csv"
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
        paternal_haplotypes = outputdir + "{example}/{chr}_output.csv"
    params:
        folder = outputdir + "{example}/",
        chromosomes = lambda wc: wc.get('chr')
    shell:
        """
        python3 code/halfsib_v2.py \
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
        linked = outputdir +  "{example}/locations_{chr}.csv"
    params:
        min_length = config['variants']['min_var_length'],
        n_fraction_max = config["variants"]["n_fraction_max"],
        fuse_adjacent = config["variants"]["fuse_adjacent"],
        fuse_adjacent_nr = config["variants"]["fuse_adjacent_nr"],
        min_markers = config["variants"]["min_markers"],
        chrs = lambda wc: wc.get('chr')
    shell:
        """
        python3 code/locations_v2.py \
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
        linked = expand(outputdir + "{{example}}/locations_{chr}.csv", chr=expand_chromosomes(config["chrs"]))
    output:
        linked = outputdir +  "{example}/locations.csv"
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
        outputdir + "{example}/IDlist.txt"
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
        input_ = outputdir +  "{example}/recoded.map",
        ids_to_remove = rules.filter_founder.output
    output:
        output = outputdir + "{example}/founder.bed"
    params:
        in_ = outputdir +  "{example}/recoded",
        out = outputdir +  "{example}/founder",
        chrs = config["chrs"]
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
        module load PLINK 
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

rule output_table:
    message:
        """Finding variants which are not homozygous in the paternal genome"""
    input:
        homozygous = rules.reformat_homozygosity.output,
        variants = rules.merge_linked.output
    output:
        table = outputdir +  "{example}/table.csv"
    shell:
        """
        python3 code/write_tables.py \
        --homozygosity {input.homozygous} \
        --variants {input.variants} \
        --output_csv {output.table}
        """

rule plot:
    message:
        """constructing chromosome map plot with homozygosity and linked haplotypes"""
    input:
        variants = rules.merge_linked.output,
        chr_map_cattle = inputdir + "chr_map.csv",
        homozyg = rules.reformat_homozygosity.output,
    output:
        plot = outputdir + "{example}/plot.pdf"
    shell:
        """
        module load matplotlib
        python3 code/plot.py \
        --chr {input.chr_map_cattle} \
        --variants {input.variants} \
        --homozyg {input.homozyg} \
        --plot {output.plot}
        """