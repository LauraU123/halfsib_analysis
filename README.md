# Paternal Halfsibling Linkage Analysis tool

This tool can be used for finding common haplotypes in paternal halfsiblings.
Chromosome number of the species of interest should be specified in the config file. 
Other parameters such as variant length and homozygosity parameters can also be specified in the config file.  

Required to run : matplotlib, PLINK, Snakemake 

Required input:
 * plink bim, bam, bed files
 * id list .txt file including all but the paternal sequence ids (IDlist.txt - should be added to the data folder.)

## Preprocess Data


### QC 

Filtering the plink-format data using --maf and --geno, followed by checking for mendelian errors --me 

### Recode Genes to .ped format

To find the paternal haplotypes, a .ped format file is required. 

### Rename Genes

Genes are renamed to a standard format. 

## Finding Common Mutations in Halfsibs

### Find paternal halfsib haplotypes

Paternal haplotypes of the half-sibs are determined here using maternal, paternal and offspring data. 

### Find markers common in all halfsib paternal haplotypes

## Annotating Homozygous Loci

To narrow down the areas of interest, homozygous loci are determined. 

### Finding Homozygous Loci using Plink

Homozygous loci are determined using plink.



## Visualising

### Reformatting the outputs for graphing

The homozygous loci output is reformatted to a ; separated csv file.
The same is done for the variant file. 

### Graph

In the final step, the data is visualised as chromosomes using matplotlib. 

