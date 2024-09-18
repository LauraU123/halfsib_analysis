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

To narrow down the areas of interest, homozygous loci are determined in this step. 

### Finding Homozygous Loci using Plink



### Reformatting the output for graphing

The homozygous loci output is reformatted to a ; separated csv file for the next step. 

## Graphing

In the final step, the data is visualised as chromosomes using matplotlib. 





### Step 3: Finding the Paternal Haplotypes and Annotating Homozygous Loci

Input: file with haplotype data for parents and offspring 

The first step includes finding the common paternal haplotypes based on

