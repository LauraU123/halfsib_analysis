# Halfsib Linkage

This tool can be used to find common haplotype regions in a dataset of half-siblings and their parents.  
It requires PLINK and Snakemake to run. 

## Inputs

The input should be in plink-friendly format - .bim, .bed and .fam.


![dag.png](https://github.com/LauraU123/halfsib_analysis/blob/parallelized_rules/config/dag.png)


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

## Finding Common Variants in Halfsibs

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

