# Paternal Halfsibling Linkage Analysis tool

This tool can be used for finding common haplotypes in paternal halfsiblings.
Chromosome number of the species of interest should be specified in the config file. 

## Preprocess Data

### Recode to Plink Format

The raw data should first be converted to a plink-friendly format.
This is done in the first step. 

### QC 

Filtering the plink-format data using --maf and --geno

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


## Graphing

Visualising the data using R





### Step 3: Finding the Paternal Haplotypes and Annotating Homozygous Loci

Input: file with haplotype data for parents and offspring 

The first step includes finding the common paternal haplotypes based on

