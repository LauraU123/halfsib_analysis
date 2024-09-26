# Workflow


## 1. Preprocessing Data


### QC 

Filtering the plink-format data using --maf and --geno, followed by checking for mendelian errors --me.
The output is a filtered set of .bim, bed and .fam files.

### Recoding Genes to .ped format

The filtered files are converted to a human-readable format with plink. 

### Renaming Genes

Genes are renamed to a standardised format: CHR_LOCATION.
This ensures locus data is retained for later steps.

## 2. Finding Common Paternally Inherited Locations in Halfsibs

### Determining Paternal Haplotypes in Halfsibs

Paternal haplotypes of the half-sibs are determined here using maternal, paternal and offspring data. 
Output files are written for each chromosome in the format "hafsibID,sequence"
e.g.
AB_123, AAAABBBBAAB...
AB_124, AAABBBBBAAB...
AB_125, AAABBBBBBBA...

### Finding Common Markers in Halfsib Paternal Haplotypes in each Chromosome + Merging

Common markers are determined in all 


## 3. Annotating Homozygous Loci

To narrow down the areas of interest, homozygous loci are determined. 

### Finding Homozygous Loci using Plink

Homozygous loci are determined using plink.



## 4. Visualising Data

### Reformatting the outputs for graphing

The homozygous loci output is reformatted to a ; separated csv file.
The same is done for the variant file. 

### Graph

In the final step, the data is visualised as chromosomes using matplotlib. 


![dag.png](https://github.com/LauraU123/halfsib_analysis/blob/parallelized_rules/config/dag.png)