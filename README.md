# Halfsib Linkage

This tool can be used to find common haplotype regions in a dataset of half-siblings and their parents.  
It requires PLINK and Snakemake to run. 

## Data Density Dependent Parameters

Different config files should be used depending on the input data density.
For WGS data, config/high_density.yaml should be used.

For less dense data, config/medium_density.yaml can be used.
These can also be manually adjusted. 

## Input Dataset


The input should be in plink-friendly format - .bim, .bed and .fam.

## Output

The output consists of a table and graph. The graph shows homozygous regions, linked regions and overlaps


