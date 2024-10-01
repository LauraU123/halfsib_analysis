#! /bin/bash

module load snakemake

snakemake --cores all --configfile config/medium_density.yaml -s run.smk --jobs 70 --latency-wait 30 -F --cluster "sbatch -ppibu_el8  --mem={resources.mem} --cpus-per-task={threads} --time={resources.time}"