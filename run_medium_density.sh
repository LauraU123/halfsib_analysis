#! /bin/bash

module load snakemake

snakemake --cores all -F  --configfile config/medium_density.yaml -s run.smk --jobs 70 --latency-wait 30 --cluster "sbatch -ppibu_el8  --output=/dev/null --mem={resources.mem} --cpus-per-task={threads} --time={resources.time}"