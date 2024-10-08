#! /bin/bash
module load snakemake

snakemake --cores all cleanup -F --configfile config/high_density.yaml -s run.smk --jobs 15 --cluster "sbatch -ppibu_el8  --mem={resources.mem} --cpus-per-task={threads} --time={resources.time}" --latency-wait 30