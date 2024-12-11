#! /bin/bash
module load snakemake

while getopts "c:" opt; do
    case ${opt} in
        c) CONFIG_FILE="$OPTARG";;
        *) echo "Usage: $0 -c configfile"; exit 1 ;;
    esac
done

snakemake --cores all -F --configfile "${CONFIGFILE:-config/medium_density.yaml}" -s run.smk --jobs 15 #--cluster "sbatch -ppibu_el8  --mem={resources.mem} --cpus-per-task={threads} --time={resources.time}" --latency-wait 30