#!/bin/bash -l
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hamsh010@umn.edu
#SBATCH --job-name calc_afreqs 
#SBATCH -o slurm_logs/%j.calc_afreqs.out
#SBATCH -e slurm_logs/%j.calc_afreqs.err
#SBATCH -A fried255
#SBATCH -p agsmall,aglarge,amdsmall,amdlarge

set -e

conda activate snakemake
cd $SLURM_SUBMIT_DIR

VCF_DIR=/panfs/jay/groups/0/fried255/fried255/working/pipeline/UU_Cfam_GSD_1.0_ROSY

snakemake -s freq_analysis.smk \
    --use-singularity \
    --singularity-args "-B $PWD,$VCF_DIR" \
    --profile slurm.afreq \
    --configfile config.yaml \
    --keep-going \
    --rerun-triggers mtime

