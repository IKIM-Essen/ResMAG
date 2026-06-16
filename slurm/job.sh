#!/usr/bin/env sh
#SBATCH --nodes 1
#SBATCH --exclusive
#SBATCH --time 48:00:00

#SBATCH --output=slurm/out/slurm-%j.out

echo "my job ran on $(date)"
hostname
echo "job id: $SLURM_JOB_ID"  # this variable is defined by Slurm
echo "working directory: $(pwd)"


snakemake --unlock

echo "starting snakemake workflow"
snakemake \
    --cores 64 \
    --profile /groups/ds/metagenomes/ResMAG/config/profile/ \
    --singularity-args "--bind /groups/ds/databases_refGenomes"

echo "cleaning snakemake temp output"
snakemake --delete-temp-output