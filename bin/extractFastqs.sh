#!/bin/bash
#SBATCH --job-name=ex_fqs    # Job name
#SBATCH --ntasks=1                    # one task
#SBATCH --cpus-per-task=128       # Run on 128 CPUs
#SBATCH --output=extract_fastqs.log   # Standard output and error log

cd $SLURM_SUBMIT_DIR
echo $$
pwd; hostname; date

bash -c /home/jarvis/bin/bcl2fastq_wrapper.sh $SLURM_CPUS_PER_TASK

date
