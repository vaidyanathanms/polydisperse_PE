#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=2g 
#SBATCH --job-name=compress
#SBATCH --mail-user=vsethura@umn.edu 
echo "job_start"
tar cvzf n_150_feb22_2022.tar.gz n_150
wait
tar cvzf n_128_feb22_2022.tar.gz n_128
wait
tar cvzf n_64_feb22_2022.tar.gz n_64
wait
tar cvzf n_32_feb22_2022.tar.gz n_32
wait
tar cvzf n_16_feb22_2022.tar.gz n_16
echo "done"

