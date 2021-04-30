#!/bin/bash -l
#SBATCH --time=py_hours:00:00
#SBATCH -N py_nodes
#SBATCH -n py_procs
#SBATCH --mem=2g 
#SBATCH --job-name=py_jobname
#SBATCH --mail-user=vsethura@umn.edu 
cd ${SLURM_SUBMIT_DIR}
echo job_start
mpirun -np py_totnp ./lmp_mesabi -in in.init -e screen
wait
mpirun -np py_totnp ./lmp_mesabi -in in.run1 -e screen
wait
mpirun -np py_totnp ./lmp_mesabi -in in.run2 -e screen
wait
mpirun -np py_totnp ./lmp_mesabi -in in.longrun -e screen
