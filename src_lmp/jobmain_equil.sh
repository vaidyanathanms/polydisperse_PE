#!/bin/bash -l
#PBS -l walltime=24:00:00,nodes=2:ppn=24,pmem=2580mb
#PBS -m abe
#PBS -N py_jobname
#PBS -M vsethura@umn.edu
cd ${PBS_O_WORKDIR}
echo job_start
mpirun -np py_totnp ./lmp_mesabi -in in.run2 -e screen
wait
mpirun -np py_totnp ./lmp_mesabi -in in.longrun -e screen
