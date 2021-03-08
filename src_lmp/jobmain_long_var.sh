#!/bin/bash -l
#PBS -l walltime=py_hours:00:00,nodes=py_nodes:ppn=py_procs,pmem=2580mb
#PBS -m abe
#PBS -N py_jobname
#PBS -M vsethura@umn.edu
cd ${PBS_O_WORKDIR}
echo job_start
mpirun -np py_totnp ./lmp_mesabi -in in.longrun -e screen
