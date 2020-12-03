#!/bin/bash -l
#PBS -l walltime=03:00:00,nodes=1:ppn=1,pmem=2580mb
#PBS -m abe
#PBS -N py_jobname
#PBS -M vsethura@umn.edu
cd ${PBS_O_WORKDIR}
echo job_start
export OMP_NUM_THREADS=1
./anainp.o anainp_pyfylval.txt


wait

mkdir results_py_freech_py_pdifree_py_caselen
mv adsfrac* results_py_freech_py_pdifree_py_caselen/
mv tether* results_py_freech_py_pdifree_py_caselen/
mv chainadsval* results_py_freech_py_pdifree_py_caselen/
mv log.config* results_py_freech_py_pdifree_py_caselen/
mv dens_* results_py_freech_py_pdifree_py_caselen/
mv chdens_* results_py_freech_py_pdifree_py_caselen/
mv PErdf_* results_py_freech_py_pdifree_py_caselen/
mv chgrpdens_* results_py_freech_py_pdifree_py_caselen/
mv grpdens_* results_py_freech_py_pdifree_py_caselen/
mv polydens_* results_py_freech_py_pdifree_py_caselen/
