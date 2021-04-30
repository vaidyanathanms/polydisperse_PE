#!/bin/bash -l
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2g
#SBATCH --tmp=2g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=vsethura@umn.edu 
cd pyoutdir
echo ${PWD}
module load intel 
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
