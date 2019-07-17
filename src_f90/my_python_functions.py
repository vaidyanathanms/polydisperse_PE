# Generic Function Definitions
# Version_1: V_Mar_26_2019

import numpy
import os
import shutil
import subprocess
import sys
import glob
import re

def cpy_main_files(maindir,destdir,fylname):

    srcfyl = maindir + '/' + fylname
    desfyl = destdir + '/' + fylname
    shutil.copy2(srcfyl, desfyl)

def init_pdi_write(freepdi, freemw,freechains,\
                   graftpdi,graftmw,graftch,destdir):

    pdi_fyl = destdir + '/init_pdi.txt'
    finit   = open(pdi_fyl,'w')
    finit.write('%s\t %s\n' %('free_data', '#pdi mw nchains'))
    finit.write('%g\t %g\t %g\n' %(freepdi,freemw,freechains))
    finit.write('%s\t %s\n' %('graft_data', '#pdi mw nchains'))
    finit.write('%g\t %g\t %g\n' %(graftpdi,graftmw,graftch))
    finit.close()


def compile_and_run_pdi(destdir):
    
    os.chdir(destdir)
    if not os.path.exists('init_pdi.txt'):
        print('init_pdi.txt not found')
        return

    if not os.path.exists('SZDist2.f90'):
        print('SZDist2.f90 not found')
        return
        
    # Generate PDI data

    print( "Running FORTRAN script for generating New Datafile")
    subprocess.call(["ifort","-r8","-check","-traceback",
                     "SZDist2.f90","-o","pdiinp.o"])

    subprocess.call(["./pdiinp.o"])


def check_pdi_files(destdir,pdi_files,flagcheck):

    os.chdir(destdir)
    for fyl in pdi_files:
        if not os.path.exists(fyl):
            print(fyl, 'not found')            
            flagcheck = -1


def create_paramfyl_for_datafyl(destdir,inpfyles,nch_free,archstr\
                                ,casenum,mw_free,nch_graft,mw_graft,\
                                nsalt,f_charge,ntail,iarch,pdi_files):
    
    os.chdir(destdir)
    datafyle = "PEdata_"+str(free_chains[ifree])+"_" \
               +fylstr+"_"+str(caselen+1)+".dat"
    logout = "log_" + dataname

    fr  = open(inpfyles[0],'r')
    fw  = open(inpfyles[1],'w')

    fid = fr.read().replace("py_dataname",datafyle).\
          replace("py_freepdi_dataname",pdi_files[0]).\
          replace("py_graft_pdi_dataname",pdi_files[1]).\
          replace("py_numfree_chains", str(nch_free)).\
          replace("py_avgmw_free",str(mw_free)).\
          replace("py_numgraft_chains",str(nch_graft)).\
          replace("py_avgmw_graft", str(mw_graft)).\
          replace("py_numtail_mons",str(ntail)).\
          replace("py_saltonespecies",str(nsalt)).\
          replace("py_chargfrac_perchain",str(f_charge)).\
          replace("py_arch",str(iarch)).\
          replace("py_log_fylename",logout)
    fw.write(fid)
    fw.close()
    fr.close()

    outfyle = launch_fyl


def compile_and_run_inpgenfyles(launch_fyle,destdir):

    if not os.path.exists('ran_numbers.f90'):
        print('ran_numbers.f90 not found')
        return
    if not os.path.exists('lmp_params.f90'):
        print('lmp_params.f90 not found')
        return
    if not os.path.exists('lammps_inp.f90'):
        print('lammps_inp.f90 not found')
        return

    os.chdir(destdir)
    subprocess.call(['ran_numbers.f90','lmp_params.f90','lammps_inp.f90'\
                     '-o','inpgen.o'])
    subprocess.call(['./inpgen.o',launch_fyle])

    
def clean_backup_initfiles(f90_files,pdi_files,par_files,destdir):
    
    initdir = destdir + '/init_files'
    if not os.path.isdir(initdir):
        os.mkdir(initdir)

    for fyl in f90_files:
        if os.path.exists(fyl):
            cpy_main_files(destdir,initdir,fyl)
            os.remove(fyl)

    for fyl in pdi_files:
        if os.path.exists(fyl):
            cpy_main_files(destdir,initdir,fyl)
            os.remove(fyl)

    for fyl in par_files:
        if os.path.exists(fyl):
            cpy_main_files(destdir,initdir,fyl)
            os.remove(fyl)

    files = glob.glob('init_*')
    for fyl in files:
        if not os.path.isdir(fyl):
            cpy_main_files(destdir,initdir,fyl)
            os.remove(fyl)

    files = glob.glob(destdir +'/*var*')
    for f in files:
        os.remove(f)

    files = glob.glob(destdir +'/*.mod')
    for f in files:
        os.remove(f)

    files = glob.glob(destdir +'/*.o')
    for f in files:
        os.remove(f)

