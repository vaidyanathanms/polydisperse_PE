# Generic Function Definitions
# Version_1: V_Mar_26_2019

import numpy
import os
import shutil
import subprocess
import sys
import glob
import re

def cpy_main_files(dum_maindir,dum_destdir,fylname):

    srcfyl = dum_maindir + '/' + fylname

    if not os.path.exists(srcfyl):
        print('ERROR', fylname, 'not found')
        return

    desfyl = dum_destdir + '/' + fylname
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

    print( "Running FORTRAN script for generating PDI files")
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
    datafyle = "PEdata_"+str(nch_free)+"_" +archstr+ \
               "_"+str(casenum+1)+".dat"
    logout = "log_" + datafyle

    fr  = open(inpfyles[0],'r')
    fw  = open(inpfyles[1],'w')

    fid = fr.read().replace("py_dataname",datafyle).\
          replace("py_freepdi_dataname",pdi_files[0]).\
          replace("py_graftpdi_dataname",pdi_files[1]).\
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

    return datafyle

def compile_and_run_inpgenfyles(launch_fyle,destdir):

    os.chdir(destdir)
    if not os.path.exists('ran_numbers.f90'):
        print('ran_numbers.f90 not found')
        return
    if not os.path.exists('lmp_params.f90'):
        print('lmp_params.f90 not found')
        return
    if not os.path.exists('lammps_inp.f90'):
        print('lammps_inp.f90 not found')
        return

    subprocess.call(['ifort','-r8','-qopenmp','-mkl','-check','-traceback',\
                     'ran_numbers.f90','lmp_params.f90','lammps_inp.f90',\
                     '-o','inpgen.o'])
    subprocess.call(['./inpgen.o',launch_fyle])


def edit_generate_input_lmp_files(lmp_infyle,lmp_datafyle):

    if not os.path.exists(lmp_infyle):
        print('ERROR: ', lmp_infyle, 'not found')
        return

    fr  = open(lmp_infyle,'r')
    fw  = open('in.init','w')
    fid = fr.read().replace("py_dataname",lmp_datafyle)
    fw.write(fid)
    fw.close()
    fr.close()
    

def run_lammps(nch_free,pdifree,casenum,dirstr,inpjob,outjob):

    if not os.path.exists(inpjob):
        print('ERROR: ', inpjob,'not found')
        return
    
    jobstr = "job_" + str(nch_free) + "_" + str(pdifree) + "_" \
             + str(casenum+1) + "_" + dirstr
    fr  = open(inpjob,'r')
    fw  = open(outjob,'w')
    fid = fr.read().replace("py_jobname",jobstr)
    fw.write(fid)
    fw.close()
    fr.close()

    subprocess.call(["qsub", outjob])
    
    
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

    files = glob.glob('*var*')
    for fyl in files:
        if os.path.exists(fyl):
            cpy_main_files(destdir,initdir,fyl)
            os.remove(fyl)

    files = glob.glob('*.o')
    for fyl in files:
        if os.path.exists(fyl):
            cpy_main_files(destdir,initdir,fyl)
            os.remove(fyl)
    
    files = glob.glob(destdir +'/*.mod')
    for fyl in files:
        os.remove(fyl)


def compile_anafiles():

    if not os.path.exists('pe_params.f90'):
        print('ERROR: pe_params.f90 not found')
        return

    if not os.path.exists('pe_analyze.f90'):
        print('ERROR: pe_analyze.f90 not found')
        return


    subprocess.call(['ifort','-r8','-qopenmp','-mkl',\
                     '-check','-traceback','pe_params.f90'\
                     'pe_analyze.f90','-o','anainp.o'])

def find_datafyle(nch_free,archstr,caseum,destdir):

    os.chdir(destdir)
    curr_dir = os.getcwd()
    datafyle = "PEdata_"+str(nch_free)+"_" +archstr+ \
               "_"+str(casenum+1)+".dat"

    if not os.path.exists(datafyle):
        print ("Trying to make datafile from restart files")
        restart_fyles = glob.glob('archival_*')
        
        if restart_fyles = []:
            print("ERROR: No restart files found")
            return

        if not os.path.exists('lmp_mesabi'):
        
            src_lmp = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
            destfyle = curr_dir + '/lmp_mesabi'
            shutil.copy2(src_lmp,destfyle)

        subprocess.call(['./lmp_mesabi','-np','48','./lmp_mesabi','-r',restart_fyles[0],datafyle])
        
    return datafyle


def find_latest_trajfyle(pref,destdir):
    
    os.chdir(destdir)
    traj_arr = glob.glob(pref)
    latest_fyle = max(traj_arr,key = os.path.getctime)
    return latest_fyle


def edit_generate_anainp_files(inpdata,inptraj,nch_tot,nch_free,\
                               nch_graft):
    
    if not os.path.exists('anainp_var.txt'):
        print('ERROR: pe_params not found')
        return

    fr  = open('anainp_var.txt','r')
    fw  = open('anainp.txt','w')

    datafyle = re.split('/',inpdata)
    dataname = datafyle[len(datafyle)-1]

    trajfyle = re.split('/',inptraj)
    trajname = trajfyle[len(trajfyle)-1]
    

    fid = fr.read().replace("py_datafyl",dataname).\
          replace("py_trajfyl",trajname).\
          replace("py_ntotchains",str(nch_tot)).\
          replace("py_nfrchains", str(nch_free)).\
          replace("py_ngrchains",str(nch_graft))
    fw.write(fid)
    fw.close()
    fr.close()

def run_lammps(nch_free,pdifree,casenum,dirstr,inpjob,outjob):

    if not os.path.exists(inpjob):
        print('ERROR: ', inpjob,'not found')
        return
    
    jobstr = "ana_" + str(nch_free) + "_" + str(pdifree) + "_" \
             + str(casenum+1) + "_" + dirstr
    fr  = open(inpjob,'r')
    fw  = open(outjob,'w')
    fid = fr.read().replace("py_jobname",jobstr)
    fw.write(fid)
    fw.close()
    fr.close()

    subprocess.call(["qsub", outjob])
    
