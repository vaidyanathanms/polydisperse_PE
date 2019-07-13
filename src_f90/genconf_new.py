import numpy
import os
import shutil
import subprocess
import sys
import glob
from subprocess import call


#---------mypython functions------------------------------

from my_python_functions import cpy_main_files
from my_python_functions import init_pdi_write
from my_python_functions import compile_and_run_pdi
from my_python_functions import check_pdi_files
from my_python_functions import compute_total_particles

#---------input flags------------------------------------------
#0-initial run  1- production
restart = 0 # For restarting from given configurations


#---------input details----------------------------------------
free_chains  = [32]#,80,32,48]
free_avg_mw  = 30
graft_chains = 64
graft_avg_mw = 32 
tail_mons    = 2
nsalt        = 510
f_charge     = 0.5
archarr      = [1]#,2,3,4]
ncases_pdi   = 5
pdi_free     = 1.2
pdi_graft    = 1.0

#--------file_lists--------------------------------------------

f90_files = ['ran_numbers.f90','lammps_inp.f90','lmp_params.f90'
             ,'SZDist2.f90','init_pdi.txt']
lmp_files = ['in.longrun','in.init']

#---------directory info---------------------------------------
maindir = os.getcwd()
scratchdir = '/scratch.global/vaidya/'
lmpdir = '/home/dorfmank/vsethura/allfiles/files_polydisperse/src_lmp'

   
#Main analysis
for ifree in range(len(free_chains)):
    
    print( "Free Number of Chains: ", free_chains[ifree])
    workdir1 = scratchdir + 'polydisp_pe'
    
    if not os.path.isdir(workdir1):
        os.mkdir(workdir1)

    workdir2 = workdir1 + '/n_' + str(free_chains[ifree])
	
    if not os.path.isdir(workdir2):
        os.mkdir(workdir2)
        
    for iarch in range(len(archarr)):
             
        if archarr[iarch] == 1:
            print( "Archval: Block_Block")
            dirstr = 'bl_bl'
            fylstr = 'block_block'
        elif archarr[iarch] == 2:
            print( "Archval: Block_Alter")
            dirstr = 'bl_al'
            fylstr = 'block_alter'
        elif archarr[iarch] == 3:
            print( "Archval: Alter_Block")
            dirstr = 'al_bl'
            fylstr = 'alter_block'
        elif archarr[iarch] == 4:
            print( "Archval: Alter_Alter")
            dirstr = 'al_al'
            fylstr = 'alter_alter'
        else:
            print( "Unknown Architecture")
            
        workdir_arch = workdir2 + '/' + dirstr

        if not os.path.isdir(workdir_arch):
            os.mkdir(workdir_arch)

        workdir_headpdi = workdir2 + '/pdi_' + str(pdi_free)
        if not os.path.isdir(workdir_headpdi):
            os.mkdir(workdir_headpdi)

        for caselen in range(ncases_pdi):

            workdir_subpdi = workdir_headpdi + '/Case_' + str(caselen+1)
            if not os.path.isdir(workdir_subpdi):
                os.mkdir(workdir_subpdi)

            if restart == 0:

                os.chdir(workdir_subpdi)
                destdir = os.getcwd()

                print( "Starting new case for ", free_chains[ifree],dirstr)

                #---Copying files------
                print( "Current Dir ", destdir)
                print( "Copying Files")
                
                for fyllist in range(len(f90_files)):
                    cpy_main_files(maindir,destdir,f90_files[fyllist])

                for fyllist in range(len(lmp_files)):
                    cpy_main_files(lmpdir,destdir,lmp_files[fyllist])

                srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
                desfyl = destdir + '/lmp_mesabi'
                shutil.copy2(srcfyl, desfyl)

                #----Run PDI code------
                print("Running PDI calculation file")
                tot_chains = free_chains[ifree] + graft_chains
                init_pdi_write(free_chains[ifree],free_avg_mw,pdi_graft,\
                               graft_chains,graft_avg_mw,pdi_graft,destdir)
                compile_and_run_pdi(destdir)

                pdiflag = 1
                check_pdi_files(destdir,pdiflag)

                if pdiflag == -1:
                    print("Check PDI files")
                    continue

                    
                nfree_mons = 0
                ngraft_mons = 0
                ntotal = 0
                compute_total_particles(destdir,free_chains[ifree],\
                                        graft_chains,tail_mons,nsalt,f_charge, 
                                        nfree_mons,ngraft_mons,ntotal)

                #----Generate input files-----

                print( "Copy Successful - Generating Input Files")

                dataname = "PEdata_new_"+str(free_chains[ifree])+\
                           "_"+fylstr+".dat"


                initdir = destdir + '/init_files'
                if not os.path.isdir(initdir):
                    os.mkdir(initdir)


                for fyl in f90_files:
                    if os.path.exists(fyl):
                        cpy_main_files(destdir,initdir,fyl)
                        os.remove(fyl)

                files = glob.glob(destdir +'/init*')
                for fyl in files:
                    cpy_main_files(destdir,initdir,fyl)
                    os.remove(f)
                

                files = glob.glob(destdir +'/*var*')
                for f in files:
                    os.remove(f)

                files = glob.glob(destdir +'/*.mod')
                for f in files:
                    os.remove(f)

                files = glob.glob(destdir +'/*.o')
                for f in files:
                    os.remove(f)

                if not os.path.isfile(newdataname):
                    print(newdataname, "not found")
                    continue


                print( "Submitting Jobs..")

                subprocess.call(["qsub","jobswitch.sh"])

                os.chdir(maindir)


            else:

                os.chdir(workdir3)
                destdir = os.getcwd()

                if not os.path.isdir(workdir3):
                    print( workdir3, "not found")
                    continue


                archfiles = destdir + '/archival*'
                list_of_files = glob.glob(archfiles)

                if not list_of_files:
                    print("No archival files found in ", destdir)
                    continue

                srcfyl = lmpdir + '/in.longrun'
                desfyl = destdir + '/in.longrun'
                shutil.copy2(srcfyl, desfyl)

                srcfyl = lmpdir + '/jobmain2.sh'
                desfyl = destdir + '/jobmain2.sh'
                shutil.copy2(srcfyl, desfyl)

                fylename = destdir + '/lmp_mesabi'

                if not fylename: 
                    srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
                    desfyl = destdir + '/lmp_mesabi'
                    shutil.copy2(srcfyl, desfyl)

                print( "Submitting Jobs..")

                subprocess.call(["qsub","jobmain2.sh"])

                os.chdir(maindir)
	 
