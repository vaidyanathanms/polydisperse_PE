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
from my_python_functions import create_paramfyl_for_datafyl
from my_python_functions import compile_and_run_inpgenfyles
from my_python_functions import edit_generate_input_lmp_files
from my_python_functions import run_lammps
from my_python_functions import clean_backup_initfiles

#---------input flags------------------------------------------
#0-initial run  1- production
restart = 0 # For restarting from given configurations


#---------input details----------------------------------------
free_chains  = [32]#,80,32,48]
free_avg_mw  = 30
graft_chains = 64
graft_avg_mw = 35 
tail_mons    = 5
nsalt        = 510
f_charge     = 0.5
archarr      = [1]#,2,3,4]
ncases_pdi   = 5
pdi_free     = 1.2
pdi_graft    = 1.0

#--------file_lists--------------------------------------------

f90_files = ['ran_numbers.f90','lammps_inp.f90','lmp_params.f90'\
             ,'SZDist2.f90','init_pdi.txt','polyinp_var.dat']
lmp_files = ['in.longrun','in.init_var','in.run1','jobmain_var.sh']
lmp_long  = ['in.longrun','jobmain_long_var.sh']
pdi_files = ['FreeChains.dat','GraftChains.dat'] #Order: freepdi,graftpdi
par_files = ['polyinp_var.dat','polyinp.txt'] #Order: varfile,oufyle

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

        workdir_freepdi = workdir_arch + '/pdi_free_' + str(pdi_free)
        if not os.path.isdir(workdir_freepdi):
            os.mkdir(workdir_freepdi)

        workdir_graftpdi = workdir_freepdi + '/pdi_graft_' + str(pdi_graft)
        if not os.path.isdir(workdir_graftpdi):
            os.mkdir(workdir_graftpdi)

        for caselen in range(ncases_pdi):

            workdir_subpdi = workdir_graftpdi + '/Case_' + str(caselen+1)
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
                init_pdi_write(pdi_free,free_avg_mw,free_chains[ifree],\
                               pdi_graft,graft_avg_mw,graft_chains,destdir)
                compile_and_run_pdi(destdir)

                pdiflag = 1
                check_pdi_files(destdir,pdi_files,pdiflag)

                if pdiflag == -1:
                    print("Check PDI files")
                    continue

                    
                #----Generate input files-----
                print( "Copy Successful - Generating Input Files")
                
                lmp_data_fyle = create_paramfyl_for_datafyl(destdir,\
                                                      par_files,free_chains[ifree]\
                                                      ,fylstr,caselen,free_avg_mw,\
                                                      graft_chains,graft_avg_mw,nsalt,\
                                                      f_charge,tail_mons,\
                                                      archarr[iarch],pdi_files)

                compile_and_run_inpgenfyles(par_files[1],destdir)            

                #---Run LAMMPS files-------------

                edit_generate_input_lmp_files('in.init_var',lmp_data_fyle)
#                run_lammps(free_chains[ifree],pdifree,caselen,fylstr,\
#                           'jobmain_var.sh','jobmain.sh')
                
                #----Copy/Backup initial files---
                clean_backup_initfiles(f90_files,pdi_files,par_files,destdir)

                os.chdir(maindir)

            else:

                os.chdir(workdir_subpdi)
                destdir = os.getcwd()

                if not os.path.isdir(workdir_subpdi):
                    print( workdir_subpdi, "not found")
                    continue


                archfiles = destdir + '/archival*'
                list_of_files = glob.glob(archfiles)

                if not list_of_files:
                    print("No archival files found in ", destdir)
                    continue


                for fyllist in range(len(lmp_long)):
                    cpy_main_files(lmp_dir,destdir,lmp_long[fyllist])

                fylename = destdir + '/lmp_mesabi'
                if not fylename: 
                    srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
                    desfyl = destdir + '/lmp_mesabi'
                    shutil.copy2(srcfyl, desfyl)

#                run_lammps(free_chains[ifree],pdifree,caselen,fylstr,\
#                           'jobmain_long_var.sh','jobmain_long.sh')


                os.chdir(maindir)
	 
