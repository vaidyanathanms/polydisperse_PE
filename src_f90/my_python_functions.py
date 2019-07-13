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

def init_pdi_write(freechains, freemw, freepdi,\
                   graftch,graftmw,graftpdi,destdir):

    pdi_fyl = destdir + '/init_pdi.txt'
    finit   = open(pdi_fyl,'w')
    finit.write('%s\t %s' %('free_data', '#pdi mw nchains'))
    finit.write('%d\t %d\t %g\n' %(freepdi,freemw,freechains))
    finit.write('%s\t %s' %('graft_data', '#pdi mw nchains'))
    finit.write('%d\t %d\t %g\n' %(graftpdi,graftmw,graftch))
    finit.close()


def compile_and_run_pdi(destdir):
    
    os.chdir(destdir)
    if not os.path.exists('init_pdi.txt'):
        print('init_pdi.txt not found')
        break

    if not os.path.exists('SZDist2.f90'):
        print('SZDist2.f00 not found')
        break
        
    # Generate PDI data

    print( "Running FORTRAN script for generating New Datafile")
    subprocess.call(["ifort","-r8","-check","-traceback",
                     "SZDist2.f90","-o","pdiinp.o"])

    subprocess.call(["./pdiinp.o"])


def check_pdifiles(destdir,flagcheck):

    os.chdir(destdir)
    if not os.path.exists('FreeChains.dat'):
        print('FreeChains.dat not found')
        flagcheck = -1
        break

    if not os.path.exists('GraftChains.dat'):
        print('GraftChains.dat not found')
        flagcheck = -1
        break
        


    


