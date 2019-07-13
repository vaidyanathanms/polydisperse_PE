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


def check_pdi_files(destdir,flagcheck):

    os.chdir(destdir)
    if not os.path.exists('FreeChains.dat'):
        print('FreeChains.dat not found')
        flagcheck = -1
        return

    if not os.path.exists('GraftChains.dat'):
        print('GraftChains.dat not found')
        flagcheck = -1
        return
        

def compute_total_particles(destdir,freech,graftch,tailmons,nsalt,\
                            fcharge,free_mons,graft_mons,totalmons):
    os.chdir(destdir)
    free_fyl = destdir + '/FreeChains.dat'
    ffree = open(free_fyl,'r')
    lyne = ffree.readline()
    dat_lyne = lyne.split()
    if int(dat_lyne[0]) != freech:
        print("Free chains do not match",int(dat_lyne[0]),freech)
    
    free_mons = int(dat_lyne[1])
    ffree.close()

    graft_fyl = destdir + '/GraftChains.dat'
    fgraft = open(graft_fyl,'r')
    lyne = fgraft.readline()
    dat_lyne = lyne.split()
    if int(dat_lyne[0]) != graftch:
        print("Graft chains do not match",int(dat_lyne[0]),graftch)
    
    graft_mons = int(dat_lyne[1])
    fgraft.close()
    
    nfree_cntr = int(fcharge*free_mons)
    ngraft_cntr = int(fcharge*(graft_mons - \
                                graftch*tailmons))
    ntotal = nfree_cntr + ngraft_cntr + free_mons + \
                         graft_mons + 2*nsalt
    
    fout_main = open('init_alldata.txt','w')
    fout_main.write('%s\t %g\n' %('Free chains', freech))
    fout_main.write('%s\t %g\n' %('Graft chains',graftch))
    fout_main.write('%s\t %g\n' %('Total freemons', free_mons))
    fout_main.write('%s\t %g\n' %('Total graftmons',graft_mons))
    fout_main.write('%s\t %g\n' %('Free counter', nfree_cntr))
    fout_main.write('%s\t %g\n' %('Graft counter',\
                                  ngraft_cntr))
    fout_main.write('%s\t %g\n' %('Net salt (sum both)',\
                                  2*nsalt))
    fout_main.write('%s\t %g\n' %('Total particle',ntotal))
    fout_main.close()

