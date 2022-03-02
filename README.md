# Charge Sequence-Specific Polydisperse Polyelectrolytes

Adsorption of polydisperse polyelectrolytes onto grafted polyelectrolyte systems with polydispersity

## Code Details

- LAMMPS codes: Aug-11-2017 version (see src_lmp) 
- FORTRAN codes: Requires ifort compiler (see src_f90)
- Python codes: Requires V3.0+ (see src_f90)
- MATLAB codes: Requires 2019b+ (see src_matlab)

NOTE1: For the Python codes, check the directory paths and rename accordingly. The scratch (where the codes will be running) will be DIFFERENT and needs to be changed accordingly. Specifically, one needs to set "jobdir", "scratchdir" and "lmpdir" whenever necessary. Please do a grep command to find the points where they are defined. 

NOTE2: The path to the executable for LAMMPS - lmp_mesabi also needs to be changed. Please do a grep for "lmp_mesabi" to see where they needs to be changed. 

## File Uses

### FORTRAN/Python Codes

- newparams_genconf2.py - To generate the input configurations. See "Inputs" within the code to see the kind of inputs required.
- ana.py - To compute the analysis of the data generated.
- copy_to_analyze.py - To copy all the analysis output files to a new directory for further analysis with MATLAB. 
- my_python_functions.py - Supplemental functions for running python codes.

NOTE1: The FORTRAN files are called within the python files to generate the SZ distribution (SZdist2.f90) and generating the LAMMPS inputs (params_inpconf.f90, inpconf_generator.f90, ran_numbers.f90). Analysis files are called using ana.py and the dependencies are pe_analyze.f90 and pe_params.f90

NOTE2: Some of the FORTRAN and input files may have a *_var* version which is used by the python codes to replace the inputs/parameters with user input versions. DO NOT delete this versions if the user is using the Python codes

### MATLAB Codes
See README.txt in src_matlab directory to see usages on MATLAB codes.

### LAMMPS Codes
- in.init/in.init_var - Initialize system and move use NVE/limit
- in.run1 - First equilibration step
- in.run2 - Second equilibration step/Production cycle
= in.longrun - Production cycle

## Order of Usage

- Use newparams_genconf2.py to create and run LAMMPS
- Use ana.py to analyze data
- Use copy_to_analyze.py to consolidate outputs
- See MATLAB folder to postprocess data.
