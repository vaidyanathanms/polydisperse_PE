!----Version - April 30 2019----------------------------------------
! Input file to generate LAMMPS file for polyelectrolyte simulations
! Use in conjunction with inpconf_generator.f90 and ran_numbers.f90
! Included polydispersity
! Change the required parameters
! arch details: "arch=number:polycation-polyanion 
! arch = 1:block-block;arch=2:block-alt;arch=3:alt-block;arch=4:alt
! -alt.
MODULE PARAMS

  USE RAN_NUMBERS

  IMPLICIT NONE

! Parameter data for creating the data file

  INTEGER :: nbulk_chains,ngraft_chains 
  INTEGER :: N_salt
  INTEGER :: default_dim = 1 !If 1=>53*53*120.Or else use  !box.dat
  REAL    :: charge_frac = 0.5
  INTEGER :: ncntr_graft!
  !N_graft*INT((M_graft-tail_graft)*charge_frac)
  INTEGER :: ncntr_free!   = N*INT(M*charge_frac)
  REAL    :: graft_dist
  REAL    :: pd_graft_val, pd_bulk_val
  INTEGER :: arch
  INTEGER :: mw_graft, mw_bulk

! Box/Particle details

  REAL :: boxl_x, boxl_y, boxl_z
  REAL :: volbox, density
  INTEGER :: totpart
  INTEGER :: n_cntr_ions, nchains, npolyatoms

! Flags for creating topology

  INTEGER, PARAMETER :: bondtype = 1
  INTEGER, PARAMETER :: angltype = 1

! Flags for creating files
  INTEGER, PARAMETER :: outfile = 300, datawrite=320, anaread=330
  INTEGER, PARAMETER :: poly_graft_read = 340, poly_bulk_read = 350

! File names

  CHARACTER(LEN = 256) :: data_fname, log_fname, ana_fname
  CHARACTER(LEN = 256) :: polydisp_graft_fname, polydisp_bulk_fname
  CHARACTER(LEN = 256) :: box_fname

! Global Arrays involved in creating data file
  
  REAL,ALLOCATABLE,DIMENSION(:,:) :: rxyz, uxyz
  REAL,ALLOCATABLE,DIMENSION(:) :: charge
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: aidvals,ixyz
  INTEGER, ALLOCATABLE, DIMENSION(:) :: mwarr_graft,mwarr_bulk

! Character Arrays for creating the data file name

  CHARACTER (LEN = 12) :: f_char
  CHARACTER (LEN = 60 ):: datafile

  TYPE (RAN_SAVE) :: X
  INTEGER(4) :: S
  
END MODULE PARAMS
