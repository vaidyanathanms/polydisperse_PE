! Input file to generate LAMMPS file for polyelectrolyte simulations
! Use in conjunction with lammps_inp.f90 and ran_numbers.f90
! Change the required parameters
! arch details: "arch=number:polycation-polyanion 
! arch = 1:block-block;arch=2:block-alt;arch=3:alt-block;arch=4:alt
! -alt.
MODULE PARAMS

  USE RAN_NUMBERS

  IMPLICIT NONE

! Parameter data for creating the data file

  INTEGER :: nch_free, avg_mon_free, nch_brush, avg_mon_brush
  INTEGER :: mon_tail_brush
  INTEGER :: n_salt
!  Default box size is set inside the code
!  INTEGER :: default_dim !If 1=>53*53*120.Or else use  !box.dat
  REAL    :: charge_frac = 0.5
  REAL    :: brush_dist
  INTEGER :: arch

! Box/Particle details

  REAL :: boxl_x, boxl_y, boxl_z
  REAL :: volbox, density
  INTEGER :: totpart
  INTEGER :: nchains, npolyatoms
  INTEGER :: mw_tot_free, mw_tot_brush
  INTEGER :: ncntr_brush
  INTEGER :: ncntr_free
  REAL    :: pdi_free, pdi_brush

! Files

  CHARACTER(LEN = 256) :: data_fname, log_fname, inp_fname
  CHARACTER(LEN = 256) :: free_pdi_fname, brush_pdi_fname,box_fname
  INTEGER, PARAMETER :: outdata = 100, logout = 110, inpread = 120
  INTEGER, PARAMETER :: freeread = 130, brushread = 140


! Flags for creating the data file

  INTEGER, PARAMETER :: numatomtypes = 8
  INTEGER, PARAMETER :: numbondtypes = 1
  INTEGER, PARAMETER :: numangltypes = 1
  INTEGER, PARAMETER :: numdihdtypes = 0
  INTEGER, PARAMETER :: bondtype = 1
  INTEGER, PARAMETER :: angltype = 1
  INTEGER, PARAMETER :: dihdtype = 0

! Global arrays involved in creating data file
  
  REAL,ALLOCATABLE,DIMENSION(:,:) :: rxyz, uxyz
  REAL,ALLOCATABLE,DIMENSION(:) :: charge
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: aidvals,ixyz
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: free_mon_ptr,brush_mon_ptr

! Random number generator

  TYPE (RAN_SAVE) :: X
  INTEGER(4) :: S
  
END MODULE PARAMS
