! To generate LAMMPS input file for polyelectrolyte simulations
! Use in conjunction with lmp_params.f90 and ran_numbers.f90
! Polydispersity is added: Ver: July-15-2019
PROGRAM LAMMPSINP

  USE PARAMS

  IMPLICIT NONE
  
  LOGICAL :: input_coor = .false.
  REAL :: bondl, bondlsq
  INTEGER :: ierror,narg

  bondlsq = 0

  CALL SYSTEM_CLOCK(S)

  narg = IARGC()
  
  CALL READ_INP_FILE()
  CALL PDI_READ_AND_ANALYZE()
  CALL COMPUTE_COUNTERIONS()
  CALL COMPUTE_POLY_ION_BOX_SALT_DETAILS()
  CALL WRITE_INIT_DETAILS()
  CALL ALLOCATE_ARRAYS()
  CALL GENERATE_INPCOR()
  CALL CREATE_ATOMTYPE_CHARGE()
  CALL CREATE_WRITE_DATAFILE()
  CALL DEALLOCATE_ARRAYS()

END PROGRAM LAMMPSINP

!--------------------------------------------------------------------

SUBROUTINE READ_INP_FILE()

  USE PARAMS
  IMPLICIT NONE
  
  INTEGER :: nargs,ierr,logflag,AllocateStatus,i,j, pdiflag
  CHARACTER(256) :: dumchar

  CALL DEFAULT_VALUES()

  nargs = IARGC()
  IF(nargs .NE. 1) STOP "Input incorrect"

  logflag = 0; pdiflag = 0

  CALL GETARG(nargs,inp_fname)

  OPEN(unit = inpread,file=trim(inp_fname),action="read",status="old"&
       &,iostat=ierr)
  
  IF(ierr /= 0) THEN

     PRINT *, trim(inp_fname), "not found"
     STOP

  END IF

  DO

     READ(inpread,*,iostat=ierr) dumchar

     IF(ierr .LT. 0) EXIT

     IF(dumchar == 'out_datafile') THEN
        
        READ(inpread,*,iostat=ierr) data_fname

     ELSEIF(dumchar == 'free_pdi_file') THEN
        
        READ(inpread,*,iostat=ierr) free_pdi_fname

     ELSEIF(dumchar == 'brush_pdi_file') THEN

        READ(inpread,*,iostat=ierr) brush_pdi_fname

     ELSEIF(dumchar == 'nfree_chains') THEN

        READ(inpread,*,iostat=ierr) nch_free

     ELSEIF(dumchar == 'free_avg_mw') THEN

        READ(inpread,*,iostat=ierr) avg_mon_free

     ELSEIF(dumchar == 'nbrush_chains') THEN
        READ(inpread,*,iostat=ierr) nch_brush

     ELSEIF(dumchar == 'brush_avg_mw') THEN

        READ(inpread,*,iostat=ierr) avg_mon_brush

     ELSEIF(dumchar == 'tail_mons') THEN

        READ(inpread,*,iostat=ierr) mon_tail_brush

     ELSEIF(dumchar == 'n_salt') THEN

        READ(inpread,*,iostat=ierr) n_salt

     ELSEIF(dumchar == 'charge_frac_chain') THEN

        READ(inpread,*,iostat=ierr) charge_frac

     ELSEIF(dumchar == 'architecture') THEN

        READ(inpread,*,iostat=ierr) arch
        logflag  = 1

     ELSEIF(dumchar == 'box_dimensions') THEN

        READ(inpread,*,iostat=ierr) boxl_x, boxl_y, boxl_z

     ELSEIF(dumchar == 'log_file') THEN

        READ(inpread,*,iostat=ierr) log_fname
        logflag  = 1

     ELSE
        
        PRINT *, "unknown keyword: ", trim(dumchar)
        STOP

     END IF

  END DO

  IF(logflag == 0) log_fname = "log."//trim(adjustl(data_fname))
  OPEN(unit = logout,file=trim(log_fname),action="write",status="repla&
       &ce",iostat=ierr)

  PRINT *, "Analysis input file read finished .."

END SUBROUTINE READ_INP_FILE

!--------------------------------------------------------------------

SUBROUTINE DEFAULT_VALUES()

  USE PARAMS
  IMPLICIT NONE

  nch_free = 0; avg_mon_free = 0; nch_brush = 0; avg_mon_brush = 0
  mon_tail_brush = 0; n_salt = 0; charge_frac = 0.0
  boxl_x = 53.0; boxl_y = 53.0; boxl_z = 120.0

END SUBROUTINE DEFAULT_VALUES

!---------------------------------------------------------------------

SUBROUTINE PDI_READ_AND_ANALYZE()

  USE PARAMS
  IMPLICIT NONE

  INTEGER :: AllocateStatus, ierr
  INTEGER :: i,j, dum_freech, dum_brushch, misq !misq is for PDI
  !calculation; not important here

! read free chains
  OPEN(unit = freeread,file=trim(free_pdi_fname),action="read",status&
       &="old",iostat=ierr)
  
  IF(ierr /= 0) THEN

     PRINT *, trim(free_pdi_fname), "not found"
     STOP

  END IF

  READ(freeread,*) dum_freech, mw_tot_free, misq, pdi_free
  IF(dum_freech .NE. nch_free) THEN
     PRINT *, dum_freech, nch_free
     STOP "Mismatch of number of free chains"
  END IF

  ALLOCATE(free_mon_ptr(1:nch_free,2), Stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "Could not allocate free_mon_ptr"

  DO i = 1, nch_free
     
     READ(freeread,*) free_mon_ptr(i,1), free_mon_ptr(i,2)

  END DO

  CLOSE(freeread)


! read brush chains
  OPEN(unit = brushread,file=trim(brush_pdi_fname),action="read"&
       &,status="old",iostat=ierr)
  
  IF(ierr /= 0) THEN

     PRINT *, trim(brush_pdi_fname), "not found"
     STOP

  END IF

  READ(brushread,*) dum_brushch, mw_tot_brush, misq, pdi_brush
  IF(dum_brushch .NE. nch_brush) THEN
     PRINT *, dum_brushch, nch_brush
     STOP "Mismatch of number of brush chains"
  END IF

  ALLOCATE(brush_mon_ptr(1:nch_brush,2), Stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "Could not allocate brush_mon_ptr"

  DO i = 1, nch_brush
     
     READ(brushread,*) brush_mon_ptr(i,1), brush_mon_ptr(i,2)

  END DO
  CLOSE(brushread)

END SUBROUTINE PDI_READ_AND_ANALYZE

!----------------------------------------------------------------------

SUBROUTINE CREATE_WRITE_DATAFILE()

  USE PARAMS
  
  IMPLICIT NONE
  
  INTEGER :: i,j, k, ierror
  INTEGER :: ntotbonds, ntotangls,ntotdihds
  REAL ::  massval
  
  i = 1
  
  PRINT *, "Writing to LAMMPS Datafile .. "

  OPEN (unit=outdata, file = data_fname, status="replace",action="write"&
       &,iostat = ierror)
  
  IF(ierror /= 0) STOP "Failed to open datafile"
     
  WRITE (outdata,*) "Data for CG-PE simulations "
  WRITE (outdata,*) 
  WRITE (outdata,'(5X,I0,2X,A)') totpart, "atoms"

  ntotbonds = 0; ntotangls = 0; ntotdihds = 0

  CALL COMPUTE_TOTAL_TOPO_DETAILS(ntotbonds,ntotangls,ntotdihds)

  WRITE (outdata,'(5X,I0,2X,A)') ntotbonds, "bonds"
  WRITE (outdata,'(5X,I0,2X,A)') ntotangls, "angles"
  WRITE (outdata,'(5X,I0,2X,A)') ntotdihds, "dihedrals"

  WRITE (outdata,'(5X,I0,2X,A)') 0, "impropers"


  WRITE (outdata,'(5X,I0,2X,A)') numatomtypes, "atom types"
  WRITE (outdata,'(5X,I0,2X,A)') numbondtypes, "bond types"
  WRITE (outdata,'(5X,I0,2X,A)') numangltypes, "angle types"
  WRITE (outdata,'(5X,I0,2X,A)') numdihdtypes, "dihedral types"
  WRITE (outdata,'(5X,I0,2X,A)') 0, "improper types"

  WRITE (outdata,*)
  WRITE (outdata,'(5X,I0,2X,F14.6,2X,A)') 0, boxl_x, "xlo xhi"
  WRITE (outdata,'(5X,I0,2X,F14.6,2X,A)') 0, boxl_y, "ylo yhi"
  WRITE (outdata,'(5X,I0,2X,F14.6,2X,A)') 0, boxl_z, "zlo zhi"
  WRITE (outdata,*)
  WRITE (outdata,*) "Masses"
  WRITE (outdata,*)
  
  ! Writing Masses

  DO i = 1,numatomtypes

     massval = 1.000
     WRITE(outdata,'(I0,1X,F14.8)') i, massval

  END DO

  CALL SET_IMGFLAGS()
  CALL WRITE_POSITIONS()
  IF(numbondtypes /= 0) CALL CREATE_WRITE_BOND_TOPO()
  IF(numangltypes /= 0) CALL CREATE_WRITE_ANGLE_TOPO()
  IF(numdihdtypes /= 0) CALL CREATE_WRITE_DIHEDRAL_TOPO()

  CLOSE(unit = outdata)  

END SUBROUTINE CREATE_WRITE_DATAFILE

!----------------------------------------------------------------------

SUBROUTINE COMPUTE_TOTAL_TOPO_DETAILS(ntotbonds,ntotangls,ntotdihds)

  USE PARAMS
  IMPLICIT NONE

  INTEGER, INTENT(OUT):: ntotbonds,ntotangls,ntotdihds
  INTEGER :: i

  ntotbonds = 0; ntotangls = 0; ntotdihds = 0  

  IF(numbondtypes == 0) THEN

     ntotbonds = 0

  ELSE

     DO i = 1,nch_brush 
        ntotbonds = ntotbonds + brush_mon_ptr(i,2) -1
     END DO

     DO i = 1,nch_free
        ntotbonds = ntotbonds + free_mon_ptr(i,2) -1
     END DO

  END IF

  IF(numangltypes == 0) THEN

     ntotangls = 0

  ELSE

     DO i = 1,nch_brush 
        ntotangls = ntotangls + brush_mon_ptr(i,2) -2
     END DO

     DO i = 1,nch_free
        ntotangls = ntotangls + free_mon_ptr(i,2) -2
     END DO

  END IF

  IF(numdihdtypes == 0) THEN

     ntotdihds = 0

  ELSE

     DO i = 1,nch_brush 
        ntotdihds = ntotdihds + brush_mon_ptr(i,2) - 3
     END DO

     DO i = 1,nch_free
        ntotdihds = ntotdihds + free_mon_ptr(i,2) - 3
     END DO

  END IF


END SUBROUTINE COMPUTE_TOTAL_TOPO_DETAILS

!----------------------------------------------------------------------

SUBROUTINE WRITE_POSITIONS()

  USE PARAMS
  IMPLICIT NONE

  INTEGER i

  ! Writing atomic corrdinates
  
  WRITE (outdata,*) 
  WRITE (outdata,*) "Atoms"
  WRITE (outdata,*)

  DO i = 1,totpart
     
     WRITE(outdata,'(3(I0,1X),4(F14.6,1X),3(I0,1X))') aidvals(i,1),&
          & aidvals(i,2), aidvals(i,3), charge(i), rxyz(i,1), rxyz(i&
          &,2),rxyz(i,3), ixyz(i,1), ixyz(i,2), ixyz(i,3)
        
  END DO


END SUBROUTINE WRITE_POSITIONS

!----------------------------------------------------------------------

SUBROUTINE CREATE_WRITE_BOND_TOPO()

  USE PARAMS
  IMPLICIT NONE
  
  INTEGER :: i,j,k
  INTEGER :: bondid, sum_tot_brush, sum_tot_free

  bondid = 0
  WRITE (outdata,*)
  WRITE (outdata,*) "Bonds"
  WRITE (outdata,*)
  
  sum_tot_free = 0; sum_tot_brush = 0

  !Brush chains
  DO i = 1,nch_brush
     
     DO j = 1,brush_mon_ptr(i,2)-1
        
        k = sum_tot_brush + j
        bondid = bondid + 1
        
        WRITE(outdata,'(4(I0,2X))') bondid, bondtype, aidvals(k,1)&
             &,aidvals(k+1,1)
        
     END DO

     sum_tot_brush = sum_tot_brush + brush_mon_ptr(i,2)
     
  END DO


  ! Free chains
  DO i = 1,nch_free
     
     DO j = 1,free_mon_ptr(i,2) - 1
        
        bondid = bondid + 1
        k = sum_tot_free + j + mw_tot_brush
        
        WRITE(outdata,'(4(I0,2X))') bondid, bondtype, aidvals(k,1)&
             &,aidvals(k+1,1)
        
     END DO

     sum_tot_free = sum_tot_free + free_mon_ptr(i,2)
     
  END DO


END SUBROUTINE CREATE_WRITE_BOND_TOPO

!----------------------------------------------------------------------

SUBROUTINE CREATE_WRITE_ANGLE_TOPO()

  USE PARAMS
  IMPLICIT NONE

  INTEGER :: i,j,k
  INTEGER :: anglid, sum_tot_free, sum_tot_brush

  anglid = 0
  WRITE (outdata,*)
  WRITE (outdata,*) "Angles"
  WRITE (outdata,*)

  sum_tot_free = 0; sum_tot_brush = 0

  !Brush chains
  DO i = 1,nch_brush
     
     DO j = 1,brush_mon_ptr(i,2)-2
     
        anglid = anglid + 1
        k = sum_tot_brush + j           
        WRITE(outdata,'(5(I0,2X))') anglid, angltype, aidvals(k,1)&
             &,aidvals(k+1,1),aidvals(k+2,1)
           
     END DO

     sum_tot_brush = sum_tot_brush + brush_mon_ptr(i,2)
     
  END DO
  
  !Free chains
  DO i = 1,nch_free
     
     DO j = 1,free_mon_ptr(i,2)-2
        
        anglid = anglid + 1
        k = sum_tot_free + j + mw_tot_brush
        WRITE(outdata,'(5(I0,2X))') anglid, angltype, aidvals(k,1)&
             &,aidvals(k+1,1),aidvals(k+2,1)
        
     END DO

     sum_tot_free = sum_tot_free + free_mon_ptr(i,2)
     
  END DO
  
END SUBROUTINE CREATE_WRITE_ANGLE_TOPO

!----------------------------------------------------------------------

SUBROUTINE CREATE_WRITE_DIHEDRAL_TOPO()

  USE PARAMS
  IMPLICIT NONE

  INTEGER :: i,j,k
  INTEGER :: dihdid, sum_tot_brush, sum_tot_free
     
  dihdid = 0
  WRITE (outdata,*)
  WRITE (outdata,*) "Dihedrals"
  WRITE (outdata,*)
  

  sum_tot_free = 0; sum_tot_brush = 0

  !Brush chains
  DO i = 1,nch_brush
     
     DO j = 1,brush_mon_ptr(i,2)-3
     
        dihdid = dihdid + 1
        k = sum_tot_brush + j           
        WRITE(outdata,'(6(I0,2X))') dihdid, dihdtype, aidvals(k,1)&
             &,aidvals(k+1,1), aidvals(k+2,1), aidvals(k+3,1)
        
     END DO

     sum_tot_brush = sum_tot_brush + brush_mon_ptr(i,2)
     
  END DO
  
  !Free chains
  DO i = 1,nch_free
     
     DO j = 1,free_mon_ptr(i,2)-3
        
        dihdid = dihdid + 1
        k = sum_tot_free + j + mw_tot_brush
        WRITE(outdata,'(6(I0,2X))') dihdid, dihdtype, aidvals(k,1)&
             &,aidvals(k+1,1), aidvals(k+2,1), aidvals(k+3,1)
        
     END DO

     sum_tot_free = sum_tot_free + free_mon_ptr(i,2)
     
  END DO


END SUBROUTINE CREATE_WRITE_DIHEDRAL_TOPO

!--------------------------------------------------------------------  

SUBROUTINE SET_IMGFLAGS()

  USE PARAMS
  IMPLICIT NONE

  INTEGER :: i,j,k,sum_mon_brush, sum_mon_free
  REAL :: rx, ry, rz

!!$Fetching image information - Brush

  ixyz = 0
  sum_mon_brush = 0

  DO i = 1,nch_brush
     
     k = sum_mon_brush + 1
     ixyz(k,1) = 0
     ixyz(k,2) = 0
     ixyz(k,3) = 0
     
     DO j = 1,brush_mon_ptr(i,2)-1
        
        k = sum_mon_brush + j

        rx = rxyz(k,1) - rxyz(k+1,1)
        ry = rxyz(k,2) - rxyz(k+1,2)
        rz = rxyz(k,3) - rxyz(k+1,3)
        
        CALL IMGFLAGS(rx,ixyz(k,1),boxl_x,ixyz(k+1,1))
        CALL IMGFLAGS(ry,ixyz(k,2),boxl_y,ixyz(k+1,2))
        CALL IMGFLAGS(rz,ixyz(k,3),boxl_z,ixyz(k+1,3))
        
     END DO

     sum_mon_brush = sum_mon_brush + brush_mon_ptr(i,2)

  END DO

!!$Fetching image information - Free

  sum_mon_free = 0

  DO i = 1,nch_free
     
     k = mw_tot_brush + sum_mon_free + 1

     ixyz(k,1) = 0
     ixyz(k,2) = 0
     ixyz(k,3) = 0
     
     DO j = 1,free_mon_ptr(i,2)-1
        
        k = mw_tot_brush + sum_mon_free + j

        rx = rxyz(k,1) - rxyz(k+1,1)
        ry = rxyz(k,2) - rxyz(k+1,2)
        rz = rxyz(k,3) - rxyz(k+1,3)
        
        CALL IMGFLAGS(rx,ixyz(k,1),boxl_x,ixyz(k+1,1))
        CALL IMGFLAGS(ry,ixyz(k,2),boxl_y,ixyz(k+1,2))
        CALL IMGFLAGS(rz,ixyz(k,3),boxl_z,ixyz(k+1,3))

     END DO
     
     sum_mon_free = sum_mon_free + free_mon_ptr(i,2)

  END DO


END SUBROUTINE SET_IMGFLAGS

!---------------------------------------------------------------------  

SUBROUTINE IMGFLAGS(dist, img, boxl, imgout)
    
  USE PARAMS
  IMPLICIT NONE
  
  REAL, INTENT(IN) :: dist,boxl
  INTEGER, INTENT(IN) :: img
  INTEGER, INTENT(OUT) :: imgout
  INTEGER :: nx
  
  IF(dist > boxl/2) THEN
     
     nx = img + 1
     
  ELSEIF(dist < -boxl/2) THEN
     
     nx = img - 1
     
  ELSE
     
     nx = img
     
  END IF
  
  imgout = nx
  
END SUBROUTINE IMGFLAGS

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_POLY_ION_BOX_SALT_DETAILS()

  USE PARAMS

  IMPLICIT NONE
  INTEGER :: ierr

  nchains = nch_free + nch_brush
  totpart = mw_tot_brush + mw_tot_free + ncntr_brush + ncntr_free +&
       & 2.0*n_salt
  volbox  = boxl_x*boxl_y*boxl_z
  density = REAL(totpart)/REAL(volbox)
  
  npolyatoms = mw_tot_brush + mw_tot_free
  brush_dist = SQRT(REAL(boxl_x*boxl_y))/SQRT(REAL(2.0*nch_brush))
  ! Extra 2.0 factor above so that diagonal is approximately equal to
  !  sqrt(area/num_chains)

END SUBROUTINE COMPUTE_POLY_ION_BOX_SALT_DETAILS
!--------------------------------------------------------------------

SUBROUTINE WRITE_INIT_DETAILS() 
  
  USE PARAMS
  IMPLICIT NONE

  WRITE(logout,*) "Simulation Details ...."
  PRINT *, "Writing Simulation Details ...."
  WRITE(logout,*) "Total particles: ", totpart
  WRITE(logout,*) "Number of atomtypes: ", numatomtypes
  WRITE(logout,*) "# of Polyanion (free) chains: ", nch_free
  WRITE(logout,*) "Total # of polyanion monomers: ", mw_tot_free
  WRITE(logout,*) "PDI of free chains: ", pdi_free
  WRITE(logout,*) "# of Polycation (brush) chains: ", nch_brush
  WRITE(logout,*) "Total # of polycation monomers: ", mw_tot_brush
  WRITE(logout,*) "PDI of brush chains: ", pdi_brush  
  WRITE(logout,*) "# of Salt: ", 2.0*n_salt
  WRITE(logout,*) "# of Brush counterions: ", ncntr_brush
  WRITE(logout,*) "# of Free counterions: ", ncntr_free
  WRITE(logout,*) "# of Polyelectrolytes: ", npolyatoms
  WRITE(logout,*) "Total Number of Particles: ", totpart
  WRITE(logout,*) "LX/LY/LZ: ", boxl_x, boxl_y, boxl_z
  WRITE(logout,*) "Density: ", density
  WRITE(logout,*) "Minimum distance between brush base: ", brush_dist

END SUBROUTINE WRITE_INIT_DETAILS

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_COUNTERIONS()

  USE PARAMS
  IMPLICIT NONE

  INTEGER :: i,neut_brush_mons,neut_free_mons
  INTEGER :: inv_charge

  inv_charge = INT(1.0/charge_frac)

  ncntr_brush = 0; ncntr_free = 0

  !IMPORTANT: THE mod value is to adjust for odd number of
  !monomers. If there are odd number of monomers (after subtracting
  !the tail monomers for the case of brush), then the extra monomer(s)
  !will be accounted as the neutral part.
  DO i = 1,nch_brush
     neut_brush_mons = mon_tail_brush + INT(charge_frac&
          &*(brush_mon_ptr(i,2)-mon_tail_brush)) +&
          & mod(brush_mon_ptr(i,2)-mon_tail_brush,inv_charge)
     ncntr_brush = ncntr_brush + brush_mon_ptr(i,2) - neut_brush_mons
  END DO

  DO i = 1,nch_free
     neut_free_mons = INT(charge_frac*free_mon_ptr(i,2)) +&
          & mod(free_mon_ptr(i,2),inv_charge)
     ncntr_free = ncntr_free + free_mon_ptr(i,2) - neut_free_mons
  END DO

END SUBROUTINE COMPUTE_COUNTERIONS

!--------------------------------------------------------------------

SUBROUTINE GENERATE_INPCOR()
  
  USE PARAMS

  IMPLICIT NONE
  
  INTEGER :: i,j,k,u,v,ierror,brushcnt, kold
  INTEGER :: sum_brush_mons, sum_free_mons
  REAL, PARAMETER :: r0init  = 0.97
  REAL, PARAMETER :: r0sq3   = r0init/sqrt(3.0)
  REAL, PARAMETER :: rmaxsq  = r0init*r0init
  REAL, PARAMETER :: math_pi = 3.14159265359
  REAL :: theta, phi
  REAL :: rx, ry, rz,rval
  LOGICAL :: in_box
  CALL RAN_INIT(S,X)
  
  WRITE(logout,*) "Random Initial Configuration : NRRW"

  i = 1; k = 0; sum_brush_mons = 0

! Create Brush

  WRITE(logout,*) "Generating Brush Configurations .. "


  DO WHILE(i .LE. nch_brush)
     
     k = sum_brush_mons + 1

     rxyz(k,1) = RAN1(X)*boxl_x
     rxyz(k,2) = RAN1(X)*boxl_y
     rxyz(k,3) = 0.50000

     !Non overlapping base

     brushcnt = 1; kold = 1

     DO WHILE(brushcnt .LT. i)
        
        rx = rxyz(k,1) - rxyz(kold,1)
        ry = rxyz(k,2) - rxyz(kold,2)
        rz = rxyz(k,3) - rxyz(kold,3)

        rx = rx - boxl_x*ANINT(rx/boxl_x)
        ry = ry - boxl_y*ANINT(ry/boxl_y)
        rz = rz - boxl_z*ANINT(rz/boxl_z)

        rval = sqrt(rx**2 + ry**2 + rz**2)

        IF(rval .GE. brush_dist) THEN
           kold = kold + brush_mon_ptr(brushcnt,2)
           brushcnt = brushcnt + 1
        ELSE
           rxyz(k,1) = RAN1(X)*boxl_x
           rxyz(k,2) = RAN1(X)*boxl_y
           brushcnt = 1
           kold = 1
        END IF

     END DO

     j = 2

     DO WHILE (j .LE. brush_mon_ptr(i,2))

        k = sum_brush_mons + j
        in_box = .false.

        DO WHILE(in_box == .false.)

           theta     = math_pi*RAN1(X)
           phi       = 2*math_pi*RAN1(X)
           rxyz(k,1) = rxyz(k-1,1) + r0init*sin(theta)*cos(phi)
           rxyz(k,2) = rxyz(k-1,2) + r0init*sin(theta)*sin(phi)
           rxyz(k,3) = rxyz(k-1,3) + r0init*cos(theta)
           IF(rxyz(k,3) .GE. 1.0 .AND. rxyz(k,3) .LE. boxl_z-1.0)&
                & in_box=.true.

        END DO

        j = j + 1
                
     END DO
     
     sum_brush_mons = sum_brush_mons + brush_mon_ptr(i,2)
     i = i + 1

        
  END DO

! Create Free Chains

  WRITE(logout,*) "Generating Free Chain Configurations .. "

  i = 1; sum_free_mons = 0

  DO WHILE(i .LE. nch_free)
     
     k = sum_free_mons + 1 + mw_tot_brush

     in_box = .false.
     DO WHILE(in_box == .false.)

        rxyz(k,1) = RAN1(X)*boxl_x
        rxyz(k,2) = RAN1(X)*boxl_y
        rxyz(k,3) = RAN1(X)*boxl_z
        IF(rxyz(k,3) .GT. 1.0 .AND. rxyz(k,3) .LT. boxl_z-1.0)&
             & in_box=.true.

     END DO

     j = 2

     DO WHILE (j .LE. free_mon_ptr(i,2))
        
        k = sum_free_mons + j + mw_tot_brush
        in_box = .false.

        DO WHILE(in_box == .false.)

           theta     = math_pi*RAN1(X)
           phi       = 2*math_pi*RAN1(X)
           rxyz(k,1) = rxyz(k-1,1) + r0init*sin(theta)*cos(phi)
           rxyz(k,2) = rxyz(k-1,2) + r0init*sin(theta)*sin(phi)
           rxyz(k,3) = rxyz(k-1,3) + r0init*cos(theta)
           IF(rxyz(k,3) .GE. 1.0 .AND. rxyz(k,3) .LE. boxl_z-1.0)&
                & in_box=.true.

        END DO

        j = j + 1
                
     END DO

     sum_free_mons = sum_free_mons + free_mon_ptr(i,2)     
     i = i + 1
        
  END DO

! Create Salt 

  WRITE(logout,*) "Generating Salt .. "

  k = 1 + mw_tot_brush + mw_tot_free; i = 1
  DO WHILE(i .LE. 2.0*N_salt)
     
     in_box = .false.
     
     DO WHILE(in_box == .false.)

        rxyz(k,1) = RAN1(X)*boxl_x
        rxyz(k,2) = RAN1(X)*boxl_y
        rxyz(k,3) = RAN1(X)*boxl_z
        IF(rxyz(k,3) .GE. 1.0 .AND. rxyz(k,3) .LE. boxl_z-1.0)&
             & in_box=.true.
        
     END DO
     
     i = i + 1; k = k + 1
        
  END DO


! Create brush counterions

  k = 1 + mw_tot_brush + mw_tot_free +2.0*N_salt
  i = 1
  DO WHILE(i .LE. ncntr_brush)
     
     in_box = .false.
     
     DO WHILE(in_box == .false.)

        rxyz(k,1) = RAN1(X)*boxl_x
        rxyz(k,2) = RAN1(X)*boxl_y
        rxyz(k,3) = RAN1(X)*boxl_z
        IF(rxyz(k,3) .GE. 1.0 .AND. rxyz(k,3) .LE. boxl_z-1.0)&
             & in_box=.true.
        
     END DO
     
     i = i + 1; k = k + 1
        
  END DO

! Create free counterions

  k = 1 + mw_tot_brush + mw_tot_free +2.0*N_salt + ncntr_brush
  i = 1
  DO WHILE(i .LE. ncntr_free)
     
     in_box = .false.
     
     DO WHILE(in_box == .false.)

        rxyz(k,1) = RAN1(X)*boxl_x
        rxyz(k,2) = RAN1(X)*boxl_y
        rxyz(k,3) = RAN1(X)*boxl_z
        IF(rxyz(k,3) .GE. 1.0 .AND. rxyz(k,3) .LE. boxl_z-1.0)&
             & in_box=.true.
        
     END DO
     
     i = i + 1; k = k + 1
        
  END DO


  ! PBC
  
  DO i = 1,totpart
     
     rxyz(i,1) = rxyz(i,1) - boxl_x*floor(rxyz(i,1)/boxl_x)
     rxyz(i,2) = rxyz(i,2) - boxl_y*floor(rxyz(i,2)/boxl_y)
     rxyz(i,3) = rxyz(i,3) - boxl_z*floor(rxyz(i,3)/boxl_z)
        
  END DO
  
END SUBROUTINE GENERATE_INPCOR

!--------------------------------------------------------------------

SUBROUTINE CREATE_ATOMTYPE_CHARGE()

  USE PARAMS
  IMPLICIT NONE

  INTEGER :: i,j,k
  INTEGER :: fin_neut_brush, fin_neut_free
  INTEGER :: sum_mon_brush, sum_mon_free
  REAL :: csum
  INTEGER :: inv_charge
  INTEGER :: charg_mon_polyan, charg_mon_polycat
  INTEGER :: n_charg_polybrush, n_charg_polyfree


  OPEN(unit=45,file="poly_charg.txt",action="write",status="replace")

  inv_charge = INT(1.0/charge_frac)
  csum = 0; charge = 0
  
  WRITE(logout,*) "Generating Atomtypes and Charges .. "
! Brush

  sum_mon_brush = 0

  WRITE(45,'(A)') trim(adjustl("Polycation brush: chID, totmons, neut,&
       & charged"))

  n_charg_polybrush = 0; n_charg_polyfree = 0

  DO i = 1, nch_brush

     fin_neut_brush = mon_tail_brush + INT(charge_frac&
          &*(brush_mon_ptr(i,2)-mon_tail_brush)) +&
          & mod(brush_mon_ptr(i,2)-mon_tail_brush,inv_charge)
     !The reminder is for a generic charge fraction case.

     WRITE(45,'(4(I0,1X))') i,brush_mon_ptr(i,2),fin_neut_brush&
          &,brush_mon_ptr(i,2)-fin_neut_brush

     DO j = 1,brush_mon_ptr(i,2)

        k = sum_mon_brush + j

        aidvals(k,1) = k
        aidvals(k,2) = i

        IF(j == 1) THEN
           aidvals(k,3) = 1
           charge(k)  = 0.0 
        ELSEIF(j .LE. mon_tail_brush) THEN
           aidvals(k,3) = 2
           charge(k) = 0.0
        ELSEIF(arch .LE. 2) THEN
           IF(j .LE. fin_neut_brush) THEN
              aidvals(k,3) = 3
              charge(k) = 0.0
           ELSE
              aidvals(k,3) = 4
              charge(k) = 1.0
              n_charg_polybrush = n_charg_polybrush + 1
           END IF
        ELSE
           IF(mod(j-mon_tail_brush,2) .NE. 0) THEN
              aidvals(k,3) = 3
              charge(k) = 0.0
           ELSE
              aidvals(k,3) = 4
              charge(k) = 1.0
              n_charg_polybrush = n_charg_polybrush + 1
           END IF
        END IF

     END DO

     sum_mon_brush = sum_mon_brush + brush_mon_ptr(i,2)

  END DO

  IF(n_charg_polybrush .NE. ncntr_brush) THEN
     PRINT *, "ERROR: Unequal brush and counterions .."
     PRINT *, ncntr_brush, n_charg_polybrush
     STOP
  END IF

! Free

  WRITE(45,'(A)') trim(adjustl("Polyanion free: chID, totmons, neut, c&
       &harged"))
  sum_mon_free = 0

  DO i = 1,nch_free

     fin_neut_free = INT(charge_frac*free_mon_ptr(i,2)) +&
          & mod(free_mon_ptr(i,2),inv_charge)

     WRITE(45,'(4(I0,1X))') i,free_mon_ptr(i,2),fin_neut_free&
          &,free_mon_ptr(i,2)-fin_neut_free

     DO j = 1,free_mon_ptr(i,2)

        k = sum_mon_free + j + mw_tot_brush

        aidvals(k,1) = k
        aidvals(k,2) = i + nch_brush

        IF(arch == 1 .OR. arch == 3) THEN
           IF(j .LE. fin_neut_free) THEN
              aidvals(k,3) = 5
              charge(k) = 0.0
           ELSE
              aidvals(k,3) = 6
              charge(k) = -1.0
              n_charg_polyfree = n_charg_polyfree + 1
           END IF
        ELSE
           IF(mod(j,2) /= 0) THEN
              aidvals(k,3) = 5
              charge(k) = 0.0
           ELSE
              aidvals(k,3) = 6
              charge(k) = -1.0
              n_charg_polyfree = n_charg_polyfree + 1
           END IF
        END IF

     END DO

     sum_mon_free = sum_mon_free + free_mon_ptr(i,2)

  END DO

  CLOSE(45)

  IF(n_charg_polyfree .NE. ncntr_free) THEN
     PRINT *, "ERROR: Unequal freepoly and counterions .."
     PRINT *, ncntr_free, n_charg_polyfree
     STOP
  END IF


! Salt

  DO i = 1,2*N_salt

     k = i + mw_tot_free + mw_tot_brush

     aidvals(k,1) = k
     aidvals(k,2) = 1 + nchains

     IF(mod(i,2) /=  0) THEN
        aidvals(k,3) = 7
        charge(k) = 1.0
     ELSE
        aidvals(k,3) = 8
        charge(k) = -1.0
     END IF
     
  END DO

! Counterions in polymers - will have same type as the salt

  DO i = 1,ncntr_brush

     k = i + mw_tot_free + mw_tot_brush + 2*N_salt
     aidvals(k,1) = k
     aidvals(k,2) = 1 + nchains
     aidvals(k,3) = 8
     charge(k) = -1.0

  END DO

  DO i = 1,ncntr_free

     k = i + mw_tot_free + mw_tot_brush + 2*N_salt + ncntr_brush
     aidvals(k,1) = k
     aidvals(k,2) = 1 + nchains
     aidvals(k,3) = 7
     charge(k) = 1.0

  END DO

! Charge Neutrality

  DO i = 1, totpart

     csum = csum + charge(i)
     
  END DO

  IF(abs(csum) .GT. 0.00000001) THEN

     OPEN(unit = 77,file="charge.txt",action="write",status="replace")

     DO i = 1,totpart
        WRITE(77,*), i,aidvals(i,3), charge(i)
     END DO

     CLOSE(77)
     WRITE(logout,*) "System not charge neutral", csum
     PRINT *, "ERROR: System not charge neutral", csum
     STOP
     
  ELSE

     OPEN(unit = 77,file="charge_init.txt",action="write",status="repl&
          &ace")

     DO i = 1,totpart
        WRITE(77,*), i,aidvals(i,3), charge(i)
     END DO
     
     CLOSE(77)
     WRITE(logout,*) "Good Charge Neutrality ", csum
     PRINT *, "Good Charge Neutrality ", csum

  END IF
     
END SUBROUTINE CREATE_ATOMTYPE_CHARGE

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_ARRAYS()

  USE PARAMS
  IMPLICIT NONE

  INTEGER :: AllocateStatus

  ALLOCATE(rxyz(1:totpart,3), Stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "Could not allocate rxyz"
  ALLOCATE(uxyz(1:totpart,3), Stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "Could not allocate uxyz"
  ALLOCATE(charge(1:totpart), Stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "Could not allocate charge"
  ALLOCATE(aidvals(1:totpart,3), Stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "Could not allocate aidvals"
  ALLOCATE(ixyz(1:totpart,3), Stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "Could not allocate ixyz"


END SUBROUTINE ALLOCATE_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE DEALLOCATE_ARRAYS()

  USE PARAMS
  IMPLICIT NONE

  DEALLOCATE(rxyz)
  DEALLOCATE(uxyz)
  DEALLOCATE(ixyz)
  DEALLOCATE(charge)
  DEALLOCATE(aidvals)

END SUBROUTINE DEALLOCATE_ARRAYS

!--------------------------------------------------------------------
