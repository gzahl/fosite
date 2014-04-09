!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_pringle.f90                                                  #
!#                                                                           #
!# Copyright (C) 2008 - 2010                                                 #
!# Bjoern Sperling  <sperling@astrophysik.uni-kiel.de>                       #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!#                                                                           #
!# This program is free software; you can redistribute it and/or modify      #
!# it under the terms of the GNU General Public License as published by      #
!# the Free Software Foundation; either version 2 of the License, or (at     #
!# your option) any later version.                                           #
!#                                                                           #
!# This program is distributed in the hope that it will be useful, but       #
!# WITHOUT ANY WARRANTY; without even the implied warranty of                #
!# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, GOOD TITLE or        #
!# NON INFRINGEMENT.  See the GNU General Public License for more            #
!# details.                                                                  #
!#                                                                           #
!# You should have received a copy of the GNU General Public License         #
!# along with this program; if not, write to the Free Software               #
!# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
!#                                                                           #
!#############################################################################

!----------------------------------------------------------------------------!
! Program and data initialization for a thin viscous ring rotating around a
! central pointmass
! References:
!  [1] Lynden-Bell, D. & Pringle, J. E. Evolution of Viscous Disks and Origin
!      of Nebular Variables, M.N.R.A.S vol. 168(3) (1974) pp. 603-637
!      http://adsabs.harvard.edu/abs/1974MNRAS.168..603L
!  [2] Pringle, J. E. Accretion discs in astrophysics
!      Annual review of astronomy and astrophysics. vol. 19 (1981) pp. 137-162  
!      DOI: 10.1146/annurev.aa.19.090181.001033
!  [3] Speith, R. and Kley, W. Stability of the viscously spreading ring
!      Astron. Astrophys. vol. 399 (2003), pp. 395-407
!      DOI: 10.1051/0004-6361:20021783
!----------------------------------------------------------------------------!

!**************************************!
!* IMPORTANT:                         *!
!* - compile with autodouble          *!
!**************************************!
 
MODULE Init
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE geometry_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE sources_generic
  USE timedisc_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! general constants
  REAL, PARAMETER :: GN      = 6.6742D-11  ! Newtons grav. constant [SI]
  ! simulation parameters
  REAL, PARAMETER :: TSIM     = 1.0        ! simulation time [TAU] see below
  ! the equations are solved in non-dimensional units, using the Keplerian
  ! velocity at the location of the initial ring at R=1 as the velocity scale
  REAL, PARAMETER :: CENTMASS = 1./GN      ! fixed for non-dim. equations
  ! these are the basic parameters; for stable solutions one requires
  ! MA >> 1 and RE > MA
  REAL, PARAMETER :: RE       = 1.0E+3     ! Reynolds number (at R=1)
  REAL, PARAMETER :: MA       = 1.0E+2     ! Mach number (at R=1)
  ! lower limit for nitial density
  REAL, PARAMETER :: RHOMIN   = 1.0E-40    ! minimal initial density
  ! viscosity prescription
!!$  INTEGER, PARAMETER :: VISTYPE = BETA     
  INTEGER, PARAMETER :: VISTYPE = PRINGLE
  REAL, PARAMETER :: TAU0     = 0.01       ! time for initial condition [TAU]
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = POLAR       ! geometry
!!$  INTEGER, PARAMETER :: MGEO = LOGPOLAR
!!$  INTEGER, PARAMETER :: MGEO = TANPOLAR
!!$  INTEGER, PARAMETER :: MGEO = SINHPOLAR
  INTEGER, PARAMETER :: XRES = 100         ! x-resolution
  INTEGER, PARAMETER :: YRES = 1           ! y-resolution
  REAL, PARAMETER    :: RMIN = 0.1         ! min radius of comp. domain
  REAL, PARAMETER    :: RMAX = 2.0         ! max radius of comp. domain
  REAL, PARAMETER    :: GPAR = 0.1         ! geometry scaling parameter
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 100         ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'pringle' 
  ! derived parameters
  REAL, PARAMETER    :: CSISO = 1./MA      ! isothermal speed of sound
  REAL               :: TAU                ! viscous time scale
  REAL               :: OMEGA0             ! angular velocity at R0
  !--------------------------------------------------------------------------!

  
  PUBLIC :: &
       ! methods
       InitProgram
  !--------------------------------------------------------------------------!


CONTAINS

  SUBROUTINE InitProgram(Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Timedisc_TYP):: Timedisc
    TYPE(FILEIO_TYP)  :: Datafile
    TYPE(FILEIO_TYP)  :: Logfile
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL              :: x1,x2,y1,y2
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!
    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER2D_ISOTHERM, &
         gamma   = 1.4, &                   ! ratio of specific heats        !
         cs      = CSISO, &                 ! isothermal speed of sound      !
         rhomin  = 1.0E-50, &
         dpmax   = 1.0)                     ! for advanced time step control !

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)                 ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         variables = PRIMITIVE, &           ! vars. to use for reconstruction!
         limiter   = MINMOD, &
         theta     = 1.2)                   ! optional parameter for limiter !

    ! mesh settings
    SELECT CASE(MGEO)
    CASE(POLAR)
       x1 = RMIN
       x2 = RMAX
       y1 = 0.0 
       y2 = 2*PI       
    CASE(LOGPOLAR)
       x1 = LOG(RMIN/GPAR)
       x2 = LOG(RMAX/GPAR)
       y1 = 0.0 
       y2 = 2*PI       
    CASE(TANPOLAR)
       x1 = ATAN(RMIN/GPAR)
       x2 = ATAN(RMAX/GPAR)
       y1 = 0.0 
       y2 = 2*PI       
    CASE(SINHPOLAR)
       x1 = RMIN/GPAR
       x1 = LOG(x1+SQRT(1.0+x1*x1))  ! = ASINH(RMIN/GPAR))
       x2 = RMAX/GPAR
       x2 = LOG(x2+SQRT(1.0+x2*x2))  ! = ASINH(RMAX/GPAR))
       y1 = 0.0 
       y2 = 2*PI       
    CASE DEFAULT
       CALL Error(Physics,"InitProgram","mesh geometry not supported for Pringle disk")
    END SELECT
    CALL InitMesh(Mesh,Fluxes, &
         geometry = MGEO, &
             inum = XRES, &
             jnum = YRES, &
             xmin = x1, &
             xmax = x2, &
             ymin = y1, &
             ymax = y2, &
           gparam = GPAR)

    ! boundary conditions
    CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = NO_GRADIENTS, &
!!$         western  = CUSTOM, &
         eastern  = NO_GRADIENTS, &
         southern = PERIODIC, &
         northern = PERIODIC)

    ! source term due to a point mass
    CALL InitSources(Physics%sources,Mesh,Fluxes,Physics,Timedisc%Boundary, &
         stype  = POINTMASS, &            ! grav. accel. of a point mass     !
         mass = CENTMASS)                 ! mass of the accreting object[kg] !
    Physics%sources%outbound = 0          ! disable accretion

    ! compute viscous time scale
    SELECT CASE(VISTYPE)
    CASE(BETA)
       TAU = 4./27. * RE
    CASE(PRINGLE)
       TAU = RE / 12.0
    CASE DEFAULT
       CALL Error(Physics,"InitProgram","viscosity type not supported")
    END SELECT

    ! source term due to viscosity
    CALL InitSources(Physics%sources,Mesh,Fluxes,Physics,Timedisc%Boundary, &
         stype    = VISCOSITY, &
         vismodel = VISTYPE, &
         cvis     = 0.5, &
         dynconst = 1./RE) 

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = TSIM * TAU, &
         dtlimit  = 1E-8 * TAU, &
         maxiter  = 100000000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc)

    ! initialize log input/output
!!$    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,&
!!$         fileformat = BINARY, &
!!$         filename   = TRIM(ODIR) // TRIM(OFNAME) // 'log', &
!!$         filecycles = 1)

    ! initialize data input/output
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc, &
!!$         fileformat = NETCDF, &
         fileformat = GNUPLOT, &
         filecycles = 0, &
         filename   = TRIM(ODIR) // TRIM(OFNAME), &
         count      = ONUM)
  END SUBROUTINE InitProgram


  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: r
    CHARACTER(LEN=64) :: value
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! distance to the origin
    r(:,:) = SQRT(Mesh%bccart(:,:,1)**2+Mesh%bccart(:,:,2)**2)
        
    ! rotational velocity
    Timedisc%pvar(:,:,Physics%YVELOCITY) = SQRT(Physics%constants%GN  &
         * CENTMASS / r(:,:))

    SELECT CASE(VISTYPE)
    CASE(BETA)
       r(:,:) = r(:,:)**0.75
       ! surface density
       Timedisc%pvar(:,:,Physics%DENSITY) = 3. / SQRT(TAU0*(4*PI*r(:,:))**3) &
            * EXP(-((1.0-r(:,:))**2)/TAU0) + RHOMIN
       ! radial velocity
       Timedisc%pvar(:,:,Physics%XVELOCITY) = 4.5/(RE*TAU0) * r(:,:)*(r(:,:)-1.0)  &
            * Timedisc%pvar(:,:,Physics%YVELOCITY)
    CASE(PRINGLE)
       ! surface density
       Timedisc%pvar(:,:,Physics%DENSITY) = 1. / SQRT(4*TAU0*(PI*SQRT(r(:,:)))**3) &
            * EXP(-((1.0-r(:,:))**2)/TAU0) + RHOMIN
       ! radial velocity
       Timedisc%pvar(:,:,Physics%XVELOCITY) = 6.0/(RE*TAU0) * (r(:,:)-1.0)
    END SELECT

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)

    ! custom boundary conditions if requested
    IF (GetType(Timedisc%boundary(WEST)).EQ.CUSTOM) THEN
       Timedisc%boundary(WEST)%cbtype(:,Physics%DENSITY) = CUSTOM_EXTRAPOL
       Timedisc%boundary(WEST)%cbtype(:,Physics%XVELOCITY) = CUSTOM_EXTRAPOL
       Timedisc%boundary(WEST)%cbtype(:,Physics%YVELOCITY) = CUSTOM_NOGRAD
    END IF

    CALL Info(Mesh, " DATA-----> initial condition: " // "pringle disk")
    WRITE(value,"(ES9.3)") TAU
    CALL Info(Mesh, "                               " // "viscous timescale:  " //TRIM(value))
    WRITE(value,"(ES9.3)") RE
    CALL Info(Mesh, "                               " // "Reynolds number:    " //TRIM(value))
    WRITE(value,"(ES9.3)") MA
    CALL Info(Mesh, "                               " // "Mach number:        " //TRIM(value))

  END SUBROUTINE InitData

END MODULE Init
