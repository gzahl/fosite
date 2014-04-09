!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_gauss2d.f90                                                  #
!#                                                                           #
!# Copyright (C) 2006-2012                                                   #
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
! Program and data initialization for 2D Gaussian pressure pulse
!----------------------------------------------------------------------------!
PROGRAM Init
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE sources_generic
  USE timedisc_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameter
  REAL, PARAMETER :: TSIM   = 0.3         ! simulation time
  REAL, PARAMETER :: GAMMA  = 1.4         ! ratio of specific heats
  REAL, PARAMETER :: CSISO  = 0.0         ! if .ne. 0.0 -> isothermal simulation
                                          !   with CSISO as sound speed
  ! initial condition (dimensionless units)
  REAL, PARAMETER :: RHO0   = 1.0         ! ambient density
  REAL, PARAMETER :: P0     = 1.0         ! ambient pressure
  REAL, PARAMETER :: AMP    = 1.0         ! amplitude of the pulse
  REAL, PARAMETER :: PWIDTH = 0.06        ! half width of the Gaussian
  REAL, PARAMETER :: ETA    = 0.0         ! dynamic viscosity (0.0 disables)
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CARTESIAN   ! geometry of the mesh
!!$  INTEGER, PARAMETER :: MGEO = POLAR
!!$  INTEGER, PARAMETER :: MGEO = LOGPOLAR
!!$  INTEGER, PARAMETER :: MGEO = TANPOLAR
!!$  INTEGER, PARAMETER :: MGEO = SINHPOLAR
  INTEGER, PARAMETER :: XRES = 100         ! resolution
  INTEGER, PARAMETER :: YRES = 100
  REAL, PARAMETER    :: RMIN = 1.0E-2      ! inner radius for polar grids
  REAL, PARAMETER    :: RMAX = 0.3         ! outer radius
  REAL, PARAMETER    :: GPAR = 0.2         ! geometry scaling parameter
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 10        ! number of output time steps
  CHARACTER(LEN=256), PARAMETER :: ODIR &! output directory
                                 = "./"
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'gauss2d' 
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

CALL InitFosite(Sim)

CALL InitProgram(Sim%Mesh, Sim%Physics, Sim%Fluxes, Sim%Timedisc, &
                 Sim%Datafile, Sim%Logfile)

! set initial condition
CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)

CALL RunFosite(Sim)

CALL CloseFosite(Sim)

CONTAINS

  SUBROUTINE InitProgram(Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Timedisc_TYP):: Timedisc
    TYPE(FileIO_TYP)  :: Datafile
    TYPE(FileIO_TYP)  :: Logfile
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(4)
    REAL              :: x1,x2,y1,y2
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!
    ! mesh settings
    SELECT CASE(MGEO)
    CASE(CARTESIAN)
       x1 = -0.5
       x2 = 0.5
       y1 = -0.5
       y2 = 0.5
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = NO_GRADIENTS
       bc(NORTH) = NO_GRADIENTS
    CASE(POLAR)
       x1 = RMIN
       x2 = 0.5*SQRT(2.0)
       y1 = 0.0
       y2 = 2*PI
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(LOGPOLAR)
       x1 = LOG(RMIN/GPAR)
       x2 = LOG(0.5*SQRT(2.0)/GPAR)
       y1 = 0.0
       y2 = 2*PI
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(TANPOLAR)
       x1 = ATAN(RMIN/GPAR)
       x2 = ATAN(0.5*SQRT(2.0)/GPAR)
       y1 = 0.0
       y2 = 2*PI
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(SINHPOLAR)
       x1 = RMIN/GPAR                ! temporary
       x1 = LOG(x1+SQRT(1.0+x1*x1))  ! = ASINH(RMIN/sc)
       x2 = 0.5*SQRT(2.0)/GPAR       ! temporary
       x2 = LOG(x2+SQRT(1.0+x2*x2))  ! = ASINH(RMAX/sc)
       y1 = 0.0
       y2 = 2*PI
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE DEFAULT
       CALL Error(Physics,"InitProgram", &
            " geometry should be one of cartesian,polar,logpolar,tanpolar,sinhpolar or bipolar")
    END SELECT

    CALL InitMesh(Mesh,&
         meshtype = MIDPOINT, &
         geometry = MGEO, &
             inum = XRES, &       ! resolution in x and            !
             jnum = YRES, &       !   y direction                  !             
             xmin = x1, &
             xmax = x2, &
             ymin = y1, &
             ymax = y2, &
           gparam = GPAR)

    ! physics settings
    IF (CSISO.GT.TINY(CSISO)) THEN
       CALL InitPhysics(Physics,Mesh, &
            problem = EULER2D_ISOTHERM, &
            cs      = CSISO)                       ! isothermal sound speed  !
    ELSE
       CALL InitPhysics(Physics,Mesh, &
            problem   = EULER2D, &
            gamma     = GAMMA, &            ! ratio of specific heats        !
            dpmax     = 1.0)                ! for advanced time step control !
    END IF

    ! flux calculation and reconstruction method
    CALL InitFluxes(Fluxes,Mesh,Physics, &
         order     = LINEAR, &
         variables = CONSERVATIVE, &        ! vars. to use for reconstruction!
         limiter   = MONOCENT, &    ! one of: minmod, monocent,...   !
         theta     = 1.2)           ! optional parameter for limiter !

    ! boundary conditions (depends on the geometry, see above)
    CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = bc(WEST), &
         eastern  = bc(EAST), &
         southern = bc(SOUTH), &
         northern = bc(NORTH))

    ! viscosity source term
    IF (ETA.GT.TINY(ETA)) THEN
       CALL InitSources(Physics%sources,Mesh,Fluxes,Physics,Timedisc%boundary, &
            stype    = VISCOSITY, &
            vismodel = MOLECULAR, &
            dynconst = ETA)
    END IF

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = TSIM, &
         tol_rel  = 0.01, &
         dtlimit  = 1.0E-5, &
         maxiter  = 1000000)

    ! initialize log input/output
!!$    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,&
!!$         fileformat = BINARY, &
!!$         filename   = TRIM(ODIR) // TRIM(OFNAME) // 'log', &
!!$         filecycles = 1)

    ! initialize data input/output
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc, &
!!$         fileformat = VTK, &
         fileformat = GNUPLOT, filecycles = 0, &
         filename   = TRIM(ODIR) // TRIM(OFNAME), &
         count      = ONUM)
  END SUBROUTINE InitProgram


  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)  :: Physics
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Timedisc_TYP) :: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Sources_TYP), POINTER :: sp
    INTEGER           :: i,j
    REAL              :: x0, y0
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: dv
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,Physics
    INTENT(INOUT)      :: Timedisc
    !------------------------------------------------------------------------!
    ! center of the pressure pulse (give in cartesian coordinates)
    x0 = 0.0
    y0 = 0.0

    ! velocities
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.

    SELECT CASE(GetType(Physics))
    CASE(EULER2D)
       ! non-isothermal setup with constant density and pressure pulse
       Timedisc%pvar(:,:,Physics%DENSITY)   = RHO0
       FORALL (i=Mesh%IGMIN:Mesh%IGMAX,j=Mesh%JGMIN:Mesh%JGMAX)
          Timedisc%pvar(i,j,Physics%PRESSURE) = P0 + AMP*EXP(-LOG(2.0) * &
               ((Mesh%bccart(i,j,1)-x0)**2+(Mesh%bccart(i,j,2)-y0)**2)/PWIDTH**2)
       END FORALL
    CASE(EULER2D_ISOTHERM)
       ! in isothermal configurations the pressure is proportional to
       ! the density; thus the pulse is applied to the density
       FORALL (i=Mesh%IGMIN:Mesh%IGMAX,j=Mesh%JGMIN:Mesh%JGMAX)
          Timedisc%pvar(i,j,Physics%DENSITY) = RHO0 + AMP*EXP(-LOG(2.0) * &
               ((Mesh%bccart(i,j,1)-x0)**2+(Mesh%bccart(i,j,2)-y0)**2)/PWIDTH**2)
       END FORALL
    END SELECT

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)

    CALL Info(Mesh, " DATA-----> initial condition: " // &
         "2D gaussian pressure pulse")

  END SUBROUTINE InitData

END PROGRAM Init
