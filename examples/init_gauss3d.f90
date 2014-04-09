!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_gauss3d.f90                                                  #
!#                                                                           #
!# Copyright (C) 2006-2010                                                   #
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
! Program and data initialization for 3D Gaussian pressure or density pulse
! with and without rotation
! [1] Illenseer, T. F., Duschl, W. J.: Two-dimensional central-upwind schemes
!     for curvilinear grids and application to gas dynamics with angular momentum,
!     Comput. Phys. Comm. 180 (2009), 2283-2302
!     DOI: 10.1016/j.cpc.2009.07.016
!----------------------------------------------------------------------------!
MODULE Init
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE timedisc_generic
  USE sources_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! simulation parameter
  REAL, PARAMETER :: TSIM   = 0.6         ! simulation time
  REAL, PARAMETER :: GAMMA  = 1.4         ! ratio of specific heats
  ! initial condition (dimensionless units)
  REAL, PARAMETER :: RHO0   = 1.0         ! ambient density
  REAL, PARAMETER :: RHO1   = 0.0         ! peak density above RHO0
  REAL, PARAMETER :: RWIDTH = 0.06        ! half width of the Gaussian
  REAL, PARAMETER :: P0     = 1.0         ! ambient pressure
  REAL, PARAMETER :: P1     = 1.0         ! peak pressure above P0
  REAL, PARAMETER :: PWIDTH = 0.06        ! half width of the Gaussian
  REAL, PARAMETER :: OMEGA0 = 0.0         ! angular velocity 
  REAL, PARAMETER :: ETA    = 0.0         ! dynamic viscosity (0.0 disables)
  ! location of the pulse in cylindrical coordinates
  REAL, PARAMETER :: R0     = 0.0         ! radial position 
  REAL, PARAMETER :: Z0     = 0.0         ! vertical position
  ! mesh settings
!!$  INTEGER, PARAMETER :: MGEO = SPHERICAL   ! geometry
  INTEGER, PARAMETER :: MGEO = CYLINDRICAL
!!$  INTEGER, PARAMETER :: MGEO = OBLATE_SPHEROIDAL
!!$  INTEGER, PARAMETER :: MGEO = SINHSPHERICAL
  INTEGER, PARAMETER :: XRES  = 100       ! x-resolution
  INTEGER, PARAMETER :: YRES  = 100       ! y-resolution
  REAL, PARAMETER    :: RMAX  = 1.0       ! width of square that fits into
                                          !   computational domain
  REAL, PARAMETER    :: GPAR  = 0.8       ! geometry scaling parameter
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 10         ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &         ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &         ! output data file name
                     :: OFNAME = 'gauss3d' 
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
    INTEGER           :: bc(4)
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!
    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER3D_ROTSYM, &
         gamma   = GAMMA, &                ! ratio of specific heats        !
         dpmax   = 1.0)                    ! for advanced time step control !

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)                ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
!!$         variables = PRIMITIVE, &       ! vars. to use for reconstruction!
         variables = CONSERVATIVE, &
         limiter   = MONOCENT, &           ! one of: minmod, monocent,...   !
         theta     = 1.2)                  ! optional parameter for limiter !

    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(SPHERICAL)
       x1 = 0.0
       x2 = SQRT(2.0)*RMAX
       y1 = 0.0
       y2 = 0.5*PI
       bc(WEST)  = REFLECTING
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = AXIS
       bc(NORTH) = REFLECTING
    CASE(CYLINDRICAL)
       x1 = 0.0
       x2 = RMAX
       y1 = 0.0
       y2 = RMAX
       bc(WEST)  = REFLECTING
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = AXIS
       bc(NORTH) = NO_GRADIENTS
    CASE(OBLATE_SPHEROIDAL)
       x1 = 0.0
       x2 = 2./GPAR**2 * (1.0 + SQRT(1.0+0.25*GPAR**4))
       x2 = 0.5*LOG(x2+SQRT(x2**2-1.0)) ! = 0.5*ACOSH(RMAX/GPAR)
       y1 = 0.0
       y2 = 0.5*PI
       bc(WEST)  = REFLECTING
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = REFLECTING
       bc(NORTH) = AXIS
    CASE(SINHSPHERICAL)
       x1 = 0.0
       x2 = SQRT(2.0)*RMAX/GPAR
       x2 = LOG(x2+SQRT(x2**2+1.0)) ! = ASINH(RMAX/GPAR)
       y1 = 0.0
       y2 = 0.5*PI
       bc(WEST)  = REFLECTING
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = AXIS
       bc(NORTH) = REFLECTING
    CASE DEFAULT
       CALL Error(Physics,"InitProgram","geometry not supported for this test")
    END SELECT
    ! mesh settings
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
         dtlimit  = 1.0E-8, &
         maxiter  = 10000000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc%pvar,Timedisc%cvar)

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


  SUBROUTINE InitData(Mesh,Physics,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)&
                      :: pvar,cvar
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(OUT)       :: pvar,cvar
    !------------------------------------------------------------------------!
    ! initial density and pressure
    FORALL (i=Mesh%IGMIN:Mesh%IGMAX,j=Mesh%JGMIN:Mesh%JGMAX)
       pvar(i,j,Physics%DENSITY) = RHO0 + RHO1*EXP(-LOG(2.0) &
            * ((Mesh%bccart(i,j,1)-R0)**2 + (Mesh%bccart(i,j,2)-Z0)**2)/RWIDTH**2)
       pvar(i,j,Physics%PRESSURE) = P0 + P1*EXP(-LOG(2.0) &
            * ((Mesh%bccart(i,j,1)-R0)**2 + (Mesh%bccart(i,j,2)-Z0)**2)/PWIDTH**2)
    END FORALL

    ! velocities in the x-y-plane
    pvar(:,:,Physics%XVELOCITY) = 0.
    pvar(:,:,Physics%YVELOCITY) = 0.

    ! rotational velocity
    pvar(:,:,Physics%ZVELOCITY) = OMEGA0 * Mesh%bhz(:,:)

    ! for specific angular momentum transport
    IF (GetType(Physics).EQ.EULER3D_ROTAMT) THEN
       ! specific angular momentum
       pvar(:,:,Physics%ZVELOCITY) = pvar(:,:,Physics%ZVELOCITY)*Mesh%bhz(:,:)
    END IF
    
    CALL Convert2Conservative(Physics,Mesh,pvar,cvar)

    CALL Info(Mesh, " DATA-----> initial condition: 3D Gaussian pulse")

  END SUBROUTINE InitData

END MODULE Init
