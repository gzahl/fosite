!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_vortex2d.f90                                                 #
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
! Program and data initialization for 2D isentropic vortex
! References:
! [1] Yee, H. C. et al.: Low-dissipative high-order shock-capturing methods
!     using characteristic-based filters, J. Comput. Phys. 150 (1999), 199-238
!     DOI: 10.1006/jcph.1998.6177
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
  ! simulation parameters
  REAL, PARAMETER    :: TSIM    = 100.0    ! simulation stop time
  REAL, PARAMETER    :: GAMMA   = 1.4      ! ratio of specific heats
  ! initial condition (dimensionless units)
  REAL, PARAMETER    :: RHOINF  = 1.       ! ambient density
  REAL, PARAMETER    :: PINF    = 1.       ! ambient pressure
  REAL, PARAMETER    :: VSTR    = 5.0      ! nondimensional vortex strength
  REAL, PARAMETER    :: UINF    = 0.1      ! cartesian components of constant
  REAL, PARAMETER    :: VINF    = 0.1      !   global velocity field
  REAL, PARAMETER    :: X0      = 0.0      ! vortex position (cart. coords.)
  REAL, PARAMETER    :: Y0      = 0.0
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CARTESIAN   ! geometry
!!$  INTEGER, PARAMETER :: MGEO = POLAR    
!!$  INTEGER, PARAMETER :: MGEO = LOGPOLAR
!!$  INTEGER, PARAMETER :: MGEO = TANPOLAR
!!$  INTEGER, PARAMETER :: MGEO = SINHPOLAR
  INTEGER, PARAMETER :: XRES = 100         ! x-resolution
  INTEGER, PARAMETER :: YRES = 100         ! y-resolution
  REAL, PARAMETER    :: RMIN = 1.0E-2      ! inner radius for polar grids
  REAL, PARAMETER    :: RMAX = 5.0         ! outer radius
  REAL, PARAMETER    :: GPAR = 1.0         ! geometry scaling parameter     !
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 100         ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'vortex2d' 
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
    INTEGER           :: bc(4)
    REAL              :: x1,x2,y1,y2
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!
    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER2D, &
         gamma   = GAMMA, &                 ! ratio of specific heats        !
         dpmax   = 1.0)                     ! for advanced time step control !

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)                 ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         variables = CONSERVATIVE, &        ! vars. to use for reconstruction!
         limiter   = MONOCENT, &            ! one of: minmod, monocent,...   !
         theta     = 1.2)                   ! optional parameter for limiter !

    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(CARTESIAN)
       x1 =-RMAX
       x2 = RMAX
       y1 =-RMAX 
       y2 = RMAX
       bc(WEST)  = PERIODIC
       bc(EAST)  = PERIODIC
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(POLAR)
       x1 = RMIN
       x2 = RMAX
       y1 = 0.0 
       y2 = 2*PI       
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(LOGPOLAR)
       x1 = LOG(RMIN/GPAR)
       x2 = LOG(RMAX/GPAR)
       y1 = 0.0 
       y2 = 2*PI       
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(TANPOLAR)
       x1 = ATAN(RMIN/GPAR)
       x2 = ATAN(RMAX/GPAR)
       y1 = 0.0 
       y2 = 2*PI       
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(SINHPOLAR)
       x1 = RMIN/GPAR
       x1 = LOG(x1+SQRT(1.0+x1*x1))  ! = ASINH(RMIN/GPAR))
       x2 = RMAX/GPAR
       x2 = LOG(x2+SQRT(1.0+x2*x2))  ! = ASINH(RMAX/GPAR))
       y1 = 0.0 
       y2 = 2*PI       
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE DEFAULT
       CALL Error(Physics,"InitProgram","mesh geometry not supported for 2D isentropic vortex")
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

    CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = bc(WEST), &
         eastern  = bc(EAST), &
         southern = bc(SOUTH), &
         northern = bc(NORTH))
 
    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = TSIM, &
         dtlimit  = 1.0E-4, &
         maxiter  = 1000000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc)

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
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: vxy
    REAL              :: x,y,r2
    REAL              :: du,dv
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! initial condition
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          x  = Mesh%bccart(i,j,1) - X0
          y  = Mesh%bccart(i,j,2) - Y0
          r2 = x*x + y*y
          du = -0.5*VSTR/PI*y*EXP(0.5*(1.-r2))
          dv = -du*x/y
          ! cartesian velocity components
          vxy(i,j,1) = UINF + du
          vxy(i,j,2) = VINF + dv
          ! density
          ! ATTENTION: there's a factor of 1/PI missing in the density
          ! formula  eq. (3.3) in [1]
          Timedisc%pvar(i,j,Physics%DENSITY) = RHOINF * (1.0 - (GAMMA-1.0)*VSTR**2 / &
               (8*PI**2*GAMMA) * EXP(1.-r2) )**(1./(GAMMA-1.))
          ! pressure
          Timedisc%pvar(i,j,Physics%PRESSURE) = PINF &
               * (Timedisc%pvar(i,j,Physics%DENSITY)/RHOINF)**GAMMA
       END DO
    END DO

    ! transform velocities
    CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,vxy(:,:,:),&
         Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY))

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh," DATA-----> initial condition: 2D isentropic vortex")

  END SUBROUTINE InitData

END MODULE Init
