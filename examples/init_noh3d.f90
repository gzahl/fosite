!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_noh3d.f90                                                    #
!#                                                                           #
!# Copyright (C) 2008-2010                                                   #
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
! Program and data initialization for 3D Noh problem
! References:
! [1] Noh, W. F.: Errors for calculations of strong shocks using an artificial
!     viscosity and an artificial heat-flux, J. Comput. Phys. 72 (1987), 78-120
!     DOI: 10.1016/0021-9991(87)90074-X
! [2] Rider, W. J.: Revisiting wall heating, J. Comput. Phys. 162 (2000), 395-410
!     DOI: 10.1006/jcph.2000.6544
!----------------------------------------------------------------------------!
MODULE Init
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE timedisc_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! simulation parameters
  REAL, PARAMETER    :: TSIM    = 1.2      ! simulation stop time
  REAL, PARAMETER    :: GAMMA   = 5./3.    ! ratio of specific heats
  ! initial condition (dimensionless units)
  REAL, PARAMETER    :: RHO0    = 1.0      ! density
  REAL, PARAMETER    :: P0      = 1.0E-5   ! pressure
  REAL, PARAMETER    :: VR0     = -1.0     ! radial velocity
  ! mesh settings
!!$  INTEGER, PARAMETER :: MGEO = SPHERICAL   ! geometry
  INTEGER, PARAMETER :: MGEO = CYLINDRICAL
!!$  INTEGER, PARAMETER :: MGEO = TANCYLINDRICAL
!!$  INTEGER, PARAMETER :: MGEO = OBLATE_SPHEROIDAL !!! ONLY WORKING for P0 >= 1e-2 !!!
!!$  INTEGER, PARAMETER :: MGEO = SINHSPHERICAL
  INTEGER, PARAMETER :: XRES  = 100        ! x-resolution
  INTEGER, PARAMETER :: YRES  = 50         ! y-resolution
  REAL, PARAMETER    :: RMAX = 0.5         ! outer radius
  REAL, PARAMETER    :: GPAR = 0.2         ! geometry scaling parameter < RMAX
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 10          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'noh3d' 
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
    INTEGER           :: j
    INTEGER           :: bc(4)
    REAL              :: x1,x2,y1,y2
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!
    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER3D_ROTSYM, &
         gamma   = GAMMA, &                 ! ratio of specific heats        !
         dpmax   = 1.0E+12)                 ! for advanced time step control !

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)                 ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         ! IMPORTANT: always use primitive reconstruction for NOH problem
         variables = PRIMITIVE, &           ! vars. to use for reconstruction!
         limiter   = MONOCENT, &            ! one of: minmod, monocent,...   !
         theta     = 1.2)                   ! optional parameter for limiter !

    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(SPHERICAL)
       x1 = 0.0
       x2 = RMAX
       y1 = 0.0
       y2 = PI
       bc(WEST)  = REFLECTING
       bc(EAST)  = NOH3D
       bc(SOUTH) = AXIS
       bc(NORTH) = AXIS
    CASE(CYLINDRICAL)
       x1 = -RMAX
       x2 = RMAX
       y1 = 0.0
       y2 = RMAX
       bc(WEST)  = NOH3D
       bc(EAST)  = NOH3D
       bc(SOUTH) = AXIS
       bc(NORTH) = NOH3D
!!! FIXME: not working for low ambient pressure 
    CASE(OBLATE_SPHEROIDAL)
       x1 = 0.0
       x2 = RMAX/GPAR
       x2 = LOG(x2+SQRT(x2*x2+1.0)) ! = ASINH(RMAX/GPAR)
       y1 = 0.0
       y2 = 0.5*PI
       bc(WEST)  = REFLECTING
       bc(EAST)  = NOH3D
       bc(SOUTH) = REFLECTING
       bc(NORTH) = AXIS
    CASE(SINHSPHERICAL)
       x1 = 0.0
       x2 = RMAX/GPAR
       x2 = LOG(x2+SQRT(x2**2+1.0)) ! = ASINH(RMAX/GPAR)
       y1 = 0.0
       y2 = PI
       bc(WEST)  = REFLECTING
       bc(EAST)  = NOH3D
       bc(SOUTH) = AXIS
       bc(NORTH) = AXIS
    CASE DEFAULT
       CALL Error(Physics,"InitProgram","geometry not supported for 3D Noh problem")
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
         dtlimit  = 1.0E-10, &
         maxiter  = 100000)

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
    TYPE(Physics_TYP)  :: Physics
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Timedisc_TYP) :: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: vcart
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! initial condition
    Timedisc%pvar(:,:,Physics%DENSITY)   = RHO0
    Timedisc%pvar(:,:,Physics%PRESSURE)  = P0
    Timedisc%pvar(:,:,Physics%ZVELOCITY) = 0.

    ! cartesian components of the velocity
    vcart(:,:,1) = VR0 * Mesh%bccart(:,:,1) / &
         SQRT(Mesh%bccart(:,:,1)**2 + Mesh%bccart(:,:,2)**2)
    vcart(:,:,2) = vcart(:,:,1) * Mesh%bccart(:,:,2) / Mesh%bccart(:,:,1)
    CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,vcart,&
         Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY))

    ! supersonic inflow boundary conditions at outer boundaries;
    ! set boundary data equal to initial values in ghost cells
    IF (GetType(Timedisc%boundary(WEST)).EQ.NOH3D) &
       Timedisc%boundary(WEST)%data(:,:,:) = &
       Timedisc%pvar(Mesh%IGMIN:Mesh%IMIN-1,Mesh%JMIN:Mesh%JMAX,:)
    IF (GetType(Timedisc%boundary(EAST)).EQ.NOH3D) &
       Timedisc%boundary(EAST)%data(:,:,:) = &
       Timedisc%pvar(Mesh%IMAX+1:Mesh%IGMAX,Mesh%JMIN:Mesh%JMAX,:)
    IF (GetType(Timedisc%boundary(SOUTH)).EQ.NOH3D) &
       Timedisc%boundary(SOUTH)%data(:,:,:) = &
       Timedisc%pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JGMIN:Mesh%JMIN-1,:)
    IF (GetType(Timedisc%boundary(NORTH)).EQ.NOH3D) &
       Timedisc%boundary(NORTH)%data(:,:,:) = &
       Timedisc%pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+1:Mesh%JGMAX,:)

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: 3D Noh problem")
  END SUBROUTINE InitData

END MODULE Init
