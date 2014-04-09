!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_poiseuille.f90                                               #
!#                                                                           #
!# Copyright (C) 2006-2010                                                   #
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
! Program and data initialization for a test of Hagen-Poiseuille equ. in a tube
!----------------------------------------------------------------------------!
MODULE Init
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE sources_generic
  USE fileio_generic
  USE timedisc_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! simulation parameters
  REAL, PARAMETER :: TSIM  = 1.0                ! simulation time            !
  REAL, PARAMETER :: RE    = 4.375              ! Reynolds Number            !
  REAL, PARAMETER :: PIN   = 24.0               ! inflow pressure            !
  REAL, PARAMETER :: POUT  = 23.0               ! outflow pressure           !
  REAL, PARAMETER :: RHO0  = 1.0                ! initial density            !
  ! mesh settings
  INTEGER, PARAMETER :: RRES = 40               ! radial resolution          !
  INTEGER, PARAMETER :: ZRES = 80               ! resolution along the tube  !
  REAL, PARAMETER    :: LTUBE = 10.0            ! length of the tube         !
  REAL, PARAMETER    :: RTUBE = 1.0             ! radius of the tube         !
  !--------------------------------------------------------------------------!
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 100         ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'poiseuille' 
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
    REAL              :: dvis, bvis
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!

    ! compute viscosity constants
    ! dynamic viscosity
    dvis = SQRT(0.25 * (PIN-POUT) * RTUBE**3 * RHO0 / LTUBE / RE)
    ! bulk viscosity
    bvis = -2./3. * dvis

    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER3D_ROTSYM, &
         gamma   = 1.4, &           ! ratio of specific heats        !
         dpmax   = 1.0)             ! for advanced time step control !

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)         ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         variables = PRIMITIVE, &! vars. to use for reconstruction!
         limiter   = MONOCENT, &    ! one of: minmod, monocent,...   !
         theta     = 1.2)           ! optional parameter for limiter !

    ! mesh settings
    CALL InitMesh(Mesh,Fluxes, &
         geometry = CYLINDRICAL, &
             inum = ZRES, &         ! resolution in x and            !
             jnum = RRES, &         !   y direction                  !             
             xmin = 0.0, &
             xmax = LTUBE, &
             ymin = 0.0, &
             ymax = RTUBE)

    ! boundary conditions
    CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = FIXED, &
         eastern  = FIXED, &
         southern = AXIS, &
         northern = NOSLIP)

    ! viscosity source term
    CALL InitSources(Physics%sources,Mesh,Fluxes,Physics,Timedisc%Boundary, &
          stype    = VISCOSITY, &
          vismodel = MOLECULAR, &
          dynconst = dvis, &
         bulkconst = bvis)


    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = TSIM, &
         dtlimit  = 1.0E-8, &
         maxiter  = 1000000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc)

    ! initialize log input/output
    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,&
         fileformat = BINARY, &
         filename   = TRIM(ODIR) // TRIM(OFNAME) // 'log', &
         filecycles = 1)

    ! initialize data input/output
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc, &
         fileformat = VTK, &
!!$         fileformat = GNUPLOT, &
!!$         filecycles = 0, &
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
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: dv
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! initial condition    
    Timedisc%pvar(:,:,Physics%DENSITY)   = RHO0
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%ZVELOCITY) = 0. 
    Timedisc%pvar(:,:,Physics%PRESSURE) = PIN + (POUT-PIN)/LTUBE &
         * Mesh%bccart(:,:,2)

    ! fixed boundary conditions 
    ! inflow
    IF (GetType(Timedisc%Boundary(WEST)).EQ.FIXED) THEN
       Timedisc%Boundary(WEST)%data(:,:,Physics%DENSITY)    = RHO0
       Timedisc%Boundary(WEST)%data(:,:,Physics%XVELOCITY)  = 0.0
       Timedisc%Boundary(WEST)%data(:,:,Physics%YVELOCITY)  = 0.0
       Timedisc%Boundary(WEST)%data(:,:,Physics%ZVELOCITY)  = 0.0
       Timedisc%Boundary(WEST)%data(:,:,Physics%PRESSURE)   = PIN
       ! imposed density, pressure and tangential velocities;
       ! extrapolated normal velocity
       Timedisc%Boundary(WEST)%fixed(:,Physics%DENSITY)   = .TRUE.
       Timedisc%Boundary(WEST)%fixed(:,Physics%XVELOCITY) = .FALSE.
       Timedisc%Boundary(WEST)%fixed(:,Physics%YVELOCITY) = .TRUE.
       Timedisc%Boundary(WEST)%fixed(:,Physics%ZVELOCITY) = .TRUE.
       Timedisc%Boundary(WEST)%fixed(:,Physics%PRESSURE)  = .TRUE.
    ENDIF
    ! outflow
    IF (GetType(Timedisc%Boundary(EAST)).EQ.FIXED) THEN
       Timedisc%Boundary(EAST)%data(:,:,Physics%DENSITY)    = RHO0
       Timedisc%Boundary(EAST)%data(:,:,Physics%XVELOCITY)  = 0.0
       Timedisc%Boundary(EAST)%data(:,:,Physics%YVELOCITY)  = 0.0
       Timedisc%Boundary(EAST)%data(:,:,Physics%ZVELOCITY)  = 0.0
       Timedisc%Boundary(EAST)%data(:,:,Physics%PRESSURE)   = POUT
       ! imposed pressure; 
       ! extrapolated density, tangential and normal velocities
       Timedisc%Boundary(EAST)%fixed(:,Physics%DENSITY)   = .FALSE.
       Timedisc%Boundary(EAST)%fixed(:,Physics%XVELOCITY) = .FALSE.
       Timedisc%Boundary(EAST)%fixed(:,Physics%YVELOCITY) = .FALSE.
       Timedisc%Boundary(EAST)%fixed(:,Physics%ZVELOCITY) = .FALSE.
       Timedisc%Boundary(EAST)%fixed(:,Physics%PRESSURE)  = .TRUE.
    ENDIF

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: " // &
         "tube with pressure gradient")

  END SUBROUTINE InitData

END MODULE Init
