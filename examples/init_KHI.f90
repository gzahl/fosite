!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_KHI.f90                                                      #
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
! Program and data initialization for Kelvin-Helmholtz instability
! References: 
! [1] Lord Kelvin (William Thomson) (1871). "Hydrokinetic solutions and 
!     observations" Philosophical Magazine 42: 362-377 
! [2] Hermann von Helmholtz (1868). "On the discontinuous movements of fluids"
!     Monthly Reports of the Royal Prussian Academy of Philosophy in Berlin 23: 215-228
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
  REAL, PARAMETER    :: TSIM  = 5.0        ! simulation time
  REAL, PARAMETER    :: GAMMA = 1.4        ! ratio of specific heats
  REAL, PARAMETER    :: ETA   = 0.0        ! dynamic viscosity (0.0 disables)           
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CARTESIAN   ! geometry of the mesh
  INTEGER, PARAMETER :: XRES = 100         ! resolution
  INTEGER, PARAMETER :: YRES = 100
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 100         ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'KHI' 
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

    ! mesh settings
    CALL InitMesh(Mesh,Fluxes, &
         geometry = MGEO, &
             inum = XRES, &                 ! resolution in x and            !
             jnum = YRES, &                 !   y direction                  !             
             xmin = -0.5, &
             xmax = 0.5, &
             ymin = -0.5, &
             ymax = 0.5)

    ! boundary conditions
    CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = PERIODIC, &
         eastern  = PERIODIC, &
         southern = PERIODIC, &
         northern = PERIODIC)

    ! viscosity source term
    IF (ETA.GT.TINY(ETA)) THEN
       CALL InitSources(Physics%sources,Mesh,Fluxes,Physics,Timedisc%boundary, &
            stype    = VISCOSITY, &
            vismodel = MOLECULAR, &
            dynconst = ETA, &
           bulkconst = -2./3.*ETA)
    END IF

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = TSIM, &
         dtlimit  = 1.0E-4, &
         maxiter  = 100000)

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
    INTEGER           :: i,j
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: dv
    REAL              :: rho0, rho1, v0, v1, P0, P1
    REAL              :: ylen
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    
    ! outer regions
    rho0 = 1.0
    v0   = -0.5
    P0   = 2.5

    ! inner region
    rho1 = 2.0
    v1   = 0.5
    P1   = P0

    ! y-velocity vanishes everywhere
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.

    ! extent along the the y-direction
    ylen = ABS(Mesh%ymax-Mesh%ymin)

    WHERE ((Mesh%bcenter(:,:,2).LT.(Mesh%ymin+0.25*ylen)).OR. &
         (Mesh%bcenter(:,:,2).GT.(Mesh%ymin+0.75*ylen)))
       Timedisc%pvar(:,:,Physics%DENSITY) = rho0
       Timedisc%pvar(:,:,Physics%XVELOCITY) = v0
       Timedisc%pvar(:,:,Physics%PRESSURE) = P0
    ELSEWHERE
       Timedisc%pvar(:,:,Physics%DENSITY) = rho1
       Timedisc%pvar(:,:,Physics%XVELOCITY) = v1
       Timedisc%pvar(:,:,Physics%PRESSURE) = P1
    END WHERE
       
    ! add velocity perturbations
    CALL RANDOM_SEED
    CALL RANDOM_NUMBER(dv)
    Timedisc%pvar(:,:,Physics%XVELOCITY) = Timedisc%pvar(:,:,Physics%XVELOCITY) &
         + (dv(:,:,1)-0.5)*0.02
    Timedisc%pvar(:,:,Physics%YVELOCITY) = Timedisc%pvar(:,:,Physics%YVELOCITY) &
         + (dv(:,:,2)-0.5)*0.02

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: " // &
         "Kelvin-Helmholtz instability")

  END SUBROUTINE InitData

END MODULE Init
