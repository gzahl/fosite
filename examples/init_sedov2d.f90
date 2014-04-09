!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_sedov2d.f90                                                  #
!#                                                                           #
!# Copyright (C) 2006-2008                                                   #
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
! Program and data initialization for 2D Sedov explosion
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
    INTEGER           :: geometry
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!

    ! set the geometry
    geometry = CARTESIAN
!    geometry = POLAR

    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER2D, &
         gamma   = 1.4, &           ! ratio of specific heats        !
         dpmax   = 100.0)           ! for advanced time step control !

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)         ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         variables = CONSERVATIVE, &! vars. to use for reconstruction!
         limiter   = MONOCENT, &    ! one of: minmod, monocent,...   !
         theta     = 1.2)           ! optional parameter for limiter !

    SELECT CASE(geometry)
    CASE(CARTESIAN)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes, &
            geometry = CARTESIAN, &
                inum = 100, &       ! resolution in x and            !
                jnum = 100, &       !   y direction                  !             
                xmin = -0.3, &
                xmax = 0.3, &
                ymin = -0.3, &
                ymax = 0.3)
       ! boundary conditions
       CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = NO_GRADIENTS, &
         eastern  = NO_GRADIENTS, &
         southern = NO_GRADIENTS, &
         northern = NO_GRADIENTS)
    CASE(POLAR)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes, &
            geometry = POLAR, &
                inum = 50, &        ! resolution in x and            !
                jnum = 60, &        !   y direction                  !             
                xmin = 0.0, &
                xmax = 0.3, &
                ymin = 0.0, &
                ymax = 2*PI)
       ! boundary conditions
       CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = REFLECTING, &
         eastern  = NO_GRADIENTS, &
         southern = PERIODIC, &
         northern = PERIODIC)
    CASE DEFAULT
       CALL Error(Physics,"InitProgram",&
            "geometry should be either cartesian or polar")
    END SELECT

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = 0.05, &
         dtlimit  = 1.0E-13, &
         maxiter  = 100000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc)

    ! initialize log input/output
    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,&
         fileformat = BINARY, &
#ifdef PARALLEL
         filename   = "/tmp/sedov2dlog", &
#else
         filename   = "sedov2dlog", &
#endif
         filecycles = 1)

    ! initialize data input/output
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc, &
         fileformat = GNUPLOT, &
#ifdef PARALLEL
         filename   = "/tmp/sedov2d", &
#else
         filename   = "sedov2d", &
#endif
         count      = 10)

  END SUBROUTINE InitProgram


  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: n
    REAL              :: dr,rho0,P0,P1,E
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: cart
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! 2D Sedov explosion parameters
    n    = 2       ! 2D
    E    = 1.0     ! total energy input
    rho0 = 1.0     ! ambient density
    P0   = 1.0E-05 ! ambient pressure

    ! spatial with of the puls should be at least 5 cells
    ! be careful, if you wish to compare with the polar
    ! results dr should be the same
    dr = 0.03
    ! peak pressure
    P1 = 3.*(Physics%gamma-1.0) * E / ((n + 1) * PI * dr**n)

    ! cartesian coordinates
    CALL Convert2Cartesian(Mesh%geometry,Mesh%bcenter,cart)

    WHERE ((cart(:,:,1)**2 + cart(:,:,2)**2).LE.dr**2)
       Timedisc%pvar(:,:,Physics%DENSITY)   = rho0
       Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
       Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.
       Timedisc%pvar(:,:,Physics%PRESSURE)  = P1
    ELSEWHERE
       Timedisc%pvar(:,:,Physics%DENSITY)   = rho0
       Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
       Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.
       Timedisc%pvar(:,:,Physics%PRESSURE)  = P0
    END WHERE
     
    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: 2D Sedov explosion")

  END SUBROUTINE InitData


END MODULE Init
