!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_riemann1d.f90                                                #
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
! Program and data initialization for 1D Riemann problems
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
  INTEGER :: testnum
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
    TYPE(FileIO_TYP)  :: Datafile
    TYPE(FileIO_TYP)  :: Logfile
    !------------------------------------------------------------------------!
    ! Local variable declaration
    CHARACTER(LEN=256):: ofname,lfname
    INTEGER           :: ierror
    REAL              :: test_stoptime, test_gamma, test_cs
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!

    ! 1D test problems (for initial conditions see InitData)
    !   1: Toro test no. 1
    !   2: Toro test no. 2
    !   3: Toro test no. 3
    !   4: Sod problem
    !   5: Noh problem
    !   6: simple isothermal Riemann problem
    testnum = 4

    SELECT CASE(testnum)
    CASE(1) ! Toro test 1
       ofname  = "toro1"
       lfname  = "toro1log"
       test_gamma    = 1.4
       test_stoptime = 0.2
    CASE(2) ! Toro test 2
       ofname  = "toro2"
       lfname  = "toro2log"
       test_gamma    = 1.4
       test_stoptime = 0.15
    CASE(3) ! Toro test 3
       ofname  = "toro3"
       lfname  = "toro3log"
       test_gamma    = 1.4
       test_stoptime = 0.012
    CASE(4) ! Sod problem
       ofname  = "sod"
       lfname  = "sodlog"
       test_gamma    = 1.4
       test_stoptime = 0.25
    CASE(5) ! Noh problem
       ofname  = "noh"
       lfname  = "nohlog"
       test_gamma    = 5./3.
       test_stoptime = 1.0
    CASE(6) ! isothermal shock tube
       ofname  = "isotherm"
       lfname  = "isothermlog"
       test_cs = 1.0       ! isothermal sound speed
       test_stoptime = 0.25
    CASE DEFAULT
       CALL Error(Mesh,"InitProgram", "Test problem number should be 1,2,3,4,5 or 6")
    END SELECT

    ! physics settings
    IF (testnum.EQ.6) THEN
       CALL InitPhysics(Physics, &
            problem = EULER2D_ISOTHERM, &
            cs      = test_cs)
    ELSE   
       CALL InitPhysics(Physics, &
            problem = EULER2D, &
            gamma   = test_gamma, &    ! ratio of specific heats        !
            dpmax   = 1.0E+06)         ! for advanced time step control !
    END IF

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)         ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         variables = CONSERVATIVE,& ! vars. to use for reconstruction!
         limiter   = MONOCENT, &    ! one of: minmod, monocent,...   !
         theta     = 1.2)           ! optional parameter for limiter !

    ! mesh settings
    CALL InitMesh(Mesh,Fluxes, &
         geometry = CARTESIAN, &
             inum = 200,&           ! resolution in x and            !
             jnum = 1, &            !   y direction                  !             
             xmin = 0.0, &
             xmax = 1.0, &
             ymin = 0.0, &
             ymax = 1.0)

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = test_stoptime, &
         dtlimit  = 1.0E-10, &
         maxiter  = 100000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc)

    ! boundary conditions
    CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = NO_GRADIENTS, &
         eastern  = NO_GRADIENTS, &
         southern = NO_GRADIENTS, &
         northern = NO_GRADIENTS)

    ! initialize log input/output
    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc, &
         fileformat = BINARY, &
#ifdef PARALLEL
         filename   = "/tmp/" // TRIM(lfname), &
#else
         filename   = TRIM(lfname), &
#endif
         filecycles = 1)

    ! initialize data input/output
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc, &
         fileformat = GNUPLOT, &
#ifdef PARALLEL
         filename   = "/tmp/" // TRIM(ofname), &
#else
         filename   = TRIM(ofname), &
#endif
         filecycles = 0, &
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
    REAL              :: x0, rho_l, rho_r, u_l, u_r, p_l, p_r
    CHARACTER(LEN=64) :: teststr
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!

    SELECT CASE(testnum)
    CASE(1)
       teststr = "1D Riemannn problem: Toro test no. 1"
       x0    = 0.3
       rho_l = 1.0
       rho_r = 0.125
       u_l   = 0.75
       u_r   = 0.0
       p_l   = 1.0
       p_r   = 0.1
    CASE(2)
       teststr = "1D Riemannn problem: Toro test no. 2"
       x0    = 0.5
       rho_l = 1.0
       rho_r = 1.0
       u_l   = -2.0
       u_r   = 2.0
       p_l   = 0.4
       p_r   = 0.4       
    CASE(3)
       teststr = "1D Riemannn problem: Toro test no. 3"
       x0    = 0.5
       rho_l = 1.0
       rho_r = 1.0
       u_l   = 0.0
       u_r   = 0.0
       p_l   = 1000.
       p_r   = 0.01
    CASE(4)
       teststr = "1D Riemannn problem: Sod problem"
       x0    = 0.5
       rho_l = 1.0
       rho_r = 0.125
       u_l   = 0.0
       u_r   = 0.0
       p_l   = 1.0
       p_r   = 0.1
    CASE(5)
       teststr = "1D Riemannn problem: Noh problem"
       x0    = 0.5
       rho_l = 1.0
       rho_r = 1.0
       u_l   = 1.0
       u_r   = -1.0
       p_l   = 1.0E-05
       p_r   = 1.0E-05
    CASE(6)
       teststr = "1D Riemannn problem: isothermal shock tube"
       x0    = 0.5
       rho_l = 1.0
       rho_r = 0.125
       u_l   = 0.0
       u_r   = 0.0
    CASE DEFAULT
       CALL Error(Mesh,"InitData", "Test problem should be 1,2,3,4, or 5")
    END SELECT
    
    IF (Mesh%INUM.GT.Mesh%JNUM) THEN
       WHERE (Mesh%bcenter(:,:,1).LT.x0)
          Timedisc%pvar(:,:,Physics%DENSITY)   = rho_l
          Timedisc%pvar(:,:,Physics%XVELOCITY) = u_l
          Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.
       ELSEWHERE
          Timedisc%pvar(:,:,Physics%DENSITY)   = rho_r
          Timedisc%pvar(:,:,Physics%XVELOCITY) = u_r
          Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.
       END WHERE
    ELSE
        WHERE (Mesh%bcenter(:,:,2).LT.x0)
          Timedisc%pvar(:,:,Physics%DENSITY)   = rho_l
          Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
          Timedisc%pvar(:,:,Physics%YVELOCITY) = u_l
       ELSEWHERE
          Timedisc%pvar(:,:,Physics%DENSITY)   = rho_r
          Timedisc%pvar(:,:,Physics%XVELOCITY) = 0. 
          Timedisc%pvar(:,:,Physics%YVELOCITY) = u_r
       END WHERE
    END IF
    
    IF (GetType(Physics).EQ.EULER2D) THEN
       IF (Mesh%INUM.GT.Mesh%JNUM) THEN
          WHERE (Mesh%bcenter(:,:,1).LT.x0)
             Timedisc%pvar(:,:,Physics%PRESSURE)  = p_l
          ELSEWHERE
             Timedisc%pvar(:,:,Physics%PRESSURE)  = p_r
          END WHERE
       ELSE
          WHERE (Mesh%bcenter(:,:,2).LT.x0)
             Timedisc%pvar(:,:,Physics%PRESSURE)  = p_l
          ELSEWHERE
             Timedisc%pvar(:,:,Physics%PRESSURE)  = p_r
          END WHERE
       END IF
    END IF

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: " // TRIM(teststr))

  END SUBROUTINE InitData

END MODULE Init
