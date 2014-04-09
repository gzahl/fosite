!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_riemann1d.f90                                                #
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
! Program and data initialization for 1D Riemann problems
! References:
! [1] Toro, E. F.: Riemann Solvers and Numerical Methods for Fluid Dynamics,
!     A Practical Introduction, Springer-Verlag 1999, 2nd ed., Chapter 4.3.3
! [2] Sod, G. A.: A survey of several finite difference methods for systems of
!     nonlinear hyperbolic conservation laws, J. Comput. Phys. 27 (1978), 1-31
!     DOI: 10.1016/0021-9991(78)90023-2
! [3] Noh, W. F.: Errors for calculations of strong shocks using an artificial
!     viscosity and an artificial heat-flux, J. Comput. Phys. 72 (1987), 78-120
!     DOI: 10.1016/0021-9991(87)90074-X
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
  ! initial condition, one of
  !   1: Sod problem, see ref. [1],[2]
  !   2: Toro test no. 2, see ref. [1]
  !   3: Toro test no. 3, see ref. [1]
  !   4: Toro test no. 4, see ref. [1]
  !   5: Toro test no. 5, see ref. [1]
  !   6: Noh problem, see ref. [3]
  !   7: simple isothermal Riemann problem
  INTEGER, PARAMETER :: ICNUM = 1
  CHARACTER(LEN=256) :: TESTSTR          ! test description
  ! mesh settings
  INTEGER, PARAMETER :: XRES = 100       ! x-resolution
  INTEGER, PARAMETER :: YRES = 1         ! y-resolution
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 1         ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &        ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256) :: OFNAME           ! output file name
  ! some global constants
  REAL               :: GAMMA            ! ratio of specific heats
  REAL               :: TSIM             ! simulation time
  REAL               :: CSISO            ! isothermal sound speed (test no. 6)
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
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!
    SELECT CASE(ICNUM)
    CASE(1)
       TESTSTR= "1D Riemannn problem: Sod shock tube"
       OFNAME = "sod"
       GAMMA  = 1.4
       TSIM   = 0.25
    CASE(2)
       TESTSTR= "1D Riemannn problem: Toro test no. 2"
       OFNAME = "toro2"
       GAMMA  = 1.4
       TSIM   = 0.15
    CASE(3)
       TESTSTR= "1D Riemannn problem: Toro test no. 3"
       OFNAME = "toro3"
       GAMMA  = 1.4
       TSIM   = 0.012
    CASE(4)
       TESTSTR= "1D Riemannn problem: Toro test no. 4"
       OFNAME = "toro4"
       GAMMA  = 1.4
       TSIM   = 0.035
    CASE(5)
       TESTSTR= "1D Riemannn problem: Toro test no. 5"
       OFNAME = "toro5"
       GAMMA  = 1.4
       TSIM   = 0.035
    CASE(6)
       TESTSTR= "1D Riemannn problem: Noh problem"
       OFNAME = "noh"
       GAMMA  = 5./3.
       TSIM   = 1.0
    CASE(7)
       TESTSTR= "1D Riemannn problem: Isothermal shock tube"
       OFNAME = "isotherm"
       CSISO  = 1.0
       TSIM   = 0.25
    CASE DEFAULT
       CALL Error(Mesh,"InitProgram", "Test problem number should be 1,2,3,4,5,6 or 7")
    END SELECT

    ! physics settings
    IF (ICNUM.EQ.7) THEN
       CALL InitPhysics(Physics, &
            problem = EULER2D_ISOTHERM, &
            cs      = CSISO)
    ELSE   
       CALL InitPhysics(Physics, &
            problem = EULER2D, &
            gamma   = GAMMA, &         ! ratio of specific heats        !
            dpmax   = 1.0E+06)         ! for advanced time step control !
    END IF

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)            ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         variables = CONSERVATIVE,&    ! vars. to use for reconstruction!
         limiter   = MONOCENT, &       ! one of: minmod, monocent,...   !
         theta     = 1.2)              ! optional parameter for limiter !

    ! mesh settings
    CALL InitMesh(Mesh,Fluxes, &
         geometry = CARTESIAN, &
             inum = XRES,&             ! resolution in x and            !
             jnum = YRES, &            !   y direction                  !             
             xmin = 0.0, &
             xmax = 1.0, &
             ymin = 0.0, &
             ymax = 1.0)

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.3, &
         stoptime = TSIM, &
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
!!$    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,&
!!$         fileformat = BINARY, &
!!$         filename   = TRIM(ODIR) // TRIM(OFNAME) // 'log', &
!!$         filecycles = 1)

    ! initialize data input/output
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc, &
         fileformat = GNUPLOT, &
         filename   = TRIM(ODIR) // TRIM(OFNAME), &
         filecycles = 0, &
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
    REAL              :: x0, rho_l, rho_r, u_l, u_r, p_l, p_r
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!

    SELECT CASE(ICNUM)
    CASE(1) ! Sod shock tube (Toro's test no. 1)
       x0    = 0.5
       rho_l = 1.0
       rho_r = 0.125
       u_l   = 0.0
       u_r   = 0.0
       p_l   = 1.0
       p_r   = 0.1
    CASE(2) ! Toro's test no. 2
       x0    = 0.5
       rho_l = 1.0
       rho_r = 1.0
       u_l   = -2.0
       u_r   = 2.0
       p_l   = 0.4
       p_r   = 0.4       
    CASE(3) ! Toro's test no. 3
       x0    = 0.5
       rho_l = 1.0
       rho_r = 1.0
       u_l   = 0.0
       u_r   = 0.0
       p_l   = 1000.
       p_r   = 0.01
    CASE(4) ! Toro's test no. 4
       x0    = 0.5
       rho_l = 1.0
       rho_r = 1.0
       u_l   = 0.0
       u_r   = 0.
       p_l   = 0.01
       p_r   = 100.0
    CASE(5) ! Toro's test no. 5
       x0    = 0.3
       rho_l = 5.99924
       rho_r = 5.99242
       u_l   = 19.5975
       u_r   = -6.19633
       p_l   = 460.894
       p_r   = 46.0950
    CASE(6) ! Noh problem
       x0    = 0.5
       rho_l = 1.0
       rho_r = 1.0
       u_l   = 1.0
       u_r   = -1.0
       p_l   = 1.0E-05
       p_r   = 1.0E-05
    CASE(7) ! isothermal shock tube
       x0    = 0.5
       rho_l = 1.0
       rho_r = 0.125
       u_l   = 0.0
       u_r   = 0.0
    CASE DEFAULT
       CALL Error(Mesh,"InitData", "Test problem should be 1,2,3,4,5,6 or 7")
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
    CALL Info(Mesh, " DATA-----> initial condition: " // TRIM(TESTSTR))

  END SUBROUTINE InitData

END MODULE Init
