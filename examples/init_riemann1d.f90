!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_riemann1d.f90                                                #
!#                                                                           #
!# Copyright (C) 2006 Tobias Illenseer <tillense@ita.uni-heidelberg.de>      #
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
  USE output_generic
  USE logio_generic
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

  SUBROUTINE InitProgram(Mesh,Physics,Fluxes,Timedisc,Output,Logio)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Timedisc_TYP):: Timedisc
    TYPE(Output_TYP)  :: Output
    TYPE(Logio_TYP)   :: Logio
    !------------------------------------------------------------------------!
    ! Local variable declaration
    CHARACTER(LEN=256):: ofname,lfname
    REAL              :: test_stoptime, test_gamma
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Output,Logio
    !------------------------------------------------------------------------!

    ! 1D test problems (for initial conditions see InitData)
    !   1: Toro test no. 1
    !   2: Toro test no. 2
    !   3: Toro test no. 3
    !   4: Sod problem
    !   5: Noh problem
    testnum = 4

    SELECT CASE(testnum)
    CASE(1) ! Toro test 1
       ofname  = "toro1.dat"
       lfname  = "toro1.log"
       test_gamma    = 1.4
       test_stoptime = 0.2
    CASE(2) ! Toro test 2
       ofname  = "toro2.dat"
       lfname  = "toro2.log"
       test_gamma    = 1.4
       test_stoptime = 0.15
    CASE(3) ! Toro test 3
       ofname  = "toro3.dat"
       lfname  = "toro3.log"
       test_gamma    = 1.4
       test_stoptime = 0.012
    CASE(4) ! Sod problem
       ofname  = "sod.dat"
       lfname  = "sod.log"
       test_gamma    = 1.4
       test_stoptime = 0.25
    CASE(5) ! Noh problem
       ofname  = "noh.dat"
       lfname  = "noh.log"
       test_gamma    = 5./3.
       test_stoptime = 1.0
    CASE DEFAULT
       PRINT *, "ERROR in InitProgram: Test problem number should be 1,2,3,4 or 5"
       STOP
    END SELECT

    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER2D, &
         gamma   = test_gamma, &    ! ratio of specific heats        !
         dpmax   = 1.0E+06)         ! for advanced time step control !

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)         ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         variables = CONSERVATIVE,& ! vars. to use for reconstruction!
         limiter   = MONOCENT, &    ! one of: minmod, monocent,...   !
         theta     = 1.3)           ! optional parameter for limiter !

    ! mesh settings
    CALL InitMesh(Mesh,Fluxes, &
         geometry = CARTESIAN, &
             inum = 100, &          ! resolution in x and            !
             jnum = 1, &            !   y direction                  !             
             xmin = 0.0, &
             xmax = 1.0, &
             ymin = -0.1, &
             ymax = 0.1)

    ! boundary conditions
    CALL InitBoundary(Mesh%boundary,Mesh,Physics,NO_GRADIENTS,WEST)
    CALL InitBoundary(Mesh%boundary,Mesh,Physics,NO_GRADIENTS,EAST)
    CALL InitBoundary(Mesh%boundary,Mesh,Physics,PERIODIC,SOUTH)
    CALL InitBoundary(Mesh%boundary,Mesh,Physics,PERIODIC,NORTH)

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = test_stoptime, &
         dtlimit  = 1.0E-10, &
         maxiter  = 10000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc%pvar,Timedisc%cvar)

    ! initialize log input/output
    CALL InitLogio(Logio,Mesh,Physics,Timedisc, &
         logformat = NOLOG, &
         filename  = lfname, &
         logdt     = 300)

    ! set parameters for data output
    CALL InitOutput(Output,Mesh,Physics,Timedisc,&
         filetype  = GNUPLOT, &
         filename  = ofname, &
         mode      = OVERWRITE, &
         starttime = 0.0, &
         stoptime  = Timedisc%stoptime, &
         count     = 1)

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
    REAL              :: rho_l, rho_r, u_l, u_r, p_l, p_r
    CHARACTER(LEN=64) :: teststr
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(OUT)       :: pvar,cvar
    !------------------------------------------------------------------------!

    pvar(:,:,:) = 0.

    SELECT CASE(testnum)
    CASE(1)
       teststr = "1D Riemannn problem: Toro test no. 1"
       j = 3 / 10 * Mesh%INUM
       rho_l = 1.0
       rho_r = 0.125
       u_l   = 0.75
       u_r   = 0.0
       p_l   = 1.0
       p_r   = 0.1
    CASE(2)
       teststr = "1D Riemannn problem: Toro test no. 2"
       j = Mesh%INUM / 2
       rho_l = 1.0
       rho_r = 1.0
       u_l   = -2.0
       u_r   = 2.0
       p_l   = 0.4
       p_r   = 0.4       
    CASE(3)
       teststr = "1D Riemannn problem: Toro test no. 3"
       j = Mesh%INUM / 2
       rho_l = 1.0
       rho_r = 1.0
       u_l   = 0.0
       u_r   = 0.0
       p_l   = 1000.
       p_r   = 0.01
    CASE(4)
       teststr = "1D Riemannn problem: Sod problem"
       j = Mesh%INUM / 2
       rho_l = 1.0
       rho_r = 0.125
       u_l   = 0.0
       u_r   = 0.0
       p_l   = 1.0
       p_r   = 0.1
    CASE(5)
       teststr = "1D Riemannn problem: Noh problem"
       j = Mesh%INUM / 2
       rho_l = 1.0
       rho_r = 1.0
       u_l   = 1.0
       u_r   = -1.0
       p_l   = 1.0E-05
       p_r   = 1.0E-05
    CASE DEFAULT
       PRINT *, "ERROR in InitData: Test problem should be 1,2,3,4, or 5"
       STOP
    END SELECT

    pvar(Mesh%IMIN-1:j,:,1)   = rho_l
    pvar(j+1:Mesh%IMAX+1,:,1) = rho_r
    pvar(Mesh%IMIN-1:j,:,2)   = u_l
    pvar(j+1:Mesh%IMAX+1,:,2) = u_r
    pvar(Mesh%IMIN-1:j,:,4)   = p_l
    pvar(j+1:Mesh%IMAX+1,:,4) = p_r

    CALL Convert2Conservative(Physics,Mesh,pvar,cvar)
    PRINT "(A,A)", " DATA-----> initial condition: ", teststr

  END SUBROUTINE InitData

END MODULE Init
