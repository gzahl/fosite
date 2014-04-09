!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_riemann2d.f90                                                #
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
! Program and data initialization for 2D Riemann problems
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
    REAL              :: test_stoptime
    INTEGER           :: geometry
    CHARACTER(LEN=256):: ofname, lfname
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Output,Logio
    !------------------------------------------------------------------------!

    ! 2D cartesian test problems  (for initial conditions see InitData)
    ! References:
    ! [1] C. W. Schulz-Rinne et al.: Numerical Solution of the Riemann Problem
    !     for Gas Dynamics, SIAM J. Sci. Comp. 14 (1993), 1394-1414
    ! [2] P. Lax, X.-D. Liu: Solution of Two-dimensional Riemann Problems of
    !     Gas Dynamics by Positive Schemes, SIAM J. Sci. Comp. 19 (1998), 319-340
    ! [3] A. Kurganov, E. Tadmor: Solution of Two-Dimensional Riemann Problems for 
    !     Gas Dynamics without Riemann Problem Solvers, NMPDE 18 (2002), 561-588
    testnum = 18

    ! set the geometry
    geometry = CARTESIAN
!    geometry = POLAR

    ! runtime of the test problem
    SELECT CASE(testnum)
    CASE(1)  ! KT test 1
       test_stoptime = 0.2
    CASE(2)  ! Riemann problem no. 2
       test_stoptime = 0.2
    CASE(3)  ! Riemann problem no. 3
       test_stoptime = 0.3
    CASE(4)  ! Riemann problem no. 4
       test_stoptime = 0.25
    CASE(5)  ! Riemann problem no. 5
       test_stoptime = 0.23
    CASE(6)  ! Riemann problem no. 6
       test_stoptime = 0.3
    CASE(12) ! Riemann problem no. 12
       test_stoptime = 0.25
    CASE(15) ! Riemann problem no. 15
       test_stoptime = 0.2
    CASE(17) ! Riemann problem no. 17
       test_stoptime = 0.3
    CASE(18) ! Riemann problem no. 18
       test_stoptime = 0.2
    CASE(19) ! Riemann problem no. 19
       test_stoptime = 0.3
    CASE DEFAULT
       PRINT *, "ERROR in InitProgram: Sorry, this 2D Riemann problem is currently not supported!"
       STOP
    END SELECT

    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER2D, &
         gamma   = 1.4, &           ! ratio of specific heats        !
         dpmax   = 1.0)             ! for advanced time step control !

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)         ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         variables = CONSERVATIVE, &! vars. to use for reconstruction!
         limiter   = MONOCENT, &    ! one of: minmod, monocent,...   !
         theta     = 1.3)           ! optional parameter for limiter !

    SELECT CASE(geometry)
    CASE(CARTESIAN)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes, &
            geometry = CARTESIAN, &
                inum = 100, &       ! resolution in x and            !
                jnum = 100, &       !   y direction                  !             
                xmin = -0.5, &
                xmax = 0.5, &
                ymin = -0.5, &
                ymax = 0.5)
       ! boundary conditions
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,NO_GRADIENTS,WEST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,NO_GRADIENTS,EAST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,NO_GRADIENTS,SOUTH)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,NO_GRADIENTS,NORTH)
    CASE(POLAR)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes, &
            geometry = POLAR, &
                inum = 50, &        ! resolution in x and            !
                jnum = 120, &       !   y direction                  !             
                xmin = 0.0, &
                xmax = 0.5, &
                ymin = 0.0, &
                ymax = 2*PI)
       ! boundary conditions
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,NO_GRADIENTS,WEST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,NO_GRADIENTS,EAST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,PERIODIC,SOUTH)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,PERIODIC,NORTH)
    CASE DEFAULT
       PRINT *, "ERROR in InitProgram: geometry should be either ", ACHAR(13), &
            "cartesian or polar"
       STOP
    END SELECT

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = test_stoptime, &
         dtlimit  = 1.0E-4, &
         maxiter  = 100000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc%pvar,Timedisc%cvar)

    ! log file name
    WRITE (lfname, '(A,I2.2,A)'), "riemann2d_", testnum, ".log"
    ! initialize log input/output
    CALL InitLogio(Logio,Mesh,Physics,Timedisc, &
         logformat = NOLOG, &
         filename  = lfname, &
         logdt     = 300)

    ! output file name
    WRITE (ofname, '(A,I2.2,A)') "riemann2d_", testnum, ".dat"
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
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: cart
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: vxy
    INTEGER           :: i,j
    REAL              :: xmin,ymin,xlength,ylength,x0,y0
    CHARACTER(LEN=64) :: teststr
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(OUT)       :: pvar,cvar
    !------------------------------------------------------------------------!

    ! convert to cartesian coordinates
    CALL Convert2Cartesian(Mesh%geometry,Mesh%bcenter,cart)

    xmin = MINVAL(cart(:,:,1))
    ymin = MINVAL(cart(:,:,2))
    xlength = MAXVAL(cart(:,:,1))-xmin
    ylength = MAXVAL(cart(:,:,2))-ymin
    x0 = xmin + 0.5*xlength
    y0 = ymin + 0.5*ylength

    ! Shock/contact discontinuity and rarefraction wave interaction;
    ! the computational domain is seperated in four quadrants
    !  -----------
    ! |  2  |  1  |
    ! |-----------|
    ! |  3  |  4  |
    !  -----------

    pvar(:,:,:) = 0.
    vxy(:,:,:) = 0.

    SELECT CASE(testnum)
    CASE(1)
       teststr = "2D Riemann problem no. 1"
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          pvar(:,:,1) = 1.
          pvar(:,:,4) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          pvar(:,:,1) = .5197
          vxy(:,:,1) = -.7259
          pvar(:,:,4) = .4
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          pvar(:,:,1) = .1072
          vxy(:,:,1) = -.7259
          vxy(:,:,2) = -1.4045
          pvar(:,:,4) = .0439
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          pvar(:,:,1) = .2579
          vxy(:,:,2) = -1.4045
          pvar(:,:,4) = .15
       END WHERE

    CASE(2)
       teststr = "2D Riemann problem no. 2" 
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          pvar(:,:,1) = 1.
          pvar(:,:,4) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          pvar(:,:,1) = .5197
          vxy(:,:,1) = -.7259
          pvar(:,:,4) = .4
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          pvar(:,:,1) = 1.
          vxy(:,:,1) = -.7259
          vxy(:,:,2) = -.7259
          pvar(:,:,4) = 1.
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          pvar(:,:,1) = .5197
          vxy(:,:,2) = -.7259
          pvar(:,:,4) = .4
       END WHERE

    CASE(3)
       teststr = "2D Riemann problem no. 3" 
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          pvar(:,:,1) = 1.5
          pvar(:,:,4) = 1.5
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          pvar(:,:,1) = .5323
          vxy(:,:,1) = 1.206
          pvar(:,:,4) = .3
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          pvar(:,:,1) = 0.138
          vxy(:,:,1) = 1.206
          vxy(:,:,2) = 1.206
          pvar(:,:,4) = 0.029
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          pvar(:,:,1) = .5323
          vxy(:,:,2) = 1.206
          pvar(:,:,4) = .3
       END WHERE

    CASE(4)
       teststr = "2D Riemann problem no. 4" 

       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          pvar(:,:,1) = 1.1
          pvar(:,:,4) = 1.1
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          pvar(:,:,1) = .5065
          vxy(:,:,1) = .8939
          pvar(:,:,4) = .35
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          pvar(:,:,1) = 1.1
          vxy(:,:,1) = .8939
          vxy(:,:,2) = .8939
          pvar(:,:,4) = 1.1
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          pvar(:,:,1) = .5065
          vxy(:,:,2) = .8939
          pvar(:,:,4) = .35
       END WHERE

    CASE(5)
       teststr = "2D Riemann problem no. 5" 
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          pvar(:,:,1) = 1.
          vxy(:,:,1) = -.75
          vxy(:,:,2) = -.5
          pvar(:,:,4) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          pvar(:,:,1) = 2.
          vxy(:,:,1) = -.75
          vxy(:,:,2) = .5
          pvar(:,:,4) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          pvar(:,:,1) = 1.
          vxy(:,:,1) = .75
          vxy(:,:,2) = .5
          pvar(:,:,4) = 1.
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          pvar(:,:,1) = 3.
          vxy(:,:,1) = .75
          vxy(:,:,2) = -.5
          pvar(:,:,4) = 1.
       END WHERE

    CASE(6)
       teststr = "2D Riemann problem no. 6" 
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          pvar(:,:,1) = 1.
          vxy(:,:,1) = 0.75
          vxy(:,:,2) = -0.5
          pvar(:,:,4) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          pvar(:,:,1) = 2.
          vxy(:,:,1) = 0.75
          vxy(:,:,2) = 0.5
          pvar(:,:,4) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          pvar(:,:,1) = 1.
          vxy(:,:,1) = -0.75
          vxy(:,:,2) = 0.5
          pvar(:,:,4) = 1.
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          pvar(:,:,1) = 3.
          vxy(:,:,1) = -0.75
          vxy(:,:,2) = -0.5
          pvar(:,:,4) = 1.
       END WHERE
       
    CASE(12)
       teststr = "2D Riemann problem no. 12" 
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          pvar(:,:,1) = 0.5313
          vxy(:,:,1) = 0.
          vxy(:,:,2) = 0.
          pvar(:,:,4) = 0.4
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          pvar(:,:,1) = 1.
          vxy(:,:,1) = 0.7276
          vxy(:,:,2) = 0.
          pvar(:,:,4) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          pvar(:,:,1) = 0.8
          vxy(:,:,1) = 0.
          vxy(:,:,2) = 0.
          pvar(:,:,4) = 1.
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          pvar(:,:,1) = 1.
          vxy(:,:,1) = 0.
          vxy(:,:,2) = 0.7276
          pvar(:,:,4) = 1.
       END WHERE
       
    CASE(15)
       teststr = "2D Riemann problem no. 15"
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          pvar(:,:,1) = 1.
          vxy(:,:,1) = 0.1
          vxy(:,:,2) = -0.3
          pvar(:,:,4) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          pvar(:,:,1) = 0.5197
          vxy(:,:,1) = -0.6259
          vxy(:,:,2) = -0.3
          pvar(:,:,4) = 0.4
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          pvar(:,:,1) = 0.8
          vxy(:,:,1) = 0.1
          vxy(:,:,2) = -0.3
          pvar(:,:,4) = 0.4
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          pvar(:,:,1) = 0.5313
          vxy(:,:,1) = 0.1
          vxy(:,:,2) = 0.4276
          pvar(:,:,4) = 0.4
       END WHERE
       
    CASE(17)
       teststr = "2D Riemann problem no. 17"
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          pvar(:,:,1) = 1.
          vxy(:,:,1) = 0.
          vxy(:,:,2) = -0.4
          pvar(:,:,4) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          pvar(:,:,1) = 2.
          vxy(:,:,1) = 0.
          vxy(:,:,2) = -0.3
          pvar(:,:,4) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          pvar(:,:,1) = 1.0625
          vxy(:,:,1) = 0.
          vxy(:,:,2) = 0.2145
          pvar(:,:,4) = 0.4
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          pvar(:,:,1) = 0.5197
          vxy(:,:,1) = 0.
          vxy(:,:,2) = -1.1259
          pvar(:,:,4) = 0.4
       END WHERE
       
    CASE(18)
       teststr = "2D Riemann problem no. 18"
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          pvar(:,:,1) = 1.
          vxy(:,:,1) = 0.
          vxy(:,:,2) = 1.
          pvar(:,:,4) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          pvar(:,:,1) = 2.
          vxy(:,:,1) = 0.
          vxy(:,:,2) = -0.3
          pvar(:,:,4) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          pvar(:,:,1) = 1.0625
          vxy(:,:,1) = 0.
          vxy(:,:,2) = 0.2145
          pvar(:,:,4) = 0.4
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          pvar(:,:,1) = 0.5197
          vxy(:,:,1) = 0.
          vxy(:,:,2) = 0.2741
          pvar(:,:,4) = 0.4
       END WHERE
       
    CASE(19)
       teststr = "2D Riemann problem no. 19"
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          pvar(:,:,1) = 1.
          vxy(:,:,1) = 0.
          vxy(:,:,2) = 0.3
          pvar(:,:,4) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          pvar(:,:,1) = 2.
          vxy(:,:,1) = 0.
          vxy(:,:,2) = -0.3
          pvar(:,:,4) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          pvar(:,:,1) = 1.0625
          vxy(:,:,1) = 0.
          vxy(:,:,2) = 0.2145
          pvar(:,:,4) = 0.4
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          pvar(:,:,1) = 0.5197
          vxy(:,:,1) = 0.
          vxy(:,:,2) = -0.4259
          pvar(:,:,4) = 0.4
       END WHERE
       
    CASE DEFAULT
       PRINT *, "ERROR in InitData: Sorry, this 2D Riemann problem is currently not supported!"
       STOP
    END SELECT
    
    CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,vxy,pvar(:,:,2:3))

    CALL Convert2Conservative(Physics,Mesh,pvar,cvar)
    PRINT "(A,A)", " DATA-----> initial condition: ", teststr

  END SUBROUTINE InitData

END MODULE Init
