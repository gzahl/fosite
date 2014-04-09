!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_riemann2d.f90                                                #
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
! Program and data initialization for 2D Riemann problems
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
    REAL              :: test_stoptime
    INTEGER           :: geometry
    CHARACTER(LEN=256):: ofname, lfname
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!

    ! 2D cartesian test problems  (for initial conditions see InitData)
    ! References:
    ! [1] C. W. Schulz-Rinne et al.: Numerical Solution of the Riemann Problem
    !     for Gas Dynamics, SIAM J. Sci. Comp. 14 (1993), 1394-1414
    ! [2] P. Lax, X.-D. Liu: Solution of Two-dimensional Riemann Problems of
    !     Gas Dynamics by Positive Schemes, SIAM J. Sci. Comp. 19 (1998), 319-340
    ! [3] A. Kurganov, E. Tadmor: Solution of Two-Dimensional Riemann Problems for 
    !     Gas Dynamics without Riemann Problem Solvers, NMPDE 18 (2002), 561-588
    testnum = 1

    ! set the geometry
    geometry = CARTESIAN
!    geometry = POLAR

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
         theta     = 1.2)           ! optional parameter for limiter !

    ! mesh settings
    SELECT CASE(geometry)
    CASE(CARTESIAN)
       CALL InitMesh(Mesh,Fluxes, &
            geometry = CARTESIAN, &
                inum = 100, &       ! resolution in x and            !
                jnum = 100, &       !   y direction                  !             
                xmin = -0.5, &
                xmax = 0.5, &
                ymin = -0.5, &
                ymax = 0.5)
    CASE(POLAR)
       CALL InitMesh(Mesh,Fluxes, &
            geometry = POLAR, &
                inum = 50, &        ! resolution in x and            !
                jnum = 120, &       !   y direction                  !             
                xmin = 0.0, &
                xmax = 0.5, &
                ymin = 0.0, &
                ymax = 2*PI)
    CASE DEFAULT
       CALL Error(Physics,"InitProgram", &
            " geometry should be either cartesian or polar")
    END SELECT

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
       CALL Error(Physics,"InitProgram", &
            "Sorry, this 2D Riemann problem is currently not supported!")
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
    CALL InitData(Mesh,Physics,Timedisc)

    ! boundary conditions
    SELECT CASE(geometry)
    CASE(CARTESIAN)
       CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
            western  = NO_GRADIENTS, &
            eastern  = NO_GRADIENTS, &
            southern = NO_GRADIENTS, &
            northern = NO_GRADIENTS)
    CASE(POLAR)
       CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
            western  = NO_GRADIENTS, &
            eastern  = NO_GRADIENTS, &
            southern = PERIODIC, &
            northern = PERIODIC)
    CASE DEFAULT
       PRINT *, "ERROR in InitProgram: geometry should be either ", ACHAR(13), &
            "cartesian or polar"
       STOP
    END SELECT

    ! initialize log input/output
    WRITE (lfname, '(A,I2.2,A)') "riemann2dlog_", testnum
    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,&
         fileformat = BINARY, &
#ifdef PARALLEL
         filename   = "/tmp/" // lfname, &
#else
         filename   = lfname, &
#endif
         dtwall     = 1800, &
         filecycles = 1)

    ! initialize data input/output
    WRITE (ofname, '(A,I2.2,A)') "riemann2d_", testnum
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc,&
         fileformat = BINARY, &
#ifdef PARALLEL
         filename   = "/tmp/" // ofname, &
#else
         filename   = ofname, &
#endif
         filecycles = 2, &
         count      = 1)
  END SUBROUTINE InitProgram


  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
#ifdef PARALLEL
    include 'mpif.h'
#endif
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: cart
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: vxy
    REAL              :: xmin,ymin,xmax,ymax,x0,y0
#ifdef PARALLEL
    REAL              :: xmin_all,xmax_all,ymin_all,ymax_all
    INTEGER           :: ierr
#endif
    CHARACTER(LEN=64) :: teststr
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!

    ! convert to cartesian coordinates
    CALL Convert2Cartesian(Mesh%geometry,Mesh%bcenter,cart)

    ! minima and maxima of _cartesian_ coordinates
    xmin = MINVAL(cart(:,:,1))
    ymin = MINVAL(cart(:,:,2))
    xmax = MAXVAL(cart(:,:,1))
    ymax = MAXVAL(cart(:,:,2))
#ifdef PARALLEL
    CALL MPI_Allreduce(xmin,xmin_all,1,DEFAULT_MPI_REAL,MPI_MIN,Mesh%comm_cart,ierr)
    xmin = xmin_all
    CALL MPI_Allreduce(ymin,ymin_all,1,DEFAULT_MPI_REAL,MPI_MIN,Mesh%comm_cart,ierr)
    ymin = ymin_all
    CALL MPI_Allreduce(xmax,xmax_all,1,DEFAULT_MPI_REAL,MPI_MAX,Mesh%comm_cart,ierr)
    xmax = xmax_all
    CALL MPI_Allreduce(ymax,ymax_all,1,DEFAULT_MPI_REAL,MPI_MAX,Mesh%comm_cart,ierr)
    ymax = ymax_all
#endif
    x0 = xmin + 0.5*ABS(xmax-xmin)
    y0 = ymin + 0.5*ABS(ymax-ymin)

    ! Shock/contact discontinuity and rarefraction wave interaction;
    ! the computational domain is seperated in four quadrants
    !  -----------
    ! |  2  |  1  |
    ! |-----------|
    ! |  3  |  4  |
    !  -----------

    Timedisc%pvar(:,:,:) = 0.
    vxy(:,:,:) = 0.

    SELECT CASE(testnum)
    CASE(1)
       teststr = "2D Riemann problem no. 1"
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          Timedisc%pvar(:,:,Physics%DENSITY) = .5197
          vxy(:,:,1) = -.7259
          Timedisc%pvar(:,:,Physics%PRESSURE) = .4
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar(:,:,Physics%DENSITY) = .1072
          vxy(:,:,1) = -.7259
          vxy(:,:,2) = -1.4045
          Timedisc%pvar(:,:,Physics%PRESSURE) = .0439
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar(:,:,Physics%DENSITY) = .2579
          vxy(:,:,2) = -1.4045
          Timedisc%pvar(:,:,Physics%PRESSURE) = .15
       END WHERE

    CASE(2)
       teststr = "2D Riemann problem no. 2" 
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          Timedisc%pvar(:,:,Physics%DENSITY) = .5197
          vxy(:,:,1) = -.7259
          Timedisc%pvar(:,:,Physics%PRESSURE) = .4
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.
          vxy(:,:,1) = -.7259
          vxy(:,:,2) = -.7259
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar(:,:,Physics%DENSITY) = .5197
          vxy(:,:,2) = -.7259
          Timedisc%pvar(:,:,Physics%PRESSURE) = .4
       END WHERE

    CASE(3)
       teststr = "2D Riemann problem no. 3" 
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.5
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.5
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          Timedisc%pvar(:,:,Physics%DENSITY) = .5323
          vxy(:,:,1) = 1.206
          Timedisc%pvar(:,:,Physics%PRESSURE) = .3
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar(:,:,Physics%DENSITY) = 0.138
          vxy(:,:,1) = 1.206
          vxy(:,:,2) = 1.206
          Timedisc%pvar(:,:,Physics%PRESSURE) = 0.029
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar(:,:,Physics%DENSITY) = .5323
          vxy(:,:,2) = 1.206
          Timedisc%pvar(:,:,Physics%PRESSURE) = .3
       END WHERE

    CASE(4)
       teststr = "2D Riemann problem no. 4" 

       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.1
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.1
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          Timedisc%pvar(:,:,Physics%DENSITY) = .5065
          vxy(:,:,1) = .8939
          Timedisc%pvar(:,:,Physics%PRESSURE) = .35
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.1
          vxy(:,:,1) = .8939
          vxy(:,:,2) = .8939
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.1
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar(:,:,Physics%DENSITY) = .5065
          vxy(:,:,2) = .8939
          Timedisc%pvar(:,:,Physics%PRESSURE) = .35
       END WHERE

    CASE(5)
       teststr = "2D Riemann problem no. 5" 
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.
          vxy(:,:,1) = -.75
          vxy(:,:,2) = -.5
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          Timedisc%pvar(:,:,Physics%DENSITY) = 2.
          vxy(:,:,1) = -.75
          vxy(:,:,2) = .5
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.
          vxy(:,:,1) = .75
          vxy(:,:,2) = .5
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar(:,:,Physics%DENSITY) = 3.
          vxy(:,:,1) = .75
          vxy(:,:,2) = -.5
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       END WHERE

    CASE(6)
       teststr = "2D Riemann problem no. 6" 
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.
          vxy(:,:,1) = 0.75
          vxy(:,:,2) = -0.5
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          Timedisc%pvar(:,:,Physics%DENSITY) = 2.
          vxy(:,:,1) = 0.75
          vxy(:,:,2) = 0.5
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.
          vxy(:,:,1) = -0.75
          vxy(:,:,2) = 0.5
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar(:,:,Physics%DENSITY) = 3.
          vxy(:,:,1) = -0.75
          vxy(:,:,2) = -0.5
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       END WHERE
       
    CASE(12)
       teststr = "2D Riemann problem no. 12" 
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          Timedisc%pvar(:,:,Physics%DENSITY) = 0.5313
          vxy(:,:,1) = 0.
          vxy(:,:,2) = 0.
          Timedisc%pvar(:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.
          vxy(:,:,1) = 0.7276
          vxy(:,:,2) = 0.
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar(:,:,Physics%DENSITY) = 0.8
          vxy(:,:,1) = 0.
          vxy(:,:,2) = 0.
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.
          vxy(:,:,1) = 0.
          vxy(:,:,2) = 0.7276
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       END WHERE
       
    CASE(15)
       teststr = "2D Riemann problem no. 15"
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.
          vxy(:,:,1) = 0.1
          vxy(:,:,2) = -0.3
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          Timedisc%pvar(:,:,Physics%DENSITY) = 0.5197
          vxy(:,:,1) = -0.6259
          vxy(:,:,2) = -0.3
          Timedisc%pvar(:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar(:,:,Physics%DENSITY) = 0.8
          vxy(:,:,1) = 0.1
          vxy(:,:,2) = -0.3
          Timedisc%pvar(:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar(:,:,Physics%DENSITY) = 0.5313
          vxy(:,:,1) = 0.1
          vxy(:,:,2) = 0.4276
          Timedisc%pvar(:,:,Physics%PRESSURE) = 0.4
       END WHERE
       
    CASE(17)
       teststr = "2D Riemann problem no. 17"
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.
          vxy(:,:,1) = 0.
          vxy(:,:,2) = -0.4
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          Timedisc%pvar(:,:,Physics%DENSITY) = 2.
          vxy(:,:,1) = 0.
          vxy(:,:,2) = -0.3
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.0625
          vxy(:,:,1) = 0.
          vxy(:,:,2) = 0.2145
          Timedisc%pvar(:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar(:,:,Physics%DENSITY) = 0.5197
          vxy(:,:,1) = 0.
          vxy(:,:,2) = -1.1259
          Timedisc%pvar(:,:,Physics%PRESSURE) = 0.4
       END WHERE
       
    CASE(18)
       teststr = "2D Riemann problem no. 18"
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.
          vxy(:,:,1) = 0.
          vxy(:,:,2) = 1.
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          Timedisc%pvar(:,:,Physics%DENSITY) = 2.
          vxy(:,:,1) = 0.
          vxy(:,:,2) = -0.3
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.0625
          vxy(:,:,1) = 0.
          vxy(:,:,2) = 0.2145
          Timedisc%pvar(:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar(:,:,Physics%DENSITY) = 0.5197
          vxy(:,:,1) = 0.
          vxy(:,:,2) = 0.2741
          Timedisc%pvar(:,:,Physics%PRESSURE) = 0.4
       END WHERE
       
    CASE(19)
       teststr = "2D Riemann problem no. 19"
       WHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 1
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.
          vxy(:,:,1) = 0.
          vxy(:,:,2) = 0.3
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).GT.y0) )
          ! no. 2
          Timedisc%pvar(:,:,Physics%DENSITY) = 2.
          vxy(:,:,1) = 0.
          vxy(:,:,2) = -0.3
          Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (cart(:,:,1).LT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar(:,:,Physics%DENSITY) = 1.0625
          vxy(:,:,1) = 0.
          vxy(:,:,2) = 0.2145
          Timedisc%pvar(:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (cart(:,:,1).GT.x0).AND.(cart(:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar(:,:,Physics%DENSITY) = 0.5197
          vxy(:,:,1) = 0.
          vxy(:,:,2) = -0.4259
          Timedisc%pvar(:,:,Physics%PRESSURE) = 0.4
       END WHERE
       
    CASE DEFAULT
       CALL Error(Mesh,"InitData", &
            "Sorry, this 2D Riemann problem is currently not supported!")
    END SELECT
    
    CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,vxy,&
         Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY))

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: " // TRIM(teststr))

  END SUBROUTINE InitData

END MODULE Init