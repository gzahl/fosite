!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_riemann3d.f90                                                #
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
! Program and data initialization for 3D Riemann problems
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
    INTEGER           :: onum
    INTEGER           :: geometry
    REAL              :: test_stoptime
    CHARACTER(LEN=256):: ofname, lfname
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Output,Logio
    !------------------------------------------------------------------------!


    ! 3D spherical Riemann problem between walls
    ! 1. setup according to J. O. Langseth and R. J. LeVeque. "A wave propagation 
    !    "method for three-dimensional hyperbolic conservation laws", J. o. Comp.
    !    Phys., 165(1), 2000;
    ! 2. setup with angular velocity
    testnum=1

    ! set the geometry
    geometry = CYLINDRICAL
!    geometry = SPHERICAL
!    geometry = OBLATE_SPHEROIDAL

    SELECT CASE(testnum)
    CASE(1) ! original Test of Langseth and LeVeque
       test_stoptime = 0.7
       onum=7                              ! number of output data sets
       ofname="leveque.dat"
       lfname="leveque.log"
    CASE(2) ! angular momentum test
       test_stoptime = 1.5
       onum=15                             ! number of output data sets
       ofname="rotspher.dat"
       lfname="rotspher.log"
     CASE DEFAULT
       PRINT *, "ERROR in InitProgram: testnum should be one of 1, 2"
       STOP
    END SELECT

    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER3D_ROTSYM, &
         gamma   = 1.4, &           ! ratio of specific heats        !
         mu      = 0.029, &         ! mean molecular weight          !
         dpmax   = 10.0)            ! for advanced time step control !

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)         ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         variables = PRIMITIVE, &   ! vars. to use for reconstruction!
         limiter   = MONOCENT, &    ! one of: minmod, monocent,...   !
         theta     = 1.3)           ! optional parameter for limiter !

    SELECT CASE(geometry)
    CASE(CYLINDRICAL)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes, &
            geometry = CYLINDRICAL, &
                inum = 80,  &       ! resolution in x and            !
                jnum = 120, &       !   y direction                  !             
                xmin = 0.0, &
                xmax = 1.0, &
                ymin = 0.0, &
                ymax = 1.5)
       ! boundary conditions
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,WEST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,EAST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,AXIS,SOUTH)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,NORTH)
    CASE(SPHERICAL)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes,SPHERICAL,50,60,0.0,1.5,0.0,0.5*PI)
       ! boundary conditions
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,WEST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,EAST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,AXIS,SOUTH)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,AXIS,NORTH)
    CASE(OBLATE_SPHEROIDAL)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes,OBLATE_SPHEROIDAL,50,60, &
            0.0,1.4,-0.5*PI,0.5*PI, &
            0.75)                   ! optional geometry parameter    !
       ! boundary conditions
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,WEST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,EAST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,AXIS,SOUTH)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,AXIS,NORTH)
    CASE DEFAULT
       PRINT *, "ERROR in InitProgram: geometry should be either ", ACHAR(13), &
            "cylindrical, spherical or oblate spheroidal"
       STOP
    END SELECT

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = test_stoptime, &
         dtlimit  = 1.0E-9, &
         maxiter  = 100000)

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
         count     = onum)

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
    INTEGER, PARAMETER :: icmax=1, jcmax=1
    REAL, DIMENSION(-icmax:icmax,-jcmax:jcmax,2) :: curv, cart
    REAL              :: dcx, dcy
    REAL              :: rmax
    REAL              :: xpos,ypos
    REAL              :: rho_in, rho_out, omega_in, omega_out, P_in, P_out
    REAL              :: ratio
    INTEGER           :: i,j,ic,jc,nc,ncmax
    CHARACTER(LEN=80) :: teststr
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(OUT)       :: pvar,cvar
    !------------------------------------------------------------------------!

    SELECT CASE(testnum)
    CASE(1) ! original Test of Langseth and LeVeque
       teststr = "Spherical pressure discontinuity between walls"
       ! radius and position
       rmax = 0.2
       xpos = 0.0
       ypos = 0.4
       ! density and pressure
       rho_in  = 1.0
       rho_out = 1.0
       P_in    = 5.0
       P_out   = 1.0
       ! angular velocity
       omega_in  = 0.0
       omega_out = 0.0
    CASE(2) ! angular momentum test
       teststr = "Homogenious sphere w/ uniform angular velocity"
       ! radius and position
       rmax = 0.2
       xpos = 0.0
       ypos = 0.4
       ! density and pressure
       rho_in  = 1.0
       rho_out = 1.0E-04
       P_in    = 1.0
       P_out   = 1.0
       ! angular velocity
       omega_in  = 5.
       omega_out = 5.
    CASE DEFAULT
       PRINT *, "ERROR in InitProgram: testnum should be one of 1, 2"
       STOP
    END SELECT

    ! set density and pressure for homogenious sphere
    dcx = Mesh%dx / (2*icmax+1)
    dcy = Mesh%dy / (2*jcmax+1)
    ncmax = (2*icmax+1)*(2*jcmax+1)
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! use subgrid for initial data 
          ! to obtain better resolution
          DO jc=-jcmax,jcmax
             DO ic=-icmax,icmax
                curv(ic,jc,1) = Mesh%center(i,j,1) + ic*dcx
                curv(ic,jc,2) = Mesh%center(i,j,2) + jc*dcy
             END DO
          END DO
          CALL Convert2Cartesian(Mesh%geometry,curv,cart)
          nc = 0
          DO jc=-jcmax,jcmax
             DO ic=-icmax,icmax
                IF (SQRT((xpos-cart(ic,jc,1))**2+(ypos-cart(ic,jc,2))**2).LE.rmax) THEN
                   nc = nc + 1
                END IF
             END DO
          END DO
          ratio = (1.0 * nc) / ncmax
          ! density
          pvar(i,j,1) = ratio * rho_in + (1.-ratio) * rho_out
          ! angular velocity
          pvar(i,j,4) = (ratio * omega_in + (1.-ratio) * omega_out) * Mesh%bhz(i,j)
          ! pressure
          pvar(i,j,5) = ratio * P_in + (1.-ratio) * P_out
       END DO
    END DO

    ! velocities
    pvar(:,:,2:3) = 0.

    ! for specific angular momentum transport
    IF (GetType(Physics).EQ.EULER3D_ROTAMT) THEN
       ! specific angular momentum
       pvar(:,:,4) = pvar(:,:,4)*Mesh%bhz(:,:)
    END IF

    CALL Convert2Conservative(Physics,Mesh,pvar,cvar)

    PRINT "(A,A)", " DATA-----> initial condition: ", TRIM(teststr)
  END SUBROUTINE InitData

END MODULE Init
