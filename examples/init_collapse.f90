!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_collapse.f90                                                 #
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
! Program and data initialization for collapse with angular momentum
!----------------------------------------------------------------------------!

!**************************************!
!* IMPORTANT:                         *!
!* - compile with autodouble          *!
!**************************************!
 
MODULE Init
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE output_generic
  USE logio_generic
  USE sources_generic
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
    INTEGER           :: geometry
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Output,Logio
    !------------------------------------------------------------------------!

    ! set the geometry
    geometry = CYLINDRICAL
!    geometry = SPHERICAL
!    geometry = OBLATE_SPHEROIDAL

    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER3D_ROTSYM, &
         gamma   = 1.4, &           ! ratio of specific heats        !
         mu      = 0.602E-03, &     ! mean molecular weight          !
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
                inum = 100, &       ! resolution in x and            !
                jnum = 50, &        !   y direction                  !             
                xmin = -2.0E+14, &
                xmax = 2.0E+14, &
                ymin = 0.0, &
                ymax = 2.0E+14)
       ! boundary conditions
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,WEST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,EAST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,AXIS,SOUTH)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,NORTH)
    CASE(SPHERICAL)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes,SPHERICAL,50,30,0.0,2.0E+14,0.0,PI)
       ! boundary conditions
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,WEST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,EAST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,AXIS,SOUTH)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,AXIS,NORTH)
    CASE(OBLATE_SPHEROIDAL)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes,OBLATE_SPHEROIDAL,50,60, &
            0.0,1.4,-0.5*PI,0.5*PI, &
            1.0E+14)                ! optional geometry parameter    !
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

    ! source term due to a point mass
    CALL InitSources(Physics%sources,Mesh,Fluxes,Physics, &
         stype  = POINTMASS, &       ! grav. accel. of a point mass   !
         sparam = 1.0E+30)           ! mass [kg]                      !

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = 1.0E+11, &
         dtlimit  = 1.0, &
         maxiter  = 1000000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc%pvar,Timedisc%cvar)

    ! initialize log input/output
    CALL InitLogio(Logio,Mesh,Physics,Timedisc, &
         logformat = NOLOG, &
         filename  = "collapse.log", &
         logdt     = 300)

    ! set parameters for data output
    CALL InitOutput(Output,Mesh,Physics,Timedisc,&
         filetype  = GNUPLOT, &
         filename  = "collapse.dat", &
         mode      = OVERWRITE, &
         starttime = 0.0, &
         stoptime  = Timedisc%stoptime, &
         count     = 10)

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
    REAL              :: rmax,omega
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(OUT)       :: pvar,cvar
    !------------------------------------------------------------------------!

    ! radius of the sphere
    rmax = 1.0D+14

    ! angular velocity of the sphere
    omega = 1.0E-13

    CALL Convert2Cartesian(Mesh%geometry,Mesh%bcenter,cart)

    ! homogenious sphere with constant angular velocity
    WHERE (SQRT(cart(:,:,1)**2+cart(:,:,2)**2).LE.rmax)
       pvar(:,:,1) = 1.0E-16
    ELSEWHERE
       pvar(:,:,1) = 1.0E-20
    END WHERE

    ! velocities
    pvar(:,:,2:3) = 0.
    ! constant angular velocity
    pvar(:,:,4) = omega * cart(:,:,1)

    IF (GetType(Physics).EQ.EULER3D_ROTAMT) THEN
       ! specific angular momentum
       pvar(:,:,4) = pvar(:,:,4)*Mesh%bhz(:,:)
    END IF

    ! constant temperature 10 K
    pvar(:,:,5) = Physics%constants%RG/Physics%mu * 10.0 * pvar(:,:,1)

    CALL Convert2Conservative(Physics,Mesh,pvar,cvar)

    PRINT "(A,A)", " DATA-----> initial condition: ", &
         "Homogenious density w/ uniform angular velocity"
  END SUBROUTINE InitData

END MODULE Init
