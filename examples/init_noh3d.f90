!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_noh3d.f90                                                    #
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
! Program and data initialization for 3D Noh problem
!----------------------------------------------------------------------------!

!**************************************!
!* IMPORTANT:                         *!
!* - use primitive reconstruction     *!
!**************************************!

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
         gamma   = 5./3., &         ! ratio of specific heats        !
         dpmax   = 1000.0)          ! for advanced time step control !

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
                jnum = 100, &       !   y direction                  !             
                xmin = 0.0, &
                xmax = 0.5, &
                ymin = 0.0, &
                ymax = 0.5)
                                    ! of output data sets            !
       ! boundary conditions
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! WARNING WARNING WARNING !!!
       !!! bad boundary conditions !!!
       !!! WARNING WARNING WARNING !!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,WEST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,FIXED,EAST)
       CALL InitBoundaryData(Mesh%boundary,Mesh,Fluxes,Physics,EAST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,SOUTH)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,FIXED,NORTH)
       CALL InitBoundaryData(Mesh%boundary,Mesh,Fluxes,Physics,NORTH)
    CASE(SPHERICAL)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes,SPHERICAL,50,30,0.0,0.5,0.0,0.5*PI)
       ! boundary conditions
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! WARNING WARNING WARNING !!!
       !!! bad boundary conditions !!!
       !!! WARNING WARNING WARNING !!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,WEST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,FIXED,EAST)
       CALL InitBoundaryData(Mesh%boundary,Mesh,Fluxes,Physics,EAST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,AXIS,SOUTH)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,NORTH)
    CASE(OBLATE_SPHEROIDAL)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! WARNING WARNING WARNING !!!
       !!!        not tested       !!!
       !!! WARNING WARNING WARNING !!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! set parameters for data output
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes,OBLATE_SPHEROIDAL,50,30, &
            0.0,0.48,0.0,0.5*PI, &
            1.0)                    ! optional geometry parameter    !
       ! boundary conditions
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,WEST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,NO_GRADIENTS,EAST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,AXIS,SOUTH)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,NORTH)
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
         stoptime = 0.3, &
         dtlimit  = 1.0E-6, &
         maxiter  = 100000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc%pvar,Timedisc%cvar)

    ! initialize log input/output
    CALL InitLogio(Logio,Mesh,Physics,Timedisc, &
         logformat = NOLOG, &
         filename  = "noh3d.log", &
         logdt     = 300)

    ! set parameters for data output
    CALL InitOutput(Output,Mesh,Physics,Timedisc,&
         filetype  = GNUPLOT, &
         filename  = "noh3d.dat", &
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
    INTEGER           :: i,j
    REAL              :: rho0, vr0, P0
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: cart
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: vcart
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(OUT)       :: pvar,cvar
    !------------------------------------------------------------------------!

    CALL Convert2Cartesian(Mesh%geometry,Mesh%bcenter,cart)
    
    rho0 = 1.0
    vr0  = -1.0
    P0   = 1.0E-06

    pvar(:,:,1) = rho0
    pvar(:,:,4) = 0.
    pvar(:,:,5) = P0

    FORALL (i=Mesh%IGMIN:Mesh%IGMAX,j=Mesh%JGMIN:Mesh%JGMAX)
       vcart(i,j,1) = vr0 * cart(i,j,1) / SQRT(cart(i,j,1)**2 + cart(i,j,2)**2)
       vcart(i,j,2) = vcart(i,j,1) * cart(i,j,2) / cart(i,j,1)
    END FORALL

    CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,vcart,pvar(:,:,2:3))

    CALL Convert2Conservative(Physics,Mesh,pvar,cvar)
    PRINT "(A,A)", " DATA-----> initial condition: ", &
         "3D Noh problem"

  END SUBROUTINE InitData


  SUBROUTINE InitBoundaryData(Boundary,Mesh,Fluxes,Physics,direction)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), DIMENSION(4) :: Boundary
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    INTEGER           :: direction
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j
    REAL              :: rho0, vr0, P0, gm1
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: cart, vcart, vcurv
    REAL, DIMENSION(2,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) :: pdata_east
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,2,Physics%vnum) :: pdata_north
    !------------------------------------------------------------------------!

    rho0 = 1.0
    vr0  = -1.0
    P0   = 1.0E-06
    gm1  = Physics%gamma - 1.

    CALL Convert2Cartesian(Mesh%geometry,Mesh%bcenter,cart)

    ! supersonic inflow -> set all variables
    Boundary(direction)%fixed = (/ .TRUE., .TRUE., .TRUE., .TRUE., .TRUE. /)
    SELECT CASE(GetType(Mesh%geometry))
    CASE(CYLINDRICAL)
       SELECT CASE(direction)
       CASE(EAST)
          pdata_east(:,:,1) = rho0
          pdata_east(:,:,4) = 0.
          pdata_east(:,:,5) = P0
          FORALL (i=1:2,j=Mesh%JGMIN:Mesh%JGMAX)
             pdata_east(i,j,2) = vr0 * cart(Mesh%IMAX+i,j,2) &
                  / SQRT(cart(Mesh%IMAX+i,j,1)**2 + cart(Mesh%IMAX+i,j,2)**2)
             pdata_east(i,j,3) = pdata_east(i,j,2) * cart(Mesh%IMAX+i,j,1) &
                  / cart(Mesh%IMAX+i,j,2)
          END FORALL
       CASE(NORTH)
          pdata_north(:,:,1) = rho0
          pdata_north(:,:,4) = 0.
          pdata_north(:,:,5) = P0
          FORALL (i=Mesh%IGMIN:Mesh%IGMAX,j=1:2)
             pdata_north(i,j,2) = vr0 * cart(i,Mesh%JMAX+j,2) &
                  / SQRT(cart(i,Mesh%JMAX+j,1)**2 + cart(i,Mesh%JMAX+j,2)**2)
             pdata_north(i,j,3) = pdata_north(i,j,2)* cart(i,Mesh%JMAX+j,1) &
                  / cart(i,Mesh%JMAX+j,2)
          END FORALL
       CASE DEFAULT
          PRINT *, "ERROR in InitBoundaryData: direction should be either EAST or NORTH"
          STOP
       END SELECT
    CASE(SPHERICAL)
       pdata_east(:,:,1) = rho0
       pdata_east(:,:,2) = vr0
       pdata_east(:,:,3) = 0.
       pdata_east(:,:,4) = 0.
       pdata_east(:,:,5) = P0
    CASE DEFAULT
       PRINT *, "ERROR in InitBoundaryData: geometry should be either ", ACHAR(13), &
            "cylindrical, spherical or oblate spheroidal"
       STOP
    END SELECT

    ! set boundary data to either primitive or conservative values
    ! depending on the reconstruction
    SELECT CASE(direction)
    CASE(WEST,EAST)
       IF (PrimRecon(Fluxes%reconstruction).EQV.PRIMITIVE) THEN
          Boundary(direction)%data = pdata_east
       ELSE
          Boundary(direction)%data(:,:,1) = pdata_east(:,:,1)
          Boundary(direction)%data(:,:,2) = pdata_east(:,:,1)*pdata_east(:,:,2)
          Boundary(direction)%data(:,:,3) = pdata_east(:,:,1)*pdata_east(:,:,3)
          Boundary(direction)%data(:,:,4) = pdata_east(:,:,1)*pdata_east(:,:,4)
          Boundary(direction)%data(:,:,5) = pdata_east(:,:,5)/gm1 &
               + 0.5*pdata_east(:,:,1)*pdata_east(:,:,2)**2
       END IF
    CASE(SOUTH,NORTH)
       IF (PrimRecon(Fluxes%reconstruction).EQV.PRIMITIVE) THEN
          Boundary(direction)%data = pdata_north
       ELSE
          Boundary(direction)%data(:,:,1) = pdata_north(:,:,1)
          Boundary(direction)%data(:,:,2) = pdata_north(:,:,1)*pdata_north(:,:,2)
          Boundary(direction)%data(:,:,3) = pdata_north(:,:,1)*pdata_north(:,:,3)
          Boundary(direction)%data(:,:,4) = pdata_north(:,:,1)*pdata_north(:,:,4)
          Boundary(direction)%data(:,:,5) = pdata_north(:,:,5)/gm1 &
               + 0.5*pdata_north(:,:,1)*pdata_north(:,:,2)**2
       END IF
    CASE DEFAULT
       PRINT *, "ERROR in InitBoundaryData: direction should be one of WEST,EAST,SOUTH,NORTH."
       STOP
    END SELECT

  END SUBROUTINE InitBoundaryData

END MODULE Init
