!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_sedov3d.f90                                                  #
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
! Program and data initialization for 3D Sedov explosion
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
         dpmax   = 100.0)           ! for advanced time step control !

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)         ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         variables = CONSERVATIVE,& ! vars. to use for reconstruction!
         limiter   = MONOCENT, &    ! one of: minmod, monocent,...   !
         theta     = 1.3)           ! optional parameter for limiter !

    SELECT CASE(geometry)
    CASE(CYLINDRICAL)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes, &
            geometry = CYLINDRICAL, &
                inum = 100, &       ! resolution in x and            !
                jnum = 50, &        !   y direction                  !             
                xmin = -0.4, &
                xmax = 0.4, &
                ymin = 0.0, &
                ymax = 0.4)
       ! boundary conditions
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,NO_GRADIENTS,WEST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,NO_GRADIENTS,EAST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,AXIS,SOUTH)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,NO_GRADIENTS,NORTH)
    CASE(SPHERICAL)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes,SPHERICAL,50,30,0.0,0.4,0.0,PI)
       ! boundary conditions
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,WEST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,NO_GRADIENTS,EAST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,AXIS,SOUTH)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,AXIS,NORTH)
    CASE(OBLATE_SPHEROIDAL)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes,OBLATE_SPHEROIDAL,50,60, &
            0.0,1.4,-0.5*PI,0.5*PI, &
            0.2)                    ! optional geometry parameter    !
       ! boundary conditions
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,FOLDED,WEST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,NO_GRADIENTS,EAST)
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
         stoptime = 0.05, &
         dtlimit  = 1.0E-13, &
         maxiter  = 100000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc%pvar,Timedisc%cvar)

    ! initialize log input/output
    CALL InitLogio(Logio,Mesh,Physics,Timedisc, &
         logformat = NOLOG, &
         filename  = "sedov3d.log", &
         logdt     = 300)

    ! set parameters for data output
    CALL InitOutput(Output,Mesh,Physics,Timedisc,&
         filetype  = GNUPLOT, &
         filename  = "sedov3d.dat", &
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
    INTEGER           :: n
    REAL              :: dr,P0,E
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: cart
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(OUT)       :: pvar,cvar
    !------------------------------------------------------------------------!

    ! 3D Sedov explosion
    n  = 3    ! 3D
    E  = 1.0  ! energy

    pvar(:,:,1) = 1.
    pvar(:,:,2:4) = 0.
    pvar(:,:,5) = 1.0D-05
     
    CALL Convert2Cartesian(Mesh%geometry,Mesh%bcenter,cart)

    ! should be at least 5 cells
    dr = 0.04
    P0 = 3.*(Physics%gamma-1.0) * E / ((n + 1) * PI * dr**n)
    WHERE (SQRT(cart(:,:,1)**2+cart(:,:,2)**2).LE.dr)
       pvar(:,:,5) = P0
    END WHERE
    
    CALL Convert2Conservative(Physics,Mesh,pvar,cvar)
    PRINT "(A,A)", " DATA-----> initial condition  ", &
         "3D Sedov explosion"

  END SUBROUTINE InitData

END MODULE Init
