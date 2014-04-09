!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_sedov2d.f90                                                  #
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
! Program and data initialization for 2D Sedov explosion
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
         theta     = 1.3)           ! optional parameter for limiter !

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
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,NO_GRADIENTS,WEST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,NO_GRADIENTS,EAST)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,NO_GRADIENTS,SOUTH)
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,NO_GRADIENTS,NORTH)
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
       CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,WEST)
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
         stoptime = 0.05, &
         dtlimit  = 1.0E-13, &
         maxiter  = 100000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc%pvar,Timedisc%cvar)

    ! initialize log input/output
    CALL InitLogio(Logio,Mesh,Physics,Timedisc, &
         logformat = NOLOG, &
         filename  = "sedov2d.log", &
         logdt     = 300)

    ! set parameters for data output
    CALL InitOutput(Output,Mesh,Physics,Timedisc,&
         filetype  = GNUPLOT, &
         filename  = "sedov2d.dat", &
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
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(OUT)       :: pvar,cvar
    !------------------------------------------------------------------------!

    ! 2D Sedov explosion
    n  = 2    ! 2D
    E  = 1.0  ! energy

    pvar(:,:,1) = 1.
    pvar(:,:,2) = 0.
    pvar(:,:,3) = 0.
    pvar(:,:,4) = 1.0D-05
     
    SELECT CASE(GetType(Mesh%geometry))
    CASE(CARTESIAN)
       ! should be at least 5 cells
!       dr = 3.5*SQRT(Mesh%dx**2+Mesh%dy**2)
       dr = 0.03
       P0 = 3.*(Physics%gamma-1.0) * E / ((n + 1) * PI * dr**n)
       WHERE ((Mesh%bcenter(:,:,1)**2 + Mesh%bcenter(:,:,2)**2).LE.dr**2)
          pvar(:,:,4) = P0
       END WHERE
    CASE(POLAR)
       ! should be at least 5 cells
       ! be careful, if you wish to compare with the cartesian
       ! results dr should be the same
       dr = 0.03
       P0 = 3.*(Physics%gamma-1.0) * E / ((n + 1) * PI * dr**n)
       WHERE (Mesh%bcenter(:,:,1).LE.dr)
          pvar(:,:,4) = P0
       END WHERE
    CASE DEFAULT
       PRINT *, "ERROR in InitData: geometry should be either cartesian or polar"
       STOP
    END SELECT

    CALL Convert2Conservative(Physics,Mesh,pvar,cvar)
    PRINT "(A,A)", " DATA-----> initial condition  ", &
         "2D Sedov explosion"

  END SUBROUTINE InitData


END MODULE Init
