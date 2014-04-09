!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_KHI.f90                                                      #
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
! Program and data initialization for Kelvin-Helmholtz instability
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
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Output,Logio
    !------------------------------------------------------------------------!

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

    ! mesh settings
    CALL InitMesh(Mesh,Fluxes, &
         geometry = CARTESIAN, &
             inum = 100, &          ! resolution in x and            !
             jnum = 100, &          !   y direction                  !             
             xmin = -0.5, &
             xmax = 0.5, &
             ymin = -0.5, &
             ymax = 0.5)

    ! boundary conditions
    CALL InitBoundary(Mesh%boundary,Mesh,Physics,PERIODIC,WEST)
    CALL InitBoundary(Mesh%boundary,Mesh,Physics,PERIODIC,EAST)
    CALL InitBoundary(Mesh%boundary,Mesh,Physics,PERIODIC,SOUTH)
    CALL InitBoundary(Mesh%boundary,Mesh,Physics,PERIODIC,NORTH)

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = 5.0, &
         dtlimit  = 1.0E-4, &
         maxiter  = 100000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc%pvar,Timedisc%cvar)

    ! initialize log input/output
    CALL InitLogio(Logio,Mesh,Physics,Timedisc, &
         logformat = NOLOG, &
         filename  = "KHI.log", &
         logdt     = 300)

    ! set parameters for data output
    CALL InitOutput(Output,Mesh,Physics,Timedisc,&
         filetype  = GNUPLOT, &
         filename  = "KHI.dat", &
         mode      = OVERWRITE, &
         starttime = 0.0, &
         stoptime  = Timedisc%stoptime, &
         count     = 50)

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
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2):: dv
    REAL              :: rho0, rho1, v0, v1, P0, P1
    REAL              :: ylen
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh, Physics
    INTENT(OUT)       :: pvar, cvar
    !------------------------------------------------------------------------!
    
    ! outer regions
    rho0 = 1.0
    v0   = -0.5
    P0   = 2.5

    ! inner region
    rho1 = 2.0
    v1   = 0.5
    P1   = P0

    ! y-velocity vanishes everywhere
    pvar(:,:,3) = 0.

    ylen = ABS(Mesh%ymax-Mesh%ymin)

    WHERE ((Mesh%bcenter(:,:,2).LT.(Mesh%ymin+0.25*ylen)).OR. &
         (Mesh%bcenter(:,:,2).GT.(Mesh%ymin+0.75*ylen)))
       pvar(:,:,1) = rho0
       pvar(:,:,2) = v0
       pvar(:,:,4) = P0
    ELSEWHERE
       pvar(:,:,1) = rho1
       pvar(:,:,2) = v1
       pvar(:,:,4) = P1
    END WHERE
       
    ! add velocity perturbations
    CALL RANDOM_SEED
    CALL RANDOM_NUMBER(dv)
    pvar(:,:,2:3) = pvar(:,:,2:3) + (dv(:,:,:)-0.5)*0.02

    CALL Convert2Conservative(Physics,Mesh,pvar,cvar)
    PRINT "(A,A)", " DATA-----> initial condition: ", &
         "Kelvin-Helmholtz instability"

  END SUBROUTINE InitData

END MODULE Init
