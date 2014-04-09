!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: timedisc_modeuler.f90                                             #
!#                                                                           #
!# Copyright (C) 2007 Tobias Illenseer <tillense@ita.uni-heidelberg.de>      #
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
! subroutines for modified Euler i.e. Runge-Kutta methods
!----------------------------------------------------------------------------!
MODULE timedisc_modeuler
  USE timedisc_common
  USE mesh_generic
  USE fluxes_generic
  USE boundary_generic
  USE physics_generic, GeometricalSources_Physics => GeometricalSources, &
       ExternalSources_Physics => ExternalSources
  USE sources_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: ODEsolver_name = "modified Euler"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Timedisc_TYP, &
       ! constants
       ! methods 
       InitTimedisc_modeuler, &
       SolveODE_modeuler, &
       SetBoundaries_modeuler, &
       GetType, &
       GetName, &
       GetOrder, &
       GetCFL, &
       CloseTimedisc_modeuler
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitTimedisc_modeuler(this,os,order,cfl,stoptime,dtlimit,maxiter)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    INTEGER            :: os, order,maxiter
    REAL               :: cfl,stoptime,dtlimit
    !------------------------------------------------------------------------!
    INTEGER            :: err
    !------------------------------------------------------------------------!
    INTENT(IN)         :: os,order,cfl,stoptime,dtlimit,maxiter
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    CALL InitTimedisc(this,os,ODEsolver_name,order,cfl,stoptime,dtlimit,maxiter)

    ! allocate memory parameter matrix
    ALLOCATE(this%eta(3,3),STAT=err) 
    IF (err.NE.0) THEN
       PRINT *, "ERROR in InitTimedisc_modeuler: Unable to allocate memory!"
       STOP
    END IF

    ! initialize parameter matrix
    this%eta(1,:) = 0.0                       ! time order=1 -> forward Euler!
    
    this%eta(2,1) = 0.0                       ! time order=2 ->              !
    this%eta(2,2) = 0.5                       !          two step mod. Euler !
    this%eta(2,3) = 0.0
    
    this%eta(3,1) = 0.0                       ! time order=3 ->              !
    this%eta(3,2) = 0.75                      !        three step mod. Euler !
    this%eta(3,3) = 1.0/3.0
  END SUBROUTINE InitTimedisc_modeuler


  SUBROUTINE SetBoundaries_modeuler(this,Mesh,Physics,Fluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    !------------------------------------------------------------------------!    
    INTENT(IN)         :: Mesh,Fluxes
    INTENT(INOUT)      :: this,Physics
    !------------------------------------------------------------------------!

    IF (PrimRecon(Fluxes%reconstruction)) THEN
       ! set center boundaries for primitive vars
       CALL Convert2Primitive(Physics,Mesh,this%cvar,this%pvar)
       CALL CenterBoundary(Mesh%boundary,Mesh,Physics,this%pvar)
    ELSE
       ! set center boundaries for conservative vars
       CALL CenterBoundary(Mesh%boundary,Mesh,Physics,this%cvar)
       CALL Convert2Primitive(Physics,Mesh,this%cvar,this%pvar)
    END IF
  END SUBROUTINE SetBoundaries_modeuler


  SUBROUTINE SolveODE_modeuler(this,Mesh,Physics,Fluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    !------------------------------------------------------------------------!
    INTEGER            :: n,i,j,k
    !------------------------------------------------------------------------!    
    INTENT(IN)         :: Mesh
    INTENT(INOUT)      :: this,Physics,Fluxes
    !------------------------------------------------------------------------!

    DO n=1,GetOrder(this)
       ! get the numerical fluxes
       CALL CalculateFluxes(Fluxes,Mesh,Physics, &
            this%pvar,this%cvar,this%xflux,this%yflux)
       
       ! get geometrical sources for non-cartesian mesh
       IF (GetType(Mesh%geometry).NE.CARTESIAN) THEN 
          CALL GeometricalSources(Physics,Mesh,Fluxes, &
               this%pvar,this%cvar,this%geo_src)
       END IF
       
       ! get source terms due to external forces if present
       IF (ASSOCIATED(Physics%sources)) THEN
          CALL ExternalSources(Physics%sources,Mesh,Fluxes,Physics, &
               this%pvar,this%cvar,this%src)
       END IF
       
       ! time step update
       FORALL (i=Mesh%IMIN:Mesh%IMAX, j=Mesh%JMIN:Mesh%JMAX, k=1:Physics%vnum)
          this%cvar(i,j,k) = this%eta(GetOrder(this),n)*this%cold(i,j,k) &
               + (1.0-this%eta(GetOrder(this),n)) * (this%cvar(i,j,k) &
               - this%dt * ( &
               (this%xflux(i,j,k) - this%xflux(i-1,j,k)) * Mesh%dydV(i,j) &
               + (this%yflux(i,j,k) - this%yflux(i,j-1,k)) * Mesh%dxdV(i,j) &
               - this%geo_src(i,j,k) - this%src(i,j,k)) &
               )
       END FORALL
       
       ! set boundary values
       CALL SetBoundaries_modeuler(this,Mesh,Physics,Fluxes)
    END DO
  END SUBROUTINE SolveODE_modeuler


  SUBROUTINE CloseTimedisc_modeuler(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP)   :: this
    !------------------------------------------------------------------------!

    DEALLOCATE(this%eta)
  END SUBROUTINE CloseTimedisc_modeuler

END MODULE timedisc_modeuler
