!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: timedisc_modeuler.f90                                             #
!#                                                                           #
!# Copyright (C) 2007-2008                                                   #
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
       ! methods 
       InitTimedisc_modeuler, &
       CloseTimedisc_modeuler, &
       SolveODE_modeuler, &
       SetBoundaries_modeuler, &
       GetType, &
       GetName, &
       GetOrder, &
       GetCFL, &
       GetRank, &
       GetNumProcs, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitTimedisc_modeuler(this,os,order,stoptime,cfl,dtlimit,maxiter)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    INTEGER            :: os, order,maxiter
    REAL               :: stoptime,cfl,dtlimit
    !------------------------------------------------------------------------!
    INTEGER            :: err
    !------------------------------------------------------------------------!
    INTENT(IN)         :: os,order,stoptime,cfl,dtlimit,maxiter
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    CALL InitTimedisc(this,os,ODEsolver_name,order,stoptime,cfl,dtlimit,maxiter)

    ! allocate memory parameter matrix
    ALLOCATE(this%eta(3,3),STAT=err) 
    IF (err.NE.0) THEN
       CALL Error(this,"InitTimedisc_modeuler", "Unable to allocate memory.")
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
    INTEGER            :: i
    INTENT(IN)         :: Mesh,Physics,Fluxes
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!

    IF (PrimRecon(Fluxes%reconstruction)) THEN
       ! set center boundaries for primitive vars
       CALL Convert2Primitive(Physics,Mesh,this%cvar,this%pvar)
       CALL CenterBoundary(this%boundary,Mesh,Physics,this%time,this%pvar)
    ELSE
       ! set center boundaries for conservative vars
       CALL CenterBoundary(this%boundary,Mesh,Physics,this%time,this%cvar)
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
    INTEGER, SAVE      :: count=0
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

       ! flux ballance in x-direction
       IF (Mesh%INUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR NODEP
             DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NOUNROLL
                DO i=Mesh%IMIN,Mesh%IMAX
                   this%dxflux(i,j,k) = this%xflux(i,j,k) - this%xflux(i-1,j,k)
                END DO
             END DO
          END DO
       ELSE
          this%dxflux(:,:,:) = 0.
       END IF

       ! flux ballance in y-direction
       IF (Mesh%JNUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR NOUNROLL
             DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
                DO i=Mesh%IMIN,Mesh%IMAX
                   this%dyflux(i,j,k) = this%yflux(i,j,k) - this%yflux(i,j-1,k)
                END DO
             END DO
          END DO
       ELSE
          this%dyflux(:,:,:) = 0.
       END IF

       ! time step update
       DO k=1,Physics%VNUM
          this%cnew(:,:,k) = this%cvar(:,:,k) - this%dt * &
              (Mesh%dydV(:,:)*this%dxflux(:,:,k) + Mesh%dxdV(:,:)*this%dyflux(:,:,k) &
            - this%geo_src(:,:,k) - this%src(:,:,k))
       END DO

       this%cvar(:,:,:) = this%eta(GetOrder(this),n)*this%cold(:,:,:) &
            + (1.0-this%eta(GetOrder(this),n)) * this%cnew(:,:,:)

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
