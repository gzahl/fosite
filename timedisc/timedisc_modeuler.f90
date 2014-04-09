!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: timedisc_modeuler.f90                                             #
!#                                                                           #
!# Copyright (C) 2007 - 2010                                                 #
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
  REAL, DIMENSION(3) :: eta
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Timedisc_TYP, &
       ! methods 
       InitTimedisc_modeuler, &
       CloseTimedisc_modeuler, &
       SolveODE_modeuler, &
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

    ! initialize parameter vector    
    SELECT CASE(GetOrder(this))
    CASE(1) ! time order=1 -> forward Euler
       eta(:) = (/ 0.0, 0.0, 0.0 /)
    CASE(2) ! time order=2 -> two step mod. Euler
       eta(:) = (/ 0.0, 0.5, 0.0 /)
    CASE(3) ! time order=3 -> three step mod. Euler
       eta(:) = (/ 0.0, 0.75, 1.0/3.0 /)
    CASE DEFAULT
       CALL Error(this,"InitTimedisc_modeuler","time order must be one of 1,2,3")
    END SELECT
  END SUBROUTINE InitTimedisc_modeuler


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

       DO k=1,Physics%VNUM
          ! flux ballance in x-direction
          IF (Mesh%INUM.GT.1) THEN
!CDIR OUTERUNROLL=8
             DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
                DO i=Mesh%IMIN,Mesh%IMAX
                   this%dxflux(i,j,k) = this%xflux(i,j,k) - this%xflux(i-1,j,k)
                END DO
                ! update western/eastern boundary fluxes (for accounting only)
                Fluxes%bxflux(j,1,k) = UpdateTimestep(eta(n),this%dt, Fluxes%bxfold(j,1,k), &
                     Fluxes%bxflux(j,1,k), Mesh%dy * this%xflux(Mesh%IMIN-1,j,k))
                Fluxes%bxflux(j,2,k) = UpdateTimestep(eta(n),this%dt, Fluxes%bxfold(j,2,k), &
                     Fluxes%bxflux(j,2,k), -Mesh%dy * this%xflux(Mesh%IMAX,j,k))
             END DO
          ELSE
             this%dxflux(:,:,k) = 0.
          END IF

          ! flux ballance in y-direction
          IF (Mesh%JNUM.GT.1) THEN
!CDIR OUTERUNROLL=8
             DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
                DO i=Mesh%IMIN,Mesh%IMAX
                   this%dyflux(i,j,k) = this%yflux(i,j,k) - this%yflux(i,j-1,k)
                END DO
             END DO
!CDIR NODEP
             DO i=Mesh%IMIN,Mesh%IMAX
                ! update southern/northern boundary fluxes (for accounting only)
                Fluxes%byflux(i,1,k) = UpdateTimestep(eta(n),this%dt,Fluxes%byfold(i,1,k), &
                     Fluxes%byflux(i,1,k), Mesh%dx * this%yflux(i,Mesh%JMIN-1,k))
                Fluxes%byflux(i,2,k) = UpdateTimestep(eta(n),this%dt,Fluxes%byfold(i,2,k), &
                     Fluxes%byflux(i,2,k), -Mesh%dx * this%yflux(i,Mesh%JMAX,k))
             END DO
          ELSE
             this%dyflux(:,:,k) = 0.
          END IF

          ! time step update of cell mean values
!CDIR OUTERUNROLL=8
          DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
             DO i=Mesh%IMIN,Mesh%IMAX
                this%cvar(i,j,k) = UpdateTimestep(eta(n),this%dt,this%cold(i,j,k),this%cvar(i,j,k), &
                       Mesh%dydV(i,j)*this%dxflux(i,j,k) + Mesh%dxdV(i,j)*this%dyflux(i,j,k) &
                     - this%geo_src(i,j,k) - this%src(i,j,k))
             END DO
          END DO
       END DO

       ! set boundary values
       CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,this%time,this%pvar,this%cvar)
    END DO

    CONTAINS

      ELEMENTAL FUNCTION UpdateTimestep(eta_n,dt,y0,yn,rhs) RESULT(y)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        REAL, INTENT(IN) :: eta_n,dt,y0,yn,rhs
        REAL :: y
        !------------------------------------------------------------------------!
        y = eta_n * y0 + (1.0 - eta_n) * (yn - dt * rhs)
      END FUNCTION UpdateTimestep

  END SUBROUTINE SolveODE_modeuler

  
  SUBROUTINE CloseTimedisc_modeuler(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP)   :: this
    !------------------------------------------------------------------------!
  END SUBROUTINE CloseTimedisc_modeuler

END MODULE timedisc_modeuler
