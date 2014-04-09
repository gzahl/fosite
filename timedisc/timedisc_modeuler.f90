!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: timedisc_modeuler.f90                                             #
!#                                                                           #
!# Copyright (C) 2007-2011                                                   #
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
! Reference: Shu & Osher, J. Comput Phys. 77(2), 439 (1988)
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
  REAL, PARAMETER :: EPS = 1.0E-1
  REAL, PARAMETER :: eta(3,3) = &
       RESHAPE((/ 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.75, 1.0/3.0 /),(/3,3/))
  REAL, PARAMETER :: zeta(3,3) = &
       RESHAPE((/ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.5 /),(/3,3/))
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Timedisc_TYP, &
       ! methods 
       InitTimedisc_modeuler, &
       CloseTimedisc_modeuler, &
       SolveODE_modeuler, &
       CloseTimedisc, &
       GetOrder, &
       GetCFL, &
       GetType, &
       GetName, &
       GetRank, &
       GetNumProcs, &
       Initialized, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitTimedisc_modeuler(this,os,order,stoptime,cfl,dtlimit,maxiter, &
       tol_rel,tol_abs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    INTEGER            :: os, order,maxiter
    REAL               :: stoptime,cfl,dtlimit
    REAL               :: tol_rel
    REAL,DIMENSION(:)  :: tol_abs
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    INTENT(IN)         :: os,order,stoptime,cfl,dtlimit,maxiter,tol_rel,tol_abs
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    CALL InitTimedisc(this,os,ODEsolver_name,order,stoptime,cfl,dtlimit,maxiter)

    ! relative and absolute tolerance for adaptive step size control
    this%tol_rel    = tol_rel
    this%tol_abs(:) = tol_abs(:)
    SELECT CASE(GetOrder(this))    
    CASE(1)
       ! set relative error tolarance to value > 1.0
       ! to disable adaptive step size control
       ! (not available for 1st order scheme)
       this%tol_rel    = 10.0
       this%tol_abs(:) = 1.0
      CALL Warning(this,"InitTimedisc_modeuler", &
            "adaptive step size control not supported in 1st order scheme")
    CASE(2,3)
       IF (tol_rel.GT.1.0) &
            CALL Warning(this,"InitTimedisc_modeuler", &
            "adaptive step size control disabled (tol_rel>1)")
    CASE DEFAULT
       CALL Error(this,"InitTimedisc_modeuler","time order must be one of 1,2,3")
    END SELECT
    IF ((this%tol_rel.LT.0.0).OR.MINVAL(this%tol_abs(:)).LT.0.0) &
         CALL Error(this,"InitTimedisc_modeuler", &
         "error tolerance levels must be greater than 0")
  END SUBROUTINE InitTimedisc_modeuler


  SUBROUTINE SolveODE_modeuler(this,Mesh,Physics,Fluxes,time,dt,maxerr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    REAL               :: time,dt,maxerr
    !------------------------------------------------------------------------!
    INTEGER            :: n,i,j,k
    REAL               :: t,dtnew,dtold
    REAL               :: rel_err(Physics%VNUM)
    !------------------------------------------------------------------------!    
    INTENT(IN)         :: Mesh,time
    INTENT(INOUT)      :: this,Physics,Fluxes,dt,maxerr
    !------------------------------------------------------------------------!
    t = time
    ! check if adaptive step size control is enabled
    IF (this%tol_rel.GT.1.0) THEN
       ! no adaptive step size control
       DO n=1,GetOrder(this)
          t = time+zeta(n,GetOrder(this))*dt
          CALL ComputeRHS_modeuler(this,Mesh,Physics,Fluxes,t, &
               this%pvar,this%cvar,this%rhs)
          ! time step update of cell mean values
          CALL ComputeCVar_modeuler(this,Mesh,Physics,eta(n,GetOrder(this)),dt,&
               this%rhs,this%cold,this%cvar,this%cvar)
          ! update boundary fluxes (for accounting only)
          CALL ComputeBFluxes_modeuler(this,Mesh,Physics,eta(n,GetOrder(this)),dt,&
               this%rhs,Fluxes%bxfold,Fluxes%byfold,Fluxes%bxflux,Fluxes%byflux)
          ! set boundary values
          CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,t,this%pvar,this%cvar)
       END DO
       maxerr = 0.0
    ELSE
       CALL ComputeRHS_modeuler(this,Mesh,Physics,Fluxes,t,this%pvar,this%cvar,this%rhs)
       ! time step update of cell mean values
!CDIR IEXPAND
       IF (GetOrder(this).EQ.2) THEN
          ! first order explicit Euler step
          CALL ComputeCVar_modeuler(this,Mesh,Physics,eta(1,2),dt,this%rhs,&
               this%cold,this%cvar,this%ctmp)
          ! update boundary fluxes (for accounting only)
          CALL ComputeBFluxes_modeuler(this,Mesh,Physics,eta(1,2),dt, &
               this%rhs,Fluxes%bxfold,Fluxes%byfold,Fluxes%bxflux,Fluxes%byflux)
          ! set boundary values
          CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,t,this%ptmp,this%ctmp)
          ! next time
          t = time+zeta(2,2)*dt
          CALL ComputeRHS_modeuler(this,Mesh,Physics,Fluxes,t, &
               this%ptmp,this%ctmp,this%rhs)
          ! second step for 2nd order modified Euler scheme
          CALL ComputeCVar_modeuler(this,Mesh,Physics,eta(2,2),dt,this%rhs,&
               this%cold,this%ctmp,this%cvar)
          ! set boundary values
          CALL ComputeBFluxes_modeuler(this,Mesh,Physics,eta(2,2),dt, &
               this%rhs,Fluxes%bxfold,Fluxes%byfold,Fluxes%bxflux,Fluxes%byflux)
       ELSE
          ! order is either 1 or 3
          ! first order explicit Euler step
          CALL ComputeCVar_modeuler(this,Mesh,Physics,eta(1,3),dt,this%rhs,&
               this%cold,this%cvar,this%cvar)
          ! update boundary fluxes (for accounting only)
          CALL ComputeBFluxes_modeuler(this,Mesh,Physics,eta(1,3),dt, &
               this%rhs,Fluxes%bxfold,Fluxes%byfold,Fluxes%bxflux,Fluxes%byflux)
!CDIR IEXPAND
          IF (GetOrder(this).EQ.3) THEN
             ! set boundary values
             CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,t, &
                  this%pvar,this%cvar)
             ! right hand side for second step of 2nd & 3rd order
             ! modified Euler scheme
             t = time+zeta(2,3)*dt ! = time+zeta(2,2)*dt
             CALL ComputeRHS_modeuler(this,Mesh,Physics,Fluxes,t, &
                  this%pvar,this%cvar,this%rhs)
             ! second step of 2nd oder scheme
             CALL ComputeCVar_modeuler(this,Mesh,Physics,eta(2,2),dt,this%rhs,&
                  this%cold,this%cvar,this%ctmp)
             ! no update of boundary fluxes for 2nd order scheme
             ! set boundary values for ctmp
             CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,t, &
                  this%ptmp,this%ctmp)
             ! second step of 3rd order scheme (uses same rhs)
             CALL ComputeCVar_modeuler(this,Mesh,Physics,eta(2,3),dt,this%rhs,&
                  this%cold,this%cvar,this%cvar)
             ! update boundary fluxes (for accounting only)
             CALL ComputeBFluxes_modeuler(this,Mesh,Physics,eta(2,3),dt, &
                  this%rhs,Fluxes%bxfold,Fluxes%byfold,Fluxes%bxflux,Fluxes%byflux)
             ! set boundary values
             CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,t, &
                  this%pvar,this%cvar)
             ! third step for 3rd order modified Euler scheme
             t = time+zeta(3,3)*dt
             CALL ComputeRHS_modeuler(this,Mesh,Physics,Fluxes,t, &
                  this%pvar,this%cvar,this%rhs)
             CALL ComputeCVar_modeuler(this,Mesh,Physics,eta(3,3),dt,this%rhs,&
                  this%cold,this%cvar,this%cvar)
             CALL ComputeBFluxes_modeuler(this,Mesh,Physics,eta(3,3),dt, &
                  this%rhs,Fluxes%bxfold,Fluxes%byfold,Fluxes%bxflux,Fluxes%byflux)
          END IF
       END IF
       ! set boundary values
       CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,t,this%pvar,this%cvar)
       
       ! maximum of truncation error
       DO k=1,Physics%VNUM
          rel_err(k) = MAXVAL(ABS(this%cvar(:,:,k)-this%ctmp(:,:,k)) &
               / (this%tol_rel*ABS(this%cvar(:,:,k)) + this%tol_abs(k)))
       END DO
       maxerr = MAXVAL(rel_err(:))

       ! compute new time step
       dtnew = 0.9*dt*(1./maxerr)**(1./(GetOrder(this)+1.0))
       IF (maxerr.LT.1.0) THEN
          dt = MIN(dtnew,4.*dt)   ! not too large
       ELSE
          dt = MAX(dtnew,0.25*dt) ! not too small
       END IF
!!$       PRINT '(7(ES14.6))',time,dt,maxerr,rel_err(:)
    END IF
  END SUBROUTINE SolveODE_modeuler

  
  SUBROUTINE ComputeCVar_modeuler(this,Mesh,Physics,eta,dt,rhs,cold,cvar,cnew)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL               :: eta,dt
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                       :: rhs,cold,cvar,cnew
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this,Mesh,Physics,eta,dt,rhs,cold
    INTENT(INOUT)      :: cvar,cnew
    !------------------------------------------------------------------------!
!CDIR NOVECTOR
    DO k=1,Physics%VNUM
!CDIR OUTERUNTOLL=8
       DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
          DO i=Mesh%IMIN,Mesh%IMAX
!CDIR IEXPAND
             cnew(i,j,k) = UpdateTimestep_modeuler(eta,dt, &
                  cold(i,j,k),cvar(i,j,k),rhs(i,j,k))
          END DO
       END DO
    END DO
  END SUBROUTINE ComputeCVar_modeuler


  SUBROUTINE ComputeBFluxes_modeuler(this,Mesh,Physics,eta,dt,rhs,bxfold,byfold, &
       bxflux,byflux)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL               :: eta,dt
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                       :: rhs
    REAL,DIMENSION(Mesh%JGMIN:Mesh%JGMAX,2,Physics%VNUM) :: bxfold,bxflux
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,2,Physics%VNUM) :: byfold,byflux
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this,Mesh,Physics,eta,dt,rhs,bxfold,byfold
    INTENT(INOUT)      :: bxflux,byflux
    !------------------------------------------------------------------------!
    DO k=1,Physics%VNUM
       ! western and eastern
!CDIR NODEP
       DO j=Mesh%JMIN,Mesh%JMAX
!CDIR IEXPAND
          bxflux(j,1,k) = UpdateTimestep_modeuler(eta,dt,bxfold(j,1,k), &
               bxflux(j,1,k),rhs(Mesh%IMIN-1,j,k))
!CDIR IEXPAND
          bxflux(j,2,k) = UpdateTimestep_modeuler(eta,dt,bxfold(j,2,k), &
               bxflux(j,2,k),rhs(Mesh%IMAX+1,j,k))
       END DO
       ! southern and northern
!CDIR NODEP
       DO i=Mesh%IMIN,Mesh%IMAX
!CDIR IEXPAND
          byflux(i,1,k) = UpdateTimestep_modeuler(eta,dt,byfold(i,1,k), &
               byflux(i,1,k),this%rhs(i,Mesh%JMIN-1,k))
!CDIR IEXPAND
          byflux(i,2,k) = UpdateTimestep_modeuler(eta,dt,byfold(i,2,k), &
               byflux(i,2,k),this%rhs(i,Mesh%JMIN-1,k))
       END DO
    END DO
  END SUBROUTINE ComputeBFluxes_modeuler


  SUBROUTINE ComputeRHS_modeuler(this,Mesh,Physics,Fluxes,time,pvar,cvar,rhs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    REAL               :: time
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                       :: pvar,cvar,rhs
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,k
    REAL               :: dyflux
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,time,pvar,cvar
    INTENT(INOUT)      :: this,Physics,Fluxes
    INTENT(OUT)        :: rhs
    !------------------------------------------------------------------------!
    ! get the numerical fluxes
    CALL CalculateFluxes(Fluxes,Mesh,Physics,pvar,cvar,this%xflux,this%yflux)

    ! get geometrical sources for non-cartesian mesh
    IF (GetType(Mesh%geometry).NE.CARTESIAN) &
       CALL GeometricalSources(Physics,Mesh,Fluxes,pvar,cvar,this%geo_src)

    ! get source terms due to external forces if present
    IF (ASSOCIATED(Physics%sources)) &
       CALL ExternalSources(Physics%sources,Mesh,Fluxes,Physics, &
            time,pvar,cvar,this%src)

    DO k=1,Physics%VNUM
       ! compute flux differences
       ! x-direction (if not 1D)
       IF (Mesh%INUM.GT.1) THEN
!CDIR OUTERUNROLL=8
          DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
             DO i=Mesh%IMIN,Mesh%IMAX
                ! temporary use rhs for flux difference in x-direction
                rhs(i,j,k) = Mesh%dydV(i,j)*( &
                     this%xflux(i,j,k) - this%xflux(i-1,j,k))
             END DO
          END DO
       ELSE
          rhs(:,:,:) = 0.0
       END IF
!CDIR COLLAPSE
       DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
          DO i=Mesh%IGMIN,Mesh%IGMAX
             ! one may exclude computation of dyflux for 1D computations
             ! but this prevents vectorization; thus we allways compute dyflux
             dyflux = Mesh%dxdV(i,j)*(this%yflux(i,j,k) - this%yflux(i,j-1,k))
             rhs(i,j,k) = rhs(i,j,k) & ! = dxflux (see above)
                  + dyflux - this%geo_src(i,j,k) - this%src(i,j,k)
          END DO
       END DO
       ! compute RHS for boundary fluxes
       ! western and eastern
!CDIR NODEP
       DO j=Mesh%JMIN,Mesh%JMAX
          rhs(Mesh%IMIN-1,j,k) = Mesh%dy * this%xflux(Mesh%IMIN-1,j,k)
          rhs(Mesh%IMAX+1,j,k) = -Mesh%dy * this%xflux(Mesh%IMAX,j,k)
       END DO
       ! southern and northern
!CDIR NODEP
       DO i=Mesh%IMIN,Mesh%IMAX
          rhs(i,Mesh%JMIN-1,k) = Mesh%dx * this%yflux(i,Mesh%JMIN-1,k)
          rhs(i,Mesh%JMAX+1,k) = -Mesh%dx * this%yflux(i,Mesh%JMAX,k)
       END DO
    END DO
  END SUBROUTINE ComputeRHS_modeuler


  SUBROUTINE CloseTimedisc_modeuler(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP)   :: this
    !------------------------------------------------------------------------!
  END SUBROUTINE CloseTimedisc_modeuler


  ELEMENTAL FUNCTION UpdateTimestep_modeuler(eta_n,dt,y0,yn,rhs) RESULT(y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: eta_n,dt,y0,yn,rhs
    REAL :: y
    !------------------------------------------------------------------------!
    y = eta_n * y0 + (1.0 - eta_n) * (yn - dt * rhs)
  END FUNCTION UpdateTimestep_modeuler


END MODULE timedisc_modeuler
