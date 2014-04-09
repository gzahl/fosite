!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: timedisc_modeuler.f90                                             #
!#                                                                           #
!# Copyright (C) 2007-2012                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
!CDIR IEXPAND
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
    INTEGER            :: n,k
    INTEGER            :: order
    REAL               :: t,dtnew
    REAL               :: rel_err(Physics%VNUM)
    TYPE var_typ
       REAL, DIMENSION(:,:,:), POINTER :: var
    END TYPE var_typ
    TYPE(var_typ)      :: p(4),c(4)
    !------------------------------------------------------------------------!    
    INTENT(IN)         :: Mesh,time
    INTENT(INOUT)      :: this,Physics,Fluxes,dt,maxerr
    !------------------------------------------------------------------------!
    t = time
!CDIR IEXPAND
    order = GetOrder(this)
    ! check if adaptive step size control is enabled
    IF (this%tol_rel.GT.1.0) THEN
       ! no adaptive step size control
!CDIR UNROLL=3
       DO n=1,order
          ! update time variable
          t = time+zeta(n,order)*dt
          ! compute right hand side for explicit Euler step
          ! time step update of cell mean values
          CALL ComputeCVar_modeuler(this,Mesh,Physics,Fluxes,eta(n,order), &
               time,dt,this%cold,this%pvar,this%cvar,this%rhs,this%cvar)
          ! set boundary values
          CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,t,this%pvar,this%cvar)
       END DO
       maxerr = 0.0
    ELSE
       p(1)%var => this%pvar
       p(2)%var => this%ptmp  ! store intermediate result for error control
       p(3)%var => this%pvar
       p(4)%var => this%pvar
       c(1)%var => this%cvar
       c(2)%var => this%ctmp  ! store intermediate result for error control
       c(3)%var => this%cvar          
       c(4)%var => this%cvar          
!CDIR UNROLL=3
       DO n=1,order
          ! update time variable
          t = time+zeta(n,order)*dt
          ! compute right hand side and update cvar and bfluxes
          CALL ComputeCVar_modeuler(this,Mesh,Physics,Fluxes,eta(n,order), &
               time,dt,this%cold,p(n)%var,c(n)%var,this%rhs,c(n+1)%var)
          ! for 3rd order scheme compute the 2nd order result with the same RHS
          ! and store it in this%ctmp, bfluxes and boundary conditions are not required
          IF (n.EQ.2.AND.order.EQ.3) &
!CDIR IEXPAND
             this%ctmp(:,:,:) = UpdateTimestep_modeuler(eta(2,2),dt,this%cold(:,:,:), &
                                this%ctmp(:,:,:),this%rhs(:,:,:))
          ! set boundary values
          CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,t,p(n+1)%var,c(n+1)%var)
       END DO
       ! maximum of truncation error
       DO k=1,Physics%VNUM
          rel_err(k) = MAXVAL(ABS(this%cvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k) &
                                 -this%ctmp(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k)) &
               / (this%tol_rel*ABS(this%cvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k)) + this%tol_abs(k)))
       END DO
       maxerr = MAXVAL(rel_err(:))

       ! compute new time step
       dtnew = 0.9*dt*exp(-log(maxerr)/order)
       IF (maxerr.LT.1.0) THEN
          dt = MIN(dtnew,4.*dt)   ! not too large
       ELSE
          dt = MAX(dtnew,0.25*dt) ! not too small
       END IF
!!$       PRINT '(7(ES14.6))',time,dt,maxerr,rel_err(:)
    END IF
  END SUBROUTINE SolveODE_modeuler

  
  SUBROUTINE ComputeCVar_modeuler(this,Mesh,Physics,Fluxes,eta,time,dt, &
                                  cold,pvar,cvar,rhs,cnew)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    REAL               :: eta,time,dt
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                       :: cold,pvar,cvar,cnew,rhs
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,eta,time,dt,cold,pvar,cvar
    INTENT(INOUT)      :: this,Physics,Fluxes
    INTENT(OUT)        :: cnew,rhs
    !------------------------------------------------------------------------!
    ! get the numerical fluxes
    CALL CalculateFluxes(Fluxes,Mesh,Physics,pvar,cvar,this%xflux,this%yflux)

    ! get geometrical sources for non-cartesian mesh
!CDIR IEXPAND
    IF (GetType(Mesh%geometry).NE.CARTESIAN) &
       CALL GeometricalSources(Physics,Mesh,Fluxes,pvar,cvar,this%geo_src)

    ! get source terms due to external forces if present
    IF (ASSOCIATED(Physics%sources)) &
       CALL ExternalSources(Physics%sources,Mesh,Fluxes,Physics, &
            time,pvar,cvar,this%src)

!CDIR NOVECTOR
    DO k=1,Physics%VNUM
!CDIR OUTERUNROLL=8
!CDIR NODEP
       DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
          DO i=Mesh%IMIN,Mesh%IMAX
             ! update right hand side of ODE
             rhs(i,j,k) = Mesh%dydV(i,j)*(this%xflux(i,j,k) - this%xflux(i-1,j,k)) &
                  + Mesh%dxdV(i,j)*(this%yflux(i,j,k) - this%yflux(i,j-1,k)) &
                  - this%geo_src(i,j,k) - this%src(i,j,k)
             ! update cvar -> cnew
!CDIR IEXPAND
             cnew(i,j,k) = UpdateTimestep_modeuler(eta,dt,cold(i,j,k),cvar(i,j,k),rhs(i,j,k))
          END DO
          ! compute RHS for boundary fluxes
          ! western and eastern
!CDIR IEXPAND
          Fluxes%bxflux(j,1,k) = UpdateTimestep_modeuler(eta,dt,Fluxes%bxfold(j,1,k), &
               Fluxes%bxflux(j,1,k),Mesh%dy*this%xflux(Mesh%IMIN-1,j,k))
!CDIR IEXPAND
          Fluxes%bxflux(j,2,k) = UpdateTimestep_modeuler(eta,dt,Fluxes%bxfold(j,2,k), &
               Fluxes%bxflux(j,2,k),-Mesh%dy*this%xflux(Mesh%IMAX,j,k))
       END DO

       ! compute RHS for boundary fluxes
       ! southern and northern
!CDIR NODEP
       DO i=Mesh%IMIN,Mesh%IMAX
 !CDIR IEXPAND
          Fluxes%byflux(i,1,k) = UpdateTimestep_modeuler(eta,dt,Fluxes%byfold(i,1,k), &
               Fluxes%byflux(i,1,k),Mesh%dx*this%yflux(i,Mesh%JMIN-1,k))
!CDIR IEXPAND
          Fluxes%byflux(i,2,k) = UpdateTimestep_modeuler(eta,dt,Fluxes%byfold(i,2,k), &
               Fluxes%byflux(i,2,k),-Mesh%dx*this%yflux(i,Mesh%JMAX,k))
       END DO
     END DO
  END SUBROUTINE ComputeCVar_modeuler




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
    ! ATTENTION:
    ! The time step update is computed according to:
    !    y = eta_n * y0 + (1.0 - eta_n) * (yn - dt * rhs)
    ! but to minimize the truncation error it is essential to sort the terms
    ! in this way:
    y = yn-dt*rhs+eta_n*(y0-yn+dt*rhs)
  END FUNCTION UpdateTimestep_modeuler


END MODULE timedisc_modeuler
