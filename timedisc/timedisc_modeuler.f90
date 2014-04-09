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
!> \author Tobias Illenseer
!! \author Björn Sperling
!!
!! \brief subroutines for modified Euler i.e. Runge-Kutta methods
!!
!! Reference: Shu & Osher, J. Comput Phys. 77(2), 439 (1988)
!!
!! \extends timedisc_common
!! \ingroup timedisc
!----------------------------------------------------------------------------!
MODULE timedisc_modeuler
  USE timedisc_common
  USE mesh_generic
  USE fluxes_generic
  USE boundary_generic
  USE physics_generic, GeometricalSources_Physics => GeometricalSources, &
       ExternalSources_Physics => ExternalSources
  USE sources_generic, CalcTimestep_sources => CalcTimestep
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: ODEsolver_name = "modified Euler"
  REAL, PARAMETER :: eta(3,3) = &
       RESHAPE((/ 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.75, 1.0/3.0 /),(/3,3/))
  REAL, PARAMETER :: zeta(3,3) = &
       RESHAPE((/ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.5 /),(/3,3/))
  INTEGER, PARAMETER :: DTCAUSE_SMALLERR = -2  ! smallest ts due to err      !
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Timedisc_TYP, &
       ! constants
       DTCAUSE_SMALLERR, &
       ! methods 
       InitTimedisc_modeuler, &
       CloseTimedisc_modeuler, &
       SolveODE_modeuler, &
       CalcTimestep_modeuler, &
       ComputeError_modeuler, &
       ComputeSources_modeuler, &
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

  SUBROUTINE InitTimedisc_modeuler(this,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Dict_TYP), POINTER &
                       :: config
    !------------------------------------------------------------------------!
    INTEGER            :: method
    !------------------------------------------------------------------------!
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    ! set default order 
    CALL RequireKey(config, "order", 3)
    CALL GetAttr(config, "order", this%order)
    
    CALL GetAttr(config, "method", method)
    CALL InitTimedisc(this,method,ODEsolver_name)

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
       IF (this%tol_rel.GT.1.0) &
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
    IF (this%tol_rel.GE.1.0) THEN
       ! no adaptive step size control
!CDIR UNROLL=3
       DO n=1,order
          ! update time variable
          t = time+zeta(n,order)*dt
          ! compute right hand side for explicit Euler step
          ! time step update of cell mean values
          CALL ComputeCVar_modeuler(this,Mesh,Physics,Fluxes,eta(n,order), &
               t,dt,this%cold,this%pvar,this%cvar,this%rhs,this%cvar)
          ! set boundary values
          CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,t,this%pvar,this%cvar)
       END DO
       maxerr = 0.0
       dt = HUGE(dt)
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
               t,dt,this%cold,p(n)%var,c(n)%var,this%rhs,c(n+1)%var)
          ! for 3rd order scheme compute the 2nd order result with the same RHS
          ! and store it in this%ctmp, bfluxes and boundary conditions are not required
          IF (n.EQ.2.AND.order.EQ.3) &
!CDIR IEXPAND
             this%ctmp(:,:,:) = UpdateTimestep_modeuler(eta(2,2),dt,this%cold(:,:,:), &
                                this%ctmp(:,:,:),this%rhs(:,:,:))
          ! set boundary values
          CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,t,p(n+1)%var,c(n+1)%var)
       END DO
       ! maxerr and dt are global values (MPI)
       CALL ComputeError_modeuler(this,Mesh,Physics,dt,maxerr)
    END IF

  END SUBROUTINE SolveODE_modeuler


  SUBROUTINE ComputeError_modeuler(this,Mesh,Physics,dt,maxerr)
#ifdef PARALLEL
#ifdef HAVE_MPI_MOD
  USE mpi
#endif
#endif
  IMPLICIT NONE
#ifdef PARALLEL
#ifdef HAVE_MPIF_H
  include 'mpif.h'
#endif
#endif
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL               :: dt,maxerr
    !------------------------------------------------------------------------!
    REAL               :: dtnew,error,alpha
    REAL               :: rel_err(Physics%VNUM)
    INTEGER            :: k,l(1),i,j
#ifdef PARALLEL
    INTEGER            :: ierror
    REAL               :: maxerr_all
#endif
    !------------------------------------------------------------------------!    
    INTENT(IN)         :: Mesh,Physics
    INTENT(INOUT)      :: this,dt,maxerr
    !------------------------------------------------------------------------!
  
    IF (this%tol_rel.LT.1.0) THEN
      ! maximum of truncation error
      IF(this%write_error) THEN
        DO k=1,Physics%VNUM
          rel_err(k) = 0.
          DO j=Mesh%JMIN,Mesh%JMAX
            DO i=Mesh%IMIN,Mesh%IMAX
              error = ABS(this%cvar(i,j,k) - this%ctmp(i,j,k)) &
                      / (this%tol_rel*ABS(this%cvar(i,j,k)) + this%tol_abs(k))
              rel_err(k) = MAX(rel_err(k),error)
              this%error(i,j,k) = MAX(this%error(i,j,k),error)
            END DO
          END DO
        END DO
      ELSE
        DO k=1,Physics%VNUM
          rel_err(k) = MAXVAL(ABS(this%cvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k) &
                                -this%ctmp(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k)) &
                            / (this%tol_rel &
                               *ABS(this%cvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k)) &
                               + this%tol_abs(k)))
        END DO
      END IF
      maxerr = MAXVAL(rel_err(:))

#ifdef PARALLEL
       CALL MPI_Allreduce(maxerr,maxerr_all,1,DEFAULT_MPI_REAL,MPI_MAX,&
            Mesh%comm_cart,ierror)
       maxerr = maxerr_all
#endif

      !compute new time step
!      dtnew = 0.9*dt*(1./maxerr)**(1./(GetOrder(this)))
!CDIR IEXPAND
      ! see E. Hairer, Solving Ordinary Differential Equ. II, 2ed, Springer (2.43c)
      alpha = 1./GetOrder(this) + this%beta
      IF(this%maxerrold.GT.0.) THEN
        dtnew = 0.9*dt*exp(-log(maxerr)*alpha+log(this%maxerrold)*this%beta)
      ELSE
        dtnew = 0.9*dt*exp(-log(maxerr)/GetOrder(this))
      END IF
      IF (maxerr.LT.1.0) THEN
         dt = MIN(dtnew,4.*dt)   ! not too large
         ! it is possible: maxerr < 1 =>  0.9*dt < dtnew => maybe dtnew < dt
         IF (dtnew .LT. dt) this%dtcause = DTCAUSE_SMALLERR
      ELSE
         !print *,maxerr,l
         dt = MAX(dtnew,0.25*dt) ! not too small
      END IF
!         PRINT '(7(ES14.6))',time,dt,maxerr,rel_err(:)
      this%maxerrold = maxerr
    ELSE
      ! no adaptive step size control
      maxerr = 0.0
      dt = HUGE(dt)
    END IF
  END SUBROUTINE ComputeError_modeuler


  SUBROUTINE ComputeSources_modeuler(this,Mesh,Physics,Fluxes,time,dt,pvar,cvar,geosrc,src)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    REAL               :: time,dt
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                       :: pvar,cvar,geosrc,src
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,time,dt,pvar,cvar
    INTENT(INOUT)      :: this,Physics,Fluxes
    INTENT(OUT)        :: geosrc,src
    !------------------------------------------------------------------------!
  
    ! get geometrical sources
    CALL GeometricalSources(Physics,Mesh,Fluxes,pvar,cvar,geosrc)

    ! get source terms due to external forces if present
    IF (ASSOCIATED(Physics%sources)) &
       CALL ExternalSources(Physics%sources,Mesh,Fluxes,Physics, &
            time,dt,pvar,cvar,src)

  END SUBROUTINE ComputeSources_modeuler
  
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
    REAL               :: dyflux
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,eta,time,dt,cold,pvar,cvar
    INTENT(INOUT)      :: this,Physics,Fluxes
    INTENT(OUT)        :: cnew,rhs
    !------------------------------------------------------------------------!
    ! get the numerical fluxes
    CALL CalculateFluxes(Fluxes,Mesh,Physics,pvar,cvar,this%xfluxdy,this%yfluxdx)

    CALL ComputeSources_modeuler(this,Mesh,Physics,Fluxes,time,dt,pvar,cvar,&
      this%geo_src,this%src)

!CDIR NOVECTOR
    DO k=1,Physics%VNUM
!CDIR OUTERUNROLL=8
       DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
          DO i=Mesh%IMIN,Mesh%IMAX
             ! update right hand side of ODE
             rhs(i,j,k) = Mesh%dydV(i,j)*(this%xfluxdy(i,j,k) - this%xfluxdy(i-1,j,k)) &
                  + Mesh%dxdV(i,j)*(this%yfluxdx(i,j,k) - this%yfluxdx(i,j-1,k)) &
                  - this%geo_src(i,j,k) - this%src(i,j,k)
             ! time step update
!CDIR IEXPAND
             cnew(i,j,k) = UpdateTimestep_modeuler(eta,dt,cold(i,j,k),cvar(i,j,k),rhs(i,j,k))
          END DO
          ! western and eastern boundary fluxes
          ! update right hand side of boundary ODE
          rhs(Mesh%IMIN-1,j,k) = Mesh%dy * this%xfluxdy(Mesh%IMIN-1,j,k)
          rhs(Mesh%IMAX+1,j,k) = -Mesh%dy * this%xfluxdy(Mesh%IMAX,j,k)
          ! time step update of boundary fluxes
!CDIR IEXPAND
          Fluxes%bxflux(j,1,k) = UpdateTimestep_modeuler(eta,dt,Fluxes%bxfold(j,1,k), &
               Fluxes%bxflux(j,1,k),rhs(Mesh%IMIN-1,j,k))
!CDIR IEXPAND
          Fluxes%bxflux(j,2,k) = UpdateTimestep_modeuler(eta,dt,Fluxes%bxfold(j,2,k), &
               Fluxes%bxflux(j,2,k),rhs(Mesh%IMAX+1,j,k))
       END DO

       ! southern and northern boundary fluxes
!CDIR NODEP
       DO i=Mesh%IMIN,Mesh%IMAX
          ! update right hand side of boundary ODE
          rhs(i,Mesh%JMIN-1,k) = Mesh%dx * this%yfluxdx(i,Mesh%JMIN-1,k)
          rhs(i,Mesh%JMAX+1,k) = -Mesh%dx * this%yfluxdx(i,Mesh%JMAX,k)
          ! time step update of boundary fluxes
!CDIR IEXPAND
          Fluxes%byflux(i,1,k) = UpdateTimestep_modeuler(eta,dt,Fluxes%byfold(i,1,k), &
                                 Fluxes%byflux(i,1,k),rhs(i,Mesh%JMIN-1,k))
!CDIR IEXPAND
          Fluxes%byflux(i,2,k) = UpdateTimestep_modeuler(eta,dt,Fluxes%byfold(i,2,k), &
                                 Fluxes%byflux(i,2,k),rhs(i,Mesh%JMAX+1,k))
       END DO

     END DO
  END SUBROUTINE ComputeCVar_modeuler


  REAL FUNCTION CalcTimestep_modeuler(this,Mesh,Physics,Fluxes,time,dtcause) RESULT(dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    REAL               :: time
    INTEGER            :: dtcause
    !------------------------------------------------------------------------!
    REAL               :: invdt_x, invdt_y
    REAL               :: dt_cfl, dt_src
    !------------------------------------------------------------------------!    
    INTENT(IN)         :: Mesh, time
    INTENT(INOUT)      :: this,Physics,dtcause
    !------------------------------------------------------------------------!   
    ! CFL condition:
    ! maximal wave speeds in each direction
    CALL MaxWaveSpeeds(Physics,Mesh,time,this%pvar,this%amax)

    ! inverse of time step in each direction
!    IF (Mesh%INUM.GT.1) THEN
!       invdt_x = MAXVAL(this%amax(:,:,1) / Mesh%dlx(:,:))
!    ELSE
!       ! set to zero, i.e. no CFL limit in x-direction
!       invdt_x = 0.0
!    END IF
!    IF (Mesh%JNUM.GT.1) THEN
!       invdt_y = MAXVAL(this%amax(:,:,2) / Mesh%dly(:,:))
!    ELSE
!       ! set to zero, i.e. no CFL limit in y-direction
!       invdt_y = 0.0
!    END IF
  
    ! largest time step due to CFL condition
    !dt_cfl = this%cfl / MAX(invdt_x, invdt_y)
    dt_cfl = this%cfl / &
    MAXVAL(this%amax(:,:,1)/Mesh%dlx(:,:)+this%amax(:,:,2)/Mesh%dly(:,:))
    ! due to cfl = 0
    dtcause = 0
    
    ! initialize this to be sure dt_src > 0
    dt_src = dt_cfl
    CALL CalcTimestep_sources(Physics%sources,Mesh,Physics,Fluxes,time, &
         this%pvar,this%cvar,dt_src,dtcause)

    dt = MIN(dt_cfl,dt_src)
  END FUNCTION CalcTimestep_modeuler


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
