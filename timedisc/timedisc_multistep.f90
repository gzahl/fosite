!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: timedisc_multistep.f90                                            #
!#                                                                           #
!# Copyright (C) 2013                                                        #
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
!> \author Björn Sperling
!!
!! \brief multistep method (adam method ) with variable step size and order
!!
!! see e.g. Hairer, Norsett, Wanner, "Solving ODE I", 3rd. 2008 Springer-Verlag
!!          ISBN: 978-3-540-56670-0
!!
!! \extends timedisc_common
!! \ingroup timedisc
!----------------------------------------------------------------------------!
MODULE timedisc_multistep
  USE timedisc_common
  USE mesh_generic
  USE fluxes_generic
  USE boundary_generic
  USE physics_generic, GeometricalSources_Physics => GeometricalSources, &
       ExternalSources_Physics => ExternalSources
  USE timedisc_modeuler, ONLY: CalcTimestep_multistep => CalcTimestep_modeuler
  USE sources_generic
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: ODEsolver_name = "multistep"
  INTEGER, PARAMETER :: DTCAUSE_SMALLERR = -2  ! smallest ts due to err      !
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Timedisc_TYP, &
       ! methods 
       InitTimedisc_multistep, &
       CloseTimedisc_multistep, &
       SolveODE_multistep, &
       CalcTimestep_multistep, &
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

  SUBROUTINE InitTimedisc_multistep(this,Mesh,Physics,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Dict_TYP), POINTER &
                       :: config
    !------------------------------------------------------------------------!
    INTEGER            :: err,order,method
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,Physics
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    ! set default order 
    CALL RequireKey(config, "order", 12)
    CALL GetAttr(config, "order", this%order)
    order=this%order

    ! enable predictor-corrector approach ("implicit")
    ! implicit is strongly recommended!!!
    CALL RequireKey(config, "implicit", 1) ! 0 = disable, 1 = enable
    CALL GetAttr(config, "implicit", this%pc)

    CALL GetAttr(config, "method", method)
    CALL InitTimedisc(this,method,ODEsolver_name)

    IF (order .LE. 0) CALL Error(this, "InitTimedisc_multistep", "max(order) < 1")
    this%m = 1
    ! allocate memory for data structures only needed in multistep 
    ALLOCATE(this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM,0:order), &
         this%newphi_s(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM,0:order), &
         this%oldphi_s(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM,0:order), &
         this%gamma(0:order), &
         this%c(0:order,1:order+1), &
         STAT = err)
    IF (err.NE.0) THEN
       CALL Error(this,"InitTimedisc_multistep", "Unable to allocate memory.")
    END IF

    this%gamma = 0.
    this%c = 0.
    this%phi = 0.
    this%oldphi_s = 0.
    this%newphi_s = 0.

    ! Init ring structure
    CALL InitRing(this%ts)


    IF (this%pc.EQ.0)&
         CALL Warning(this,"InitTimedisc_multistep", &
         "Explicit multistep is experimental. Use 'implicit' instead.")

    IF ((this%tol_rel.LT.0.0).OR.MINVAL(this%tol_abs(:)).LT.0.0) &
         CALL Error(this,"InitTimedisc_multistep", &
         "error tolerance levels must be greater than 0")
  END SUBROUTINE InitTimedisc_multistep


  SUBROUTINE SolveODE_multistep(this,Mesh,Physics,Fluxes,time,dt,maxerr)
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
    TYPE(Fluxes_TYP)   :: Fluxes
    REAL               :: time,dt,maxerr
    !------------------------------------------------------------------------!
    REAL,DIMENSION(:,:,:,:),POINTER :: ptemp
    REAL               :: rel_err(Physics%VNUM),err(-1:2),dtnew
    INTEGER            :: k,j,o,nextorder,maxorder
#ifdef PARALLEL
    INTEGER            :: ierror
    REAL               :: maxerr_all,err_all(-1:2)
#endif
    !------------------------------------------------------------------------!    
    INTENT(IN)         :: Mesh,time
    INTENT(INOUT)      :: this,Physics,Fluxes,dt,maxerr
    !------------------------------------------------------------------------!
    maxorder = GetOrder(this)    
    o = this%m

    ! compute Phi_i(time), 0<= i <=o-1
    CALL CalcPhi(this,Mesh,Physics,Fluxes,time,dt,o,this%oldphi_s,this%phi)
    ! compute PhiS_i, 0<= i <=o-1
    CALL CalcPhi_s(this,Mesh,Physics,time,dt,o,this%phi,this%gamma,this%newphi_s)
    ! compute c_ij, 0<= i <=o , 1<= j <= o-1
    CALL CalcC(this,time,dt,o,this%c)

    ! starting process
IF (o .le. 1) THEN
     CALL SolveODE_RK2(this,Mesh,Physics,Fluxes,time,dt,maxerr)
     nextorder = o+1
#ifdef PARALLEL
     CALL MPI_Allreduce(maxerr,maxerr_all,1,DEFAULT_MPI_REAL,MPI_MAX,&
            Mesh%comm_cart,ierror)
     maxerr = maxerr_all
#endif

 ELSE

    !explicit adam step (P)
    this%ctmp(:,:,:) = 0.0
    ! Note: max(0,o-2) special treatment for o==1 
    DO j=0,max(0,o-2)
      ! temp
      this%ctmp(:,:,:) = this%ctmp(:,:,:) + this%c(j,1)*this%newphi_s(:,:,:,j)
    END DO
    this%cvar = this%cvar - dt*this%ctmp
    
    ! time step update of boundary fluxes
    Fluxes%bxflux(:,1,:) = Fluxes%bxfold(:,1,:) - dt*this%ctmp(Mesh%IMIN-1,:,:)
    Fluxes%bxflux(:,2,:) = Fluxes%bxfold(:,2,:) - dt*this%ctmp(Mesh%IMAX+1,:,:)
    Fluxes%byflux(:,1,:) = Fluxes%byfold(:,1,:) - dt*this%ctmp(:,Mesh%JMIN-1,:)
    Fluxes%byflux(:,2,:) = Fluxes%byfold(:,2,:) - dt*this%ctmp(:,Mesh%JMAX+1,:)

   ! set boundary values
    CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,time+dt,this%pvar,this%cvar)
   
   IF (this%pc .eq. 1) THEN
     ! implicit adam (E)
     ! compute Phi_i(time+dt), 0<= i <=o
     CALL CalcPhi(this,Mesh,Physics,Fluxes,time+dt,dt,o+1,this%newphi_s,this%phi)
     ! (C)
     ! Note: max(1,o-1) special treatment for o==1 
     this%cvar = this%cvar-dt*this%c(max(1,o-1),1)*this%phi(:,:,:,max(1,o-1))

     ! time step update of boundary fluxes
     Fluxes%bxflux(:,1,:) = Fluxes%bxflux(:,1,:) &
         - dt*this%c(max(1,o-1),1)*this%phi(Mesh%IMIN-1,:,:,max(1,o-1))
     Fluxes%bxflux(:,2,:) = Fluxes%bxflux(:,2,:) &
         - dt*this%c(max(1,o-1),1)*this%phi(Mesh%IMAX+1,:,:,max(1,o-1))
     Fluxes%byflux(:,1,:) = Fluxes%byflux(:,1,:) &
         - dt*this%c(max(1,o-1),1)*this%phi(:,Mesh%JMIN-1,:,max(1,o-1))
     Fluxes%byflux(:,2,:) = Fluxes%byflux(:,2,:) &
         - dt*this%c(max(1,o-1),1)*this%phi(:,Mesh%JMAX+1,:,max(1,o-1))

     ! set boundary values
     CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,time+dt,this%pvar,this%cvar)
     ! new computation of phi at next step will be done (E) => (PECE)
   END IF
 
    ! error estimate: err(0) = error of current order, err(-1) = error of higher order
    ! and err(1:2) = errors of lower orders
    err(:) = huge(err(0))
    IF (this%pc .eq. 1 .AND. o .GT. 1 .AND. this%tol_rel.LT.1.0) THEN
      ! "implicit"
      ! compute err_o-1, err_o-2
      DO j=0,min(3,o-1)
        DO k=1,Physics%VNUM
          rel_err(k) = MAXVAL(ABS(dt*(this%c(o-j,1)-this%c(o-j-1,1))* & 
            this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,o-j)) &
            / (this%tol_rel* &
          ABS(this%cvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k)) &
                 + this%tol_abs(k)))
        END DO
        err(j-1) = MAXVAL(rel_err(:))
      END DO
    ELSE IF (this%pc .eq. 0 .AND. o .GT. 1 .AND. this%tol_rel.LT.1.0) THEN
      ! "explicit"
      ! compute err_o-1, err_o-2
      DO j=0,min(3,o-1)
        DO k=1,Physics%VNUM
          rel_err(k) = MAXVAL(ABS(dt*this%c(o-j-1,1)* & 
            this%newphi_s(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,o-j-1)) &
            / (this%tol_rel* &
          ABS(this%cvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k)) &
                 + this%tol_abs(k)))
        END DO
        err(j-1) = MAXVAL(rel_err(:))
      END DO
    ELSE
      ! no error control
      err(0) = 0.0
      err(-1) = -1.0
    END IF

#ifdef PARALLEL
    CALL MPI_Allreduce(err,err_all,4,DEFAULT_MPI_REAL,MPI_MAX,&
            Mesh%comm_cart,ierror)
    err = err_all
#endif

    IF ((MAXVAL(err(1:2)) .LE. err(0))) THEN
!    IF (err(1) .LE. err(0)) THEN
       nextorder = o-1
       CALL RemoveLast(this%ts)
    ELSE IF ((err(0).LT.1.0) .AND. err(-1) .LT. err(0)) THEN 
       nextorder = o+1
    ELSE
       nextorder = o
    END IF

    ! But... 
    nextorder = min(nextorder,maxorder)

    !return value
    maxerr = err(0)
END IF

    ! otherwise restart
    IF (maxerr.LT.1.0) THEN
      IF (nextorder .GT. o) THEN
        CALL AddFirst(this%ts)
      ELSE
        CALL Rotate(this%ts)
      END IF
      CALL Set(this%ts,0,time+dt)
      this%m = nextorder

      ! switch phi_s (n -> n+1)
      ptemp => this%oldphi_s
      this%oldphi_s => this%newphi_s
      this%newphi_s => ptemp
    END IF
    dtnew = 0.9*dt*exp(-log(maxerr)/(this%m+this%pc))
    IF (maxerr.LT.1.0) THEN
       dt = MIN(dtnew,4.*dt)   ! not too large
       ! it is possible: maxerr < 1 =>  0.9*dt < dtnew => maybe dtnew < dt
       IF (dtnew .LT. dt) this%dtcause = DTCAUSE_SMALLERR
    ELSE
       dt = MAX(dtnew,0.25*dt) ! not too small
    END IF
  END SUBROUTINE SolveODE_multistep

 SUBROUTINE SolveODE_RK2(this,Mesh,Physics,Fluxes,time,dt,maxerr)
  IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    REAL               :: time,dt,maxerr
    !------------------------------------------------------------------------!
    INTEGER            :: k
    REAL               :: rel_err(Physics%VNUM)
    !------------------------------------------------------------------------!    
    INTENT(IN)         :: Mesh,time
    INTENT(INOUT)      :: this,Physics,Fluxes,dt,maxerr
    !------------------------------------------------------------------------!
    ! euler step (first order)
    this%ctmp(:,:,:) = this%cvar(:,:,:) - dt * this%phi(:,:,:,0)
    ! set boundary values and convert2primitive ctmp => ptmp
    CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,time+dt,this%ptmp,this%ctmp)


    CALL ComputeRHS(this,Mesh,Physics,Fluxes,time+dt,dt,this%ptmp,this%ctmp,this%rhs)
    ! heun method (mod. euler, second order)
    this%cvar(:,:,:) = this%cvar(:,:,:) - 0.5*dt* &
                        (this%phi(:,:,:,0) + this%rhs(:,:,:))
    CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,time+dt,this%pvar,this%cvar)

    ! time step update of boundary fluxes
    Fluxes%bxflux(:,1,:) = Fluxes%bxfold(:,1,:) - 0.5*dt* &
                        (this%phi(Mesh%IMIN-1,:,:,0) + this%rhs(Mesh%IMIN-1,:,:))
    Fluxes%bxflux(:,2,:) = Fluxes%bxfold(:,2,:) - 0.5*dt* &
                        (this%phi(Mesh%IMAX+1,:,:,0) + this%rhs(Mesh%IMAX+1,:,:))
    Fluxes%byflux(:,1,:) = Fluxes%byfold(:,1,:) - 0.5*dt* &
                        (this%phi(:,Mesh%JMIN-1,:,0) + this%rhs(:,Mesh%JMIN-1,:))
    Fluxes%byflux(:,2,:) = Fluxes%byfold(:,2,:) - 0.5*dt* &
                        (this%phi(:,Mesh%JMAX+1,:,0) + this%rhs(:,Mesh%JMAX+1,:))

   
    ! compute error
    DO k=1,Physics%VNUM
      rel_err(k) = MAXVAL(ABS(this%cvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k) &
                            -this%ctmp(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k)) &
                            / (this%tol_rel *&
                         ABS(this%cvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k)) &
                            + this%tol_abs(k)))
    END DO
    maxerr = MAXVAL(rel_err(:))
  END SUBROUTINE SolveODE_RK2


  SUBROUTINE CalcC(this,time,dt,order,c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    REAL               :: time,dt
    REAL,DIMENSION(0:this%order,1:this%order+1) &
                       :: c
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,q,order
    !------------------------------------------------------------------------!    
    INTENT(IN)         :: time,dt,order
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    ! compute c_j,q coefficients because g_j = c_j,1
    ! c_0,q = 1/q, c_1,q = 1/(q(q+1)) and 
    ! c_j,q = c_j-1,q - c_j-1,q+1 * dt/(t_n+1-t_n-j+1)
    ! see F.T. Krogh (1974)
    DO q=1,order+1
      c(0,q) = 1.0 / q
      c(1,q) = 1.0 / (q*(q+1.0))
    END DO
    ! do not interchange the loops!
    DO j=2,order
      DO q=1,order-j+1
        c(j,q) =  c(j-1,q) - c(j-1,q+1)*dt &
                             /(time+dt-Get_t(this%ts,j-1))
      END DO
    END DO
  END SUBROUTINE CalcC

  SUBROUTINE CalcPhi_s(this,Mesh,Physics,time,dt,order,phi,gamma,phi_s)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL               :: time,dt
    INTEGER            :: order
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM,0:this%order) &
                       :: phi,phi_s
    REAL,DIMENSION(0:this%order) &
                       :: gamma
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,q
    !------------------------------------------------------------------------!    
    INTENT(IN)         :: Mesh,Physics,time,dt,order,phi
    INTENT(INOUT)      :: this,gamma,phi_s
    !------------------------------------------------------------------------!
    ! compute gamma and phi_s
    phi_s(:,:,:,0) = phi(:,:,:,0)
    gamma(0) = 1.0
    DO i=1,order-1
      gamma(i) = gamma(i-1) * (time+dt-Get_t(this%ts,i-1))/&
                            (Get_t(this%ts,0)-Get_t(this%ts,i))
      phi_s(:,:,:,i) = gamma(i)*phi(:,:,:,i)
    END DO
  END SUBROUTINE CalcPhi_s

  SUBROUTINE CalcPhi(this,Mesh,Physics,Fluxes,time,dt,order,oldphi_s,phi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    REAL               :: time,dt
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM,0:this%order) &
                       :: oldphi_s,phi
    INTEGER            :: order
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,q
    !------------------------------------------------------------------------!    
    INTENT(IN)         :: Mesh,time,dt,oldphi_s,order
    INTENT(INOUT)      :: this,Physics,Fluxes,phi
    !------------------------------------------------------------------------!
    ! compute phi
    CALL ComputeRHS(this,Mesh,Physics,Fluxes,time,dt, &
                 this%pvar,this%cvar,phi(:,:,:,0))
    DO i=1,order-1
      phi(:,:,:,i) = phi(:,:,:,i-1)-oldphi_s(:,:,:,i-1)
    END DO
  END SUBROUTINE CalcPhi

  SUBROUTINE ComputeRHS(this,Mesh,Physics,Fluxes,time,dt, &
                                  pvar,cvar,rhs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    REAL               :: time,dt
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                       :: pvar,cvar,rhs
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,time,dt,pvar,cvar
    INTENT(INOUT)      :: this,Physics,Fluxes
    INTENT(OUT)        :: rhs
    !------------------------------------------------------------------------!
    ! get the numerical fluxes
    CALL CalculateFluxes(Fluxes,Mesh,Physics,pvar,cvar,this%xfluxdy,this%yfluxdx)

    ! get geometrical sources
!CDIR IEXPAND
    CALL GeometricalSources(Physics,Mesh,Fluxes,pvar,cvar,this%geo_src)

    ! get source terms due to external forces if present
    IF (ASSOCIATED(Physics%sources)) &
       CALL ExternalSources(Physics%sources,Mesh,Fluxes,Physics, &
            time,dt,pvar,cvar,this%src)

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
          END DO
          ! western and eastern boundary fluxes
          ! update right hand side of boundary ODE
          rhs(Mesh%IMIN-1,j,k) = Mesh%dy * this%xfluxdy(Mesh%IMIN-1,j,k)
          rhs(Mesh%IMAX+1,j,k) = -Mesh%dy * this%xfluxdy(Mesh%IMAX,j,k)
       END DO
       ! southern and northern boundary fluxes
       DO i=Mesh%IMIN,Mesh%IMAX
          ! update right hand side of boundary ODE
          rhs(i,Mesh%JMIN-1,k) = Mesh%dx * this%yfluxdx(i,Mesh%JMIN-1,k)
          rhs(i,Mesh%JMAX+1,k) = -Mesh%dx * this%yfluxdx(i,Mesh%JMAX,k)
       END DO
     END DO
  END SUBROUTINE ComputeRHS


  SUBROUTINE CloseTimedisc_multistep(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP)   :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%phi,this%oldphi_s,this%newphi_s,this%gamma,this%c)
  END SUBROUTINE CloseTimedisc_multistep


END MODULE timedisc_multistep
