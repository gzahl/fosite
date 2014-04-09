!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: timedisc_rkfehlberg .f90                                          #
!#                                                                           #
!# Copyright (C) 2011                                                        #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
! subroutines for Runge-Kutta Fehlberg method
! Reference: G.Engeln-Müllges & F.Reutter; .....
!----------------------------------------------------------------------------!
MODULE timedisc_rkfehlberg
  
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
  CHARACTER(LEN=32), PARAMETER :: ODEsolver_name = "Runge-Kutta Fehlberg"

  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Timedisc_TYP, &
       ! methods 
       InitTimedisc_rkfehlberg, &
       CloseTimedisc_rkfehlberg, &
       SolveODE_rkfehlberg, &
       ComputeCVar_rkfehlberg, &
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

  SUBROUTINE InitTimedisc_rkfehlberg(this,Mesh,Physics,os,order,stoptime,cfl,dtlimit,maxiter, &
       tol_rel,tol_abs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: os, order,maxiter
    REAL               :: stoptime,cfl,dtlimit
    REAL               :: tol_rel
    REAL,DIMENSION(:)  :: tol_abs
    !------------------------------------------------------------------------!
    INTEGER            :: err
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
    CASE(3)
       !set number of coefficients
       this%m = 3 
       ! allocate memory 
       ALLOCATE(this%coeff(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM,this%m), &
                this%A1(this%m),this%A2(this%m),this%a(this%m),this%b(this%m,this%m), &
         STAT = err)
       IF (err.NE.0) THEN
          CALL Error(this,"timedisc_rkfehlberg", "Unable to allocate memory.")
       END IF
       !set coefficient scheme of RK-Fehlberg rkf23
       this%A1 = (/ 1.0/6.0, 2.0/3.0, 1.0/6.0 /)
       this%A2 = (/ 0.0, 1.0, 0.0 /)
       this%a  = (/ 0.0, 0.5, 1.0 /)
       this%b  = RESHAPE((/ 0.0, 0.0, 0.0, &
                            0.5, 0.0, 0.0, &
                            -1.0, 2.0, 0.0 /),(/this%m,this%m/))
    CASE(5)
       !set number of coefficients
       this%m = 6 
       ! allocate memory 
       ALLOCATE(this%coeff(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM,this%m), &
                this%A1(this%m),this%A2(this%m),this%a(this%m),this%b(this%m,this%m), &
         STAT = err)
       IF (err.NE.0) THEN
          CALL Error(this,"timedisc_rkfehlberg", "Unable to allocate memory.")
       END IF
       !set coefficient scheme of RK-Fehlberg rkf45
       this%A1 = (/ 16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0 /)
       this%A2 = (/ 25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -0.2, 0.0 /)
       this%a  = (/ 0.0, 0.25, 3.0/8.0, 12.0/13.0, 1.0, 0.5 /)
       this%b  = RESHAPE((/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                            0.25, 0.0, 0.0, 0.0, 0.0, 0.0, &
                            3.0/32.0, 9.0/32.0, 0.0, 0.0, 0.0, 0.0, &
                            1932.0/2197.0, -7200.0/2179.0, 7296.0/2197.0, 0.0, 0.0, 0.0, &
                            489.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0, 0.0, 0.0, &
                            -8.0/28.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0, 0.0/),(/this%m,this%m/))
!*****************************************************
!IMPORTANT: CHECK b(1,6) = -8/28 OR -8/27 ???????????
!*****************************************************
    CASE DEFAULT
       CALL Error(this,"timedisc_rkfehlberg","time order must be 3 or 5")
    END SELECT
    IF ((this%tol_rel.LT.0.0).OR.MINVAL(this%tol_abs(:)).LT.0.0) &
         CALL Error(this,"timedisc_rkfehlberg", &
         "error tolerance levels must be greater than 0")
    IF (tol_rel.GT.1.0) THEN
         CALL Warning(this,"timedisc_rkfehlberg", &
            "adaptive step size control disabled (tol_rel>1)")
    ELSE IF(tol_rel.GE.0.01 .AND. this%order .GE. 5) THEN
         CALL Warning(this,"timedisc_rkfehlberg", &
             "You chose a relatively high tol_rel (in comparison to order)")
    END IF
  END SUBROUTINE InitTimedisc_rkfehlberg


  SUBROUTINE SolveODE_rkfehlberg(this,Mesh,Physics,Fluxes,time,dt,maxerr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    REAL               :: time,dt,maxerr
    !------------------------------------------------------------------------!
    INTEGER            :: n,i,j,k,m
    REAL               :: t,dtnew,dtold
    REAL               :: rel_err(Physics%VNUM)
    !------------------------------------------------------------------------!    
    INTENT(IN)         :: Mesh,time
    INTENT(INOUT)      :: this,Physics,Fluxes,dt,maxerr
    !------------------------------------------------------------------------!
    t = time
    ! compute right-hand-side
    CALL ComputeRHS_rkfehlberg (this,Mesh,Physics,Fluxes,t,this%pvar,this%cvar,this%coeff(:,:,:,1))
    DO m=2,this%m
       ! time step update of cell mean values
!CDIR IEXPAND
       CALL ComputeCVar_rkfehlberg(this,Mesh,Physics,dt,m,this%coeff,this%cvar,this%ctmp)
       ! set boundary values and convert2primitive ctmp => ptmp
       CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,t+this%a(m)*dt,this%ptmp,this%ctmp)
       CALL ComputeRHS_rkfehlberg(this,Mesh,Physics,Fluxes,t+this%a(m)*dt,this%ptmp,this%ctmp,this%coeff(:,:,:,m))
    END DO
   
    !reset ctmp
    this%ctmp(:,:,:) = this%cvar(:,:,:) 
!CDIR NOVECTOR   
    DO m=1,this%m
!CDIR NOVECTOR
      DO k=1,Physics%VNUM
!CDIR COLLAPSE
        DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
          DO i=Mesh%IGMIN,Mesh%IGMAX
             ! compute two solutions with different numerical order
             this%ctmp(i,j,k) = this%ctmp(i,j,k) - dt*this%A2(m)*this%coeff(i,j,k,m)
             this%cvar(i,j,k) = this%cvar(i,j,k) - dt*this%A1(m)*this%coeff(i,j,k,m)
          END DO
        END DO
      END DO
    END DO
    ! set boundary values and convert2primitive cvar => pvar
    CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,t+dt,this%pvar,this%cvar)

    !at boundary the coeff() contains fluxes and NOT rhs! (see SubRoutine ComputeRHS_rkfehlberg )
    DO m=1,this%m
      DO k=1,Physics%VNUM
        ! western and eastern
!CDIR NODEP
        DO j=Mesh%JMIN,Mesh%JMAX
           Fluxes%bxflux(j,1,k) = Fluxes%bxflux(j,1,k) - dt*this%A1(m)*this%coeff(Mesh%IMIN-1,j,k,m)
           Fluxes%bxflux(j,2,k) = Fluxes%bxflux(j,2,k) - dt*this%A1(m)*this%coeff(Mesh%IMAX+1,j,k,m)
        END DO
        ! southern and northern
!CDIR NODEP
        DO i=Mesh%IMIN,Mesh%IMAX
          Fluxes%byflux(i,1,k) = Fluxes%byflux(i,1,k) - dt*this%A1(m)*this%coeff(i,Mesh%JMIN-1,k,m)
          Fluxes%byflux(i,2,k) = Fluxes%byflux(i,2,k) - dt*this%A1(m)*this%coeff(i,Mesh%JMAX+1,k,m)
        END DO
      END DO
    END DO

    IF (this%tol_rel.LT.1.0) THEN
      ! maximum of truncation error
      DO k=1,Physics%VNUM
          rel_err(k) = MAXVAL(ABS(this%cvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k) &
                                 -this%ctmp(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k)) &
               / (this%tol_rel*ABS(this%cvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k)) + this%tol_abs(k)))
       END DO
      maxerr = MAXVAL(rel_err(:))

      !compute new time step
!      dtnew = 0.9*dt*(1./maxerr)**(1./(GetOrder(this)))
!CDIR IEXPAND
      dtnew = 0.9*dt*exp(-log(maxerr)/GetOrder(this))
      IF (maxerr.LT.1.0) THEN
         dt = MIN(dtnew,4.*dt)   ! not too large
      ELSE
         dt = MAX(dtnew,0.25*dt) ! not too small
      END IF
!         PRINT '(7(ES14.6))',time,dt,maxerr,rel_err(:)
    ELSE
       maxerr = 0.0
    END IF
  END SUBROUTINE SolveODE_rkfehlberg

  
  SUBROUTINE ComputeCVar_rkfehlberg(this,Mesh,Physics,dt,m,coeff,cvar,cnew)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL               :: dt
    INTEGER            :: m
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                       :: cvar,cnew
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM,this%m) &
                       :: coeff
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,k,mm
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this,Mesh,Physics,dt,m,coeff
    INTENT(INOUT)      :: cvar,cnew
    !------------------------------------------------------------------------!
    cnew(:,:,:) = cvar(:,:,:)
!CDIR NOVECTOR   
    DO mm=1,m-1
!CDIR NOVECTOR
       DO k=1,Physics%VNUM
!CDIR OUTERUNTOLL=8
         DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
            DO i=Mesh%IMIN,Mesh%IMAX
               cnew(i,j,k) = cnew(i,j,k) - dt*this%b(mm,m)*coeff(i,j,k,mm)
            END DO
         END DO
       END DO
    END DO
  END SUBROUTINE ComputeCVar_rkfehlberg


  SUBROUTINE ComputeRHS_rkfehlberg(this,Mesh,Physics,Fluxes,time,pvar,cvar,rhs)
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
!CDIR IEXPAND
    IF (GetType(Mesh%geometry).NE.CARTESIAN) &
       CALL GeometricalSources(Physics,Mesh,Fluxes,pvar,cvar,this%geo_src)

    ! get source terms due to external forces if present
    IF (ASSOCIATED(Physics%sources)) &
       CALL ExternalSources(Physics%sources,Mesh,Fluxes,Physics, &
            time,pvar,cvar,this%src)

    DO k=1,Physics%VNUM
       ! compute flux differences
       ! x-direction
!CDIR OUTERUNROLL=8
       DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
          DO i=Mesh%IMIN,Mesh%IMAX
             ! temporary use rhs for flux difference in x-direction
             rhs(i,j,k) = Mesh%dydV(i,j)*( &
                  this%xflux(i,j,k) - this%xflux(i-1,j,k))
          END DO
       END DO

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
  END SUBROUTINE ComputeRHS_rkfehlberg


  SUBROUTINE CloseTimedisc_rkfehlberg(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP)   :: this
    !------------------------------------------------------------------------!
       DEALLOCATE(this%coeff,this%A1,this%A2,this%a,this%b)
  END SUBROUTINE CloseTimedisc_rkfehlberg

END MODULE timedisc_rkfehlberg
