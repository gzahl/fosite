!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: timedisc_dumka.f90                                                #
!#                                                                           #
!# Copyright (C) 2012                                                        #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
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
!> \author Manuel Jung
!!
!! \brief DUMKA3 - a Chebyshev explicit Runge-Kutta method
!!
!! A Chebyshev explicit Runge-Kutta method of order three for
!! mildly stiff ODEs arising from semidiscretization of parabolic PDE
!! Reference: Alexei A. Medovikiv (1998) "High Order Explicit Methods For
!! Parabolic Equations"
!!
!! \extends timedisc_common
!! \ingroup timedisc
!----------------------------------------------------------------------------!

MODULE timedisc_dumka
  USE timedisc_common
  USE mesh_generic
  USE fluxes_generic
  USE boundary_generic
  USE physics_generic, GeometricalSources_Physics => GeometricalSources, &
       ExternalSources_Physics => ExternalSources
  USE timedisc_rkfehlberg, ONLY: ComputeRHS_rkfehlberg
  USE timedisc_modeuler, ONLY: CalcTimestep_modeuler
  USE dumka_constants
  USE sources_generic
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: ODEsolver_name = "Dumka 3"
  INTEGER, PARAMETER           :: order = 3
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Timedisc_TYP, &
       ! methods 
       InitTimedisc_dumka, &
       CloseTimedisc_dumka, &
       SolveODE_dumka, &
       CalcTimestep_dumka, &
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

  SUBROUTINE InitTimedisc_dumka(this,Mesh,Physics,config)
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
    CALL RequireKey(config, "order", 3)
    CALL GetAttr(config, "order", this%order)
    order = this%order

    ! time step maximum in units of [CFL timestep] 
    CALL RequireKey(config, "dtmax", 5.0)
    CALL GetAttr(config, "dtmax", this%dtmax)

    CALL GetAttr(config, "method", method)
    CALL InitTimedisc(this,method,ODEsolver_name)


    IF ((this%tol_rel.LT.0.0).OR.MINVAL(this%tol_abs(:)).LT.0.0) &
         CALL Error(this,"InitTimedisc_dumka", &
         "error tolerance levels must be greater than 0")
    
    this%H_N = 0.
    this%ERR_N = 0.

    ALLOCATE(this%coeff(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM,order+1),&
         STAT = err)
    IF (err.NE.0) THEN
       CALL Error(this,"timedisc_dumka", "Unable to allocate memory.")
    END IF


  END SUBROUTINE InitTimedisc_dumka


  SUBROUTINE SolveODE_dumka(this,Mesh,Physics,Fluxes,time,dt,maxerr)
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
    INTEGER, DIMENSION(12)&
                       :: iwork
    REAL               :: C_2,C_3,C_4,A_21,A_31,A_41,A_32,A_42,A_43,T,T2,T3,&
                          T4,R,dtnew
#ifdef PARALLEL
    INTEGER            :: ierror
    REAL               :: maxerr_all,dt_all
#endif
    REAL, DIMENSION(order) &
                       :: A, B, C
    INTEGER            :: K, N_POL, N_POL_DEGREE
    INTEGER            :: m, i, j,v
    !------------------------------------------------------------------------!    
    INTENT(IN)         :: Mesh,time
    INTENT(INOUT)      :: this,Physics,Fluxes,dt
    INTENT(OUT)        :: maxerr
    !------------------------------------------------------------------------!
    
    T = time

    IF(this%n_adj.EQ.0) THEN
      !CALL FUN(N,T,Y,Z)
      !CALL ComputeRHS_dumka(this,Mesh,Physics,Fluxes,T,Y,Z)
      ! at boundary the this%rhs contains fluxes and NOT rhs! (see SubRoutine ComputeRHS_rkfehlberg )
      CALL ComputeRHS_rkfehlberg(this,Mesh,Physics,Fluxes,T,dt,this%pvar,this%cvar,this%coeff(:,:,:,1))
    END IF
    

    N_POL_DEGREE=N_DEG(this%degree)
  
    N_POL=0
    DO K=INDEX_FIRST(this%degree),INDEX_LAST(this%degree)
      IF(K.EQ.INDEX_FIRST(this%degree)) THEN
        this%rhsold = this%coeff(:,:,:,1)
      ENDIF
      
      N_POL = N_POL + 3
      A_21  = dt * COEF(K,1)
      !C_2   = A_21
      A_31  = dt * COEF(K,2)
      A_32  = dt * COEF(K,3)
      !C_3   = A_31 + A_32
      A_41  = dt * COEF(K,4)
      A_42  = dt * COEF(K,5)
      A_43  = dt * COEF(K,6)
      !C_4   = A_41 + A_42 + A_43

      C(:) = (/ A_21, A_31 + A_32, A_41 + A_42 + A_43 /)
      A(:) = (/ A_21, A_32, A_43 /)

      DO m=1,order
!CDIR IEXPAND
        CALL ComputeCVar_dumka(this,Mesh,Physics,A(m),this%cvar,this%coeff(:,:,:,m))
!CDIR IEXPAND
        CALL ComputeFluxes_dumka(this,Mesh,Physics,A(m),Fluxes%bxflux,Fluxes%byflux,&
          this%coeff(:,:,:,m))

        T = T + C(m)
        CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,T,this%ptmp,this%cvar)
        CALL ComputeRHS_rkfehlberg(this,Mesh,Physics,Fluxes,T,dt,this%ptmp,this%cvar,this%coeff(:,:,:,m+1))
      END DO

      B(:) = (/ dt*(COEF(K,2)-COEF(K,1)), 0.5 * C(2), 0.5*(C(3)-C(2)) /)
      this%coeff(:,:,:,2) = B(2) * this%coeff(:,:,:,2) &
                          - B(3) * this%coeff(:,:,:,1)
      IF(N_POL.EQ.N_DEG(this%degree)) THEN
        CALL ComputeCVar_dumka(this,Mesh,Physics,B(1),this%cvar,this%coeff(:,:,:,1))
        CALL ComputeCVar_dumka(this,Mesh,Physics,-B(3),this%coeff(:,:,:,2),this%coeff(:,:,:,3))
      ENDIF

      this%coeff(:,:,:,1) = this%coeff(:,:,:,4)

    END DO

!CDIR NOVECTOR
    DO v=1,Physics%VNUM
!CDIR OUTERUNROLL=8
      DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
        DO i=Mesh%IMIN,Mesh%IMAX
          this%coeff(i,j,v,2) = (B(2) * this%coeff(i,j,v,1) - this%coeff(i,j,v,2))&
            /(this%tol_rel * MAX(ABS(this%cvar(i,j,v)),ABS(this%cold(i,j,v)))&
            +this%tol_abs(v))
        END DO
      END DO
    END DO

    CALL ST(this,Mesh,Physics,Fluxes,T,this%cvar,this%coeff(:,:,:,2),dt,N_POL_DEGREE,maxerr)

#ifdef PARALLEL
       CALL MPI_Allreduce(maxerr,maxerr_all,1,DEFAULT_MPI_REAL,MPI_MAX,&
            Mesh%comm_cart,ierror)
       maxerr = maxerr_all
#endif
    
    IF(maxerr.GE.1.0) THEN
      this%coeff(:,:,:,1) = this%rhsold
    END IF
    
    dt = CheckStability_dumka(this,Mesh,Physics,Fluxes,T,dt,this%dtcause)

#ifdef PARALLEL
       CALL MPI_Allreduce(dt,dt_all,1,DEFAULT_MPI_REAL,MPI_MIN,&
            Mesh%comm_cart,ierror)
       dt = dt_all
#endif

    ! set boundary values and convert2primitive cvar => pvar
    CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,T,this%pvar,this%cvar)

  END SUBROUTINE SolveODE_dumka

  SUBROUTINE ComputeCVar_dumka(this,Mesh,Physics,a,cvar,rhs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL               :: a
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                       :: cvar,rhs
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this,Mesh,Physics,a,rhs
    INTENT(INOUT)      :: cvar
    !------------------------------------------------------------------------!
!CDIR NOVECTOR
    DO k=1,Physics%VNUM
!CDIR OUTERUNTOLL=8
      DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
        DO i=Mesh%IMIN,Mesh%IMAX
          cvar(i,j,k) = cvar(i,j,k) - a*rhs(i,j,k)
        END DO
      END DO
    END DO
  END SUBROUTINE ComputeCVar_dumka

  SUBROUTINE ComputeFluxes_dumka(this,Mesh,Physics,a,bxflux,byflux,rhs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL               :: a
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                       :: rhs
    REAL,DIMENSION(Mesh%JMIN:Mesh%JMAX,2,Physics%VNUM) &
                       :: bxflux
    REAL,DIMENSION(Mesh%IMIN:Mesh%IMAX,2,Physics%VNUM) &
                       :: byflux
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this,Mesh,Physics,a,rhs
    INTENT(INOUT)      :: bxflux, byflux
    !------------------------------------------------------------------------!
!CDIR NOVECTOR
    DO k=1,Physics%VNUM
!CDIR NODEP
      DO j=Mesh%JMIN,Mesh%JMAX
        bxflux(j,1,k) = bxflux(j,1,k) - a * rhs(Mesh%IMIN-1,j,k)
        bxflux(j,2,k) = bxflux(j,2,k) - a * rhs(Mesh%IMAX+1,j,k)
      END DO
!CDIR NODEP
      DO i=Mesh%IMIN,Mesh%IMAX
        byflux(i,1,k) = byflux(i,1,k) - a * rhs(i,Mesh%JMIN-1,k)
        byflux(i,2,k) = byflux(i,2,k) - a * rhs(i,Mesh%JMAX+1,k)
      END DO
    END DO
  END SUBROUTINE ComputeFluxes_dumka

! returns maxerr and calculates new timestep H_NEW
  SUBROUTINE ST(this,Mesh,Physics,Fluxes,T,Y,Z1,H_NEW,N_POL_DEGREE,maxerr)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    INTEGER :: N_POL_DEGREE,ISTEP_ID
    REAL :: T,H_NEW,maxerr
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                       :: Y,Z1
    !----------------------------------------------------------------------!
    REAL :: EPS1, FRAC, FRAC2, FRACMAX, FRACMIN, H_OLD, dteuler, EPS
    INTEGER :: N, I, J, M
    INTEGER :: dtcause
    !----------------------------------------------------------------------!
    INTENT(IN) :: T,Y,Z1,N_POL_DEGREE, Mesh
    INTENT(INOUT) :: this, Physics, Fluxes, H_NEW, maxerr
    !----------------------------------------------------------------------!
    N = (Mesh%IMAX-Mesh%IMIN+1)&
       *(Mesh%JMAX-Mesh%JMIN+1)&
       *Physics%VNUM
    
! FIXME: dtcause is a dummy
    dteuler = CalcTimestep_modeuler(this,Mesh,Physics,Fluxes,T,dtcause)

!    EPS1=0.0
!    DO M=1,Physics%VNUM
!      DO J=Mesh%JMIN,Mesh%JMAX
!        DO I=Mesh%IMIN,Mesh%IMAX
!          EPS1 = EPS1 + Z1(I,J,M) * Z1(I,J,M)
!        END DO
!      END DO
!    END DO
!    print *,"EPS1: ", EPS1
!    EPS=SQRT( EPS1/(1.*N) )
    EPS = MAXVAL(Z1(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:))
    FRACMIN=0.1
    !IF(EPS.EQ.0.E0) &
    !  EPS=1.E-14       
    !EPS = MAX(EPS,TINY(EPS))
    !IF(EPS.LT.1.E-14) &
    !  print *,"Small EPS: ", EPS
    !print *,"EPS: ", EPS
    EPS = MAX(EPS,1.E-14)
    FRAC=(1./EPS)**(1./3.)
    IF(EPS.LE.1.E0) THEN
      ! step accepted
      IF ((this%ERR_N.GT.0.).AND.&
          (this%H_N.GT.0.)) THEN
        FRAC2=this%ERR_N**(1./3.) * FRAC * FRAC * (H_NEW/this%H_N)
        FRAC=MIN(FRAC,FRAC2)
      END IF
      FRACMAX=2.E0
      FRAC=MIN(FRACMAX,MAX(FRACMIN,0.8E0*FRAC))
      this%H_N = H_NEW
      H_OLD = H_NEW
      H_NEW=FRAC*H_NEW
      maxerr=EPS
      this%ERR_N = EPS
!        WRITE(*,102) EPS, T
        !WRITE(*,101) N_POL_DEGREE, H_NEW, H_OLD      
    ELSE
      ! step rejected
      FRACMAX=1.E0
      FRAC=0.8E0*MIN(FRACMAX,MAX(FRACMIN,0.8E0*FRAC))
      H_OLD = H_NEW
      IF(this%n_adj.EQ.1) THEN
        H_NEW=FRACMIN*H_NEW
      ELSE
        H_NEW=FRAC*H_NEW
      ENDIF
      maxerr=EPS
!      WRITE(*,100) EPS, T
      !WRITE(*,101) N_POL_DEGREE, H_NEW, H_OLD 
    END IF

!    WRITE(*,101) N_POL_DEGREE, H_NEW/dteuler, H_OLD/dteuler

100   FORMAT('STEP REJECTED, BECAUSE ERROR=',1E18.10,&
             ' TIME=',1E18.10)
102   FORMAT('STEP ACCEPTED, BECAUSE ERROR=',1E18.10,&
             ' TIME=',1E18.10)
101   FORMAT('  POLYNOMIAL DEGREE =', I3,&
             ' NEW STEP =',1E18.10,' OLD STEP',1E18.10)
    RETURN
    END SUBROUTINE ST

!  SUBROUTINE ComputeRHS_dumka(this,Mesh,Physics,Fluxes,time,cvar,rhs)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Timedisc_TYP) :: this
!    TYPE(Mesh_TYP)     :: Mesh
!    TYPE(Physics_TYP)  :: Physics
!    TYPE(Fluxes_TYP)   :: Fluxes
!    REAL               :: time
!    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
!                       :: cvar,rhs
!    !------------------------------------------------------------------------!
!    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
!                       :: pvar
!    !------------------------------------------------------------------------!
!    INTENT(IN)         :: Mesh,time
!    INTENT(INOUT)      :: this,Physics,Fluxes,cvar
!    INTENT(OUT)        :: rhs
!    !------------------------------------------------------------------------!
!    ! set boundary values and convert2primitive ctmp => ptmp
!  END SUBROUTINE ComputeRHS_dumka 


  REAL FUNCTION CheckStability_dumka(this,Mesh,Physics,Fluxes,time,dtold,dtcause) RESULT(dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    REAL               :: time, dtold
    INTEGER            :: dtcause
    !------------------------------------------------------------------------!    
    REAL               :: dteuler
    !------------------------------------------------------------------------!    
    INTENT(IN)         :: Mesh, time, dtold
    INTENT(INOUT)      :: this,Physics,dtcause
    !------------------------------------------------------------------------!   
!FIXME: that's not right! see CalcTimestep_dumka "double call" of this subr.
    dt = dtold
    dteuler = CalcTimestep_modeuler(this,Mesh,Physics,Fluxes,time,dtcause)

    ! Choose the right polynom degree to stay in the stable zone and
    ! adjust timestep if necessary
    this%degree=1
    DO WHILE((this%degree.NE.N_P).AND.((2.*dt/dteuler).GT.STAB_REG(this%degree)))
      this%degree = this%degree + 1
    END DO
  
    IF(this%degree.GT.1) THEN
      IF(((STAB_REG(this%degree-1)/(N_DEG(this%degree-1)*1.))*N_DEG(this%degree)*1.)&
         .GT.(2.*dt/dteuler))THEN
        this%degree = this%degree - 1
      END IF
    END IF
  

    dt = MIN((STAB_REG(this%degree)*dteuler*0.5),dt)

    ! Do not use a bigger timestep than the user defined this%dtmax * dteuler
    dt = MIN(dt/dteuler,this%dtmax) * dteuler

!    print *,"dt/dteuler: ", dt/dteuler
!    WRITE(*,300)  N_POL_DEGREE,H/(N_POL_DEGREE*COU),T
!300 FORMAT (1X,'SPEED UP=',I3,' TIMES',' H/COU=',F14.7,'  TIME',F14.7)
!    print *,"stability: H/COU < stab_reg: ", H/COU, " < ", STAB_REG(INDEX)
 
  END FUNCTION CheckStability_dumka


  REAL FUNCTION CalcTimestep_dumka(this,Mesh,Physics,Fluxes,time,dtcause) RESULT(dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    REAL               :: time
    !------------------------------------------------------------------------!    
    REAL               :: dteuler
    INTEGER            :: dtcause
    !------------------------------------------------------------------------!    
    INTENT(IN)         :: Mesh, time
    INTENT(INOUT)      :: this,Physics,dtcause
    !------------------------------------------------------------------------!   

    IF(this%dtold.EQ.this%stoptime) THEN
      ! Use Euler timestep as initial guess
      dt = CalcTimestep_modeuler(this,Mesh,Physics,Fluxes,time,dtcause)
      dt = CheckStability_dumka(this,Mesh,Physics,Fluxes,time,dt,dtcause) 
    ELSE
      ! Use guess from the last timestep
      dt = this%dtold
    END IF

  END FUNCTION CalcTimestep_dumka

  SUBROUTINE CloseTimedisc_dumka(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP)   :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%coeff)
  END SUBROUTINE CloseTimedisc_dumka

END MODULE timedisc_dumka
