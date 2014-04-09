!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_forcing.f90                                               #
!#                                                                           #
!# Copyright (C) 2012                                                        #
!# Bjoern Sperling <sperling@astrophysik.uni-kiel.de>                        #
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
!! \brief forcing term
!!
!! \extends sources_c_accel
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_forcing
  USE common_types, ONLY : Common_TYP, InitCommon
  USE sources_c_accel
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE common_dict
  USE functions
#ifdef HAVE_FFTW
  USE fftw
#endif
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  REAL, PARAMETER :: TINY = 1.0E-30              ! to avoid division by 0    !
  CHARACTER(LEN=32), PARAMETER  :: source_name = "forcing term"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! constants
       ! methods
       InitSources_forcing, &
       InfoSources_forcing, &
       CalcTimestep_forcing, &
       ExternalSources_forcing, &
       CloseSources_forcing
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources_forcing(this,Mesh,Physics,Fluxes,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Dict_TYP),POINTER :: config,IO
    INTEGER           :: stype
    !------------------------------------------------------------------------!
    REAL              :: c_velo
    INTEGER           :: err,l,m,lenwrk,lensav
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Fluxes
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "stype", stype)
    IF (.NOT.Initialized(Fluxes)) &
         CALL Error(this,"InitSources_forcing","fluxes module uninitialized")

    IF (.NOT.Initialized(Physics)) &
         CALL Error(this,"InitSources_forcing","physics module uninitialized")

    CALL InitSources(this,stype,source_name)

    ! check mesh
    IF (GetType(Fluxes).NE.MIDPOINT) &
         CALL Error(this,"InitSources_forcing","only midpoint rule is currently supported")


#ifndef HAVE_FFTW
    CALL Error(this,"InitSources_forcing", &
         "Mandatory requirement fftw has not been enabled. "//&
         "Please add --with-fftw=$FFTWDIR or similar to configure call.")
#endif
#ifndef HAVE_ISO_C_BINDING
    CALL Error(this,"InitSources_forcing", &
         "No ISO_C_BINDINGs are available for this compiler, but they are "//&
         "a Mandatory requirement for fftw.")
#endif

#ifdef HAVE_FFTW
    CALL RequireKey(config, "L", 1.0) 
    CALL GetAttr(config, "L", this%L)
    IF (this%L .LT. MINVAL(Mesh%dlx(:,:))) &
      CALL Warning(this, "InitSources_forcing", "Your length scale is smaller then your grid!")
    CALL RequireKey(config, "V", 1.0) 
    CALL GetAttr(config, "V", c_velo)
    CALL RequireKey(config, "CHI", 1.0) 
    CALL GetAttr(config, "CHI", this%CHI)
    IF (this%CHI .GT. 1.0 .OR. this%CHI .LT. 0.0) &
       CALL Error(this, "InitSources_forcing", "chi must be between 0.0 and 1.0")
    this%K0 = 2*PI/this%L
    this%F0 = c_velo**2/this%L
    this%T  = this%L/c_velo    !=sqrt(L/F0)

    IF (HasKey(config, "stoptime")) THEN 
       CALL GetAttr(config, "stoptime", this%stoptime)
    ELSE
       this%stoptime = HUGE(this%stoptime)
    END IF
       
    CALL RequireKey(config, "cvis", .1) 
    CALL GetAttr(config, "cvis", this%cvis)
        
    this%invsqrtN = 0.5/sqrt(1.0*Mesh%INUM*Mesh%JNUM)

    ALLOCATE(this%fk(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM), &
!             fftback(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM), &
             this%rand(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%DIM), &
             this%Ftemp_c(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX), &
             this%temp_c(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX), &
             this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM), &
             STAT=err)
    IF (err.NE.0) CALL Error(this,"InitSources_forcing", "Unable to allocate memory.")
    this%fk(:,:,:) = 0.0
    this%accel(:,:,:) = 0.0

    this%Ftemp_c = 0.0
    this%temp_c = 0.0

    !  Make a plan for the backward FFT, and recover the original data.
    ! Important reversed array index!!!! JNUM,INUM!
    this%plan_r2r = fftw_plan_r2r_2d(Mesh%JNUM, Mesh%INUM, this%Ftemp_c, this%temp_c, &
        FFTW_REDFT01,FFTW_REDFT01,FFTW_MEASURE)
!    testplan = fftw_plan_r2r_2d(Mesh%JNUM, Mesh%INUM, this%temp_c, this%Ftemp_c, &
!        FFTW_REDFT10,FFTW_REDFT10,FFTW_MEASURE)
   
    CALL SetOutput(this,Mesh,Physics,config,IO)
#endif
  END SUBROUTINE InitSources_forcing

 SUBROUTINE InfoSources_forcing(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: L_str,chi_str,V_str,T_str
    !------------------------------------------------------------------------!
#ifdef HAVE_FFTW
    WRITE (L_str,'(ES9.2)') this%L
    WRITE (chi_str,'(ES9.2)') this%chi
    WRITE (V_str,'(ES9.2)') this%L/this%T
    WRITE (T_str,'(ES9.2)') this%T
       CALL Info(this,"            length scale :    " // TRIM(L_str) // &
           ACHAR(10)//"            timescale:        " // TRIM(T_str) //&
           ACHAR(10)//"            spectral weight:  " // TRIM(chi_str) //&
           ACHAR(10)//"            chara. velocity:  " //TRIM(V_str))
#endif
  END SUBROUTINE InfoSources_forcing

  PURE SUBROUTINE CalcTimestep_forcing(this,Mesh,Physics,pvar,cvar,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar,cvar
    REAL              :: dt, dlx_min, dly_min
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,pvar,cvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: dt
    !------------------------------------------------------------------------!
#ifdef HAVE_FFTW
    ! timestep control
    ! x-direction
    IF (Mesh%INUM.GT.1) THEN
       dlx_min = MINVAL(Mesh%dlx(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX))
    ELSE
       ! set to zero, i.e. no limit in x-direction due to pointmass
       dlx_min = HUGE(dlx_min)
    END IF
    ! y-direction
    IF (Mesh%JNUM.GT.1) THEN
       dly_min = MINVAL(Mesh%dly(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX))
    ELSE
       ! set to zero, i.e. no limit in x-direction due to pointmass
       dly_min = HUGE(dly_min)
    END IF

    dt = this%cvis * this%T/this%L*MIN(dlx_min,dly_min)
#endif    
  END SUBROUTINE CalcTimestep_forcing

 SUBROUTINE SetOutput(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    TYPE(Physics_TYP)    :: Physics
    TYPE(Dict_TYP),POINTER  :: config,IO
    !------------------------------------------------------------------------!
    INTEGER              :: valwrite
    !------------------------------------------------------------------------!
    INTENT(IN)           :: Mesh,Physics
    INTENT(INOUT)        :: this
    !------------------------------------------------------------------------!
#ifdef HAVE_FFTW 
    valwrite = 0
    IF (HasKey(config, "output/accel")) CALL GetAttr(config, "output/accel", valwrite)
    IF (valwrite .EQ. 1) THEN
      
       CALL AddField(IO, &
               "accel", &
               this%accel, &
               Dict("name" / "accel"))
    END IF
    valwrite = 0
    IF (HasKey(config, "output/fk")) CALL GetAttr(config, "output/fk", valwrite)
    IF (valwrite .EQ. 1) THEN
      
       CALL AddField(IO, &
               "fk", &
               this%fk, &
               Dict("name" / "fk"))
    END IF

    valwrite = 0
    IF (HasKey(config, "output/f")) CALL GetAttr(config, "output/f", valwrite)
    IF (valwrite .EQ. 1) THEN
      
       CALL AddField(IO, &
               "f", &
               this%temp_c, &
               Dict("name" / "f"))
    END IF
   
! CALL AddField(IO, &
!               "fftback", &
!               fftback, &
!               Dict("name" / "fftbacks"))
#endif
  END SUBROUTINE SetOutput
 
 
  SUBROUTINE ExternalSources_forcing(this,Mesh,Physics,time,dt,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time,dt
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: cvar,pvar,sterm 
    !------------------------------------------------------------------------!
    INTEGER           :: err,i,j,k
    REAL,DIMENSION(Physics%DIM,Physics%DIM) :: B
    REAL,DIMENSION(Physics%DIM) :: Wb
    REAL              :: ki,kj,kk,k2
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,time,dt,pvar,cvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    sterm(:,:,:) = 0.0
#ifdef HAVE_FFTW
    ! after stoptime => decaying turbluence
    IF (time .GT. this%stoptime)  THEN
       this%accel(:,:,:) = 0.0
       return
    END IF

    ! uniform random numbers
    CALL RANDOM_NUMBER(this%rand(:,:,1:2))
    ! normal distributed random numbers
    CALL NormalRand(this%rand(:,:,1),this%rand(:,:,2), &
                    this%rand(:,:,1),this%rand(:,:,2))

    IF (Physics%DIM .GE. 3) THEN
       !2 wird überschrieben... macht aber nix
       CALL RANDOM_NUMBER(this%rand(:,:,2:3))
       CALL NormalRand(this%rand(:,:,2),this%rand(:,:,3), &
                       this%rand(:,:,2),this%rand(:,:,3))
    END IF

    kk = 0.0 !1 or 0????
    DO i=Mesh%IMIN,Mesh%IMAX
      DO j=Mesh%JMIN,Mesh%JMAX
!TODO: hier evtl. MATMUL auf gesamten ARRAY ausführen (effizienter?)
!FIXME: k richtig berechnen
        !k2 = 1.0*i**2+j**2+k**2
        ki = 2.0*PI/(Mesh%xmax-Mesh%xmin)*i
        kj = 2.0*PI/(Mesh%ymax-Mesh%ymin)*j
        k2 = ki**2+kj**2+kk**2
        B(1,1) = this%chi + (1.0-2.0*this%chi)*ki**2/k2
        B(1,2) =            (1.0-2.0*this%chi)*ki*kj/k2
        B(2,2) = this%chi + (1.0-2.0*this%chi)*kj**2/k2
        B(2,1) =            (1.0-2.0*this%chi)*kj*ki/k2
        IF (Physics%DIM .GE. 3) THEN
          B(1,3) =            (1.0-2.0*this%chi)*ki*kk/k2    
          B(2,3) =            (1.0-2.0*this%chi)*kj*kk/k2    
          B(3,1) =            (1.0-2.0*this%chi)*kk*ki/k2    
          B(3,2) =            (1.0-2.0*this%chi)*kk*kj/k2    
          B(3,3) = this%chi + (1.0-2.0*this%chi)*kk**2/k2 
        END IF        
        Wb(:)=MATMUL(B(:,:),this%rand(i,j,:))
        this%fk(i,j,:) = this%fk(i,j,:)*(1.0 - dt/this%T) &
                + this%F0*ABS(sigma(Mesh,this%K0,i,j))*SQRT(2.0*dt/this%T)*Wb(:)
      END DO
    END DO

    DO i=1,Physics%DIM
      this%Ftemp_c(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) = &
           this%fk(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,i)
      CALL fftw_execute_r2r(this%plan_r2r,this%Ftemp_c, this%temp_c)
      this%temp_c = this%temp_c*this%invsqrtN

!this%Ftemp_c = 0.0
!CALL fftw_execute_r2r(testplan, this%temp_c, this%Ftemp_c)
!fftback(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,i) = &
!   this%invsqrtN*this%Ftemp_c(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)


       this%accel(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,i) = &
          this%temp_c(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) &
          /pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY)
    END DO


    ! compute source terms due to constant acceleration
    CALL ExternalSources(Physics,Mesh,this%accel,pvar,cvar,sterm)
#endif
  END SUBROUTINE ExternalSources_forcing
  
  FUNCTION sigma(Mesh,k0,i,j) 
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP):: Mesh
    INTEGER       :: i,j
    !------------------------------------------------------------------------!
    REAL          :: k0,k,ki,kj,a,sigma
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,k0,i,j
    !------------------------------------------------------------------------!
!FIXME: Hier korrektes k berechnen
    !k = SQRT(1.0*i**2+j**2)
    ki = 2.0*PI/(Mesh%xmax-Mesh%xmin)*i
    kj = 2.0*PI/(Mesh%ymax-Mesh%ymin)*j
    k = SQRT(ki**2+kj**2)

    a = -0.75/K0**3
    IF (k .GE. 2.0*K0) THEN 
!    IF (k .GE. 1.5*K0 .OR. k .LE. 0.5*K0) THEN 
      sigma = 0.0
    ELSEIF (i .GE. Mesh%INUM/2) THEN
      sigma = 0.0
    ELSEIF (j .GE. Mesh%JNUM/2) THEN
      sigma = 0.0
    ELSE
      sigma = a*k*(k-2.0*K0)
    END IF
  END FUNCTION sigma

 
  SUBROUTINE CloseSources_forcing(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
#ifdef HAVE_FFTW
     CALL fftw_destroy_plan(this%plan_r2r)
     DEALLOCATE(this%fk,this%rand,this%accel,this%Ftemp_c,this%temp_c)
#endif
    CALL CloseSources(this)
  END SUBROUTINE CloseSources_forcing
 

END MODULE sources_forcing
