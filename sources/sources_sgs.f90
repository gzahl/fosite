!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_sgs.f90                                                   #
!#                                                                           #
!# Copyright (C) 2010-2012                                                   #
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
!# along with Mesh program; if not, write to the Free Software               #
!# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
!#                                                                           #
!#############################################################################
!> \addtogroup sources
!! - parameters of sgs source as key-values
!! \key{model,INTEGER,type/model of used filter,1}
!! \key{output/stress,INTEGER,enable(=1) output of the stress tensor,0}
!! \key{output/diffusion,INTEGER,enable(=1) output of diffusion term,0}
!! \key{output/rhoeps,INTEGER,enable(=1) output of dissipation term,0}
!! \key{output/sigma,INTEGER,enable(=1) output of production term,0}
!! \key{output/dynvis,INTEGER,enable(=1) output of eddy viscosity (dynamic),0}
!! \key{output/kinvis,INTEGER,enable(=1) output of eddy viscosity (kinematic),0}
!----------------------------------------------------------------------------!
!> \author Björn Sperling
!!
!! \brief sgs closure model
!!
!! \extends sources_c_accel
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_sgs
  USE common_types, ONLY : Common_TYP, InitCommon
  USE sources_c_accel
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  REAL, PARAMETER :: TINY = 1.0E-30              ! to avoid division by 0    !
  CHARACTER(LEN=32), PARAMETER :: source_name = "sgs closure model"
  ! flags for viscosity model
  !W. Schmidt and C. Federrath, A&A 528, A106 (2011)
  INTEGER, PARAMETER :: SCHMIDT2011 = 1
  REAL, PARAMETER    :: Cd          = 0.65     !fixed
  REAL, PARAMETER    :: C1          = 0.102    !eddy-viscosity
  REAL, PARAMETER    :: Cdiff       = 5.0E-0
!   REAL, PARAMETER    :: C1          = 0.833  !nonlinear
!   REAL, PARAMETER    :: C1          = 0.02   !nonlinear/both  (C1 and C2)
!   REAL, PARAMETER    :: C2          = 0.07   !both (C1 and C2)
  REAL, PARAMETER    :: Ceps        = 1.5
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! constants
       ! methods
       InitSources_sgs, &
       CalcTimestep_sgs, &
       ExternalSources_sgs, &
       CloseSources_sgs
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources_sgs(this,Mesh,Physics,Fluxes,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Dict_TYP),POINTER :: config,IO
    INTEGER           :: stype
    !------------------------------------------------------------------------!
    INTEGER           :: err,m
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Fluxes
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "stype", stype)
    IF (.NOT.Initialized(Fluxes)) &
         CALL Error(this,"InitSources_sgs","fluxes module uninitialized")

    IF (.NOT.Initialized(Physics)) &
         CALL Error(this,"InitSources_sgs","physics module uninitialized")

    CALL InitSources(this,stype,source_name)

    ! check mesh
    IF (GetType(Fluxes).NE.MIDPOINT) &
         CALL Error(this,"InitSources_sgs","only midpoint rule is currently supported")

    SELECT CASE(GetType(Physics))
    CASE(EULER3D_ROTSYMSGS)
       !do nothing
    CASE(EULER3D_ROTAMTSGS,EULER2D_SGS)
       CALL Warning(this,"InitSources_sgs","this physics module is in development")
    CASE DEFAULT
       CALL Error(this,"InitSources_sgs","use physics module with sgs support")
    END SELECT

    ALLOCATE(this%btxx(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btyy(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btzz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btxy(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btxz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btyz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%Sxx(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%Syy(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%Szz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%Sxy(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%Sxz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%Syz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%delta(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%diff(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%sigma(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%rhoeps(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%dynvis(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%kinvis(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%tmp(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%tmp2(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%tmp3(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         STAT=err)
    IF (err.NE.0) CALL Error(this,"InitSources_viscosity", "Unable to allocate memory.")

    this%Sxx(:,:)=0.0
    this%Syy(:,:)=0.0
    this%Szz(:,:)=0.0
    this%Sxy(:,:)=0.0
    this%Sxz(:,:)=0.0
    this%Syz(:,:)=0.0


    m = 1
    IF (HasKey(config, "model")) CALL GetAttr(config, "model", m)
    IF (m .eq. 2) THEN
       this%delta(:,:) = (Mesh%dlx(:,:)*Mesh%dly(:,:)*2.0*PI*Mesh%bhz(:,:))**(1.0/3.0)
       CALL Warning(this,"experimental! SGS_Model", "special Filter; only working in rotsym") 
    ELSE IF (m .eq. 3) THEN
       this%delta(:,:) = (Mesh%dlx(:,:)*Mesh%dly(:,:)*&
            sqrt(Mesh%bccart(:,:,1)**2+Mesh%bccart(:,:,2)**2)/10.0)**(1.0/3.0)
       CALL Warning(this,"experimental! SGS_Model", "special Filter; only working in flat geo")
    ELSE IF (m .eq. 4) THEN
       this%delta(:,:) = (Mesh%dlx(:,:)*Mesh%dly(:,:)*Mesh%bhz(:,:))**(1.0/3.0)
       CALL Warning(this,"experimental! SGS_Model", "special Filter; only working in rotsym NO 2pi")
    ELSE
       this%delta(:,:) = sqrt(Mesh%dlx(:,:)*Mesh%dly(:,:))
    END IF
    

    CALL SetOutput(this,Mesh,Physics,config,IO)

  END SUBROUTINE InitSources_sgs

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
    valwrite = 0
    IF (HasKey(config, "output/stress")) CALL GetAttr(config, "output/stress", valwrite)
    IF (valwrite .EQ. 1) THEN
       CALL AddField(IO, &
               "SGS-stress:xx", &
               this%btxx, &
               Dict("name" / "SGS-stress:xx"))

       CALL AddField(IO, &
               "SGS-stress:xy", &
               this%btxy, &
               Dict("name" / "SGS-stress:xy"))

       CALL AddField(IO, &
               "SGS-stress:yy", &
               this%btyy, &
               Dict("name" / "SGS-stress:yy"))
     SELECT CASE(GetType(Physics))
     CASE (EULER3D_ROTSYMSGS)
       CALL AddField(IO, &
               "SGS-stress:xz", &
               this%btxz, &
               Dict("name" / "SGS-stress:xz"))

       CALL AddField(IO, &
               "SGS-stress:yz", &
               this%btyz, &
               Dict("name" / "SGS-stress:yz"))

       CALL AddField(IO, &
               "SGS-stress:zz", &
               this%btzz, &
               Dict("name" / "SGS-stress:zz"))
     END SELECT
    END IF

    valwrite = 0
    IF (HasKey(config, "output/diffusion")) CALL GetAttr(config, "output/diffusion", valwrite)
    IF (valwrite .EQ. 1) &
       CALL AddField(IO, &
               "SGS-diffusion", &
               this%diff, &
               Dict("name" / "SGS-diffusion"))

    valwrite = 0
    IF (HasKey(config, "output/rhoeps")) CALL GetAttr(config, "output/rhoeps", valwrite)
    IF (valwrite .EQ. 1) &
       CALL AddField(IO, &
               "SGS-rhoeps", &
               this%rhoeps, &
               Dict("name" / "SGS-rhoeps"))

    valwrite = 0
    IF (HasKey(config, "output/sigma")) CALL GetAttr(config, "output/sigma", valwrite)
    IF (valwrite .EQ. 1) &
       CALL AddField(IO, &
               "SGS-sigma", &
               this%sigma, &
               Dict("name" / "SGS-sigma"))

    valwrite = 0
    IF (HasKey(config, "output/dynvis")) CALL GetAttr(config, "output/dynvis", valwrite)
    IF (valwrite .EQ. 1) &
       CALL AddField(IO, &
               "SGS-dynvis", &
               this%dynvis, &
               Dict("name" / "SGS-dynvis"))

    valwrite = 0
    IF (HasKey(config, "output/kinvis")) CALL GetAttr(config, "output/kinvis", valwrite)
    IF (valwrite .EQ. 1) &
       CALL AddField(IO, &
               "SGS-kinvis", &
               this%kinvis, &
               Dict("name" / "SGS-kinvis"))

  END SUBROUTINE SetOutput
 
 
  PURE SUBROUTINE ExternalSources_sgs(this,Mesh,Physics,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: cvar,pvar,sterm 
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,cvar,pvar
    INTENT(INOUT)     :: this,Physics
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!

    DO i=Mesh%IMIN,Mesh%IMAX
       DO j=Mesh%JMIN,Mesh%JMAX
          this%dynvis(i,j) = C1SCHMIDT(this%delta(i,j),&
                             cvar(i,j,Physics%DENSITY),cvar(i,j,Physics%SGSENERGY))
       END DO
    END DO 

   ! calc SGS tensor field
    CALL CalculateSGSTensor(Physics,Mesh,this,this%dynvis,pvar)

    ! rho*epsilon (SGS energy dissipation)
    this%rhoeps(:,:) = SGSdissipation(this%delta(:,:),cvar(:,:,Physics%SGSENERGY),&
                                    cvar(:,:,Physics%DENSITY))
    ! Sigma (SGS energy production)
    this%sigma(:,:) = SGSproduction(this%btxx(:,:),this%btxy(:,:),this%btxz(:,:),&
                                   this%btyy(:,:),this%btyz(:,:),this%btzz(:,:),&
                                   this%Sxx(:,:), this%Sxy(:,:), this%Sxz(:,:),&
                                   this%Syy(:,:), this%Syz(:,:), this%Szz(:,:))

    ! D (transport term in ksgs)
    CALL SGStransport(this,Mesh,this%delta(:,:),cvar(:,:,Physics%DENSITY),&
                      cvar(:,:,Physics%SGSENERGY),this%diff(:,:)) 

    CALL SGSSources(Physics,Mesh,this,pvar,cvar,sterm)
  
  END SUBROUTINE ExternalSources_sgs

  SUBROUTINE CalcTimestep_sgs(this,Mesh,Physics,time,pvar,cvar,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar,cvar
    REAL              :: dt
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    REAL              :: invdtmax, invdt
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,time,pvar,cvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: dt
    !------------------------------------------------------------------------!
    !compute coefficient
    DO i=Mesh%IMIN,Mesh%IMAX
      DO j=Mesh%JMIN,Mesh%JMAX
        this%kinvis(i,j) = C1SCHMIDT(this%delta(i,j),&
                        cvar(i,j,Physics%DENSITY),cvar(i,j,Physics%SGSENERGY)) &
                        /cvar(i,j,Physics%DENSITY)
      END DO
    END DO

    ! x-direction
    IF (Mesh%INUM.GT.1) THEN
       invdtmax = MAXVAL(this%kinvis(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) &
            / Mesh%dlx(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)**2) 
    ELSE
       ! set to zero, i.e. no limit in x-direction due to diffusion
       invdtmax = 0.0
    END IF
    ! y-direction
    IF (Mesh%JNUM.GT.1) THEN
       invdt = MAXVAL(this%kinvis(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) &
            / Mesh%dly(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)**2) 
    ELSE
       ! set to zero, i.e. no limit in y-direction due to diffusion
       invdt = 0.0
    END IF
    ! largest time step due to diffusion
    invdtmax = MAX(invdt, invdtmax)

!part two
    DO i=Mesh%IMIN,Mesh%IMAX
      DO j=Mesh%JMIN,Mesh%JMAX
        this%tmp(i,j) = SGSdiffusion(this%delta(i,j),&
                        cvar(i,j,Physics%SGSENERGY),cvar(i,j,Physics%DENSITY)) &
                        /cvar(i,j,Physics%DENSITY)
      END DO
    END DO

    ! x-direction
    IF (Mesh%INUM.GT.1) THEN
       invdt = MAXVAL(this%tmp(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) &
            / Mesh%dlx(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)**2)
    ELSE
       ! set to zero, i.e. no limit in x-direction due to diffusion
       invdt = 0.0
    END IF

    invdtmax = MAX(invdt, invdtmax)
    ! y-direction
    IF (Mesh%JNUM.GT.1) THEN
       invdt = MAXVAL(this%tmp(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) &
            / Mesh%dly(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)**2)
    ELSE
       ! set to zero, i.e. no limit in y-direction due to diffusion
       invdt = 0.0
    END IF
    ! largest time step due to diffusion
    invdtmax = MAX(invdt, invdtmax)

!What's wrong in this statement?
!     IF (invdt_x .GT. TINY(invdt_x) )THEN
!   print *, dt, Cdiff / invdt_x
    IF (invdtmax .GT. 1.0E-30 )THEN
       dt = Cdiff / invdtmax
!  print *, "Zeitskala sgs      :", dt
    ELSE
       dt = HUGE(dt)
    END IF
  END SUBROUTINE CalcTimestep_sgs


  PURE SUBROUTINE SGStransport(this,Mesh,delta,rho,Ksgs,diff) 
     IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: delta,rho,Ksgs,diff
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    REAL              :: k
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,delta,rho,Ksgs 
    INTENT(INOUT)     :: this                        ! provide temp storage  !
    INTENT(OUT)       :: diff
    !------------------------------------------------------------------------!
    DO i=Mesh%IMIN,Mesh%IMAX
      DO j=Mesh%JMIN,Mesh%JMAX
        k = SGSdiffusion(delta(i,j),Ksgs(i,j),rho(i,j))
        this%tmp2(i,j) = k*(Ksgs(i+1,j)/rho(i+1,j)-Ksgs(i-1,j)/rho(i-1,j))*0.5/Mesh%dlx(i,j)
        this%tmp3(i,j) = k*(Ksgs(i,j+1)/rho(i,j+1)-Ksgs(i,j-1)/rho(i,j-1))*0.5/Mesh%dly(i,j)
      END DO
    END DO 
    
    CALL Divergence(Mesh,this%tmp2(:,:),this%tmp3(:,:),diff(:,:))

  END SUBROUTINE SGStransport

  ELEMENTAL FUNCTION SGSdissipation(delta,Ksgs,rho) RESULT(SGSdiss)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
     REAL, INTENT(IN) :: delta,Ksgs,rho
     REAL             :: SGSdiss
    !------------------------------------------------------------------------!
     SGSdiss = Ceps/delta*sqrt(Ksgs**3/rho)
  END FUNCTION SGSdissipation


  ELEMENTAL FUNCTION SGSdiffusion(delta,Ksgs,rho) RESULT(SGSdiff)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
     REAL, INTENT(IN) :: delta,Ksgs,rho
     REAL             :: SGSdiff
    !------------------------------------------------------------------------!
     SGSdiff = Cd*delta*sqrt(rho*Ksgs)
  END FUNCTION SGSdiffusion 


  ELEMENTAL FUNCTION SGSproduction(txx,txy,txz,tyy,tyz,tzz,Sxx,Sxy,Sxz,Syy,Syz,Szz) RESULT(Sigma)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
     REAL, INTENT(IN) :: txx,txy,txz,tyy,tyz,tzz,Sxx,Sxy,Sxz,Syy,Syz,Szz
     REAL             :: Sigma 
    !------------------------------------------------------------------------!
    ! works also in 2D (Sxz,Syz,Szz = 0)
    Sigma = txx*Sxx+tyy*Syy+tzz*Szz+2.0*(txy*Sxy+txz*Sxz+tyz*Syz) 
  END FUNCTION SGSproduction

  ELEMENTAL FUNCTION C1SCHMIDT(delta,rho,ksgs) RESULT(nu)
    !------------------------------------------------------------------------!
     REAL, INTENT(IN) :: delta,rho,ksgs
     REAL             :: nu
    !------------------------------------------------------------------------!
      nu = C1*delta*sqrt(2.0*rho*ksgs) 
  END FUNCTION C1SCHMIDT
 
  SUBROUTINE CloseSources_sgs(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
     DEALLOCATE(this%btxx, &
         this%btyy, &
         this%btzz, &
         this%btxy, &
         this%btxz, &
         this%btyz, &
         this%Sxx, &
         this%Syy, &
         this%Szz, &
         this%Sxy, &
         this%Sxz, &
         this%Syz, &
         this%delta, & 
         this%sigma, &
         this%diff, &
         this%rhoeps, &
         this%tmp, &
         this%tmp2, &
         this%tmp3)
    CALL CloseSources(this)
  END SUBROUTINE CloseSources_sgs
 

END MODULE sources_sgs
