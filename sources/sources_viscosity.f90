!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_viscosity.f90                                             #
!#                                                                           #
!# Copyright (C) 2008-2012                                                   #
!# Bjoern Sperling <sperling@astrophysik.uni-kiel.de>                        #
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
!# along with Mesh program; if not, write to the Free Software               #
!# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
!#                                                                           #
!#############################################################################
!> \addtogroup sources
!! - parameters of viscosity source as key-values
!! \key{vismodel,INTEGER,viscosity model}
!! \key{dynconst,REAL,dynamic viscosity constant,0.1}
!! \key{bulkconst,REAL,bulk viscosity constant,-2/3*dynconst}
!! \key{cvis,REAL,viscous courant number,0.5}
!! \key{output/stress,INTEGER,enable(=1) output of the stress tensor,0}
!! \key{output/dynvis,INTEGER,enable(=1) output of dynamic viscosity,0}
!! \key{output/kinvis,INTEGER,enable(=1) output of kinematic viscosity,0}
!! \key{output/bulkvis,INTEGER,enable(=1) output of bulk viscosity,0}
!----------------------------------------------------------------------------!
!> \author Björn Sperling
!! \author Tobias Illenseer
!!
!! \brief viscosity of Newtonian fluid
!!
!! \extends sources_c_accel
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_viscosity
  USE common_types, ONLY : Common_TYP, InitCommon
  USE sources_c_accel
  USE gravity_generic
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: source_name = "viscosity of Newtonian fluid"
  ! flags for viscosity model
  INTEGER, PARAMETER :: MOLECULAR = 1     ! constant viscosity
  INTEGER, PARAMETER :: ALPHA     = 2     ! Shakura-Sunyaev prescription
  INTEGER, PARAMETER :: BETA      = 3     ! Duschl prescription
  INTEGER, PARAMETER :: PRINGLE   = 4     ! constant kinematic viscosity
  INTEGER, PARAMETER :: ALPHA_ALT = 5     ! alternative Shakura-Sunyaev
  CHARACTER(LEN=32), PARAMETER, DIMENSION(5) :: viscosity_name = (/ &
                                     "constant viscosity              ", &
                                     "turbulent Shakura-Sunyaev       ", &
                                     "turbulent Duschl                ", &
                                     "const. kinematic viscosity      ", &
                                     "alternative Shakura-Sunyaev     "/)
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! constants
       MOLECULAR, ALPHA, BETA, PRINGLE, ALPHA_ALT, &
       ! methods
       InitSources_viscosity, &
       InfoSources_viscosity, &
       ExternalSources_viscosity, &
       CalcTimestep_viscosity, &
       CloseSources_viscosity
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources_viscosity(this,Mesh,Physics,Fluxes,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Dict_TYP),POINTER :: config,IO
    INTEGER           :: stype
    INTEGER           :: model
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Fluxes
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "stype", stype)
    IF (.NOT.Initialized(Fluxes)) &
         CALL Error(this,"InitSources_viscosity","fluxes module uninitialized")

    CALL InitSources(this,stype,source_name)

    ! check mesh
    IF (GetType(Fluxes).NE.MIDPOINT) &
         CALL Error(this,"InitSources_viscosity", &
         "only midpoint rule is currently supported")

    ! viscosity model
    CALL RequireKey(config, "vismodel") ! no default!
    CALL GetAttr(config, "vismodel", model)
    CALL InitCommon(this%viscosity,model,viscosity_name(model))

    ! dynamic viscosity constant
    CALL RequireKey(config, "dynconst", 0.1)
    CALL GetAttr(config, "dynconst", this%dynconst)
    
    ! bulk viscosity constant (disabled by default)
    CALL RequireKey(config, "bulkconst", -2./3.*this%dynconst)
    CALL GetAttr(config, "bulkconst", this%bulkconst)

    ! set viscous courant number
    CALL RequireKey(config, "cvis", 0.5)
    CALL GetAttr(config, "cvis", this%cvis)

    ALLOCATE(this%dynvis(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%kinvis(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%bulkvis(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btxx(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btyy(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btzz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btxy(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btxz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btyz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%height(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         STAT=err)
    IF (err.NE.0) &
         CALL Error(this,"InitSources_viscosity","Memory allocation failed.")

    ! do the model specific initialization
    SELECT CASE(GetType(this%viscosity))
    CASE(MOLECULAR,PRINGLE)
       ! do nothing
    CASE(ALPHA)
       ALLOCATE(this%invr(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
            STAT=err)
       IF (err.NE.0) CALL Error(this,"InitSources_viscosity",&
               "Unable to allocate memory.")
       ! check geometry
       SELECT CASE(GetType(Mesh%geometry))
       CASE(POLAR,LOGPOLAR,TANPOLAR,SINHPOLAR)
          ! compute inverse of distance to origin
          this%invr(:,:) = 1./(TINY(1.0)+Physics%bcradius(:,:))
       CASE(CYLINDRICAL,TANCYLINDRICAL,SPHERICAL,OBLATE_SPHEROIDAL)
          ! compute inverse of distance to axis
          this%invr(:,:) = 1./(TINY(1.0)+Mesh%bhz(:,:))
       CASE DEFAULT
          CALL Error(this,"InitSources_viscosity",&
               "Geometry not supported for alpha-viscosity")
       END SELECT
    CASE(BETA)
       SELECT CASE(GetType(Physics))
       CASE (EULER2D,EULER2D_ISOTHERM,&
             EULER2D_IAMT,EULER2D_ISOIAMT,&
             EULER3D_ROTSYM,EULER3D_ROTAMT,EULER3D_ROTSYMSGS)
          ! do nothing
       CASE DEFAULT
          CALL Error(this,"InitSources_viscosity",&
               "Physics not supported for beta-viscosity")
       END SELECT
    CASE(ALPHA_ALT)
       IF (Physics%DIM.NE.2) &
          CALL Error(this,"InitSources_viscosity",&
               "alternative alpha-viscosity works only for flat disks")
    CASE DEFAULT
       CALL Error(this,"InitSources_viscosity",&
            "Viscosity prescription not supported.")
    END SELECT

    ! set initial time to negative value;
    ! this guarantees update of viscosity for time t=0.0
    this%time = -1.0
    ! initialize viscosity arrays
    this%dynvis(:,:)  = this%dynconst
    this%bulkvis(:,:) = this%bulkconst
    ! set stress tensor arrays to zero
    this%btxx(:,:) = 0.0
    this%btyy(:,:) = 0.0
    this%btzz(:,:) = 0.0
    this%btxy(:,:) = 0.0
    this%btxz(:,:) = 0.0
    this%btyz(:,:) = 0.0

    CALL SetOutput(this,Mesh,Physics,config,IO)
  END SUBROUTINE InitSources_viscosity

  SUBROUTINE SetOutput(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP),POINTER    :: this
    TYPE(Mesh_TYP)       :: Mesh
    TYPE(Physics_TYP)    :: Physics
    TYPE(Dict_TYP),POINTER  :: config,IO
    !------------------------------------------------------------------------!
    INTEGER              :: valwrite
    !------------------------------------------------------------------------!
    INTENT(IN)           :: Mesh,Physics
    !------------------------------------------------------------------------! 
    valwrite = 0
    IF (HasKey(config, "output/stress")) CALL GetAttr(config, "output/stress", valwrite)
    IF (valwrite .EQ. 1) THEN
       CALL AddField(IO, &
               "stress:xx", &
               this%btxx, &
               Dict("name" / "stress:xx"))

       CALL AddField(IO, &
               "stress:xy", &
               this%btxy, &
               Dict("name" / "stress:xy"))

       CALL AddField(IO, &
               "stress:yy", &
               this%btyy, &
               Dict("name" / "stress:yy"))

     SELECT CASE(GetType(Physics))
       CASE (EULER3D_ROTSYM,EULER3D_ROTAMT,EULER3D_ROTSYMSGS)
         CALL AddField(IO, &
               "stress:xz", &
               this%btxz, &
               Dict("name" / "stress:xz"))

         CALL AddField(IO, &
               "stress:yz", &
               this%btyz, &
               Dict("name" / "stress:yz"))

         CALL AddField(IO, &
               "stress:zz", &
               this%btzz, &
               Dict("name" / "stress:zz"))

      END SELECT
    END IF

    valwrite = 0 
    IF (HasKey(config, "output/dynvis")) CALL GetAttr(config, "output/dynvis", valwrite)
    IF (valwrite .EQ. 1) THEN
         CALL AddField(IO, &
               "dynvis", &
               this%dynvis, &
               Dict("name" / "dynvis"))
    END IF 


    valwrite = 0
    IF (HasKey(config, "output/kinvis")) CALL GetAttr(config, "output/kinvis", valwrite)
    IF (valwrite .EQ. 1) THEN
         CALL AddField(IO, &
               "kinvis", &
               this%kinvis, &
               Dict("name" / "kinvis"))
    END IF 

    valwrite = 0 
    IF (HasKey(config, "output/bulkvis")) CALL GetAttr(config, "output/bulkvis", valwrite)
    IF (valwrite .EQ. 1) THEN
         CALL AddField(IO, &
               "bulkvis", &
               this%bulkvis, &
               Dict("name" / "bulkvis"))
    END IF 

  END SUBROUTINE SetOutput


  SUBROUTINE InfoSources_viscosity(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: dynconst_str,bulkconst_str
    !------------------------------------------------------------------------!
    WRITE (dynconst_str,'(ES9.2)') this%dynconst
    WRITE (bulkconst_str,'(ES9.2)') this%bulkconst
    CALL Info(this,"            viscosity model:   " // GetName(this%viscosity))
    SELECT CASE(GetType(this%viscosity))
    CASE(MOLECULAR)
       CALL Info(this,"            dynamic viscosity: " // TRIM(dynconst_str))
       CALL Info(this,"            bulk viscosity:    " // TRIM(bulkconst_str))
    CASE(ALPHA,ALPHA_ALT)
       CALL Info(this,"            alpha:             " // TRIM(dynconst_str))
    CASE(BETA)
       CALL Info(this,"            beta:              " // TRIM(dynconst_str))
    CASE(PRINGLE)
       CALL Info(this,"            kinemat. viscosity:" // TRIM(dynconst_str))
    END SELECT
  END SUBROUTINE InfoSources_viscosity


  SUBROUTINE UpdateViscosity(this,Mesh,Physics,Fluxes,time,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,kp=0,kv=0
    REAL              :: P,Omega,rotOmega,cs2=1.0,height
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,time,pvar,cvar
    INTENT(INOUT)     :: Physics
    !------------------------------------------------------------------------!
    ! only update viscosity if time has changed
    IF ((time.NE.this%time) .OR. (time .EQ. 0.)) THEN
       SELECT CASE(GetType(this%viscosity))
       CASE(MOLECULAR)
          ! do nothing, since it has already been initialized
          ! in InitSources_viscosity
       CASE(ALPHA)
          ! Shakura-Sunyaev type alpha viscosity
          ! standard alpha prescription: nu = alpha*cs*h
          ! or nu = alpha / A * cs**2 / omega with A = -d ln(omega)/d ln(r)
          ! (see Kato, Fukue & Minishige: Black Hole Accretion Disks, 2008; 
          ! equation (3.46))
          ! 
          ! this is a rough estimation assuming that the logarithmic derivative
          ! of the angular velocity (d ln(omega)/d ln(r)) is of the order of -1
          !
          ! check if PRESSURE exists
          IF (Physics%PRESSURE.EQ.0) THEN
             ! isothermal Physics
             kp  = Physics%DENSITY
             cs2 = Physics%csiso**2
          ELSE
             ! non-isothermal Physics
             kp  = Physics%PRESSURE
             cs2 = 1.0
          END IF
          ! check for 2D / 3D
          SELECT CASE(Physics%DIM)
          CASE(2)
             kv  = Physics%YVELOCITY
          CASE(3)
             kv  = Physics%ZVELOCITY
          END SELECT
!CDIR COLLAPSE
          DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
             DO i=Mesh%IGMIN,Mesh%IGMAX
                ! get/compute pressure and angular velocity
                P = cs2*pvar(i,j,kp)
                ! consider Omega of rotating frame => Physics%Omega 
                Omega = pvar(i,j,kv)*this%invr(i,j) + Physics%Omega
                ! compute alpha viscosity
                this%dynvis(i,j) = etafkt_alpha(this%dynconst,P,Omega)
                ! trace-free condition!
                this%bulkvis(i,j) = -2./3.*this%dynvis(i,j)
             END DO
          END DO
       CASE(BETA)
          ! Duschl type beta viscosity
!CDIR IEXPAND
          SELECT CASE(GetType(Physics))
          CASE (EULER2D,EULER2D_ISOTHERM)
             this%dynvis(:,:) = etafkt_beta(this%dynconst, &
                  Mesh%bhy(:,:)*cvar(:,:,Physics%YMOMENTUM))
          CASE (EULER2D_IAMT,EULER2D_ISOIAMT)
             this%dynvis(:,:) = etafkt_beta(this%dynconst, &
                  cvar(:,:,Physics%YMOMENTUM))
          CASE (EULER3D_ROTSYM,EULER3D_ROTSYMSGS)
             this%dynvis(:,:) = etafkt_beta(this%dynconst, &
                  Mesh%bhz(:,:)*cvar(:,:,Physics%ZMOMENTUM))
          CASE (EULER3D_ROTAMT)
             this%dynvis(:,:) = etafkt_beta(this%dynconst, &
                  cvar(:,:,Physics%ZMOMENTUM))
          END SELECT
          ! trace-free condition!
          this%bulkvis(:,:) = -2./3.*this%dynvis(:,:)
       CASE(PRINGLE) ! constant kinematic viscosity
!CDIR COLLAPSE
          DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
             DO i=Mesh%IGMIN,Mesh%IGMAX
                this%dynvis(i,j) = etafkt_pringle(this%dynconst, &
                     pvar(i,j,Physics%DENSITY))
                ! trace-free condition!
                this%bulkvis(i,j) = -2./3.*this%dynvis(i,j)
             END DO
          END DO
       CASE(ALPHA_ALT)
          ! eta = nu*rho = alpha*cs*h * rho
          ! get disk height and compute alpha viscosity
          ! disk height
!CDIR IEXPAND
!FIXME: 1 => GRAVITY!!!!
          this%height(:,:) = GetDiskHeight(GetSourcesPointer(Physics%Sources,1)&
                                        ,Mesh,Physics,Fluxes,time,pvar)
          CALL UpdateSoundSpeed(Physics,Mesh,time,pvar)
!CDIR COLLAPSE
          DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
             DO i=Mesh%IGMIN,Mesh%IGMAX
                ! dynamic viscosity
                this%dynvis(i,j) = this%dynconst * Physics%bccsound(i,j) &
                      * this%height(i,j) * pvar(i,j,Physics%DENSITY)
             END DO
          END DO

          ! trace-free condition!
          this%bulkvis(:,:) = -2./3.*this%dynvis(:,:)
       END SELECT

       ! update time
       this%time = time
    END IF
  CONTAINS
    ! some elemental functions for computation of the viscosity

    ! Alpha viscosity
    ELEMENTAL FUNCTION etafkt_alpha(alpha,P,omega) RESULT(eta)
      IMPLICIT NONE
      !----------------------------------------------------------------------!
      REAL, INTENT(IN) :: alpha,P,omega
      REAL :: eta
      !----------------------------------------------------------------------!
      eta = alpha * P / ABS(omega)
    END FUNCTION etafkt_alpha

    ! Beta viscosity
    ELEMENTAL FUNCTION etafkt_beta(beta,L) RESULT(eta)
      IMPLICIT NONE
      !----------------------------------------------------------------------!
      REAL, INTENT(IN) :: beta,L
      REAL :: eta
      !----------------------------------------------------------------------!
      ! L = r * rho * v_phi is the angular momentum
      eta = beta * abs(L)
    END FUNCTION etafkt_beta
    
    ! Pringle disk (kinematic viscosity is constant)
    ELEMENTAL FUNCTION etafkt_pringle(nu,rho) RESULT(eta)
      IMPLICIT NONE
      !----------------------------------------------------------------------!
      REAL, INTENT(IN) :: nu,rho
      REAL :: eta
      !----------------------------------------------------------------------!
      eta = nu * rho   !nu is the kinematic viscosity
    END FUNCTION etafkt_pringle
    
  END SUBROUTINE UpdateViscosity

 
  SUBROUTINE ExternalSources_viscosity(this,Mesh,Physics,Fluxes,time,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP),POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: cvar,pvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,time,cvar,pvar
    INTENT(INOUT)     :: Physics
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    CALL UpdateViscosity(this,Mesh,Physics,Fluxes,time,pvar,cvar)
    CALL CalculateStresses(Physics,Mesh,pvar,this%dynvis,this%bulkvis, &
             this%btxx,this%btxy,this%btxz,this%btyy,this%btyz,this%btzz)
    CALL ViscositySources(Physics,Mesh,pvar,this%btxx,this%btxy,this%btxz, &
             this%btyy,this%btyz,this%btzz,sterm)
  END SUBROUTINE ExternalSources_viscosity

  
  SUBROUTINE CalcTimestep_viscosity(this,Mesh,Physics,Fluxes,time,pvar,cvar,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar,cvar
    REAL              :: dt
    !------------------------------------------------------------------------!
    REAL              :: invdt_x, invdt_y
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,time,pvar,cvar
    INTENT(INOUT)     :: Physics
    INTENT(OUT)       :: dt
    !------------------------------------------------------------------------!
    ! update dynamic and bulk viscosity
    CALL UpdateViscosity(this,Mesh,Physics,Fluxes,time,pvar,cvar)

    ! compute kinematic viscosity
    this%kinvis(:,:) = this%dynvis(:,:) / pvar(:,:,Physics%DENSITY)

    ! x-direction
    IF (Mesh%INUM.GT.1) THEN
       invdt_x = MAXVAL(this%kinvis(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) &
            / Mesh%dlx(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)**2)
    ELSE
       ! set to zero, i.e. no limit in x-direction due to diffusion
       invdt_x = 0.0
    END IF
    ! y-direction
    IF (Mesh%JNUM.GT.1) THEN
       invdt_y = MAXVAL(this%kinvis(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) &
            / Mesh%dly(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)**2)
    ELSE
       ! set to zero, i.e. no limit in y-direction due to diffusion
       invdt_y = 0.0
    END IF
    ! largest time step due to diffusion
    invdt_x = MAX(invdt_x, invdt_y)
    IF (invdt_x.GT.TINY(invdt_x)) THEN
       dt = this%cvis / invdt_x
    ELSE
       dt = HUGE(dt)
    END IF
  END SUBROUTINE CalcTimestep_viscosity


  SUBROUTINE CloseSources_viscosity(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%dynvis,this%kinvis,this%bulkvis,this%height &
               ,this%btxx,this%btyy,this%btzz,this%btxy,this%btxz,this%btyz)
    CALL CloseSources(this)
  END SUBROUTINE CloseSources_viscosity
 

END MODULE sources_viscosity
