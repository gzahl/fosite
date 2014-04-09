!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_viscosity.f90                                             #
!#                                                                           #
!# Copyright (C) 2008-2011                                                   #
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
!----------------------------------------------------------------------------!
! viscosity of Newtonian fluid
!----------------------------------------------------------------------------!
MODULE sources_viscosity
  USE common_types, ONLY : Common_TYP, InitCommon
  USE sources_pointmass
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: source_name = "viscosity of Newtonian fluid"
  ! flags for viscosity model
  INTEGER, PARAMETER :: MOLECULAR = 1     ! constant viscosity
  INTEGER, PARAMETER :: ALPHA     = 2     ! Shakura-Sunyaev prescription
  INTEGER, PARAMETER :: BETA      = 3     ! Duschl prescription
  INTEGER, PARAMETER :: PRINGLE   = 4     ! constant kinematic viscosity
  CHARACTER(LEN=32), PARAMETER, DIMENSION(4) :: viscosity_name = (/ &
                                     "constant viscosity              ", &
                                     "turbulent Shakura-Sunyaev       ", &
                                     "turbulent Duschl                ", &
                                     "const. kinematic viscosity      " /)
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! constants
       MOLECULAR, ALPHA, BETA, PRINGLE, &
       ! methods
       InitSources_viscosity, &
       InfoSources_viscosity, &
       ExternalSources_viscosity, &
       CalcTimestep_viscosity, &
       CloseSources_viscosity
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources_viscosity(this,Mesh,Physics,Fluxes,stype,model,dynconst, &
       bulkconst,cvis)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    INTEGER           :: stype
    INTEGER,OPTIONAL  :: model
    REAL,OPTIONAL     :: dynconst,bulkconst,cvis
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Fluxes,stype,model,dynconst,bulkconst,cvis
    !------------------------------------------------------------------------!
    IF (.NOT.Initialized(Fluxes)) &
         CALL Error(this,"InitSources_viscosity","fluxes module uninitialized")

    CALL InitSources(this,stype,source_name)

    ! check mesh
    IF (GetType(Fluxes).NE.MIDPOINT) &
         CALL Error(this,"InitSources_viscosity", &
         "only midpoint rule is currently supported")

    ! viscosity model
    IF (PRESENT(model)) THEN
       CALL InitCommon(this%viscosity,model,viscosity_name(model))
    ELSE
       ! default: molecular viscosity
       CALL InitCommon(this%viscosity,MOLECULAR,viscosity_name(MOLECULAR))
    END IF

    ! dynamic viscosity constant
    IF (PRESENT(dynconst)) THEN
       this%dynconst = dynconst
    ELSE
       this%dynconst  = 0.1
    END IF
    
    ! bulk viscosity constant (disabled by default)
    IF (PRESENT(bulkconst)) THEN
       this%bulkconst =  bulkconst
    ELSE
       this%bulkconst =  0.0
    END IF

    ! set viscous courant number
    IF (PRESENT(cvis)) THEN
       this%cvis = cvis
    ELSE
       this%cvis = 0.5
    END IF

    ALLOCATE(this%dynvis(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%kinvis(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%bulkvis(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btxx(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btyy(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btzz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btxy(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btxz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btyz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         STAT=err)
    IF (err.NE.0) &
         CALL Error(this,"InitSources_viscosity","Memory allocation failed.")

    ! do the model specific initialization
    SELECT CASE(GetType(this%viscosity))
    CASE(MOLECULAR,PRINGLE)
       ! do nothing
    CASE(ALPHA)
       ! check physics
       SELECT CASE(GetType(Physics))
       CASE(EULER2D,EULER2D_ISOTHERM,EULER3D_ROTSYM)
          ! do nothing
       CASE DEFAULT
          CALL Error(this,"InitSources_viscosity",&
               "Physics not supported for alpha-viscosity")
       END SELECT
       ALLOCATE(this%invr(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
            STAT=err)
       IF (err.NE.0) CALL Error(this,"InitSources_viscosity",&
               "Unable to allocate memory.")
       ! check geometry
       SELECT CASE(GetType(Mesh%geometry))
       CASE(POLAR,LOGPOLAR,TANPOLAR,SINHPOLAR)
          ! compute inverse of distance to origin
          this%invr(:,:) = 1./(TINY(1.0)+SQRT(Mesh%bccart(:,:,1)**2 &
               + Mesh%bccart(:,:,2)**2))
       CASE(CYLINDRICAL,TANCYLINDRICAL,SPHERICAL,OBLATE_SPHEROIDAL)
          ! compute inverse of distance to axis
          this%invr(:,:) = 1./(TINY(1.0)+ABS(Mesh%bccart(:,:,1)))
       CASE DEFAULT
          CALL Error(this,"InitSources_viscosity",&
               "Geometry not supported for alpha-viscosity")
       END SELECT
    CASE(BETA)
       SELECT CASE(GetType(Physics))
       CASE (EULER2D,EULER2D_ISOTHERM,EULER3D_ROTSYM)
          ! do nothing
       CASE DEFAULT
          CALL Error(this,"InitSources_viscosity",&
               "Physics not supported for beta-viscosity")
       END SELECT
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
  END SUBROUTINE InitSources_viscosity


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
       CALL Info(this,"            dynamic viscosity: " // TRIM(dynconst_str) // &
           ACHAR(10)//"            bulk viscosity:    " // TRIM(bulkconst_str))
    CASE(ALPHA)
       CALL Info(this,"            alpha:             " // TRIM(dynconst_str))
    CASE(BETA)
       CALL Info(this,"            beta:              " // TRIM(dynconst_str))
    CASE(PRINGLE)
       CALL Info(this,"            kinemat. viscosity:" // TRIM(dynconst_str))
    END SELECT
  END SUBROUTINE InfoSources_viscosity


  PURE SUBROUTINE UpdateViscosity(this,Mesh,Physics,time,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,kp,kv
    REAL              :: P,Omega,cs2
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,time,pvar,cvar
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! only update viscosity if time has changed
    IF (time.NE.this%time) THEN
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
          ! set some constants depending on physics
!CDIR IEXPAND
          SELECT CASE(GetType(Physics))
          CASE(EULER2D)
             kp  = Physics%PRESSURE
             kv  = Physics%YVELOCITY
             cs2 = 1.0
          CASE(EULER2D_ISOTHERM)
             kp  = Physics%DENSITY
             kv  = Physics%YVELOCITY
             cs2 = Physics%csiso**2
          CASE(EULER3D_ROTSYM)
             kp  = Physics%PRESSURE
             kv  = Physics%ZVELOCITY
             cs2 = 1.0
          END SELECT
!CDIR COLLAPSE
          DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
             DO i=Mesh%IGMIN,Mesh%IGMAX
                ! get pressure and angular velocity
                P = cs2*pvar(i,j,kp)
                Omega = pvar(i,j,kv)*this%invr(i,j)
                ! compute alpha viscosity
                this%dynvis(i,j) = etafkt_alpha(this%dynconst,P,Omega)
             END DO
          END DO
       CASE(BETA)
          ! Duschl type beta viscosity
!CDIR IEXPAND
          SELECT CASE(GetType(Physics))
          CASE (EULER2D,EULER2D_ISOTHERM)
             this%dynvis(:,:) = etafkt_beta(this%dynconst, &
                  Mesh%bhy(:,:)*cvar(:,:,Physics%YMOMENTUM))
          CASE (EULER3D_ROTSYM)
             this%dynvis(:,:) = etafkt_beta(this%dynconst, &
                  Mesh%bhz(:,:)*cvar(:,:,Physics%ZMOMENTUM))
          CASE (EULER3D_ROTAMT)
             this%dynvis(:,:) = etafkt_beta(this%dynconst, &
                  cvar(:,:,Physics%ZMOMENTUM))
          END SELECT
       CASE(PRINGLE) ! constant kinematic viscosity
!CDIR COLLAPSE
          DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
             DO i=Mesh%IGMIN,Mesh%IGMAX
                this%dynvis(i,j) = etafkt_pringle(this%dynconst, &
                     pvar(i,j,Physics%DENSITY))
             END DO
          END DO
       END SELECT
       ! set viscosity in ghost cells to the viscosity
       ! in adjacent cells inside the computational domain
       IF (GetType(this%viscosity).NE.MOLECULAR) THEN
          DO i=1,Mesh%GNUM
             this%dynvis(Mesh%IMIN-i,:) = this%dynvis(Mesh%IMIN,:)
             this%dynvis(Mesh%IMAX+i,:) = this%dynvis(Mesh%IMAX,:)
             this%dynvis(:,Mesh%JMIN-i) = this%dynvis(:,Mesh%JMIN)
             this%dynvis(:,Mesh%JMAX+i) = this%dynvis(:,Mesh%JMAX)
          END DO
       END IF
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
      eta = alpha * P / omega
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

 
  PURE SUBROUTINE ExternalSources_viscosity(this,Mesh,Physics,time,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP),POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: cvar,pvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,time,cvar,pvar
    INTENT(INOUT)     :: Physics
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    CALL UpdateViscosity(this,Mesh,Physics,time,pvar,cvar)
    CALL CalculateStresses(Physics,Mesh,pvar,this%dynvis,this%bulkvis, &
             this%btxx,this%btxy,this%btxz,this%btyy,this%btyz,this%btzz)
    CALL ViscositySources(Physics,Mesh,pvar,this%btxx,this%btxy,this%btxz, &
             this%btyy,this%btyz,this%btzz,sterm)
  END SUBROUTINE ExternalSources_viscosity

  
  PURE SUBROUTINE CalcTimestep_viscosity(this,Mesh,Physics,time,pvar,cvar,dt)
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
    REAL              :: invdt_x, invdt_y
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,time,pvar,cvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: dt
    !------------------------------------------------------------------------!
    ! update dynamic and bulk viscosity
    CALL UpdateViscosity(this,Mesh,Physics,time,pvar,cvar)

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
    DEALLOCATE(this%dynvis,this%kinvis,this%bulkvis &
               ,this%btxx,this%btyy,this%btzz,this%btxy,this%btxz,this%btyz)
    CALL CloseSources(this)
  END SUBROUTINE CloseSources_viscosity
 

END MODULE sources_viscosity
