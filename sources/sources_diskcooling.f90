!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_diskcooling.f90                                           #
!#                                                                           #
!# Copyright (C) 2011                                                        #
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
!> \author Anna Feiler
!!
!! \brief source terms module for gray cooling of geometrically thin
!! accretion disks
!!
!! uses opacities taken from Bell & Lin, ApJ, 427, 1994 to
!! compute Rosseland mean using the interpolation formula of Gail 2003
!! (private communication)
!!
!! \warning use SI units
!!
!! \extends sources_c_accel
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_diskcooling
  USE constants_common, ONLY : KE
  USE sources_c_accel
  USE gravity_generic
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: source_name = "gray accretion disk cooling"
!!$  LOGICAL, PARAMETER :: DEBUG = .TRUE.
!!  LOGICAL, PARAMETER :: DEBUG = .FALSE.
  REAL, PARAMETER :: SQRT_THREE = 1.73205080757
  REAL, PARAMETER :: SQRT_TWOPI = 2.50662827463
  REAL, PARAMETER :: T0 = 3000      ! temperature constant (opacity interpolation)

  ! Rosseland mean opacity constants in SI units;
  ! taken from Bell & Lin, ApJ, 427, 1994
  ! kappa_i= kappa_0i * rho**rexp(i) * T**Texp(i)

!!$  DOUBLE PRECISION, PARAMETER :: kappa0(8) = (/ &
!!$       2.00D-05, & ! ice grains                     [m^2/kg/K^2]
!!$       2.00D+15, & ! evaporation of ice grains      [m^2/kg*K^7]
!!$       1.00D-02, & ! metal grains                 [m^2/kg/K^0.5]
!!$       2.00D+77, & ! evaporation of metal grains [m^5/kg^2*K^24]
!!$       1.00D-11, & ! molecules                [m^4/kg^(5/3)/K^3]
!!$       1.00D-38, & ! H-scattering            [m^3/kg^(4/3)/K^10]
!!$       1.50D+16, & ! bound-free and free-free [m^5/kg^2*K^(5/2)]
!!$       KE /)       ! electron scattering                [m^2/kg]

  REAL, PARAMETER :: logkappa0(8) = (/ &
       -10.8197782844, & ! ice grains                     [m^2/kg/K^2]
       35.2319235755, &  ! evaporation of ice grains      [m^2/kg*K^7]
       -4.60517018599, & ! metal grains                 [m^2/kg/K^0.5]
       177.992199341, &  ! evaporation of metal grains [m^5/kg^2*K^24]
       -25.3284360229, & ! molecules                [m^4/kg^(5/3)/K^3]
       -87.4982335338, & ! H-scattering            [m^3/kg^(4/3)/K^10]
       37.246826596, &   ! bound-free and free-free [m^5/kg^2*K^(5/2)]
       -3.3581378922 /)  ! electron scattering                [m^2/kg]
  REAL, PARAMETER :: Texp(8) = (/ 2.0, -7.0, 0.5, -24.0, 3.0, 10.0, -2.5, 0.0 /)
  REAL, PARAMETER :: rexp(8) = (/ 0.0, 0.0, 0.0, 1.0, 2./3., 1./3., 1.0, 0.0 /)
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! methods
       InitSources_diskcooling, &
       ExternalSources_diskcooling, &
       CalcTimestep_diskcooling, &
       CloseSources_diskcooling
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources_diskcooling(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Dict_TYP),POINTER :: config,IO
    INTEGER           :: stype
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "stype", stype)
    CALL InitSources(this,stype,source_name)
    ! some sanity checks
    IF (GetType(Physics).NE.EULER2D) &
         CALL Error(this,"InitSources_diskcooling","physics not supported")
    IF (GetType(Physics%constants).NE.SI) &
         CALL Error(this,"InitSources_diskcooling","only SI units supported")

    ALLOCATE(this%tmp(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX),&
         this%Qcool(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         STAT=err)
    IF (err.NE.0) CALL Error(this,"InitSources_diskcooling","memory allocation failed")

    ! Courant number, i.e. safety factor for numerical stability
    CALL RequireKey(config, "cvis", 0.9)
    CALL GetAttr(config, "cvis", this%cvis)

    ! set initial time < 0
    this%time = -1.0

    ! initialize arrays
    this%tmp(:,:)   = 0.0
    this%Qcool(:,:) = 0.0

    ! get speed of sound pointer
    this%cs => GetSoundSpeeds(Physics)

    !initialise output
    CALL SetOutput(this,Mesh,Physics,config,IO)

  END SUBROUTINE InitSources_diskcooling

  
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

    !cooling energy source term
    valwrite = 0 
    IF (HasKey(config, "output/Qcool")) CALL GetAttr(config, "output/Qcool", valwrite)
    IF (valwrite .EQ. 1) THEN
         CALL AddField(IO, &
               "Qcool", &
               this%Qcool,&
               Dict("name" / "Qcool"))
    END IF 

  END SUBROUTINE SetOutput


  SUBROUTINE ExternalSources_diskcooling(this,Mesh,Physics,Fluxes,time,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,time,pvar,cvar
    INTENT(INOUT)     :: Physics
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    sterm(:,:,Physics%DENSITY) = 0.0
    sterm(:,:,Physics%XMOMENTUM) = 0.0
    sterm(:,:,Physics%YMOMENTUM) = 0.0

    CALL UpdateCooling(this,Mesh,Physics,Fluxes,time,pvar)
    ! energy loss due to radiation processes
    sterm(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%ENERGY) = &
         -this%Qcool(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)

  END SUBROUTINE ExternalSources_diskcooling

 
  SUBROUTINE CalcTimestep_diskcooling(this,Mesh,Physics,Fluxes,time,pvar,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                      :: pvar
    REAL              :: time,dt
    !------------------------------------------------------------------------!
    REAL              :: invdt
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,time
    INTENT(INOUT)     :: Physics
    INTENT(OUT)       :: dt
    !------------------------------------------------------------------------!
    CALL UpdateCooling(this,Mesh,Physics,Fluxes,time,pvar)
    ! maximum of inverse cooling timescale t_cool ~ P/Q_cool
    invdt = MAXVAL(ABS(this%Qcool(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) &
         / pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE)))
    IF (invdt.GT.TINY(invdt)) THEN
       dt = this%cvis / invdt
    ELSE
       dt = HUGE(invdt)
    END IF
  END SUBROUTINE CalcTimestep_diskcooling


  SUBROUTINE UpdateCooling(this,Mesh,Physics,Fluxes,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP),POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP) :: Fluxes
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    !------------------------------------------------------------------------!
    INTEGER, SAVE     :: counter = 0
    INTEGER           :: i,j
    REAL              :: cs,kappa,tau,tau_eff,muRgamma,sqrtgamma,logrho,logT
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,pvar
    INTENT(INOUT)     :: Physics
    !------------------------------------------------------------------------!
    ! enery loss due to radiation processes
    IF ((time.NE.this%time) .OR. (time .EQ. 0.)) THEN
        ! make sure that the sound spped is up to date
       CALL UpdateSoundSpeed(Physics,Mesh,time,pvar)
       ! compute mu/R_G/gamma
       muRgamma = Physics%mu / (Physics%Constants%RG * Physics%gamma)
       ! compute sqrt of the ratio of specific heats
       sqrtgamma = SQRT(Physics%gamma)

       ! disk height
!FIXME!!!!! not 1 => Gravity!!!!
       this%tmp(:,:) = GetDiskHeight(GetSourcesPointer(&
                         Physics%Sources,1),Mesh,Physics,Fluxes,time,pvar)
       
       ! compute gray cooling term
!CDIR OUTERUNROLL=8
       DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
          DO i=Mesh%IMIN,Mesh%IMAX

             ! logarithm of midplane density
             ! SQRT(2*Pi)^-1 * Sigma / H    
             logrho = LOG(pvar(i,j,Physics%DENSITY) / &
                  (TINY(SQRT_TWOPI) + SQRT_TWOPI*this%tmp(i,j)))

             ! logarithm of midplane temperature
             logT = LOG(muRgamma*this%cs(i,j)*this%cs(i,j))

             ! compute Rosseland mean absorption coefficient
!CDIR IEXPAND
             kappa = RosselandMeanOpacity_new(logrho,logT)
             ! optical depth
             tau = 0.5*kappa*pvar(i,j,Physics%DENSITY)
             ! effective optical depth
             tau_eff = 0.5*tau + 1./(3.*tau) + 1./SQRT_THREE
             ! cooling term
             this%Qcool(i,j) = (8.*Physics%Constants%SB)/(3.*tau_eff) * &
                  EXP(4.0*logT) ! = T**4
          END DO
       END DO

       this%time=time
    END IF
  END SUBROUTINE UpdateCooling


  SUBROUTINE CloseSources_diskcooling(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%tmp,this%Qcool)
    CALL CloseSources(this)
  END SUBROUTINE CloseSources_diskcooling
 


  ELEMENTAL FUNCTION RosselandMeanOpacity_new(logrho,logT) RESULT(kappa_R)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: logrho,logT
    REAL             :: kappa_R
    !------------------------------------------------------------------------!
    REAL :: kappa4(8)
    REAL :: Tfactor
    !------------------------------------------------------------------------!
    ! compute kappa_i^4 for all cooling mechanisms
!CDIR UNROLL=8
    kappa4(:) = EXP(MAX(-40.0,MIN(40.0, &
         4.*(logkappa0(:)+rexp(:)*logrho+Texp(:)*logT))))
    ! compute (T/T0)**10
    Tfactor = EXP(MAX(-40.0,MIN(40.0,10.*(logT-LOG(T0)))))
    ! compute Rosseland mean using Gails interpolation formula
    kappa_R = 1. /(TINY(kappa_R) + &
           (1./kappa4(1) + 1./(1.+Tfactor)/(kappa4(2)+kappa4(3)))**0.25 &
         + (1./(kappa4(4)+kappa4(5)+kappa4(6)) + 1./(kappa4(7)+kappa4(8)))**0.25)
  END FUNCTION RosselandMeanOpacity_new

END MODULE sources_diskcooling
