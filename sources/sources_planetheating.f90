!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_planet_heating.f90                                        #
!#                                                                           #
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
!> \author Jannes Klee
!!
!! \brief heating of a planet
!!
!! \warning use SI units
!!
!! \extends sources_common
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_planetheating
  USE constants_common,ONLY :SB
  USE mesh_generic
  USE gravity_binary
  USE gravity_pointmass
  USE sources_common
  USE common_dict
  USE physics_generic
  USE geometry_generic
  USE sources_planetcooling
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: source_name = "planetheating"
!!$  LOGICAL, PARAMETER :: DEBUG = .TRUE.
!!  LOGICAL, PARAMETER :: DEBUG = .FALSE.
 REAL                        :: MDISK
 REAL, DIMENSION(:,:),POINTER::eff_Alpha,Te,H_tau1
! REAL, DIMENSION(:,:),POINTER::Te,H_tau1,H_tau2
 !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! methods
       InitSources_planetheating, &
       ExternalSources_planetheating, &
       CloseSources_planetheating,&
       UpdatePlanetHeating,&
       CalcTimestep_PlanetHeating,&
        CalcTimestep_HeatingCooling
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources_planetheating(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Dict_TYP),POINTER :: config,IO
    INTEGER           :: stype,i,j
    REAL              :: rho_dust
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "stype", stype)
    CALL InitSources(this,stype,source_name)
    ! some sanity checks
    IF (GetType(Physics).NE.EULER2D) &
         CALL Error(this,"InitSources_planetheating","physics not supported")
    IF (GetType(Physics%constants).NE.SI) &
         CALL Error(this,"InitSources_planetheating","only SI units supported")

    ALLOCATE(&
         this%Qstar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%T_s(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%cos1(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
         this%sin1(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
         STAT=err)
    IF (err.NE.0) CALL Error(this,"InitSources_planetheating","memory allocation failed")

    ! Courant number, i.e. safety factor for numerical stability
    CALL RequireKey(config, "cvis", 0.1)
    CALL GetAttr(config, "cvis", this%cvis)
   
    ! intensity of the planets star at 1 AU
    CALL RequireKey(config, "intensity")
    CALL GetAttr(config, "intensity", this%intensity)
    ! albedo of the planet
    CALL RequireKey(config, "albedo")
    CALL GetAttr(config, "albedo", this%albedo)
    ! distance planet-star 
    CALL RequireKey(config, "distance")
    CALL GetAttr(config, "distance", this%distance)
    ! molar mass
    CALL RequireKey(config, "mu")
    CALL GetAttr(config, "mu", this%mu)
    ! day-night omega
    CALL RequireKey(config, "omegasun",0.0)
    CALL GetAttr(config, "omegasun", this%omegasun)
    ! year
    CALL RequireKey(config, "year")
    CALL GetAttr(config, "year", this%year)
    ! gravitational acceleration
    CALL RequireKey(config, "gacc")
    CALL GetAttr(config, "gacc", this%gacc)
    ! ratio of specific heats
    CALL RequireKey(config, "gamma")
    CALL GetAttr(config, "gamma", this%gamma)
    ! beginning points 
    CALL RequireKey(config, "theta0", 0.0)
    CALL RequireKey(config, "phi0", 0.0)
    CALL GetAttr(config, "theta0", this%theta0)
    CALL GetAttr(config, "phi0", this%phi0)
 
    IF (GetType(Mesh%geometry).EQ.BIANGLESPHERICAL) THEN
    this%cos1(:,:,1) = COS(Mesh%bcenter(:,:,1))
!    this%cos1(:,:,2) = COS(Mesh%bcenter(:,:,2))
    this%sin1(:,:,1) = SIN(Mesh%bcenter(:,:,1))
!    this%sin1(:,:,2) = SIN(Mesh%bcenter(:,:,2))
    ELSE
    ! do nothing
    END IF

   ! Radius of the Planet
!     CALL RequireKey(config, "R_planet")
!     CALL GetAttr(config, "R_planet", this%R_planet)
    
!     CALL GetAttr(config, "dz", Mesh%dz)
    
    ! set initial time < 0
    this%time = -1.0

    ! initialize arrays
    this%Qstar(:,:)     = 0.0
    this%T_s(:,:)  = 0.0

    !initialise output
     CALL SetOutput(this,Mesh,Physics,config,IO)
  END SUBROUTINE InitSources_planetheating


  SUBROUTINE ExternalSources_planetheating(this,Mesh,Physics,time,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar,cvar,sterm
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,time,pvar,cvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    sterm(:,:,Physics%DENSITY) = 0.0
    sterm(:,:,Physics%XMOMENTUM) = 0.0
    sterm(:,:,Physics%YMOMENTUM) = 0.0


    CALL UpdatePlanetHeating(this,Mesh,Physics,time,pvar)
    !radiative heating by the central stars
    sterm(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%ENERGY) = &
         +this%Qstar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)
  END SUBROUTINE ExternalSources_planetheating

 
  SUBROUTINE CalcTimestep_planetheating(this,Mesh,Physics,time,pvar,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                      :: pvar
    REAL              :: time,dt, dt_new
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    REAL              :: invdt
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,time
    INTENT(INOUT)     :: this,Physics
    INTENT(OUT)       :: dt
    !------------------------------------------------------------------------!
    CALL  UpdatePlanetHeating(this,Mesh,Physics,time,pvar)
    ! maximum of inverse heating timescale t_star ~ P/Q_star
    invdt = MAXVAL(ABS(this%Qstar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) &
         / pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE)))  
  
    
    IF (invdt.GT.TINY(invdt)) THEN
       dt = this%cvis / invdt
    ELSE
       dt = HUGE(invdt)
    END IF
  END SUBROUTINE CalcTimestep_planetheating

  
 SUBROUTINE CalcTimestep_HeatingCooling(this,Mesh,Physics,time,pvar,dt,PLANET_COOLING)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Sources_TYP) :: this
   TYPE(Mesh_TYP)    :: Mesh
   TYPE(Physics_TYP) :: Physics
   REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                     :: pvar
   REAL              :: time,dt, dt_new
   INTEGER           :: i,j,PLANET_COOLING
   !------------------------------------------------------------------------!
   REAL              :: invdt
   TYPE(Sources_TYP), POINTER :: cooling,heating
   !------------------------------------------------------------------------!
   INTENT(IN)        :: Mesh,time
   INTENT(INOUT)     :: this,Physics
   INTENT(OUT)       :: dt
   !------------------------------------------------------------------------!
   CALL  UpdatePlanetHeating(this,Mesh,Physics,time,pvar)
   cooling=>GetSourcesPointer(Physics%sources,PLANET_COOLING)
   CALL  UpdatePlanetCooling(cooling,Mesh,Physics,time,pvar)
   ! -energy loss due to radiation processes
   ! +energy gain due to planet irradiation
   ! calculated in one source term to avoid small timesteps
   invdt = MAXVAL(ABS((this%Qstar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) &
          -this%Qcool(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX))&
        / pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE)))  
 
   
   IF (invdt.GT.TINY(invdt)) THEN
      dt = this%cvis / invdt
   ELSE
      dt = HUGE(invdt)
   END IF
   
 END SUBROUTINE CalcTimestep_HeatingCooling


  SUBROUTINE UpdatePlanetHeating(this,Mesh,Physics,time,pvar)
    USE functions
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    !------------------------------------------------------------------------!
    INTEGER, SAVE     :: counter = 0
    INTEGER           :: i,j,l
    REAL              :: Qstar,R_planet
    REAL              :: theta1,theta2,phi1,phi2,phi0,theta0
    REAL              :: T_0,T_s,distance_change
    REAL              :: da
    !------------------------------------------------------------------------!
    REAL, PARAMETER   :: AU      = 1.49597870691E+11    ! astronomical unit [m]
    REAL, PARAMETER    :: DAY      = 8.6400E+4        ! Day [sec]                    !
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,pvar
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! heating by the star
    IF (time.NE.this%time) THEN
      DO i=Mesh%IMIN,Mesh%IMAX
         DO j=Mesh%JMIN,Mesh%JMAX
              !--------------------------------------------------------------!
              phi1=Mesh%bcenter(i,j,2) + 2.*PI*this%omegasun*this%time
              ! transformation
              theta2 = ACOS(this%cos1(i,j,1)*COS(this%theta0*&
                       COS(2*PI*this%time/this%year)) + &
                       this%sin1(i,j,1)*SIN(this%theta0*&
                       COS(2*PI*this%time/this%year))*COS(phi1-this%phi0))
              
              ! modulo function to prevent phi2>2*PI,phi2<0
              phi2   = MODULO(ATAN2(this%sin1(i,j,1)*SIN(phi1-this%phi0),&
                       (-this%cos1(i,j,1)*SIN(this%theta0*COS(2.*PI*this%time/&
                       this%year)) + this%sin1(i,j,1)*COS(phi1-this%phi0)*&
                       COS(this%theta0*COS(2*PI*this%time/this%year)))),2*PI)

!             theta2 =  ACOS(this%cos1(i,j,1)*COS(this%theta0*&
!                       COS(2*PI*this%time/this%year)) - &
!                       this%sin1(i,j,1)*SIN(this%theta0*&
!                       COS(2*PI*this%time/this%year))*COS(phi1))
!                             
!                       ! modulo function to prevent phi2>2*PI,phi2<0
!             phi2    = MODULO(this%phi0 + &
!                       ATAN2(this%sin1(i,j,1)*SIN(phi1),&
!                       (this%cos1(i,j,1)*SIN(this%theta0)*COS(2*PI*this%time/&
!                       (this%year))+this%sin1(i,j,1)*COS(phi1)*&
!                       COS(this%theta0*COS(2*PI*this%time/this%year)))),2*PI)        


              ! calculate heating source                                     !
              !--------------------------------------------------------------!
              IF (phi2.LE.PI/2..OR.phi2.GE.3./2.*PI) THEN
              ! for sin,cos have a look at spherical trigonometry
              ! the albedo can contain things like scattering etc. (not implemented)

                   ! distance to the star              
                   distance_change=this%distance

                   Qstar = this%intensity*&
                   (AU/distance_change)**(2.)*(1.-(this%albedo))*SIN(theta2)*COS(phi2)
              ELSE
              Qstar = 0.0
              END IF
              this%Qstar(i,j) = Qstar  

           END DO         
        END DO 
       this%time=time
    END IF
!     print *, sum(this%Energy(:,:)), sum(this%Qstar(:,:))*this%time, &
!                 sum(this%Energy(:,:))- sum(this%Qstar(:,:)*this%time)
  END SUBROUTINE UpdatePlanetHeating
 
   

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

    !haeting source term
    valwrite = 0 
    IF (HasKey(config, "output/Qstar")) CALL GetAttr(config, "output/Qstar", valwrite)
    IF (valwrite .EQ. 1) THEN
         CALL AddField(IO, &
               "Qstar", &
               this%Qstar,&
               Dict("name" / "Qstar"))
    END IF 
  END SUBROUTINE SetOutput

  SUBROUTINE CloseSources_planetheating(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%Qstar,this%T_s,this%cos1,this%sin1)

    CALL CloseSources(this)
  END SUBROUTINE CloseSources_planetheating

END MODULE sources_planetheating
