!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_planetcooling.f90                                         #
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
!> \author Tobias Illenseer
!!
!! \brief source terms module for gray cooling of thin planet atmospheres
!!
!! \warning use SI units
!!
!! \extends sources_common
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_planetcooling
  USE constants_common, ONLY : KE
  USE gravity_pointmass
  USE gravity_binary
  USE sources_common
  USE physics_generic
  USE mesh_generic
  USE common_dict
  USE functions
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: source_name = "planetcooling"

  REAL, PARAMETER :: Texp(8) = (/ 2.0, -7.0, -0.5, -24.0, 3.0, 10.0, -2.5, 0.0 /)
  REAL, PARAMETER :: rexp(8) = (/ 0.0, 0.0, 0.0, 1.0, 2./3., 1./3., 1.0, 0.0 /)
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! methods
       InitSources_planetcooling, &
       ExternalSources_planetcooling, &
       CalcTimestep_planetcooling, &
       CloseSources_planetcooling, &
       UpdatePlanetCooling
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources_planetcooling(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Dict_TYP),POINTER :: config,IO
    INTEGER           :: stype
    REAL              :: intensity 
    REAL, PARAMETER   :: AU      = 1.49597870691E+11    ! astronomical unit [m]
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "stype", stype)
    CALL InitSources(this,stype,source_name)
    ! some sanity checks
    IF (GetType(Physics).NE.EULER2D) &
         CALL Error(this,"InitSources_planetcooling","physics not supported")
    IF (GetType(Physics%constants).NE.SI) &
         CALL Error(this,"InitSources_planetcooling","only SI units supported")

    ALLOCATE(this%cs(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX),&
         this%tmp(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX),&
         this%Qcool(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%T_s(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX),&
         this%RHO_s(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%P_s(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         STAT=err)
    IF (err.NE.0) CALL Error(this,"InitSources_planetcooling","memory allocation failed")

    ! Courant number, i.e. safety factor for numerical stability
    CALL RequireKey(config, "cvis", 0.9)
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
    ! equilibrium temperature
    CALL RequireKey(config, "T_0")
    CALL GetAttr(config,"T_0", this%T_0)
    ! gravitational acceleration
    CALL RequireKey(config, "gacc")
    CALL GetAttr(config, "gacc", this%gacc)
    ! ratio of specific heats
    CALL RequireKey(config, "gamma")
    CALL GetAttr(config, "gamma", this%gamma)

    ! set initial time < 0
    this%time = -1.0

    ! initialize arrays
    this%cs(:,:)    = 1.0
    this%tmp(:,:)   = 0.0
    this%Qcool(:,:) = 0.0
    this%T_s(:,:)   = this%T_0
    this%RHO_s(:,:) = 0.0
    this%P_s(:,:)   = 0.0
    
    ! initial calculation of optical thickness
    intensity =  (this%intensity/4.)*(AU/this%distance)**(2.)
    this%tau_inf= (((1.-this%albedo)*intensity/Physics%Constants%SB)&
             **(-0.25)*EXP(LnGamma(4.*(Physics%Constants%RG/Physics%Constants%RG)&
             *(this%gamma-1.)/this%gamma))**(0.25)*this%T_0)&
             **(this%gamma/(this%gamma-1.)) 

    print *, this%tau_inf

    CALL SetOutput(this,Mesh,Physics,config,IO)
  END SUBROUTINE InitSources_planetcooling


  SUBROUTINE ExternalSources_planetcooling(this,Mesh,Physics,time,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,time,pvar,cvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    sterm(:,:,Physics%DENSITY) = 0.0
    sterm(:,:,Physics%XMOMENTUM) = 0.0
    sterm(:,:,Physics%YMOMENTUM) = 0.0

    CALL UpdatePlanetCooling(this,Mesh,Physics,time,pvar)
    ! enery loss due to radiation processes
    sterm(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%ENERGY) = &
         -this%Qcool(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)
  END SUBROUTINE ExternalSources_planetcooling

 
  SUBROUTINE CalcTimestep_planetcooling(this,Mesh,Physics,time,pvar,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                      :: pvar
    REAL              :: time,dt
    !------------------------------------------------------------------------!
    REAL              :: invdt
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,time
    INTENT(INOUT)     :: this,Physics
    INTENT(OUT)       :: dt
    !------------------------------------------------------------------------!
    CALL UpdatePlanetCooling(this,Mesh,Physics,time,pvar)
    ! maximum of inverse cooling timescale t_cool ~ P/Q_cool
    invdt = MAXVAL(ABS(this%Qcool(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) &
         / pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE)))
    IF (invdt.GT.TINY(invdt)) THEN
       dt = this%cvis / invdt
    ELSE
       dt = HUGE(invdt)
    END IF
  END SUBROUTINE CalcTimestep_planetcooling


  SUBROUTINE UpdatePlanetCooling(this,Mesh,Physics,time,pvar)
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
    INTEGER           :: i,j
    REAL              :: T_0,T_s,Qcool,abs_vel
    REAL, PARAMETER   :: AU      = 1.49597870691E+11    ! astronomical unit [m]
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,pvar
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    abs_vel = 0.0
    ! calculation of cooling source  
    IF (time.NE.this%time) THEN
       DO j=Mesh%JMIN,Mesh%JMAX
          DO i=Mesh%IMIN,Mesh%IMAX
             ! calculate surface density and pressure
             ! (irrelevant for calculation)
             this%P_s(i,j)   = pvar(i,j,Physics%DENSITY)*this%gacc
             this%RHO_s(i,j) = this%P_s(i,j)*this%mu/(Physics%Constants%RG*this%T_s(i,j))
             
             !--------------------------------------------------------------!
             ! calculate cooling source                                     ! 
             !--------------------------------------------------------------!
             ! basic concept of the model in "Principles of Planetary
             ! Climate" from Raymond T. Pierrehumbert 199ff

             ! calculate planet-surface temperature from initial conditions
             T_s   = pvar(i,j,Physics%PRESSURE)*this%mu*&
                     (1.+(this%gamma-1.)/this%gamma)&  
                     /(Physics%Constants%RG*pvar(i,j,Physics%DENSITY)) 
             this%T_s(i,j) = T_s
             
             ! cooling source
             Qcool = Physics%Constants%SB*this%tau_inf**&
                     (-4.*(this%gamma-1.)/this%gamma)&
                     *EXP(LnGamma(4.*(Physics%Constants%RG/Physics%Constants%RG)*(this%gamma-1.)/this%gamma))&
                     *T_s**(4.)
             
             this%Qcool(i,j)=Qcool 
             
             abs_vel=abs_vel + SQRT(pvar(i,j,Physics%XVELOCITY)**2. + &
                                    pvar(i,j,Physics%YVELOCITY)**2.)
          END DO
       END DO
!       print *, abs_vel
       this%time=time
    END IF
  END SUBROUTINE UpdatePlanetCooling
  
  
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

    !cooling source term
    valwrite = 0 
    IF (HasKey(config, "output/Qcool")) CALL GetAttr(config, "output/Qcool", valwrite)
    IF (valwrite .EQ. 1) THEN
         CALL AddField(IO, &
               "Qcool", &
               this%Qcool,&
               Dict("name" / "Qcool"))
    END IF 
    
    ! temperature
    valwrite = 0 
    IF (HasKey(config, "output/T_s")) CALL GetAttr(config, "output/T_s", valwrite)
    IF (valwrite .EQ. 1) THEN
         CALL AddField(IO, &
               "T_s", &
               this%T_s,&
               Dict("name" / "T_s"))
    END IF 
    
    ! surface pressure
    valwrite = 0 
    IF (HasKey(config, "output/P_s")) CALL GetAttr(config, "output/P_s", valwrite)
    IF (valwrite .EQ. 1) THEN
         CALL AddField(IO, &
               "P_s", &
               this%P_s,&
               Dict("name" / "P_s"))
    END IF 

    ! surface density
    valwrite = 0 
    IF (HasKey(config, "output/RHO_s")) CALL GetAttr(config, "output/RHO_s", valwrite)
    IF (valwrite .EQ. 1) THEN
         CALL AddField(IO, &
               "RHO_s", &
               this%RHO_s,&
               Dict("name" / "RHO_s"))
    END IF 
  END SUBROUTINE SetOutput

  SUBROUTINE CloseSources_planetcooling(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%cs,this%tmp,this%Qcool,this%T_s,this%RHO_s,this%P_s)
    CALL CloseSources(this)
  END SUBROUTINE CloseSources_planetcooling


END MODULE sources_planetcooling
