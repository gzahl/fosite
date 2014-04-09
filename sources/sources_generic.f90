!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_generic.f90                                               #
!#                                                                           #
!# Copyright (C) 2007-2013                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
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
!> \addtogroup sources
!! - general parameters of sources group as key-values
!! \key{stype,INTEGER,Type of source}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Björn Sperling
!!
!! \brief generic source terms module providing functionaly common
!! to all source terms
!!
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_generic
  USE mesh_common, ONLY : Mesh_TYP
  USE boundary_common, ONLY : Boundary_TYP
  USE timedisc_common, ONLY : Timedisc_TYP
USE sources_c_accel, InitSources_common => InitSources, &
       CloseSources_common => CloseSources
  USE sources_diskthomson
  USE sources_viscosity
  USE sources_wave_damping
  USE sources_cooling
  USE sources_rotframe
  USE sources_sgs
  USE sources_diskcooling
  USE sources_planetheating
  USE sources_planetcooling
  USE sources_forcing
  USE gravity_generic
  USE physics_generic, GeometricalSources_Physics => GeometricalSources, &
       ExternalSources_Physics => ExternalSources
  USE fluxes_generic
  USE mesh_generic, ONLY : PI
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! tempory storage for source terms
  REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: temp_sterm
  ! flags for source terms
  INTEGER, PARAMETER :: GRAVITY          = 1
  INTEGER, PARAMETER :: DISK_THOMSON     = 2
  INTEGER, PARAMETER :: VISCOSITY        = 3
  INTEGER, PARAMETER :: C_ACCEL          = 4
  INTEGER, PARAMETER :: COOLING          = 5
  INTEGER, PARAMETER :: ROTATING_FRAME   = 20
  INTEGER, PARAMETER :: SGS              = 23
  INTEGER, PARAMETER :: DISK_COOLING     = 24
  INTEGER, PARAMETER :: WAVE_DAMPING     = 25
  INTEGER, PARAMETER :: FORCING          = 26
  INTEGER, PARAMETER :: PLANET_HEATING   = 27
  INTEGER, PARAMETER :: PLANET_COOLING   = 28
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! constants
       GRAVITY, DISK_THOMSON, VISCOSITY, C_ACCEL, COOLING, &
       ROTATING_FRAME, SGS, DISK_COOLING, WAVE_DAMPING, FORCING, &
       POINTMASS, POINTMASS_BINARY, MONOPOL, &
       NEWTON, WIITA,  &
       MOLECULAR, ALPHA, BETA, PRINGLE, ALPHA_ALT, &
       MULTIGRID, SPECTRAL, POTENTIAL, &
       RED_BLACK_GAUSS_SEIDEL,BLOCK_GAUSS_SEIDEL,GAUSS_SEIDEL, &
       SPHERMULTEXPAN, CYLINMULTEXPAN, &
       PLANET_HEATING, PLANET_COOLING, & 
       ! methods
       InitSources, &
       MallocSources, &
       CloseSources, &
       GeometricalSources, &
       ExternalSources, &
       CalcTimestep, &
       GetSourcesPointer, &
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

  SUBROUTINE InitSources(list,Mesh,Fluxes,Physics,Boundary,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: list
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Physics_TYP) :: Physics
    TYPE(Boundary_TYP), DIMENSION(4) :: Boundary
    TYPE(Dict_TYP),POINTER :: config,IO
    INTEGER           :: stype
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: sp
    TYPE(Dict_TYP),POINTER :: dir,src,IOsrc
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes
    INTENT(INOUT)     :: Boundary,Physics
    !------------------------------------------------------------------------!
    IF (.NOT.Initialized(Physics).OR..NOT.Initialized(Mesh)) &
         CALL Error(list,"InitSources","physics and/or mesh module uninitialized")
    ! allocate common memory for all sources
    IF (.NOT.ALLOCATED(temp_sterm)) THEN
       CALL MallocSources(list,Mesh,Physics)
    END IF
    dir => config
    DO WHILE(ASSOCIATED(dir))
        CALL GetAttr(config, GetKey(dir), src)
        CALL GetAttr(IO, GetKey(dir), IOsrc)
        dir => GetNext(dir)

        CALL RequireKey(src, "stype")
        CALL GetAttr(src, "stype", stype)

        SELECT CASE(stype)
        CASE(GRAVITY)
           ! gravity modules
           CALL InitGravity(list,Mesh,Fluxes,Physics,Boundary,stype,src,IOsrc)
        CASE(DISK_THOMSON)
           ! radiational acceleration due to Thomson scattering
           ! of accretion disk radiation
           CALL InitSources_diskthomson(list,Mesh,Physics,src,IOsrc)
        CASE(VISCOSITY)
           ! viscous diffusion and heating
           CALL InitSources_viscosity(list,Mesh,Physics,Fluxes,src,IOsrc)
        CASE(C_ACCEL)
           ! constant acceleration in x- and y-direction
           CALL InitSources_c_accel(list,Mesh,Physics,src)
        CASE(WAVE_DAMPING)
           ! wave damping for planet eu experiments
           CALL InitSources_wave_damping(list,Mesh,Physics,src)
        CASE(COOLING)
           ! simple cooling function
           CALL InitSources_cooling(list,Mesh,Physics,src)
        CASE(ROTATING_FRAME)
           ! inertial forces due to rotating reference frame
           CALL InitSources_rotframe(list,Mesh,Physics,src)
        CASE(SGS)
           CALL InitSources_sgs(list,Mesh,Physics,Fluxes,src,IOsrc)
        CASE(DISK_COOLING)
           CALL InitSources_diskcooling(list,Mesh,Physics,src,IOsrc)
        CASE(FORCING)
           CALL InitSources_forcing(list,Mesh,Physics,Fluxes,src,IOsrc)
        CASE(PLANET_COOLING)
          CALL InitSources_planetcooling(list,Mesh,Physics,src,IOsrc)
        CASE(PLANET_HEATING)
          CALL InitSources_planetheating(list,Mesh,Physics,src,IOsrc)
        CASE DEFAULT
           CALL Error(list,"InitSources", "unknown source term")
        END SELECT
        ! print some information
        IF (ASSOCIATED(list)) THEN
           CALL Info(list, " SOURCES--> source term:       " // GetName(list))
           ! print setup information of the individual source terms
           CALL InfoSources(list)
        END IF
    END DO
  END SUBROUTINE InitSources


  SUBROUTINE MallocSources(list,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: list
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    !------------------------------------------------------------------------!
    ! temporay storage
    ALLOCATE(temp_sterm(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
         STAT=err)
    IF (err.NE.0) CALL Error(list, "MallocSources_generic", "Unable allocate memory!")
  END SUBROUTINE MallocSources


  SUBROUTINE InfoSources(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    !------------------------------------------------------------------------!
    IF (ASSOCIATED(this)) THEN
!CDIR IEXPAND
       SELECT CASE(GetType(this))
       CASE(C_ACCEL,COOLING,SGS,DISK_COOLING,GRAVITY,&
         PLANET_COOLING,PLANET_HEATING)
          ! do nothing
       CASE(DISK_THOMSON)
          CALL InfoSources_diskthomson(this)
       CASE(VISCOSITY)
          CALL InfoSources_viscosity(this)
       CASE(ROTATING_FRAME)
          CALL InfoSources_rotframe(this)
       CASE(FORCING)
          CALL InfoSources_forcing(this)
       END SELECT
    END IF
  END SUBROUTINE InfoSources


  SUBROUTINE GeometricalSources(Physics,Mesh,Fluxes,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)  :: Physics
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Fluxes_TYP)   :: Fluxes
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
         :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,pvar,cvar
    INTENT(INOUT)     :: Physics
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! calculate geometrical sources depending on the integration rule
!CDIR IEXPAND
    SELECT CASE(GetType(Fluxes))
    CASE(MIDPOINT)
       ! use center values for midpoint rule
       CALL GeometricalSources_physics(Physics,Mesh,pvar,cvar,sterm)
    CASE(TRAPEZOIDAL)
       ! use reconstructed corner values for trapezoidal rule
       CALL GeometricalSources_physics(Physics,Mesh,Fluxes%prim,Fluxes%cons,sterm)
    END SELECT
  END SUBROUTINE GeometricalSources


  SUBROUTINE ExternalSources(this,Mesh,Fluxes,Physics,time,dt,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Physics_TYP) :: Physics
    REAL              :: time,dt
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: cvar,pvar,sterm
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: srcptr
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,time,pvar,cvar
    INTENT(INOUT)     :: Physics
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! reset sterm
    sterm(:,:,:) = 0.
    ! go through all source terms in the list
    srcptr => this
    DO WHILE (ASSOCIATED(srcptr))
       ! call specific subroutine

!CDIR IEXPAND
       SELECT CASE(GetType(srcptr))
       CASE(GRAVITY)
          CALL GravitySources(srcptr,Mesh,Physics,Fluxes,time,dt,pvar,cvar,temp_sterm)
       CASE(DISK_THOMSON)
          CALL ExternalSources_diskthomson(srcptr,Mesh,Physics,pvar,cvar,temp_sterm)
       CASE(VISCOSITY)
          CALL ExternalSources_viscosity(srcptr,Mesh,Physics,Fluxes,time,pvar,cvar,temp_sterm)
       CASE(C_ACCEL)
          CALL ExternalSources_c_accel(srcptr,Mesh,Physics,pvar,cvar,temp_sterm)
       CASE(WAVE_DAMPING)
          CALL ExternalSources_wave_damping(srcptr,Mesh,Physics,time,dt,pvar,cvar,temp_sterm)
       CASE(COOLING)
          CALL ExternalSources_cooling(srcptr,Mesh,Physics,time,pvar,cvar,temp_sterm)
       CASE(ROTATING_FRAME)
          CALL ExternalSources_rotframe(srcptr,Mesh,Physics,pvar,cvar,temp_sterm)
       CASE(SGS)
          CALL ExternalSources_sgs(srcptr,Mesh,Physics,pvar,cvar,temp_sterm)
       CASE(DISK_COOLING)
          CALL ExternalSources_diskcooling(srcptr,Mesh,Physics,Fluxes,time,pvar,cvar,temp_sterm)
       CASE(FORCING)
          CALL ExternalSources_forcing(srcptr,Mesh,Physics,time,dt,pvar,cvar,temp_sterm)
       CASE(PLANET_COOLING)
          CALL ExternalSources_planetcooling(srcptr,Mesh,Physics,time,pvar,cvar,temp_sterm)
       CASE(PLANET_HEATING)
          CALL ExternalSources_planetheating(srcptr,Mesh,Physics,time,pvar,cvar,temp_sterm)
       CASE DEFAULT
          CALL Error(srcptr,"ExternalSources", "unknown source term")
       END SELECT

       ! add to the sources
       sterm(:,:,:) = sterm(:,:,:) + temp_sterm(:,:,:)
       ! next source term
       srcptr => srcptr%next
    END DO
    ! reset ghost cell data
    sterm(Mesh%IGMIN:Mesh%IMIN-1,:,:) = 0.0
    sterm(Mesh%IMAX+1:Mesh%IGMAX,:,:) = 0.0
    sterm(:,Mesh%JGMIN:Mesh%JMIN-1,:) = 0.0
    sterm(:,Mesh%JMAX+1:Mesh%JGMAX,:) = 0.0
  END SUBROUTINE ExternalSources


  SUBROUTINE CalcTimestep(this,Mesh,Physics,Fluxes,time,pvar,cvar,dt,dtcause)
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
    INTEGER           :: dtcause
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: srcptr
    REAL              :: dt_new
    INTEGER           :: hc
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,time,pvar,cvar
    INTENT(INOUT)     :: dt,dtcause,Physics
    !------------------------------------------------------------------------!
    
    ! go through all source terms in the list
    srcptr => this
    DO WHILE(ASSOCIATED(srcptr))
       ! call specific subroutine
!CDIR IEXPAND
       SELECT CASE(GetType(srcptr))
       CASE(DISK_THOMSON,C_ACCEL,ROTATING_FRAME,WAVE_DAMPING,GRAVITY,&
         PLANET_HEATING)
          ! do nothing
          dt_new = dt
       CASE(VISCOSITY)
          CALL CalcTimestep_viscosity(srcptr,Mesh,Physics,Fluxes,time,pvar,cvar,dt_new)
       CASE(COOLING)
          CALL CalcTimestep_cooling(srcptr,Mesh,Physics,time,pvar,dt_new)
       CASE(DISK_COOLING)
          CALL CalcTimestep_diskcooling(srcptr,Mesh,Physics,Fluxes,time,pvar,dt_new)
       CASE(SGS)
          CALL CalcTimestep_sgs(srcptr,Mesh,Physics,time,pvar,cvar,dt_new)
       CASE(FORCING)
          CALL CalcTimestep_forcing(srcptr,Mesh,Physics,pvar,cvar,dt_new)
       CASE(PLANET_COOLING)
          CALL CalcTimestep_planetcooling(srcptr,Mesh,Physics,time,pvar,dt_new)   
       CASE DEFAULT
          CALL Error(srcptr,"CalcTimestep", "unknown source term")
       END SELECT
       ! who was it?
       IF (dt_new .LT. dt) dtcause=GetType(srcptr)
       dt = MIN(dt,dt_new)
!        print*, 'dt', dt
       ! next source term
       srcptr => srcptr%next
    END DO    
  END SUBROUTINE CalcTimestep


  SUBROUTINE CloseSources(this,Fluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Fluxes_TYP)  :: Fluxes
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: srcptr
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Fluxes
    !------------------------------------------------------------------------!
    ! call deallocation procedures for all source terms
    DO
       srcptr => this
       IF (.NOT.ASSOCIATED(srcptr)) EXIT
       this => srcptr%next
       IF (.NOT.Initialized(srcptr)) &
            CALL Error(this,"CloseSources","not initialized")
       ! call specific deconstructor
!CDIR IEXPAND
       SELECT CASE(GetType(srcptr))
       CASE(GRAVITY)
          CALL CloseGravity(srcptr)
       CASE(DISK_THOMSON)
          CALL CloseSources_diskthomson(srcptr)
       CASE(VISCOSITY)
          CALL CloseSources_viscosity(srcptr)
       CASE(C_ACCEL)
          CALL CloseSources_c_accel(srcptr,Fluxes)
       CASE(WAVE_DAMPING)
          CALL CloseSources_wave_damping(srcptr,Fluxes)
       CASE(COOLING)
          CALL CloseSources_cooling(srcptr)
       CASE(ROTATING_FRAME)
          CALL CloseSources_rotframe(srcptr)
       CASE(SGS)
          CALL CloseSources_sgs(srcptr)
       CASE(DISK_COOLING)
          CALL CloseSources_diskcooling(srcptr)
       CASE(FORCING)
          CALL CloseSources_forcing(srcptr)
       CASE(PLANET_COOLING)
          CALL CloseSources_planetcooling(srcptr)
       CASE(PLANET_HEATING)
          CALL CloseSources_planetheating(srcptr)
       END SELECT
       ! deallocate source term structure
       DEALLOCATE(srcptr)
    END DO
    ! release temporary storage
    DEALLOCATE(temp_sterm)
  END SUBROUTINE CloseSources

END MODULE sources_generic
