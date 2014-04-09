!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_generic.f90                                               #
!#                                                                           #
!# Copyright (C) 2007-2012                                                   #
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
!----------------------------------------------------------------------------!
! generic source terms module providing functionaly common
! to all source terms
!----------------------------------------------------------------------------!
MODULE sources_generic
  USE mesh_common, ONLY : Mesh_TYP
  USE boundary_common, ONLY : Boundary_TYP
  USE timedisc_common, ONLY : Timedisc_TYP
  USE sources_pointmass, InitSources_common => InitSources, &
       CloseSources_common => CloseSources
  USE sources_diskthomson
  USE sources_viscosity
  USE sources_c_accel
  USE sources_cooling
  USE poisson_generic
  USE physics_generic, GeometricalSources_Physics => GeometricalSources, &
       ExternalSources_Physics => ExternalSources
  USE fluxes_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! tempory storage for source terms
  REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: temp_sterm
  ! flags for source terms
  INTEGER, PARAMETER :: POINTMASS        = 1
  INTEGER, PARAMETER :: DISK_THOMSON     = 2
  INTEGER, PARAMETER :: VISCOSITY        = 3
  INTEGER, PARAMETER :: C_ACCEL          = 4
  INTEGER, PARAMETER :: COOLING          = 5
  INTEGER, PARAMETER :: POISSON          = 6
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! constants
       POINTMASS, DISK_THOMSON, VISCOSITY, C_ACCEL, COOLING, POISSON, &
       NEWTON, WIITA, &
       MOLECULAR, ALPHA, BETA, PRINGLE, &
       MULTIGRID, &
       RED_BLACK_GAUSS_SEIDEL,BLOCK_GAUSS_SEIDEL,GAUSS_SEIDEL, &
       SPHERMULTEXPAN, CYLINMULTEXPAN, &
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

  SUBROUTINE InitSources(list,Mesh,Fluxes,Physics,Boundary,stype,potential,vismodel, &
       mass,mdot,rin,rout,dynconst,bulkconst,cvis,xaccel,yaccel,solver,maxresidnorm, &
       maxmult,bndrytype,relaxtype,npre,npost,minres,nmaxcycle)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: list
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Physics_TYP) :: Physics
    TYPE(Boundary_TYP), DIMENSION(4) :: Boundary
    INTEGER           :: stype
    INTEGER, OPTIONAL :: potential,vismodel,solver,maxmult,bndrytype,relaxtype,&
                         npre,npost,minres,nmaxcycle
    REAL, OPTIONAL    :: mass,mdot,rin,rout,dynconst,bulkconst,cvis, &
                         xaccel,yaccel,maxresidnorm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,Physics,Boundary,stype,potential,vismodel, &
                         mass,mdot,rin,rout,dynconst,bulkconst,cvis,xaccel,yaccel, &
                         solver,maxresidnorm,maxmult,bndrytype,relaxtype,&
                         npre,npost,minres,nmaxcycle
    !------------------------------------------------------------------------!
    IF (.NOT.Initialized(Physics).OR..NOT.Initialized(Mesh)) &
         CALL Error(list,"InitSources","physics and/or mesh module uninitialized")
    ! allocate common memory for all sources
    IF (.NOT.ALLOCATED(temp_sterm)) THEN
       CALL MallocSources(list,Mesh,Physics)
    END IF

    SELECT CASE(stype)
    CASE(POINTMASS)
       ! gravitational acceleration due to point mass
       CALL InitSources_pointmass(list,Mesh,Physics,stype,potential,mass,cvis)
    CASE(DISK_THOMSON)
       ! radiational acceleration due to Thomson scattering
       ! of accretion disk radiation
       CALL InitSources_diskthomson(list,Mesh,Physics,stype,mass,mdot,rin,rout)
    CASE(VISCOSITY)
       ! viscous diffusion and heating
       CALL InitSources_viscosity(list,Mesh,Physics,Fluxes,stype,vismodel, &
            dynconst,bulkconst,cvis)
    CASE(C_ACCEL)
       ! constant acceleration in x- and y-direction
       CALL InitSources_c_accel(list,Mesh,Physics,stype,xaccel,yaccel)
    CASE(COOLING)
       ! simple cooling function
       CALL InitSources_cooling(list,Mesh,Physics,stype,cvis)
    CASE(POISSON)
       ! poisson solver to compute gravitational potenial of
       ! the density distribution
       CALL InitPoisson(list,Mesh,Physics,Boundary,stype,solver,maxresidnorm,&
            maxmult,bndrytype,relaxtype,npre,npost,minres,nmaxcycle)
    CASE DEFAULT
       CALL Error(list,"InitSources", "unknown source term")
    END SELECT

    ! print some information
    IF (ASSOCIATED(list)) THEN
       CALL Info(list, " SOURCES--> source term:       " // GetName(list))
       ! print setup information of the individual source terms
       CALL InfoSources(list)
    END IF
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
       CASE(C_ACCEL,COOLING,POISSON)
          ! do nothing
       CASE(POINTMASS)
          CALL InfoSources_pointmass(this)
       CASE(DISK_THOMSON)
          CALL InfoSources_diskthomson(this)
       CASE(VISCOSITY)
          CALL InfoSources_viscosity(this)
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


  SUBROUTINE ExternalSources(this,Mesh,Fluxes,Physics,time,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
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
       CASE(POINTMASS)
          CALL ExternalSources_pointmass(srcptr,Mesh,Physics,Fluxes,pvar,cvar,temp_sterm)
       CASE(DISK_THOMSON)
          CALL ExternalSources_diskthomson(srcptr,Mesh,Physics,pvar,cvar,temp_sterm)
       CASE(VISCOSITY)
          CALL ExternalSources_viscosity(srcptr,Mesh,Physics,time,pvar,cvar,temp_sterm)
       CASE(C_ACCEL)
          CALL ExternalSources_c_accel(srcptr,Mesh,Physics,pvar,cvar,temp_sterm)
       CASE(COOLING)
          CALL ExternalSources_cooling(srcptr,Mesh,Physics,time,pvar,cvar,temp_sterm)
       CASE(POISSON)
          CALL PoissonSource(srcptr%poisson,Mesh,Physics,pvar,cvar,temp_sterm)
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


  SUBROUTINE CalcTimestep(this,Mesh,Physics,time,pvar,cvar,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar,cvar
    REAL              :: dt
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: srcptr
    REAL              :: dt_new
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,time,pvar,cvar
    INTENT(INOUT)     :: dt,Physics
    !------------------------------------------------------------------------!
    ! go through all source terms in the list
    srcptr => this
    DO
       IF (.NOT.ASSOCIATED(srcptr)) EXIT
       ! call specific subroutine
!CDIR IEXPAND
       SELECT CASE(GetType(srcptr))
       CASE(DISK_THOMSON,C_ACCEL,POISSON)
          ! do nothing
          dt_new = dt
       CASE(POINTMASS)
          CALL CalcTimestep_pointmass(srcptr,Mesh,Physics,pvar,cvar,dt_new)
       CASE(VISCOSITY)
          CALL CalcTimestep_viscosity(srcptr,Mesh,Physics,time,pvar,cvar,dt_new)
       CASE(COOLING)
          CALL CalcTimestep_cooling(srcptr,Mesh,Physics,time,pvar,dt_new)
       CASE DEFAULT
          CALL Error(srcptr,"CalcTimestep", "unknown source term")
       END SELECT
       dt = MIN(dt,dt_new)
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
       CASE(POINTMASS)
          CALL CloseSources_pointmass(srcptr)
       CASE(DISK_THOMSON)
          CALL CloseSources_diskthomson(srcptr)
       CASE(VISCOSITY)
          CALL CloseSources_viscosity(srcptr)
       CASE(C_ACCEL)
          CALL CloseSources_c_accel(srcptr,Fluxes)
       CASE(COOLING)
          CALL CloseSources_cooling(srcptr)
       CASE(POISSON)
          CALL ClosePoisson(srcptr%poisson)
       END SELECT
       ! deallocate source term structure
       DEALLOCATE(srcptr)
    END DO
    ! release temporary storage
    DEALLOCATE(temp_sterm)
  END SUBROUTINE CloseSources

END MODULE sources_generic
