!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_generic.f90                                               #
!#                                                                           #
!# Copyright (C) 2007-2008                                                   #
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
! generic source terms module providing functionaly common
! to all source terms
!----------------------------------------------------------------------------!
MODULE sources_generic
  USE sources_pointmass, InitSources_all => InitSources
  USE sources_diskthomson
  USE sources_viscosity
  USE physics_generic, GeometricalSources_Physics => GeometricalSources, &
       ExternalSources_Physics => ExternalSources
  USE fluxes_generic
  USE mesh_common, ONLY : Mesh_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! tempory storage for source terms
  REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: temp_sterm
  ! flags for source terms
  INTEGER, PARAMETER :: POINTMASS    = 1
  INTEGER, PARAMETER :: DISK_THOMSON = 2
  INTEGER, PARAMETER :: VISCOSITY    = 3
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! constants
       POINTMASS, DISK_THOMSON, VISCOSITY, &
       NEWTON, WIITA, &
       MOLECULAR, ALPHA, BETA, &
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
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources(list,Mesh,Fluxes,Physics,stype,potential,vismodel, &
       mass,mdot,rin,rout,dynconst,bulkconst,cvis)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: list
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Physics_TYP) :: Physics
    INTEGER           :: stype
    INTEGER, OPTIONAL :: potential,vismodel
    REAL, OPTIONAL    :: mass,mdot,rin,rout,dynconst,bulkconst,cvis
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: srcptr
    INTEGER           :: potential_def,vismodel_def
    REAL              :: mass_def,mdot_def,rin_def,rout_def,dynconst_def, &
                         bulkconst_def,cvis_def
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,Physics,stype,potential,vismodel,mass, &
                         mdot,rin,rout,dynconst,bulkconst,cvis
    !------------------------------------------------------------------------!
    ! allocate common memory for all sources
    IF (.NOT.ALLOCATED(temp_sterm)) THEN
       CALL MallocSources(list,Mesh,Physics)
    END IF

    SELECT CASE(stype)
    CASE(POINTMASS)
       ! default central mass
       IF (PRESENT(mass)) THEN
          mass_def = mass
       ELSE
          mass_def = 1.0
       END IF
       ! type of the potential
       IF (PRESENT(potential)) THEN
          potential_def=potential
       ELSE
          potential_def=NEWTON
       END IF
       CALL InitSources_pointmass(list,Mesh,Physics,stype,potential_def,mass_def)
    CASE(DISK_THOMSON)
       ! accretion rate
       IF (PRESENT(mdot)) THEN
          mdot_def = mdot
       ELSE
          mdot_def = 1.0
       END IF
       ! inner and outer disk radius
       IF (PRESENT(rin)) THEN
          rin_def = rin
       ELSE
          rin_def = 1.0
       END IF
       IF (PRESENT(rout)) THEN
          rout_def = rout
       ELSE
          rout_def = 2.0
       END IF
       CALL InitSources_diskthomson(list,Mesh,Physics,stype,mass_def,mdot_def, &
            rin_def,rout_def)
    CASE(VISCOSITY)
       ! viscosity model
       IF (PRESENT(vismodel)) THEN
          vismodel_def = vismodel
       ELSE
          vismodel_def = MOLECULAR
       END IF
       ! dynamic viscosity constant
       IF (PRESENT(dynconst)) THEN
          dynconst_def=dynconst
       ELSE
          dynconst_def=0.1
       END IF
       ! bulk viscosity constant (disabled by default)
       IF (PRESENT(bulkconst)) THEN
          bulkconst_def=bulkconst
       ELSE
          bulkconst_def=0.0
       END IF
       ! viscous Courant number
       IF (PRESENT(cvis)) THEN
          cvis_def=cvis
       ELSE
          cvis_def=0.5
       END IF
       ! this is necessary for alpha viscosity
       srcptr => GetSourcesPointer(list,POINTMASS)
       CALL InitSources_viscosity(list,Mesh,Physics,Fluxes,stype,vismodel_def, &
            dynconst_def,bulkconst_def,cvis_def,srcptr)
    CASE DEFAULT
       CALL Error(list,"InitSources", "unknown source term")
    END SELECT

    ! print some information
    IF (ASSOCIATED(list)) THEN
       CALL Info(list, " SOURCES--> source term:       " // GetName(list))
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
    SELECT CASE(GetType(Fluxes))
    CASE(MIDPOINT)
       ! use center values for midpoint rule
       CALL GeometricalSources_physics(Physics,Mesh,pvar,cvar,sterm)
    CASE(TRAPEZOIDAL)
       ! use reconstructed corner values for trapezoidal rule
       CALL GeometricalSources_physics(Physics,Mesh,Fluxes%prim,Fluxes%cons,sterm)
    END SELECT
  END SUBROUTINE GeometricalSources


  SUBROUTINE ExternalSources(this,Mesh,Fluxes,Physics,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: cvar,pvar,sterm
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: srcptr
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,Physics,pvar,cvar
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! reset sterm
    sterm(:,:,:) = 0.
    ! go through all source terms in the list
    srcptr => this
    DO
       IF (.NOT.ASSOCIATED(srcptr)) EXIT
       ! call specific subroutine
       SELECT CASE(GetType(srcptr))
       CASE(POINTMASS)
          CALL ExternalSources_pointmass(srcptr,Mesh,Physics,pvar,cvar,temp_sterm)
       CASE(DISK_THOMSON)
          CALL ExternalSources_diskthomson(srcptr,Mesh,Physics,pvar,cvar,temp_sterm)
       CASE(VISCOSITY)
          CALL ExternalSources_viscosity(srcptr,Mesh,Physics,pvar,cvar,temp_sterm)
       CASE DEFAULT
          CALL Error(srcptr,"ExternalSources", "unknown source term")
       END SELECT
       ! add to the sources
       sterm(:,:,:) = sterm(:,:,:) + temp_sterm(:,:,:)
       ! next source term
       srcptr => srcptr%next
    END DO    
  END SUBROUTINE ExternalSources


  SUBROUTINE CalcTimestep(this,Mesh,Physics,pvar,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    REAL              :: dt
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: srcptr
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,pvar
    INTENT(OUT)       :: dt
    !------------------------------------------------------------------------!
    ! go through all source terms in the list
    srcptr => this
    DO
       IF (.NOT.ASSOCIATED(srcptr)) EXIT
       ! call specific subroutine
       SELECT CASE(GetType(srcptr))
       CASE(POINTMASS,DISK_THOMSON)
          ! do nothing
       CASE(VISCOSITY)
          CALL CalcTimestep_viscosity(srcptr,Mesh,Physics,pvar,dt)
       CASE DEFAULT
          CALL Error(srcptr,"CalcTimestep", "unknown source term")
       END SELECT
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
       ! call specific deconstructor
       SELECT CASE(GetType(srcptr))
       CASE(POINTMASS)
          CALL CloseSources_pointmass(srcptr,Fluxes)
       CASE(VISCOSITY)
          CALL CloseSources_viscosity(srcptr)
       END SELECT
       ! deallocate source term structure
       DEALLOCATE(srcptr)
    END DO
    ! release temporary storage
    DEALLOCATE(temp_sterm)
  END SUBROUTINE CloseSources

END MODULE sources_generic
