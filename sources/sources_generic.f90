!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_generic.f90                                               #
!#                                                                           #
!# Copyright (C) 2007 Tobias Illenseer <tillense@ita.uni-heidelberg.de>      #
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
  USE physics_generic, GeometricalSources_Physics => GeometricalSources, &
       ExternalSources_Physics => ExternalSources
  USE fluxes_generic
  USE mesh_common, ONLY : Mesh_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE InitSources
     MODULE PROCEDURE InitSources_param1, InitSources_param4
  END INTERFACE
  ! tempory storage for source terms
  REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: temp_sterm
  ! flags for source terms
  INTEGER, PARAMETER :: POINTMASS    = 1
  INTEGER, PARAMETER :: DISK_THOMSON = 2
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! constants
       POINTMASS, DISK_THOMSON, &
       ! methods
       InitSources, &
       MallocSources, &
       GeometricalSources, &
       ExternalSources, &
       GetType, &
       GetName, &
       GetSourcesTimescale, &
       GetSourcesPointer, &
       CloseSources
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources_param1(list,Mesh,Fluxes,Physics,stype,sparam)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: list
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Physics_TYP) :: Physics
    INTEGER           :: stype
    REAL              :: sparam
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,Physics,stype,sparam
    !------------------------------------------------------------------------!
    ! allocate common memory for all sources
    IF (.NOT.ALLOCATED(temp_sterm)) THEN
       CALL MallocSources(list,Mesh,Physics)
    END IF

    SELECT CASE(stype)
    CASE(POINTMASS)
       CALL InitSources_pointmass(list,Mesh,Physics,stype,sparam)
    CASE DEFAULT
       PRINT *, "ERROR in InitSources: unknown source term"
       STOP
    END SELECT

    ! print some information
    IF (ASSOCIATED(list)) THEN
       PRINT "(A,A)", " SOURCES--> source term:       ", TRIM(GetName(list))
    END IF
  END SUBROUTINE InitSources_param1


  SUBROUTINE InitSources_param4(list,Mesh,Fluxes,Physics,stype,sp1,sp2,sp3,sp4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: list
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Physics_TYP) :: Physics
    INTEGER           :: stype
    REAL              :: sp1,sp2,sp3,sp4
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,Physics,stype,sp1,sp2,sp3,sp4
    !------------------------------------------------------------------------!
    ! allocate common memory for all sources
    IF (.NOT.ALLOCATED(temp_sterm)) THEN
       CALL MallocSources(list,Mesh,Physics)
    END IF

    SELECT CASE(stype)
    CASE(DISK_THOMSON)
       CALL InitSources_diskthomson(list,Mesh,Physics,stype,sp1,sp2,sp3,sp4)
    CASE DEFAULT
       PRINT *, "ERROR in InitSources: unknown source term"
       STOP
    END SELECT

    ! print some information
    IF (ASSOCIATED(list)) THEN
       PRINT "(A,A)", " SOURCES--> source term:       ", TRIM(GetName(list))
    END IF
  END SUBROUTINE InitSources_param4


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
    IF (err.NE.0) THEN
       PRINT *, "ERROR in MallocSources_generic: Unable allocate memory!"
       STOP
    END IF
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
    CASE DEFAULT
       PRINT *, "ERROR in GeometricalSources: unknown integration rule"
       STOP
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
       CASE DEFAULT
          PRINT *, "ERROR in CalculateSources: unknown source term"
          STOP
       END SELECT
       ! add to the sources
       sterm(:,:,:) = sterm(:,:,:) + temp_sterm(:,:,:)
       ! next source term
       srcptr => srcptr%next
    END DO    
  END SUBROUTINE ExternalSources


  FUNCTION GetSourcesTimescale(this,dtin) RESULT(dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    REAL, INTENT(IN) :: dtin
    REAL :: dt
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: srcptr
    !------------------------------------------------------------------------!
    srcptr => this
    dt = dtin
    DO
       IF (.NOT.ASSOCIATED(srcptr)) EXIT
       SELECT CASE(GetType(srcptr))
       CASE(POINTMASS,DISK_THOMSON)
          dt = MIN(dt,srcptr%dtmin)
       END SELECT
       srcptr => srcptr%next
    END DO    
  END FUNCTION GetSourcesTimescale
  

  FUNCTION GetSourcesPointer(this,stype) RESULT(srcptr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    INTEGER, INTENT(IN) :: stype
    TYPE(Sources_TYP), POINTER :: srcptr
    !------------------------------------------------------------------------!
    srcptr => this
    DO
       IF (.NOT.ASSOCIATED(srcptr)) EXIT
       IF (GetType(srcptr).EQ.stype) RETURN
       srcptr => srcptr%next
    END DO
  END FUNCTION GetSourcesPointer
  

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
       END SELECT
       ! deallocate source term structure
       DEALLOCATE(srcptr)
    END DO
    ! release temporary storage
    DEALLOCATE(temp_sterm)
  END SUBROUTINE CloseSources

END MODULE sources_generic
