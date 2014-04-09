!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_c_accel.f90                                               #
!#                                                                           #
!# Copyright (C) 2009,2011                                                   #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
!> \author Björn Sperling
!! \author Tobias Illenseer
!!
!! \brief source terms module for constant acceleration
!!
!! \extends sources_common
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_c_accel
  USE sources_common
  USE fluxes_common, ONLY : Fluxes_TYP
  USE physics_generic
  USE mesh_generic
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: source_name = "constant acceleration"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! methods
       InitSources, &
       CloseSources, &
       InitSources_c_accel, &
       ExternalSources_c_accel, &
       CloseSources_c_accel, &
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


  SUBROUTINE InitSources_c_accel(this,Mesh,Physics,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)          :: Mesh
    TYPE(Fluxes_TYP)        :: Fluxes
    TYPE(Physics_TYP)       :: Physics
    TYPE(Dict_TYP),POINTER  :: config
    INTEGER                 :: stype
    REAL                    :: xaccel, yaccel, zaccel
    !------------------------------------------------------------------------!
    INTEGER                 :: err
    !------------------------------------------------------------------------!
    INTENT(IN)              :: Mesh,Physics
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "stype", stype)
    CALL InitSources(this,stype,source_name)
  
    ALLOCATE(this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM), &
         STAT = err)
    IF (err.NE.0) &
         CALL Error(this,"InitSources_c_accel","memory allocation failed")

    ! initialize constant acceleration
    CALL RequireKey(config, "xaccel", 0.0)
    CALL GetAttr(config, "xaccel", xaccel)
    this%accel(:,:,1) = xaccel

    CALL RequireKey(config, "yaccel", 0.0)
    CALL GetAttr(config, "yaccel", yaccel)
    this%accel(:,:,2) = yaccel

    IF (Physics%DIM .GE. 3) THEN
      CALL RequireKey(config, "zaccel", 0.0)
      CALL GetAttr(config, "zaccel", zaccel)
      this%accel(:,:,3) = zaccel
    END IF
  END SUBROUTINE InitSources_c_accel


  PURE SUBROUTINE ExternalSources_c_accel(this,Mesh,Physics,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: cvar,pvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,cvar,pvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! compute source terms due to constant acceleration
    CALL ExternalSources(Physics,Mesh,this%accel,pvar,cvar,sterm)
  END SUBROUTINE ExternalSources_c_accel

 
  SUBROUTINE CloseSources_c_accel(this,Fluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Fluxes_TYP)  :: Fluxes
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Fluxes
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%accel)
    CALL CloseSources(this)
  END SUBROUTINE CloseSources_c_accel

END MODULE sources_c_accel
