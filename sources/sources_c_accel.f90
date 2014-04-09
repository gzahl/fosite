!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_c_accel.f90                                               #
!#                                                                           #
!# Copyright (C) 2008                                                        #
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
! source terms module for constant acceleration 
!----------------------------------------------------------------------------!
MODULE sources_c_accel
  USE sources_common, InitSources_common => InitSources
  USE fluxes_common, ONLY : Fluxes_TYP
  USE physics_generic
  USE mesh_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  REAL, PARAMETER :: TINY = 1.0E-30              ! to avoid division by 0    !
  CHARACTER(LEN=32), PARAMETER :: source_name = "constant acceleration"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! methods
       InitSources2, &
       InitSources_c_accel, &
       GetType, &
       GetName, &
       ExternalSources_c_accel, &
       CloseSources_c_accel
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources2(this,stype,sname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    INTEGER           :: stype
    CHARACTER(LEN=32) :: sname
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: newsrc, tmpsrc
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: stype,sname
    !------------------------------------------------------------------------!
    ! allocate memory for new source term
    ALLOCATE(newsrc,STAT=err)
    IF (err.NE.0) THEN
       PRINT *, "ERROR in InitSources: Unable to allocate memory"
       STOP
    END IF
    
    ! basic initialization
    CALL InitSources_common(newsrc,stype,sname)

    ! add new source term to beginning of
    ! list of source terms
    IF (ASSOCIATED(this).EQV..FALSE.) THEN
       this => newsrc
    ELSE
       tmpsrc => this
       this => newsrc
       this%next => tmpsrc
    END IF

  END SUBROUTINE InitSources2


  SUBROUTINE InitSources_c_accel(this,Mesh,Physics,stype,xaccel,yaccel)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Physics_TYP) :: Physics
    INTEGER           :: stype
    REAL              :: xaccel, yaccel
    !------------------------------------------------------------------------!
    INTEGER           :: err
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,stype
    !------------------------------------------------------------------------!
    CALL InitSources2(this,stype,source_name)

  
    ALLOCATE(this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
         STAT = err)
    IF (err.NE.0) THEN
       PRINT *, "ERROR in InitSources_c_accel: Unable allocate memory!"
       STOP
    END IF
    ! initialize gravitational acceleration
    this%accel(:,:,1) = xaccel
    this%accel(:,:,2) = yaccel
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
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,cvar,pvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! gravitational source terms
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
  END SUBROUTINE CloseSources_c_accel

END MODULE sources_c_accel
