!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_common.f90                                                #
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
! basic sources module
!----------------------------------------------------------------------------!
MODULE sources_common
  USE common_types, GetType_common => GetType, GetName_common => GetName
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE GetType
     MODULE PROCEDURE GetSourceType, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetSourceTypeName, GetName_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  TYPE Sources_TYP
     TYPE(Common_TYP)                :: sourcetype   ! type of source term   !
     TYPE(Sources_TYP), POINTER      :: next => null() ! next source in list !
     REAL                            :: dtmin        ! minimal time step     !
     REAL                            :: mass         ! mass of point source  !
     REAL                            :: mdot         ! disk accretion rate   !
     REAL, DIMENSION(:,:,:), POINTER :: accel        ! acceleration          
     REAL, DIMENSION(:,:,:), POINTER :: temp1,temp2  ! temp. storage         !
  END TYPE SOURCES_TYP
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! methods
       InitSources, &
       GetType, &
       GetName
  !--------------------------------------------------------------------------!

CONTAINS

  PURE SUBROUTINE InitSources(this,stype,sname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    INTEGER           :: stype
    CHARACTER(LEN=32) :: sname
    !------------------------------------------------------------------------!
    INTENT(IN)        :: stype,sname
    INTENT(OUT)       :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%sourcetype,stype,sname)
    NULLIFY(this%next)
  END SUBROUTINE InitSources


  PURE FUNCTION GetSourceType(this) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    INTEGER :: st
    !------------------------------------------------------------------------!
    st = GetType_common(this%sourcetype)
  END FUNCTION GetSourceType


  PURE FUNCTION GetSourceTypeName(this) RESULT(sn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: sn
    !------------------------------------------------------------------------!
    sn = GetName_common(this%sourcetype)
  END FUNCTION GetSourceTypeName


END MODULE sources_common
