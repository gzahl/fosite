!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_common.f90                                               #
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
! basic geometry module
!----------------------------------------------------------------------------!
MODULE geometry_common
  USE common_types, GetType_common => GetType, GetName_common => GetName
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE GetType
     MODULE PROCEDURE GetCoordsys, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetCoordsysName, GetName_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  TYPE Geometry_TYP
     TYPE(Common_TYP) :: coordsys                   ! cartesian, polar, etc. !
     REAL             :: scalefactor !scale factor, e.g. for oblate-spheroi. !
  END TYPE Geometry_TYP
  SAVE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Geometry_TYP, &
       ! methods
       InitGeometry, &
       GetType, &
       GetName
  !--------------------------------------------------------------------------!

CONTAINS

  PURE SUBROUTINE InitGeometry(this,cs,cn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP) :: this
    INTEGER            :: cs
    CHARACTER(LEN=32)  :: cn
    !------------------------------------------------------------------------!
    INTENT(IN)         :: cs,cn
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%coordsys,cs,cn)
  END SUBROUTINE InitGeometry


  PURE FUNCTION GetCoordsys(this) RESULT(cs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    INTEGER :: cs
    !------------------------------------------------------------------------!
    cs = GetType_common(this%coordsys)
  END FUNCTION GetCoordsys


  PURE FUNCTION GetCoordsysName(this) RESULT(cn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: cn    
    !------------------------------------------------------------------------!
    cn = GetName_common(this%coordsys)
  END FUNCTION GetCoordsysName


END MODULE geometry_common
