!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_common.f90                                               #
!#                                                                           #
!# Copyright (C) 2006-2008                                                   #
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
! basic geometry module
!----------------------------------------------------------------------------!
MODULE geometry_common
  USE common_types, GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, Info_common => Info, &
       Warning_common => Warning, Error_common => Error
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE GetType
     MODULE PROCEDURE GetCoordsys, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetCoordsysName, GetName_common
  END INTERFACE
  INTERFACE GetRank
     MODULE PROCEDURE GetGeometryRank, GetRank_common
  END INTERFACE
  INTERFACE Info
     MODULE PROCEDURE GeometryInfo, Info_common
  END INTERFACE
  INTERFACE Warning
     MODULE PROCEDURE GeometryWarning, Warning_common
  END INTERFACE
  INTERFACE Error
     MODULE PROCEDURE GeometryError, Error_common
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
       GetName, &
       GetRank, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry(this,cs,cn)
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


  PURE FUNCTION GetGeometryRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%coordsys)
  END FUNCTION GetGeometryRank


  SUBROUTINE GeometryInfo(this,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: msg
    !------------------------------------------------------------------------!
    CALL Info_common(this%coordsys,msg)
  END SUBROUTINE GeometryInfo


  SUBROUTINE GeometryWarning(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Warning_common(this%coordsys,modproc,msg)
  END SUBROUTINE GeometryWarning


  SUBROUTINE GeometryError(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Error_common(this%coordsys,modproc,msg)
  END SUBROUTINE GeometryError


END MODULE geometry_common
