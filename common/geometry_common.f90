!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_common.f90                                               #
!#                                                                           #
!# Copyright (C) 2006-2010                                                   #
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
!> \defgroup geometry geometry
!! \{
!! \brief Family of geometry modules
!! \}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!!
!! \brief basic geometry module
!!
!! \extends common_types
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_common
  USE common_types, &
       GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Initialized_common => Initialized, Info_common => Info, &
       Warning_common => Warning, Error_common => Error
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! global constants
  REAL, PARAMETER :: PI = 3.1415926535897932384626433832795028842
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE GetType
     MODULE PROCEDURE GetCoordsys, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetCoordsysName, GetName_common
  END INTERFACE
  INTERFACE GetRank
     MODULE PROCEDURE GetGeometryRank, GetRank_common
  END INTERFACE
  INTERFACE GetNumProcs
     MODULE PROCEDURE GetGeometryNumProcs, GetNumProcs_common
  END INTERFACE
  INTERFACE GetScale
     MODULE PROCEDURE GetScale1, GetScale2
  END INTERFACE
  INTERFACE SetScale
     MODULE PROCEDURE SetScale1, SetScale2, SetScale3
  END INTERFACE
  INTERFACE Initialized
     MODULE PROCEDURE GeometryInitialized, Initialized_common
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
  !> \endcond
  !--------------------------------------------------------------------------!
  TYPE Geometry_TYP
     TYPE(Common_TYP)  :: coordsys                  !< cartesian, polar, etc.
     REAL,DIMENSION(3) :: geoparam                  !< geometry parameter
  END TYPE Geometry_TYP
  SAVE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Geometry_TYP, &
       ! constants
       PI, &
       ! methods
       InitGeometry, &
       CloseGeometry, &
       GetType, &
       GetName, &
       GetRank, &
       GetNumProcs, &
       GetScale, &
       SetScale, &
       Initialized, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public
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
    ! set geometry parameter to default value
    CALL SetScale(this,1.0,1.0,1.0)
  END SUBROUTINE InitGeometry


  !> \public
  SUBROUTINE CloseGeometry(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL CloseCommon(this%coordsys)
  END SUBROUTINE CloseGeometry


  !> \public
  PURE FUNCTION GetCoordsys(this) RESULT(cs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    INTEGER :: cs
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    cs = GetType_common(this%coordsys)
  END FUNCTION GetCoordsys


  !> \public
  PURE FUNCTION GetCoordsysName(this) RESULT(cn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: cn    
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    cn = GetName_common(this%coordsys)
  END FUNCTION GetCoordsysName


  !> \public
  PURE FUNCTION GetGeometryRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%coordsys)
  END FUNCTION GetGeometryRank


  !> \public
  PURE FUNCTION GetGeometryNumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    INTEGER :: p
    !------------------------------------------------------------------------!
    p = GetNumProcs_common(this%coordsys)
  END FUNCTION GetGeometryNumProcs


  PURE FUNCTION GetScale1(this) RESULT(gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    REAL :: gp
    !------------------------------------------------------------------------!
    gp = this%geoparam(1)
  END FUNCTION GetScale1

  PURE FUNCTION GetScale2(this,i) RESULT(gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    INTEGER, INTENT(IN)            :: i
    REAL :: gp
    !------------------------------------------------------------------------!
    gp = this%geoparam(i)
  END FUNCTION GetScale2

  PURE SUBROUTINE SetScale1(this,gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    REAL, INTENT(IN) :: gp
    !------------------------------------------------------------------------!
    this%geoparam(1) = gp
  END SUBROUTINE SetScale1

  PURE SUBROUTINE SetScale2(this,gp,gp2)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    REAL, INTENT(IN) :: gp,gp2
    !------------------------------------------------------------------------!
    this%geoparam(1) = gp
    this%geoparam(2) = gp2
  END SUBROUTINE SetScale2

  PURE SUBROUTINE SetScale3(this,gp,gp2,gp3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    REAL, INTENT(IN) :: gp,gp2,gp3
    !------------------------------------------------------------------------!
    this%geoparam(1) = gp
    this%geoparam(2) = gp2
    this%geoparam(3) = gp3
  END SUBROUTINE SetScale3

  !> \public
  PURE FUNCTION GeometryInitialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    LOGICAL :: i
    !------------------------------------------------------------------------!
    i = Initialized_common(this%coordsys)
  END FUNCTION GeometryInitialized


  !> \public
  SUBROUTINE GeometryInfo(this,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: msg
    !------------------------------------------------------------------------!
    CALL Info_common(this%coordsys,msg)
  END SUBROUTINE GeometryInfo


  !> \public
  SUBROUTINE GeometryWarning(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Warning_common(this%coordsys,modproc,msg)
  END SUBROUTINE GeometryWarning


  !> \public
  SUBROUTINE GeometryError(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Error_common(this%coordsys,modproc,msg)
  END SUBROUTINE GeometryError


END MODULE geometry_common
