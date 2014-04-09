!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_common.f90                                               #
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
! basic boundary module
!----------------------------------------------------------------------------!
MODULE boundary_common
  USE common_types, GetType_common => GetType, GetName_common => GetName
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE GetType
     MODULE PROCEDURE GetCondition, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetConditionName, GetName_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  TYPE Boundary_TYP
     TYPE(Common_TYP)  :: condition         ! outflow, reflect, periodic, .. !
     TYPE(Common_TYP)  :: direction         ! west, east, south, north       !
     REAL, DIMENSION(:,:,:), POINTER &
                       :: data              ! boundary data for fixed bc     !
     LOGICAL, DIMENSION(:), POINTER  &      ! masks for reflecting or fixed  !
                       :: reflX,reflY,fixed !   boundaries                   !
     INTEGER           :: IMID, JMID        ! indices of cells in the middle !
  END TYPE Boundary_TYP
  !--------------------------------------------------------------------------!
  INTEGER, PARAMETER :: WEST  = 1
  INTEGER, PARAMETER :: EAST  = 2
  INTEGER, PARAMETER :: SOUTH = 3
  INTEGER, PARAMETER :: NORTH = 4
  CHARACTER(LEN=32), DIMENSION(4), PARAMETER :: direction_name = (/ &
       'west ', 'east ', 'south', 'north' /) 
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Boundary_TYP, &
       ! constants
       WEST, EAST, SOUTH, NORTH, &
       ! methods
       InitBoundary, &
       GetType, &
       GetName, &
       GetDirection, &
       GetDirectionName
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitBoundary(this,bctype,bcname,direction)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    INTEGER            :: bctype,direction
    CHARACTER(LEN=32)  :: bcname
    !------------------------------------------------------------------------!
    INTENT(IN)    :: bctype,bcname,direction
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    ! set boundary condition
    CALL InitCommon(this%condition,bctype,bcname)
    ! check for wrong direction
    SELECT CASE(direction)
    CASE(WEST,EAST,NORTH,SOUTH)
       ! ok
    CASE DEFAULT
       PRINT *, "ERROR in InitBoundary: unknown direction"
       STOP
    END SELECT
    ! set direction
    CALL InitCommon(this%direction,direction,direction_name(direction))
  END SUBROUTINE InitBoundary


  PURE FUNCTION GetCondition(this) RESULT(bt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this
    INTEGER :: bt
    !------------------------------------------------------------------------!
    bt=GetType_common(this%condition)
  END FUNCTION GetCondition


  PURE FUNCTION GetConditionName(this) RESULT(bn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: bn
    !------------------------------------------------------------------------!
    bn=GetName_common(this%condition)
  END FUNCTION GetConditionName


  PURE FUNCTION GetDirection(this) RESULT(dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this
    INTEGER :: dir
    !------------------------------------------------------------------------!
    dir=GetType_common(this%direction)
  END FUNCTION GetDirection


  PURE FUNCTION GetDirectionName(this) RESULT(dn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: dn
    !------------------------------------------------------------------------!
    dn=GetName_common(this%direction)
  END FUNCTION GetDirectionName

END MODULE boundary_common
