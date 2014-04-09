!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_common.f90                                               #
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
! basic boundary module
!----------------------------------------------------------------------------!
MODULE boundary_common
  USE common_types, &
       GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Info_common => Info, Warning_common => Warning, Error_common => Error
  IMPLICIT NONE
#ifdef PARALLEL
    include 'mpif.h'
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE GetType
     MODULE PROCEDURE GetCondition, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetConditionName, GetName_common
  END INTERFACE
  INTERFACE GetRank
     MODULE PROCEDURE GetBoundaryRank, GetRank_common
  END INTERFACE
  INTERFACE GetNumProcs
     MODULE PROCEDURE GetBoundaryNumProcs, GetNumProcs_common
  END INTERFACE
  INTERFACE Info
     MODULE PROCEDURE BoundaryInfo_rank0, BoundaryInfo_rankX, Info_common
  END INTERFACE
  INTERFACE Warning
     MODULE PROCEDURE BoundaryWarning, Warning_common
  END INTERFACE
  INTERFACE Error
     MODULE PROCEDURE BoundaryError_rank0, BoundaryError_rankX, Error_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  TYPE Boundary_TYP
     TYPE(Common_TYP)  :: condition         ! outflow, reflect, periodic, .. !
     TYPE(Common_TYP)  :: direction         ! west, east, south, north       !
     INTEGER           :: IMID, JMID        ! indices of cells in the middle !
     INTEGER           :: nohdim            ! dimension of Noh problem       !
     REAL, DIMENSION(:,:,:), POINTER &
                       :: data              ! boundary data for fixed bc     !
     REAL, DIMENSION(:,:), POINTER &        ! inverse distance to center     !
                       :: invr              !   used for Noh boundary        !
     LOGICAL, DIMENSION(:), POINTER  &      ! masks for reflecting or fixed  !
                       :: reflX,reflY,fixed !   boundaries                   !
#ifdef PARALLEL
     REAL, DIMENSION(:,:,:), POINTER &      ! send and receive buffer for    !
          :: sendbuf,recvbuf                ! boundary data                  !
#endif
  END TYPE Boundary_TYP
  !--------------------------------------------------------------------------!
  INTEGER, PARAMETER :: WEST  = 1
  INTEGER, PARAMETER :: EAST  = 2
  INTEGER, PARAMETER :: SOUTH = 3
  INTEGER, PARAMETER :: NORTH = 4
  CHARACTER(LEN=32), DIMENSION(4), PARAMETER :: direction_name = (/ &
       ' west', ' east', 'south', 'north' /) 
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Boundary_TYP, &
       ! constants
#ifdef PARALLEL
       DEFAULT_MPI_REAL, &
#endif
       WEST, EAST, SOUTH, NORTH, &
       ! methods
       InitBoundary, &
       GetType, &
       GetName, &
       GetDirection, &
       GetDirectionName, &
       GetRank, &
       GetNumProcs, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitBoundary(this,bctype,bcname,direction)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    INTEGER            :: bctype,direction
    CHARACTER(LEN=*)   :: bcname
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
       CALL Error(this,"InitBoundary_common", "Unknown direction")
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


  PURE FUNCTION GetBoundaryRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%condition)
  END FUNCTION GetBoundaryRank


  PURE FUNCTION GetBoundaryNumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this
    INTEGER :: p
    !------------------------------------------------------------------------!
    p = GetNumProcs_common(this%condition)
  END FUNCTION GetBoundaryNumProcs


  SUBROUTINE BoundaryInfo_rank0(this,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: msg
    !------------------------------------------------------------------------!
    CALL Info_common(this%condition,msg)
  END SUBROUTINE BoundaryInfo_rank0


  SUBROUTINE BoundaryInfo_rankX(this,msg,rank)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN)  :: msg
    INTEGER, INTENT(IN)            :: rank
    !------------------------------------------------------------------------!
    CALL Info_common(this%condition,msg,rank)
  END SUBROUTINE BoundaryInfo_rankX


  SUBROUTINE BoundaryWarning(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN)  :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Warning_common(this%condition,modproc,msg)
  END SUBROUTINE BoundaryWarning


  SUBROUTINE BoundaryError_rank0(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN)  :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Error_common(this%condition,modproc,msg)
  END SUBROUTINE BoundaryError_rank0


  SUBROUTINE BoundaryError_rankX(this,modproc,msg,rank)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN)  :: modproc,msg
    INTEGER, INTENT(IN)            :: rank
    !------------------------------------------------------------------------!
    CALL Error_common(this%condition,modproc,msg,rank)
  END SUBROUTINE BoundaryError_rankX

END MODULE boundary_common
