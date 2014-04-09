!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_common.f90                                                #
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
! basic physics module
!----------------------------------------------------------------------------!
MODULE physics_common
  USE common_types, GetType_common => GetType, GetName_common => GetName
  USE sources_common, ONLY : Sources_TYP
  USE constants_common, ONLY : Constants_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE GetType
     MODULE PROCEDURE GetAdvProblem, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetAdvProblemName, GetName_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  TYPE Physics_TYP
     TYPE(Common_TYP)       :: advproblem            ! advection problem     !
     TYPE(Sources_TYP), POINTER &
                            :: sources               ! list of source terms  !
     TYPE(Constants_TYP)    :: constants             ! physical constants    !
     REAL, DIMENSION(:,:,:), POINTER &
                            :: csound                ! sound speed           !
     REAL, DIMENSION(:,:), POINTER &
                            :: amin, amax, bmin, bmax! wave speeds           !
     REAL, DIMENSION(:,:,:), POINTER &
                            :: tmin, tmax            ! temporary storage     !
     REAL, DIMENSION(:,:,:,:), POINTER &
                            :: fcent                 ! centrifugal force     !
     REAL                   :: gamma                 ! ratio of spec. heats  !
     REAL                   :: mu                    ! mean molecular weight !
     REAL                   :: dpmax                 ! for time step control !
     INTEGER                :: vnum                  ! number of variables   !
  END TYPE Physics_TYP
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Physics_TYP, &
       ! methods
       InitPhysics, &
       GetType, &
       GetName
  !--------------------------------------------------------------------------!

CONTAINS

  PURE SUBROUTINE InitPhysics(this,atype,aname,vnum)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    INTEGER           :: atype
    CHARACTER(LEN=32) :: aname
    INTEGER           :: vnum
    !------------------------------------------------------------------------!
    INTENT(IN)        :: atype,aname,vnum
    INTENT(OUT)       :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%advproblem,atype,aname)
    this%vnum=vnum
  END SUBROUTINE InitPhysics


  PURE FUNCTION GetAdvProblem(this) RESULT(ap)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP), INTENT(IN) :: this
    INTEGER :: ap
    !------------------------------------------------------------------------!
    ap = GetType_common(this%advproblem)
  END FUNCTION GetAdvProblem


  PURE FUNCTION GetAdvProblemName(this) RESULT(an)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: an
    !------------------------------------------------------------------------!
    an = GetName_common(this%advproblem)
  END FUNCTION GetAdvProblemName


END MODULE physics_common
