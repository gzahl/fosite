!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_common.f90                                                #
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
! basic physics module
!----------------------------------------------------------------------------!
MODULE physics_common
  USE common_types, &
       GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, Info_common => Info, &
       Warning_common => Warning, Error_common => Error
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
  INTERFACE GetRank
     MODULE PROCEDURE GetPhysicsRank, GetRank_common
  END INTERFACE
  INTERFACE Info
     MODULE PROCEDURE PhysicsInfo, Info_common
  END INTERFACE
  INTERFACE Warning
     MODULE PROCEDURE PhysicsWarning, Warning_common
  END INTERFACE
  INTERFACE Error
     MODULE PROCEDURE PhysicsError, Error_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  TYPE Physics_TYP
     TYPE(Common_TYP)       :: advproblem            ! advection problem     !
     TYPE(Constants_TYP)    :: constants             ! physical constants    !
     TYPE(Sources_TYP), POINTER &
                            :: sources               ! list of source terms  !
     REAL                   :: gamma                 ! ratio of spec. heats  !
     REAL                   :: mu                    ! mean molecular weight !
     REAL                   :: rhomin                ! density minimum       !
     REAL                   :: pmin                  ! pressure minimum      !
     REAL                   :: dpmax                 ! for time step control !
     INTEGER                :: VNUM                  ! number of variables   !
     INTEGER                :: DENSITY               ! array indices for     !
     INTEGER                :: PRESSURE, ENERGY      !    primitive and      !
     INTEGER                :: XVELOCITY, XMOMENTUM  !    conservative       !
     INTEGER                :: YVELOCITY, YMOMENTUM  !    variables          !
     INTEGER                :: ZVELOCITY, ZMOMENTUM  !                       !
     CHARACTER(LEN=16), DIMENSION(:), POINTER &
                            :: pvarname,cvarname     ! names of variables    !
     REAL, DIMENSION(:,:,:), POINTER &
                            :: csound                ! sound speed           !
     REAL, DIMENSION(:,:), POINTER &
                            :: amin, amax, bmin, bmax! wave speeds           !
     REAL, DIMENSION(:,:,:), POINTER &
                            :: tmin, tmax            ! temporary storage     !
     REAL, DIMENSION(:,:,:,:), POINTER &
                            :: fcent                 ! centrifugal force     !
  END TYPE Physics_TYP
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Physics_TYP, &
       ! methods
       InitPhysics, &
       ClosePhysics, &
       GetType, &
       GetName, &
       GetRank, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitPhysics(this,atype,aname,vnum)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    INTEGER           :: atype
    CHARACTER(LEN=32) :: aname
    INTEGER           :: vnum
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: atype,aname,vnum
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%advproblem,atype,aname)
    this%vnum=vnum
    ALLOCATE(this%pvarname(this%vnum),this%cvarname(this%vnum),STAT=err)
    IF (err.NE.0) CALL Error(this,"InitPhysics", &
         "unable to allocate memory")
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


  PURE FUNCTION GetPhysicsRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP), INTENT(IN) :: this
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%advproblem)
  END FUNCTION GetPhysicsRank


  SUBROUTINE PhysicsInfo(this,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: msg
    !------------------------------------------------------------------------!
    CALL Info_common(this%advproblem,msg)
  END SUBROUTINE PhysicsInfo


  SUBROUTINE PhysicsWarning(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Warning_common(this%advproblem,modproc,msg)
  END SUBROUTINE PhysicsWarning


  SUBROUTINE PhysicsError(this,modproc,msg,rank)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    INTEGER, OPTIONAL, INTENT(IN) :: rank
    !------------------------------------------------------------------------!
    IF (PRESENT(rank)) THEN
       CALL Error_common(this%advproblem,modproc,msg,rank)
    ELSE
       CALL Error_common(this%advproblem,modproc,msg)
    END IF
  END SUBROUTINE PhysicsError


  SUBROUTINE ClosePhysics(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%pvarname,this%cvarname)
  END SUBROUTINE ClosePhysics


END MODULE physics_common
