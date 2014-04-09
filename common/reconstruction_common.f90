!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: reconstruction_common.f90                                         #
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
! basic reconstruction module
!----------------------------------------------------------------------------!
MODULE reconstruction_common
  USE common_types, GetType_common => GetType, GetName_common => GetName
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE GetType
     MODULE PROCEDURE GetOrder, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetOrderName, GetName_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  ! Reconstruction data structure
  TYPE Reconstruction_TYP
     TYPE(Common_TYP)  :: order                    ! constant, linear, ...   !
     TYPE(Common_TYP)  :: limiter                  ! limiter function        !
     LOGICAL           :: primcons                 ! true if primitive       !
     REAL              :: limiter_param            ! limiter parameter       !
     REAL, DIMENSION(:,:,:), &
          POINTER               :: xslopes,yslopes ! limited slopes          !
     REAL, DIMENSION(:,:,:), &
          POINTER  :: temp1,temp2,temp3,temp4 ! temporary storage for slopes !
  END TYPE Reconstruction_TYP
  !--------------------------------------------------------------------------!
  LOGICAL, PARAMETER :: PRIMITIVE    = .TRUE.
  LOGICAL, PARAMETER :: CONSERVATIVE = .FALSE.
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Reconstruction_TYP, &
       ! constants
       PRIMITIVE, CONSERVATIVE, &
       ! methods
       InitReconstruction, &
       GetType, &
       GetName, &
       PrimRecon
  !--------------------------------------------------------------------------!

CONTAINS

  PURE SUBROUTINE InitReconstruction(this,rtype,rname,pc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP) :: this
    INTEGER                  :: rtype
    CHARACTER(LEN=32)        :: rname
    LOGICAL                  :: pc
    !------------------------------------------------------------------------!
    INTENT(IN)               :: rtype,rname,pc
    INTENT(INOUT)            :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%order,rtype,rname)
    this%primcons  = pc
  END SUBROUTINE InitReconstruction


  PURE FUNCTION GetOrder(this) RESULT(rt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP), INTENT(IN) :: this
    INTEGER :: rt
    !------------------------------------------------------------------------!
    rt = GetType_common(this%order)
  END FUNCTION GetOrder


  PURE FUNCTION GetOrderName(this) RESULT(rn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: rn
    !------------------------------------------------------------------------!
    rn = GetName(this%order)
  END FUNCTION GetOrderName


  PURE FUNCTION PrimRecon(this) RESULT(pc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP), INTENT(IN) :: this
    LOGICAL :: pc
    !------------------------------------------------------------------------!
    pc = this%primcons
  END FUNCTION PrimRecon
  
END MODULE reconstruction_common
