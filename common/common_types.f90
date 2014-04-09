!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: common_types.f90                                                  #
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
! basic data and methods common to all objects
!----------------------------------------------------------------------------!
MODULE common_types
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! common data structure
  TYPE Common_TYP
     PRIVATE
     INTEGER           :: type
     CHARACTER(LEN=32) :: name
  END TYPE Common_TYP
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Common_TYP, &
       ! methods
       InitCommon, &
       GetType, &
       GetName
  !--------------------------------------------------------------------------!

CONTAINS

  PURE SUBROUTINE InitCommon(this,t,n)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Common_TYP)  :: this
    INTEGER           :: t
    CHARACTER(LEN=32) :: n
    !------------------------------------------------------------------------!
    INTENT(IN)        :: t,n
    INTENT(OUT)       :: this
    !------------------------------------------------------------------------!
    this%type = t
    this%name = n
  END SUBROUTINE InitCommon


  PURE FUNCTION GetType(this) RESULT(t)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Common_TYP), INTENT(IN) :: this
    INTEGER :: t
    !------------------------------------------------------------------------!
    t = this%type
  END FUNCTION GetType


  PURE FUNCTION GetName(this) RESULT(n)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Common_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: n
    !------------------------------------------------------------------------!
    n = this%name
  END FUNCTION GetName

END MODULE common_types
