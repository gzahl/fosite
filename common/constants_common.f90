!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: constants_common.f90                                              #
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
! basic physical constants module
!----------------------------------------------------------------------------!
MODULE constants_common
  USE common_types, GetType_common => GetType, GetName_common => GetName
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE GetType
     MODULE PROCEDURE GetUnits, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetUnitsName, GetName_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  TYPE Constants_TYP
     TYPE(Common_TYP) :: units                      ! cartesian, polar, etc. !
     ! some physical constants
     DOUBLE PRECISION :: C                          ! light speed            !
     DOUBLE PRECISION :: GN                         ! Newtons grav. constant !
     DOUBLE PRECISION :: KB                         ! Boltzmann constant     !
     DOUBLE PRECISION :: KE                         ! electr. scat. opacity  !
     DOUBLE PRECISION :: NA                         ! Avogadro constant      !
     DOUBLE PRECISION :: RG                         ! gas constant           !
     ! factors for convertion from SI to other units
     DOUBLE PRECISION :: cf_time                    ! time scale             !
     DOUBLE PRECISION :: cf_mass                    ! mass scale             !
     DOUBLE PRECISION :: cf_momentum                ! momentum scale         !
     DOUBLE PRECISION :: cf_energy                  ! energy scale           !
     DOUBLE PRECISION :: cf_power                   ! power scale            !
     DOUBLE PRECISION :: cf_temperature             ! temperature scale      !
     DOUBLE PRECISION :: cf_density                 ! density scale          !
  END TYPE Constants_TYP
  SAVE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Constants_TYP, &
       ! methods
       InitConstants, &
       GetType, &
       GetName
  !--------------------------------------------------------------------------!

CONTAINS

  PURE SUBROUTINE InitConstants(this,ut,un)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Constants_TYP) :: this
    INTEGER             :: ut
    CHARACTER(LEN=32)   :: un
    !------------------------------------------------------------------------!
    INTENT(IN)          :: ut,un
    INTENT(INOUT)       :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%units,ut,un)
  END SUBROUTINE InitConstants


  PURE FUNCTION GetUnits(this) RESULT(ut)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Constants_TYP), INTENT(IN) :: this
    INTEGER :: ut
    !------------------------------------------------------------------------!
    ut = GetType_common(this%units)
  END FUNCTION GetUnits


  PURE FUNCTION GetUnitsName(this) RESULT(un)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Constants_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: un
    !------------------------------------------------------------------------!
    un = GetName_common(this%units)
  END FUNCTION GetUnitsName


END MODULE constants_common
