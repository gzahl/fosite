!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: constants_generic.f90                                             #
!#                                                                           #
!# Copyright (C) 2007-2008                                                   #
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
! generic module for units and physical constants
!----------------------------------------------------------------------------!
MODULE constants_generic
  USE constants_SI
  USE constants_geometrical
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: SI          = 1
  INTEGER, PARAMETER :: CGS         = 2
  INTEGER, PARAMETER :: GEOMETRICAL = 3
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Constants_TYP, &
       ! constant flags
       SI, CGS, GEOMETRICAL, &
       ! methods
       InitConstants
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitConstants(this,units)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Constants_TYP) :: this
    INTEGER             :: units
    !------------------------------------------------------------------------!
    INTENT(IN)          :: units
    INTENT(INOUT)       :: this
    !------------------------------------------------------------------------!
    SELECT CASE(units)
    CASE(SI)
       CALL InitConstants_SI(this,units)
    CASE(GEOMETRICAL)
       CALL InitConstants_geometrical(this,units)
    CASE DEFAULT
       CALL Error(this, "InitConstants", "Unknown physical units.")
    END SELECT

    ! derived constants
    this%RG = this%KB * this%NA         ![J/mol/K]    universal gas constant !

    ! print some information
    CALL Info(this, " CONSTANTS> physical units:    " // TRIM(GetName(this)))
  END SUBROUTINE InitConstants

END MODULE constants_generic
