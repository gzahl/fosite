!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: constants_SI.f90                                                  #
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
! module for SI units and physical constants
!----------------------------------------------------------------------------!
MODULE constants_SI
  USE constants_common
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: units_name = 'SI'
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Constants_TYP, &
       ! methods
       InitConstants_SI, &
       GetType, &
       GetName, &
       GetRank, &
       GetNumProcs, &
       Initialized, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitConstants_SI(this,units)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Constants_TYP) :: this
    INTEGER             :: units
    !------------------------------------------------------------------------!
    INTENT(IN)          :: units
    INTENT(INOUT)       :: this
    !------------------------------------------------------------------------!
    CALL InitConstants(this,units,units_name)
    ! assign numerical values of physical constants in SI units;
    ! (C, GN, etc. are defined in constants_common)
    this%C  = C
    this%GN = GN
    this%KB = KB
    this%NA = NA
    this%SB = SB
    this%KE = KE
    ! conversion factors to SI units are unity
    this%cf_time = 1.0
    this%cf_mass = 1.0
    this%cf_momentum = 1.0
    this%cf_energy = 1.0
    this%cf_power = 1.0
    this%cf_temperature = 1.0
    this%cf_density = 1.0
    this%cf_opacity = 1.0
  END SUBROUTINE InitConstants_SI

END MODULE constants_SI
