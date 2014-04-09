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
    ! numerical values of physical constants in SI units
    this%C  = 2.99792458E+08            ![m/s]        lightspeed             !
    this%GN = 6.6742E-11                ![m^3/kg/s^2] Newtons grav. constant !
    this%KB = 1.3806505E-23             ![J/K]        Boltzmann constant     !
    this%KE = 3.48E-02                  ![m^2/kg]     electron scat. opacity !
    this%NA = 6.022E+23                 ![1/mol]      Avogadro constant      !
  END SUBROUTINE InitConstants_SI

END MODULE constants_SI
