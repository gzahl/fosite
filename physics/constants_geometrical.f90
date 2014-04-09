!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: constants_geometrical.f90                                         #
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
! module for geometrical units and physical constants
!----------------------------------------------------------------------------!
MODULE constants_geometrical
  USE constants_common
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: units_name = 'geometrical'
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Constants_TYP, &
       ! methods
       InitConstants_geometrical
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitConstants_geometrical(this,units)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Constants_TYP) :: this
    INTEGER             :: units
    !------------------------------------------------------------------------!
    REAL                :: C,GN,KB
    !------------------------------------------------------------------------!
    INTENT(IN)          :: units
    INTENT(INOUT)       :: this
    !------------------------------------------------------------------------!
    CALL InitConstants(this,units,units_name)
    ! numerical values of physical constants in geometrical units (c = G = 1)
    this%C  = 1.0                                   ! lightspeed             !
    this%GN = 1.0                                   ! Newtons grav. constant !
    this%KB = 1.0                                   ! Boltzmann constant     !
    this%KE = 3.48E-02                  ![m^2/kg]     electron scat. opacity !
    this%NA = 6.022E+23                 ![1/mol]      Avogadro constant      !
    ! basic constants in SI units
    C  = 2.99792458E+08                 ![m/s]        lightspeed             !
    GN = 6.6742E-11                     ![m^3/kg/s^2] Newtons grav. constant !
    KB = 1.3806505E-23                  ![J/K]        Boltzmann constant     !
    ! factors for conversion to SI units
    ! time and mass have unit of length i.e. metre
    this%cf_time = C                 ! i.e. 1 sec [SI] = 3e8 m [geometrical] !
    this%cf_mass = (GN/C)/C
    this%cf_momentum = this%cf_mass/C
    this%cf_energy = this%cf_momentum/C
    this%cf_power = this%cf_energy/C
    this%cf_temperature = KB * this%cf_energy
    this%cf_density = this%cf_mass
  END SUBROUTINE InitConstants_geometrical

END MODULE constants_geometrical
