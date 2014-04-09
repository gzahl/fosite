!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_sinhspherical.f90                                        #
!#                                                                           #
!# Copyright (C) 2010                                                        #
!# Björn Sperling <sperling@astrophysik.uni-kiel.de>                         #
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
! define properties of a hyperbolic 2D spherical mesh with dimensionless radial
! coordinate rho (with 0 < rho < +inf) according to:
!    x = r0 * sinh(rho) * sin(theta)  
!    y = r0 * sinh(rho) * cos(theta)
!----------------------------------------------------------------------------!
MODULE geometry_sinhspherical
  USE geometry_cartesian, ONLY : Geometry_TYP, PI, InitGeometry
  USE geometry_spherical
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTERFACE Convert2Cartesian_sinhspher
     MODULE PROCEDURE Sinhspherical2Cartesian_coords, Spherical2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear_sinhspher
     MODULE PROCEDURE Cartesian2Sinhspherical_coords, Cartesian2Spherical_vectors
  END INTERFACE
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "sinhspherical"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitGeometry_sinhspher, &
       ScaleFactors_sinhspher, &
       Convert2Cartesian_sinhspher, &
       Convert2Curvilinear_sinhspher, &
       Sinhspherical2Cartesian_coords, &
       Cartesian2Sinhspherical_coords
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_sinhspher(this,gt,gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    REAL, INTENT(IN) :: gp
    !------------------------------------------------------------------------!
    CALL InitGeometry(this,gt,geometry_name)
    this%geoparam = gp
  END SUBROUTINE InitGeometry_sinhspher
    

  ELEMENTAL SUBROUTINE ScaleFactors_sinhspher(gp,rho,theta,hrho,htheta,hphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,rho,theta
    REAL, INTENT(OUT) :: hrho,htheta,hphi
    !------------------------------------------------------------------------!
    hrho   = gp*COSH(rho)
    htheta = gp*SINH(rho)
    hphi   = htheta*SIN(theta)
  END SUBROUTINE ScaleFactors_sinhspher


  ! coordinate transformations
  ELEMENTAL SUBROUTINE Sinhspherical2Cartesian_coords(gp,rho,theta,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,rho,theta
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    REAL :: r
    !------------------------------------------------------------------------!
    r = gp*SINH(rho)
    x = r*SIN(theta)
    y = r*COS(theta)
  END SUBROUTINE Sinhspherical2Cartesian_coords

  ELEMENTAL SUBROUTINE Cartesian2Sinhspherical_coords(gp,x,y,rho,theta)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,x,y
    REAL, INTENT(OUT) :: rho,theta
    !------------------------------------------------------------------------!
    REAL :: x1,y1,r1
    !------------------------------------------------------------------------!
    x1 = x/gp
    y1 = y/gp
    r1 = SQRT(x1*x1+y1*y1)
    rho = LOG(r1+SQRT(1.+r1*r1))  ! = ASINH(r1))
    theta=ACOS(y1/r1)
  END SUBROUTINE Cartesian2Sinhspherical_coords

END MODULE geometry_sinhspherical
