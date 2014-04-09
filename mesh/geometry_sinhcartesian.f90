!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_sinhcartesian.f90                                        #
!#                                                                           #
!# Copyright (C) 2013                                                        #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
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
!> \author Manuel Jung
!!
!! \brief define properties of a sinh cartesian mesh
!!
!! dimensionless coordinates xi and eta are as following:
!!    x = gp1 * sinh(xi)
!!    y = gp1 * sinh(eta)
!!
!! \extends geometry_cartesian
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_sinhcartesian
  USE geometry_cartesian
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE Convert2Cartesian_sinhcart
     MODULE PROCEDURE Sinhcartesian2Cartesian_coords
  END INTERFACE
  INTERFACE Convert2Curvilinear_sinhcart
     MODULE PROCEDURE Cartesian2Sinhcartesian_coords
  END INTERFACE
  INTERFACE PositionVector_sinhcartesian
     MODULE PROCEDURE Sinhcartesian2Cartesian_coords
  END INTERFACE
  !> \endcond
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "sinhcartesian"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitGeometry_sinhcartesian, &
       ScaleFactors_sinhcartesian, &
       Radius_sinhcartesian, &
       PositionVector_sinhcartesian, &
       Convert2Cartesian_sinhcart, &
       Convert2Curvilinear_sinhcart, &
       Sinhcartesian2Cartesian_coords, &
       Cartesian2Sinhcartesian_coords
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_sinhcartesian(this,gt,gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    REAL, INTENT(IN) :: gp
    !------------------------------------------------------------------------!
    CALL InitGeometry(this,gt,geometry_name)
    CALL SetScale(this,gp)
  END SUBROUTINE InitGeometry_sinhcartesian
    

  ELEMENTAL SUBROUTINE ScaleFactors_sinhcartesian(gp,xi,eta,hxi,heta,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,xi,eta
    REAL, INTENT(OUT) :: hxi,heta,hz
    !------------------------------------------------------------------------!
    hxi = gp*COSH(xi)
    heta = gp*COSH(eta)
    hz   = 1.
  END SUBROUTINE ScaleFactors_sinhcartesian

  ELEMENTAL FUNCTION Radius_sinhcartesian(gp,xi,eta) RESULT(radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,xi,eta
    REAL :: radius
    !------------------------------------------------------------------------!
    radius = gp*SQRT(SINH(xi)**2+SINH(eta)**2)
  END FUNCTION Radius_sinhcartesian

  ! coordinate transformations
  ELEMENTAL SUBROUTINE Sinhcartesian2Cartesian_coords(gp,xi,eta,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,xi,eta
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    x = gp*SINH(xi)
    y = gp*SINH(eta)
  END SUBROUTINE Sinhcartesian2Cartesian_coords

  ELEMENTAL SUBROUTINE Cartesian2Sinhcartesian_coords(gp,x,y,xi,eta)
    USE functions, ONLY : Asinh
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,x,y
    REAL, INTENT(OUT) :: xi,eta
    !------------------------------------------------------------------------!
    xi = Asinh(x/gp)
    eta = Asinh(y/gp)
  END SUBROUTINE Cartesian2Sinhcartesian_coords
    
END MODULE geometry_sinhcartesian
