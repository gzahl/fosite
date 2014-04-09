!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_polar.f90                                                #
!#                                                                           #
!# Copyright (C) 2007 Tobias Illenseer <tillense@astrophysik.uni-kiel.de>    #
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
!> \author Tobias Illenseer
!!
!! \brief define properties of a 2D polar mesh
!!
!! \extends geometry_cartesian
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_polar
  USE geometry_cartesian
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE Convert2Cartesian_polar
     MODULE PROCEDURE Polar2Cartesian_coords, Polar2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear_polar
     MODULE PROCEDURE Cartesian2Polar_coords, Cartesian2Polar_vectors
  END INTERFACE
  !> \endcond
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "polar"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitGeometry_polar, &
       ScaleFactors_polar, &
       Radius_polar, &
       PositionVector_polar, &
       Convert2Cartesian_polar, &
       Convert2Curvilinear_polar, &
       Polar2Cartesian_coords, &
       Polar2Cartesian_vectors, &
       Cartesian2Polar_coords, &
       Cartesian2Polar_vectors
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_polar(this,gt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    !------------------------------------------------------------------------!
    CALL InitGeometry(this,gt,geometry_name)
  END SUBROUTINE InitGeometry_polar
    

  ELEMENTAL SUBROUTINE ScaleFactors_polar(r,hr,hphi,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r
    REAL, INTENT(OUT) :: hr,hphi,hz
    !------------------------------------------------------------------------!
    hr   = 1.
    hphi = r
    hz   = 1.
  END SUBROUTINE ScaleFactors_polar

  ELEMENTAL FUNCTION Radius_polar(r) RESULT(radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r
    REAL :: radius
    !------------------------------------------------------------------------!
    radius = r
  END FUNCTION Radius_polar

  ELEMENTAL SUBROUTINE PositionVector_polar(r,rx,ry)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r
    REAL, INTENT(OUT) :: rx,ry
    !------------------------------------------------------------------------!
    rx = Radius_polar(r)
    ry = 0.0
  END SUBROUTINE PositionVector_polar

  ! coordinate transformations
  ELEMENTAL SUBROUTINE Polar2Cartesian_coords(r,phi,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r,phi
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    x = r*COS(phi)
    y = r*SIN(phi)
  END SUBROUTINE Polar2Cartesian_coords

  ELEMENTAL SUBROUTINE Cartesian2Polar_coords(x,y,r,phi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: x,y
    REAL, INTENT(OUT) :: r,phi
    !------------------------------------------------------------------------!
    r = SQRT(x*x+y*y)
    phi = ATAN2(y,x)
    IF(phi.LT.0.0) THEN
        phi = phi + 2.0*PI
    END IF
  END SUBROUTINE Cartesian2Polar_coords


  ! vector transformations
  ELEMENTAL SUBROUTINE Polar2Cartesian_vectors(phi,vr,vphi,vx,vy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: phi,vr,vphi
    REAL, INTENT(OUT) :: vx,vy
    !------------------------------------------------------------------------!
    vx = vr * COS(phi) - vphi * SIN(phi)
    vy = vr * SIN(phi) + vphi * COS(phi)
  END SUBROUTINE Polar2Cartesian_vectors

  ELEMENTAL SUBROUTINE Cartesian2Polar_vectors(phi,vx,vy,vr,vphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: phi,vx,vy
    REAL, INTENT(OUT) :: vr,vphi
    !------------------------------------------------------------------------!
    vr   = vx * COS(phi) + vy * SIN(phi)
    vphi = -vx * SIN(phi) + vy * COS(phi)
  END SUBROUTINE Cartesian2Polar_vectors
  
END MODULE geometry_polar
