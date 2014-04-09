!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_polypolar.f90                                            #
!#                                                                           #
!# Copyright (C) 2012 Manuel Jung <mjung@astrophysik.uni-kiel.de>            #
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
!! \brief define properties of a 2D polar mesh with a quadratic scaled radius
!!
!! x = gp * r**2 * cos(phi)
!! y = gp * r**2 * sin(phi)
!!
!! \extends geometry_cartesian
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_polypolar
  USE geometry_cartesian
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE Convert2Cartesian_polypolar
     MODULE PROCEDURE Polypolar2Cartesian_coords, Polypolar2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear_polypolar
     MODULE PROCEDURE Cartesian2Polypolar_coords, Cartesian2Polypolar_vectors
  END INTERFACE
  !> \endcond
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "polypolar"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitGeometry_polypolar, &
       ScaleFactors_polypolar, &
       Radius_polypolar, &
       PositionVector_polypolar, &
       Convert2Cartesian_polypolar, &
       Convert2Curvilinear_polypolar, &
       Polypolar2Cartesian_coords, &
       Polypolar2Cartesian_vectors, &
       Cartesian2Polypolar_coords, &
       Cartesian2Polypolar_vectors
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_polypolar(this,gt,gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    REAL, INTENT(IN)    :: gp
    !------------------------------------------------------------------------!
    IF(gp.LE.0.0) THEN
        CALL Error(this, "InitGeometry_polypolar", "geometry parameter 1 has "&
        //"to be greater than 0.0.")
    END IF
    CALL InitGeometry(this,gt,geometry_name)
    CALL SetScale(this,gp)
  END SUBROUTINE InitGeometry_polypolar
    

  ELEMENTAL SUBROUTINE ScaleFactors_polypolar(gp,rho,hr,hphi,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,rho
    REAL, INTENT(OUT) :: hr,hphi,hz
    !------------------------------------------------------------------------!
    hr   = 2.0*gp*rho
    hphi = Radius_polypolar(gp,rho)
    hz   = 1.
  END SUBROUTINE ScaleFactors_polypolar

  ELEMENTAL FUNCTION Radius_polypolar(gp,rho) RESULT(radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,rho
    REAL :: radius
    !------------------------------------------------------------------------!
    radius = gp * rho * rho
  END FUNCTION Radius_polypolar

  ELEMENTAL SUBROUTINE PositionVector_polypolar(gp,rho,rx,ry)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,rho
    REAL, INTENT(OUT) :: rx,ry
    !------------------------------------------------------------------------!
    rx = Radius_polypolar(gp,rho)
    ry = 0.0
  END SUBROUTINE PositionVector_polypolar


  ! coordinate transformations
  ELEMENTAL SUBROUTINE Polypolar2Cartesian_coords(gp,r,phi,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r,phi,gp
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    REAL              :: rs
    !------------------------------------------------------------------------!
    rs = gp*r*r
    x = rs*COS(phi)
    y = rs*SIN(phi)
  END SUBROUTINE Polypolar2Cartesian_coords

  ELEMENTAL SUBROUTINE Cartesian2Polypolar_coords(gp,x,y,r,phi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: x,y,gp
    REAL, INTENT(OUT) :: r,phi
    !------------------------------------------------------------------------!
    r = gp**(-0.5) * (x*x+y*y)**0.25
    phi = ATAN2(y,x)
    IF(phi.LT.0.0) THEN
        phi = phi + 2.0*PI
    END IF
  END SUBROUTINE Cartesian2Polypolar_coords

  ! vector transformations
  ELEMENTAL SUBROUTINE Polypolar2Cartesian_vectors(gp,phi,vr,vphi,vx,vy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: phi,vr,vphi,gp
    REAL, INTENT(OUT) :: vx,vy
    !------------------------------------------------------------------------!
    vx = vr * COS(phi) - vphi * SIN(phi)
    vy = vr * SIN(phi) + vphi * COS(phi)
  END SUBROUTINE Polypolar2Cartesian_vectors

  ELEMENTAL SUBROUTINE Cartesian2Polypolar_vectors(gp,phi,vx,vy,vr,vphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: phi,vx,vy,gp
    REAL, INTENT(OUT) :: vr,vphi
    !------------------------------------------------------------------------!
    vr   = vx * COS(phi) + vy * SIN(phi)
    vphi = -vx * SIN(phi) + vy * COS(phi)
  END SUBROUTINE Cartesian2Polypolar_vectors
  
END MODULE geometry_polypolar
