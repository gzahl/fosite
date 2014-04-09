!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_sinhpolar.f90                                            #
!#                                                                           #
!# Copyright (C) 2008 Tobias Illenseer <tillense@astrophysik.uni-kiel.de>    #
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
!! \brief define properties of a hyperbolic 2D polar mesh
!!
!! dimensionless radial coordinate rho (with 0 < rho < +inf) according to:
!!    x = (gp * sinh(rho-gp3) + gp2) * cos(phi)
!!    y = (gp * sinh(rho-gp3) + gp2) * sin(phi)
!!
!! \extends geometry_cartesian
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_sinhpolar
  USE geometry_cartesian
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE Convert2Cartesian_sinhpolar
     MODULE PROCEDURE Sinhpolar2Cartesian_coords, Sinhpolar2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear_sinhpolar
     MODULE PROCEDURE Cartesian2Sinhpolar_coords, Cartesian2Sinhpolar_vectors
  END INTERFACE
  !> \endcond
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "sinhpolar"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitGeometry_sinhpolar, &
       ScaleFactors_sinhpolar, &
       Radius_sinhpolar, &
       PositionVector_sinhpolar, &
       Convert2Cartesian_sinhpolar, &
       Convert2Curvilinear_sinhpolar, &
       Sinhpolar2Cartesian_coords, &
       Sinhpolar2Cartesian_vectors, &
       Cartesian2Sinhpolar_coords, &
       Cartesian2Sinhpolar_vectors
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_sinhpolar(this,gt,gp,gp2,gp3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    REAL, INTENT(IN) :: gp,gp2,gp3
    !------------------------------------------------------------------------!
    CALL InitGeometry(this,gt,geometry_name)
    CALL SetScale(this,gp,gp2,gp3)
  END SUBROUTINE InitGeometry_sinhpolar
    

  ELEMENTAL SUBROUTINE ScaleFactors_sinhpolar(gp,gp2,gp3,rho,hrho,hphi,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,gp2,gp3,rho
    REAL, INTENT(OUT) :: hrho,hphi,hz
    !------------------------------------------------------------------------!
    hrho = ABS(gp*COSH(rho-gp3))
    hphi = Radius_sinhpolar(gp,gp2,gp3,rho)
    hz   = 1.
  END SUBROUTINE ScaleFactors_sinhpolar


  ELEMENTAL FUNCTION Radius_sinhpolar(gp,gp2,gp3,rho) RESULT(radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,gp2,gp3,rho
    REAL :: radius
    !------------------------------------------------------------------------!
    radius = ABS(gp*SINH(rho-gp3)+gp2)
  END FUNCTION Radius_sinhpolar

  ELEMENTAL SUBROUTINE PositionVector_sinhpolar(gp,gp2,gp3,rho,rx,ry)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,gp2,gp3,rho
    REAL, INTENT(OUT) :: rx,ry
    !------------------------------------------------------------------------!
    rx = Radius_sinhpolar(gp,gp2,gp3,rho)
    ry = 0.0
  END SUBROUTINE PositionVector_sinhpolar

  ! coordinate transformations
  ELEMENTAL SUBROUTINE Sinhpolar2Cartesian_coords(gp,gp2,gp3,rho,phi,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,gp2,gp3,rho,phi
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    REAL :: r
    !------------------------------------------------------------------------!
    r = gp*SINH(rho-gp3)+gp2
    x = r*COS(phi)
    y = r*SIN(phi)
  END SUBROUTINE Sinhpolar2Cartesian_coords

  ELEMENTAL SUBROUTINE Cartesian2Sinhpolar_coords(gp,gp2,gp3,x,y,rho,phi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,gp2,gp3,x,y
    REAL, INTENT(OUT) :: rho,phi
    !------------------------------------------------------------------------!
    REAL :: x1,y1,r1
    !------------------------------------------------------------------------!
    x1 = x
    y1 = y
    r1 = (SQRT(x1*x1+y1*y1)-gp2)/gp
    rho = LOG(r1+SQRT(1.+r1*r1)) & ! = ASINH(r1))
          + gp3
    phi = ATAN2(y,x)
    IF(phi.LT.0.0) THEN
        phi = phi + 2.0*PI
    END IF
  END SUBROUTINE Cartesian2Sinhpolar_coords


  ! vector transformations
  ELEMENTAL SUBROUTINE Sinhpolar2Cartesian_vectors(phi,vrho,vphi,vx,vy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: phi,vrho,vphi
    REAL, INTENT(OUT) :: vx,vy
    !------------------------------------------------------------------------!
    vx = vrho * COS(phi) - vphi * SIN(phi)
    vy = vrho * SIN(phi) + vphi * COS(phi)
  END SUBROUTINE Sinhpolar2Cartesian_vectors
  
  ELEMENTAL SUBROUTINE Cartesian2Sinhpolar_vectors(phi,vx,vy,vrho,vphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: phi,vx,vy
    REAL, INTENT(OUT) :: vrho,vphi
    !------------------------------------------------------------------------!
    vrho = vx * COS(phi) + vy * SIN(phi)
    vphi = -vx * SIN(phi) + vy * COS(phi)
  END SUBROUTINE Cartesian2Sinhpolar_vectors
  
END MODULE geometry_sinhpolar
