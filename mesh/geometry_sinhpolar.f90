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
! define properties of a hyperbolic 2D polar mesh with dimensionless radial
! coordinate rho (with 0 < rho < +inf) according to:
!    x = r0 * sinh(rho) * cos(phi)  
!    y = r0 * sinh(rho) * sin(phi)
!----------------------------------------------------------------------------!
MODULE geometry_sinhpolar
  USE geometry_cartesian
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTERFACE Convert2Cartesian_sinhpolar
     MODULE PROCEDURE Sinhpolar2Cartesian_coords, Sinhpolar2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear_sinhpolar
     MODULE PROCEDURE Cartesian2Sinhpolar_coords, Cartesian2Sinhpolar_vectors
  END INTERFACE
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "sinhpolar"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitGeometry_sinhpolar, &
       ScaleFactors_sinhpolar, &
       Convert2Cartesian_sinhpolar, &
       Convert2Curvilinear_sinhpolar, &
       Sinhpolar2Cartesian_coords, &
       Sinhpolar2Cartesian_vectors, &
       Cartesian2Sinhpolar_coords, &
       Cartesian2Sinhpolar_vectors
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_sinhpolar(this,gt,gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    REAL, INTENT(IN) :: gp
    !------------------------------------------------------------------------!
    CALL InitGeometry(this,gt,geometry_name)
    this%geoparam = gp
  END SUBROUTINE InitGeometry_sinhpolar
    

  ELEMENTAL SUBROUTINE ScaleFactors_sinhpolar(gp,rho,hrho,hphi,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,rho
    REAL, INTENT(OUT) :: hrho,hphi,hz
    !------------------------------------------------------------------------!
    hrho = gp*COSH(rho)
    hphi = gp*SINH(rho)
    hz   = 1.
  END SUBROUTINE ScaleFactors_sinhpolar


  ! coordinate transformations
  ELEMENTAL SUBROUTINE Sinhpolar2Cartesian_coords(gp,rho,phi,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,rho,phi
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    REAL :: r
    !------------------------------------------------------------------------!
    r = gp*SINH(rho)
    x = r*COS(phi)
    y = r*SIN(phi)
  END SUBROUTINE Sinhpolar2Cartesian_coords

  ELEMENTAL SUBROUTINE Cartesian2Sinhpolar_coords(gp,x,y,rho,phi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,x,y
    REAL, INTENT(OUT) :: rho,phi
    !------------------------------------------------------------------------!
    REAL :: x1,y1,r1
    !------------------------------------------------------------------------!
    x1 = x/gp
    y1 = y/gp
    r1 = SQRT(x1*x1+y1*y1)
    rho = LOG(r1+SQRT(1.+r1*r1))  ! = ASINH(r1))
    IF (x1.GT.0.0) THEN
       IF (y1.GE.0.0) THEN
          phi = ATAN(y1/x1)  ! = ATAN(y/x)
       ELSE
          phi = ATAN(y1/x1) + 2*PI
       END IF
    ELSE
       phi = ATAN(y1/x1) + PI
    END IF
  END SUBROUTINE Cartesian2Sinhpolar_coords


  ! vector transformations
  ELEMENTAL SUBROUTINE Sinhpolar2Cartesian_vectors(gp,phi,vrho,vphi,vx,vy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,phi,vrho,vphi
    REAL, INTENT(OUT) :: vx,vy
    !------------------------------------------------------------------------!
    vx = vrho * COS(phi) - vphi * SIN(phi)
    vy = vrho * SIN(phi) + vphi * COS(phi)
  END SUBROUTINE Sinhpolar2Cartesian_vectors
  
  ELEMENTAL SUBROUTINE Cartesian2Sinhpolar_vectors(gp,phi,vx,vy,vrho,vphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,phi,vx,vy
    REAL, INTENT(OUT) :: vrho,vphi
    !------------------------------------------------------------------------!
    vrho = vx * COS(phi) + vy * SIN(phi)
    vphi = -vx * SIN(phi) + vy * COS(phi)
  END SUBROUTINE Cartesian2Sinhpolar_vectors
  
END MODULE geometry_sinhpolar
