!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_logpolar.f90                                             #
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
! define properties of a logarithmic 2D polar mesh with dimensionless radial
! coordinate rho according to:
!    x = r0 * exp(rho) * cos(phi)  
!    y = r0 * exp(rho) * sin(phi)
!----------------------------------------------------------------------------!
MODULE geometry_logpolar
  USE geometry_cartesian
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTERFACE Convert2Cartesian_logpolar
     MODULE PROCEDURE Logpolar2Cartesian_coords, Logpolar2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear_logpolar
     MODULE PROCEDURE Cartesian2Logpolar_coords, Cartesian2Logpolar_vectors
  END INTERFACE
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "logpolar"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitGeometry_logpolar, &
       ScaleFactors_logpolar, &
       Convert2Cartesian_logpolar, &
       Convert2Curvilinear_logpolar, &
       Logpolar2Cartesian_coords, &
       Logpolar2Cartesian_vectors, &
       Cartesian2Logpolar_coords, &
       Cartesian2Logpolar_vectors
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_logpolar(this,gt,gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    REAL, INTENT(IN) :: gp
    !------------------------------------------------------------------------!
    CALL InitGeometry(this,gt,geometry_name)
    this%geoparam = gp
  END SUBROUTINE InitGeometry_logpolar
    

  ELEMENTAL SUBROUTINE ScaleFactors_logpolar(gp,rho,hrho,hphi,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,rho
    REAL, INTENT(OUT) :: hrho,hphi,hz
    !------------------------------------------------------------------------!
    hrho = gp*EXP(rho)
    hphi = hrho
    hz   = 1.
  END SUBROUTINE ScaleFactors_logpolar


  ! coordinate transformations
  ELEMENTAL SUBROUTINE Logpolar2Cartesian_coords(gp,rho,phi,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,rho,phi
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    REAL :: r
    !------------------------------------------------------------------------!
    r = gp*EXP(rho)
    x = r*COS(phi)
    y = r*SIN(phi)
  END SUBROUTINE Logpolar2Cartesian_coords

  ELEMENTAL SUBROUTINE Cartesian2Logpolar_coords(gp,x,y,rho,phi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,x,y
    REAL, INTENT(OUT) :: rho,phi
    !------------------------------------------------------------------------!
    REAL :: x1,y1
    !------------------------------------------------------------------------!
    x1 = x/gp
    y1 = y/gp
    rho = 0.5*LOG(x1*x1+y1*y1)
    phi = ATAN2(y1,x1)+pi
  END SUBROUTINE Cartesian2Logpolar_coords


  ! vector transformations
  ELEMENTAL SUBROUTINE Logpolar2Cartesian_vectors(gp,phi,vrho,vphi,vx,vy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,phi,vrho,vphi
    REAL, INTENT(OUT) :: vx,vy
    !------------------------------------------------------------------------!
    vx = vrho * COS(phi) - vphi * SIN(phi)
    vy = vrho * SIN(phi) + vphi * COS(phi)
  END SUBROUTINE Logpolar2Cartesian_vectors
  
  ELEMENTAL SUBROUTINE Cartesian2Logpolar_vectors(gp,phi,vx,vy,vrho,vphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,phi,vx,vy
    REAL, INTENT(OUT) :: vrho,vphi
    !------------------------------------------------------------------------!
    vrho = vx * COS(phi) + vy * SIN(phi)
    vphi = -vx * SIN(phi) + vy * COS(phi)
  END SUBROUTINE Cartesian2Logpolar_vectors
  
END MODULE geometry_logpolar
