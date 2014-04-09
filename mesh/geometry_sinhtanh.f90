!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_sinhtanhpolar.f90                                        #
!#                                                                           #
!# Copyright (C) 2008-2011                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
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
!> \author Tobias Illenseer
!! \author Manuel Jung
!!
!! \brief define properties of a hyperbolic 2D polar mesh
!!
!! dimensionless radial coordinate rho (with 0 < rho < +inf) according to:
!!    f(rho) = r0 * sinh(r-gp2) + gp3
!!    g(phi) = alpha * pi * tanh(phi) + pi
!!    x = f(rho) * cos(g(phi))
!!    y = f(rho) * sin(g(phi))
!! here alpha is scaling constant > 1, because we need a mapping of tanh to
!! -pi..pi, which would otherwise only be acchieved for -oo..oo as arguments.
!!
!! \extends geometry_cartesian
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_sinhtanh
  USE geometry_cartesian
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE Convert2Cartesian_sinhtanh
     MODULE PROCEDURE Sinhtanh2Cartesian_coords, &
                      Sinhtanh2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear_sinhtanh
     MODULE PROCEDURE Cartesian2Sinhtanh_coords, &
                      Cartesian2Sinhtanh_vectors
  END INTERFACE
  !> \endcond
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "sinhtanhpolar"
  REAL, PARAMETER              :: alpha = 1.01
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitGeometry_sinhtanh, &
       ScaleFactors_sinhtanh, &
       Radius_sinhtanh, &
       PositionVector_sinhtanh, &
       Convert2Cartesian_sinhtanh, &
       Convert2Curvilinear_sinhtanh, &
       Sinhtanh2Cartesian_coords, &
       Sinhtanh2Cartesian_vectors, &
       Cartesian2Sinhtanh_coords, &
       Cartesian2Sinhtanh_vectors
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_sinhtanh(this,gt,gp,gp2,gp3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    REAL, INTENT(IN) :: gp,gp2,gp3
    !------------------------------------------------------------------------!
    CALL InitGeometry(this,gt,geometry_name)
    CALL SetScale(this,gp,gp2,gp3)
  END SUBROUTINE InitGeometry_sinhtanh
    

  ELEMENTAL SUBROUTINE ScaleFactors_sinhtanh(gp,gp2,gp3,rho,phi,hrho,hphi,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,gp2,gp3,rho,phi
    REAL, INTENT(OUT) :: hrho,hphi,hz
    !------------------------------------------------------------------------!
    hrho = gp*COSH(rho-gp3)
    hphi = (gp*SINH(rho-gp3)+gp2)*alpha*PI/(COSH(phi))**2
    hz   = 1.
  END SUBROUTINE ScaleFactors_sinhtanh


  ELEMENTAL FUNCTION Radius_sinhtanh(gp,gp2,gp3,rho) RESULT(radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,gp2,gp3,rho
    REAL :: radius
    !------------------------------------------------------------------------!
    radius = ABS(gp*SINH(rho-gp3)+gp2)
  END FUNCTION Radius_sinhtanh

  ELEMENTAL SUBROUTINE PositionVector_sinhtanh(gp,gp2,gp3,rho,rx,ry)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,gp2,gp3,rho
    REAL, INTENT(OUT) :: rx,ry
    !------------------------------------------------------------------------!
    rx = Radius_sinhtanh(gp,gp2,gp3,rho)
    ry = 0.0
  END SUBROUTINE PositionVector_sinhtanh

  ! coordinate transformations
  ELEMENTAL SUBROUTINE Sinhtanh2Cartesian_coords(gp,gp2,gp3,rho,phi,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,gp2,gp3,rho,phi
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    REAL :: r, sigma
    !------------------------------------------------------------------------!
    r = gp*SINH(rho-gp3)+gp2
    sigma = alpha * PI * TANH(phi) + PI
    x = r*COS(sigma)
    y = r*SIN(sigma)
  END SUBROUTINE Sinhtanh2Cartesian_coords

  ELEMENTAL SUBROUTINE Cartesian2Sinhtanh_coords(gp,gp2,gp3,x,y,rho,phi)
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
    rho = LOG(r1+SQRT(1.+r1*r1)) +gp3 ! = ASINH(r1))
    phi = ATAN2(y,x)
    IF(phi.LT.0.0) THEN
        ! phi = ATANH(phi/PI - 1.0)  !ATANH is only available in gfortran
        ! for x<1: atanh(x) = 0.5 * ln((1+x)/(1-x))
        x1 = phi/PI + 1.0
        phi = 0.5*LOG((1.0+x1)/(1.0-x1))
    END IF
  END SUBROUTINE Cartesian2Sinhtanh_coords


  ! vector transformations
  ELEMENTAL SUBROUTINE Sinhtanh2Cartesian_vectors(gp,gp2,gp3,phi,vrho,vphi,vx,vy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,gp2,gp3,phi,vrho,vphi
    REAL, INTENT(OUT) :: vx,vy
    !------------------------------------------------------------------------!
    REAL              :: sigma
    !------------------------------------------------------------------------!
    sigma = alpha * PI * TANH(phi) + PI
    vx = vrho * COS(sigma) - vphi * SIN(sigma)
    vy = vrho * SIN(sigma) + vphi * COS(sigma)
  END SUBROUTINE Sinhtanh2Cartesian_vectors
  
  ELEMENTAL SUBROUTINE Cartesian2Sinhtanh_vectors(gp,phi,vx,vy,vrho,vphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,phi,vx,vy
    REAL, INTENT(OUT) :: vrho,vphi
    !------------------------------------------------------------------------!
    REAL              :: sigma
    !------------------------------------------------------------------------!
    sigma = alpha * PI * TANH(phi) + PI
    vrho = vx * COS(sigma) + vy * SIN(sigma)
    vphi = -vx * SIN(sigma) + vy * COS(sigma)
  END SUBROUTINE Cartesian2Sinhtanh_vectors
  
END MODULE geometry_sinhtanh
