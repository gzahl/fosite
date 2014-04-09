!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_elliptic.f90                                             #
!#                                                                           #
!# Copyright (C) 2010 Tobias Illenseer <tillense@astrophysik.uni-kiel.de>    #
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
!! \brief define properties of a 2D bipolar mesh
!!
!! \extends geometry_cartesian
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_elliptic
  USE geometry_cartesian
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE Convert2Cartesian_elliptic
     MODULE PROCEDURE Elliptic2Cartesian_coords, Elliptic2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear_elliptic
     MODULE PROCEDURE Cartesian2Elliptic_coords, Cartesian2Elliptic_vectors
  END INTERFACE
  !> \endcond
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "elliptic"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitGeometry_elliptic, &
       ScaleFactors_elliptic, &
       Radius_elliptic, &
       PositionVector_elliptic, &
       Convert2Cartesian_elliptic, &
       Convert2Curvilinear_elliptic, &
       Elliptic2Cartesian_coords, &
       Elliptic2Cartesian_vectors, &
       Cartesian2Elliptic_coords, &
       Cartesian2Elliptic_vectors
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_elliptic(this,gt,gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    REAL, INTENT(IN) :: gp
    !------------------------------------------------------------------------!
    CALL InitGeometry(this,gt,geometry_name)
    CALL SetScale(this,gp)
  END SUBROUTINE InitGeometry_elliptic
    

  ELEMENTAL SUBROUTINE ScaleFactors_elliptic(gp,mu,nu,hmu,hnu,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,mu,nu
    REAL, INTENT(OUT) :: hmu
    REAL, INTENT(OUT), OPTIONAL :: hnu,hz
    !------------------------------------------------------------------------!
    hmu = gp*SQRT(0.5*(COSH(2*mu)-COS(2*nu)))
    IF (PRESENT(hnu)) hnu = hmu
    IF (PRESENT(hz))  hz = 1.
  END SUBROUTINE ScaleFactors_elliptic

  ELEMENTAL FUNCTION Radius_elliptic(gp,mu,nu) RESULT(radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,mu,nu
    REAL :: radius
    !------------------------------------------------------------------------!
    radius = gp*SQRT(0.5*(COSH(2*mu)+COS(2*nu)))
  END FUNCTION Radius_elliptic

  ELEMENTAL SUBROUTINE PositionVector_elliptic(gp,mu,nu,rx,ry)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,mu,nu
    REAL, INTENT(OUT) :: rx,ry
    REAL              :: hmu
    !------------------------------------------------------------------------!
    CALL ScaleFactors_elliptic(gp,mu,nu,hmu)
    rx = 0.5*gp*gp*SINH(2*mu) / hmu
    ry = 0.5*gp*gp*SIN(2*nu) / hmu
  END SUBROUTINE PositionVector_elliptic


  ! coordinate transformations
  ELEMENTAL SUBROUTINE Elliptic2Cartesian_coords(gp,mu,nu,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,mu,nu
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    x = gp * COSH(mu) * COS(nu)
    y = gp * SINH(mu) * SIN(nu)
  END SUBROUTINE Elliptic2Cartesian_coords

  ELEMENTAL SUBROUTINE Cartesian2Elliptic_coords(gp,x,y,mu,nu)
    USE functions, ONLY : ASINH
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,x,y
    REAL, INTENT(OUT) :: mu,nu
    !------------------------------------------------------------------------!
    REAL :: xpgp,xmgp,tmp1,tmp2
    !------------------------------------------------------------------------!
    xpgp = x+gp
    xmgp = x-gp
    tmp1 = SQRT((y*y+xpgp*xpgp)*(y*y+xmgp*xmgp))
    tmp2 = y*y + xpgp*xmgp
    mu = ASINH(SQRT(0.5*(tmp1+tmp2))/gp)
    ! compute nu (correct only in the first quadrant (x>0,y>0))
    nu = ASIN(SQRT(0.5*(tmp1-tmp2))/gp)
    ! correct for other quadrants to map nu to the interval [0,2*PI)
    IF (y.GT.0.0) THEN
       IF (x.LT.0.0) THEN
          ! second quadrant
          nu = PI - nu
       END IF
    ELSE
       IF (x.LT.0.0) THEN
          ! third quadrant
          nu = PI + nu
       ELSE
          ! fourth quadrant
          nu = 2*PI - nu
       END IF
    END IF
  END SUBROUTINE Cartesian2Elliptic_coords


  ! vector transformations
  ELEMENTAL SUBROUTINE Elliptic2Cartesian_vectors(mu,nu,vmu,vnu,vx,vy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)    :: mu,nu,vmu,vnu
    REAL, INTENT(INOUT) :: vx,vy
    !------------------------------------------------------------------------!
    REAL :: h,invh,a,b
    !------------------------------------------------------------------------!
    ! we can set the geometry parameter to unity, because the  transformation
    ! matrix is multiplied by gp/h and the factor gp cancels out
!CDIR IEXPAND
    CALL ScaleFactors_elliptic(1.0,mu,nu,h)
    a = SINH(mu)*COS(nu)
    b = COSH(mu)*SIN(nu)
    invh = 1./h
    vx   = (a*vmu - b*vnu) * invh
    vy   = (b*vmu + a*vnu) * invh
  END SUBROUTINE Elliptic2Cartesian_vectors

  ELEMENTAL SUBROUTINE Cartesian2Elliptic_vectors(mu,nu,vx,vy,vmu,vnu)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: mu,nu,vx,vy
    REAL, INTENT(OUT) :: vmu,vnu
    !------------------------------------------------------------------------!
    REAL :: h,invh,a,b
    !------------------------------------------------------------------------!
    ! we can set the geometry parameter to unity, because the  transformation
    ! matrix is multiplied by gp/h and the factor gp cancels out
!CDIR IEXPAND
    CALL ScaleFactors_elliptic(1.0,mu,nu,h)
    a = SINH(mu)*COS(nu)
    b = COSH(mu)*SIN(nu)
    invh = 1./h
    vmu  =  (a*vx + b*vy) * invh
    vnu  = (-b*vx + a*vy) * invh
  END SUBROUTINE Cartesian2Elliptic_vectors
  
END MODULE geometry_elliptic
