!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_spherical.f90                                            #
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
! define properties of a 2.5D spherical mesh
!----------------------------------------------------------------------------!
MODULE geometry_spherical
  USE geometry_cartesian
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTERFACE Convert2Cartesian_spherical
     MODULE PROCEDURE Spherical2Cartesian_coords, Spherical2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear_spherical
     MODULE PROCEDURE Cartesian2Spherical_coords, Cartesian2Spherical_vectors
  END INTERFACE
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "spherical"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! methods
       InitGeometry_spherical, &
       ScaleFactors_spherical, &
       Convert2Cartesian_spherical, &
       Convert2Curvilinear_spherical, &
       Spherical2Cartesian_coords, &
       Spherical2Cartesian_vectors, &
       Cartesian2Spherical_coords, &
       Cartesian2Spherical_vectors
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_spherical(this,gt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    !------------------------------------------------------------------------!
    CALL InitGeometry(this,gt,geometry_name)
    this%spherical_like = .TRUE.
  END SUBROUTINE InitGeometry_spherical
    

  ELEMENTAL SUBROUTINE ScaleFactors_spherical(r,theta,hr,htheta,hphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r,theta
    REAL, INTENT(OUT) :: hr,htheta,hphi
    !------------------------------------------------------------------------!
    hr     = 1.
    htheta = r
    hphi   = r*sin(theta)
  END SUBROUTINE ScaleFactors_spherical

  
  ! coordinate transformations
  ELEMENTAL SUBROUTINE Spherical2Cartesian_coords(r,theta,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r,theta
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    x = r*SIN(theta)
    y = r*COS(theta)
  END SUBROUTINE Spherical2Cartesian_coords

  ELEMENTAL SUBROUTINE Cartesian2Spherical_coords(x,y,r,theta)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: x,y
    REAL, INTENT(OUT) :: r,theta
    !------------------------------------------------------------------------!
    r = SQRT(x*x+y*y)
    theta = ACOS(y/r)
  END SUBROUTINE Cartesian2Spherical_coords


  ! vector transformations
  ELEMENTAL SUBROUTINE Spherical2Cartesian_vectors(gp,theta,vr,vtheta,vx,vy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,theta,vr,vtheta
    REAL, INTENT(OUT) :: vx,vy
    !------------------------------------------------------------------------!
    vx = vr*SIN(theta) + vtheta*COS(theta)
    vy = vr*COS(theta) - vtheta*SIN(theta)
  END SUBROUTINE Spherical2Cartesian_vectors

  ELEMENTAL SUBROUTINE Cartesian2Spherical_vectors(gp,theta,vx,vy,vr,vtheta)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,theta,vx,vy
    REAL, INTENT(OUT) :: vr,vtheta
    !------------------------------------------------------------------------!
    vr = vx*SIN(theta) + vy*COS(theta)
    vtheta = vx*COS(theta) - vy*SIN(theta)
  END SUBROUTINE Cartesian2Spherical_vectors

END MODULE geometry_spherical
