!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_bianglespherical.f90                                     #
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
!> \author Björn Sperling
!! \author Tobias Illenseer
!!
!! \brief define properties of spherical mesh with two angles
!!
!! \extends geometry_cartesian
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_bianglespherical
  USE geometry_cartesian
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE Convert2Cartesian_bianglespher
     MODULE PROCEDURE Bianglespherical2Cartesian_coords,&
    Bianglespherical2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear_bianglespher
     MODULE PROCEDURE Cartesian2Bianglespherical_coords,&
    Cartesian2Bianglespherical_vectors
  END INTERFACE
  !> \endcond
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "bianglespherical"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitGeometry_bianglespher, &
       ScaleFactors_bianglespher, &
       Radius_bianglespher, &
       PositionVector_bianglespher, &
       Convert2Cartesian_bianglespher, &
       Convert2Curvilinear_bianglespher, &
       Bianglespherical2Cartesian_coords, &
       Bianglespherical2Cartesian_vectors, &
       Cartesian2Bianglespherical_coords, &
       Cartesian2Bianglespherical_vectors
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_bianglespher(this,gt,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    REAL, INTENT(IN) :: r ! r 
    !------------------------------------------------------------------------!
    CALL InitGeometry(this,gt,geometry_name)
    CALL SetScale(this,r)
  END SUBROUTINE InitGeometry_bianglespher
  
  ELEMENTAL SUBROUTINE ScaleFactors_bianglespher(r,theta,phi,htheta,hphi,hr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r,theta,phi
    REAL, INTENT(OUT) :: hr,htheta,hphi
    !------------------------------------------------------------------------!
    hr     = 1.
    htheta = r     
    hphi   = r*sin(theta)
  END SUBROUTINE ScaleFactors_bianglespher

  ELEMENTAL FUNCTION Radius_bianglespher(r) RESULT(radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r
    REAL :: radius
    !------------------------------------------------------------------------!
    radius = r
  END FUNCTION Radius_bianglespher

  ELEMENTAL SUBROUTINE PositionVector_bianglespher(r,rx,ry,rz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r
    REAL, INTENT(OUT) :: rx,ry,rz
    !------------------------------------------------------------------------!
    rx = Radius_bianglespher(r)
    ry = 0.0
    rz = 0.0
  END SUBROUTINE PositionVector_bianglespher

  ! coordinate transformations
  ELEMENTAL SUBROUTINE Bianglespherical2Cartesian_coords(r,theta,phi,x,y,z)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r,theta,phi
    REAL, INTENT(OUT) :: x,y,z
    !------------------------------------------------------------------------!
    x = r*SIN(theta)*COS(phi)
    y = r*SIN(theta)*SIN(phi)
    z = r*COS(theta)  
  END SUBROUTINE Bianglespherical2Cartesian_coords

  ELEMENTAL SUBROUTINE Cartesian2Bianglespherical_coords(r,x,y,z,theta,phi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r,x,y,z
    REAL, INTENT(OUT) :: theta,phi
    !------------------------------------------------------------------------!
    theta = ACOS(z/r)
    phi = ATAN2(y,x)
  END SUBROUTINE Cartesian2Bianglespherical_coords
 
  ! vector transformations
  ELEMENTAL SUBROUTINE Bianglespherical2Cartesian_vectors(r,theta,phi,&
  vtheta,vphi,vr,vx,vy,vz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r,theta,phi,vr,vtheta,vphi
    REAL, INTENT(OUT) :: vx,vy,vz
    !------------------------------------------------------------------------!
    vx = vr*SIN(theta)*COS(phi) + vtheta*COS(theta)*COS(phi) - vphi*SIN(PHI)
    vy = vr*SIN(theta)*SIN(phi) + vtheta*COS(theta)*SIN(phi) + vphi*COS(PHI)
    vz = vr*COS(theta)          - vphi*SIN(theta)
  END SUBROUTINE Bianglespherical2Cartesian_vectors

  ELEMENTAL SUBROUTINE Cartesian2Bianglespherical_vectors(r,theta,phi,&
  vx,vy,vz,vtheta,vphi,vr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r,theta,phi,vx,vy,vz
    REAL, INTENT(OUT) :: vr,vtheta,vphi
    !------------------------------------------------------------------------!
    vr     = vx*SIN(theta)*COS(phi) + vy*SIN(theta)*SIN(phi) + vz*COS(theta)
    vtheta = vx*COS(theta)*COS(phi) + vy*COS(theta)*SIN(phi) - vz*SIN(theta)
    vphi   = -vx*SIN(phi)           + vy*COS(phi)
  END SUBROUTINE Cartesian2Bianglespherical_vectors

END MODULE geometry_bianglespherical
