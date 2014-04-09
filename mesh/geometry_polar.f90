!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_polar.f90                                                #
!#                                                                           #
!# Copyright (C) 2007 Tobias Illenseer <tillense@ita.uni-heidelberg.de>      #
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
! define properties of a 2D polar mesh
!----------------------------------------------------------------------------!
MODULE geometry_polar
  USE geometry_common
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTERFACE Convert2Cartesian_polar
     MODULE PROCEDURE Convert2Cartesian_coords, Convert2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear_polar
     MODULE PROCEDURE Convert2Curvilinear_coords, Convert2Curvilinear_vectors
  END INTERFACE
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "polar"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitGeometry_polar, &
       ScaleFactors_polar, &
       Convert2Cartesian_polar, &
       Convert2Curvilinear_polar
  !--------------------------------------------------------------------------!

CONTAINS

  PURE SUBROUTINE InitGeometry_polar(this,gt)
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


  ! coordinate transformation
  ! polar -> cartesian
  ELEMENTAL SUBROUTINE Convert2Cartesian_coords(r,phi,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r,phi
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    x = r*COS(phi)
    y = r*SIN(phi)
  END SUBROUTINE Convert2Cartesian_coords


  ! cartesian -> polar
  ELEMENTAL SUBROUTINE Convert2Curvilinear_coords(x,y,r,phi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: x,y
    REAL, INTENT(OUT) :: r,phi
    !------------------------------------------------------------------------!
    r = SQRT(x*x+y*y)
    phi = ATAN(y/x)
  END SUBROUTINE Convert2Curvilinear_coords


  ! vector transformation
  ! polar -> cartesian  
  ELEMENTAL SUBROUTINE Convert2Cartesian_vectors(phi,vr,vphi,vx,vy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: phi,vr,vphi
    REAL, INTENT(OUT) :: vx,vy
    !------------------------------------------------------------------------!
    vx = vr * COS(phi) - vphi * SIN(phi)
    vy = vr * SIN(phi) + vphi * COS(phi)
  END SUBROUTINE Convert2Cartesian_vectors


  ! cartesian -> polar
  ELEMENTAL SUBROUTINE Convert2Curvilinear_vectors(phi,vx,vy,vr,vphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: phi,vx,vy
    REAL, INTENT(OUT) :: vr,vphi
    !------------------------------------------------------------------------!
    vr   = vx * COS(phi) + vy * SIN(phi)
    vphi = -vx * SIN(phi) + vy * COS(phi)
  END SUBROUTINE Convert2Curvilinear_vectors
  
END MODULE geometry_polar
