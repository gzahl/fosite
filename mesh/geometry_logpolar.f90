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
! define properties of a logarithmic 2D polar mesh
!    x = r0 * exp(r) * cos(phi)  
!    y = r0 * exp(r) * sin(phi)
!----------------------------------------------------------------------------!
MODULE geometry_logpolar
  USE geometry_cartesian
  USE geometry_oblatespheroidal
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTERFACE Convert2Cartesian_logpolar
     MODULE PROCEDURE Convert2Cartesian_coords, Convert2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear_logpolar
     MODULE PROCEDURE Convert2Curvilinear_coords, Convert2Curvilinear_vectors
  END INTERFACE
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "logpolar"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitGeometry_logpolar, &
       ScaleFactors_logpolar, &
       Convert2Cartesian_logpolar, &
       Convert2Curvilinear_logpolar
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_logpolar(this,gt,gs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    REAL, INTENT(IN) :: gs
    !------------------------------------------------------------------------!
    CALL InitGeometry(this,gt,geometry_name)
    this%scalefactor = gs
  END SUBROUTINE InitGeometry_logpolar
    

  ELEMENTAL SUBROUTINE ScaleFactors_logpolar(gs,r,hr,hphi,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gs,r
    REAL, INTENT(OUT) :: hr,hphi,hz
    !------------------------------------------------------------------------!
    hr   = gs*EXP(r)
    hphi = hr
    hz   = 1.
  END SUBROUTINE ScaleFactors_logpolar


  ! coordinate transformation
  ! logpolar -> cartesian
  ELEMENTAL SUBROUTINE Convert2Cartesian_coords(gs,r,phi,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gs,r,phi
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    REAL :: r0
    !------------------------------------------------------------------------!
    r0 = gs*EXP(r)
    x = r0*COS(phi)
    y = r0*SIN(phi)
  END SUBROUTINE Convert2Cartesian_coords


  ! cartesian -> logpolar
  ELEMENTAL SUBROUTINE Convert2Curvilinear_coords(gs,x,y,r,phi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gs,x,y
    REAL, INTENT(OUT) :: r,phi
    !------------------------------------------------------------------------!
    REAL :: x1,y1
    !------------------------------------------------------------------------!
    x1 = x/gs
    y1 = y/gs
    r = 0.5*LOG(x1*x1+y1*y1)
    phi = ATAN(y1/x1) ! = ATAN(y/x)
  END SUBROUTINE Convert2Curvilinear_coords


  ! vector transformation
  ! logpolar -> cartesian  
  ELEMENTAL SUBROUTINE Convert2Cartesian_vectors(gs,phi,vr,vphi,vx,vy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gs,phi,vr,vphi
    REAL, INTENT(OUT) :: vx,vy
    !------------------------------------------------------------------------!
    vx = vr * COS(phi) - vphi * SIN(phi)
    vy = vr * SIN(phi) + vphi * COS(phi)
  END SUBROUTINE Convert2Cartesian_vectors


  ! cartesian -> logpolar
  ELEMENTAL SUBROUTINE Convert2Curvilinear_vectors(gs,phi,vx,vy,vr,vphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gs,phi,vx,vy
    REAL, INTENT(OUT) :: vr,vphi
    !------------------------------------------------------------------------!
    vr   = vx * COS(phi) + vy * SIN(phi)
    vphi = -vx * SIN(phi) + vy * COS(phi)
  END SUBROUTINE Convert2Curvilinear_vectors
  
END MODULE geometry_logpolar
