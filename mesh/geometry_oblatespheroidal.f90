!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_oblatespheroidal.f90                                     #
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
! define properties of a 2.5D oblate spheroidal mesh
!----------------------------------------------------------------------------!
MODULE geometry_oblatespheroidal
  USE geometry_common
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTERFACE Convert2Cartesian_oblatespheroidal
     MODULE PROCEDURE Convert2Cartesian_coords, Convert2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear_oblatespheroidal
     MODULE PROCEDURE Convert2Curvilinear_coords, Convert2Curvilinear_vectors
  END INTERFACE
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "oblate spheroidal"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitGeometry_oblatespheroidal, &
       GetScale, &
       ScaleFactors_oblatespheroidal, &
       Convert2Cartesian_oblatespheroidal, &
       Convert2Curvilinear_oblatespheroidal
  !--------------------------------------------------------------------------!

CONTAINS

  PURE SUBROUTINE InitGeometry_oblatespheroidal(this,gt,gs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    REAL, INTENT(IN) :: gs
    !------------------------------------------------------------------------!
    CALL InitGeometry(this,gt,geometry_name)
    this%scalefactor = gs
  END SUBROUTINE InitGeometry_oblatespheroidal
    

  PURE FUNCTION GetScale(this) RESULT(gs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    REAL :: gs
    !------------------------------------------------------------------------!
    gs = this%scalefactor
  END FUNCTION GetScale


  ELEMENTAL SUBROUTINE ScaleFactors_oblatespheroidal(gs,xi,eta,hxi,heta,hphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gs,xi,eta
    REAL, INTENT(OUT) :: hxi,heta,hphi
    !------------------------------------------------------------------------!

    hxi  = gs*SQRT(SINH(xi)**2+SIN(eta)**2)
    heta = hxi
    hphi = gs*COSH(xi)*COS(eta)
  END SUBROUTINE ScaleFactors_oblatespheroidal

  
  ! coordinate transformation
  ! oblate spheroidal -> cartesian
  ELEMENTAL SUBROUTINE Convert2Cartesian_coords(gs,xi,eta,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gs,xi,eta
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    x = gs*COSH(xi)*COS(eta)
    y = gs*SINH(xi)*SIN(eta)
  END SUBROUTINE Convert2Cartesian_coords


  ELEMENTAL SUBROUTINE Convert2Cartesian_vectors(gs,xi,eta,vxi,veta,vx,vy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gs,xi,eta,vxi,veta
    REAL, INTENT(OUT) :: vx,vy
    !------------------------------------------------------------------------!
    ! local variables
    REAL :: hh
    !------------------------------------------------------------------------!
    hh = SQRT(SINH(xi)**2+SIN(eta)**2)
    vx = (vxi*SINH(xi)*COS(eta) - veta*COSH(xi)*SIN(eta))/hh
    vy = (vxi*COSH(xi)*SIN(eta) + veta*SINH(xi)*COS(eta))/hh
  END SUBROUTINE Convert2Cartesian_vectors


  ! vector transformation
  ! cartesian -> oblate spheroidal
  ELEMENTAL SUBROUTINE Convert2Curvilinear_coords(gs,x,y,xi,eta)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gs,x,y
    REAL, INTENT(OUT) :: xi,eta
    !------------------------------------------------------------------------!
    ! local variables
    REAL :: r2,a2,tmp,sinhxi,sinhxi2
    !------------------------------------------------------------------------!
    a2=gs*gs
    r2=.5*(a2 - x*x - y*y)
    sinhxi2=(SQRT(1.+(gs*y/r2)**2)-1.)*r2/a2
    sinhxi=SQRT(sinhxi2)
    xi=LOG(sinhxi+SQRT(sinhxi2+1.))             ! = ASINH(SQRT(sinhxi2))
    eta=ASIN(y/(gs*sinhxi))
  END SUBROUTINE Convert2Curvilinear_coords


  ELEMENTAL SUBROUTINE Convert2Curvilinear_vectors(gs,xi,eta,vx,vy,vxi,veta)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gs,xi,eta,vx,vy
    REAL, INTENT(OUT) :: vxi,veta
    !------------------------------------------------------------------------!
    ! local variables
    REAL :: hh
    !------------------------------------------------------------------------!
    hh = SQRT(SINH(xi)**2+SIN(eta)**2) 
    vxi = (vx*SINH(xi)*COS(eta) + vy*COSH(xi)*SIN(eta))/hh
    veta = (-vx*COSH(xi)*SIN(eta) + vy*SINH(xi)*COS(eta))/hh 
  END SUBROUTINE Convert2Curvilinear_vectors

END MODULE geometry_oblatespheroidal
