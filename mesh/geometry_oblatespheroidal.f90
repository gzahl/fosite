!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_oblatespheroidal.f90                                     #
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
! define properties of a 2.5D oblate spheroidal mesh
!----------------------------------------------------------------------------!
MODULE geometry_oblatespheroidal
  USE geometry_cartesian
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTERFACE Convert2Cartesian_oblatespher
     MODULE PROCEDURE Oblatespher2Cartesian_coords, Oblatespher2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear_oblatespher
     MODULE PROCEDURE Cartesian2Oblatespher_coords, Cartesian2Oblatespher_vectors
  END INTERFACE
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "oblate spheroidal"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitGeometry_oblatespher, &
       ScaleFactors_oblatespher, &
       Convert2Cartesian_oblatespher, &
       Convert2Curvilinear_oblatespher, & 
       Oblatespher2Cartesian_coords, &
       Oblatespher2Cartesian_vectors, &
       Cartesian2Oblatespher_coords, &
       Cartesian2Oblatespher_vectors
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_oblatespher(this,gt,gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    REAL, INTENT(IN) :: gp
    !------------------------------------------------------------------------!
    CALL InitGeometry(this,gt,geometry_name)
    this%geoparam = gp
  END SUBROUTINE InitGeometry_oblatespher
    

  ELEMENTAL SUBROUTINE ScaleFactors_oblatespher(gp,xi,eta,hxi,heta,hphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,xi,eta
    REAL, INTENT(OUT) :: hxi,heta,hphi
    !------------------------------------------------------------------------!
    hxi  = gp*SQRT(SINH(xi)**2+SIN(eta)**2)
    heta = hxi
    hphi = gp*COSH(xi)*COS(eta)
  END SUBROUTINE ScaleFactors_oblatespher

  
  ! coordinate transformations
  ELEMENTAL SUBROUTINE Oblatespher2Cartesian_coords(gp,xi,eta,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,xi,eta
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    x = gp*COSH(xi)*COS(eta)
    y = gp*SINH(xi)*SIN(eta)
  END SUBROUTINE Oblatespher2Cartesian_coords

  ELEMENTAL SUBROUTINE Cartesian2Oblatespher_coords(gp,x,y,xi,eta)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,x,y
    REAL, INTENT(OUT) :: xi,eta
    !------------------------------------------------------------------------!
    ! local variables
    REAL :: xa,ya,coshxi
    !------------------------------------------------------------------------!
    xa = x/gp
    ya = y/gp
    coshxi = 0.5*(SQRT((1.0+xa)**2+ya*ya) + SQRT((1.0-xa)**2+ya*ya))
    xi = LOG(coshxi+SQRT(coshxi*coshxi-1.0))    ! = ACOSH(coshxi)
    eta=ACOS(xa/coshxi)
  END SUBROUTINE Cartesian2Oblatespher_coords

  ! vector transformations
  ELEMENTAL SUBROUTINE Oblatespher2Cartesian_vectors(gp,xi,eta,vxi,veta,vx,vy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,xi,eta,vxi,veta
    REAL, INTENT(OUT) :: vx,vy
    !------------------------------------------------------------------------!
    ! local variables
    REAL :: hh
    !------------------------------------------------------------------------!
    hh = SQRT(SINH(xi)**2+SIN(eta)**2)
    vx = (vxi*SINH(xi)*COS(eta) - veta*COSH(xi)*SIN(eta))/hh
    vy = (vxi*COSH(xi)*SIN(eta) + veta*SINH(xi)*COS(eta))/hh
  END SUBROUTINE Oblatespher2Cartesian_vectors

  ELEMENTAL SUBROUTINE Cartesian2Oblatespher_vectors(xi,eta,vx,vy,vxi,veta)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: xi,eta,vx,vy
    REAL, INTENT(OUT) :: vxi,veta
    !------------------------------------------------------------------------!
    ! local variables
    REAL :: hh
    !------------------------------------------------------------------------!
    hh = SQRT(SINH(xi)**2+SIN(eta)**2) 
    vxi = (vx*SINH(xi)*COS(eta) + vy*COSH(xi)*SIN(eta))/hh
    veta = (-vx*COSH(xi)*SIN(eta) + vy*SINH(xi)*COS(eta))/hh 
  END SUBROUTINE Cartesian2Oblatespher_vectors

END MODULE geometry_oblatespheroidal
