!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_tancylindrical.f90                                       #
!#                                                                           #
!# Copyright (C) 2009 Tobias Illenseer <tillense@astrophysik.uni-kiel.de>    #
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
!! \brief define properties of a 2.5D tancylindrical mesh
!!
!! dimensionless vertical coordinate zeta (with -pi/2 < zeta < +pi/2)
!! according to:
!!   x = r
!!   z = z0 * tan(zeta)
!!
!! \extends geometry_cartesian
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_tancylindrical
  USE geometry_cartesian
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE Convert2Cartesian_tancyl
     MODULE PROCEDURE Tancyl2Cartesian_coords, Tancyl2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear_tancyl
     MODULE PROCEDURE Cartesian2Tancyl_coords, Cartesian2Tancyl_vectors
  END INTERFACE
  !> \endcond
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "tancylindrical"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
      InitGeometry_tancyl, &
      GetScale, &
      ScaleFactors_tancyl, &
      Radius_tancyl, &
      PositionVector_tancyl, &
      Convert2Cartesian_tancyl, &
      Convert2Curvilinear_tancyl, &
      Tancyl2Cartesian_coords, &
      Tancyl2Cartesian_vectors, &
      Cartesian2Tancyl_coords, &
      Cartesian2Tancyl_vectors
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_tancyl(this,gt,gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    REAL, INTENT(IN) :: gp
    !------------------------------------------------------------------------!
    CALL InitGeometry(this,gt,geometry_name)
    CALL SetScale(this,gp)
  END SUBROUTINE InitGeometry_tancyl
    

  ELEMENTAL SUBROUTINE ScaleFactors_tancyl(gp,zeta,r,hzeta,hr,hphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,zeta,r
    REAL, INTENT(OUT) :: hzeta,hr,hphi
    !------------------------------------------------------------------------!
    hzeta = gp/COS(zeta)**2
    hr = 1.
    hphi = r
  END SUBROUTINE ScaleFactors_tancyl


  ELEMENTAL FUNCTION Radius_tancyl(gp,zeta,r) RESULT(radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,zeta,r
    REAL              :: radius
    !------------------------------------------------------------------------!
    radius = SQRT((gp*TAN(zeta))**2+r**2)
  END FUNCTION Radius_tancyl


  ELEMENTAL SUBROUTINE PositionVector_tancyl(gp,zeta,r,rx,ry)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,zeta,r
    REAL, INTENT(OUT) :: rx,ry
    !------------------------------------------------------------------------!
    rx = gp*TAN(zeta)
    ry = r
  END SUBROUTINE PositionVector_tancyl


  ! coordinate transformations
  ELEMENTAL SUBROUTINE Tancyl2Cartesian_coords(gp,zeta,r,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,zeta,r
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    x = r
    y = gp*TAN(zeta)
  END SUBROUTINE Tancyl2Cartesian_coords

  ELEMENTAL SUBROUTINE Cartesian2Tancyl_coords(gp,x,y,zeta,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,x,y
    REAL, INTENT(OUT) :: zeta,r
    !------------------------------------------------------------------------!
    zeta = ATAN(y/gp)
    r = x
  END SUBROUTINE Cartesian2Tancyl_coords


  ! vector transformations
  ELEMENTAL SUBROUTINE Tancyl2Cartesian_vectors(vzeta,vr,vx,vy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: vzeta,vr
    REAL, INTENT(OUT) :: vx,vy
    !------------------------------------------------------------------------!
    vx = vr
    vy = vzeta
  END SUBROUTINE Tancyl2Cartesian_vectors

  ! cartesian -> tancylindrical
  ELEMENTAL SUBROUTINE Cartesian2Tancyl_vectors(vx,vy,vzeta,vr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: vx,vy
    REAL, INTENT(OUT) :: vzeta,vr
    !------------------------------------------------------------------------!
    vzeta = vy
    vr = vx
  END SUBROUTINE Cartesian2Tancyl_vectors

END MODULE geometry_tancylindrical
