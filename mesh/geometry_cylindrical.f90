!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_cylindrical.f90                                          #
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
!> \author Tobias Illenseer
!!
!! \brief define properties of a 2.5D cylindrical mesh
!!
!! \extends geometry_cartesian
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_cylindrical
  USE geometry_cartesian, Radius_cylindrical => Radius_cartesian
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "cylindrical"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
      InitGeometry_cylindrical, &
      ScaleFactors_cylindrical, &
      Radius_cylindrical, &
      Convert2Cartesian_cylindrical, &
      Convert2Curvilinear_cylindrical
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_cylindrical(this,gt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    !------------------------------------------------------------------------!
    CALL InitGeometry(this,gt,geometry_name)
  END SUBROUTINE InitGeometry_cylindrical
    

  ELEMENTAL SUBROUTINE ScaleFactors_cylindrical(r,hz,hr,hphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r
    REAL, INTENT(OUT) :: hz,hr,hphi
    !------------------------------------------------------------------------!
    hz   = 1.
    hr   = 1.
    hphi = r
  END SUBROUTINE ScaleFactors_cylindrical


  ! coordinate/vector transformation
  ! cylindrical -> cartesian
  ! this works for both, vectors and coords
  ELEMENTAL SUBROUTINE Convert2Cartesian_cylindrical(z,r,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: z,r
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    x = r
    y = z
  END SUBROUTINE Convert2Cartesian_cylindrical


  ! coordinate/vector transformation
  ! cartesian -> cylindrical
  ! this works for both, vectors and coords
  ELEMENTAL SUBROUTINE Convert2Curvilinear_cylindrical(x,y,z,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: x,y
    REAL, INTENT(OUT) :: z,r
    !------------------------------------------------------------------------!
    z=y
    r=x
  END SUBROUTINE Convert2Curvilinear_cylindrical

END MODULE geometry_cylindrical
