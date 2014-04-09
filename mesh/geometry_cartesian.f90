!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_cartesian.f90                                            #
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
! define properties of a 2D cartesian mesh
!----------------------------------------------------------------------------!
MODULE geometry_cartesian
  USE geometry_common
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "cartesian"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Geometry_TYP, &
       ! constants
       PI, &
       ! methods
       InitGeometry, &
       InitGeometry_cartesian, &
       ScaleFactors_cartesian, &
       GetType, &
       GetName, &
       GetRank, &
       GetScale, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_cartesian(this,gt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    !------------------------------------------------------------------------!
    CALL InitGeometry(this,gt,geometry_name)
  END SUBROUTINE InitGeometry_cartesian
    

  ELEMENTAL SUBROUTINE ScaleFactors_cartesian(hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(OUT) :: hx,hy,hz
    !------------------------------------------------------------------------!
    ! scale factors are unity
    hx = 1.
    hy = 1.
    hz = 1.
  END SUBROUTINE ScaleFactors_cartesian

END MODULE geometry_cartesian
