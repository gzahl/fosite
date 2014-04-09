!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_generic.f90                                              #
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
! generic module for geometrical properties
!----------------------------------------------------------------------------!
MODULE geometry_generic
  USE geometry_common, ONLY : Geometry_TYP
  USE geometry_cartesian
  USE geometry_polar
  USE geometry_cylindrical
  USE geometry_spherical
  USE geometry_oblatespheroidal
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE ScaleFactors
     MODULE PROCEDURE ScaleFactors_1, ScaleFactors_2
  END INTERFACE
  INTERFACE Convert2Cartesian
     MODULE PROCEDURE Convert2Cartesian_coords
     MODULE PROCEDURE Convert2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear
     MODULE PROCEDURE Convert2Curvilinear_coords
     MODULE PROCEDURE Convert2Curvilinear_vectors
  END INTERFACE
  !--------------------------------------------------------------------------!
  INTEGER, PARAMETER :: CARTESIAN         = 1
  INTEGER, PARAMETER :: POLAR             = 2
  INTEGER, PARAMETER :: CYLINDRICAL       = 3
  INTEGER, PARAMETER :: SPHERICAL         = 4
  INTEGER, PARAMETER :: OBLATE_SPHEROIDAL = 5
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Geometry_TYP, &
       ! constants
       CARTESIAN, POLAR, CYLINDRICAL, SPHERICAL, OBLATE_SPHEROIDAL, &
       ! methods
       InitGeometry, &
       GetType, &
       GetName, &
       GetScale, &
       ScaleFactors, &
       Convert2Cartesian, &
       Convert2Curvilinear
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry(this,gt,gs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    REAL, INTENT(IN), OPTIONAL :: gs
    !------------------------------------------------------------------------!

    SELECT CASE(gt)
    CASE(CARTESIAN)
       CALL InitGeometry_cartesian(this,gt)
    CASE(POLAR)
       CALL InitGeometry_polar(this,gt)
    CASE(CYLINDRICAL)
       CALL InitGeometry_cylindrical(this,gt)
    CASE(SPHERICAL)
       CALL InitGeometry_spherical(this,gt)
    CASE(OBLATE_SPHEROIDAL)
       IF (PRESENT(gs)) THEN
          CALL InitGeometry_oblatespheroidal(this,gt,gs)
       ELSE
          CALL InitGeometry_oblatespheroidal(this,gt,1.0)
       END IF
    CASE DEFAULT
       PRINT *, "ERROR in InitGeometry: unknown geometry"
       STOP
    END SELECT

    ! print some information
    PRINT "(A,A)", " GEOMETRY-> coordinates:       ", TRIM(GetName(this))
  END SUBROUTINE InitGeometry


  PURE SUBROUTINE ScaleFactors_1(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    REAL, INTENT(IN), DIMENSION(:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:) :: hx,hy,hz
    !------------------------------------------------------------------------!
    
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       CALL ScaleFactors_cartesian(hx,hy,hz)
    CASE(POLAR)
       CALL ScaleFactors_polar(coords(:,:,1),hx,hy,hz)
    CASE(CYLINDRICAL)
       CALL ScaleFactors_cylindrical(coords(:,:,2),hx,hy,hz)
    CASE(SPHERICAL)
       CALL ScaleFactors_spherical(coords(:,:,1),coords(:,:,2),hx,hy,hz)
    CASE(OBLATE_SPHEROIDAL)
       CALL ScaleFactors_oblatespheroidal(GetScale(this),coords(:,:,1), &
            coords(:,:,2),hx,hy,hz)
    END SELECT
  END SUBROUTINE ScaleFactors_1


  PURE SUBROUTINE ScaleFactors_2(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    REAL, INTENT(IN), DIMENSION(:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:) :: hx,hy,hz
    !------------------------------------------------------------------------!
    
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       CALL ScaleFactors_cartesian(hx,hy,hz)
    CASE(POLAR)
       CALL ScaleFactors_polar(coords(:,:,:,1),hx,hy,hz)
    CASE(CYLINDRICAL)
       CALL ScaleFactors_cylindrical(coords(:,:,:,2),hx,hy,hz)
    CASE(SPHERICAL)
       CALL ScaleFactors_spherical(coords(:,:,:,1),coords(:,:,:,2),hx,hy,hz)
    CASE(OBLATE_SPHEROIDAL)
       CALL ScaleFactors_oblatespheroidal(GetScale(this),coords(:,:,:,1), &
            coords(:,:,:,2),hx,hy,hz)
    END SELECT
  END SUBROUTINE ScaleFactors_2


  ! coordinate conversion
  PURE SUBROUTINE Convert2Cartesian_coords(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP)     :: this
    REAL, DIMENSION(:,:,:) :: curv, cart
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,curv
    INTENT(OUT)   :: cart
    !------------------------------------------------------------------------!

    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       cart(:,:,:) = curv(:,:,:)
    CASE(POLAR)
       CALL Convert2Cartesian_polar(curv(:,:,1),curv(:,:,2),cart(:,:,1),cart(:,:,2))
    CASE(CYLINDRICAL)
       CALL Convert2Cartesian_cylindrical(curv(:,:,1),curv(:,:,2),cart(:,:,1),cart(:,:,2))
    CASE(SPHERICAL)
       CALL Convert2Cartesian_spherical(curv(:,:,1),curv(:,:,2),cart(:,:,1),cart(:,:,2))
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Cartesian_oblatespheroidal(GetScale(this),curv(:,:,1),curv(:,:,2), &
            cart(:,:,1),cart(:,:,2))
    END SELECT
    
  END SUBROUTINE Convert2Cartesian_coords


  PURE SUBROUTINE Convert2Curvilinear_coords(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP)     :: this
    REAL, DIMENSION(:,:,:) :: cart, curv
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,cart
    INTENT(OUT)   :: curv
    !------------------------------------------------------------------------!

    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       curv(:,:,:) = cart(:,:,:)
    CASE(POLAR)
       CALL Convert2Curvilinear_polar(cart(:,:,1),cart(:,:,2),curv(:,:,1),curv(:,:,2))
    CASE(CYLINDRICAL)
       CALL Convert2Curvilinear_cylindrical(cart(:,:,1),cart(:,:,2),curv(:,:,1),curv(:,:,2))
    CASE(SPHERICAL)
       CALL Convert2Curvilinear_spherical(cart(:,:,1),cart(:,:,2),curv(:,:,1),curv(:,:,2))
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Curvilinear_spherical(GetScale(this),cart(:,:,1),cart(:,:,2), &
            curv(:,:,1),curv(:,:,2))
    END SELECT
    
  END SUBROUTINE Convert2Curvilinear_coords


  ! vector conversion
  PURE SUBROUTINE Convert2Cartesian_vectors(this,curv,v_curv,v_cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP)     :: this
    REAL, DIMENSION(:,:,:) :: curv, v_curv, v_cart
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,curv,v_curv
    INTENT(OUT)   :: v_cart
    !------------------------------------------------------------------------!

    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       v_cart(:,:,:) = v_curv(:,:,:)
    CASE(POLAR)
       CALL Convert2Cartesian_polar(curv(:,:,2),v_curv(:,:,1),v_curv(:,:,2), &
            v_cart(:,:,1),v_cart(:,:,2))
    CASE(CYLINDRICAL)
       CALL Convert2Cartesian_cylindrical(v_curv(:,:,1),v_curv(:,:,2), &
            v_cart(:,:,1),v_cart(:,:,2))
    CASE(SPHERICAL)
       CALL Convert2Cartesian_spherical(curv(:,:,2),v_curv(:,:,1),v_curv(:,:,2), &
            v_cart(:,:,1),v_cart(:,:,2))
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Cartesian_oblatespheroidal(GetScale(this),curv(:,:,1),curv(:,:,2), &
            v_curv(:,:,1),v_curv(:,:,2),v_cart(:,:,1),v_cart(:,:,2))
    END SELECT
    
  END SUBROUTINE Convert2Cartesian_vectors


  PURE SUBROUTINE Convert2Curvilinear_vectors(this,curv,v_cart,v_curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP)     :: this
    REAL, DIMENSION(:,:,:) :: curv, v_curv, v_cart
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,curv,v_cart
    INTENT(OUT)   :: v_curv
    !------------------------------------------------------------------------!

    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       v_curv(:,:,:) = v_cart(:,:,:)
    CASE(POLAR)
       CALL Convert2Curvilinear_polar(curv(:,:,2),v_cart(:,:,1),v_cart(:,:,2), &
            v_curv(:,:,1),v_curv(:,:,2))
    CASE(CYLINDRICAL)
       CALL Convert2Curvilinear_cylindrical(v_cart(:,:,1),v_cart(:,:,2), &
            v_curv(:,:,1),v_curv(:,:,2))
    CASE(SPHERICAL)
       CALL Convert2Curvilinear_spherical(curv(:,:,2),v_cart(:,:,1),v_cart(:,:,2), &
            v_curv(:,:,1),v_curv(:,:,2))
     CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Curvilinear_oblatespheroidal(GetScale(this),curv(:,:,1),curv(:,:,2), &
            v_cart(:,:,1),v_cart(:,:,2),v_curv(:,:,1),v_curv(:,:,2))
    END SELECT
    
  END SUBROUTINE Convert2Curvilinear_vectors

END MODULE geometry_generic
