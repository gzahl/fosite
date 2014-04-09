!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_generic.f90                                              #
!#                                                                           #
!# Copyright (C) 2007-2010                                                   #
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
! generic module for geometrical properties
!----------------------------------------------------------------------------!
MODULE geometry_generic
  USE geometry_cartesian, InitGeometry_common => InitGeometry
  USE geometry_polar
  USE geometry_logpolar
  USE geometry_cylindrical
  USE geometry_spherical
  USE geometry_oblatespheroidal
  USE geometry_tancylindrical
  USE geometry_tanpolar
  USE geometry_sinhpolar
  USE geometry_sinhspherical
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE ScaleFactors
     MODULE PROCEDURE ScaleFactors_1, ScaleFactors_2
  END INTERFACE
  INTERFACE Convert2Cartesian
     MODULE PROCEDURE Convert2Cartesian_coords_1, Convert2Cartesian_coords_2
     MODULE PROCEDURE Convert2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear
     MODULE PROCEDURE Convert2Curvilinear_coords_1, Convert2Curvilinear_coords_2
     MODULE PROCEDURE Convert2Curvilinear_vectors
  END INTERFACE
  !--------------------------------------------------------------------------!
  INTEGER, PARAMETER :: CARTESIAN         = 1
  INTEGER, PARAMETER :: POLAR             = 20
  INTEGER, PARAMETER :: LOGPOLAR          = 21
  INTEGER, PARAMETER :: TANPOLAR          = 22
  INTEGER, PARAMETER :: SINHPOLAR         = 23
  INTEGER, PARAMETER :: CYLINDRICAL       = 30
  INTEGER, PARAMETER :: TANCYLINDRICAL    = 31
  INTEGER, PARAMETER :: SPHERICAL         = 40
  INTEGER, PARAMETER :: SINHSPHERICAL     = 41
  INTEGER, PARAMETER :: OBLATE_SPHEROIDAL = 50
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Geometry_TYP, &
       ! constants
       CARTESIAN, POLAR, LOGPOLAR, TANPOLAR, SINHPOLAR, &
       CYLINDRICAL, TANCYLINDRICAL, SPHERICAL, SINHSPHERICAL, &
       OBLATE_SPHEROIDAL, &
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
    CHARACTER(LEN=8) :: gs_str
    REAL :: gs_def
    !------------------------------------------------------------------------!
    ! default scale
    IF (PRESENT(gs)) THEN
       gs_def = gs
    ELSE
       gs_def = 1.0
    END IF

    ! default: do not assume a spherical geometry (<e_eta,e_r> /= 0 )
    this%spherical_like = .FALSE.

    SELECT CASE(gt)
    CASE(CARTESIAN)
       CALL InitGeometry_cartesian(this,gt)
    CASE(POLAR)
       CALL InitGeometry_polar(this,gt)
    CASE(LOGPOLAR)
       CALL InitGeometry_logpolar(this,gt,gs_def)
    CASE(TANPOLAR)
       CALL InitGeometry_tanpolar(this,gt,gs_def)
    CASE(SINHPOLAR)
       CALL InitGeometry_sinhpolar(this,gt,gs_def)
    CASE(CYLINDRICAL)
       CALL InitGeometry_cylindrical(this,gt)
    CASE(TANCYLINDRICAL)
       CALL InitGeometry_tancyl(this,gt,gs_def)
    CASE(SPHERICAL)
       CALL InitGeometry_spherical(this,gt)
    CASE(SINHSPHERICAL)
       CALL InitGeometry_sinhspherical(this,gt,gs_def)
    CASE(OBLATE_SPHEROIDAL)
       CALL InitGeometry_oblatespher(this,gt,gs_def)
    CASE DEFAULT
       CALL Error(this,"InitGeometry",  "Unknown geometry.")
    END SELECT

    ! print some information
    CALL Info(this, " GEOMETRY-> coordinates:       " // TRIM(GetName(this)))
    SELECT CASE(gt)
    CASE(LOGPOLAR,TANPOLAR,SINHPOLAR,TANCYLINDRICAL,SINHSPHERICAL,OBLATE_SPHEROIDAL)
       WRITE (gs_str,'(ES8.1)') GetScale(this)
       CALL Info(this, "            geometry scale:    " // TRIM(gs_str))
    END SELECT
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
    CASE(LOGPOLAR)
       CALL ScaleFactors_logpolar(GetScale(this),coords(:,:,1),hx,hy,hz)
    CASE(TANPOLAR)
       CALL ScaleFactors_tanpolar(GetScale(this),coords(:,:,1),hx,hy,hz)
    CASE(SINHPOLAR)
       CALL ScaleFactors_sinhpolar(GetScale(this),coords(:,:,1),hx,hy,hz)
    CASE(CYLINDRICAL)
       CALL ScaleFactors_cylindrical(coords(:,:,2),hx,hy,hz)
    CASE(TANCYLINDRICAL)
       CALL ScaleFactors_tancyl(GetScale(this),coords(:,:,1), &
            coords(:,:,2),hx,hy,hz)
    CASE(SPHERICAL)
       CALL ScaleFactors_spherical(coords(:,:,1),coords(:,:,2),hx,hy,hz)
    CASE(SINHSPHERICAL)
       CALL ScaleFactors_sinhspherical(GetScale(this),coords(:,:,1),coords(:,:,2),&
            hx,hy,hz)
    CASE(OBLATE_SPHEROIDAL)
       CALL ScaleFactors_oblatespher(GetScale(this),coords(:,:,1), &
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
    CASE(LOGPOLAR)
       CALL ScaleFactors_logpolar(GetScale(this),coords(:,:,:,1),hx,hy,hz)
    CASE(TANPOLAR)
       CALL ScaleFactors_tanpolar(GetScale(this),coords(:,:,:,1),hx,hy,hz)
    CASE(SINHPOLAR)
       CALL ScaleFactors_sinhpolar(GetScale(this),coords(:,:,:,1),hx,hy,hz)
    CASE(CYLINDRICAL)
       CALL ScaleFactors_cylindrical(coords(:,:,:,2),hx,hy,hz)
    CASE(TANCYLINDRICAL)
       CALL ScaleFactors_tancyl(GetScale(this),coords(:,:,:,1), &
            coords(:,:,:,2),hx,hy,hz)
    CASE(SPHERICAL)
       CALL ScaleFactors_spherical(coords(:,:,:,1),coords(:,:,:,2),hx,hy,hz)
    CASE(SINHSPHERICAL)
       CALL ScaleFactors_sinhspherical(GetScale(this),coords(:,:,:,1),coords(:,:,:,2),&
            hx,hy,hz)
    CASE(OBLATE_SPHEROIDAL)
       CALL ScaleFactors_oblatespher(GetScale(this),coords(:,:,:,1), &
            coords(:,:,:,2),hx,hy,hz)
    END SELECT
  END SUBROUTINE ScaleFactors_2


  ! coordinate conversion
  PURE SUBROUTINE Convert2Cartesian_coords_1(this,curv,cart)
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
    CASE(LOGPOLAR)
       CALL Convert2Cartesian_logpolar(GetScale(this),curv(:,:,1),curv(:,:,2),&
            cart(:,:,1),cart(:,:,2))
    CASE(TANPOLAR)
       CALL Convert2Cartesian_tanpolar(GetScale(this),curv(:,:,1),curv(:,:,2),&
            cart(:,:,1),cart(:,:,2))
    CASE(SINHPOLAR)
       CALL Convert2Cartesian_sinhpolar(GetScale(this),curv(:,:,1),curv(:,:,2),&
            cart(:,:,1),cart(:,:,2))
    CASE(CYLINDRICAL)
       CALL Convert2Cartesian_cylindrical(curv(:,:,1),curv(:,:,2),cart(:,:,1),&
            cart(:,:,2))
    CASE(TANCYLINDRICAL)
       CALL Convert2Cartesian_tancyl(GetScale(this),curv(:,:,1),curv(:,:,2),&
            cart(:,:,1),cart(:,:,2))
    CASE(SPHERICAL)
       CALL Convert2Cartesian_spherical(curv(:,:,1),curv(:,:,2),cart(:,:,1),&
            cart(:,:,2))
    CASE(SINHSPHERICAL)
       CALL Convert2Cartesian_sinhspherical(GetScale(this),curv(:,:,1),curv(:,:,2),&
            cart(:,:,1),cart(:,:,2))
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Cartesian_oblatespher(GetScale(this),curv(:,:,1),curv(:,:,2),&
            cart(:,:,1),cart(:,:,2))
    END SELECT
    
  END SUBROUTINE Convert2Cartesian_coords_1


  PURE SUBROUTINE Convert2Cartesian_coords_2(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP)     :: this
    REAL, DIMENSION(:,:,:,:) :: curv, cart
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,curv
    INTENT(OUT)   :: cart
    !------------------------------------------------------------------------!

    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       cart(:,:,:,:) = curv(:,:,:,:)
    CASE(POLAR)
       CALL Convert2Cartesian_polar(curv(:,:,:,1),curv(:,:,:,2),cart(:,:,:,1),cart(:,:,:,2))
    CASE(LOGPOLAR)
       CALL Convert2Cartesian_logpolar(GetScale(this),curv(:,:,:,1),curv(:,:,:,2),&
            cart(:,:,:,1),cart(:,:,:,2))
    CASE(TANPOLAR)
       CALL Convert2Cartesian_tanpolar(GetScale(this),curv(:,:,:,1),curv(:,:,:,2),&
            cart(:,:,:,1),cart(:,:,:,2))
    CASE(SINHPOLAR)
       CALL Convert2Cartesian_sinhpolar(GetScale(this),curv(:,:,:,1),curv(:,:,:,2),&
            cart(:,:,:,1),cart(:,:,:,2))
    CASE(CYLINDRICAL)
       CALL Convert2Cartesian_cylindrical(curv(:,:,:,1),curv(:,:,:,2),cart(:,:,:,1),&
            cart(:,:,:,2))
    CASE(TANCYLINDRICAL)
       CALL Convert2Cartesian_tancyl(GetScale(this),curv(:,:,:,1),curv(:,:,:,2),&
            cart(:,:,:,1),cart(:,:,:,2))
    CASE(SPHERICAL)
       CALL Convert2Cartesian_spherical(curv(:,:,:,1),curv(:,:,:,2),cart(:,:,:,1),&
            cart(:,:,:,2))
    CASE(SINHSPHERICAL)
       CALL Convert2Cartesian_sinhspherical(GetScale(this),curv(:,:,:,1),curv(:,:,:,2), &
            cart(:,:,:,1),cart(:,:,:,2))
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Cartesian_oblatespher(GetScale(this),curv(:,:,:,1),curv(:,:,:,2),&
            cart(:,:,:,1),cart(:,:,:,2))
    END SELECT
  
  END SUBROUTINE Convert2Cartesian_coords_2


  PURE SUBROUTINE Convert2Curvilinear_coords_1(this,cart,curv)
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
    CASE(LOGPOLAR)
       CALL Convert2Curvilinear_logpolar(GetScale(this),cart(:,:,1),cart(:,:,2),&
            curv(:,:,1),curv(:,:,2))
    CASE(TANPOLAR)
       CALL Convert2Curvilinear_tanpolar(GetScale(this),cart(:,:,1),cart(:,:,2),&
            curv(:,:,1),curv(:,:,2))
    CASE(SINHPOLAR)
       CALL Convert2Curvilinear_sinhpolar(GetScale(this),cart(:,:,1),cart(:,:,2),&
            curv(:,:,1),curv(:,:,2))
    CASE(CYLINDRICAL)
       CALL Convert2Curvilinear_cylindrical(cart(:,:,1),cart(:,:,2),curv(:,:,1),curv(:,:,2))
    CASE(TANCYLINDRICAL)
       CALL Convert2Curvilinear_tancyl(GetScale(this),cart(:,:,1),cart(:,:,2), &
            curv(:,:,1),curv(:,:,2))
    CASE(SPHERICAL)
       CALL Convert2Curvilinear_spherical(cart(:,:,1),cart(:,:,2),curv(:,:,1),curv(:,:,2))
    CASE(SINHSPHERICAL)
       CALL Convert2Curvilinear_sinhspherical(GetScale(this),cart(:,:,1),cart(:,:,2), &
            curv(:,:,1),curv(:,:,2))
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Curvilinear_oblatespher(GetScale(this),cart(:,:,1),cart(:,:,2), &
            curv(:,:,1),curv(:,:,2))
    END SELECT
    
  END SUBROUTINE Convert2Curvilinear_coords_1


  PURE SUBROUTINE Convert2Curvilinear_coords_2(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP)     :: this
    REAL, DIMENSION(:,:,:,:) :: cart, curv
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,cart
    INTENT(OUT)   :: curv
    !------------------------------------------------------------------------!

    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       curv(:,:,:,:) = cart(:,:,:,:)
    CASE(POLAR)
       CALL Convert2Curvilinear_polar(cart(:,:,:,1),cart(:,:,:,2),curv(:,:,:,1),curv(:,:,:,2))
    CASE(LOGPOLAR)
       CALL Convert2Curvilinear_logpolar(GetScale(this),cart(:,:,:,1),cart(:,:,:,2),&
            curv(:,:,:,1),curv(:,:,:,2))
    CASE(TANPOLAR)
       CALL Convert2Curvilinear_tanpolar(GetScale(this),cart(:,:,:,1),cart(:,:,:,2),&
            curv(:,:,:,1),curv(:,:,:,2))
    CASE(SINHPOLAR)
       CALL Convert2Curvilinear_sinhpolar(GetScale(this),cart(:,:,:,1),cart(:,:,:,2),&
            curv(:,:,:,1),curv(:,:,:,2))
    CASE(CYLINDRICAL)
       CALL Convert2Curvilinear_cylindrical(cart(:,:,:,1),cart(:,:,:,2),curv(:,:,:,1),curv(:,:,:,2))
    CASE(TANCYLINDRICAL)
       CALL Convert2Curvilinear_tancyl(GetScale(this),cart(:,:,:,1),cart(:,:,:,2), &
            curv(:,:,:,1),curv(:,:,:,2))
    CASE(SPHERICAL)
       CALL Convert2Curvilinear_spherical(cart(:,:,:,1),cart(:,:,:,2),curv(:,:,:,1),curv(:,:,:,2))
    CASE(SINHSPHERICAL)
       CALL Convert2Curvilinear_sinhspherical(GetScale(this),cart(:,:,:,1),cart(:,:,:,2), &
            curv(:,:,:,1),curv(:,:,:,2))
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Curvilinear_oblatespher(GetScale(this),cart(:,:,:,1),cart(:,:,:,2), &
            curv(:,:,:,1),curv(:,:,:,2))
    END SELECT
  
  END SUBROUTINE Convert2Curvilinear_coords_2


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
    CASE(LOGPOLAR)
       CALL Convert2Cartesian_logpolar(GetScale(this),curv(:,:,2),v_curv(:,:,1),&
            v_curv(:,:,2),v_cart(:,:,1),v_cart(:,:,2))
    CASE(TANPOLAR)
       CALL Convert2Cartesian_tanpolar(GetScale(this),curv(:,:,2),v_curv(:,:,1),&
            v_curv(:,:,2),v_cart(:,:,1),v_cart(:,:,2))
    CASE(SINHPOLAR)
       CALL Convert2Cartesian_sinhpolar(GetScale(this),curv(:,:,2),v_curv(:,:,1),&
            v_curv(:,:,2),v_cart(:,:,1),v_cart(:,:,2))
    CASE(CYLINDRICAL)
       CALL Convert2Cartesian_cylindrical(v_curv(:,:,1),v_curv(:,:,2), &
            v_cart(:,:,1),v_cart(:,:,2))
    CASE(TANCYLINDRICAL)
       CALL Convert2Cartesian_tancyl(v_curv(:,:,1),v_curv(:,:,2), &
            v_cart(:,:,1),v_cart(:,:,2))
    CASE(SPHERICAL)
       CALL Convert2Cartesian_spherical(GetScale(this),curv(:,:,2),v_curv(:,:,1), &
            v_curv(:,:,2),v_cart(:,:,1),v_cart(:,:,2))
    CASE(SINHSPHERICAL)
       CALL Convert2Cartesian_sinhspherical(GetScale(this),curv(:,:,2),v_curv(:,:,1), &
            v_curv(:,:,2),v_cart(:,:,1),v_cart(:,:,2))
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Cartesian_oblatespher(GetScale(this),curv(:,:,1),curv(:,:,2), &
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
    CASE(LOGPOLAR)
       CALL Convert2Curvilinear_logpolar(GetScale(this),curv(:,:,2),v_cart(:,:,1),&
            v_cart(:,:,2),v_curv(:,:,1),v_curv(:,:,2))
    CASE(TANPOLAR)
       CALL Convert2Curvilinear_tanpolar(GetScale(this),curv(:,:,2),v_cart(:,:,1),&
            v_cart(:,:,2),v_curv(:,:,1),v_curv(:,:,2))
    CASE(SINHPOLAR)
       CALL Convert2Curvilinear_sinhpolar(GetScale(this),curv(:,:,2),v_cart(:,:,1),&
            v_cart(:,:,2),v_curv(:,:,1),v_curv(:,:,2))
    CASE(CYLINDRICAL)
       CALL Convert2Curvilinear_cylindrical(v_cart(:,:,1),v_cart(:,:,2), &
            v_curv(:,:,1),v_curv(:,:,2))
    CASE(TANCYLINDRICAL)
       CALL Convert2Curvilinear_tancyl(v_cart(:,:,1),v_cart(:,:,2), &
            v_curv(:,:,1),v_curv(:,:,2))
    CASE(SPHERICAL)
       CALL Convert2Curvilinear_spherical(GetScale(this),curv(:,:,2),v_cart(:,:,1), &
            v_cart(:,:,2),v_curv(:,:,1),v_curv(:,:,2))
    CASE(SINHSPHERICAL)
       CALL Convert2Curvilinear_sinhspherical(GetScale(this),curv(:,:,2),v_cart(:,:,1), &
            v_cart(:,:,2),v_curv(:,:,1),v_curv(:,:,2))
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Curvilinear_oblatespher(curv(:,:,1),curv(:,:,2), &
            v_cart(:,:,1),v_cart(:,:,2),v_curv(:,:,1),v_curv(:,:,2))
    END SELECT
    
  END SUBROUTINE Convert2Curvilinear_vectors

END MODULE geometry_generic
