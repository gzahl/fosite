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
!> \author Tobias Illenseer
!!
!! \brief generic module for geometrical properties
!!
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_generic
  USE geometry_cartesian, InitGeometry_common => InitGeometry, &
       CloseGeometry_common => CloseGeometry
  USE geometry_sinhcartesian
  USE geometry_polar
  USE geometry_logpolar
  USE geometry_polypolar
  USE geometry_cylindrical
  USE geometry_spherical
  USE geometry_oblatespheroidal
  USE geometry_tancylindrical
  USE geometry_tanpolar
  USE geometry_sinhpolar
  USE geometry_sinhtanh
  USE geometry_bianglespherical
  USE geometry_sinhspherical
  USE geometry_channel
  USE geometry_lncoshcylindrical
  USE geometry_elliptic
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE ScaleFactors
     MODULE PROCEDURE ScaleFactors_1, ScaleFactors_2
  END INTERFACE
  INTERFACE Convert2Cartesian
     MODULE PROCEDURE Convert2Cartesian_point
     MODULE PROCEDURE Convert2Cartesian_coords_1, Convert2Cartesian_coords_2
     MODULE PROCEDURE Convert2Cartesian_vectors_1, Convert2Cartesian_vectors_2,&
                      Convert2Cartesian_vectors_3
  END INTERFACE
  INTERFACE Convert2Curvilinear
     MODULE PROCEDURE Convert2Curvilinear_point
     MODULE PROCEDURE Convert2Curvilinear_coords_1, Convert2Curvilinear_coords_2
     MODULE PROCEDURE Convert2Curvilinear_vectors_1, Convert2Curvilinear_vectors_2,&
                      Convert2Curvilinear_vectors_3
  END INTERFACE
  INTERFACE Radius
     MODULE PROCEDURE Radius_1, Radius_2
  END INTERFACE
  INTERFACE PositionVector
     MODULE PROCEDURE PositionVector_1, PositionVector_2
  END INTERFACE
  !> \endcond
  !> \name Public Attributes
  !! #### geometries
  !--------------------------------------------------------------------------!
  INTEGER, PARAMETER :: CARTESIAN         = 1
  INTEGER, PARAMETER :: SINHCARTESIAN     = 2
  INTEGER, PARAMETER :: POLAR             = 20
  INTEGER, PARAMETER :: LOGPOLAR          = 21
  INTEGER, PARAMETER :: TANPOLAR          = 22
  INTEGER, PARAMETER :: SINHPOLAR         = 23
  INTEGER, PARAMETER :: SINHTANHPOLAR     = 24
  INTEGER, PARAMETER :: POLYPOLAR         = 25
  INTEGER, PARAMETER :: ELLIPTIC          = 27
  INTEGER, PARAMETER :: CYLINDRICAL       = 30
  INTEGER, PARAMETER :: TANCYLINDRICAL    = 31
  INTEGER, PARAMETER :: LNCOSHCYLINDRICAL = 32
  INTEGER, PARAMETER :: SPHERICAL         = 40
  INTEGER, PARAMETER :: SINHSPHERICAL     = 41
  INTEGER, PARAMETER :: BIANGLESPHERICAL  = 42
  INTEGER, PARAMETER :: OBLATE_SPHEROIDAL = 50
  INTEGER, PARAMETER :: CHANNEL           = 60
  !> \}
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Geometry_TYP, &
       ! constants
       PI, &
       CARTESIAN, POLAR, LOGPOLAR, TANPOLAR, SINHPOLAR, SINHTANHPOLAR, &
       CYLINDRICAL, TANCYLINDRICAL, LNCOSHCYLINDRICAL, SPHERICAL, SINHSPHERICAL, &
       BIANGLESPHERICAL, OBLATE_SPHEROIDAL, CHANNEL, POLYPOLAR, ELLIPTIC, &
       SINHCARTESIAN, &
       ! methods
       InitGeometry, &
       CloseGeometry, &
       Convert2Cartesian, &
       Convert2Curvilinear, &
       ScaleFactors, &
       Radius, &
       PositionVector, &
       GetScale, &
       SetScale, &
       GetType, &
       GetName, &
       GetRank, &
       GetNumProcs, &
       Initialized, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor of generic geometry module
  SUBROUTINE InitGeometry(this,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP)      :: this
    TYPE(Dict_TYP),POINTER  :: config
    INTEGER                 :: gt
    REAL                    :: gs_def,gs_def2,gs_def3, dr
    !------------------------------------------------------------------------!
    CHARACTER(LEN=8)        :: gs_str
    !------------------------------------------------------------------------!
    INTENT(INOUT)           :: this
    !------------------------------------------------------------------------!
    CALL RequireKey(config, "geometry")
    CALL GetAttr(config, "geometry", gt)
    
    ! check if geometry parameters were given
    ! and set to defaults if not
    SELECT CASE(gt)
    CASE(CARTESIAN,POLAR,CYLINDRICAL,SPHERICAL)
       ! do nothing (no parameters needed)
    CASE DEFAULT
       ! geometries with at least one parameter
       ! gparam defaults to 1 except for CHANNEL
       IF (gt.EQ.CHANNEL) THEN
          CALL RequireKey(config, "gparam", 0.5)
       ELSE
          CALL RequireKey(config, "gparam", 1.0)
       END IF
       CALL GetAttr(config, "gparam", gs_def)
       ! geometries with three parameters
       SELECT CASE(gt)
       CASE(SINHPOLAR)
          ! gp2,gp3 defaults to 0
          CALL RequireKey(config, "gparam2", 0.0)
          CALL GetAttr(config, "gparam2", gs_def2)
          CALL RequireKey(config, "gparam3", 0.0)
          CALL GetAttr(config, "gparam3", gs_def3)
       CASE(SINHTANHPOLAR,LNCOSHCYLINDRICAL)
          ! gp2,gp3 defaults to 1
          CALL RequireKey(config, "gparam2", 1.0)
          CALL GetAttr(config, "gparam2", gs_def2)
          CALL RequireKey(config, "gparam3", 1.0)
          CALL GetAttr(config, "gparam3", gs_def3)
       END SELECT
    END SELECT


    SELECT CASE(gt)
    CASE(CARTESIAN)
       CALL InitGeometry_cartesian(this,gt)
       CALL RequireKey(config, "dz", 1.0)
    CASE(SINHCARTESIAN)
       CALL InitGeometry_sinhcartesian(this,gt,gs_def)
       CALL RequireKey(config, "dz", 1.0)
    CASE(POLAR)
       CALL InitGeometry_polar(this,gt)
       CALL RequireKey(config, "dz", 1.0)
    CASE(LOGPOLAR)
       CALL InitGeometry_logpolar(this,gt,gs_def)
       CALL RequireKey(config, "dz", 1.0)
    CASE(TANPOLAR)
       CALL InitGeometry_tanpolar(this,gt,gs_def)
       CALL RequireKey(config, "dz", 1.0)
    CASE(SINHPOLAR)
       CALL InitGeometry_sinhpolar(this,gt,gs_def,gs_def2,gs_def3)
       CALL RequireKey(config, "dz", 1.0)
    CASE(SINHTANHPOLAR)
       CALL InitGeometry_sinhtanh(this,gt,gs_def,gs_def2,gs_def3)
       CALL RequireKey(config, "dz", 1.0)
    CASE(POLYPOLAR)
       CALL InitGeometry_polypolar(this,gt,gs_def)
       CALL RequireKey(config, "dz", 1.0)
    CASE(ELLIPTIC)
       CALL InitGeometry_elliptic(this,gt,gs_def)
       CALL RequireKey(config, "dz", 1.0)
    CASE(CYLINDRICAL)
       CALL InitGeometry_cylindrical(this,gt)
       CALL RequireKey(config, "dz", 2.0*PI)
    CASE(TANCYLINDRICAL)
       CALL InitGeometry_tancyl(this,gt,gs_def)
       CALL RequireKey(config, "dz", 2.0*PI)
    CASE(LNCOSHCYLINDRICAL)
       CALL InitGeometry_lncoshcyl(this,gt,gs_def,gs_def2,gs_def3)
       CALL RequireKey(config, "dz", 2.0*PI)
    CASE(SPHERICAL)
       CALL InitGeometry_spherical(this,gt)
       CALL RequireKey(config, "dz", 2.0*PI)
    CASE(SINHSPHERICAL)
       CALL InitGeometry_sinhspher(this,gt,gs_def)
       CALL RequireKey(config, "dz", 2.0*PI)
    CASE(BIANGLESPHERICAL)
       CALL InitGeometry_bianglespher(this,gt,gs_def)
       CALL RequireKey(config, "dz", 1.0)
    CASE(OBLATE_SPHEROIDAL)
       CALL InitGeometry_oblatespher(this,gt,gs_def)
       CALL RequireKey(config, "dz", 2.0*PI)
    CASE(CHANNEL)
       CALL InitGeometry_channel(this,gt,gs_def)
       CALL RequireKey(config, "dz", 1.0)
    CASE DEFAULT
       CALL Error(this,"InitGeometry",  "Unknown geometry.")
    END SELECT

    ! print some information
    CALL Info(this, " GEOMETRY-> coordinates:       " // TRIM(GetName(this)))
    SELECT CASE(gt)
    CASE(LOGPOLAR,TANPOLAR,SINHPOLAR,TANCYLINDRICAL,SINHSPHERICAL,&
         BIANGLESPHERICAL,OBLATE_SPHEROIDAL,CHANNEL,POLYPOLAR,ELLIPTIC,&
         SINHCARTESIAN)
       WRITE (gs_str,'(ES8.1)') GetScale(this)
       CALL Info(this, "            geometry scale:    " // TRIM(gs_str))
    CASE(LNCOSHCYLINDRICAL,SINHTANHPOLAR)
       WRITE (gs_str,'(ES8.1)') GetScale(this)
       CALL Info(this, "            geometry scale:    " // TRIM(gs_str))
       WRITE (gs_str,'(ES8.1)') GetScale(this,2)
       CALL Info(this, "            geometry scale:    " // TRIM(gs_str))
       WRITE (gs_str,'(ES8.1)') GetScale(this,3)
       CALL Info(this, "            geometry scale:    " // TRIM(gs_str))
    
    END SELECT
  END SUBROUTINE InitGeometry


  !> \public Compute scale factors
  !!
  !! Computes the scale factors of the given geometry.
  !! \f[
  !!   h_x = \left|\frac{\partial \vec{r}}{\partial x}\right|,\quad
  !!   h_y = \left|\frac{\partial \vec{r}}{\partial y}\right|,\quad
  !!   h_z = \left|\frac{\partial \vec{r}}{\partial z}\right|
  !! \f]
  PURE SUBROUTINE ScaleFactors_1(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    REAL, INTENT(IN), DIMENSION(:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:) :: hx,hy,hz
    !------------------------------------------------------------------------!
!CDIR IEXPAND    
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       CALL ScaleFactors_cartesian(hx,hy,hz)
    CASE(SINHCARTESIAN)
       CALL ScaleFactors_sinhcartesian(GetScale(this),&
          coords(:,:,1),coords(:,:,2),hx,hy,hz)
    CASE(POLAR)
       CALL ScaleFactors_polar(coords(:,:,1),hx,hy,hz)
    CASE(LOGPOLAR)
       CALL ScaleFactors_logpolar(GetScale(this),coords(:,:,1),hx,hy,hz)
    CASE(TANPOLAR)
       CALL ScaleFactors_tanpolar(GetScale(this),coords(:,:,1),hx,hy,hz)
    CASE(SINHPOLAR)
       CALL ScaleFactors_sinhpolar(GetScale(this),GetScale(this,2),&
            GetScale(this,3),coords(:,:,1),hx,hy,hz)
    CASE(SINHTANHPOLAR)
       CALL ScaleFactors_sinhtanh(GetScale(this),GetScale(this,2),&
                                  GetScale(this,3),coords(:,:,1),&
                                  coords(:,:,2),hx,hy,hz)
    CASE(POLYPOLAR)
       CALL ScaleFactors_polypolar(GetScale(this),coords(:,:,1),hx,hy,hz)
    CASE(ELLIPTIC)
       CALL ScaleFactors_elliptic(GetScale(this),coords(:,:,1), &
            coords(:,:,2),hx,hy,hz)
    CASE(CYLINDRICAL)
       CALL ScaleFactors_cylindrical(coords(:,:,2),hx,hy,hz)
    CASE(TANCYLINDRICAL)
       CALL ScaleFactors_tancyl(GetScale(this),coords(:,:,1), &
            coords(:,:,2),hx,hy,hz)
    CASE(LNCOSHCYLINDRICAL)
       CALL ScaleFactors_lncoshcyl(GetScale(this),GetScale(this,2),GetScale(this,3),&
            coords(:,:,2),hx,hy,hz)
    CASE(SPHERICAL)
       CALL ScaleFactors_spherical(coords(:,:,1),coords(:,:,2),hx,hy,hz)
    CASE(SINHSPHERICAL)
       CALL ScaleFactors_sinhspher(GetScale(this),coords(:,:,1),coords(:,:,2),&
            hx,hy,hz)
    CASE(BIANGLESPHERICAL)
       CALL ScaleFactors_bianglespher(GetScale(this),coords(:,:,1),coords(:,:,2), &
            hx,hy,hz)
    CASE(OBLATE_SPHEROIDAL)
       CALL ScaleFactors_oblatespher(GetScale(this),coords(:,:,1), &
            coords(:,:,2),hx,hy,hz)
    CASE(CHANNEL)
       CALL ScaleFactors_channel(GetScale(this),coords(:,:,2),hx,hy,hz)
    END SELECT
  END SUBROUTINE ScaleFactors_1


  !> \public Compute scale factors
  !!
  !! Computes the scale factors of the given geometry.
  !! \f[
  !!   h_x = \left|\frac{\partial \vec{r}}{\partial x}\right|,\quad
  !!   h_y = \left|\frac{\partial \vec{r}}{\partial y}\right|,\quad
  !!   h_z = \left|\frac{\partial \vec{r}}{\partial z}\right|
  !! \f]
  PURE SUBROUTINE ScaleFactors_2(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    REAL, INTENT(IN), DIMENSION(:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:) :: hx,hy,hz
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       CALL ScaleFactors_cartesian(hx,hy,hz)
    CASE(SINHCARTESIAN)
       CALL ScaleFactors_sinhcartesian(GetScale(this),&
            coords(:,:,:,1),coords(:,:,:,2),hx,hy,hz)
    CASE(POLAR)
       CALL ScaleFactors_polar(coords(:,:,:,1),hx,hy,hz)
    CASE(LOGPOLAR)
       CALL ScaleFactors_logpolar(GetScale(this),coords(:,:,:,1),hx,hy,hz)
    CASE(TANPOLAR)
       CALL ScaleFactors_tanpolar(GetScale(this),coords(:,:,:,1),hx,hy,hz)
    CASE(SINHPOLAR)
       CALL ScaleFactors_sinhpolar(GetScale(this),GetScale(this,2),&
            GetScale(this,3),coords(:,:,:,1),hx,hy,hz)
    CASE(SINHTANHPOLAR)
       CALL ScaleFactors_sinhtanh(GetScale(this),GetScale(this,2),&
                                  GetScale(this,3),coords(:,:,:,1),&
                                  coords(:,:,:,2),hx,hy,hz)
    CASE(POLYPOLAR)
       CALL ScaleFactors_polypolar(GetScale(this),coords(:,:,:,1),hx,hy,hz)
    CASE(ELLIPTIC)
       CALL ScaleFactors_elliptic(GetScale(this),coords(:,:,:,1), &
            coords(:,:,:,2),hx,hy,hz)
    CASE(CYLINDRICAL)
       CALL ScaleFactors_cylindrical(coords(:,:,:,2),hx,hy,hz)
    CASE(TANCYLINDRICAL)
       CALL ScaleFactors_tancyl(GetScale(this),coords(:,:,:,1), &
            coords(:,:,:,2),hx,hy,hz)
    CASE(LNCOSHCYLINDRICAL)
       CALL ScaleFactors_lncoshcyl(GetScale(this),GetScale(this,2),GetScale(this,3),&
            coords(:,:,:,2),hx,hy,hz)
    CASE(SPHERICAL)
       CALL ScaleFactors_spherical(coords(:,:,:,1),coords(:,:,:,2),hx,hy,hz)
    CASE(SINHSPHERICAL)
       CALL ScaleFactors_sinhspher(GetScale(this),coords(:,:,:,1),coords(:,:,:,2),&
            hx,hy,hz)
	  CASE(BIANGLESPHERICAL)
       CALL ScaleFactors_bianglespher(GetScale(this),coords(:,:,:,1),coords(:,:,:,2),&
            hx,hy,hz)
    CASE(OBLATE_SPHEROIDAL)
       CALL ScaleFactors_oblatespher(GetScale(this),coords(:,:,:,1), &
            coords(:,:,:,2),hx,hy,hz)
    CASE(CHANNEL)
       CALL ScaleFactors_channel(GetScale(this),coords(:,:,:,2),hx,hy,hz)
    END SELECT
  END SUBROUTINE ScaleFactors_2


  !> \public Compute radial distances to the origin
  !!
  !! Computes the radial distances to the origin for the given
  !! coordinates depending on the geometry.
  PURE SUBROUTINE Radius_1(this,coords,radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    REAL, INTENT(IN), DIMENSION(:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:) :: radius
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
      radius = Radius_cartesian(coords(:,:,1),coords(:,:,2))
    CASE(SINHCARTESIAN)
      radius = Radius_sinhcartesian(GetScale(this),coords(:,:,1),coords(:,:,2))
    CASE(POLAR)
      radius = Radius_polar(coords(:,:,1))
    CASE(LOGPOLAR)
      radius = Radius_logpolar(GetScale(this),coords(:,:,1))
    CASE(TANPOLAR)
      radius = Radius_tanpolar(GetScale(this),coords(:,:,1))
    CASE(SINHPOLAR)
      radius = Radius_sinhpolar(GetScale(this),GetScale(this,2),GetScale(this,3), &
                               coords(:,:,1))
    CASE(SINHTANHPOLAR)
      radius = Radius_sinhtanh(GetScale(this),GetScale(this,2),GetScale(this,3), &
                               coords(:,:,1))
    CASE(POLYPOLAR)
      radius = Radius_polypolar(GetScale(this),coords(:,:,1))
    CASE(ELLIPTIC)
      radius = Radius_elliptic(GetScale(this),coords(:,:,1),coords(:,:,2))
    CASE(CYLINDRICAL)
      radius = Radius_cylindrical(coords(:,:,1),coords(:,:,2))
    CASE(TANCYLINDRICAL)
      radius = Radius_tancyl(GetScale(this),coords(:,:,1),coords(:,:,2))
    CASE(LNCOSHCYLINDRICAL)
      radius = Radius_lncoshcyl(GetScale(this),GetScale(this,2),GetScale(this,3), &
                                coords(:,:,1))
    CASE(SPHERICAL)
      radius = Radius_spherical(coords(:,:,1))
    CASE(SINHSPHERICAL)
      radius = Radius_sinhspher(GetScale(this),coords(:,:,1))
    CASE(BIANGLESPHERICAL)
      radius = Radius_bianglespher(coords(:,:,1))
    CASE(OBLATE_SPHEROIDAL)
      radius = Radius_oblatespher(GetScale(this),coords(:,:,1),coords(:,:,2))
    CASE(CHANNEL)
      radius = Radius_channel(GetScale(this),coords(:,:,1),coords(:,:,2))
    END SELECT
  END SUBROUTINE Radius_1


  !> \public Compute radial distances to the origin
  !!
  !! Computes the radial distances to the origin for the given
  !! coordinates depending on the geometry.
  SUBROUTINE Radius_2(this,coords,radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    REAL, INTENT(IN), DIMENSION(:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:) :: radius
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
      radius = Radius_cartesian(coords(:,:,:,1),coords(:,:,:,2))
     CASE(SINHCARTESIAN)
      radius = Radius_sinhcartesian(GetScale(this),coords(:,:,:,1),coords(:,:,:,2))
   CASE(POLAR)
      radius = Radius_polar(coords(:,:,:,1))
    CASE(LOGPOLAR)
      radius = Radius_logpolar(GetScale(this),coords(:,:,:,1))
    CASE(TANPOLAR)
      radius = Radius_tanpolar(GetScale(this),coords(:,:,:,1))
    CASE(SINHPOLAR)
      radius = Radius_sinhpolar(GetScale(this),GetScale(this,2),GetScale(this,3), &
                               coords(:,:,:,1))
    CASE(SINHTANHPOLAR)
      radius = Radius_sinhtanh(GetScale(this),GetScale(this,2),GetScale(this,3), &
                               coords(:,:,:,1))
    CASE(POLYPOLAR)
      radius = Radius_polypolar(GetScale(this),coords(:,:,:,1))
    CASE(ELLIPTIC)
      radius = Radius_elliptic(GetScale(this),coords(:,:,:,1),coords(:,:,:,2))
    CASE(CYLINDRICAL)
      radius = Radius_cylindrical(coords(:,:,:,1),coords(:,:,:,2))
    CASE(TANCYLINDRICAL)
      radius = Radius_tancyl(GetScale(this),coords(:,:,:,1),coords(:,:,:,2))
    CASE(LNCOSHCYLINDRICAL)
      radius = Radius_lncoshcyl(GetScale(this),GetScale(this,2),GetScale(this,3), &
                                coords(:,:,:,1))
    CASE(SPHERICAL)
      radius = Radius_spherical(coords(:,:,:,1))
    CASE(SINHSPHERICAL)
      radius = Radius_sinhspher(GetScale(this),coords(:,:,:,1))
    CASE(BIANGLESPHERICAL)
      radius = Radius_bianglespher(coords(:,:,:,1))
    CASE(OBLATE_SPHEROIDAL)
      radius = Radius_oblatespher(GetScale(this),coords(:,:,:,1),coords(:,:,:,2))
    CASE(CHANNEL)
      radius = Radius_channel(GetScale(this),coords(:,:,:,1),coords(:,:,:,2))
    END SELECT
  END SUBROUTINE Radius_2


  !> \public Compute position vector components
  !!
  !! Computes the curvilinear position vector components with respect to
  !! the given geometry.
  PURE SUBROUTINE PositionVector_1(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    REAL, INTENT(IN), DIMENSION(:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:) :: posvec
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
      posvec(:,:,:) = coords(:,:,:)
    CASE(SINHCARTESIAN)
      CALL PositionVector_sinhcartesian(GetScale(this),coords(:,:,1),coords(:,:,2), &
                                        posvec(:,:,1),posvec(:,:,2)) 
    CASE(POLAR)
      CALL PositionVector_polar(coords(:,:,1),posvec(:,:,1),posvec(:,:,2)) 
    CASE(LOGPOLAR)
      CALL PositionVector_logpolar(GetScale(this),coords(:,:,1),posvec(:,:,1),posvec(:,:,2)) 
    CASE(TANPOLAR)
      CALL PositionVector_tanpolar(GetScale(this),coords(:,:,1),posvec(:,:,1),posvec(:,:,2)) 
    CASE(SINHPOLAR)
      CALL PositionVector_sinhpolar(GetScale(this),GetScale(this,2),GetScale(this,3), &
                                   coords(:,:,1),posvec(:,:,1),posvec(:,:,2)) 
    CASE(SINHTANHPOLAR)
      CALL PositionVector_sinhtanh(GetScale(this),GetScale(this,2),GetScale(this,3), &
                                   coords(:,:,1),posvec(:,:,1),posvec(:,:,2)) 
    CASE(POLYPOLAR)
      CALL PositionVector_polypolar(GetScale(this),coords(:,:,1),posvec(:,:,1),posvec(:,:,2)) 
    CASE(ELLIPTIC)
      CALL PositionVector_elliptic(GetScale(this),coords(:,:,1),coords(:,:,2), &
                                   posvec(:,:,1),posvec(:,:,2)) 
    CASE(CYLINDRICAL)
      posvec(:,:,:) = coords(:,:,:)
    CASE(TANCYLINDRICAL)
      CALL PositionVector_tancyl(GetScale(this),coords(:,:,1),coords(:,:,2), &
                                 posvec(:,:,1),posvec(:,:,2))
    CASE(LNCOSHCYLINDRICAL)
      CALL PositionVector_lncoshcyl(GetScale(this),GetScale(this,2),GetScale(this,3), &
                               coords(:,:,1),coords(:,:,2),posvec(:,:,1),posvec(:,:,2)) 
    CASE(SPHERICAL)
      CALL PositionVector_spherical(coords(:,:,1),posvec(:,:,1),posvec(:,:,2))
    CASE(SINHSPHERICAL)
      CALL PositionVector_sinhspher(GetScale(this),coords(:,:,1), &
                                   posvec(:,:,1),posvec(:,:,2)) 
    CASE(BIANGLESPHERICAL)
       CALL PositionVector_bianglespher(coords(:,:,1),posvec(:,:,1),posvec(:,:,2), &
                                        posvec(:,:,3))
    CASE(OBLATE_SPHEROIDAL)
      CALL PositionVector_oblatespher(GetScale(this),coords(:,:,1),coords(:,:,2), &
                                   posvec(:,:,1),posvec(:,:,2)) 
    CASE(CHANNEL)
      CALL PositionVector_channel(GetScale(this),coords(:,:,1),coords(:,:,2), &
                                        posvec(:,:,1),posvec(:,:,2)) 
    END SELECT    
  END SUBROUTINE PositionVector_1


  !> \public Compute position vector components
  !!
  !! Computes the curvilinear position vector components with respect to
  !! the given geometry.
  PURE SUBROUTINE PositionVector_2(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(IN) :: this
    REAL, INTENT(IN), DIMENSION(:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:,:) :: posvec
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
      posvec(:,:,:,:) = coords(:,:,:,:)
    CASE(SINHCARTESIAN)
      CALL PositionVector_sinhcartesian(GetScale(this),coords(:,:,:,1),coords(:,:,:,2), &
                                        posvec(:,:,:,1),posvec(:,:,:,2)) 
    CASE(POLAR)
      CALL PositionVector_polar(coords(:,:,:,1),posvec(:,:,:,1),posvec(:,:,:,2)) 
    CASE(LOGPOLAR)
      CALL PositionVector_logpolar(GetScale(this),coords(:,:,:,1),posvec(:,:,:,1), &
                                   posvec(:,:,:,2)) 
    CASE(TANPOLAR)
      CALL PositionVector_tanpolar(GetScale(this),coords(:,:,:,1),posvec(:,:,:,1), &
                                   posvec(:,:,:,2)) 
    CASE(SINHPOLAR)
      CALL PositionVector_sinhpolar(GetScale(this),GetScale(this,2),GetScale(this,3), &
                                   coords(:,:,:,1),posvec(:,:,:,1),posvec(:,:,:,2)) 
    CASE(SINHTANHPOLAR)
      CALL PositionVector_sinhtanh(GetScale(this),GetScale(this,2),GetScale(this,3), &
                                   coords(:,:,:,1),posvec(:,:,:,1),posvec(:,:,:,2)) 
    CASE(POLYPOLAR)
      CALL PositionVector_polypolar(GetScale(this),coords(:,:,:,1),posvec(:,:,:,1), &
                                    posvec(:,:,:,2)) 
    CASE(ELLIPTIC)
      CALL PositionVector_elliptic(GetScale(this),coords(:,:,:,1),coords(:,:,:,2), &
                                   posvec(:,:,:,1),posvec(:,:,:,2)) 

    CASE(CYLINDRICAL)
      posvec(:,:,:,:) = coords(:,:,:,:)
    CASE(TANCYLINDRICAL)
      CALL PositionVector_tancyl(GetScale(this),coords(:,:,:,1),coords(:,:,:,2), &
                                 posvec(:,:,:,1),posvec(:,:,:,2)) 
    CASE(LNCOSHCYLINDRICAL)
      CALL PositionVector_lncoshcyl(GetScale(this),GetScale(this,2),GetScale(this,3), &
                       coords(:,:,:,1),coords(:,:,:,2),posvec(:,:,:,1),posvec(:,:,:,2))
    CASE(SPHERICAL)
      CALL PositionVector_spherical(coords(:,:,:,1),posvec(:,:,:,1),posvec(:,:,:,2))
    CASE(SINHSPHERICAL)
      CALL PositionVector_sinhspher(GetScale(this),coords(:,:,:,1), &
                                   posvec(:,:,:,1),posvec(:,:,:,2)) 
    CASE(BIANGLESPHERICAL)
       CALL PositionVector_bianglespher(coords(:,:,:,1),posvec(:,:,:,1), &
                                        posvec(:,:,:,2),posvec(:,:,:,3))
    CASE(OBLATE_SPHEROIDAL)
      CALL PositionVector_oblatespher(GetScale(this),coords(:,:,:,1),coords(:,:,:,2), &
                                   posvec(:,:,:,1),posvec(:,:,:,2)) 
    CASE(CHANNEL)
      CALL PositionVector_channel(GetScale(this),coords(:,:,:,1),coords(:,:,:,2), &
                                        posvec(:,:,:,1),posvec(:,:,:,2)) 
    END SELECT    
  END SUBROUTINE PositionVector_2


  !> \public Convert curvilinear to cartesian coordinates
  PURE SUBROUTINE Convert2Cartesian_point(this,xi,eta,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP) :: this
    REAL        :: xi,eta,x,y
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,xi,eta
    INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       x = xi
       y = eta
    CASE(SINHCARTESIAN)
       CALL Convert2Cartesian_sinhcart(GetScale(this),xi,eta,x,y)
    CASE(POLAR)
       CALL Convert2Cartesian_polar(xi,eta,x,y)
    CASE(LOGPOLAR)
       CALL Convert2Cartesian_logpolar(GetScale(this),xi,eta,x,y)
    CASE(TANPOLAR)
       CALL Convert2Cartesian_tanpolar(GetScale(this),xi,eta,x,y)
    CASE(SINHPOLAR)
       CALL Convert2Cartesian_sinhpolar(GetScale(this),GetScale(this,2), &
            GetScale(this,3),xi,eta,x,y)
    CASE(SINHTANHPOLAR)
       CALL Convert2Cartesian_sinhtanh(GetScale(this),GetScale(this,2), &
            GetScale(this,3),xi,eta,x,y)
    CASE(POLYPOLAR)
       CALL Convert2Cartesian_polypolar(GetScale(this),xi,eta,x,y)
    CASE(ELLIPTIC)
       CALL Convert2Cartesian_elliptic(GetScale(this),xi,eta,x,y)
    CASE(CYLINDRICAL)
       CALL Convert2Cartesian_cylindrical(xi,eta,x,y)
    CASE(TANCYLINDRICAL)
       CALL Convert2Cartesian_tancyl(GetScale(this),xi,eta,x,y)
    CASE(LNCOSHCYLINDRICAL)
       CALL Convert2Cartesian_lncoshcyl(GetScale(this),GetScale(this,2), &
            GetScale(this,3),xi,eta,x,y)
    CASE(SPHERICAL)
       CALL Convert2Cartesian_spherical(xi,eta,x,y)
    CASE(SINHSPHERICAL)
       CALL Convert2Cartesian_sinhspher(GetScale(this),xi,eta,x,y)
!> \todo: Convert2Cartesian_bianglespher
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Cartesian_oblatespher(GetScale(this),xi,eta,x,y)
    CASE(CHANNEL)
       CALL Convert2Cartesian_channel(GetScale(this),xi,eta,x,y)
    END SELECT
  END SUBROUTINE Convert2Cartesian_point


  !> \public Convert curvilinear to cartesian coordinates
  PURE SUBROUTINE Convert2Cartesian_coords_1(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP)     :: this
    REAL, DIMENSION(:,:,:) :: curv, cart
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,curv
    INTENT(OUT)   :: cart
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       cart(:,:,:) = curv(:,:,:)
    CASE(SINHCARTESIAN)
       CALL Convert2Cartesian_sinhcart(GetScale(this),curv(:,:,1),curv(:,:,2),&
            cart(:,:,1),cart(:,:,2))
    CASE(POLAR)
       CALL Convert2Cartesian_polar(curv(:,:,1),curv(:,:,2),cart(:,:,1),cart(:,:,2))
    CASE(LOGPOLAR)
       CALL Convert2Cartesian_logpolar(GetScale(this),curv(:,:,1),curv(:,:,2),&
            cart(:,:,1),cart(:,:,2))
    CASE(TANPOLAR)
       CALL Convert2Cartesian_tanpolar(GetScale(this),curv(:,:,1),curv(:,:,2),&
            cart(:,:,1),cart(:,:,2))
    CASE(SINHPOLAR)
       CALL Convert2Cartesian_sinhpolar(GetScale(this),GetScale(this,2),GetScale(this,3),&
            curv(:,:,1),curv(:,:,2),cart(:,:,1),cart(:,:,2))
    CASE(SINHTANHPOLAR)
       CALL Convert2Cartesian_sinhtanh(GetScale(this),GetScale(this,2),GetScale(this,3),&
            curv(:,:,1),curv(:,:,2),cart(:,:,1),cart(:,:,2))
    CASE(POLYPOLAR)
       CALL Convert2Cartesian_polypolar(GetScale(this),&
            curv(:,:,1),curv(:,:,2),cart(:,:,1),cart(:,:,2))
    CASE(ELLIPTIC)
       CALL Convert2Cartesian_elliptic(GetScale(this),curv(:,:,1),curv(:,:,2),&
            cart(:,:,1),cart(:,:,2))
    CASE(CYLINDRICAL)
       CALL Convert2Cartesian_cylindrical(curv(:,:,1),curv(:,:,2),cart(:,:,1),&
            cart(:,:,2))
    CASE(TANCYLINDRICAL)
       CALL Convert2Cartesian_tancyl(GetScale(this),curv(:,:,1),curv(:,:,2),&
            cart(:,:,1),cart(:,:,2))
    CASE(LNCOSHCYLINDRICAL)
       CALL Convert2Cartesian_lncoshcyl(GetScale(this),GetScale(this,2),GetScale(this,3),&
            curv(:,:,1),curv(:,:,2),cart(:,:,1),cart(:,:,2))
    CASE(SPHERICAL)
       CALL Convert2Cartesian_spherical(curv(:,:,1),curv(:,:,2),cart(:,:,1),&
            cart(:,:,2))
    CASE(SINHSPHERICAL)
       CALL Convert2Cartesian_sinhspher(GetScale(this),curv(:,:,1),curv(:,:,2),&
            cart(:,:,1),cart(:,:,2))
    CASE(BIANGLESPHERICAL)
       CALL Convert2Cartesian_bianglespher(GetScale(this),curv(:,:,1),curv(:,:,2),&
            cart(:,:,1),cart(:,:,2),cart(:,:,3))
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Cartesian_oblatespher(GetScale(this),curv(:,:,1),curv(:,:,2),&
            cart(:,:,1),cart(:,:,2))
    CASE(CHANNEL)
       CALL Convert2Cartesian_channel(GetScale(this),curv(:,:,1),curv(:,:,2),&
            cart(:,:,1),cart(:,:,2))
    END SELECT
    
  END SUBROUTINE Convert2Cartesian_coords_1


  !> \public Convert curvilinear to cartesian coordinates
  PURE SUBROUTINE Convert2Cartesian_coords_2(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP)     :: this
    REAL, DIMENSION(:,:,:,:) :: curv, cart
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,curv
    INTENT(OUT)   :: cart
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       cart(:,:,:,:) = curv(:,:,:,:)
    CASE(SINHCARTESIAN)
       CALL Convert2Cartesian_sinhcart(GetScale(this),curv(:,:,:,1),curv(:,:,:,2),&
            cart(:,:,:,1),cart(:,:,:,2))
    CASE(POLAR)
       CALL Convert2Cartesian_polar(curv(:,:,:,1),curv(:,:,:,2),cart(:,:,:,1),cart(:,:,:,2))
    CASE(LOGPOLAR)
       CALL Convert2Cartesian_logpolar(GetScale(this),curv(:,:,:,1),curv(:,:,:,2),&
            cart(:,:,:,1),cart(:,:,:,2))
    CASE(TANPOLAR)
       CALL Convert2Cartesian_tanpolar(GetScale(this),curv(:,:,:,1),curv(:,:,:,2),&
            cart(:,:,:,1),cart(:,:,:,2))
    CASE(SINHPOLAR)
       CALL Convert2Cartesian_sinhpolar(GetScale(this),GetScale(this,2),GetScale(this,3),&
            curv(:,:,:,1),curv(:,:,:,2),cart(:,:,:,1),cart(:,:,:,2))
    CASE(SINHTANHPOLAR)
       CALL Convert2Cartesian_sinhtanh(GetScale(this),GetScale(this,2),&
            GetScale(this,3),curv(:,:,:,1),curv(:,:,:,2),cart(:,:,:,1),cart(:,:,:,2))
    CASE(POLYPOLAR)
       CALL Convert2Cartesian_polypolar(GetScale(this),curv(:,:,:,1),curv(:,:,:,2),&
            cart(:,:,:,1), cart(:,:,:,2))
    CASE(ELLIPTIC)
       CALL Convert2Cartesian_elliptic(GetScale(this),curv(:,:,:,1),curv(:,:,:,2),&
            cart(:,:,:,1),cart(:,:,:,2))
    CASE(CYLINDRICAL)
       CALL Convert2Cartesian_cylindrical(curv(:,:,:,1),curv(:,:,:,2),cart(:,:,:,1),&
            cart(:,:,:,2))
    CASE(TANCYLINDRICAL)
       CALL Convert2Cartesian_tancyl(GetScale(this),curv(:,:,:,1),curv(:,:,:,2),&
            cart(:,:,:,1),cart(:,:,:,2))
    CASE(LNCOSHCYLINDRICAL)
       CALL Convert2Cartesian_lncoshcyl(GetScale(this),GetScale(this,2),GetScale(this,3),&
            curv(:,:,:,1),curv(:,:,:,2),cart(:,:,:,1),cart(:,:,:,2))
    CASE(SPHERICAL)
       CALL Convert2Cartesian_spherical(curv(:,:,:,1),curv(:,:,:,2),cart(:,:,:,1),&
            cart(:,:,:,2))
    CASE(SINHSPHERICAL)
       CALL Convert2Cartesian_sinhspher(GetScale(this),curv(:,:,:,1),curv(:,:,:,2), &
            cart(:,:,:,1),cart(:,:,:,2))
    CASE(BIANGLESPHERICAL)
       CALL Convert2Cartesian_bianglespher(GetScale(this),curv(:,:,:,1),curv(:,:,:,2),&
            cart(:,:,:,1),cart(:,:,:,2),cart(:,:,:,3))
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Cartesian_oblatespher(GetScale(this),curv(:,:,:,1),curv(:,:,:,2),&
            cart(:,:,:,1),cart(:,:,:,2))
    CASE(CHANNEL)
       CALL Convert2Cartesian_channel(GetScale(this),curv(:,:,:,1),curv(:,:,:,2),&
            cart(:,:,:,1),cart(:,:,:,2))
    END SELECT
  
  END SUBROUTINE Convert2Cartesian_coords_2


  !> \public Convert cartesian to curvilinear coordinates
  PURE SUBROUTINE Convert2Curvilinear_point(this,x,y,xi,eta)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP) :: this
    REAL        :: x,y,xi,eta
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,x,y
    INTENT(OUT) :: xi,eta
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       xi  = x
       eta = y
    CASE(SINHCARTESIAN)
       CALL Convert2Curvilinear_sinhcart(GetScale(this),x,y,xi,eta)
    CASE(POLAR)
       CALL Convert2Curvilinear_polar(x,y,xi,eta)
    CASE(LOGPOLAR)
       CALL Convert2Curvilinear_logpolar(GetScale(this),x,y,xi,eta)
    CASE(TANPOLAR)
       CALL Convert2Curvilinear_tanpolar(GetScale(this),x,y,xi,eta)
    CASE(SINHPOLAR)
       CALL Convert2Curvilinear_sinhpolar(GetScale(this),GetScale(this,2), &
            GetScale(this,3),x,y,xi,eta)
    CASE(SINHTANHPOLAR)
       CALL Convert2Curvilinear_sinhtanh(GetScale(this),GetScale(this,2), &
            GetScale(this,3),x,y,xi,eta)
    CASE(POLYPOLAR)
       CALL Convert2Curvilinear_polypolar(GetScale(this),x,y,xi,eta)
    CASE(ELLIPTIC)
       CALL Convert2Curvilinear_elliptic(GetScale(this),x,y,xi,eta)
    CASE(CYLINDRICAL)
       CALL Convert2Curvilinear_cylindrical(x,y,xi,eta)
    CASE(TANCYLINDRICAL)
       CALL Convert2Curvilinear_tancyl(GetScale(this),x,y,xi,eta)
    CASE(LNCOSHCYLINDRICAL)
       CALL Convert2Curvilinear_lncoshcyl(GetScale(this),GetScale(this,2), &
            GetScale(this,3),x,y,xi,eta)
    CASE(SPHERICAL)
       CALL Convert2Curvilinear_spherical(x,y,xi,eta)
    CASE(SINHSPHERICAL)
       CALL Convert2Curvilinear_sinhspher(GetScale(this),x,y,xi,eta)
!> \todo: Convert2Curvilinear_bianglespher
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Curvilinear_oblatespher(GetScale(this),x,y,xi,eta)
    CASE(CHANNEL)
       CALL Convert2Curvilinear_channel(GetScale(this),x,y,xi,eta)
    END SELECT
  END SUBROUTINE Convert2Curvilinear_point


  !> \public Convert cartesian to curvilinear coordinates
  PURE SUBROUTINE Convert2Curvilinear_coords_1(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP)     :: this
    REAL, DIMENSION(:,:,:) :: cart, curv
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,cart
    INTENT(OUT)   :: curv
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       curv(:,:,:) = cart(:,:,:)
    CASE(SINHCARTESIAN)
       CALL Convert2Curvilinear_sinhcart(GetScale(this),cart(:,:,1),cart(:,:,2),&
            curv(:,:,1),curv(:,:,2))
    CASE(POLAR)
       CALL Convert2Curvilinear_polar(cart(:,:,1),cart(:,:,2),curv(:,:,1),curv(:,:,2))
    CASE(LOGPOLAR)
       CALL Convert2Curvilinear_logpolar(GetScale(this),cart(:,:,1),cart(:,:,2),&
            curv(:,:,1),curv(:,:,2))
    CASE(TANPOLAR)
       CALL Convert2Curvilinear_tanpolar(GetScale(this),cart(:,:,1),cart(:,:,2),&
            curv(:,:,1),curv(:,:,2))
    CASE(SINHPOLAR)
       CALL Convert2Curvilinear_sinhpolar(GetScale(this),GetScale(this,2),GetScale(this,3),&
            cart(:,:,1),cart(:,:,2),curv(:,:,1),curv(:,:,2))
    CASE(SINHTANHPOLAR)
       CALL Convert2Curvilinear_sinhtanh(GetScale(this),GetScale(this,2),&
            GetScale(this,3),cart(:,:,1),cart(:,:,2),curv(:,:,1),curv(:,:,2))
    CASE(POLYPOLAR)
       CALL Convert2Curvilinear_polypolar(GetScale(this),cart(:,:,1),cart(:,:,2),&
            curv(:,:,1),curv(:,:,2))
    CASE(ELLIPTIC)
       CALL Convert2Curvilinear_elliptic(GetScale(this),cart(:,:,1),cart(:,:,2), &
            curv(:,:,1),curv(:,:,2))
    CASE(CYLINDRICAL)
       CALL Convert2Curvilinear_cylindrical(cart(:,:,1),cart(:,:,2),curv(:,:,1),curv(:,:,2))
    CASE(TANCYLINDRICAL)
       CALL Convert2Curvilinear_tancyl(GetScale(this),cart(:,:,1),cart(:,:,2), &
            curv(:,:,1),curv(:,:,2))
    CASE(LNCOSHCYLINDRICAL)
       CALL Convert2Curvilinear_lncoshcyl(GetScale(this),GetScale(this,2),GetScale(this,3),&
            cart(:,:,1),cart(:,:,2),curv(:,:,1),curv(:,:,2))
    CASE(SPHERICAL)
       CALL Convert2Curvilinear_spherical(cart(:,:,1),cart(:,:,2),curv(:,:,1),curv(:,:,2))
    CASE(SINHSPHERICAL)
       CALL Convert2Curvilinear_sinhspher(GetScale(this),cart(:,:,1),cart(:,:,2), &
            curv(:,:,1),curv(:,:,2))
    CASE(BIANGLESPHERICAL)
       CALL Convert2Curv_bianglespher(GetScale(this),cart(:,:,1),cart(:,:,2),cart(:,:,3),&
            curv(:,:,1),curv(:,:,2))
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Curvilinear_oblatespher(GetScale(this),cart(:,:,1),cart(:,:,2), &
            curv(:,:,1),curv(:,:,2))
    CASE(CHANNEL)
       CALL Convert2Curvilinear_channel(GetScale(this),cart(:,:,1),cart(:,:,2), &
            curv(:,:,1),curv(:,:,2))
    END SELECT
    
  END SUBROUTINE Convert2Curvilinear_coords_1


  !> \public Convert cartesian to curvilinear coordinates
  PURE SUBROUTINE Convert2Curvilinear_coords_2(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP)     :: this
    REAL, DIMENSION(:,:,:,:) :: cart, curv
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,cart
    INTENT(OUT)   :: curv
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       curv(:,:,:,:) = cart(:,:,:,:)
    CASE(SINHCARTESIAN)
       CALL Convert2Curvilinear_sinhcart(GetScale(this),cart(:,:,:,1),cart(:,:,:,2),&
            curv(:,:,:,1),curv(:,:,:,2))
    CASE(POLAR)
       CALL Convert2Curvilinear_polar(cart(:,:,:,1),cart(:,:,:,2),curv(:,:,:,1),curv(:,:,:,2))
    CASE(LOGPOLAR)
       CALL Convert2Curvilinear_logpolar(GetScale(this),cart(:,:,:,1),cart(:,:,:,2),&
            curv(:,:,:,1),curv(:,:,:,2))
    CASE(TANPOLAR)
       CALL Convert2Curvilinear_tanpolar(GetScale(this),cart(:,:,:,1),cart(:,:,:,2),&
            curv(:,:,:,1),curv(:,:,:,2))
    CASE(SINHPOLAR)
       CALL Convert2Curvilinear_sinhpolar(GetScale(this),GetScale(this,2),GetScale(this,3),&
            cart(:,:,:,1),cart(:,:,:,2),curv(:,:,:,1),curv(:,:,:,2))
    CASE(SINHTANHPOLAR)
       CALL Convert2Curvilinear_sinhtanh(GetScale(this),GetScale(this,2),&
            GetScale(this,3),cart(:,:,:,1),cart(:,:,:,2),curv(:,:,:,1),curv(:,:,:,2))
    CASE(POLYPOLAR)
       CALL Convert2Curvilinear_polypolar(GetScale(this),cart(:,:,:,1),cart(:,:,:,2),&
            curv(:,:,:,1),curv(:,:,:,2))
    CASE(ELLIPTIC)
       CALL Convert2Curvilinear_elliptic(GetScale(this),cart(:,:,:,1),cart(:,:,:,2), &
            curv(:,:,:,1),curv(:,:,:,2))
    CASE(CYLINDRICAL)
       CALL Convert2Curvilinear_cylindrical(cart(:,:,:,1),cart(:,:,:,2),curv(:,:,:,1),curv(:,:,:,2))
    CASE(TANCYLINDRICAL)
       CALL Convert2Curvilinear_tancyl(GetScale(this),cart(:,:,:,1),cart(:,:,:,2), &
            curv(:,:,:,1),curv(:,:,:,2))
    CASE(LNCOSHCYLINDRICAL)
       CALL Convert2Curvilinear_lncoshcyl(GetScale(this),GetScale(this,2),GetScale(this,3),&
            cart(:,:,:,1),cart(:,:,:,2),curv(:,:,:,1),curv(:,:,:,2))
    CASE(SPHERICAL)
       CALL Convert2Curvilinear_spherical(cart(:,:,:,1),cart(:,:,:,2),curv(:,:,:,1),curv(:,:,:,2))
    CASE(SINHSPHERICAL)
       CALL Convert2Curvilinear_sinhspher(GetScale(this),cart(:,:,:,1),cart(:,:,:,2), &
            curv(:,:,:,1),curv(:,:,:,2))
    CASE(BIANGLESPHERICAL)
       CALL Convert2Curv_bianglespher(GetScale(this),cart(:,:,:,1),cart(:,:,:,2),cart(:,:,:,3),&
       curv(:,:,:,1),curv(:,:,:,2) )
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Curvilinear_oblatespher(GetScale(this),cart(:,:,:,1),cart(:,:,:,2), &
            curv(:,:,:,1),curv(:,:,:,2))
    CASE(CHANNEL)
       CALL Convert2Curvilinear_channel(GetScale(this),cart(:,:,:,1),cart(:,:,:,2), &
            curv(:,:,:,1),curv(:,:,:,2))
    END SELECT
  
  END SUBROUTINE Convert2Curvilinear_coords_2


  !> \public Convert curvilinear vector components to cartesian vector components
  PURE SUBROUTINE Convert2Cartesian_vectors_1(this,curv,v_curv,v_cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP)     :: this
    REAL, DIMENSION(:,:,:) :: curv, v_curv, v_cart
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,curv,v_curv
    INTENT(OUT)   :: v_cart
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       v_cart(:,:,:) = v_curv(:,:,:)
    CASE(SINHCARTESIAN)
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
       CALL Convert2Cartesian_sinhpolar(curv(:,:,2),v_curv(:,:,1),v_curv(:,:,2),&
            v_cart(:,:,1),v_cart(:,:,2))
    CASE(SINHTANHPOLAR)
       CALL Convert2Cartesian_sinhtanh(GetScale(this),GetScale(this,2),&
            GetScale(this,3),curv(:,:,2),v_curv(:,:,1),v_curv(:,:,2),&
            v_cart(:,:,1),v_cart(:,:,2))
    CASE(POLYPOLAR)
       CALL Convert2Cartesian_polypolar(GetScale(this),curv(:,:,2),v_curv(:,:,1),&
            v_curv(:,:,2),v_cart(:,:,1),v_cart(:,:,2))
    CASE(ELLIPTIC)
       CALL Convert2Cartesian_elliptic(curv(:,:,1),curv(:,:,2),v_curv(:,:,1), &
            v_curv(:,:,2),v_cart(:,:,1),v_cart(:,:,2))
    CASE(CYLINDRICAL)
       CALL Convert2Cartesian_cylindrical(v_curv(:,:,1),v_curv(:,:,2), &
            v_cart(:,:,1),v_cart(:,:,2))
    CASE(TANCYLINDRICAL)
       CALL Convert2Cartesian_tancyl(v_curv(:,:,1),v_curv(:,:,2), &
            v_cart(:,:,1),v_cart(:,:,2))
    CASE(LNCOSHCYLINDRICAL)
       CALL Convert2Cartesian_lncoshcyl(v_curv(:,:,1),v_curv(:,:,2), &
            v_cart(:,:,1),v_cart(:,:,2))
    CASE(SPHERICAL)
       CALL Convert2Cartesian_spherical(GetScale(this),curv(:,:,2),v_curv(:,:,1), &
            v_curv(:,:,2),v_cart(:,:,1),v_cart(:,:,2))
    CASE(SINHSPHERICAL)
       CALL Convert2Cartesian_sinhspher(GetScale(this),curv(:,:,2),v_curv(:,:,1), &
            v_curv(:,:,2),v_cart(:,:,1),v_cart(:,:,2))
    CASE(BIANGLESPHERICAL)
       CALL Convert2Cartesian_bianglespher(GetScale(this),curv(:,:,1),curv(:,:,2), &
            v_curv(:,:,1),v_curv(:,:,2),v_curv(:,:,3), &
            v_cart(:,:,1),v_cart(:,:,2),v_cart(:,:,3))
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Cartesian_oblatespher(GetScale(this),curv(:,:,1),curv(:,:,2), &
            v_curv(:,:,1),v_curv(:,:,2),v_cart(:,:,1),v_cart(:,:,2))
    CASE(CHANNEL)
       CALL Convert2Cartesian_channel(v_curv(:,:,1),v_curv(:,:,2), &
            v_cart(:,:,1),v_cart(:,:,2))
    END SELECT
    
  END SUBROUTINE Convert2Cartesian_vectors_1


  !> \public Convert curvilinear vector components to cartesian vector components
  PURE SUBROUTINE Convert2Cartesian_vectors_2(this,curv,v_curv,v_cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP)     :: this
    REAL, DIMENSION(:,:,:,:) :: curv, v_curv, v_cart
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,curv,v_curv
    INTENT(OUT)   :: v_cart
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       v_cart(:,:,:,:) = v_curv(:,:,:,:)
    CASE(SINHCARTESIAN)
       v_cart(:,:,:,:) = v_curv(:,:,:,:)
    CASE(POLAR)
       CALL Convert2Cartesian_polar(curv(:,:,:,2),v_curv(:,:,:,1),v_curv(:,:,:,2), &
            v_cart(:,:,:,1),v_cart(:,:,:,2))
    CASE(LOGPOLAR)
       CALL Convert2Cartesian_logpolar(GetScale(this),curv(:,:,:,2),v_curv(:,:,:,1),&
            v_curv(:,:,:,2),v_cart(:,:,:,1),v_cart(:,:,:,2))
    CASE(TANPOLAR)
       CALL Convert2Cartesian_tanpolar(GetScale(this),curv(:,:,:,2),v_curv(:,:,:,1),&
            v_curv(:,:,:,2),v_cart(:,:,:,1),v_cart(:,:,:,2))
    CASE(SINHPOLAR)
       CALL Convert2Cartesian_sinhpolar(curv(:,:,:,2),v_curv(:,:,:,1),v_curv(:,:,:,2),&
            v_cart(:,:,:,1),v_cart(:,:,:,2))
    CASE(SINHTANHPOLAR)
       CALL Convert2Cartesian_sinhtanh(GetScale(this),GetScale(this,2),&
            GetScale(this,3),curv(:,:,:,2),v_curv(:,:,:,1),v_curv(:,:,:,2),&
            v_cart(:,:,:,1),v_cart(:,:,:,2))
    CASE(POLYPOLAR)
       CALL Convert2Cartesian_polypolar(GetScale(this),curv(:,:,:,2),v_curv(:,:,:,1),&
            v_curv(:,:,:,2),v_cart(:,:,:,1),v_cart(:,:,:,2))
    CASE(ELLIPTIC)
       CALL Convert2Cartesian_elliptic(curv(:,:,:,1),curv(:,:,:,2),v_curv(:,:,:,1), &
            v_curv(:,:,:,2),v_cart(:,:,:,1),v_cart(:,:,:,2))
    CASE(CYLINDRICAL)
       CALL Convert2Cartesian_cylindrical(v_curv(:,:,:,1),v_curv(:,:,:,2), &
            v_cart(:,:,:,1),v_cart(:,:,:,2))
    CASE(TANCYLINDRICAL)
       CALL Convert2Cartesian_tancyl(v_curv(:,:,:,1),v_curv(:,:,:,2), &
            v_cart(:,:,:,1),v_cart(:,:,:,2))
    CASE(LNCOSHCYLINDRICAL)
       CALL Convert2Cartesian_lncoshcyl(v_curv(:,:,:,1),v_curv(:,:,:,2), &
            v_cart(:,:,:,1),v_cart(:,:,:,2))
    CASE(SPHERICAL)
       CALL Convert2Cartesian_spherical(GetScale(this),curv(:,:,:,2),v_curv(:,:,:,1), &
            v_curv(:,:,:,2),v_cart(:,:,:,1),v_cart(:,:,:,2))
    CASE(SINHSPHERICAL)
       CALL Convert2Cartesian_sinhspher(GetScale(this),curv(:,:,:,2),v_curv(:,:,:,1), &
            v_curv(:,:,:,2),v_cart(:,:,:,1),v_cart(:,:,:,2))
    CASE(BIANGLESPHERICAL)
       CALL Convert2Cartesian_bianglespher(GetScale(this),curv(:,:,:,1), &
            curv(:,:,:,2),v_curv(:,:,:,1),v_curv(:,:,:,2),v_curv(:,:,:,3), &
            v_cart(:,:,:,1),v_cart(:,:,:,2),v_cart(:,:,:,3))
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Cartesian_oblatespher(GetScale(this),curv(:,:,:,1),curv(:,:,:,2), &
            v_curv(:,:,:,1),v_curv(:,:,:,2),v_cart(:,:,:,1),v_cart(:,:,:,2))
    CASE(CHANNEL)
       CALL Convert2Cartesian_channel(v_curv(:,:,:,1),v_curv(:,:,:,2), &
            v_cart(:,:,:,1),v_cart(:,:,:,2))
    END SELECT
    
  END SUBROUTINE Convert2Cartesian_vectors_2

  !> \public Convert curvilinear vector components to cartesian vector components
  PURE SUBROUTINE Convert2Cartesian_vectors_3(this,curv,v_curv,v_cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP)     :: this
    REAL, DIMENSION(:) :: curv, v_curv, v_cart
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,curv,v_curv
    INTENT(OUT)   :: v_cart
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       v_cart(:) = v_curv(:)
    CASE(POLAR)
       CALL Convert2Cartesian_polar(curv(2),v_curv(1),v_curv(2), &
            v_cart(1),v_cart(2))
    CASE(LOGPOLAR)
       CALL Convert2Cartesian_logpolar(GetScale(this),curv(2),v_curv(1),&
            v_curv(2),v_cart(1),v_cart(2))
    CASE(TANPOLAR)
       CALL Convert2Cartesian_tanpolar(GetScale(this),curv(2),v_curv(1),&
            v_curv(2),v_cart(1),v_cart(2))
    CASE(SINHPOLAR)
       CALL Convert2Cartesian_sinhpolar(curv(2),v_curv(1),v_curv(2),&
            v_cart(1),v_cart(2))
    CASE(SINHTANHPOLAR)
       CALL Convert2Cartesian_sinhtanh(GetScale(this),GetScale(this,2),&
            GetScale(this,3),curv(2),v_curv(1),v_curv(2),&
            v_cart(1),v_cart(2))
    CASE(POLYPOLAR)
       CALL Convert2Cartesian_polypolar(GetScale(this),curv(2),v_curv(1),&
            v_curv(2),v_cart(1),v_cart(2))
    CASE(CYLINDRICAL)
       CALL Convert2Cartesian_cylindrical(v_curv(1),v_curv(2), &
            v_cart(1),v_cart(2))
    CASE(TANCYLINDRICAL)
       CALL Convert2Cartesian_tancyl(v_curv(1),v_curv(2), &
            v_cart(1),v_cart(2))
    CASE(LNCOSHCYLINDRICAL)
       CALL Convert2Cartesian_lncoshcyl(v_curv(1),v_curv(2), &
            v_cart(1),v_cart(2))
    CASE(SPHERICAL)
       CALL Convert2Cartesian_spherical(GetScale(this),curv(2),v_curv(1), &
            v_curv(2),v_cart(1),v_cart(2))
    CASE(SINHSPHERICAL)
       CALL Convert2Cartesian_sinhspher(GetScale(this),curv(2),v_curv(1), &
            v_curv(2),v_cart(1),v_cart(2))
    CASE(BIANGLESPHERICAL)
       CALL Convert2Cartesian_bianglespher(GetScale(this),curv(1), &
            curv(2),v_curv(1),v_curv(2),v_curv(3),v_cart(1), &
            v_cart(2),v_cart(3))
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Cartesian_oblatespher(GetScale(this),curv(1),curv(2), &
            v_curv(1),v_curv(2),v_cart(1),v_cart(2))
    CASE(CHANNEL)
       CALL Convert2Cartesian_channel(v_curv(1),v_curv(2), &
            v_cart(1),v_cart(2))
    END SELECT
    
  END SUBROUTINE Convert2Cartesian_vectors_3

  !> \public Convert cartesian vector components to curvilinear vector components
  PURE SUBROUTINE Convert2Curvilinear_vectors_1(this,curv,v_cart,v_curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP)     :: this
    REAL, DIMENSION(:,:,:) :: curv, v_curv, v_cart
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,curv,v_cart
    INTENT(OUT)   :: v_curv
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       v_curv(:,:,:) = v_cart(:,:,:)
    CASE(SINHCARTESIAN)
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
       CALL Convert2Curvilinear_sinhpolar(curv(:,:,2),v_cart(:,:,1),v_cart(:,:,2),&
            v_curv(:,:,1),v_curv(:,:,2))
    CASE(SINHTANHPOLAR)
       CALL Convert2Curvilinear_sinhtanh(GetScale(this),curv(:,:,2),&
            v_cart(:,:,1),v_cart(:,:,2),v_curv(:,:,1),v_curv(:,:,2))
    CASE(POLYPOLAR)
       CALL Convert2Curvilinear_polypolar(GetScale(this),curv(:,:,2),&
            v_cart(:,:,1),v_cart(:,:,2),v_curv(:,:,1),v_curv(:,:,2))
    CASE(ELLIPTIC)
       CALL Convert2Curvilinear_elliptic(curv(:,:,1),curv(:,:,2),v_cart(:,:,1), &
            v_cart(:,:,2),v_curv(:,:,1),v_curv(:,:,2))
    CASE(CYLINDRICAL)
       CALL Convert2Curvilinear_cylindrical(v_cart(:,:,1),v_cart(:,:,2), &
            v_curv(:,:,1),v_curv(:,:,2))
    CASE(TANCYLINDRICAL)
       CALL Convert2Curvilinear_tancyl(v_cart(:,:,1),v_cart(:,:,2), &
            v_curv(:,:,1),v_curv(:,:,2))
    CASE(LNCOSHCYLINDRICAL)
       CALL Convert2Curvilinear_lncoshcyl(v_cart(:,:,1),v_cart(:,:,2), &
            v_curv(:,:,1),v_curv(:,:,2))
    CASE(SPHERICAL)
       CALL Convert2Curvilinear_spherical(GetScale(this),curv(:,:,2),v_cart(:,:,1), &
            v_cart(:,:,2),v_curv(:,:,1),v_curv(:,:,2))
    CASE(SINHSPHERICAL)
       CALL Convert2Curvilinear_sinhspher(GetScale(this),curv(:,:,2),v_cart(:,:,1), &
            v_cart(:,:,2),v_curv(:,:,1),v_curv(:,:,2))
    CASE(BIANGLESPHERICAL)
       CALL Convert2Curv_bianglespher(GetScale(this),curv(:,:,1),curv(:,:,2),&
            v_cart(:,:,1),v_cart(:,:,2),v_cart(:,:,3),&
            v_curv(:,:,1),v_curv(:,:,2),v_curv(:,:,3))
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Curvilinear_oblatespher(curv(:,:,1),curv(:,:,2), &
            v_cart(:,:,1),v_cart(:,:,2),v_curv(:,:,1),v_curv(:,:,2))
    CASE(CHANNEL)
       CALL Convert2Curvilinear_channel(v_cart(:,:,1),v_cart(:,:,2), &
            v_curv(:,:,1),v_curv(:,:,2))
    END SELECT
    
  END SUBROUTINE Convert2Curvilinear_vectors_1

  !> \public Convert cartesian vector components to curvilinear vector components
  PURE SUBROUTINE Convert2Curvilinear_vectors_2(this,curv,v_cart,v_curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP)     :: this
    REAL, DIMENSION(:,:,:,:) :: curv, v_curv, v_cart
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,curv,v_cart
    INTENT(OUT)   :: v_curv
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       v_curv(:,:,:,:) = v_cart(:,:,:,:)
    CASE(SINHCARTESIAN)
       v_curv(:,:,:,:) = v_cart(:,:,:,:)
    CASE(POLAR)
       CALL Convert2Curvilinear_polar(curv(:,:,:,2),v_cart(:,:,:,1),v_cart(:,:,:,2), &
            v_curv(:,:,:,1),v_curv(:,:,:,2))
    CASE(LOGPOLAR)
       CALL Convert2Curvilinear_logpolar(GetScale(this),curv(:,:,:,2),v_cart(:,:,:,1),&
            v_cart(:,:,:,2),v_curv(:,:,:,1),v_curv(:,:,:,2))
    CASE(TANPOLAR)
       CALL Convert2Curvilinear_tanpolar(GetScale(this),curv(:,:,:,2),v_cart(:,:,:,1),&
            v_cart(:,:,:,2),v_curv(:,:,:,1),v_curv(:,:,:,2))
    CASE(SINHPOLAR)
       CALL Convert2Curvilinear_sinhpolar(curv(:,:,:,2),v_cart(:,:,:,1),v_cart(:,:,:,2),&
            v_curv(:,:,:,1),v_curv(:,:,:,2))
    CASE(SINHTANHPOLAR)
       CALL Convert2Curvilinear_sinhtanh(GetScale(this),curv(:,:,:,2),&
            v_cart(:,:,:,1),v_cart(:,:,:,2),v_curv(:,:,:,1),v_curv(:,:,:,2))
    CASE(POLYPOLAR)
       CALL Convert2Curvilinear_polypolar(GetScale(this),curv(:,:,:,2),&
            v_cart(:,:,:,1),v_cart(:,:,:,2),v_curv(:,:,:,1),v_curv(:,:,:,2))
    CASE(ELLIPTIC)
       CALL Convert2Curvilinear_elliptic(curv(:,:,:,1),curv(:,:,:,2),v_cart(:,:,:,1), &
            v_cart(:,:,:,2),v_curv(:,:,:,1),v_curv(:,:,:,2))
    CASE(CYLINDRICAL)
       CALL Convert2Curvilinear_cylindrical(v_cart(:,:,:,1),v_cart(:,:,:,2), &
            v_curv(:,:,:,1),v_curv(:,:,:,2))
    CASE(TANCYLINDRICAL)
       CALL Convert2Curvilinear_tancyl(v_cart(:,:,:,1),v_cart(:,:,:,2), &
            v_curv(:,:,:,1),v_curv(:,:,:,2))
    CASE(LNCOSHCYLINDRICAL)
       CALL Convert2Curvilinear_lncoshcyl(v_cart(:,:,:,1),v_cart(:,:,:,2), &
            v_curv(:,:,:,1),v_curv(:,:,:,2))
    CASE(SPHERICAL)
       CALL Convert2Curvilinear_spherical(GetScale(this),curv(:,:,:,2),v_cart(:,:,:,1), &
            v_cart(:,:,:,2),v_curv(:,:,:,1),v_curv(:,:,:,2))
    CASE(SINHSPHERICAL)
       CALL Convert2Curvilinear_sinhspher(GetScale(this),curv(:,:,:,2),v_cart(:,:,:,1), &
            v_cart(:,:,:,2),v_curv(:,:,:,1),v_curv(:,:,:,2))
    CASE(BIANGLESPHERICAL)
       CALL Convert2Curv_bianglespher(GetScale(this),curv(:,:,:,1),&
            curv(:,:,:,2),v_cart(:,:,:,1),v_cart(:,:,:,2),v_cart(:,:,:,3),&
            v_curv(:,:,:,1),v_curv(:,:,:,2),v_curv(:,:,:,3))
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Curvilinear_oblatespher(curv(:,:,:,1),curv(:,:,:,2), &
            v_cart(:,:,:,1),v_cart(:,:,:,2),v_curv(:,:,:,1),v_curv(:,:,:,2))
    CASE(CHANNEL)
       CALL Convert2Curvilinear_channel(v_cart(:,:,:,1),v_cart(:,:,:,2), &
            v_curv(:,:,:,1),v_curv(:,:,:,2))
    END SELECT
    
  END SUBROUTINE Convert2Curvilinear_vectors_2

  !> \public Convert cartesian vector components to curvilinear vector components
  PURE SUBROUTINE Convert2Curvilinear_vectors_3(this,curv,v_cart,v_curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP)     :: this
    REAL, DIMENSION(:) :: curv, v_curv, v_cart
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,curv,v_cart
    INTENT(OUT)   :: v_curv
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(CARTESIAN)
       v_curv(:) = v_cart(:)
    CASE(POLAR)
       CALL Convert2Curvilinear_polar(curv(2),v_cart(1),v_cart(2), &
            v_curv(1),v_curv(2))
    CASE(LOGPOLAR)
       CALL Convert2Curvilinear_logpolar(GetScale(this),curv(2),v_cart(1),&
            v_cart(2),v_curv(1),v_curv(2))
    CASE(TANPOLAR)
       CALL Convert2Curvilinear_tanpolar(GetScale(this),curv(2),v_cart(1),&
            v_cart(2),v_curv(1),v_curv(2))
    CASE(SINHPOLAR)
       CALL Convert2Curvilinear_sinhpolar(curv(2),v_cart(1),v_cart(2),&
            v_curv(1),v_curv(2))
    CASE(SINHTANHPOLAR)
       CALL Convert2Curvilinear_sinhtanh(GetScale(this),curv(2),&
            v_cart(1),v_cart(2),v_curv(1),v_curv(2))
    CASE(POLYPOLAR)
       CALL Convert2Curvilinear_polypolar(GetScale(this),curv(2),&
            v_cart(1),v_cart(2),v_curv(1),v_curv(2))
    CASE(CYLINDRICAL)
       CALL Convert2Curvilinear_cylindrical(v_cart(1),v_cart(2), &
            v_curv(1),v_curv(2))
    CASE(TANCYLINDRICAL)
       CALL Convert2Curvilinear_tancyl(v_cart(1),v_cart(2), &
            v_curv(1),v_curv(2))
    CASE(LNCOSHCYLINDRICAL)
       CALL Convert2Curvilinear_lncoshcyl(v_cart(1),v_cart(2), &
            v_curv(1),v_curv(2))
    CASE(SPHERICAL)
       CALL Convert2Curvilinear_spherical(GetScale(this),curv(2),v_cart(1), &
            v_cart(2),v_curv(1),v_curv(2))
    CASE(SINHSPHERICAL)
       CALL Convert2Curvilinear_sinhspher(GetScale(this),curv(2),v_cart(1), &
            v_cart(2),v_curv(1),v_curv(2))
    CASE(BIANGLESPHERICAL)
       CALL Convert2Curv_bianglespher(GetScale(this),curv(1),curv(2),&
            v_cart(1),v_cart(2),v_cart(3),&
            v_curv(1),v_curv(2),v_curv(3))
    CASE(OBLATE_SPHEROIDAL)
       CALL Convert2Curvilinear_oblatespher(curv(1),curv(2), &
            v_cart(1),v_cart(2),v_curv(1),v_curv(2))
    CASE(CHANNEL)
       CALL Convert2Curvilinear_channel(v_cart(1),v_cart(2), &
            v_curv(1),v_curv(2))
    END SELECT
    
  END SUBROUTINE Convert2Curvilinear_vectors_3

  !> \public Destructor of generic geometry module
  SUBROUTINE CloseGeometry(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF (.NOT.Initialized(this)) &
         CALL Error(this,"CloseGeometry","not initialized")
    CALL CloseGeometry_common(this)
  END SUBROUTINE CloseGeometry

END MODULE geometry_generic
