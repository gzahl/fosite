!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_diskthomson.f90                                           #
!#                                                                           #
!# Copyright (C) 2007-2008                                                   #
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
! source terms module for gravitational acceleration due to
! a point mass at the center of the coordinate system
!----------------------------------------------------------------------------!
MODULE sources_diskthomson
  USE sources_pointmass, ExternalSources_diskthomson => Externalsources_pointmass, &
       CloseSources_diskthomson => CloseSources_pointmass
  USE fluxes_common, ONLY : Fluxes_TYP
  USE physics_generic
  USE mesh_generic
  USE integration
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  REAL, PARAMETER :: EPS  = 1.0D-04              ! precision for integration !
  REAL, PARAMETER :: TINY = 1.0E-30              ! to avoid division by 0    !
  CHARACTER(LEN=32), PARAMETER :: source_name = "thomoson scat. of disk rad."
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! methods
       InitSources_diskthomson, &
       ExternalSources_diskthomson, &
       CloseSources_diskthomson
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources_diskthomson(this,Mesh,Physics,stype,mass,mdot,s0,s1)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    INTEGER           :: stype
    REAL              :: mass,mdot,s0,s1
    !------------------------------------------------------------------------!
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
                      :: curv, cart, accel
    REAL              :: params(2)
    REAL              :: rs, s0rs, x0, xm, factor
    INTEGER           :: err
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,stype,mass,mdot,s0,s1
    !------------------------------------------------------------------------!
    CALL InitSources(this,stype,source_name)

    ! mass of central object (black hole)
    this%mass = mass
    ! mass accretion rate of the disk
    this%mdot = mdot

    ! some constants
    rs = 2*Physics%constants%GN * &          ! Schwarzschildradius of the BH !
         (this%mass / (Physics%constants%C**2 + TINY))
    s0rs = s0 / rs                           ! inner radius in terms of R_s  !
    factor = 3.*Physics%constants%KE / &     ! constant factor               !
         (16*PI*s0rs**3 + TINY) * (this%mdot/rs) * (Physics%constants%C/rs)

    ! reserve memory for radiational acceleration
    ALLOCATE(this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
         STAT = err)
    IF (err.NE.0) CALL Error(this,"InitSources_diskthomson", "Unable allocate memory!")

    ! convert to cartesian coordinates
    curv(:,:,:) = Mesh%bcenter(:,:,:)
    CALL Convert2Cartesian(Mesh%geometry,curv,cart)

    ! initialize radiational acceleration
    x0 = s0/s1
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! calculate cartesian components of radiational acceleration
          ! due to Thomson scattering
          params(1) = ABS(cart(i,j,1) / (s0+TINY))          ! = r / s0
          params(2) = ABS(cart(i,j,2) / (s0+TINY))          ! = z / s0
          xm = MIN(1./(params(1)+TINY),1.0)
          ! integrate around the pseudo-singularity at xm 
          accel(i,j,1) = SIGN(1.0,cart(i,j,2)) * factor * (&
               integrate(integrand_rade_z,x0,xm,EPS,params,method=1) &
             + integrate(integrand_rade_z,xm,1.0,EPS,params,method=1))
          accel(i,j,2) = factor * (&
               integrate(integrand_rade_r,x0,xm,EPS,params,method=1) &
             + integrate(integrand_rade_r,xm,1.0,EPS,params,method=1))
       END DO
    END DO

    ! convert to curvilinear vector components
    CALL Convert2Curvilinear(Mesh%geometry,curv,accel,this%accel)
  END SUBROUTINE InitSources_diskthomson


  FUNCTION integrand_rade_r(x,plist) RESULT(fx)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x
    REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
    REAL :: fx
    REAL :: s,s2,z,z2,x2,tmp

    s  = plist(1)
    z  = plist(2)
    x2 = x*x
    s2 = s*s
    z2 = z*z
    tmp = 1.0 + (z2+s2)*x2
    fx = (1.-SQRT(x)) * z*s*x2*x2 * (tmp - 2.) / (SQRT(tmp*tmp-4*x2*s2)**3)
  END FUNCTION integrand_rade_r


  FUNCTION integrand_rade_z(x,plist) RESULT(fx)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x
    REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
    REAL :: fx
    REAL :: s2,z2,x2, tmp

    x2 = x*x
    s2 = plist(1)*plist(1)
    z2 = plist(2)*plist(2)
    tmp = 1.0 + (z2+s2)*x2
    fx = (1.-SQRT(x)) * z2*x2*x2 * tmp / (SQRT(tmp*tmp-4*x2*s2)**3)
  END FUNCTION integrand_rade_z


END MODULE sources_diskthomson
