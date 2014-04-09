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
! source terms module for radiational acceleration due to
! Thomson scattering of accretion disk radiation
!----------------------------------------------------------------------------!
MODULE sources_diskthomson
  USE sources_pointmass
  USE fluxes_common, ONLY : Fluxes_TYP
  USE physics_generic
  USE mesh_generic
  USE integration
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  REAL, PARAMETER :: EPS  = 1.0D-04              ! precision for integration !
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
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: accel
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
         (this%mass / (Physics%constants%C**2 + TINY(1.0)))
    s0rs = s0 / rs                           ! inner radius in terms of R_s  !
    factor = 3.*Physics%constants%KE / &     ! constant factor               !
         (16*PI*s0rs**3 + TINY(1.0)) * (this%mdot/rs) * (Physics%constants%C/rs)

    ! reserve memory for radiational acceleration
    ALLOCATE(this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
         STAT = err)
    IF (err.NE.0) CALL Error(this,"InitSources_diskthomson", "Unable allocate memory!")

    ! initialize radiational acceleration
    x0 = s0/s1
    DO j=Mesh%JMIN,Mesh%JMAX
       DO i=Mesh%IMIN,Mesh%IMAX
          ! calculate cartesian components of radiational acceleration
          ! due to Thomson scattering
          params(1) = ABS(Mesh%bccart(i,j,1) / (s0+TINY(1.0)))          ! = r / s0
          params(2) = ABS(Mesh%bccart(i,j,2) / (s0+TINY(1.0)))          ! = z / s0
          xm = MIN(1./(params(1)+TINY(1.0)),1.0)
          ! integrate around the pseudo-singularity at xm 
          accel(i,j,1) = SIGN(1.0,Mesh%bccart(i,j,2)) * factor * (&
               integrate(integrand_rade_z,x0,xm,EPS,params,method=1) &
             + integrate(integrand_rade_z,xm,1.0,EPS,params,method=1))
          accel(i,j,2) = factor * (&
               integrate(integrand_rade_r,x0,xm,EPS,params,method=1) &
             + integrate(integrand_rade_r,xm,1.0,EPS,params,method=1))
       END DO
    END DO
    accel(Mesh%IGMIN:Mesh%IMIN-1,:,:) = 0.0
    accel(Mesh%IMAX+1:Mesh%IGMAX,:,:) = 0.0
    accel(:,Mesh%JGMIN:Mesh%JMIN-1,:) = 0.0
    accel(:,Mesh%JMAX+1:Mesh%JGMAX,:) = 0.0

    ! convert to curvilinear vector components
    CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,accel,this%accel)
!!$    PRINT *, SQRT(MIN(MINVAL(ABS(Mesh%dlx(:,:) / this%accel(:,:,1))), &
!!$         MINVAL(ABS(Mesh%dly(:,:) / this%accel(:,:,2)))))
!!$    DO i=Mesh%IMIN,Mesh%IMAX
!!$       DO j=Mesh%JMIN,Mesh%JMAX
!!$          PRINT '(4(ES14.6))', Mesh%bccart(i,j,:),accel(i,j,:) 
!!$       END DO
!!$       PRINT '(A)', ""
!!$    END DO
  END SUBROUTINE InitSources_diskthomson


  PURE SUBROUTINE ExternalSources_diskthomson(this,Mesh,Physics,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: cvar,pvar,sterm
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    REAL, DIMENSION(Physics%VNUM) :: bflux
    REAL              :: oldmass
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,cvar,pvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! gravitational source terms
    CALL ExternalSources(Physics,Mesh,this%accel,pvar,cvar,sterm)
  END SUBROUTINE ExternalSources_diskthomson
  
 
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


  SUBROUTINE CloseSources_diskthomson(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%accel)
  END SUBROUTINE CloseSources_diskthomson


END MODULE sources_diskthomson
