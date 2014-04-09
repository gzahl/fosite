!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_diskthomson.f90                                           #
!#                                                                           #
!# Copyright (C) 2007-2008,2011                                              #
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
!! \brief source terms module for radiational acceleration due to Thomson
!! scattering of accretion disk radiation
!!
!! \extends sources_c_accel
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_diskthomson
  USE sources_c_accel
  USE fluxes_common, ONLY : Fluxes_TYP
  USE physics_generic
  USE mesh_generic
  USE integration
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  REAL, PARAMETER :: EPS  = 1.0D-05              ! precision for integration !
  CHARACTER(LEN=32), PARAMETER :: source_name = "thomoson scat. of disk rad."
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! methods
       InitSources_diskthomson, &
       InfoSources_diskthomson, &
       ExternalSources_diskthomson, &
       CloseSources_diskthomson
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources_diskthomson(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Dict_TYP),POINTER :: config,IO
    INTEGER           :: stype
    !------------------------------------------------------------------------!
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM) :: accel
    REAL              :: params(2)
    REAL              :: Ldisk,r0,r1,rs,r0rs,x0,xm,factor
    INTEGER           :: err,valwrite
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "stype", stype)
    CALL InitSources(this,stype,source_name)

    ! central mass
    CALL RequireKey(config, "mass", 1.0)
    CALL GetAttr(config, "mass", this%mass)
    
    ! accretion rate
    CALL RequireKey(config, "mdot", 1.0)
    CALL GetAttr(config, "mdot", this%mdot)
    
    ! inner and outer disk radius
    CALL RequireKey(config, "rin", 1.0)
    CALL GetAttr(config, "rin", r0)
    
    CALL RequireKey(config, "rout", 2.0)
    CALL GetAttr(config, "rout", r1)

    ! some constants
    rs = 2*Physics%constants%GN * &          ! Schwarzschildradius of the BH
         (this%mass / (Physics%constants%C**2 + TINY(1.0)))
    r0rs = r0 / rs                           ! inner radius in terms of R_s
    Ldisk = 0.5*this%mdot*Physics%constants%GN*this%mass/r0  ! disk luminosity
    factor = 1.5*Ldisk / r0**2 &                       ! constant factor
          * Physics%constants%KE / Physics%constants%C

    ! reserve memory for radiational acceleration
    ALLOCATE(this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM), &
         STAT = err)
    IF (err.NE.0) CALL Error(this,"InitSources_diskthomson", "Unable allocate memory!")
   
    this%accel(:,:,:)  = 0.0

    ! initialize radiational acceleration
    x0 = r0/r1
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! calculate cartesian components of radiational acceleration
          ! due to Thomson scattering
          params(1) = ABS(Mesh%bccart(i,j,1) / (r0+TINY(1.0)))          ! = r / r0
          params(2) = ABS(Mesh%bccart(i,j,2) / (r0+TINY(1.0)))          ! = z / r0
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

    ! convert to curvilinear vector components
    CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,accel(:,:,1:2),this%accel(:,:,1:2))

    ! check if output is requested
    IF (HasKey(config, "output/accel")) CALL GetAttr(config, "output/accel", valwrite)
    IF (valwrite .EQ. 1) THEN
       CALL AddField(IO, &
               "accel", &
               this%accel, &
               Dict("name" / "radaccel"))
    END IF

  END SUBROUTINE InitSources_diskthomson


  SUBROUTINE InfoSources_diskthomson(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: mass_str,mdot_str
    !------------------------------------------------------------------------!
    WRITE (mass_str,'(ES8.2)') this%mass
    WRITE (mdot_str,'(ES8.2)') this%mdot
    CALL Info(this,"            mass:              " // TRIM(mass_str) // &
        ACHAR(10)//"            accretion rate:    " // TRIM(mdot_str))
  END SUBROUTINE InfoSources_diskthomson


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
    REAL :: r,r2,z,z2,x2,tmp

    r  = plist(1)
    z  = plist(2)
    x2 = x*x
    r2 = r*r
    z2 = z*z
    tmp = 1.0 + (z2+r2)*x2
    fx = (1.-SQRT(x)) * z*r*x2*x2 * (tmp - 2.) / (SQRT(tmp*tmp-4*x2*r2)**3)
  END FUNCTION integrand_rade_r


  FUNCTION integrand_rade_z(x,plist) RESULT(fx)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x
    REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
    REAL :: fx
    REAL :: r2,z2,x2, tmp

    x2 = x*x
    r2 = plist(1)*plist(1)
    z2 = plist(2)*plist(2)
    tmp = 1.0 + (z2+r2)*x2
    fx = (1.-SQRT(x)) * z2*x2*x2 * tmp / (SQRT(tmp*tmp-4*x2*r2)**3)
  END FUNCTION integrand_rade_z


  SUBROUTINE CloseSources_diskthomson(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%accel)
    CALL CloseSources(this)
  END SUBROUTINE CloseSources_diskthomson


END MODULE sources_diskthomson
