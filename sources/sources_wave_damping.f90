!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_wave_damping.f90                                          #
!#                                                                           #
!# Copyright (C) 2012                                                        #
!# Manuel Jung      <mjung@astrophysik.uni-kiel.de>
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
!> \author Manuel Jung
!!
!! \brief source terms module for wave damping as required by the planet eu
!! experiment
!!
!! \extends sources_c_accel
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_wave_damping
  USE sources_c_accel
  USE fluxes_common, ONLY : Fluxes_TYP
  USE physics_generic
  USE mesh_generic
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: source_name = "wave damping"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! methods
       InitSources_wave_damping, &
       ExternalSources_wave_damping, &
       CloseSources_wave_damping
  !--------------------------------------------------------------------------!

CONTAINS


  SUBROUTINE InitSources_wave_damping(this,Mesh,Physics,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Physics_TYP) :: Physics
    TYPE(Dict_TYP),POINTER :: config
    INTEGER           :: stype
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "stype", stype)
    CALL InitSources(this,stype,source_name)
 
    ALLOCATE(this%init_pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,&
                            Physics%VNUM), &
         STAT = err)

    IF (err.NE.0) &
         CALL Error(this,"InitSources_wave_damping","memory allocation failed")

    CALL RequireKey(config, "r0")
    CALL RequireKey(config, "r1")
    CALL RequireKey(config, "tau0")
    CALL RequireKey(config, "tau1")

    CALL GetAttr(config, "r0", this%r(1))
    CALL GetAttr(config, "r1", this%r(2))
    CALL GetAttr(config, "tau0", this%tau(1))
    CALL GetAttr(config, "tau1", this%tau(2))
!    print *,"r: ", this%r
!    print *,"tau: ", this%tau
!    print *,"Rinnen: ", Mesh%xmin
!    print *,"Raussen: ", Mesh%xmax

  END SUBROUTINE InitSources_wave_damping


  SUBROUTINE ExternalSources_wave_damping(this,Mesh,Physics,time,dt,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL              :: time,dt
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: cvar,pvar,sterm
    INTEGER           :: i,j,k,m
    REAL              :: b, r, radius
    REAL, DIMENSION(2):: rbound
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,cvar,pvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!

    IF(time.LE.0.) THEN
        this%init_pvar(:,:,:) = cvar(:,:,:)
    END IF

    rbound(1) = Mesh%xmin
    rbound(2) = Mesh%xmax 
    sterm(:,:,:) = 0.0
    DO k=1,2
        r = this%r(k)
        b = 1./(this%tau(k) * (rbound(k)-r)**2)
        DO j=Mesh%JMIN,Mesh%JMAX
            DO i=Mesh%IMIN,Mesh%IMAX
                radius = SQRT(Mesh%bccart(i,j,1)**2+Mesh%bccart(i,j,2)**2)
                IF(((k==1).AND.(radius.LE.r))&
                   .OR.((k==2).AND.(radius.GE.r))) THEN
                    sterm(i,j,Physics%DENSITY) = &
                        -(cvar(i,j,Physics%DENSITY) &
                          - this%init_pvar(i,j,Physics%DENSITY)) &
                         * (radius-r)**2 * b

                    !sterm(i,j,Physics%XMOMENTUM) = &
                    !    -(2.*cvar(i,j,Physics%XMOMENTUM) &
                    !      -pvar(i,j,Physics%DENSITY)*this%init_pvar(i,j,Physics%XVELOCITY) &
                    !      -this%init_pvar(i,j,Physics%DENSITY)*pvar(i,j,Physics%XVELOCITY)) &
                    !     * (radius-r)**2 * b

                    !sterm(i,j,Physics%XMOMENTUM) = &
                    !    -(cvar(i,j,Physics%XMOMENTUM) &
                    !      -pvar(i,j,Physics%DENSITY)*this%init_pvar(i,j,Physics%XVELOCITY)) &
                    !     * (radius-r)**2 * b

                    !sterm(i,j,Physics%XMOMENTUM) = &
                    !    cvar(i,j,Physics%XMOMENTUM)/pvar(i,j,Physics%DENSITY)*sterm(i,j,Physics%DENSITY)
                    sterm(i,j,Physics%YMOMENTUM) = &
                        cvar(i,j,Physics%YMOMENTUM)/pvar(i,j,Physics%DENSITY)*sterm(i,j,Physics%DENSITY)

!                    sterm(i,j,Physics%YMOMENTUM) = &
!                        cvar(i,j,Physics%YMOMENTUM)/pvar(i,j,Physics%DENSITY)*sterm(i,j,Physics%DENSITY) &
!                        -(pvar(i,j,Physics%DENSITY)*Mesh%bhy(i,j)&
!                          *(pvar(i,j,Physics%YVELOCITY)-this%init_pvar(i,j,Physics%YVELOCITY)))&
!                         * (radius-r)**2 * b

                    !sterm(i,j,Physics%YMOMENTUM) = &
                    !    -(2.*pvar(i,j,Physics%DENSITY)*pvar(i,j,Physics%YVELOCITY)&
                    !      -pvar(i,j,Physics%DENSITY)*this%init_pvar(i,j,Physics%YVELOCITY)&
                    !      -this%init_pvar(i,j,Physics%DENSITY)*pvar(i,j,Physics%YVELOCITY))&
                    !     / tau &
                    !     * (ABS(radius-r)/ABS(rbound(k)-r))**2
                    DO m=Physics%XMOMENTUM,Physics%YMOMENTUM
                            sterm(i,j,m) = &
                                - (cvar(i,j,m)-this%init_pvar(i,j,m)) &
                                * (radius-r)**2 * b
                    END DO
                    IF(Physics%ENERGY.GT.0) THEN
                        sterm(i,j,Physics%ENERGY) = &
                            -(cvar(i,j,Physics%ENERGY) &
                              - this%init_pvar(i,j,Physics%ENERGY)) &
                             * (radius-r)**2 * b
                    END IF
                END IF
            END DO
        END DO
    END DO
    !CALL ExternalSources(Physics,Mesh,this%accel,pvar,cvar,sterm)
  END SUBROUTINE ExternalSources_wave_damping

 
  SUBROUTINE CloseSources_wave_damping(this,Fluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Fluxes_TYP)  :: Fluxes
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Fluxes
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%init_pvar)
    CALL CloseSources(this)
  END SUBROUTINE CloseSources_wave_damping

END MODULE sources_wave_damping
