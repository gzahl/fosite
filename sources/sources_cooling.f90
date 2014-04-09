!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_cooling.f90                                               #
!#                                                                           #
!# Copyright (C) 2009,2011                                                   #
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
!! \brief source terms module for simple optically thin cooling
!!
!! \extends sources_c_accel
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_cooling
  USE sources_c_accel
  USE fluxes_common, ONLY : Fluxes_TYP
  USE physics_generic
  USE mesh_generic
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: source_name = "optically thin cooling"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! methods
       InitSources_cooling, &
       ExternalSources_cooling, &
       CalcTimestep_cooling, &
       CloseSources_cooling
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources_cooling(this,Mesh,Physics,config)
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
    ! Courant number, i.e. safety factor for numerical stability
    CALL RequireKey(config, "cvis", 0.9)
    CALL GetAttr(config, "cvis", this%cvis)
    ALLOCATE(this%Qcool(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         STAT=err)
    IF (err.NE.0) CALL Error(this,"InitSources_cooling","memory allocation failed")
    this%Qcool(:,:) = 0.0
    this%time = -1.0
  END SUBROUTINE InitSources_cooling


  PURE SUBROUTINE ExternalSources_cooling(this,Mesh,Physics,time,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                      :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,time,pvar,cvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! update cooling function
    CALL UpdateCooling(this,Mesh,Physics,time,pvar)
    ! compute energy source terms due to optically thin cooling
    IF (Physics%ENERGY.NE.0) THEN
       sterm(:,:,1:Physics%ENERGY-1) = 0.0
       sterm(:,:,Physics%ENERGY) = this%Qcool(:,:)
       sterm(:,:,Physics%ENERGY:Physics%VNUM) = 0.0
    ELSE
       sterm(:,:,:) = 0.0
    END IF
  END SUBROUTINE ExternalSources_cooling

 
  PURE SUBROUTINE CalcTimestep_cooling(this,Mesh,Physics,time,pvar,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                      :: pvar
    REAL              :: time,dt
    !------------------------------------------------------------------------!
    REAL              :: invdt
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,time,pvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: dt
    !------------------------------------------------------------------------!
    CALL UpdateCooling(this,Mesh,Physics,time,pvar)
    ! maximum of inverse cooling timescale t_cool ~ P/Q_cool
    invdt = MAXVAL(ABS(this%Qcool(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) &
         / pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE)))
    IF (invdt.GT.TINY(invdt)) THEN
       dt = this%cvis / invdt
    ELSE
       dt = HUGE(invdt)
    END IF
  END SUBROUTINE CalcTimestep_cooling


  PURE SUBROUTINE UpdateCooling(this,Mesh,Physics,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    REAL              :: N,T
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,time,pvar
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! update the cooling function if time has changed
    IF ((time.NE.this%time) .OR. (time .EQ. 0.)) THEN
!CDIR COLLAPSE
       DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
          DO i=Mesh%IGMIN,Mesh%IGMAX
             ! compute particle number N = rho/m and temperature
             ! with mean particle mass m = mu/NA,
             ! with mean molecular mass mu in kg/mol,
             ! with the Avogadro number NA = 6.022e23 [particles/mol] and
             N = pvar(i,j,Physics%DENSITY) * Physics%Constants%NA/Physics%mu
             T = Physics%mu*pvar(i,j,Physics%PRESSURE) &
                  / (Physics%Constants%RG*pvar(i,j,Physics%DENSITY))
             this%Qcool(i,j) = -N*N * lambda(T)
          END DO
       END DO
       this%time=time
    END IF
  END SUBROUTINE UpdateCooling

  SUBROUTINE CloseSources_cooling(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%Qcool)
    CALL CloseSources(this)
  END SUBROUTINE CloseSources_cooling

  ! simple optically thin cooling function
  ELEMENTAL FUNCTION lambda(T) RESULT(L)
    IMPLICIT NONE
    REAL, INTENT(IN) :: T
    REAL :: L
    !!! SI units !!!
    IF (T.LT.1.2e4) THEN
       ! disable cooling for T < 1.2e4 K, i.e. set L to a small value
       ! this is a very very rough approximation
       ! dont't trust this function for T < 1.2e4 !!!
       L = 1.0D-40
    ELSE IF ((T.GE.1.2e4).AND.(T.LT.2.0e5)) THEN
       L = 7.1D-38 * SQRT(T)
    ELSE IF ((T.GE.2.0e5).AND.(T.LT.5.0e7)) THEN
       L = 2.3D-32 * EXP(-0.54*LOG(T))
    ELSE
       L = 1.3D-40 * SQRT(T)
    END IF
  END FUNCTION lambda

END MODULE sources_cooling
