!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_cooling.f90                                               #
!#                                                                           #
!# Copyright (C) 2009                                                        #
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
! source terms module for simple optically thin cooling
!----------------------------------------------------------------------------!
MODULE sources_cooling
  USE sources_pointmass
  USE fluxes_common, ONLY : Fluxes_TYP
  USE physics_generic
  USE mesh_generic
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
       GetType, &
       GetName
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources_cooling(this,Mesh,Physics,stype,cvis)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Physics_TYP) :: Physics
    INTEGER           :: stype
    REAL              :: cvis
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,stype,cvis
    !------------------------------------------------------------------------!
    CALL InitSources(this,stype,source_name)
    this%cvis    = cvis
  END SUBROUTINE InitSources_cooling


  PURE SUBROUTINE ExternalSources_cooling(this,Mesh,Physics,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    REAL              :: N,T
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,pvar,cvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    sterm(:,:,Physics%DENSITY) = 0.0
    sterm(:,:,Physics%XMOMENTUM) = 0.0
    sterm(:,:,Physics%YMOMENTUM) = 0.0
    IF (Physics%ZMOMENTUM.NE.0) sterm(:,:,Physics%ZMOMENTUM) = 0.0

    ! enery loss due to radiation processes: -N**2 * Lambda(T) [J/m^3/s]
!CDIR NOUNROLL
     DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
       DO i=Mesh%IMIN,Mesh%IMAX
          ! compute particle number N = rho/m and temperature
          ! with mean particle mass m = mu/NA,
          ! with mean molecular mass mu in kg/mol,
          ! with the Avogadro number NA = 6.022e23 [particles/mol] and
          N = pvar(i,j,Physics%DENSITY) * Physics%Constants%NA/Physics%mu
          T = Physics%mu*pvar(i,j,Physics%DENSITY) &
              / (Physics%Constants%RG*pvar(i,j,Physics%DENSITY))
          sterm(i,j,Physics%ENERGY) = -N*N * lambda(T)
       END DO
    END DO
  END SUBROUTINE ExternalSources_cooling

 
  SUBROUTINE CalcTimestep_cooling(this,Mesh,Physics,pvar,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                      :: pvar
    REAL              :: dt
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    INTENT(INOUT)     :: this,Physics
    INTENT(OUT)       :: dt
    !------------------------------------------------------------------------!
    ! compute temperature and particle number using Physics%tmin as temporary storage
    Physics%tmin(:,:,1) = Physics%mu*pvar(:,:,Physics%PRESSURE) &
       / (Physics%Constants%RG*pvar(:,:,Physics%DENSITY))
    Physics%tmin(:,:,2) = pvar(:,:,Physics%DENSITY) * Physics%Constants%NA/Physics%mu

    ! cooling timescale t_cool = P/(rho*n^2)
    dt = this%cvis * MINVAL( &
         pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) &
       / (Physics%tmin(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2) &
         *Physics%tmin(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2) & 
         *lambda(Physics%tmin(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1))))
  END SUBROUTINE CalcTimestep_cooling


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
