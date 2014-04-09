!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fluxes_generic.f90                                                #
!#                                                                           #
!# Copyright (C) 2007-2012                                                   #
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
! generic module for numerical flux functions
!----------------------------------------------------------------------------!
MODULE fluxes_generic
  USE physics_common, ONLY : Physics_TYP
  USE fluxes_midpoint, InitFluxes_common => InitFluxes, &
       CloseFluxes_common => CloseFluxes
  USE fluxes_trapezoidal
  USE mesh_generic
  USE reconstruction_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Fluxes_TYP, &
       ! constants
       MIDPOINT, TRAPEZOIDAL, &
       ! methods
       InitFluxes, &
       CloseFluxes, &
       CalculateFluxes, &
       PrimRecon, &
       GetBoundaryFlux, &
       GetType, &
       GetName, &
       GetRank, &
       Initialized, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitFluxes(this,Mesh,Physics,order,variables,limiter,theta)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    INTEGER, OPTIONAL :: order
    LOGICAL, OPTIONAL :: variables
    INTEGER, OPTIONAL :: limiter
    REAL, OPTIONAL    :: theta
    !------------------------------------------------------------------------!
    INTEGER                  :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,order,variables,limiter,theta
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! check initialization of Mesh and Physics
    IF (.NOT.Initialized(Mesh).OR..NOT.Initialized(Physics)) &
         CALL Error(this,"InitFluxes","mesh and/or physics module uninitialized")

    ! call specific flux initialization routines
    ! flux module type depends on mesh module type, see mesh_generic
!CDIR IEXPAND
    SELECT CASE(GetType(Mesh))
    CASE(MIDPOINT)
       CALL InitFluxes_midpoint(this,Mesh,MIDPOINT)
    CASE(TRAPEZOIDAL)
       CALL InitFluxes_trapezoidal(this,Mesh,TRAPEZOIDAL)
    CASE DEFAULT
       CALL Error(this, "InitFluxes", "Unknown mesh type.")
    END SELECT

    ! allocate memory for all arrays used in fluxes
    ALLOCATE(this%cons(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,Physics%VNUM), &
         this%prim(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,Physics%VNUM), &
         this%pfluxes(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,Physics%VNUM), &
         this%qfluxes(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,Physics%VNUM), &
         this%bxflux(Mesh%JGMIN:Mesh%JGMAX,2,Physics%VNUM), &
         this%byflux(Mesh%IGMIN:Mesh%IGMAX,2,Physics%VNUM), &
         this%bxfold(Mesh%JGMIN:Mesh%JGMAX,2,Physics%VNUM), &
         this%byfold(Mesh%IGMIN:Mesh%IGMAX,2,Physics%VNUM), &
         STAT = err)
    IF (err.NE.0) THEN
       CALL Error(this, "InitFluxes", "Unable to allocate memory.")
    END IF

    ! print some information
    CALL Info(this, " FLUXES---> quadrature rule    " // TRIM(GetName(this)))

    ! initialize reconstruction modules
    CALL InitReconstruction(this%reconstruction,Mesh,Physics,order,variables,limiter,theta)

    ! set reconstruction pointer
!CDIR IEXPAND
    IF (PrimRecon(this%Reconstruction)) THEN
       this%rstates => this%prim
    ELSE
       this%rstates => this%cons
    END IF

    ! initialize boundary fluxes
    this%bxflux(:,:,:) = 0.
    this%byflux(:,:,:) = 0.
    this%bxfold(:,:,:) = 0.
    this%byfold(:,:,:) = 0.
  END SUBROUTINE InitFluxes
  

  PURE SUBROUTINE CalculateFluxes(this,Mesh,Physics,pvar,cvar,xflux,yflux)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP)   :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) :: pvar,cvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) :: xflux,yflux
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,cvar
    INTENT(INOUT)     :: this,Physics
    INTENT(OUT)       :: xflux,yflux
    !------------------------------------------------------------------------!
    ! calculate numerical fluxes depending on the integration rule
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(MIDPOINT)
       CALL CalculateFluxes_midpoint(this,Mesh,Physics,pvar,cvar,xflux,yflux)
    CASE(TRAPEZOIDAL)
       CALL CalculateFluxes_trapezoidal(this,Mesh,Physics,pvar,cvar,xflux,yflux)
    END SELECT
  END SUBROUTINE CalculateFluxes

  
  SUBROUTINE CloseFluxes(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP)  :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    IF (.NOT.Initialized(this)) &
        CALL Error(this,"CloseFluxes","not initialized")
    DEALLOCATE(this%cons,this%prim,this%pfluxes,this%qfluxes, &
         this%bxflux,this%byflux,this%bxfold,this%byfold)
    SELECT CASE(GetType(this))
    CASE(MIDPOINT)
       CALL CloseFluxes_midpoint(this)
    CASE(TRAPEZOIDAL)
       CALL CloseFluxes_trapezoidal(this)
    END SELECT
    CALL CloseReconstruction(this%Reconstruction)
    CALL CloseFluxes_common(this)
  END SUBROUTINE CloseFluxes


END MODULE fluxes_generic
