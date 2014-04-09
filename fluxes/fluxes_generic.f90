!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fluxes_generic.f90                                                #
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
! generic module for numerical flux functions
!----------------------------------------------------------------------------!
MODULE fluxes_generic
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  USE fluxes_midpoint, InitFluxes_all => InitFluxes
  USE fluxes_trapezoidal
  USE reconstruction_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: MIDPOINT     = 1
  INTEGER, PARAMETER :: TRAPEZOIDAL  = 2
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Fluxes_TYP, &
       ! constants
       MIDPOINT, TRAPEZOIDAL, &
       ! methods
       InitFluxes, &
       MallocFluxes, &
       GetType, &
       GetName, &
       GetRank, &
       Info, &
       Warning, &
       Error, &
       PrimRecon, &
       CalculateFluxes, &
       CloseFluxes
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitFluxes(this,scheme)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP)  :: this
    INTEGER           :: scheme
    !------------------------------------------------------------------------!
    INTENT(IN)        :: scheme
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!

    ! set flux function properties
    SELECT CASE(scheme)
    CASE(MIDPOINT)
       CALL InitFluxes_midpoint(this,scheme)
    CASE(TRAPEZOIDAL)
       CALL InitFluxes_trapezoidal(this,scheme)
    CASE DEFAULT
       CALL Error(this, "InitFluxes", "Unknown flux type.")
    END SELECT

    ! print some information
    CALL Info(this, " FLUXES---> quadrature rule    " // TRIM(GetName(this)))
  END SUBROUTINE InitFluxes
  

  SUBROUTINE MallocFluxes(this,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! allocate memory for reconstruction object
    CALL MallocReconstruction(this%Reconstruction,Mesh,Physics)

    ! allocate memory for all arrays used in fluxes
    ALLOCATE(this%cons(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,Physics%vnum), &
         this%prim(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,Physics%vnum), &
         this%pfluxes(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,Physics%vnum), &
         this%qfluxes(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,Physics%vnum), &
         STAT = err)
    IF (err.NE.0) THEN
       CALL Error(this, "MallocFluxes", "Unable to allocate memory.")
    END IF

    ! set reconstruction pointer
    IF (PrimRecon(this%Reconstruction)) THEN
       this%rstates => this%prim
    ELSE
       this%rstates => this%cons
    END IF
  END SUBROUTINE MallocFluxes


  PURE SUBROUTINE CalculateFluxes(this,Mesh,Physics,pvar,cvar,xflux,yflux)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP)   :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) :: pvar,cvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) :: xflux,yflux
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,cvar
    INTENT(INOUT)     :: this,Physics
    INTENT(OUT)       :: xflux,yflux
    !------------------------------------------------------------------------!
    ! calculate numerical fluxes depending on the integration rule
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
    DEALLOCATE(this%cons,this%prim,this%pfluxes,this%qfluxes)
    CALL CloseReconstruction(this%Reconstruction)
  END SUBROUTINE CloseFluxes


END MODULE fluxes_generic
