!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fluxes_midpoint.f90                                               #
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
! numerical fluxes module for midpoint quadrature rule
!----------------------------------------------------------------------------!
MODULE fluxes_midpoint
  USE fluxes_common
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_generic
  USE reconstruction_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! name for the numerical flux function
  CHARACTER(LEN=32), PARAMETER :: qrule_name = "midpoint"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Fluxes_TYP, &
       ! methods
       InitFluxes, &
       InitFluxes_midpoint, &
       GetType, &
       GetName, &
       GetRank, &
       Info, &
       Warning, &
       Error, &
       PrimRecon, &
       CalculateFaceData, &
       CalculateFluxes_midpoint
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitFluxes_midpoint(this,qrule)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP)  :: this
    INTEGER           :: qrule
    !------------------------------------------------------------------------!
    INTENT(IN)        :: qrule
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL InitFluxes(this,qrule,qrule_name)
  END SUBROUTINE InitFluxes_midpoint


  PURE SUBROUTINE CalculateFaceData(this,Mesh,Physics,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar,cvar
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,cvar
    INTENT(INOUT)     :: this,Physics
    !------------------------------------------------------------------------!

    ! reconstruct data on cell faces
    IF (PrimRecon(this%reconstruction)) THEN
       CALL CalculateSlopes(this%Reconstruction,Mesh,Physics,pvar)
       CALL CalculateStates(this%Reconstruction,Mesh,Physics,4,Mesh%bcenter,&
            Mesh%fpos,pvar,this%rstates)
       CALL Convert2Conservative(Physics,Mesh,this%rstates,this%cons)
    ELSE
       CALL CalculateSlopes(this%Reconstruction,Mesh,Physics,cvar)
       CALL CalculateStates(this%Reconstruction,Mesh,Physics,4,Mesh%bcenter,&
            Mesh%fpos,cvar,this%rstates)
       CALL Convert2Primitive(Physics,Mesh,this%rstates,this%prim)
    END IF

    ! get minimal & maximal wave speeds on cell interfaces
    CALL CalculateWaveSpeeds(Physics,Mesh,this%prim)
  END SUBROUTINE CalculateFaceData


  PURE SUBROUTINE CalculateFluxes_midpoint(this,Mesh,Physics,pvar,cvar,xflux,yflux)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar,cvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: xflux,yflux
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,cvar
    INTENT(INOUT)     :: this,Physics
    INTENT(OUT)       :: xflux,yflux
    !------------------------------------------------------------------------!

    ! execute generic tasks common to all flux types
    CALL CalculateFaceData(this,Mesh,Physics,pvar,cvar)

    ! west and east
    IF (Mesh%INUM.GT.1) THEN
       ! physical fluxes
       CALL CalculateFluxesX(Physics,Mesh,1,2,this%prim,this%cons,this%pfluxes)
       ! numerical fluxes
       FORALL (i=Mesh%IMIN-1:Mesh%IMAX,j=Mesh%JGMIN:Mesh%JGMAX,k=1:Physics%vnum)
          ! x-direction
          xflux(i,j,k) = Mesh%dAxdy(i+1,j,1) / &
            (Physics%amax(i,j) - Physics%amin(i,j)) * &
            (Physics%amax(i,j)*this%pfluxes(i,j,2,k) - &
            Physics%amin(i,j)*this%pfluxes(i+1,j,1,k) + &
            Physics%amin(i,j)*Physics%amax(i,j) * &
            (this%cons(i+1,j,1,k) - this%cons(i,j,2,k)))
       END FORALL
    ELSE
       xflux(:,:,:) = 0.0
    END IF

    ! south and north
    IF (Mesh%JNUM.GT.1) THEN
       ! physical fluxes
       CALL CalculateFluxesY(Physics,Mesh,3,4,this%prim,this%cons,this%pfluxes)
       ! numerical fluxes
       FORALL (i=Mesh%IGMIN:Mesh%IGMAX,j=Mesh%JMIN-1:Mesh%JMAX,k=1:Physics%vnum)
         ! y-direction (nonsense for i=IMIN-1)
         yflux(i,j,k) = Mesh%dAydx(i,j+1,1) / &
            (Physics%bmax(i,j) - Physics%bmin(i,j)) * &
            (Physics%bmax(i,j)*this%pfluxes(i,j,4,k) - &
            Physics%bmin(i,j)*this%pfluxes(i,j+1,3,k) + &
            Physics%bmin(i,j)*Physics%bmax(i,j) * &
            (this%cons(i,j+1,3,k) - this%cons(i,j,4,k)))
       END FORALL
    ELSE
       yflux(:,:,:) = 0.0
    END IF

  END SUBROUTINE CalculateFluxes_midpoint

END MODULE fluxes_midpoint
