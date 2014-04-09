!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fluxes_trapezoidal.f90                                            #
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
! numerical fluxes module for trapezoidal quadrature rule
!----------------------------------------------------------------------------!
MODULE fluxes_trapezoidal
  USE mesh_common, ONLY : Mesh_TYP
  USE fluxes_midpoint, CloseFluxes_trapezoidal => CloseFluxes_midpoint
  USE physics_generic
  USE reconstruction_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitFluxes_trapezoidal, &
       CalculateFluxes_trapezoidal, &
       BilinearInterpolation, &
       CloseFluxes_trapezoidal
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitFluxes_trapezoidal(this,Mesh,qrule)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: qrule
    !------------------------------------------------------------------------!
    INTEGER           :: err,n
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,qrule
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL InitFluxes(this,qrule,GetName(Mesh))
    ! allocate arrays used in fluxes_trapezoidal
    ALLOCATE(this%dx(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4), &
         this%dy(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4), &
         STAT=err)
    IF (err.NE.0) &
       CALL Error(this, "InitFluxes_midpoint","Unable to allocate memory.")
    ! set relative positions for reconstruction points:
    ! cell corner positions
    DO n=1,4
       this%dx(:,:,n) = Mesh%cpos(:,:,n,1) - Mesh%bcenter(:,:,1)
       this%dy(:,:,n) = Mesh%cpos(:,:,n,2) - Mesh%bcenter(:,:,2)
    END DO
  END SUBROUTINE InitFluxes_trapezoidal

  PURE SUBROUTINE CalculateFluxes_trapezoidal(this,Mesh,Physics,pvar,cvar,xflux,yflux)
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
    ! get slopes and wave speeds
    CALL CalculateFaceData(this,Mesh,Physics,pvar,cvar) 

    ! west and east
    IF (Mesh%INUM.GT.1) THEN
       ! physical fluxes
       CALL CalculateFluxesX(Physics,Mesh,1,4,this%prim,this%cons,this%pfluxes)
       ! numerical fluxes
       FORALL (i=Mesh%IMIN-1:Mesh%IMAX,j=Mesh%JGMIN:Mesh%JGMAX,k=1:Physics%vnum)
          ! x-direction
          xflux(i,j,k) = 0.5/(Physics%amax(i,j) - Physics%amin(i,j)) * ( &
               Mesh%dAxdy(i+1,j,1) * ( &
                 Physics%amax(i,j)*this%pfluxes(i,j,2,k) - &
                 Physics%amin(i,j)*this%pfluxes(i+1,j,1,k) + &
                 Physics%amin(i,j)*Physics%amax(i,j)* &
               (this%cons(i+1,j,1,k) - this%cons(i,j,2,k))) + &
               Mesh%dAxdy(i+1,j,2) * ( &
                 Physics%amax(i,j)*this%pfluxes(i,j,4,k) - &
                 Physics%amin(i,j)*this%pfluxes(i+1,j,3,k) + &
                 Physics%amin(i,j)*Physics%amax(i,j)* &
               (this%cons(i+1,j,3,k) - this%cons(i,j,4,k))))
       END FORALL
    ELSE
       xflux(:,:,:) = 0.0
    END IF

    ! south and north
    IF (Mesh%JNUM.GT.1) THEN
       ! physical fluxes
       CALL CalculateFluxesY(Physics,Mesh,1,4,this%prim,this%cons,this%qfluxes)
       ! numerical fluxes
       FORALL (i=Mesh%IGMIN:Mesh%IGMAX,j=Mesh%JMIN-1:Mesh%JMAX,k=1:Physics%vnum)
          ! y-direction
          yflux(i,j,k) = 0.5 / (Physics%bmax(i,j) - Physics%bmin(i,j)) * ( &
             Mesh%dAydx(i,j+1,1) * ( &
               Physics%bmax(i,j)*this%qfluxes(i,j,3,k) - &
               Physics%bmin(i,j)*this%qfluxes(i,j+1,1,k) + &
               Physics%bmin(i,j)*Physics%bmax(i,j)* &
             (this%cons(i,j+1,1,k)-this%cons(i,j,3,k))) + &
             Mesh%dAydx(i,j+1,2) * ( &
               Physics%bmax(i,j)*this%qfluxes(i,j,4,k) - &
               Physics%bmin(i,j)*this%qfluxes(i,j+1,2,k) + &
               Physics%bmin(i,j)*Physics%bmax(i,j)* &
             (this%cons(i,j+1,2,k)-this%cons(i,j,4,k))))
       END FORALL
    ELSE
       yflux(:,:,:) = 0.0
    END IF
  END SUBROUTINE CalculateFluxes_trapezoidal


  SUBROUTINE BilinearInterpolation(this,Mesh,Physics,rvar,rint)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP)   :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
         :: rvar,rint
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,rvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: rint
    !------------------------------------------------------------------------!
    ! interpolate primitive vars in upper right cell corner
    ! with respect to cell (i,j) for inner points
    FORALL (i=Mesh%IMIN:Mesh%IMAX-1,j=Mesh%JMIN:Mesh%JMAX-1,k=1:Physics%vnum)
       rint(i,j,k) = Mesh%weights(i,j,2,2) * rvar(i,j,k) &
            + Mesh%weights(i+1,j,1,2) * rvar(i+1,j,k) &
            + Mesh%weights(i,j+1,2,1) * rvar(i,j+1,k) &
            + Mesh%weights(i+1,j+1,1,1) * rvar(i+1,j+1,k)
    END FORALL

    ! simple extrapolation for boundary values
    ! south-west corner
    rint(Mesh%IMIN-1,Mesh%JMIN-1,:) = rvar(Mesh%IMIN,Mesh%JMIN,:)
    ! south-east corner
    rint(Mesh%IMAX,Mesh%JMIN-1,:)   = rvar(Mesh%IMAX,Mesh%JMIN,:)
    ! north-west corner
    rint(Mesh%IMIN-1,Mesh%JMAX,:)   = rvar(Mesh%IMIN,Mesh%JMAX,:)
    ! north-east corner
    rint(Mesh%IMAX,Mesh%JMAX,:)     = rvar(Mesh%IMAX,Mesh%JMAX,:)

    ! western boundary
    rint(Mesh%IMIN-1,Mesh%JMIN:Mesh%JMAX-1,:) = 0.5*(rvar(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX-1,:) &
         + rvar(Mesh%IMIN,Mesh%JMIN+1:Mesh%JMAX,:))
    ! eastern boundary
    rint(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX-1,:) = 0.5*(rvar(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX-1,:) &
         + rvar(Mesh%IMAX,Mesh%JMIN+1:Mesh%JMAX,:))
    ! southern boundary
    rint(Mesh%IMIN:Mesh%IMAX-1,Mesh%JMIN-1,:) = 0.5*(rvar(Mesh%IMIN:Mesh%IMAX-1,Mesh%JMIN,:) &
         + rvar(Mesh%IMIN+1:Mesh%IMAX,Mesh%JMIN,:))
    ! northern boundary
    rint(Mesh%IMIN:Mesh%IMAX-1,Mesh%JMAX,:) = 0.5*(rvar(Mesh%IMIN:Mesh%IMAX-1,Mesh%JMAX,:) &
         + rvar(Mesh%IMIN+1:Mesh%IMAX,Mesh%JMAX,:))
  END SUBROUTINE BilinearInterpolation
  

END MODULE fluxes_trapezoidal
