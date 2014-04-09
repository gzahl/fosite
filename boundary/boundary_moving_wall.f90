!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_moving_wall.f90                                          #
!#                                                                           #
!# Copyright (C) 2006-2009                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Bjoern Sperling  <sperling@astrophysik.uni-kiel.de>                       #
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
! boundary module for moving wall 
!----------------------------------------------------------------------------!
MODULE boundary_moving_wall
  USE mesh_common, ONLY : Mesh_TYP
  USE fluxes_common
  USE reconstruction_common
  USE boundary_nogradients
  USE physics_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "moving wall"  
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Boundary_TYP, &
       ! constants
       WEST, EAST, SOUTH, NORTH, &
       ! methods
       InitBoundary, &
       InitBoundary_moving_wall, &
       CenterBoundary_moving_wall, &
       FaceBoundary_moving_wall, &
       CloseBoundary_moving_wall, &
       GetType, &
       GetName, &
       GetDirection, &
       GetDirectionName, &
       GetRank, &
       GetNumProcs, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitBoundary_moving_wall(this,Mesh,Physics,btype,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: btype,dir
    !------------------------------------------------------------------------!
    INTEGER       :: err
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,btype,dir
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL InitBoundary(this,btype,boundcond_name,dir)
    
    ALLOCATE(this%reflX(Physics%vnum), &
         this%reflY(Physics%vnum), &
         STAT=err)

    IF (err.NE.0) THEN
       CALL Error(this, "InitBoundary_moving_wall", "Unable to allocate memory.")
    END IF

    ! allocate memory for boundary data
    SELECT CASE(dir)
    CASE(WEST,EAST)
       ALLOCATE(this%data(1,Mesh%JGMIN:Mesh%JGMAX,Physics%XVELOCITY:Physics%vnum), &
            STAT=err)
       this%data(:,:,:)=0.
    CASE(SOUTH,NORTH)
       ALLOCATE(this%data(Mesh%IGMIN:Mesh%IGMAX,1,Physics%XVELOCITY:Physics%vnum), &
            STAT=err)
       this%data(:,:,:)=0.
    END SELECT
    
    IF (err.NE.0) THEN
       CALL Error(this, "InitBoundary_moving_wall", "Unable to allocate memory.")
    END IF
    ! this tells us which vars get the opposite sign/vanish at cell faces;
    ! e.g. vertical velocities (depends on the underlying physics)
    CALL ReflectionMasks(Physics,this%reflX,this%reflY)
  END SUBROUTINE InitBoundary_moving_wall

  PURE SUBROUTINE CenterBoundary_moving_wall(this,Mesh,Fluxes,Physics,rvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Fluxes_TYP)   :: Fluxes
    TYPE(Physics_TYP)  :: Physics
    REAL :: rvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)
    !------------------------------------------------------------------------!
    INTEGER       :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,Mesh,Fluxes,Physics
    INTENT(INOUT) :: rvar   
    !------------------------------------------------------------------------!
    SELECT CASE(GetDirection(this))
    CASE(WEST)
       FORALL (j=Mesh%JGMIN:Mesh%JGMAX)
          WHERE (this%reflX)
             rvar(Mesh%IMIN-1,j,:) = -rvar(Mesh%IMIN,j,:)
             rvar(Mesh%IMIN-2,j,:) = -rvar(Mesh%IMIN+1,j,:)
          ELSEWHERE
             rvar(Mesh%IMIN-1,j,:) = rvar(Mesh%IMIN,j,:)
             rvar(Mesh%IMIN-2,j,:) = rvar(Mesh%IMIN+1,j,:)
          END WHERE
       END FORALL
     IF (PrimRecon(Fluxes%reconstruction).EQV.PRIMITIVE) THEN
     !set tangential velocity of boundary (y,z)!
       IF (Physics%YVELOCITY /= 0) THEN
          rvar(Mesh%IMIN-1,:,Physics%YVELOCITY) = &
             this%data(1,:,Physics%YVELOCITY)
          rvar(Mesh%IMIN-2,:,Physics%YVELOCITY) = &
             this%data(1,:,Physics%YVELOCITY)
       END IF
       IF (Physics%ZVELOCITY /= 0) THEN
             rvar(Mesh%IMIN-1,:,Physics%ZVELOCITY) = &
                 this%data(1,:,Physics%ZVELOCITY)
             rvar(Mesh%IMIN-2,:,Physics%ZVELOCITY) = &
                 this%data(1,:,Physics%ZVELOCITY)
       END IF
       !set temperature (if pressure exit and temp != 0)
       !here PRESSURE == temperature
       IF (Physics%PRESSURE /= 0) THEN
          WHERE (this%data(1,:,Physics%PRESSURE) > 0)
             rvar(Mesh%IMIN-1,:,Physics%PRESSURE) = &
                 this%data(1,:,Physics%PRESSURE) * &
                 rvar(Mesh%IMIN-1,:,Physics%DENSITY) &
                 * Physics%constants%RG / Physics%mu
             rvar(Mesh%IMIN-2,:,Physics%PRESSURE) = &
                 this%data(1,:,Physics%PRESSURE) * &
                 rvar(Mesh%IMIN-2,:,Physics%DENSITY) &
                 * Physics%constants%RG / Physics%mu
          END WHERE
       END IF
     ELSE
      IF (Physics%YMOMENTUM /= 0) THEN
          rvar(Mesh%IMIN-1,:,Physics%YMOMENTUM) = &
             this%data(1,:,Physics%YVELOCITY)&
             *rvar(Mesh%IMIN-1,:,Physics%DENSITY)
          rvar(Mesh%IMIN-2,:,Physics%YMOMENTUM) = &
             this%data(1,:,Physics%YVELOCITY)&
             *rvar(Mesh%IMIN-2,:,Physics%DENSITY)
       END IF
       IF (Physics%ZMOMENTUM /= 0) THEN
             rvar(Mesh%IMIN-1,:,Physics%ZMOMENTUM) = &
                 this%data(1,:,Physics%ZVELOCITY)&
                 *rvar(Mesh%IMIN-1,:,Physics%DENSITY)
             rvar(Mesh%IMIN-2,:,Physics%ZMOMENTUM) = &
                 this%data(1,:,Physics%ZVELOCITY)&
                 *rvar(Mesh%IMIN-2,:,Physics%DENSITY)
       END IF
       !set temperature (if pressure exit and temp != 0)
       !here PRESSURE == temperature
       IF (Physics%PRESSURE /= 0) THEN
          WHERE (this%data(1,:,Physics%PRESSURE) > 0)
             rvar(Mesh%IMIN-1,:,Physics%ENERGY) = &
                 this%data(1,:,Physics%PRESSURE) &
                 * Physics%constants%RG / (Physics%mu*(Physics%gamma-1))&
                 * rvar(Mesh%IMIN-1,:,Physics%DENSITY)
             rvar(Mesh%IMIN-2,:,Physics%ENERGY) = &
                 this%data(1,:,Physics%PRESSURE) &
                 * Physics%constants%RG / (Physics%mu*(Physics%gamma-1))&
                 * rvar(Mesh%IMIN-2,:,Physics%DENSITY)
          END WHERE
       END IF
     END IF
    CASE(EAST)
       FORALL (j=Mesh%JGMIN:Mesh%JGMAX)
          WHERE (this%reflX)
             rvar(Mesh%IMAX+1,j,:) = -rvar(Mesh%IMAX,j,:)
             rvar(Mesh%IMAX+2,j,:) = -rvar(Mesh%IMAX-1,j,:)
          ELSEWHERE
             rvar(Mesh%IMAX+1,j,:) = rvar(Mesh%IMAX,j,:)
             rvar(Mesh%IMAX+2,j,:) = rvar(Mesh%IMAX-1,j,:)
          END WHERE
       END FORALL
    IF (PrimRecon(Fluxes%reconstruction).EQV.PRIMITIVE) THEN
       !set tangential velocity of boundary (y,z)!
       IF (Physics%YVELOCITY /= 0) THEN
          rvar(Mesh%IMAX+1,:,Physics%YVELOCITY) = &
             this%data(1,:,Physics%YVELOCITY)
          rvar(Mesh%IMAX+2,:,Physics%YVELOCITY) = &
             this%data(1,:,Physics%YVELOCITY)
       END IF
       IF (Physics%ZVELOCITY /= 0) THEN
             rvar(Mesh%IMAX+1,:,Physics%ZVELOCITY) = &
                 this%data(1,:,Physics%ZVELOCITY)
             rvar(Mesh%IMAX+2,:,Physics%ZVELOCITY) = &
                 this%data(1,:,Physics%ZVELOCITY)
       END IF
       !set temperature (if pressure exit and temp != 0)
       !here PRESSURE == temperature
       IF (Physics%PRESSURE /= 0) THEN
          WHERE (this%data(1,:,Physics%PRESSURE) > 0)
             rvar(Mesh%IMAX+1,:,Physics%PRESSURE) = &
                 this%data(1,:,Physics%PRESSURE) * &
                 rvar(Mesh%IMAX+1,:,Physics%DENSITY) &
                 * Physics%constants%RG / Physics%mu
             rvar(Mesh%IMAX+2,:,Physics%PRESSURE) = &
                 this%data(1,:,Physics%PRESSURE) * &
                 rvar(Mesh%IMAX+2,:,Physics%DENSITY) &
                 * Physics%constants%RG / Physics%mu
          END WHERE
       END IF
     ELSE
      IF (Physics%YMOMENTUM /= 0) THEN
          rvar(Mesh%IMAX+1,:,Physics%YMOMENTUM) = &
             this%data(1,:,Physics%YVELOCITY)&
             *rvar(Mesh%IMAX+1,:,Physics%DENSITY)
          rvar(Mesh%IMAX+2,:,Physics%YMOMENTUM) = &
             this%data(1,:,Physics%YVELOCITY)&
             *rvar(Mesh%IMAX+2,:,Physics%DENSITY)
       END IF
       IF (Physics%ZMOMENTUM /= 0) THEN
             rvar(Mesh%IMAX+1,:,Physics%ZMOMENTUM) = &
                 this%data(1,:,Physics%ZVELOCITY)&
                 *rvar(Mesh%IMAX+1,:,Physics%DENSITY)
             rvar(Mesh%IMAX+2,:,Physics%ZMOMENTUM) = &
                 this%data(1,:,Physics%ZVELOCITY)&
                 *rvar(Mesh%IMAX+2,:,Physics%DENSITY)
       END IF
       !set temperature (if pressure exit and temp != 0)
       !here PRESSURE == temperature
       IF (Physics%PRESSURE /= 0) THEN
          WHERE (this%data(1,:,Physics%PRESSURE) > 0)
             rvar(Mesh%IMAX+1,:,Physics%ENERGY) = &
                 this%data(1,:,Physics%PRESSURE) &
                 * Physics%constants%RG / (Physics%mu*(Physics%gamma-1))&
                 * rvar(Mesh%IMAX+1,:,Physics%DENSITY)
             rvar(Mesh%IMAX+2,:,Physics%ENERGY) = &
                 this%data(1,:,Physics%PRESSURE) &
                 * Physics%constants%RG / (Physics%mu*(Physics%gamma-1))&
                 * rvar(Mesh%IMAX+2,:,Physics%DENSITY)
          END WHERE
       END IF
     END IF
   CASE(SOUTH)
       FORALL (i=Mesh%IGMIN:Mesh%IGMAX)
          WHERE (this%reflY)
             rvar(i,Mesh%JMIN-1,:) = -rvar(i,Mesh%JMIN,:)
             rvar(i,Mesh%JMIN-2,:) = -rvar(i,Mesh%JMIN+1,:)
          ELSEWHERE
             rvar(i,Mesh%JMIN-1,:) = rvar(i,Mesh%JMIN,:)
             rvar(i,Mesh%JMIN-2,:) = rvar(i,Mesh%JMIN+1,:)
          END WHERE
        END FORALL
      IF (PrimRecon(Fluxes%reconstruction).EQV.PRIMITIVE) THEN
      !set tangential velocity of boundary (x,z)!
        IF (Physics%XVELOCITY /= 0) THEN
           rvar(:,Mesh%JMIN-1,Physics%XVELOCITY) = &
              this%data(:,1,Physics%XVELOCITY)
           rvar(:,Mesh%JMIN-2,Physics%XVELOCITY) = &
              this%data(:,1,Physics%XVELOCITY)
        END IF
        IF (Physics%ZVELOCITY /= 0) THEN
           rvar(:,Mesh%JMIN-1,Physics%ZVELOCITY) = &
              this%data(:,1,Physics%ZVELOCITY)
           rvar(:,Mesh%JMIN-2,Physics%ZVELOCITY) = &
              this%data(:,1,Physics%ZVELOCITY)
        END IF
        !set temperature (if pressure exit and temp != 0)
        IF (Physics%PRESSURE /= 0) THEN
           WHERE (this%data(:,1,Physics%PRESSURE) > 0)
              rvar(:,Mesh%JMIN-1,Physics%PRESSURE) = &
                 this%data(:,1,Physics%PRESSURE) * &
                 rvar(:,Mesh%JMIN-1,Physics%DENSITY) &
                 * Physics%constants%RG / Physics%mu
              rvar(:,Mesh%JMIN-2,Physics%PRESSURE) = &
                 this%data(:,1,Physics%PRESSURE) * &
                 rvar(:,Mesh%JMIN-2,Physics%DENSITY) &
                 * Physics%constants%RG / Physics%mu
           END WHERE
        END IF
     ELSE
        IF (Physics%XVELOCITY /= 0) THEN
           rvar(:,Mesh%JMIN-1,Physics%XMOMENTUM) = &
              this%data(:,1,Physics%XVELOCITY)&
              *rvar(:,Mesh%JMIN-1,Physics%DENSITY)
           rvar(:,Mesh%JMIN-2,Physics%XMOMENTUM) = &
              this%data(:,1,Physics%XVELOCITY)&
              *rvar(:,Mesh%JMIN-2,Physics%DENSITY)
        END IF
        IF (Physics%ZVELOCITY /= 0) THEN
           rvar(:,Mesh%JMIN-1,Physics%ZVELOCITY) = &
              this%data(:,1,Physics%ZMOMENTUM)&
              *rvar(:,Mesh%JMIN-1,Physics%DENSITY)
           rvar(:,Mesh%JMIN-2,Physics%ZMOMENTUM) = &
              this%data(:,1,Physics%ZVELOCITY)&
              *rvar(:,Mesh%JMIN-2,Physics%DENSITY)
        END IF
        !set temperature (if pressure exit and temp != 0)
        IF (Physics%PRESSURE /= 0) THEN
           WHERE (this%data(:,1,Physics%PRESSURE) > 0)
              rvar(:,Mesh%JMIN-1,Physics%ENERGY) = &
                 this%data(:,1,Physics%PRESSURE) &
                 * Physics%constants%RG / (Physics%mu*(Physics%gamma-1))&
                 * rvar(:,Mesh%JMIN-1,Physics%DENSITY)
              rvar(:,Mesh%JMIN-2,Physics%ENERGY) = &
                 this%data(:,1,Physics%PRESSURE) &
                 * Physics%constants%RG / (Physics%mu*(Physics%gamma-1))&
                 * rvar(:,Mesh%JMIN-2,Physics%DENSITY)
           END WHERE
        END IF
     END IF
    CASE(NORTH)
       FORALL (i=Mesh%IGMIN:Mesh%IGMAX)
          WHERE (this%reflY)
             rvar(i,Mesh%JMAX+1,:) = -rvar(i,Mesh%JMAX,:)
             rvar(i,Mesh%JMAX+2,:) = -rvar(i,Mesh%JMAX-1,:)
          ELSEWHERE
             rvar(i,Mesh%JMAX+1,:) = rvar(i,Mesh%JMAX,:)
             rvar(i,Mesh%JMAX+2,:) = rvar(i,Mesh%JMAX-1,:)
          END WHERE
       END FORALL

       IF (PrimRecon(Fluxes%reconstruction).EQV.PRIMITIVE) THEN
          !set tangential velocity of boundary (x,z)!
       IF (Physics%XVELOCITY /= 0) THEN
          rvar(:,Mesh%JMAX+1,Physics%XVELOCITY) = &
             this%data(:,1,Physics%XVELOCITY)
          rvar(:,Mesh%JMAX+2,Physics%XVELOCITY) = &
             this%data(:,1,Physics%XVELOCITY)
       END IF
       IF (Physics%ZVELOCITY /= 0) THEN
          rvar(:,Mesh%JMAX+1,Physics%ZVELOCITY) = &
             this%data(:,1,Physics%ZVELOCITY)
          rvar(:,Mesh%JMAX+2,Physics%ZVELOCITY) = &
             this%data(:,1,Physics%ZVELOCITY)
       END IF
       !set temperature (if pressure exit and temp != 0)
       IF (Physics%PRESSURE /= 0) THEN
          WHERE (this%data(:,1,Physics%PRESSURE) > 0)
             rvar(:,Mesh%JMAX+1,Physics%PRESSURE) = &
                this%data(:,1,Physics%PRESSURE)* &
                rvar(:,Mesh%JMAX+1,Physics%DENSITY) &
                * Physics%constants%RG / Physics%mu
             rvar(:,Mesh%JMAX+2,Physics%PRESSURE) = &
                this%data(:,1,Physics%PRESSURE)* &
                rvar(:,Mesh%JMAX+2,Physics%DENSITY) &
                * Physics%constants%RG / Physics%mu
          END WHERE
       END IF
     ELSE
        IF (Physics%XVELOCITY /= 0) THEN
           rvar(:,Mesh%JMAX+1,Physics%XMOMENTUM) = &
              this%data(:,1,Physics%XVELOCITY)&
              *rvar(:,Mesh%JMAX+1,Physics%DENSITY)
           rvar(:,Mesh%JMAX+2,Physics%XMOMENTUM) = &
              this%data(:,1,Physics%XVELOCITY)&
              *rvar(:,Mesh%JMAX+2,Physics%DENSITY)
        END IF
        IF (Physics%ZVELOCITY /= 0) THEN
           rvar(:,Mesh%JMAX+1,Physics%ZVELOCITY) = &
              this%data(:,1,Physics%ZMOMENTUM)&
              *rvar(:,Mesh%JMAX+1,Physics%DENSITY)
           rvar(:,Mesh%JMAX+2,Physics%ZMOMENTUM) = &
              this%data(:,1,Physics%ZVELOCITY)&
              *rvar(:,Mesh%JMAX+2,Physics%DENSITY)
        END IF
        !set temperature (if pressure exit and temp != 0)
        IF (Physics%PRESSURE /= 0) THEN
           WHERE (this%data(:,1,Physics%PRESSURE) > 0)
              rvar(:,Mesh%JMAX+1,Physics%ENERGY) = &
                 this%data(:,1,Physics%PRESSURE) &
                 * Physics%constants%RG / (Physics%mu*(Physics%gamma-1))&
                 * rvar(:,Mesh%JMAX+1,Physics%DENSITY)
              rvar(:,Mesh%JMAX+2,Physics%ENERGY) = &
                 this%data(:,1,Physics%PRESSURE) &
                 * Physics%constants%RG / (Physics%mu*(Physics%gamma-1))&
                 * rvar(:,Mesh%JMAX+2,Physics%DENSITY)
           END WHERE
        END IF
     END IF
    END SELECT
  END SUBROUTINE CenterBoundary_moving_wall


!***************************************************
! until now: no velo for face boundary !!!!!!!!!!!!
!***************************************************

  PURE SUBROUTINE FaceBoundary_moving_wall(this,Mesh,Physics,we,ea,so,no,rstates)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: we,ea,so,no
    REAL :: rstates(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,Physics%vnum)
    !------------------------------------------------------------------------!
    INTEGER       :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,Mesh,Physics,we,ea,so,no
    INTENT(INOUT) :: rstates
    !------------------------------------------------------------------------!
    !************************************************************************!
    ! Be careful! There is a problem with trapezoidal rule, because          !
    ! SetFaceBoundary is called twice with different pairs of boundary values!
    ! (1st call:sw/se,sw/nw; 2nd call:nw/ne,se/ne).                          !
    !************************************************************************!
    SELECT CASE(GetDirection(this))
    CASE(WEST)
       FORALL (j=Mesh%JGMIN:Mesh%JGMAX)
          WHERE (this%reflX)
             rstates(Mesh%IMIN,j,we,:) = 0.
          END WHERE
          rstates(Mesh%IMIN-1,j,ea,:) = rstates(Mesh%IMIN,j,we,:)
       END FORALL
    CASE(EAST)
       FORALL (j=Mesh%JGMIN:Mesh%JGMAX)
          WHERE (this%reflX)
             rstates(Mesh%IMAX,j,ea,:) = 0.
          END WHERE
          rstates(Mesh%IMAX+1,j,we,:) = rstates(Mesh%IMAX,j,ea,:)
       END FORALL
    CASE(SOUTH)
       FORALL (i=Mesh%IGMIN:Mesh%IGMAX)
          WHERE (this%reflY)
             rstates(i,Mesh%JMIN,so,:) = 0.
          END WHERE
          rstates(i,Mesh%JMIN-1,no,:) = rstates(i,Mesh%JMIN,so,:) 
       END FORALL
    CASE(NORTH)
       FORALL (i=Mesh%IGMIN:Mesh%IGMAX)
          WHERE (this%reflY)
             rstates(i,Mesh%JMAX,no,:) = 0.
          END WHERE
          rstates(i,Mesh%JMAX+1,so,:) = rstates(i,Mesh%JMAX,no,:) 
       END FORALL
    END SELECT
  END SUBROUTINE FaceBoundary_moving_wall


  SUBROUTINE CloseBoundary_moving_wall(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%data,this%reflX,this%reflY)
  END SUBROUTINE CloseBoundary_moving_wall

END MODULE boundary_moving_wall
