!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_reflecting.f90                                           #
!#                                                                           #
!# Copyright (C) 2006-2008                                                   #
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
! boundary module for refelcting boundaries
!----------------------------------------------------------------------------!
MODULE boundary_reflecting
  USE mesh_common, ONLY : Mesh_TYP
  USE fluxes_common, ONLY : Fluxes_TYP
  USE reconstruction_common, ONLY : Reconstruction_TYP, PrimRecon
  USE boundary_nogradients
  USE physics_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "reflecting"  
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Boundary_TYP, &
       ! constants
       WEST, EAST, SOUTH, NORTH, &
       ! methods
       InitBoundary, &
       InitBoundary_reflecting, &
       CenterBoundary_reflecting, &
       CloseBoundary_reflecting, &
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

  SUBROUTINE InitBoundary_reflecting(this,Physics,btype,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: btype,dir
    !------------------------------------------------------------------------!
    INTEGER       :: err
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Physics,btype,dir
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL InitBoundary(this,btype,boundcond_name,dir)
    
    ALLOCATE(this%reflX(Physics%vnum), &
         this%reflY(Physics%vnum), &
         STAT=err)
    IF (err.NE.0) THEN
       CALL Error(this, "InitBoundary_reflecting", "Unable to allocate memory.")
    END IF
    ! this tells us which vars get the opposite sign/vanish at cell faces;
    ! e.g. vertical velocities (depends on the underlying physics)
    CALL ReflectionMasks(Physics,this%reflX,this%reflY)
  END SUBROUTINE InitBoundary_reflecting


  PURE SUBROUTINE CenterBoundary_reflecting(this,Mesh,Physics,Fluxes,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    REAL :: pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)
    !------------------------------------------------------------------------!
    INTEGER       :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,Mesh,Physics,Fluxes
    INTENT(INOUT) :: pvar 
    !------------------------------------------------------------------------!
    SELECT CASE(GetDirection(this))
    CASE(WEST)
       FORALL (j=Mesh%JMIN:Mesh%JMAX)
          WHERE (this%reflX)
             pvar(Mesh%IMIN-1,j,:) = -pvar(Mesh%IMIN,j,:)
             pvar(Mesh%IMIN-2,j,:) = -pvar(Mesh%IMIN+1,j,:)
          ELSEWHERE
             pvar(Mesh%IMIN-1,j,:) = pvar(Mesh%IMIN,j,:)
             pvar(Mesh%IMIN-2,j,:) = pvar(Mesh%IMIN+1,j,:)
          END WHERE
       END FORALL
    CASE(EAST)
       FORALL (j=Mesh%JMIN:Mesh%JMAX)
          WHERE (this%reflX)
             pvar(Mesh%IMAX+1,j,:) = -pvar(Mesh%IMAX,j,:)
             pvar(Mesh%IMAX+2,j,:) = -pvar(Mesh%IMAX-1,j,:)
          ELSEWHERE
             pvar(Mesh%IMAX+1,j,:) = pvar(Mesh%IMAX,j,:)
             pvar(Mesh%IMAX+2,j,:) = pvar(Mesh%IMAX-1,j,:)
          END WHERE
       END FORALL
    CASE(SOUTH)
       FORALL (i=Mesh%IMIN:Mesh%IMAX)
          WHERE (this%reflY)
             pvar(i,Mesh%JMIN-1,:) = -pvar(i,Mesh%JMIN,:)
             pvar(i,Mesh%JMIN-2,:) = -pvar(i,Mesh%JMIN+1,:)
          ELSEWHERE
             pvar(i,Mesh%JMIN-1,:) = pvar(i,Mesh%JMIN,:)
             pvar(i,Mesh%JMIN-2,:) = pvar(i,Mesh%JMIN+1,:)
          END WHERE
       END FORALL
    CASE(NORTH)
       FORALL (i=Mesh%IMIN:Mesh%IMAX)
          WHERE (this%reflY)
             pvar(i,Mesh%JMAX+1,:) = -pvar(i,Mesh%JMAX,:)
             pvar(i,Mesh%JMAX+2,:) = -pvar(i,Mesh%JMAX-1,:)
          ELSEWHERE
             pvar(i,Mesh%JMAX+1,:) = pvar(i,Mesh%JMAX,:)
             pvar(i,Mesh%JMAX+2,:) = pvar(i,Mesh%JMAX-1,:)
          END WHERE
       END FORALL
    END SELECT 
  END SUBROUTINE CenterBoundary_reflecting

  SUBROUTINE CloseBoundary_reflecting(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%reflX,this%reflY)
  END SUBROUTINE CloseBoundary_reflecting

END MODULE boundary_reflecting
