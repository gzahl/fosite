!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_reflecting.f90                                           #
!#                                                                           #
!# Copyright (C) 2007 Tobias Illenseer <tillense@ita.uni-heidelberg.de>      #
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
  USE boundary_common
  USE mesh_common, ONLY : Mesh_TYP
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
       GetType, &
       GetName, &
       GetDirection, &
       GetDirectionName, &
       CenterBoundary_reflecting, &
       FaceBoundary_reflecting, &
       CloseBoundary_reflecting
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
       PRINT *, "ERROR in InitBoundary_reflecting: Can't allocate memory!"
       STOP
    END IF
    ! this tells us which vars get the opposite sign/vanish at cell faces;
    ! e.g. vertical velocities (depends on the underlying physics)
    CALL GetReflectionMasks(Physics,this%reflX,this%reflY)
  END SUBROUTINE InitBoundary_reflecting


  PURE SUBROUTINE CenterBoundary_reflecting(this,Mesh,Physics,rvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL :: rvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)
    !------------------------------------------------------------------------!
    INTEGER       :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,Mesh,Physics
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
    END SELECT
  END SUBROUTINE CenterBoundary_reflecting


  PURE SUBROUTINE FaceBoundary_reflecting(this,Mesh,Physics,we,ea,so,no,rstates)
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
  END SUBROUTINE FaceBoundary_reflecting


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
