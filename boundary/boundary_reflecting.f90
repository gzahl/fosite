!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_reflecting.f90                                           #
!#                                                                           #
!# Copyright (C) 2006-2014                                                   #
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
!! \brief Boundary module for refelcting boundaries
!!
!! \extends boundary_nogradients
!! \ingroup boundary
!----------------------------------------------------------------------------!
MODULE boundary_reflecting
  USE mesh_common, ONLY : Mesh_TYP
  USE boundary_nogradients
  USE physics_generic, ONLY : Physics_TYP, ReflectionMasks
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
       CloseBoundary, &
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

  !> \public Constructor for reflecting boundary conditions
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

  !> \public Applies the reflecting boundary condition
  PURE SUBROUTINE CenterBoundary_reflecting(this,Mesh,Physics,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL :: pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)
    !------------------------------------------------------------------------!
    INTEGER       :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,Mesh,Physics
    INTENT(INOUT) :: pvar 
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetDirection(this))
    CASE(WEST)
!CDIR UNROLL=8
       DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
          DO i=1,Mesh%GNUM
             WHERE (this%reflX)
                pvar(Mesh%IMIN-i,j,:) = -pvar(Mesh%IMIN+i-1,j,:)
             ELSEWHERE
                pvar(Mesh%IMIN-i,j,:) = pvar(Mesh%IMIN+i-1,j,:)
             END WHERE
          END DO
       END DO
    CASE(EAST)
!CDIR UNROLL=8
       DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
          DO i=1,Mesh%GNUM
             WHERE (this%reflX)
                pvar(Mesh%IMAX+i,j,:) = -pvar(Mesh%IMAX-i+1,j,:)
             ELSEWHERE
                pvar(Mesh%IMAX+i,j,:) = pvar(Mesh%IMAX-i+1,j,:)
             END WHERE
          END DO
       END DO
    CASE(SOUTH)
!CDIR UNROLL=4
       DO j=1,Mesh%GNUM
!CDIR NODEP
          DO i=Mesh%IMIN,Mesh%IMAX
             WHERE (this%reflY)
                pvar(i,Mesh%JMIN-j,:) = -pvar(i,Mesh%JMIN+j-1,:)
             ELSEWHERE
                pvar(i,Mesh%JMIN-j,:) = pvar(i,Mesh%JMIN+j-1,:)
             END WHERE
          END DO
       END DO
    CASE(NORTH)
!CDIR UNROLL=4
       DO j=1,Mesh%GNUM
!CDIR NODEP
          DO i=Mesh%IMIN,Mesh%IMAX
             WHERE (this%reflY)
                pvar(i,Mesh%JMAX+j,:) = -pvar(i,Mesh%JMAX-j+1,:)
             ELSEWHERE
                pvar(i,Mesh%JMAX+j,:) = pvar(i,Mesh%JMAX-j+1,:)
             END WHERE
          END DO
       END DO
    END SELECT 
  END SUBROUTINE CenterBoundary_reflecting

  !> \public Destructor for reflecting boundary conditions
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
