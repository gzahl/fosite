!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_folded.f90                                               #
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
! boundary module for folded boundary
! i.e. copy value(IMIN)<->value(IMAX), value(IMIN+1)<->value(IMAX-1) ...
! use this for oblate spheroidal coordinates
!----------------------------------------------------------------------------!
MODULE boundary_folded
  USE mesh_common, ONLY : Mesh_TYP
  USE boundary_reflecting, CloseBoundary_folded => CloseBoundary_reflecting
  USE physics_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "folded"  
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Boundary_TYP, &
       ! constants
       WEST, EAST, SOUTH, NORTH, &
       ! methods
       InitBoundary_folded, &
       CenterBoundary_folded, &
       CloseBoundary_folded
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitBoundary_folded(this,Mesh,Physics,btype,dir)
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
       CALL Error(this,"InitBoundary_folded", "Unable to allocate memory.")
    END IF
    ! this tells us which vars get the opposite sign/vanish at cell faces;
    ! e.g. vertical velocities (depends on the underlying physics)
    CALL ReflectionMasks(Physics,this%reflX,this%reflY)
    ! odd cell numbers: IMID+1 and JMID+1 are in the middle 
    this%IMID = Mesh%INUM / 2 + Mesh%IMIN - 1
    this%JMID = Mesh%JNUM / 2 + Mesh%JMIN - 1
  END SUBROUTINE InitBoundary_folded


  PURE SUBROUTINE CenterBoundary_folded(this,Mesh,Physics,rvar)
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
       ! middle cell transmissive (no gradient) boundaries
       rvar(Mesh%IMIN-1,this%JMID+1,:) = rvar(Mesh%IMIN,this%JMID+1,:)
       rvar(Mesh%IMIN-2,this%JMID+1,:) = rvar(Mesh%IMIN,this%JMID+1,:)
       ! JGMAX -> JGMIN, JGMAX-1 -> JGMIN+1, ..., JGMAX-JMID -> JMID
       ! and JGMIN -> JGMAX, JGMIN+1 -> JGMAX-1, ...
       FORALL (j=Mesh%JGMIN:this%JMID)
          WHERE (this%reflX)
             ! switch sign for vertical vector components
             rvar(Mesh%IMIN-1,j,:) = -rvar(Mesh%IMIN,Mesh%JGMAX+(Mesh%JGMIN-j),:)
             rvar(Mesh%IMIN-2,j,:) = -rvar(Mesh%IMIN+1,Mesh%JGMAX+(Mesh%JGMIN-j),:)
             rvar(Mesh%IMIN-1,Mesh%JGMAX+(Mesh%JGMIN-j),:) = -rvar(Mesh%IMIN,j,:)
             rvar(Mesh%IMIN-2,Mesh%JGMAX+(Mesh%JGMIN-j),:) = -rvar(Mesh%IMIN+1,j,:)
          ELSEWHERE
             rvar(Mesh%IMIN-1,j,:) = rvar(Mesh%IMIN,Mesh%JGMAX+(Mesh%JGMIN-j),:)
             rvar(Mesh%IMIN-2,j,:) = rvar(Mesh%IMIN+1,Mesh%JGMAX+(Mesh%JGMIN-j),:)
             rvar(Mesh%IMIN-1,Mesh%JGMAX+(Mesh%JGMIN-j),:) = rvar(Mesh%IMIN,j,:)
             rvar(Mesh%IMIN-2,Mesh%JGMAX+(Mesh%JGMIN-j),:) = rvar(Mesh%IMIN+1,j,:)
          END WHERE
       END FORALL
    CASE(EAST)
       ! middle cell transmissive (no gradient) boundaries
       rvar(Mesh%IMAX+1,this%JMID+1,:) = rvar(Mesh%IMAX,this%JMID+1,:)
       rvar(Mesh%IMAX+2,this%JMID+1,:) = rvar(Mesh%IMAX,this%JMID+1,:)
       ! JGMAX -> JGMIN, JGMAX-1 -> JGMIN+1, ..., JGMAX-JMID -> JMID
       ! and JGMIN -> JGMAX, JGMIN+1 -> JGMAX-1, ...
       FORALL (j=Mesh%JGMIN:this%JMID)
          WHERE (this%reflX)
             ! switch sign for vertical vector components
             rvar(Mesh%IMAX+1,j,:) = -rvar(Mesh%IMAX,Mesh%JGMAX+(Mesh%JGMIN-j),:)
             rvar(Mesh%IMAX+2,j,:) = -rvar(Mesh%IMAX-1,Mesh%JGMAX+(Mesh%JGMIN-j),:)
             rvar(Mesh%IMAX+1,Mesh%JGMAX+(Mesh%JGMIN-j),:) = -rvar(Mesh%IMAX,j,:)
             rvar(Mesh%IMAX+2,Mesh%JGMAX+(Mesh%JGMIN-j),:) = -rvar(Mesh%IMAX-1,j,:)
          ELSEWHERE
             rvar(Mesh%IMAX+1,j,:) = rvar(Mesh%IMAX,Mesh%JGMAX+(Mesh%JGMIN-j),:)
             rvar(Mesh%IMAX+2,j,:) = rvar(Mesh%IMAX-1,Mesh%JGMAX+(Mesh%JGMIN-j),:)
             rvar(Mesh%IMAX+1,Mesh%JGMAX+(Mesh%JGMIN-j),:) = rvar(Mesh%IMAX,j,:)
             rvar(Mesh%IMAX+2,Mesh%JGMAX+(Mesh%JGMIN-j),:) = rvar(Mesh%IMAX-1,j,:)
          END WHERE
       END FORALL
    CASE(SOUTH)
       ! middle cell transmissive boundaries
       rvar(this%IMID+1,Mesh%JMIN-1,:) = rvar(this%IMID+1,Mesh%JMIN,:)
       rvar(this%IMID+1,Mesh%JMIN-2,:) = rvar(this%IMID+1,Mesh%JMIN,:)
       ! IGMAX -> IGMIN, IGMAX-1 -> IGMIN+1, ..., IGMAX-IMID -> IMID
       ! and IGMIN -> IGMAX, IGMIN+1 -> IGMAX-1, ...
       FORALL (i=Mesh%IGMIN:this%IMID)
          WHERE (this%reflY)
             ! switch sign for vertical vector components
             rvar(i,Mesh%JMIN-1,:) = - rvar(Mesh%IMAX+(Mesh%IGMIN-i),Mesh%JMIN,:)
             rvar(i,Mesh%JMIN-2,:) = - rvar(Mesh%IMAX+(Mesh%IGMIN-i),Mesh%JMIN+1,:)
             rvar(Mesh%IMAX+(Mesh%IGMIN-i),Mesh%JMIN-1,:) = - rvar(i,Mesh%JMIN,:)
             rvar(Mesh%IMAX+(Mesh%IGMIN-i),Mesh%JMIN-2,:) = - rvar(i,Mesh%JMIN+1,:)
          ELSEWHERE
             rvar(i,Mesh%JMIN-1,:) = rvar(Mesh%IMAX+(Mesh%IGMIN-i),Mesh%JMIN,:)
             rvar(i,Mesh%JMIN-2,:) = rvar(Mesh%IMAX+(Mesh%IGMIN-i),Mesh%JMIN+1,:)
             rvar(Mesh%IMAX+(Mesh%IGMIN-i),Mesh%JMIN-1,:) = rvar(i,Mesh%JMIN,:)
             rvar(Mesh%IMAX+(Mesh%IGMIN-i),Mesh%JMIN-2,:) = rvar(i,Mesh%JMIN+1,:)
         END WHERE
      END FORALL
    CASE(NORTH)
       ! middle cell transmissive boundaries
       rvar(this%IMID+1,Mesh%JMAX+1,:) = rvar(this%IMID+1,Mesh%JMAX,:)
       rvar(this%IMID+1,Mesh%JMAX+2,:) = rvar(this%IMID+1,Mesh%JMAX,:)
       ! IGMAX -> IGMIN, IGMAX-1 -> IGMIN+1, ..., IGMAX-IMID -> IMID
       ! and IGMIN -> IGMAX, IGMIN+1 -> IGMAX-1, ...
       FORALL (i=Mesh%IGMIN:this%IMID)
          WHERE (this%reflY)
             ! switch sign for vertical vector components
             rvar(i,Mesh%JMAX+1,:) = - rvar(Mesh%IMAX+(Mesh%IGMIN-i),Mesh%JMAX,:)
             rvar(i,Mesh%JMAX+2,:) = - rvar(Mesh%IMAX+(Mesh%IGMIN-i),Mesh%JMAX-1,:)
             rvar(Mesh%IMAX+(Mesh%IGMIN-i),Mesh%JMAX+1,:) = - rvar(i,Mesh%JMAX,:)
             rvar(Mesh%IMAX+(Mesh%IGMIN-i),Mesh%JMAX+2,:) = - rvar(i,Mesh%JMAX-1,:)
          ELSEWHERE
             rvar(i,Mesh%JMAX+1,:) = rvar(Mesh%IMAX+(Mesh%IGMIN-i),Mesh%JMAX,:)
             rvar(i,Mesh%JMAX+2,:) = rvar(Mesh%IMAX+(Mesh%IGMIN-i),Mesh%JMAX-1,:)
             rvar(Mesh%IMAX+(Mesh%IGMIN-i),Mesh%JMAX+1,:) = rvar(i,Mesh%JMAX,:)
             rvar(Mesh%IMAX+(Mesh%IGMIN-i),Mesh%JMAX+2,:) = rvar(i,Mesh%JMAX-1,:)
         END WHERE
      END FORALL
    END SELECT
  END SUBROUTINE CenterBoundary_folded

END MODULE boundary_folded
