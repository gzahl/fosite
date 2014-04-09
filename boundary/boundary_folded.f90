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
  USE fluxes_common, ONLY : Fluxes_TYP
  USE reconstruction_common, ONLY : Reconstruction_TYP, PrimRecon
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
#ifdef PARALLEL
    CALL Error(this,"InitBoundary_folded", "Boundary condition not supported in parallel mode.")
#endif
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


  PURE SUBROUTINE CenterBoundary_folded(this,Mesh,Physics,Fluxes,pvar)
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
       ! middle cell transmissive (no gradient) boundaries
       pvar(Mesh%IMIN-1,this%JMID+1,:) = pvar(Mesh%IMIN,this%JMID+1,:)
       pvar(Mesh%IMIN-2,this%JMID+1,:) = pvar(Mesh%IMIN,this%JMID+1,:)
       ! JGMAX -> JGMIN, JGMAX-1 -> JGMIN+1, ..., JGMAX-JMID -> JMID
       ! and JGMIN -> JGMAX, JGMIN+1 -> JGMAX-1, ...
       FORALL (j=Mesh%JMIN:this%JMID)
          WHERE (this%reflX)
             ! switch sign for vertical vector components
             pvar(Mesh%IMIN-1,j,:) = -pvar(Mesh%IMIN,Mesh%JMAX+(Mesh%JMIN-j),:)
             pvar(Mesh%IMIN-2,j,:) = -pvar(Mesh%IMIN+1,Mesh%JMAX+(Mesh%JMIN-j),:)
             pvar(Mesh%IMIN-1,Mesh%JMAX+(Mesh%JMIN-j),:) = -pvar(Mesh%IMIN,j,:)
             pvar(Mesh%IMIN-2,Mesh%JMAX+(Mesh%JMIN-j),:) = -pvar(Mesh%IMIN+1,j,:)
          ELSEWHERE
             pvar(Mesh%IMIN-1,j,:) = pvar(Mesh%IMIN,Mesh%JMAX+(Mesh%JMIN-j),:)
             pvar(Mesh%IMIN-2,j,:) = pvar(Mesh%IMIN+1,Mesh%JMAX+(Mesh%JMIN-j),:)
             pvar(Mesh%IMIN-1,Mesh%JMAX+(Mesh%JMIN-j),:) = pvar(Mesh%IMIN,j,:)
             pvar(Mesh%IMIN-2,Mesh%JMAX+(Mesh%JMIN-j),:) = pvar(Mesh%IMIN+1,j,:)
          END WHERE
       END FORALL
    CASE(EAST)
       ! middle cell transmissive (no gradient) boundaries
       pvar(Mesh%IMAX+1,this%JMID+1,:) = pvar(Mesh%IMAX,this%JMID+1,:)
       pvar(Mesh%IMAX+2,this%JMID+1,:) = pvar(Mesh%IMAX,this%JMID+1,:)
       ! JGMAX -> JGMIN, JGMAX-1 -> JGMIN+1, ..., JGMAX-JMID -> JMID
       ! and JGMIN -> JGMAX, JGMIN+1 -> JGMAX-1, ...
       FORALL (j=Mesh%JMIN:this%JMID)
          WHERE (this%reflX)
             ! switch sign for vertical vector components
             pvar(Mesh%IMAX+1,j,:) = -pvar(Mesh%IMAX,Mesh%JMAX+(Mesh%JMIN-j),:)
             pvar(Mesh%IMAX+2,j,:) = -pvar(Mesh%IMAX-1,Mesh%JMAX+(Mesh%JMIN-j),:)
             pvar(Mesh%IMAX+1,Mesh%JMAX+(Mesh%JMIN-j),:) = -pvar(Mesh%IMAX,j,:)
             pvar(Mesh%IMAX+2,Mesh%JMAX+(Mesh%JMIN-j),:) = -pvar(Mesh%IMAX-1,j,:)
          ELSEWHERE
             pvar(Mesh%IMAX+1,j,:) = pvar(Mesh%IMAX,Mesh%JMAX+(Mesh%JMIN-j),:)
             pvar(Mesh%IMAX+2,j,:) = pvar(Mesh%IMAX-1,Mesh%JMAX+(Mesh%JMIN-j),:)
             pvar(Mesh%IMAX+1,Mesh%JMAX+(Mesh%JMIN-j),:) = pvar(Mesh%IMAX,j,:)
             pvar(Mesh%IMAX+2,Mesh%JMAX+(Mesh%JMIN-j),:) = pvar(Mesh%IMAX-1,j,:)
          END WHERE
       END FORALL
    CASE(SOUTH)
       ! middle cell transmissive boundaries
       pvar(this%IMID+1,Mesh%JMIN-1,:) = pvar(this%IMID+1,Mesh%JMIN,:)
       pvar(this%IMID+1,Mesh%JMIN-2,:) = pvar(this%IMID+1,Mesh%JMIN,:)
       ! IGMAX -> IGMIN, IGMAX-1 -> IGMIN+1, ..., IGMAX-IMID -> IMID
       ! and IGMIN -> IGMAX, IGMIN+1 -> IGMAX-1, ...
       FORALL (i=Mesh%IMIN:this%IMID)
          WHERE (this%reflY)
             ! switch sign for vertical vector components
             pvar(i,Mesh%JMIN-1,:) = - pvar(Mesh%IMAX+(Mesh%IMIN-i),Mesh%JMIN,:)
             pvar(i,Mesh%JMIN-2,:) = - pvar(Mesh%IMAX+(Mesh%IMIN-i),Mesh%JMIN+1,:)
             pvar(Mesh%IMAX+(Mesh%IMIN-i),Mesh%JMIN-1,:) = - pvar(i,Mesh%JMIN,:)
             pvar(Mesh%IMAX+(Mesh%IMIN-i),Mesh%JMIN-2,:) = - pvar(i,Mesh%JMIN+1,:)
          ELSEWHERE
             pvar(i,Mesh%JMIN-1,:) = pvar(Mesh%IMAX+(Mesh%IMIN-i),Mesh%JMIN,:)
             pvar(i,Mesh%JMIN-2,:) = pvar(Mesh%IMAX+(Mesh%IMIN-i),Mesh%JMIN+1,:)
             pvar(Mesh%IMAX+(Mesh%IMIN-i),Mesh%JMIN-1,:) = pvar(i,Mesh%JMIN,:)
             pvar(Mesh%IMAX+(Mesh%IMIN-i),Mesh%JMIN-2,:) = pvar(i,Mesh%JMIN+1,:)
         END WHERE
      END FORALL
    CASE(NORTH)
       ! middle cell transmissive boundaries
       pvar(this%IMID+1,Mesh%JMAX+1,:) = pvar(this%IMID+1,Mesh%JMAX,:)
       pvar(this%IMID+1,Mesh%JMAX+2,:) = pvar(this%IMID+1,Mesh%JMAX,:)
       ! IGMAX -> IGMIN, IGMAX-1 -> IGMIN+1, ..., IGMAX-IMID -> IMID
       ! and IGMIN -> IGMAX, IGMIN+1 -> IGMAX-1, ...
       FORALL (i=Mesh%IMIN:this%IMID)
          WHERE (this%reflY)
             ! switch sign for vertical vector components
             pvar(i,Mesh%JMAX+1,:) = - pvar(Mesh%IMAX+(Mesh%IMIN-i),Mesh%JMAX,:)
             pvar(i,Mesh%JMAX+2,:) = - pvar(Mesh%IMAX+(Mesh%IMIN-i),Mesh%JMAX-1,:)
             pvar(Mesh%IMAX+(Mesh%IMIN-i),Mesh%JMAX+1,:) = - pvar(i,Mesh%JMAX,:)
             pvar(Mesh%IMAX+(Mesh%IMIN-i),Mesh%JMAX+2,:) = - pvar(i,Mesh%JMAX-1,:)
          ELSEWHERE
             pvar(i,Mesh%JMAX+1,:) = pvar(Mesh%IMAX+(Mesh%IMIN-i),Mesh%JMAX,:)
             pvar(i,Mesh%JMAX+2,:) = pvar(Mesh%IMAX+(Mesh%IMIN-i),Mesh%JMAX-1,:)
             pvar(Mesh%IMAX+(Mesh%IMIN-i),Mesh%JMAX+1,:) = pvar(i,Mesh%JMAX,:)
             pvar(Mesh%IMAX+(Mesh%IMIN-i),Mesh%JMAX+2,:) = pvar(i,Mesh%JMAX-1,:)
         END WHERE
      END FORALL
    END SELECT
  END SUBROUTINE CenterBoundary_folded

END MODULE boundary_folded
