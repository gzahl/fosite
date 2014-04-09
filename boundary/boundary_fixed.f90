!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_fixed.f90                                                #
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
! boundary module for sub/supersonic in/outflow with
! fixed boundary data
!----------------------------------------------------------------------------!
MODULE boundary_fixed
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  USE boundary_nogradients
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "fixed in/outflow"  
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! methods
       InitBoundary_fixed, &
       CenterBoundary_fixed, &
       CloseBoundary_fixed
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitBoundary_fixed(this,Mesh,Physics,btype,dir,bcname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: btype,dir
    CHARACTER(LEN=*), OPTIONAL :: bcname
    !------------------------------------------------------------------------!
    INTEGER            :: err
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,btype,dir,bcname
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF (PRESENT(bcname)) THEN
       CALL InitBoundary(this,btype,bcname,dir)
    ELSE
       CALL InitBoundary(this,btype,boundcond_name,dir)
    END IF
    ! allocate memory for boundary data and mask
    SELECT CASE(GetDirection(this))
    CASE(WEST,EAST)
       ALLOCATE(this%data(2,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
            this%fixed(Physics%vnum), &
            STAT=err)
    CASE(SOUTH,NORTH)
       ALLOCATE(this%data(Mesh%IGMIN:Mesh%IGMAX,2,Physics%vnum), &
            this%fixed(Physics%vnum), &
            STAT=err)
    END SELECT
    IF (err.NE.0) THEN
       CALL Error(this,"InitBoundary_fixed", "Unable to allocate memory.")
    END IF
    ! fixed(:) defaults to NO_GRADIENTS everywhere, so that fixed boundaries
    ! work even if the boundary data remains undefined
    this%fixed(:) = .FALSE.
  END SUBROUTINE InitBoundary_fixed


  PURE SUBROUTINE CenterBoundary_fixed(this,Mesh,Physics,rvar)
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
       DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=1,Mesh%GNUM
             WHERE(this%fixed)
                ! set fixed boundary data
                rvar(Mesh%IMIN-i,j,:) = this%data(i,j,:)
             ELSEWHERE
                ! first order extrapolation
                rvar(Mesh%IMIN-i,j,:) = 2.0*rvar(Mesh%IMIN-i+1,j,:) - rvar(Mesh%IMIN-i+2,j,:)
             END WHERE
          END DO
       END DO
    CASE(EAST)
       DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=1,Mesh%GNUM
             WHERE(this%fixed)
                ! set fixed boundary data
                rvar(Mesh%IMAX+i,j,:) = this%data(i,j,:)
             ELSEWHERE
                ! first order extrapolation
                rvar(Mesh%IMAX+i,j,:) = 2.0*rvar(Mesh%IMAX+i-1,j,:) - rvar(Mesh%IMAX+i-2,j,:)
             END WHERE
          END DO
       END DO
    CASE(SOUTH)
       DO j=1,Mesh%GNUM
          DO i=Mesh%IGMIN,Mesh%IGMAX
             WHERE(this%fixed)
                ! set fixed boundary data
                rvar(i,Mesh%JMIN-j,:) = this%data(i,j,:)
             ELSEWHERE
                ! first order extrapolation
                rvar(i,Mesh%JMIN-j,:) = 2.0*rvar(i,Mesh%JMIN-j+1,:) - rvar(i,Mesh%JMIN-j+2,:)
             END WHERE
          END DO
       ENd DO
    CASE(NORTH)
       DO j=1,Mesh%GNUM
          DO i=Mesh%IGMIN,Mesh%IGMAX
             WHERE(this%fixed)
                ! set fixed boundary data
                rvar(i,Mesh%JMAX+j,:) = this%data(i,j,:)
             ELSEWHERE
                ! first order extrapolation
                rvar(i,Mesh%JMAX+j,:) = 2.0*rvar(i,Mesh%JMAX+j-1,:) - rvar(i,Mesh%JMAX+j-2,:)
             END WHERE
          END DO
       END DO
    END SELECT
  END SUBROUTINE CenterBoundary_fixed


  SUBROUTINE CloseBoundary_fixed(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%data,this%fixed)
  END SUBROUTINE CloseBoundary_fixed

END MODULE boundary_fixed
