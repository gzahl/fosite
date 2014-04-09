!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_fixed.f90                                                #
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
! boundary module for sub/supersonic in/outflow with
! fixed boundary data
!----------------------------------------------------------------------------!
MODULE boundary_fixed
  USE boundary_common
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "fixed in/outflow"  
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Boundary_TYP, &
       ! constants
       WEST, EAST, SOUTH, NORTH, &
       ! methods
       InitBoundary_fixed, &
       GetType, &
       GetName, &
       GetDirection, &
       GetDirectionName, &
       CenterBoundary_fixed, &
       FaceBoundary_fixed, &
       CloseBoundary_fixed
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitBoundary_fixed(this,Mesh,Physics,btype,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: btype,dir
    !------------------------------------------------------------------------!
    INTEGER            :: err
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,btype,dir
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL InitBoundary(this,btype,boundcond_name,dir)
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
       PRINT *, "ERROR in InitBoundary_fixed: Unable to allocate memory"
       STOP
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
       FORALL (i=1:2,j=Mesh%JGMIN:Mesh%JGMAX)
          WHERE(this%fixed)
             ! set fixed boundary data
             rvar(Mesh%IMIN-i,j,:) = this%data(i,j,:)
          ELSEWHERE
             ! first order extrapolation
             rvar(Mesh%IMIN-i,j,:) = 2.0*rvar(Mesh%IMIN-i+1,j,:) - rvar(Mesh%IMIN-i+2,j,:)
          END WHERE
       END FORALL
    CASE(EAST)
       FORALL (i=1:2,j=Mesh%JGMIN:Mesh%JGMAX)
          WHERE(this%fixed)
             ! set fixed boundary data
             rvar(Mesh%IMAX+i,j,:) = this%data(i,j,:)
          ELSEWHERE
             ! first order extrapolation
             rvar(Mesh%IMAX+i,j,:) = 2.0*rvar(Mesh%IMAX+i-1,j,:) - rvar(Mesh%IMAX+i-2,j,:)
          END WHERE
       END FORALL
    CASE(SOUTH)
       FORALL (i=Mesh%IGMIN:Mesh%IGMAX,j=1:2)
          WHERE(this%fixed)
             ! set fixed boundary data
             rvar(i,Mesh%JMIN-j,:) = this%data(i,j,:)
          ELSEWHERE
             ! first order extrapolation
             rvar(i,Mesh%JMIN-j,:) = 2.0*rvar(i,Mesh%JMIN-j+1,:) - rvar(i,Mesh%JMIN-j+2,:)
          END WHERE
       END FORALL
    CASE(NORTH)
       FORALL (i=Mesh%IGMIN:Mesh%IGMAX,j=1:2)
          WHERE(this%fixed)
             ! set fixed boundary data
             rvar(i,Mesh%JMAX+j,:) = this%data(i,j,:)
          ELSEWHERE
             ! first order extrapolation
             rvar(i,Mesh%JMAX+j,:) = 2.0*rvar(i,Mesh%JMAX+j-1,:) - rvar(i,Mesh%JMAX+j-2,:)
          END WHERE
       END FORALL
    END SELECT
  END SUBROUTINE CenterBoundary_fixed


  PURE SUBROUTINE FaceBoundary_fixed(this,Mesh,Physics,we,ea,so,no,rstates)
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
          WHERE(this%fixed)
             ! set fixed boundary data
             rstates(Mesh%IMIN-1,j,ea,:) = this%data(1,j,:)
          ELSEWHERE
             ! first order extrapolation
             rstates(Mesh%IMIN-1,j,ea,:) = 2.0*rstates(Mesh%IMIN,j,ea,:) - rstates(Mesh%IMIN+1,j,ea,:)
          END WHERE
       END FORALL
    CASE(EAST)
       FORALL (j=Mesh%JGMIN:Mesh%JGMAX)
          WHERE(this%fixed)
             ! set fixed boundary data
             rstates(Mesh%IMAX+1,j,we,:) = this%data(1,j,:)
          ELSEWHERE
             ! zero order extrapolation
             rstates(Mesh%IMAX+1,j,we,:) = 2.0*rstates(Mesh%IMAX,j,we,:) - rstates(Mesh%IMAX-1,j,we,:)
          END WHERE
       END FORALL
    CASE(SOUTH)
       FORALL (i=Mesh%IGMIN:Mesh%IGMAX)
          WHERE(this%fixed)
             ! set fixed boundary data
             rstates(i,Mesh%JMIN-1,no,:) = this%data(i,1,:)
          ELSEWHERE
             ! zero order extrapolation
             rstates(i,Mesh%JMIN-1,no,:) = 2.0*rstates(i,Mesh%JMIN,no,:) - rstates(i,Mesh%JMIN+1,no,:) 
          END WHERE
       END FORALL
    CASE(NORTH)
       FORALL (i=Mesh%IGMIN:Mesh%IGMAX)
          WHERE(this%fixed)
             ! set fixed boundary data
             rstates(i,Mesh%JMAX+1,so,:) = this%data(i,1,:)
          ELSEWHERE
             ! zero order extrapolation
             rstates(i,Mesh%JMAX+1,so,:) = 2.0*rstates(i,Mesh%JMAX,so,:) - rstates(i,Mesh%JMAX-1,so,:)
          END WHERE
       END FORALL
    END SELECT
  END SUBROUTINE FaceBoundary_fixed


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
