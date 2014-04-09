!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_periodic.f90                                             #
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
! module for periodic boundary conditions
!----------------------------------------------------------------------------!
MODULE boundary_periodic
  USE boundary_common
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "periodic"  
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Boundary_TYP, &
       ! constants
       WEST, EAST, SOUTH, NORTH, &
       ! methods
       InitBoundary_periodic, &
       GetType, &
       GetName, &
       GetDirection, &
       GetDirectionName, &
       CenterBoundary_periodic, &
       FaceBoundary_periodic
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitBoundary_periodic(this,Physics,btype,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: btype,dir
    !------------------------------------------------------------------------!
    INTENT(IN)    :: btype,dir
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL InitBoundary(this,btype,boundcond_name,dir)
  END SUBROUTINE InitBoundary_periodic


  PURE SUBROUTINE CenterBoundary_periodic(this,Mesh,Physics,rvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL :: rvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,Mesh,Physics
    INTENT(INOUT) :: rvar   
    !------------------------------------------------------------------------!
    SELECT CASE(GetDirection(this))
    CASE(WEST)
       rvar(Mesh%IMIN-1,:,:) = rvar(Mesh%IMAX,:,:)
       rvar(Mesh%IMIN-2,:,:) = rvar(Mesh%IMAX-1,:,:)
    CASE(EAST)
       rvar(Mesh%IMAX+1,:,:) = rvar(Mesh%IMIN,:,:)
       rvar(Mesh%IMAX+2,:,:) = rvar(Mesh%IMIN+1,:,:)
    CASE(SOUTH)
       rvar(:,Mesh%JMIN-1,:) = rvar(:,Mesh%JMAX,:)
       rvar(:,Mesh%JMIN-2,:) = rvar(:,Mesh%JMAX-1,:)
    CASE(NORTH)
       rvar(:,Mesh%JMAX+1,:) = rvar(:,Mesh%JMIN,:)
       rvar(:,Mesh%JMAX+2,:) = rvar(:,Mesh%JMIN+1,:)
    END SELECT
  END SUBROUTINE CenterBoundary_periodic


  PURE SUBROUTINE FaceBoundary_periodic(this,Mesh,Physics,we,ea,so,no,rstates)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: we,ea,so,no
    REAL :: rstates(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,Physics%vnum)
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
       rstates(Mesh%IMIN-1,:,ea,:) = rstates(Mesh%IMAX,:,ea,:)
       rstates(Mesh%IMIN-1,:,we,:) = rstates(Mesh%IMAX,:,we,:)
    CASE(EAST)
       rstates(Mesh%IMAX+1,:,we,:) = rstates(Mesh%IMIN,:,we,:)
       rstates(Mesh%IMAX+1,:,ea,:) = rstates(Mesh%IMIN,:,ea,:)
    CASE(SOUTH)
       rstates(:,Mesh%JMIN-1,no,:) = rstates(:,Mesh%JMAX,no,:)
       rstates(:,Mesh%JMIN-1,so,:) = rstates(:,Mesh%JMAX,so,:)
    CASE(NORTH)
       rstates(:,Mesh%JMAX+1,so,:) = rstates(:,Mesh%JMIN,so,:)
       rstates(:,Mesh%JMAX+1,no,:) = rstates(:,Mesh%JMIN,no,:)
    END SELECT
  END SUBROUTINE FaceBoundary_periodic

END MODULE boundary_periodic
