!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_nogradients.f90                                          #
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
! boundary module for vanishing gradients (zero order extrapolation)
!----------------------------------------------------------------------------!
MODULE boundary_nogradients
  USE boundary_common
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "vanishing gradients"  
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Boundary_TYP, &
       ! constants
#ifdef PARALLEL
       DEFAULT_MPI_REAL, &
#endif
       WEST, EAST, SOUTH, NORTH, &
       ! methods
       InitBoundary, &
       InitBoundary_nogradients, &
       CenterBoundary_nogradients, &
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

  SUBROUTINE InitBoundary_nogradients(this,btype,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    INTEGER            :: btype,dir
    !------------------------------------------------------------------------!
    INTENT(IN)    :: btype,dir
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL InitBoundary(this,btype,boundcond_name,dir)
  END SUBROUTINE InitBoundary_nogradients


  PURE SUBROUTINE CenterBoundary_nogradients(this,Mesh,Physics,rvar)
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
       rvar(Mesh%IMIN-1,:,:) = rvar(Mesh%IMIN,:,:)
       rvar(Mesh%IMIN-2,:,:) = rvar(Mesh%IMIN,:,:)
    CASE(EAST)
       rvar(Mesh%IMAX+1,:,:) = rvar(Mesh%IMAX,:,:)
       rvar(Mesh%IMAX+2,:,:) = rvar(Mesh%IMAX,:,:)
    CASE(SOUTH)
       rvar(:,Mesh%JMIN-1,:) = rvar(:,Mesh%JMIN,:)
       rvar(:,Mesh%JMIN-2,:) = rvar(:,Mesh%JMIN,:)
    CASE(NORTH)
       rvar(:,Mesh%JMAX+1,:) = rvar(:,Mesh%JMAX,:)
       rvar(:,Mesh%JMAX+2,:) = rvar(:,Mesh%JMAX,:)
    END SELECT
  END SUBROUTINE CenterBoundary_nogradients

END MODULE boundary_nogradients
