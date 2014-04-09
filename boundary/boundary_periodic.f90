!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_periodic.f90                                             #
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
! module for periodic boundary conditions
!----------------------------------------------------------------------------!
MODULE boundary_periodic
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  USE boundary_nogradients
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
       CenterBoundary_periodic
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitBoundary_periodic(this,btype,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    INTEGER            :: btype,dir
    !------------------------------------------------------------------------!
    INTENT(IN)    :: btype,dir
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL InitBoundary(this,btype,boundcond_name,dir)
  END SUBROUTINE InitBoundary_periodic


  PURE SUBROUTINE CenterBoundary_periodic(this,Mesh,Physics,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                       :: pvar
    !------------------------------------------------------------------------!
    INTEGER            :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this,Mesh,Physics
    INTENT(INOUT)      :: pvar   
    !------------------------------------------------------------------------!
    SELECT CASE(GetDirection(this))
    CASE(WEST)
!CDIR NODEP
       DO i=1,Mesh%GNUM
          pvar(Mesh%IMIN-i,:,:) = pvar(Mesh%IMAX-i+1,:,:)
       END DO
    CASE(EAST)
!CDIR NODEP
       DO i=1,Mesh%GNUM
          pvar(Mesh%IMAX+i,:,:) = pvar(Mesh%IMIN+i-1,:,:)
       END DO
    CASE(SOUTH)
!CDIR NODEP
       DO j=1,Mesh%GNUM
          pvar(:,Mesh%JMIN-j,:) = pvar(:,Mesh%JMAX-j+1,:)
       END DO
    CASE(NORTH)
!CDIR NODEP
       DO j=1,Mesh%GNUM
          pvar(:,Mesh%JMAX+j,:) = pvar(:,Mesh%JMIN+j-1,:)
       END DO
    END SELECT
  END SUBROUTINE CenterBoundary_periodic

END MODULE boundary_periodic
