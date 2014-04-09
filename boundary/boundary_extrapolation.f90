!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_extrapolation.f90                                        #
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
!! \brief Boundary module for first order extrapolation
!!
!! \extends boundary_nogradients
!! \ingroup boundary
!----------------------------------------------------------------------------!
MODULE boundary_extrapolation
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  USE boundary_nogradients
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "space extrapolation"  
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Boundary_TYP, &
       ! constants
       WEST, EAST, SOUTH, NORTH, &
       ! methods
       InitBoundary_extrapolation, &
       CenterBoundary_extrapolation
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor for extrapolation boundary conditions
  SUBROUTINE InitBoundary_extrapolation(this,Physics,btype,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: btype,dir
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Physics,btype,dir
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL InitBoundary(this,btype,boundcond_name,dir)
  END SUBROUTINE InitBoundary_extrapolation


  !> \public Applies the extrapolation boundary condition
  PURE SUBROUTINE CenterBoundary_extrapolation(this,Mesh,Physics,pvar)
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
       ! UNROLL=Mesh%GNUM would be sufficient, but the compiler does
       ! not know the value of Mesh%GNUM, hence we set UNROLL=4 and
       ! hope that nobody sets Mesh%GNUM to a value greater than 4
!CDIR UNROLL=4
       DO i=1,Mesh%GNUM
          pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,:) = (i+1)*pvar(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,:) &
               - i*pvar(Mesh%IMIN+1,Mesh%JMIN:Mesh%JMAX,:)
       END DO
    CASE(EAST)
!CDIR UNROLL=4
       DO i=1,Mesh%GNUM
          pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,:) = (i+1)*pvar(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:) &
               - i*pvar(Mesh%IMAX-1,Mesh%JMIN:Mesh%JMAX,:)
       END DO
    CASE(SOUTH)
!CDIR UNROLL=4
       DO j=1,Mesh%GNUM
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,:) = (j+1)*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,:) &
               - j*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+1,:)
       END DO
    CASE(NORTH)
!CDIR UNROLL=4
       DO j=1,Mesh%GNUM
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,:) = (j+1)*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX,:) &
               - j*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-1,:)
       END DO
    END SELECT
  END SUBROUTINE CenterBoundary_extrapolation

END MODULE boundary_extrapolation
