!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_noh.f90                                                  #
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
! boundary module for inflow boundary conditions of the Noh problem
! all values should be set to the fix auxiliary data except for the time
! dependend density
!----------------------------------------------------------------------------!
MODULE boundary_noh
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  USE boundary_nogradients
  USE boundary_fixed
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "Noh problem supersonic inflow"  
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! methods
       InitBoundary_noh, &
       CenterBoundary_noh, &
       CloseBoundary_noh
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitBoundary_noh(this,Mesh,Physics,btype,dir,dim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: btype,dir
    INTEGER, OPTIONAL  :: dim
    !------------------------------------------------------------------------!
    INTEGER            :: err
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,btype,dir,dim
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL InitBoundary_fixed(this,Mesh,Physics,btype,dir,boundcond_name)
    ! allocate memory for inverse distances of boundary cells to the origin
    SELECT CASE(GetDirection(this))
    CASE(WEST,EAST)
       ALLOCATE(this%invr(2,Mesh%JGMIN:Mesh%JGMAX), &
            STAT=err)
    CASE(SOUTH,NORTH)
       ALLOCATE(this%invr(Mesh%IGMIN:Mesh%IGMAX,2), &
            STAT=err)
    END SELECT
    IF (err.NE.0) THEN
       CALL Error(this,"InitBoundary_noh", "Unable to allocate memory.")
    END IF
    ! compute the inverse distances using cartesian coordinates
    SELECT CASE(GetDirection(this))
    CASE(WEST)
       this%invr(:,:) = 1./SQRT(Mesh%bccart(Mesh%IMIN-2:Mesh%IMIN-1,:,1)**2 &
            + Mesh%bccart(Mesh%IMIN-2:Mesh%IMIN-1,:,2)**2)
    CASE(EAST)
       this%invr(:,:) = 1./SQRT(Mesh%bccart(Mesh%IMAX+1:Mesh%IMAX+2,:,1)**2 &
            + Mesh%bccart(Mesh%IMAX+1:Mesh%IMAX+2,:,2)**2)
    CASE(SOUTH)
       this%invr(:,:) = 1./SQRT(Mesh%bccart(:,Mesh%JMIN-2:Mesh%JMIN-1,1)**2 &
            + Mesh%bccart(:,Mesh%JMIN-2:Mesh%JMIN-1,2)**2)
    CASE(NORTH)
       this%invr(:,:) = 1./SQRT(Mesh%bccart(:,Mesh%JMAX+1:Mesh%JMAX+2,1)**2 &
            + Mesh%bccart(:,Mesh%JMAX+1:Mesh%JMAX+2,2)**2)
    END SELECT
    ! dimensional constant of the Noh problem;
    ! for 1D, 2 for 2D, 3 for 3D
    IF (PRESENT(dim)) THEN
       this%nohdim = dim
    ELSE ! default
       this%nohdim = 3
    END IF
   END SUBROUTINE InitBoundary_noh


  PURE SUBROUTINE CenterBoundary_noh(this,Mesh,Physics,time,rvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL               :: time
    REAL :: rvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM)
    !------------------------------------------------------------------------!
    INTEGER       :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,Mesh,Physics,time
    INTENT(INOUT) :: rvar   
    !------------------------------------------------------------------------!
    SELECT CASE(GetDirection(this))
    CASE(WEST)
       FORALL (i=1:2,j=Mesh%JGMIN:Mesh%JGMAX)
          ! the field normally reserved for the density in the
          ! auxiliary data array contains the inverse of the distance to the
          ! origin; we use this to set the time dependend density
          rvar(Mesh%IMIN-i,j,Physics%DENSITY)  = this%data(3-i,j,Physics%DENSITY) &
               * (1.0 + time*this%invr(3-i,j))**(this%nohdim-1)
          ! set fixed boundary data for velocities and pressure
          rvar(Mesh%IMIN-i,j,Physics%XVELOCITY:Physics%VNUM) &
               = this%data(3-i,j,Physics%XVELOCITY:Physics%VNUM)
       END FORALL
    CASE(EAST)
       FORALL (i=1:2,j=Mesh%JGMIN:Mesh%JGMAX)
          ! the field normally reserved for the density in the
          ! auxiliary data array contains the inverse of the distance to the
          ! origin; we use this to set the time dependend density
          rvar(Mesh%IMAX+i,j,Physics%DENSITY) = this%data(i,j,Physics%DENSITY) &
               * (1.0 + time*this%invr(i,j))**(this%nohdim-1)
          ! set fixed boundary data for velocities and pressure
          rvar(Mesh%IMAX+i,j,Physics%XVELOCITY:Physics%VNUM) &
               = this%data(i,j,Physics%XVELOCITY:Physics%VNUM)
       END FORALL
    CASE(SOUTH)
       FORALL (i=Mesh%IGMIN:Mesh%IGMAX,j=1:2)
          ! the field normally reserved for the density in the
          ! auxiliary data array contains the inverse of the distance to the
          ! origin; we use this to set the time dependend density
          rvar(i,Mesh%JMIN-j,Physics%DENSITY) = this%data(i,3-j,Physics%DENSITY) &
               * (1.0 + time*this%invr(i,3-j))**(this%nohdim-1)
          ! set fixed boundary data for velocities and pressure
          rvar(i,Mesh%JMIN-j,Physics%XVELOCITY:Physics%VNUM) &
               = this%data(i,3-j,Physics%XVELOCITY:Physics%VNUM)
       END FORALL
    CASE(NORTH)
       FORALL (i=Mesh%IGMIN:Mesh%IGMAX,j=1:2)
          ! the field normally reserved for the density in the
          ! auxiliary data array contains the inverse of the distance to the
          ! origin; we use this to set the time dependend density
          rvar(i,Mesh%JMAX+j,Physics%DENSITY) = this%data(i,j,Physics%DENSITY) &
               * (1.0 + time*this%invr(i,j))**(this%nohdim-1)
          ! set fixed boundary data for velocities and pressure
          rvar(i,Mesh%JMAX+j,Physics%XVELOCITY:Physics%VNUM) &
               = this%data(i,j,Physics%XVELOCITY:Physics%VNUM)
       END FORALL
    END SELECT
  END SUBROUTINE CenterBoundary_noh


  SUBROUTINE CloseBoundary_noh(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%invr)
    CALL CloseBoundary_fixed(this)
  END SUBROUTINE CloseBoundary_noh

END MODULE boundary_noh
