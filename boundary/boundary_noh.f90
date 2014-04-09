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
  USE fluxes_common, ONLY : Fluxes_TYP
  USE reconstruction_common, ONLY : Reconstruction_TYP, PrimRecon
  USE physics_common, ONLY : Physics_TYP
  USE boundary_nogradients
  USE boundary_fixed
  USE physics_generic
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
       ALLOCATE(this%invr(Mesh%GNUM,Mesh%JMIN:Mesh%JMAX), &
            STAT=err)
    CASE(SOUTH,NORTH)
       ALLOCATE(this%invr(Mesh%IMIN:Mesh%IMAX,Mesh%GNUM), &
            STAT=err)
    END SELECT
    IF (err.NE.0) THEN
       CALL Error(this,"InitBoundary_noh", "Unable to allocate memory.")
    END IF
    ! compute the inverse distances using cartesian coordinates
    SELECT CASE(GetDirection(this))
    CASE(WEST)
       this%invr(:,Mesh%JMIN:Mesh%JMAX) =&
            1./SQRT(Mesh%bccart(Mesh%IGMIN:Mesh%IMIN-1,Mesh%JMIN:Mesh%JMAX,1)**2 &
            + Mesh%bccart(Mesh%IGMIN:Mesh%IMIN-1,Mesh%JMIN:Mesh%JMAX,2)**2)
    CASE(EAST)
       this%invr(:,Mesh%JMIN:Mesh%JMAX) =&
            1./SQRT(Mesh%bccart(Mesh%IMAX+1:Mesh%IGMAX,Mesh%JMIN:Mesh%JMAX,1)**2 &
            + Mesh%bccart(Mesh%IMAX+1:Mesh%IGMAX,Mesh%JMIN:Mesh%JMAX,2)**2)
    CASE(SOUTH)
       this%invr(Mesh%IMIN:Mesh%IMAX,:) =&
            1./SQRT(Mesh%bccart(Mesh%IMIN:Mesh%IMAX,Mesh%JGMIN:Mesh%JMIN-1,1)**2 &
            + Mesh%bccart(Mesh%IMIN:Mesh%IMAX,Mesh%JGMIN:Mesh%JMIN-1,2)**2)
    CASE(NORTH)
       this%invr(Mesh%IMIN:Mesh%IMAX,:) =&
            1./SQRT(Mesh%bccart(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+1:Mesh%JGMAX,1)**2 &
            + Mesh%bccart(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+1:Mesh%JGMAX,2)**2)
    END SELECT
    ! dimensional constant of the Noh problem;
    ! for 1D, 2 for 2D, 3 for 3D
    IF (PRESENT(dim)) THEN
       this%nohdim = dim
    ELSE ! default
       this%nohdim = 3
    END IF
   END SUBROUTINE InitBoundary_noh


  PURE SUBROUTINE CenterBoundary_noh(this,Mesh,Physics,Fluxes,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    REAL               :: time
    REAL :: pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM)
     !------------------------------------------------------------------------!
    INTEGER       :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,Mesh,Physics,Fluxes,time
    INTENT(INOUT) :: pvar 
    !------------------------------------------------------------------------!
    SELECT CASE(GetDirection(this))
    CASE(WEST)
       ! UNROLL=Mesh%GNUM would be sufficient, but the compiler does
       ! not know the value of Mesh%GNUM, hence we set UNROLL=4 and
       ! hope that nobody sets Mesh%GNUM to a value greater than 4
!CDIR UNROLL=4
       DO i=1,Mesh%GNUM
          ! the field normally reserved for the density in the
          ! auxiliary data array contains the inverse of the distance to the
          ! origin; we use this to set the time dependend density
          pvar(Mesh%IGMIN+i-1,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) &
               = this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) &
               * (1.0 + time*this%invr(i,Mesh%JMIN:Mesh%JMAX))**(this%nohdim-1)
          ! set fixed boundary data for velocities and pressure
          pvar(Mesh%IGMIN+i-1,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY:Physics%VNUM) &
               = this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY:Physics%VNUM)
       END DO
    CASE(EAST)
       ! UNROLL=Mesh%GNUM would be sufficient, but the compiler does
       ! not know the value of Mesh%GNUM, hence we set UNROLL=4 and
       ! hope that nobody sets Mesh%GNUM to a value greater than 4
!CDIR UNROLL=4
       DO i=1,Mesh%GNUM
          ! the field normally reserved for the density in the
          ! auxiliary data array contains the inverse of the distance to the
          ! origin; we use this to set the time dependend density
          pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) &
               = this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) &
               * (1.0 + time*this%invr(i,Mesh%JMIN:Mesh%JMAX))**(this%nohdim-1)
          ! set fixed boundary data for velocities and pressure
          pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY:Physics%VNUM) &
               = this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY:Physics%VNUM)
       END DO
    CASE(SOUTH)
       ! UNROLL=Mesh%GNUM would be sufficient, but the compiler does
       ! not know the value of Mesh%GNUM, hence we set UNROLL=4 and
       ! hope that nobody sets Mesh%GNUM to a value greater than 4
!CDIR UNROLL=4
       DO j=1,Mesh%GNUM
          ! the field normally reserved for the density in the
          ! auxiliary data array contains the inverse of the distance to the
          ! origin; we use this to set the time dependend density
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JGMIN+j-1,Physics%DENSITY) &
               = this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%DENSITY) &
               * (1.0 + time*this%invr(Mesh%IMIN:Mesh%IMAX,j))**(this%nohdim-1)
          ! set fixed boundary data for velocities and pressure
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JGMIN+j-1,Physics%XVELOCITY:Physics%VNUM) &
               = this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%XVELOCITY:Physics%VNUM)
       END DO
    CASE(NORTH)
       ! UNROLL=Mesh%GNUM would be sufficient, but the compiler does
       ! not know the value of Mesh%GNUM, hence we set UNROLL=4 and
       ! hope that nobody sets Mesh%GNUM to a value greater than 4
!CDIR UNROLL=4
       DO j=1,Mesh%GNUM
          ! the field normally reserved for the density in the
          ! auxiliary data array contains the inverse of the distance to the
          ! origin; we use this to set the time dependend density
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%DENSITY) &
               = this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%DENSITY) &
               * (1.0 + time*this%invr(Mesh%IMIN:Mesh%IMAX,j))**(this%nohdim-1)
          ! set fixed boundary data for velocities and pressure
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%XVELOCITY:Physics%VNUM) &
               = this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%XVELOCITY:Physics%VNUM)
       END DO
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
