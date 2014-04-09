!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_fixed.f90                                                #
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
!! \brief Boundary module for fixed in/outflow conditions
!! 
!! Implementation of sub/supersonic in/outflow conditions with
!! user defined fixed data.
!!
!! \extends boundary_nogradients
!! \ingroup boundary
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

  !> \public Constructor for fixed boundary conditions
  SUBROUTINE InitBoundary_fixed(this,Mesh,Physics,btype,dir,bcname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: btype,dir
    CHARACTER(LEN=*), OPTIONAL :: bcname
    !------------------------------------------------------------------------!
    INTEGER            :: err = 0
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
       ALLOCATE(this%data(Mesh%GNUM,Mesh%JMIN:Mesh%JMAX,Physics%VNUM), &
            this%fixed(Mesh%JMIN:Mesh%JMAX,Physics%VNUM), &
            STAT=err)
    CASE(SOUTH,NORTH)
       ALLOCATE(this%data(Mesh%IMIN:Mesh%IMAX,Mesh%GNUM,Physics%VNUM), &
            this%fixed(Mesh%IMIN:Mesh%IMAX,Physics%VNUM), &
            STAT=err)
    END SELECT
    IF (err.NE.0) THEN
       CALL Error(this,"InitBoundary_fixed", "Unable to allocate memory.")
    END IF
    ! fixed(:,:) defaults to EXTRAPOLATION everywhere
    this%fixed(:,:) = .FALSE.
  END SUBROUTINE InitBoundary_fixed


  !> \public Applies the fixed boundary condition
  PURE SUBROUTINE CenterBoundary_fixed(this,Mesh,Physics,pvar)
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
          WHERE(this%fixed)
             ! set fixed boundary data
             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,:) = this%data(i,Mesh%JMIN:Mesh%JMAX,:)
          ELSEWHERE
             ! first order extrapolation
             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,:) = (i+1)*pvar(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,:) &
                  - i*pvar(Mesh%IMIN+1,Mesh%JMIN:Mesh%JMAX,:)
          END WHERE
       END DO
    CASE(EAST)
!CDIR UNROLL=4
       DO i=1,Mesh%GNUM
          WHERE(this%fixed)
             ! set fixed boundary data
             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,:) = this%data(i,Mesh%JMIN:Mesh%JMAX,:)
          ELSEWHERE
             ! first order extrapolation
             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,:) = (i+1)*pvar(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:) &
                  - i*pvar(Mesh%IMAX-1,Mesh%JMIN:Mesh%JMAX,:)
          END WHERE
       END DO
    CASE(SOUTH)
!CDIR UNROLL=4
       DO j=1,Mesh%GNUM
          WHERE(this%fixed)
             ! set fixed boundary data
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,:) = this%data(Mesh%IMIN:Mesh%IMAX,j,:)
          ELSEWHERE
             ! first order extrapolation
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,:) = (j+1)*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,:) &
                  - j*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+1,:)
          END WHERE
       END DO
    CASE(NORTH)
!CDIR UNROLL=4
       DO j=1,Mesh%GNUM
          WHERE(this%fixed)
             ! set fixed boundary data
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,:) = this%data(Mesh%IMIN:Mesh%IMAX,j,:)
          ELSEWHERE
             ! first order extrapolation
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,:) = (j+1)*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX,:) &
                  - j*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-1,:)
          END WHERE
       END DO
    END SELECT
  END SUBROUTINE CenterBoundary_fixed

  !> \public Destructor for fixed boundary conditions
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
