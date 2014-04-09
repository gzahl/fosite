!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_custom.f90                                               #
!#                                                                           #
!# Copyright (C) 2006-2012                                                   #
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
! custom boundary conditions
!----------------------------------------------------------------------------!
MODULE boundary_custom
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  USE boundary_nogradients
  USE boundary_fixed
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "custom"
  INTEGER, PARAMETER :: CUSTOM_NOGRAD   = 1
  INTEGER, PARAMETER :: CUSTOM_PERIOD   = 2
  INTEGER, PARAMETER :: CUSTOM_REFLECT  = 3
  INTEGER, PARAMETER :: CUSTOM_REFLNEG  = 4
  INTEGER, PARAMETER :: CUSTOM_EXTRAPOL = 5
  INTEGER, PARAMETER :: CUSTOM_FIXED    = 6
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! constants
       CUSTOM_NOGRAD, CUSTOM_PERIOD, CUSTOM_REFLECT, CUSTOM_REFLNEG, &
       CUSTOM_EXTRAPOL, CUSTOM_FIXED, &
       ! methods
       InitBoundary_custom, &
       CenterBoundary_custom, &
       CloseBoundary_custom
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitBoundary_custom(this,Mesh,Physics,btype,dir)
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
    CALL InitBoundary_fixed(this,Mesh,Physics,btype,dir,boundcond_name)
    ! allocate memory for boundary data and mask
    SELECT CASE(GetDirection(this))
    CASE(WEST,EAST)
       ALLOCATE(this%cbtype(Mesh%JMIN:Mesh%JMAX,Physics%VNUM), &
            STAT=err)
    CASE(SOUTH,NORTH)
       ALLOCATE(this%cbtype(Mesh%IMIN:Mesh%IMAX,Physics%VNUM), &
            STAT=err)
    END SELECT
    IF (err.NE.0) THEN
       CALL Error(this,"InitBoundary_custom", "Unable to allocate memory.")
    END IF
    ! this array contains the boundary condition for each primitive variable;
    ! the default setting is NO_GRADIENTS; the user has to assign reasonable
    ! values after initialization of the boundary module, e.g. set
    ! this%cbtype(1..Physics%VNUM) = {CUSTOM_NOGRAD | CUSTOM_PERIOD | ...}
    ! for each physical variable at each custom boundary
    this%cbtype(:,:) = CUSTOM_NOGRAD
  END SUBROUTINE InitBoundary_custom


  PURE SUBROUTINE CenterBoundary_custom(this,Mesh,Physics,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL :: pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)
    !------------------------------------------------------------------------!
    INTEGER       :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,Mesh,Physics
    INTENT(INOUT) :: pvar  
    !------------------------------------------------------------------------!
    ! ATTENTION: If this%bctype(:) contains bogus values, this boundary module
    !            behaves like the NO_GRADIENTS boundary condition.
!CDIR IEXPAND
    SELECT CASE(GetDirection(this))
    CASE(WEST)
       DO i=1,Mesh%GNUM
          DO j=Mesh%JMIN,Mesh%JMAX
             WHERE (this%cbtype(j,:).EQ.CUSTOM_PERIOD)
                pvar(Mesh%IMIN-i,j,:) = pvar(Mesh%IMAX-i+1,j,:)
             ELSEWHERE(this%cbtype(j,:).EQ.CUSTOM_REFLECT)
                pvar(Mesh%IMIN-i,j,:) = pvar(Mesh%IMIN+i-1,j,:)
             ELSEWHERE(this%cbtype(j,:).EQ.CUSTOM_REFLNEG)
                pvar(Mesh%IMIN-i,j,:) = -pvar(Mesh%IMIN+i-1,j,:)
             ELSEWHERE(this%cbtype(j,:).EQ.CUSTOM_EXTRAPOL)
                pvar(Mesh%IMIN-i,j,:) = (i+1)*pvar(Mesh%IMIN,j,:) - i*pvar(Mesh%IMIN+1,j,:)
             ELSEWHERE(this%cbtype(j,:).EQ.CUSTOM_FIXED)
                pvar(Mesh%IMIN-i,j,:) = this%data(i,j,:)
             ELSEWHERE
                ! defaults to NO_GRADIENTS
                pvar(Mesh%IMIN-i,j,:) = pvar(Mesh%IMIN,j,:)               
             END WHERE
          END DO
       END DO
    CASE(EAST)
       DO i=1,Mesh%GNUM
          DO j=Mesh%JMIN,Mesh%JMAX
             WHERE (this%cbtype(j,:).EQ.CUSTOM_PERIOD)
                pvar(Mesh%IMAX+i,j,:) = pvar(Mesh%IMIN+i-1,j,:)
             ELSEWHERE (this%cbtype(j,:).EQ.CUSTOM_REFLECT)
                pvar(Mesh%IMAX+i,j,:) = pvar(Mesh%IMAX-i+1,j,:)
             ELSEWHERE (this%cbtype(j,:).EQ.CUSTOM_REFLNEG)
                pvar(Mesh%IMAX+i,j,:) = -pvar(Mesh%IMAX-i+1,j,:)
             ELSEWHERE (this%cbtype(j,:).EQ.CUSTOM_EXTRAPOL)
                pvar(Mesh%IMAX+i,j,:) = (i+1)*pvar(Mesh%IMAX,j,:) - i*pvar(Mesh%IMAX-1,j,:)
             ELSEWHERE (this%cbtype(j,:).EQ.CUSTOM_FIXED)
                pvar(Mesh%IMAX+i,j,:) = this%data(i,j,:)
             ELSEWHERE
                ! defaults to NO_GRADIENTS
                pvar(Mesh%IMAX+i,j,:) = pvar(Mesh%IMAX,j,:)               
             END WHERE
          END DO
       END DO
    CASE(SOUTH)
       DO j=1,Mesh%GNUM
          DO i=Mesh%IMIN,Mesh%IMAX
             WHERE (this%cbtype(i,:).EQ.CUSTOM_PERIOD)
                pvar(i,Mesh%JMIN-j,:) = pvar(i,Mesh%JMAX-j+1,:)
             ELSEWHERE (this%cbtype(i,:).EQ.CUSTOM_REFLECT)
                pvar(i,Mesh%JMIN-j,:) = pvar(i,Mesh%JMIN+j-1,:)
             ELSEWHERE (this%cbtype(i,:).EQ.CUSTOM_REFLNEG)
                pvar(i,Mesh%JMIN-j,:) = -pvar(i,Mesh%JMIN+j-1,:)
             ELSEWHERE (this%cbtype(i,:).EQ.CUSTOM_EXTRAPOL)
                pvar(i,Mesh%JMIN-j,:) = (j+1)*pvar(i,Mesh%JMIN,:) - j*pvar(i,Mesh%JMIN+1,:)
             ELSEWHERE (this%cbtype(i,:).EQ.CUSTOM_FIXED)
                pvar(i,Mesh%JMIN-j,:) = this%data(i,j,:)
             ELSEWHERE
                ! defaults to NO_GRADIENTS
                pvar(i,Mesh%JMIN-j,:) = pvar(i,Mesh%JMIN,:)
             END WHERE
          END DO
       END DO
    CASE(NORTH)
       DO j=1,Mesh%GNUM
           DO i=Mesh%IMIN,Mesh%IMAX
              WHERE (this%cbtype(i,:).EQ.CUSTOM_PERIOD)
                 pvar(i,Mesh%JMAX+j,:) = pvar(i,Mesh%JMIN+j-1,:)
              ELSEWHERE (this%cbtype(i,:).EQ.CUSTOM_REFLECT)
                 pvar(i,Mesh%JMAX+j,:) = pvar(i,Mesh%JMAX-j+1,:)
              ELSEWHERE (this%cbtype(i,:).EQ.CUSTOM_REFLNEG)
                 pvar(i,Mesh%JMAX+j,:) = -pvar(i,Mesh%JMAX-j+1,:)
              ELSEWHERE (this%cbtype(i,:).EQ.CUSTOM_EXTRAPOL)
                 pvar(i,Mesh%JMAX+j,:) = (j+1)*pvar(i,Mesh%JMAX,:) - j*pvar(i,Mesh%JMAX-1,:)
              ELSEWHERE (this%cbtype(i,:).EQ.CUSTOM_FIXED)
                 pvar(i,Mesh%JMAX+j,:) = this%data(i,j,:)
              ELSEWHERE
                 ! defaults to NO_GRADIENTS
                 pvar(i,Mesh%JMAX+j,:) = pvar(i,Mesh%JMAX,:)
              END WHERE
           END DO
        END DO
    END SELECT
  END SUBROUTINE CenterBoundary_custom


  SUBROUTINE CloseBoundary_custom(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL CloseBoundary_fixed(this)
    DEALLOCATE(this%cbtype)
  END SUBROUTINE CloseBoundary_custom

END MODULE boundary_custom
