!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_custom.f90                                               #
!#                                                                           #
!# Copyright (C) 2006-2010                                                   #
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
  USE fluxes_common, ONLY : Fluxes_TYP
  USE reconstruction_common, ONLY : Reconstruction_TYP, PrimRecon
  USE boundary_nogradients
  USE boundary_fixed
  USE physics_generic
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
       ALLOCATE(this%data(Mesh%GNUM,Mesh%JMIN:Mesh%JMAX,Physics%VNUM), &
            this%cbtype(Mesh%JMIN:Mesh%JMAX,Physics%VNUM), &
            STAT=err)
    CASE(SOUTH,NORTH)
       ALLOCATE(this%data(Mesh%IMIN:Mesh%IMAX,Mesh%GNUM,Physics%VNUM), &
            this%cbtype(Mesh%IMIN:Mesh%IMAX,Physics%VNUM), &
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


  PURE SUBROUTINE CenterBoundary_custom(this,Mesh,Physics,Fluxes,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    REAL :: pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)
    !------------------------------------------------------------------------!
    INTEGER       :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,Mesh,Physics,Fluxes
    INTENT(INOUT) :: pvar  
    !------------------------------------------------------------------------!
    ! ATTENTION: If this%bctype(:) contains bogus values, this boundary module
    !            behaves like the NO_GRADIENTS boundary condition.
    SELECT CASE(GetDirection(this))
    CASE(WEST)
       DO i=1,Mesh%GNUM
          WHERE (this%cbtype(:,:).EQ.CUSTOM_PERIOD)
             pvar(Mesh%IMIN-i,:,:) = pvar(Mesh%IMAX-i+1,:,:)
          ELSEWHERE(this%cbtype(:,:).EQ.CUSTOM_REFLECT)
             pvar(Mesh%IMIN-i,:,:) = pvar(Mesh%IMIN+i-1,:,:)
          ELSEWHERE(this%cbtype(:,:).EQ.CUSTOM_REFLNEG)
             pvar(Mesh%IMIN-i,:,:) = -pvar(Mesh%IMIN+i-1,:,:)
          ELSEWHERE(this%cbtype(:,:).EQ.CUSTOM_EXTRAPOL)
             pvar(Mesh%IMIN-i,:,:) = (i+1)*pvar(Mesh%IMIN,:,:) - i*pvar(Mesh%IMIN+1,:,:)
          ELSEWHERE(this%cbtype(:,:).EQ.CUSTOM_FIXED)
             pvar(Mesh%IMIN-i,:,:) = this%data(i,:,:)
          ELSEWHERE
             ! defaults to NO_GRADIENTS
             pvar(Mesh%IMIN-i,:,:) = pvar(Mesh%IMIN,:,:)               
          END WHERE
       END DO
    CASE(EAST)
       DO i=1,Mesh%GNUM
          WHERE (this%cbtype(:,:).EQ.CUSTOM_PERIOD)
             pvar(Mesh%IMAX+i,:,:) = pvar(Mesh%IMIN+i-1,:,:)
          ELSEWHERE (this%cbtype(:,:).EQ.CUSTOM_REFLECT)
             pvar(Mesh%IMAX+i,:,:) = pvar(Mesh%IMAX-i+1,:,:)
          ELSEWHERE (this%cbtype(:,:).EQ.CUSTOM_REFLNEG)
             pvar(Mesh%IMAX+i,:,:) = -pvar(Mesh%IMAX-i+1,:,:)
          ELSEWHERE (this%cbtype(:,:).EQ.CUSTOM_EXTRAPOL)
             pvar(Mesh%IMAX+i,:,:) = (i+1)*pvar(Mesh%IMAX,:,:) - i*pvar(Mesh%IMAX-1,:,:)
          ELSEWHERE (this%cbtype(:,:).EQ.CUSTOM_FIXED)
             pvar(Mesh%IMAX+i,:,:) = this%data(i,:,:)
          ELSEWHERE
             ! defaults to NO_GRADIENTS
             pvar(Mesh%IMAX+i,:,:) = pvar(Mesh%IMAX,:,:)               
          END WHERE
       END DO
    CASE(SOUTH)
       DO j=1,Mesh%GNUM
          WHERE (this%cbtype(:,:).EQ.CUSTOM_PERIOD)
             pvar(:,Mesh%JMIN-j,:) = pvar(:,Mesh%JMAX-j+1,:)
          ELSEWHERE (this%cbtype(:,:).EQ.CUSTOM_REFLECT)
             pvar(:,Mesh%JMIN-j,:) = pvar(:,Mesh%JMIN+j-1,:)
          ELSEWHERE (this%cbtype(:,:).EQ.CUSTOM_REFLNEG)
             pvar(:,Mesh%JMIN-j,:) = -pvar(:,Mesh%JMIN+j-1,:)
          ELSEWHERE (this%cbtype(:,:).EQ.CUSTOM_EXTRAPOL)
             pvar(:,Mesh%JMIN-j,:) = (j+1)*pvar(:,Mesh%JMIN,:) - j*pvar(:,Mesh%JMIN+1,:)
          ELSEWHERE (this%cbtype(:,:).EQ.CUSTOM_FIXED)
             pvar(:,Mesh%JMIN-j,:) = this%data(:,j,:)
          ELSEWHERE
             ! defaults to NO_GRADIENTS
             pvar(:,Mesh%JMIN-j,:) = pvar(:,Mesh%JMIN,:)
          END WHERE
       END DO
    CASE(NORTH)
       DO j=1,Mesh%GNUM
          WHERE (this%cbtype(:,:).EQ.CUSTOM_PERIOD)
             pvar(:,Mesh%JMAX+j,:) = pvar(:,Mesh%JMIN+j-1,:)
          ELSEWHERE (this%cbtype(:,:).EQ.CUSTOM_REFLECT)
             pvar(:,Mesh%JMAX+j,:) = pvar(:,Mesh%JMAX-j+1,:)
          ELSEWHERE (this%cbtype(:,:).EQ.CUSTOM_REFLNEG)
             pvar(:,Mesh%JMAX+j,:) = -pvar(:,Mesh%JMAX-j+1,:)
          ELSEWHERE (this%cbtype(:,:).EQ.CUSTOM_EXTRAPOL)
             pvar(:,Mesh%JMAX+j,:) = (j+1)*pvar(:,Mesh%JMAX,:) - j*pvar(:,Mesh%JMAX-1,:)
          ELSEWHERE (this%cbtype(:,:).EQ.CUSTOM_FIXED)
             pvar(:,Mesh%JMAX+j,:) = this%data(:,j,:)
          ELSEWHERE
             ! defaults to NO_GRADIENTS
             pvar(:,Mesh%JMAX+j,:) = pvar(:,Mesh%JMAX,:)
          END WHERE
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
    DEALLOCATE(this%data,this%cbtype)
  END SUBROUTINE CloseBoundary_custom

END MODULE boundary_custom
