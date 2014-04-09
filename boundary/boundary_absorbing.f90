!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_absorbing.f90                                            #
!#                                                                           #
!# Copyright (C) 2009-2014                                                   #
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
!! \brief Boundary module for absorbing (non-reflecting) conditions
!!
!! This module uses characteristic variables and wave speeds at the
!! boundary to determine the state of the flow. Depending on this it
!! damps oszillations by setting the characteristic variables to zero
!! for imcomming waves.
!!
!! \extends boundary_nogradients
!! \ingroup boundary
!----------------------------------------------------------------------------!
MODULE boundary_absorbing
  USE mesh_common, ONLY : Mesh_TYP
  USE boundary_nogradients
  USE physics_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "absorbing"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! methods
       InitBoundary_absorbing, &
       CenterBoundary_absorbing, &
       CloseBoundary_absorbing
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor for absorbing boundary conditions
  !!
  !! Initilizes the boundary condition type and direction and allocates
  !! memory for characteristic variables and wave speeds.
  SUBROUTINE InitBoundary_absorbing(this,Mesh,Physics,btype,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: btype,dir
    !------------------------------------------------------------------------!
    INTEGER            :: err = 0
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,btype,dir
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL InitBoundary(this,btype,boundcond_name,dir)
    ! check if physics supports absorbing boundary conditions
    IF (.NOT.Physics%supports_absorbing) &
       CALL Error(this,"InitBoundary_absorbing", &
                  "boundary condition not supported for this type of physics")
    SELECT CASE(GetDirection(this))
    CASE(WEST,EAST)
       ALLOCATE(this%xvar(Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
            this%lambda(Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
            STAT=err)
    CASE(SOUTH,NORTH)
       ALLOCATE(this%xvar(Mesh%IGMIN:Mesh%IGMAX,Physics%VNUM), &
            this%lambda(Mesh%IGMIN:Mesh%IGMAX,Physics%VNUM), &
            STAT=err)
    END SELECT
    IF (err.NE.0) & 
         CALL Error(this,"InitBoundary_absorbing", "Unable to allocate memory.")
    ! initialize the data
    this%xvar(:,:) = 0.0
    this%lambda(:,:) = 0.0
  END SUBROUTINE InitBoundary_absorbing


  !> \public Applies the absorbing boundary condition
  !!
  !! This is an implementation of characteristic variable extrapolation.
  !! The algorithm first computes the characteristic (pseudo-) variables at
  !! the boundary and then sets them to zero for incomming waves, i. e.
  !! for those associates with positive (western/southern) or negative
  !! (eastern/northern) wave speeds. After that it transforms the 
  !! new set of characteristic variables back to primitive variables.
  PURE SUBROUTINE CenterBoundary_absorbing(this,Mesh,Physics,cvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
         :: cvar,pvar
    !------------------------------------------------------------------------!
    INTEGER       :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics
    INTENT(INOUT) :: this,cvar,pvar 
    !------------------------------------------------------------------------!
    SELECT CASE(GetDirection(this))
    CASE(WEST)
       DO i=1,Mesh%GNUM
          CALL CalculateCharSystemX(Physics,Mesh,Mesh%IMIN-i+1,+1,pvar,this%lambda,this%xvar)
          ! set characteristic variables to zero for all incomming waves
          WHERE (this%lambda(:,:).GE.0.0)
             this%xvar(:,:) = 0.0
          END WHERE
          ! transform back to primitive variables at the boundary
          CALL CalculateBoundaryDataX(Physics,Mesh,Mesh%IMIN-i+1,-1,this%xvar,pvar)
       END DO
    CASE(EAST)
       ! characteristic variables
       DO i=1,Mesh%GNUM
          CALL CalculateCharSystemX(Physics,Mesh,Mesh%IMAX+i-1,-1,pvar,this%lambda,this%xvar)
          ! set characteristic variables to zero for all incomming waves
          WHERE (this%lambda(:,:).LE.0.0)
             this%xvar(:,:) = 0.0
          END WHERE
          ! transform back to primitive variables at the boundary
          CALL CalculateBoundaryDataX(Physics,Mesh,Mesh%IMAX+i-1,+1,this%xvar,pvar)
       END DO
    CASE(SOUTH)
       DO j=1,Mesh%GNUM
          CALL CalculateCharSystemY(Physics,Mesh,Mesh%JMIN-j+1,+1,pvar,this%lambda,this%xvar)
          ! set characteristic variables to zero for all incomming waves
          WHERE (this%lambda(:,:).GE.0.0)
             this%xvar(:,:) = 0.0
          END WHERE
          ! transform back to primitive variables at the boundary
          CALL CalculateBoundaryDataY(Physics,Mesh,Mesh%JMIN-j+1,-1,this%xvar,pvar)
       END DO
    CASE(NORTH)
       DO j=1,Mesh%GNUM
          CALL CalculateCharSystemY(Physics,Mesh,Mesh%JMAX+j-1,-1,pvar,this%lambda,this%xvar)
          ! set characteristic variables to zero for all incomming waves
          WHERE (this%lambda(:,:).LE.0.0)
             this%xvar(:,:) = 0.0
          END WHERE
          ! transform back to primitive variables at the boundary
          CALL CalculateBoundaryDataY(Physics,Mesh,Mesh%JMAX+j-1,+1,this%xvar,pvar)
       END DO
     END SELECT
  END SUBROUTINE CenterBoundary_absorbing


  !> \public Destructor for absorbing boundary conditions
  SUBROUTINE CloseBoundary_absorbing(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%xvar,this%lambda)
  END SUBROUTINE CloseBoundary_absorbing

END MODULE boundary_absorbing
