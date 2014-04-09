!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_farfield.f90                                             #
!#                                                                           #
!# Copyright (C) 2006-2014                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
!! \brief Boundary module for far field conditions
!! 
!! Implementation of inflow/outflow boundary conditions using Riemann invariants.
!!
!! \extends boundary_fixed 
!! \ingroup boundary
!----------------------------------------------------------------------------!
MODULE boundary_farfield
  USE mesh_common, ONLY : Mesh_TYP
  USE boundary_nogradients
  USE boundary_fixed
  USE physics_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "far-field in-/ouflow"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! methods
       InitBoundary_farfield, &
       CenterBoundary_farfield, &
       CloseBoundary_farfield
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor for farfield boundary conditions
  SUBROUTINE InitBoundary_farfield(this,Mesh,Physics,btype,dir)
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
    CALL InitBoundary_fixed(this,Mesh,Physics,btype,dir,boundcond_name)
    ! check if physics supports absorbing boundary conditions
    IF (.NOT.Physics%supports_farfield) &
       CALL Error(this,"InitBoundary_farfield", &
                  "boundary condition not supported for this type of physics")
    
    ! allocate memory for boundary data and mask
!CDIR IEXPAND
    SELECT CASE(GetDirection(this))
    CASE(WEST,EAST)
       ALLOCATE(this%data(Mesh%GNUM,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
            this%Rinv(Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
            this%RinvInf(Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
            this%lambda(Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
            STAT=err)
    CASE(SOUTH,NORTH)
       ALLOCATE(this%data(Mesh%IGMIN:Mesh%IGMAX,Mesh%GNUM,Physics%VNUM), &
            this%Rinv(Mesh%IGMIN:Mesh%IGMAX,Physics%VNUM), &
            this%RinvInf(Mesh%IGMIN:Mesh%IGMAX,Physics%VNUM), &
            this%lambda(Mesh%IGMIN:Mesh%IGMAX,Physics%VNUM), &
            STAT=err)
    END SELECT
    IF (err.NE.0) THEN
       CALL Error(this,"InitBoundary_farfield", "Unable to allocate memory.")
    END IF
    this%first_call = .TRUE.
  END SUBROUTINE InitBoundary_farfield


  !> \public Applies the farfield boundary condition
  PURE SUBROUTINE CenterBoundary_farfield(this,Mesh,Physics,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL :: pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)
    !------------------------------------------------------------------------!
    INTEGER            :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics
    INTENT(INOUT) :: this,pvar  
    !------------------------------------------------------------------------!
   SELECT CASE(GetDirection(this))
   CASE(WEST)
     IF (this%first_call) THEN
        ! compute Riemann invariants at boundary (infinity)
        pvar(Mesh%IMIN-1,:,:) = this%data(1,:,:) 
        CALL CalculatePrim2RiemannX(Physics,Mesh,Mesh%IMIN-1,&
                      pvar,this%lambda,this%RinvInf)
        this%first_call = .FALSE.
     END IF

     DO i=1,Mesh%GNUM
       ! compute Riemann invariants
       CALL CalculatePrim2RiemannX(Physics,Mesh,Mesh%IMIN-i+1,&
                                  pvar,this%lambda,this%Rinv)
       ! set infinity Riemanns for inflow 
       WHERE (this%lambda(:,:).GE.0.0)
             this%Rinv(:,:) = this%RinvInf(:,:)
       END WHERE
       ! transform back to primitive variables at the boundary
       CALL CalculateRiemann2PrimX(Physics,Mesh,Mesh%IMIN-i,this%Rinv,pvar) 
     END DO
   CASE(EAST)
     IF (this%first_call) THEN
        ! compute Riemann invariants at boundary (infinity)
        pvar(Mesh%IMAX+1,:,:) = this%data(1,:,:) 
        CALL CalculatePrim2RiemannX(Physics,Mesh,Mesh%IMAX+1,&
                      pvar,this%lambda,this%RinvInf)
        this%first_call = .FALSE.
     END IF

     DO i=1,Mesh%GNUM
       ! compute Riemann invariants
       CALL CalculatePrim2RiemannX(Physics,Mesh,Mesh%IMAX+i-1,&
                                  pvar,this%lambda,this%Rinv)
       ! set infinity Riemanns for inflow 
       WHERE (this%lambda(:,:).LE.0.0)
             this%Rinv(:,:) = this%RinvInf(:,:)
       END WHERE
       ! transform back to primitive variables at the boundary
       CALL CalculateRiemann2PrimX(Physics,Mesh,Mesh%IMAX+i,this%Rinv,pvar) 
     END DO
   CASE(SOUTH)
     IF (this%first_call) THEN
        ! compute Riemann invariants at boundary (infinity)
        pvar(:,Mesh%JMIN-1,:) = this%data(:,1,:) 
        CALL CalculatePrim2RiemannY(Physics,Mesh,Mesh%JMIN-1,&
                      pvar,this%lambda,this%RinvInf)
        this%first_call = .FALSE.
     END IF

     DO j=1,Mesh%GNUM
       ! compute Riemann invariants
       CALL CalculatePrim2RiemannY(Physics,Mesh,Mesh%JMIN-j+1,&
                                  pvar,this%lambda,this%Rinv)
       ! set infinity Riemanns for inflow 
       WHERE (this%lambda(:,:).GE.0.0)
             this%Rinv(:,:) = this%RinvInf(:,:)
       END WHERE
       ! transform back to primitive variables at the boundary
       CALL CalculateRiemann2PrimY(Physics,Mesh,Mesh%JMIN-j,this%Rinv,pvar) 
     END DO
   CASE(NORTH)
     IF (this%first_call) THEN
        ! compute Riemann invariants at boundary (infinity)
        pvar(:,Mesh%JMAX+1,:) = this%data(:,1,:) 
        CALL CalculatePrim2RiemannY(Physics,Mesh,Mesh%JMAX+1,&
                      pvar,this%lambda,this%RinvInf)
        this%first_call = .FALSE.
     END IF

     DO j=1,Mesh%GNUM
       ! compute Riemann invariants
       CALL CalculatePrim2RiemannY(Physics,Mesh,Mesh%JMAX+j-1,&
                                  pvar,this%lambda,this%Rinv)
       ! set infinity Riemanns for inflow 
       WHERE (this%lambda(:,:).LE.0.0)
             this%Rinv(:,:) = this%RinvInf(:,:)
       END WHERE
       ! transform back to primitive variables at the boundary
       CALL CalculateRiemann2PrimY(Physics,Mesh,Mesh%JMAX+j,this%Rinv,pvar) 
     END DO
    END SELECT 
  END SUBROUTINE CenterBoundary_farfield

  !> \public Destructor for farfield boundary conditions
  SUBROUTINE CloseBoundary_farfield(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%data,this%Rinv,this%RinvInf,this%lambda)
  END SUBROUTINE CloseBoundary_farfield

END MODULE boundary_farfield
