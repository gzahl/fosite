!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_farfield.f90                                             #
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
! inflow/outflow boundary conditions using Riemann invariants
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

  SUBROUTINE InitBoundary_farfield(this,Mesh,Physics,btype,dir)
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
    ! not working for isothermal Euler equations
    IF ((GetType(Physics).EQ.EULER2D).OR.(GetType(Physics).EQ.EULER2D_ISOTHERM)) &
         CALL Error(this,"InitBoundary_farfield", "Physics module is not " // &
         "supported for this kind of boundary conditions.")
    ! allocate memory for boundary data and mask
!CDIR IEXPAND
    SELECT CASE(GetDirection(this))
    CASE(WEST,EAST)
       ALLOCATE(this%data(Mesh%GNUM,Mesh%JMIN:Mesh%JMAX,Physics%VNUM), &
            this%Rinv(Mesh%GNUM,Mesh%JMIN:Mesh%JMAX,2), &
            this%Rtmp(Mesh%JMIN:Mesh%JMAX,2), &
            this%cs(Mesh%JMIN:Mesh%JMAX), &
            this%cs2gam(Mesh%JMIN:Mesh%JMAX), &
            this%vn(Mesh%JMIN:Mesh%JMAX), &
            this%s(Mesh%JMIN:Mesh%JMAX), &
            STAT=err)
    CASE(SOUTH,NORTH)
       ALLOCATE(this%data(Mesh%IMIN:Mesh%IMAX,Mesh%GNUM,Physics%VNUM), &
            this%Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%GNUM,2), &
            this%Rtmp(Mesh%IMIN:Mesh%IMAX,2), &
            this%cs(Mesh%IMIN:Mesh%IMAX), &
            this%cs2gam(Mesh%IMIN:Mesh%IMAX), &
            this%vn(Mesh%IMIN:Mesh%IMAX), &
            this%s(Mesh%IMIN:Mesh%IMAX), &
            STAT=err)
    END SELECT
    IF (err.NE.0) THEN
       CALL Error(this,"InitBoundary_farfield", "Unable to allocate memory.")
    END IF
    this%first_call = .TRUE.
  END SUBROUTINE InitBoundary_farfield


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
!CDIR IEXPAND
    SELECT CASE(GetDirection(this))
    CASE(WEST)
       IF (this%first_call) THEN
          ! compute Riemann invariant (R+) associated with the incomming wave
          ! using far field data; this has to be done only once per simulation
          this%Rinv(:,:,1) = this%data(:,:,Physics%XVELOCITY)+2./(Physics%gamma-1.0) &
               * GetSoundSpeed_adiabatic(Physics%gamma,this%data(:,:,Physics%DENSITY), &
               this%data(:,:,Physics%PRESSURE))
          this%first_call = .FALSE.
       END IF
       ! compute 2nd Riemann invariant (R-) associated with the outgoing wave
       ! at internal points
!CDIR UNROLL=2
       DO i=0,1
          this%Rtmp(:,i+1) = pvar(Mesh%IMIN+i,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY) &
               - 2./(Physics%gamma-1.0)*GetSoundSpeed_adiabatic(Physics%gamma, &
               pvar(Mesh%IMIN+i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY), &
               pvar(Mesh%IMIN+i,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE))
       END DO
       ! extrapolate R-
       this%Rinv(1,:,2) = 1.5*this%Rtmp(:,1)-0.5*this%Rtmp(:,2)
       DO i=2,Mesh%GNUM
          this%Rinv(i,:,2) = this%Rinv(1,:,2)
       END DO
       ! speed of sound at the boundary
       this%cs(:) = 0.25*(Physics%gamma-1.0)*(this%Rinv(1,:,1)-this%Rinv(1,:,2))
       ! normal velocity at the boundary
       this%vn(:) = 0.5*(this%Rinv(1,:,1)+this%Rinv(1,:,2))
       ! UNROLL=Mesh%GNUM would be sufficient, but the compiler does
       ! not know the value of Mesh%GNUM, hence we set UNROLL=4 and
       ! hope that nobody sets Mesh%GNUM to a value greater than 4
!CDIR UNROLL=4
       DO i=1,Mesh%GNUM
          WHERE (this%vn(:).LT.-this%cs(:))
             ! supersonic outflow (extrapolation)
!!$             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) = &
!!$                  (i+1)*pvar(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) &
!!$                  - i*pvar(Mesh%IMIN+1,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY)
!!$             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY) = &
!!$                  (i+1)*pvar(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY) &
!!$                  - i*pvar(Mesh%IMIN+1,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY)
!!$             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY) = &
!!$                  (i+1)*pvar(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY) &
!!$                  - i*pvar(Mesh%IMIN+1,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY)
!!$             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) = &
!!$                  (i+1)*pvar(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) &
!!$                  - i*pvar(Mesh%IMIN+1,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY)
!!$             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) = &
!!$                  (i+1)*pvar(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) &
!!$                  - i*pvar(Mesh%IMIN+1,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE)
             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) = &
                  pvar(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY)
             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY) = &
                  pvar(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY)
             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY) = &
                  pvar(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY)
             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) = &
                  pvar(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY)
             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) = &
                  pvar(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE)
          ELSEWHERE (this%vn(:).GT.this%cs(:))
             ! supersonic inflow (copy data)
             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) = &
                  this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY)
             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY) = &
                  this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY)
             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY) = &
                  this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY)
             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) = &
                  this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY)
             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) = &
                  this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE)
          ELSEWHERE
             ! subsonic flow
             WHERE (this%vn(:).LT.0.0)
                ! outflow
                ! tangential velocities
                pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY) = &
                     pvar(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY)
                pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) = &
                     pvar(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) &
                     * Mesh%bhz(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX) &
                     / Mesh%bhz(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX)
                ! entropy
                this%s(:) = pvar(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) &
                     /pvar(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY)**Physics%gamma
             ELSEWHERE
                ! inflow
                ! tangential velocities
                pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY) = &
                     this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY)
                pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) = &
                     this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY)
                ! entropy
                this%s(:) = this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) &
                     /this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY)**Physics%gamma
             END WHERE
             ! normal velocity
             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY) = &
                  0.5*(this%Rinv(i,:,1)+this%Rinv(i,:,2))
             ! cs**2 / gamma
             this%cs2gam(:) = (0.25*(Physics%gamma-1.0)*(this%Rinv(i,:,1) &
                  -this%Rinv(i,:,2)))**2 / Physics%gamma
             ! density
             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) = (this%cs2gam(:) &
                  /this%s(:))**(1./(Physics%gamma-1.0))
             ! pressure
             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) = this%cs2gam(:)  &
                  * pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY)
          END WHERE
       END DO
    CASE(EAST)
       IF (this%first_call) THEN
          ! compute Riemann invariant (R-) associated with the incomming wave
          ! using far field data; this has to be done only once per simulation
          this%Rinv(:,:,2) = this%data(:,:,Physics%XVELOCITY)-2./(Physics%gamma-1.0) &
               * GetSoundSpeed_adiabatic(Physics%gamma,this%data(:,:,Physics%DENSITY), &
               this%data(:,:,Physics%PRESSURE))
          this%first_call = .FALSE.
       END IF
       ! compute 1st Riemann invariant (R+) associated with the outgoing wave
       ! at internal points
!CDIR UNROLL=2
       DO i=0,1
          this%Rtmp(:,i+1) = pvar(Mesh%IMAX-i,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY) &
            + 2./(Physics%gamma-1.0)*GetSoundSpeed_adiabatic(Physics%gamma, &
            pvar(Mesh%IMAX-i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY), &
            pvar(Mesh%IMAX-i,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE))
       END DO
       ! extrapolate R+
       this%Rinv(1,:,1) = 1.5*this%Rtmp(:,1)-0.5*this%Rtmp(:,2)
       DO i=2,Mesh%GNUM
          this%Rinv(i,:,1) = this%Rinv(1,:,1)
       END DO
       ! sound speed at the boundary
       this%cs(:) = 0.25*(Physics%gamma-1.0)*(this%Rinv(1,:,1)-this%Rinv(1,:,2))
       ! normal velocity at the boundary
       this%vn(:) = 0.5*(this%Rinv(1,:,1)+this%Rinv(1,:,2))
       ! UNROLL=Mesh%GNUM would be sufficient, but the compiler does
       ! not know the value of Mesh%GNUM, hence we set UNROLL=4 and
       ! hope that nobody sets Mesh%GNUM to a value greater than 4
       ! speed of sound in boundary cells
!CDIR UNROLL=4
       DO i=1,Mesh%GNUM
          WHERE (this%vn(:).GT.this%cs(:))
             ! supersonic outflow (extrapolation)
!!$             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) = &
!!$                  (i+1)*pvar(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) &
!!$                  - i*pvar(Mesh%IMAX-1,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY)
!!$             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY) = &
!!$                  (i+1)*pvar(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY) &
!!$                  - i*pvar(Mesh%IMAX-1,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY)
!!$             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY) = &
!!$                  (i+1)*pvar(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY) &
!!$                  - i*pvar(Mesh%IMAX-1,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY)
!!$             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) = &
!!$                  (i+1)*pvar(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) &
!!$                  - i*pvar(Mesh%IMAX-1,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY)
!!$             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) = &
!!$                  (i+1)*pvar(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) &
!!$                  - i*pvar(Mesh%IMAX-1,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE)
             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) = &
                  pvar(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY)
             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY) = &
                  pvar(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY)
             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY) = &
                  pvar(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY)
             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) = &
                  pvar(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY)
             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) = &
                  pvar(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE)
          ELSEWHERE (this%vn(:).LT.-this%cs(:))
             ! supersonic inflow (copy data)
             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) = &
                  this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY)
             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY) = &
                  this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY)
             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY) = &
                  this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY)
             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) = &
                  this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY)
             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) = &
                  this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE)
          ELSEWHERE
             ! subsonic flow
             WHERE (this%vn(:).GT.0.0)
                ! outflow
                ! tangential velocities
                pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY) = &
                     pvar(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY)
                pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) = &
                     pvar(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) &
                     * Mesh%bhz(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) &
                     / Mesh%bhz(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX)
                ! entropy
                this%s(:) = pvar(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) &
                     /pvar(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY)**Physics%gamma
             ELSEWHERE
                ! inflow
                ! tangential velocities
                pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY) = &
                     this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY)
                pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) = &
                     this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY)
                ! entropy
                this%s(:) = this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) &
                     /this%data(i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY)**Physics%gamma
             END WHERE
             ! normal velocity
             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY) =  &
                  0.5*(this%Rinv(i,:,1)+this%Rinv(i,:,2))
             ! cs**2 / gamma
             this%cs2gam(:) = (0.25*(Physics%gamma-1.0)*(this%Rinv(i,:,1) &
                  -this%Rinv(i,:,2)))**2 / Physics%gamma
             ! density
             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) = (this%cs2gam(:) &
                  /this%s(:))**(1./(Physics%gamma-1.0))
             ! pressure
             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) = this%cs2gam(:)  &
                  * pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY)
          END WHERE
       END DO
    CASE(SOUTH)
       IF (this%first_call) THEN
          ! compute Riemann invariant (R+) associated with the incomming wave
          ! using far field data
          this%Rinv(:,:,1) = this%data(:,:,Physics%YVELOCITY)+2./(Physics%gamma-1.0) &
               * GetSoundSpeed_adiabatic(Physics%gamma,this%data(:,:,Physics%DENSITY), &
               this%data(:,:,Physics%PRESSURE))
          this%first_call = .FALSE.
       END IF
       ! compute 2nd Riemann invariant (R-) associated with the outgoing wave
       ! at internal points
!CDIR UNROLL=2
       DO j=0,1
          this%Rtmp(Mesh%IMIN:Mesh%IMAX,j+1) = &
               pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+j,Physics%YVELOCITY) &
               - 2./(Physics%gamma-1.0)*GetSoundSpeed_adiabatic(Physics%gamma, &
               pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+j,Physics%DENSITY), &
               pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+j,Physics%PRESSURE))
       END DO
       ! extrapolate R-
       this%Rinv(:,1,2) = 1.5*this%Rtmp(:,1)-0.5*this%Rtmp(:,2)
       DO j=2,Mesh%GNUM
          this%Rinv(:,j,2) = this%Rinv(:,1,2)
       END DO
       ! speed of sound at the boundary
       this%cs(:) = 0.25*(Physics%gamma-1.0)*(this%Rinv(:,1,1)-this%Rinv(:,1,2))
       ! normal velocity at the boundary
       this%vn(:) = 0.5*(this%Rinv(:,1,1)+this%Rinv(:,1,2))
       ! UNROLL=Mesh%GNUM would be sufficient, but the compiler does
       ! not know the value of Mesh%GNUM, hence we set UNROLL=4 and
       ! hope that nobody sets Mesh%GNUM to a value greater than 4
!CDIR UNROLL=4
       DO j=1,Mesh%GNUM
          WHERE (this%vn(:).LT.-this%cs(:))
             ! supersonic outflow (extrapolation)
!!$             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%DENSITY) = &
!!$                  (j+1)*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,Physics%DENSITY) &
!!$                  - j*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+1,Physics%DENSITY)
!!$             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%XVELOCITY) = &
!!$                  (j+1)*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,Physics%XVELOCITY) &
!!$                  - j*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+1,Physics%XVELOCITY)
!!$             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%YVELOCITY) = &
!!$                  (j+1)*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,Physics%YVELOCITY) &
!!$                  - j*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+1,Physics%YVELOCITY)
!!$             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%ZVELOCITY) = &
!!$                  (j+1)*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,Physics%ZVELOCITY) &
!!$                  - j*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+1,Physics%ZVELOCITY)
!!$             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%PRESSURE) = &
!!$                  (j+1)*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,Physics%PRESSURE) &
!!$                  - j*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+1,Physics%PRESSURE)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%DENSITY) = &
                  pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,Physics%DENSITY)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%XVELOCITY) = &
                  pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,Physics%XVELOCITY)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%YVELOCITY) = &
                  pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,Physics%YVELOCITY)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%ZVELOCITY) = &
                  pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,Physics%ZVELOCITY)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%PRESSURE) = &
                  pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,Physics%PRESSURE)
          ELSEWHERE (this%vn(:).GT.this%cs(:))
             ! supersonic inflow (copy data)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%DENSITY) = &
                  this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%DENSITY)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%XVELOCITY) = &
                  this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%XVELOCITY)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%YVELOCITY) = &
                  this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%YVELOCITY)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%ZVELOCITY) = &
                  this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%ZVELOCITY)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%PRESSURE) = &
                  this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%PRESSURE)
          ELSEWHERE
             ! subsonic flow
             WHERE (this%vn(:).LT.0.0)
                ! outflow
                ! tangential velocities
                pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%XVELOCITY) = &
                     pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,Physics%XVELOCITY)
                pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%ZVELOCITY) = &
                     pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,Physics%ZVELOCITY) &
                     * Mesh%bhz(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN) &
                     / Mesh%bhz(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j)
                ! entropy
                this%s(:) = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,Physics%PRESSURE) &
                     /pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,Physics%DENSITY)**Physics%gamma
             ELSEWHERE
                ! inflow
                ! tangential velocities
                pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%XVELOCITY) = &
                     this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%XVELOCITY)
                pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%ZVELOCITY) = &
                     this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%ZVELOCITY)
                ! entropy
                this%s(:) = this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%PRESSURE) &
                     /this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%DENSITY)**Physics%gamma
             END WHERE
             ! normal velocity
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%YVELOCITY) = &
                  0.5*(this%Rinv(:,j,1)+this%Rinv(:,j,2))
             ! cs**2 / gamma
             this%cs2gam(:) = (0.25*(Physics%gamma-1.0)*(this%Rinv(:,j,1) &
                  -this%Rinv(:,j,2)))**2 / Physics%gamma
             ! density
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%DENSITY) = (this%cs2gam(:) &
                  /this%s(:))**(1./(Physics%gamma-1.0))
             ! pressure
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%PRESSURE) = this%cs2gam(:)  &
                  * pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%DENSITY)
          END WHERE
       END DO
    CASE(NORTH)
       IF (this%first_call) THEN
          ! compute Riemann invariant (R-) associated with the incomming wave
          ! using far field data
          this%Rinv(:,:,2) = this%data(:,:,Physics%YVELOCITY)-2./(Physics%gamma-1.0) &
               * GetSoundSpeed_adiabatic(Physics%gamma,this%data(:,:,Physics%DENSITY), &
               this%data(:,:,Physics%PRESSURE))
          this%first_call = .FALSE.
       END IF
       ! compute 2nd Riemann invariant (R+) associated with the outgoing wave
       ! at internal points
!CDIR UNROLL=2
       DO j=0,1
          this%Rtmp(:,j+1) = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-j,Physics%YVELOCITY) &
            + 2./(Physics%gamma-1.0)*GetSoundSpeed_adiabatic(Physics%gamma, &
            pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-j,Physics%DENSITY), &
            pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-j,Physics%PRESSURE))
       END DO
       ! extrapolate R+
       this%Rinv(:,1,1) = 1.5*this%Rtmp(:,1)-0.5*this%Rtmp(:,2)
       DO j=2,Mesh%GNUM
          this%Rinv(:,j,1) = this%Rinv(:,1,1)
       END DO
       ! sound speed at the boundary
       this%cs(:) = 0.25*(Physics%gamma-1.0)*(this%Rinv(:,1,1)-this%Rinv(:,1,2))
       ! normal velocity at the boundary
       this%vn(:) = 0.5*(this%Rinv(:,1,1)+this%Rinv(:,1,2))
       ! UNROLL=Mesh%GNUM would be sufficient, but the compiler does
       ! not know the value of Mesh%GNUM, hence we set UNROLL=4 and
       ! hope that nobody sets Mesh%GNUM to a value greater than 4
!CDIR UNROLL=4
       DO j=1,Mesh%GNUM
          WHERE (this%vn(:).GT.this%cs(:))
             ! supersonic outflow (extrapolation)
!!$             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%DENSITY) = &
!!$                  (j+1)*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX,Physics%DENSITY) &
!!$                  - j*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-1,Physics%DENSITY)
!!$             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%XVELOCITY) = &
!!$                  (j+1)*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX,Physics%XVELOCITY) &
!!$                  - j*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-1,Physics%XVELOCITY)
!!$             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%YVELOCITY) = &
!!$                  (j+1)*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX,Physics%YVELOCITY) &
!!$                  - j*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-1,Physics%YVELOCITY)
!!$             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%ZVELOCITY) = &
!!$                  (j+1)*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX,Physics%ZVELOCITY) &
!!$                  - j*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-1,Physics%ZVELOCITY)
!!$             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%PRESSURE) = &
!!$                  (j+1)*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX,Physics%PRESSURE) &
!!$                  - j*pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-1,Physics%PRESSURE)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%DENSITY) = &
                  pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX,Physics%DENSITY)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%XVELOCITY) = &
                  pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX,Physics%XVELOCITY)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%YVELOCITY) = &
                  pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX,Physics%YVELOCITY)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%ZVELOCITY) = &
                  pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX,Physics%ZVELOCITY)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%PRESSURE) = &
                  pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX,Physics%PRESSURE)
          ELSEWHERE (this%vn(:).LT.-this%cs(:))
             ! supersonic inflow (copy data)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%DENSITY) = &
                  this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%DENSITY)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%XVELOCITY) = &
                  this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%XVELOCITY)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%YVELOCITY) = &
                  this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%YVELOCITY)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%ZVELOCITY) = &
                  this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%ZVELOCITY)
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%PRESSURE) = &
                  this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%PRESSURE)
          ELSEWHERE
             ! subsonic flow
             WHERE (this%vn(:).GT.0.0)
                ! outflow
                ! tangential velocities
                pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%XVELOCITY) = &
                     pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX,Physics%XVELOCITY)
                pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%ZVELOCITY) = &
                     pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX,Physics%ZVELOCITY) &
                     * Mesh%bhz(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX) &
                     / Mesh%bhz(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j)
                ! entropy
                this%s(:) = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX,Physics%PRESSURE) &
                     /pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX,Physics%DENSITY)**Physics%gamma
             ELSEWHERE
                ! inflow
                ! tangential velocities
                pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%XVELOCITY) = &
                     this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%XVELOCITY)
                pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%ZVELOCITY) = &
                     this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%ZVELOCITY)
                ! entropy
                this%s(:) = this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%PRESSURE) &
                     /this%data(Mesh%IMIN:Mesh%IMAX,j,Physics%DENSITY)**Physics%gamma
             END WHERE
             ! normal velocity
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%YVELOCITY) = &
                  0.5*(this%Rinv(:,j,1)+this%Rinv(:,j,2))
             ! cs**2 / gamma
             this%cs2gam(:) = (0.25*(Physics%gamma-1.0)*(this%Rinv(:,j,1) &
                  -this%Rinv(:,j,2)))**2 / Physics%gamma
             ! density
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%DENSITY) = (this%cs2gam(:) &
                  /this%s(:))**(1./(Physics%gamma-1.0))
             ! pressure
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%PRESSURE) = this%cs2gam(:)  &
                  * pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%DENSITY)
          END WHERE
       END DO
    END SELECT
  END SUBROUTINE CenterBoundary_farfield


  SUBROUTINE CloseBoundary_farfield(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%data,this%Rinv,this%Rtmp,this%cs,this%cs2gam,this%vn,this%s)
  END SUBROUTINE CloseBoundary_farfield

END MODULE boundary_farfield
