!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_custom.f90                                               #
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
!! \brief Boundary module for custom conditions
!! 
!! \extends boundary_fixed
!! \ingroup boundary
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
  INTEGER, PARAMETER :: CUSTOM_LOGEXPOL = 7 ! extrapolation of log values
  INTEGER, PARAMETER :: CUSTOM_OUTFLOW  = 8
  INTEGER, PARAMETER :: CUSTOM_KEPLER   = 9
  INTEGER, PARAMETER :: CUSTOM_ANGKEPLER= 10
  INTEGER, PARAMETER :: CUSTOM_POISSON  = 11
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! constants
       CUSTOM_NOGRAD, CUSTOM_PERIOD, CUSTOM_REFLECT, CUSTOM_REFLNEG, &
       CUSTOM_EXTRAPOL, CUSTOM_FIXED, CUSTOM_LOGEXPOL, &
       CUSTOM_OUTFLOW, CUSTOM_KEPLER, CUSTOM_ANGKEPLER, CUSTOM_POISSON, &
       ! methods
       InitBoundary_custom, &
       CenterBoundary_custom, &
       CloseBoundary_custom
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor for custom boundary conditions
  SUBROUTINE InitBoundary_custom(this,Mesh,Physics,btype,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: btype,dir
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,err = 0
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,btype,dir
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL InitBoundary_fixed(this,Mesh,Physics,btype,dir,boundcond_name)
    ! allocate memory for boundary data and mask
    SELECT CASE(GetDirection(this))
    CASE(WEST,EAST)
       ALLOCATE(this%cbtype(Mesh%JMIN:Mesh%JMAX,Physics%VNUM), &
            this%Rscale(Mesh%GNUM,Mesh%JMIN:Mesh%JMAX), &
            this%invRscale(Mesh%GNUM,Mesh%JMIN:Mesh%JMAX), &
            STAT=err)
    CASE(SOUTH,NORTH)
       ALLOCATE(this%cbtype(Mesh%IMIN:Mesh%IMAX,Physics%VNUM), &
            this%Rscale(Mesh%IMIN:Mesh%IMAX,Mesh%GNUM), &
            this%invRscale(Mesh%IMIN:Mesh%IMAX,Mesh%GNUM), &
            STAT=err)
    END SELECT
    IF (err.NE.0) THEN
       CALL Error(this,"InitBoundary_custom", "Unable to allocate memory.")
    END IF

    SELECT CASE(GetDirection(this))
    CASE(WEST)
       DO i=1,Mesh%GNUM
          DO j=Mesh%JMIN,Mesh%JMAX
             this%Rscale(i,j) = Mesh%bradius(Mesh%IMIN-i,j) / Mesh%bradius(Mesh%IMIN,j)
          END DO
       END DO
    CASE(EAST)
       DO i=1,Mesh%GNUM
          DO j=Mesh%JMIN,Mesh%JMAX
             this%Rscale(i,j) = Mesh%bradius(Mesh%IMAX+i,j) / Mesh%bradius(Mesh%IMAX,j)
          END DO
       END DO     
    CASE(SOUTH)
       DO i= Mesh%IMIN,Mesh%IMAX
          DO j=1,Mesh%GNUM
             this%Rscale(i,j) = Mesh%bradius(i,Mesh%JMIN-j) / Mesh%bradius(i,Mesh%JMIN)
          END DO
       END DO
    CASE(NORTH)
       DO i= Mesh%IMIN,Mesh%IMAX
          DO j=1,Mesh%GNUM
             this%Rscale(i,j) = Mesh%bradius(i,Mesh%JMAX+j) / Mesh%bradius(i,Mesh%JMAX)
          END DO
        END DO
    END SELECT
    this%Rscale(:,:) = sqrt(this%Rscale(:,:))
    this%invRscale(:,:) = 1.0 / (this%Rscale(:,:) + TINY(1.0))
    ! this array contains the boundary condition for each primitive variable;
    ! the default setting is NO_GRADIENTS; the user has to assign reasonable
    ! values after initialization of the boundary module, e.g. set
    ! this%cbtype(1..Physics%VNUM) = {CUSTOM_NOGRAD | CUSTOM_PERIOD | ...}
    ! for each physical variable at each custom boundary
    this%cbtype(:,:) = CUSTOM_NOGRAD
  END SUBROUTINE InitBoundary_custom


  !> \public Applies the custom boundary conditions
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
             ELSEWHERE(this%cbtype(j,:).EQ.CUSTOM_LOGEXPOL)
                pvar(Mesh%IMIN-i,j,:) = pvar(Mesh%IMIN,j,:) &
                   * (pvar(Mesh%IMIN,j,:) / pvar(Mesh%IMIN+1,j,:))**(i+1)
             ELSEWHERE(this%cbtype(j,:).EQ.CUSTOM_OUTFLOW &
                     .AND. pvar(Mesh%IMIN,j,:) .GE. 0.0 )
                !REFLNEG, else default (NO_GRADIENTS)
                pvar(Mesh%IMIN-i,j,:) = -pvar(Mesh%IMIN+i-1,j,:)
             ELSEWHERE(this%cbtype(j,:).EQ.CUSTOM_KEPLER)
                pvar(Mesh%IMIN-i,j,:) = (pvar(Mesh%IMIN,j,:) + Mesh%bradius(Mesh%IMIN,j)*Physics%Omega)&
                   * this%invRscale(i,j) - Mesh%bradius(Mesh%IMIN-i,j)*Physics%Omega
             ELSEWHERE(this%cbtype(j,:).EQ.CUSTOM_ANGKEPLER)
                pvar(Mesh%IMIN-i,j,:) = pvar(Mesh%IMIN,j,:) &
                   * this%Rscale(i,j) 
!             ELSEWHERE(this%cbtype(j,:).EQ.CUSTOM_POISSON.AND.this%accel(Mesh%IMIN,j,1).LT.TINY(1.e0))
!                 pvar(Mesh%IMIN-i,j,:) = pvar(Mesh%IMIN,j,:) &
!                   + (Mesh%bhy(Mesh%IMIN-i,j)-Mesh%bhy(Mesh%IMIN,j)) &
!                     * 0.5*(SQRT(Mesh%bhy(Mesh%IMIN+1,j)*(-this%accel(Mesh%IMIN+1,j,1))) &
!                        -SQRT(Mesh%bhy(Mesh%IMIN-1,j)*(-this%accel(Mesh%IMIN-1,j,1))))/Mesh%dlx(Mesh%IMIN,j)
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
             ELSEWHERE (this%cbtype(j,:).EQ.CUSTOM_LOGEXPOL)
                pvar(Mesh%IMAX+i,j,:) = pvar(Mesh%IMAX,j,:) &
                   * (pvar(Mesh%IMAX,j,:) / pvar(Mesh%IMAX-1,j,:))**(i+1)
             ELSEWHERE(this%cbtype(j,:).EQ.CUSTOM_OUTFLOW &
                     .AND. pvar(Mesh%IMAX,j,:) .LE. 0.0 )
                !REFLNEG, else default (NO_GRADIENTS)
                pvar(Mesh%IMAX+i,j,:) = -pvar(Mesh%IMAX-i+1,j,:)
             ELSEWHERE(this%cbtype(j,:).EQ.CUSTOM_KEPLER)
                pvar(Mesh%IMAX+i,j,:) = (pvar(Mesh%IMAX,j,:) + Mesh%bradius(Mesh%IMAX,j)*Physics%Omega)&
                    * this%invRscale(i,j) - Mesh%bradius(Mesh%IMAX+i,j)*Physics%Omega
             ELSEWHERE(this%cbtype(j,:).EQ.CUSTOM_ANGKEPLER)
                pvar(Mesh%IMAX+i,j,:) = pvar(Mesh%IMAX,j,:) &
                   * this%Rscale(i,j)
!             ELSEWHERE(this%cbtype(j,:).EQ.CUSTOM_POISSON.AND.this%accel(Mesh%IMAX,j,1).LT.TINY(1.e0))
!                pvar(Mesh%IMAX+i,j,:) = pvar(Mesh%IMAX,j,:) &
!                   + (Mesh%bhy(Mesh%IMAX+i,j)-Mesh%bhy(Mesh%IMAX,j)) &
!                     * 0.5*(SQRT(Mesh%bhy(Mesh%IMAX+1,j)*(-this%accel(Mesh%IMAX+1,j,1))) &
!                        -SQRT(Mesh%bhy(Mesh%IMAX-1,j)*(-this%accel(Mesh%IMAX-1,j,1))))/Mesh%dlx(Mesh%IMAX,j)
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
             ELSEWHERE (this%cbtype(i,:).EQ.CUSTOM_LOGEXPOL)
                pvar(i,Mesh%JMIN-j,:) = pvar(i,Mesh%JMIN,:) &
                   * (pvar(i,Mesh%JMIN,:) / pvar(i,Mesh%JMIN+1,:))**(j+1)
             ELSEWHERE(this%cbtype(i,:).EQ.CUSTOM_OUTFLOW &
                     .AND. pvar(i,Mesh%JMIN,:) .GE. 0.0 )
                !REFLNEG, else default (NO_GRADIENTS)
                pvar(i,Mesh%JMIN-j,:) = -pvar(i,Mesh%JMIN+j-1,:)
             ELSEWHERE(this%cbtype(i,:).EQ.CUSTOM_KEPLER)
                pvar(i,Mesh%JMIN-j,:) = pvar(i,Mesh%JMIN,:) &
                   * this%invRscale(i,j)
             ELSEWHERE(this%cbtype(i,:).EQ.CUSTOM_ANGKEPLER)
                pvar(i,Mesh%JMIN-j,:) = pvar(i,Mesh%JMIN,:) &
                   * this%Rscale(i,j)
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
             ELSEWHERE (this%cbtype(i,:).EQ.CUSTOM_LOGEXPOL)
                 pvar(i,Mesh%JMAX+j,:) = pvar(i,Mesh%JMAX,:) &
                    * (pvar(i,Mesh%JMAX,:) / pvar(i,Mesh%JMAX-1,:))**(j+1)
             ELSEWHERE(this%cbtype(i,:).EQ.CUSTOM_OUTFLOW &
                     .AND. pvar(i,Mesh%JMAX,:) .LE. 0.0 )
                !REFLNEG, else default (NO_GRADIENTS)
                pvar(i,Mesh%JMAX+j,:) = -pvar(i,Mesh%JMAX-j+1,:)
             ELSEWHERE(this%cbtype(i,:).EQ.CUSTOM_KEPLER)
                pvar(i,Mesh%JMAX+j,:) = pvar(i,Mesh%JMAX,:) &
                   * this%invRscale(i,j)
             ELSEWHERE(this%cbtype(i,:).EQ.CUSTOM_ANGKEPLER)
                pvar(i,Mesh%JMAX+j,:) = pvar(i,Mesh%JMAX,:) &
                   * this%Rscale(i,j)
             ELSEWHERE
                 ! defaults to NO_GRADIENTS
                 pvar(i,Mesh%JMAX+j,:) = pvar(i,Mesh%JMAX,:)
             END WHERE
           END DO
        END DO
    END SELECT
  END SUBROUTINE CenterBoundary_custom


  !> \public Destructor for custom boundary conditions
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
