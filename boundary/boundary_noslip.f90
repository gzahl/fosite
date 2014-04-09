!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_noslip.f90                                               #
!#                                                                           #
!# Copyright (C) 2006-2009                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Bjoern Sperling  <sperling@astrophysik.uni-kiel.de>                       #
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
! boundary module for noslip 
!----------------------------------------------------------------------------!
MODULE boundary_noslip
  USE mesh_common, ONLY : Mesh_TYP
  USE fluxes_common, ONLY : Fluxes_TYP
  USE reconstruction_common, ONLY : Reconstruction_TYP, PrimRecon
  USE boundary_nogradients
  USE physics_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "noslip and heating"  
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Boundary_TYP, &
       ! constants
       WEST, EAST, SOUTH, NORTH, &
       ! methods
       InitBoundary_noslip, &
       CenterBoundary_noslip, &
       CloseBoundary_noslip, &
       GetType, &
       GetName, &
       GetDirection, &
       GetDirectionName, &
       GetRank, &
       GetNumProcs, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitBoundary_noslip(this,Mesh,Physics,btype,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: btype,dir
    !------------------------------------------------------------------------!
    INTEGER       :: err
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,btype,dir
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL InitBoundary(this,btype,boundcond_name,dir)

    ! allocate memory for boundary data
    SELECT CASE(dir)
    CASE(WEST,EAST)
       ALLOCATE(this%data(Mesh%GNUM,Mesh%JMIN:Mesh%JMAX,Physics%VNUM), &
            STAT=err)
       this%data(:,:,:)=0.
    CASE(SOUTH,NORTH)
       ALLOCATE(this%data(Mesh%IMIN:Mesh%IMAX,Mesh%GNUM,Physics%VNUM), &
            STAT=err)
       this%data(:,:,:)=0.
    END SELECT
    
    IF (err.NE.0) THEN
       CALL Error(this, "InitBoundary_noslip", "Unable to allocate memory.")
    END IF
  END SUBROUTINE InitBoundary_noslip

  PURE SUBROUTINE CenterBoundary_noslip(this,Mesh,Physics,Fluxes,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Fluxes_TYP)   :: Fluxes
    TYPE(Physics_TYP)  :: Physics
    REAL :: pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)
    !------------------------------------------------------------------------!
    INTEGER       :: i,j
    INTEGER       :: TEMPERATURE
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,Mesh,Physics,Fluxes
    INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!

    ! boundary temperature is saved in data(Physics%PRESSURE)
    TEMPERATURE = Physics%PRESSURE

    SELECT CASE(GetDirection(this))
    CASE(WEST)
       ! UNROLL=Mesh%GNUM would be sufficient, but the compiler does
       ! not know the value of Mesh%GNUM, hence we set UNROLL=4 and
       ! hope that nobody sets Mesh%GNUM to a value greater than 4
!CDIR UNROLL=4
       DO i=1,Mesh%GNUM
          ! vanishing density gradient at the boundary
          pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) &
               = pvar(Mesh%IMIN+i-1,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY)
          ! normal velocity
          pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY) &
               = -pvar(Mesh%IMIN+i-1,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY)
          ! tangential velocities
          pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY) = &
               - pvar(Mesh%IMIN+i-1,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY) &
               + 2.0*this%data(1,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY)          
          IF (Physics%ZVELOCITY /= 0) THEN
             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) = &
                  -pvar(Mesh%IMIN+i-1,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY)&
                  +2.0*this%data(1,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY)
          END IF
          ! pressure or temperature (for non-isothermal simulations)
          IF (Physics%PRESSURE /= 0) THEN
             WHERE (this%data(1,Mesh%JMIN:Mesh%JMAX,TEMPERATURE) > 0.0)
                ! set pressure using fixed temperature data
                pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) = &
                     this%data(1,Mesh%JMIN:Mesh%JMAX,TEMPERATURE) * &
                     pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) &
                     * Physics%constants%RG / Physics%mu
             ELSEWHERE
                ! vanishing pressure gradient at the boundary
                pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) = &
                     pvar(Mesh%IMIN+i-1,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE)
             END WHERE
          END IF
       END DO
    CASE(EAST)
!CDIR UNROLL=4
       DO i=1,Mesh%GNUM
          ! vanishing density gradient at the boundary
          pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) = &
               pvar(Mesh%IMAX-i+1,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY)
          ! normal velocity
          pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY) &
               = -pvar(Mesh%IMAX-i+1,Mesh%JMIN:Mesh%JMAX,Physics%XVELOCITY)
          ! tangential velocities
          pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY) = &
               - pvar(Mesh%IMAX-i+1,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY) &
               + 2.0*this%data(1,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY)          
          IF (Physics%ZVELOCITY /= 0) THEN
             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) = &
                  - pvar(Mesh%IMAX-i+1,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) &
                  + 2.0*this%data(1,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY)          
          END IF
          ! pressure or temperature (for non-isothermal simulations)
          IF (Physics%PRESSURE /= 0) THEN
             WHERE (this%data(1,Mesh%JMIN:Mesh%JMAX,TEMPERATURE) > 0.0)
                ! set pressure using fixed temperature data
                pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) = &
                     this%data(1,Mesh%JMIN:Mesh%JMAX,TEMPERATURE) * &
                     pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) &
                     * Physics%constants%RG / Physics%mu
             ELSEWHERE
                ! vanishing pressure gradient at the boundary
                pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE) = &
                     pvar(Mesh%IMAX-i+1,Mesh%JMIN:Mesh%JMAX,Physics%PRESSURE)
             END WHERE
          END IF
       END DO
    CASE(SOUTH)
!CDIR UNROLL=4
       DO j=1,Mesh%GNUM
          ! vanishing density gradient at the boundary
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%DENSITY) &
               = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+j-1,Physics%DENSITY)
          ! normal velocity
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%YVELOCITY) &
               = -pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+j-1,Physics%YVELOCITY)
          ! tangential velocities
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%XVELOCITY) = &
               - pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+j-1,Physics%XVELOCITY) &
               + 2.0*this%data(Mesh%IMIN:Mesh%IMAX,1,Physics%XVELOCITY)
          IF (Physics%ZVELOCITY /= 0) THEN
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%ZVELOCITY) = &
                  - pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+j-1,Physics%ZVELOCITY) &
                  + 2.0*this%data(Mesh%IMIN:Mesh%IMAX,1,Physics%ZVELOCITY)
          END IF
          ! pressure or temperature (for non-isothermal simulations)
          IF (Physics%PRESSURE /= 0) THEN
             WHERE (this%data(Mesh%IMIN:Mesh%IMAX,1,TEMPERATURE) > 0.0)
                ! set pressure using fixed temperature data
                pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%PRESSURE) = &
                     this%data(Mesh%IMIN:Mesh%IMAX,1,TEMPERATURE) * &
                     pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%DENSITY) &
                     * Physics%constants%RG / Physics%mu
             ELSEWHERE
                ! vanishing pressure gradient at the boundary
                pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Physics%PRESSURE) &
                     = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+j-1,Physics%PRESSURE)
             END WHERE
          END IF
       END DO
    CASE(NORTH)
!CDIR UNROLL=4
       DO j=1,Mesh%GNUM
          ! vanishing density gradient at the boundary
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%DENSITY) &
               = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-j+1,Physics%DENSITY)
          ! normal velocity
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%YVELOCITY) &
               = -pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-j+1,Physics%YVELOCITY)
          ! tangential velocities
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%XVELOCITY) = &
               - pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-j+1,Physics%XVELOCITY) &
               + 2.0*this%data(Mesh%IMIN:Mesh%IMAX,1,Physics%XVELOCITY)
          IF (Physics%ZVELOCITY /= 0) THEN
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%ZVELOCITY) = &
                  - pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-j+1,Physics%ZVELOCITY) &
                  + 2.0*this%data(Mesh%IMIN:Mesh%IMAX,1,Physics%ZVELOCITY)
          END IF
          ! pressure or temperature (for non-isothermal simulations)
          IF (Physics%PRESSURE /= 0) THEN
             WHERE (this%data(Mesh%IMIN:Mesh%IMAX,1,TEMPERATURE) > 0.0)
                ! set pressure using fixed temperature data
                pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%PRESSURE) = &
                     this%data(Mesh%IMIN:Mesh%IMAX,1,TEMPERATURE) * &
                     pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%DENSITY) &
                     * Physics%constants%RG / Physics%mu
             ELSEWHERE
                ! vanishing pressure gradient at the boundary
                pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Physics%PRESSURE) &
                     = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-j+1,Physics%PRESSURE)
             END WHERE
          END IF
       END DO
    END SELECT
  END SUBROUTINE CenterBoundary_noslip


  SUBROUTINE CloseBoundary_noslip(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%data)
  END SUBROUTINE CloseBoundary_noslip

END MODULE boundary_noslip