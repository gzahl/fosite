!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_generic.f90                                              #
!#                                                                           #
!# Copyright (C) 2006 Tobias Illenseer <tillense@ita.uni-heidelberg.de>      #
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
! Subroutines for boundary conditions
!----------------------------------------------------------------------------!
MODULE boundary_generic
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  USE boundary_nogradients
  USE boundary_periodic
  USE boundary_reflecting, InitBoundary_common => InitBoundary
  USE boundary_axis
  USE boundary_folded
  USE boundary_fixed
  USE boundary_extrapolation
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE InitBoundary
     MODULE PROCEDURE InitBoundary_one, InitBoundary_all
  END INTERFACE
  INTERFACE CloseBoundary
     MODULE PROCEDURE CloseBoundary_one, CloseBoundary_all
  END INTERFACE
  !--------------------------------------------------------------------------!
  INTEGER, PARAMETER :: NO_GRADIENTS  = 1
  INTEGER, PARAMETER :: PERIODIC      = 2
  INTEGER, PARAMETER :: REFLECTING    = 3
  INTEGER, PARAMETER :: AXIS          = 4
  INTEGER, PARAMETER :: FOLDED        = 5
  INTEGER, PARAMETER :: FIXED         = 6
  INTEGER, PARAMETER :: EXTRAPOLATION = 7
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Boundary_TYP, &
       ! constants
       WEST, EAST, SOUTH, NORTH, &
       NO_GRADIENTS, PERIODIC, REFLECTING, AXIS, FOLDED, FIXED, EXTRAPOLATION, &
       ! methods
       InitBoundary, &
       GetType, &
       GetName, &
       CenterBoundary, &
       FaceBoundary, &
       CloseBoundary
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitBoundary_one(this,Mesh,Physics,btype,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: btype,dir
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,btype,dir
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    ! set boundary properties
    SELECT CASE(btype)
    CASE(NO_GRADIENTS)
       CALL InitBoundary_nogradients(this,Physics,btype,dir)
    CASE(PERIODIC)
       CALL InitBoundary_periodic(this,Physics,btype,dir)
    CASE(REFLECTING)
       CALL InitBoundary_reflecting(this,Physics,btype,dir)
    CASE(AXIS)
       CALL InitBoundary_axis(this,Physics,btype,dir)
    CASE(FOLDED)
       CALL InitBoundary_folded(this,Mesh,Physics,btype,dir)
    CASE(FIXED)
       CALL InitBoundary_fixed(this,Mesh,Physics,btype,dir)
    CASE(EXTRAPOLATION)
       CALL InitBoundary_extrapolation(this,Physics,btype,dir)
    CASE DEFAULT
       PRINT *, "ERROR in InitBoundary_one: unknown boundary condition"
       STOP
    END SELECT

    ! print some information
    PRINT *, "BOUNDARY-> condition:         ", TRIM(GetDirectionName(this)), &
         " ", TRIM(GetName(this))
  END SUBROUTINE InitBoundary_one


  SUBROUTINE InitBoundary_all(this,Mesh,Physics,btype,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), DIMENSION(4) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: btype,dir
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,btype,dir
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    SELECT CASE(dir)
    CASE(WEST,EAST,SOUTH,NORTH)
       CALL InitBoundary_one(this(dir),Mesh,Physics,btype,dir)
    CASE DEFAULT
       PRINT *, "ERROR in InitBoundary_all: unknown direction" 
       STOP
    END SELECT
  END SUBROUTINE InitBoundary_all


  PURE SUBROUTINE CenterBoundary(this,Mesh,Physics,rvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), DIMENSION(4) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL :: rvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)
    !------------------------------------------------------------------------!
    INTEGER       :: i
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,Mesh,Physics
    INTENT(INOUT) :: rvar   
    !------------------------------------------------------------------------!
    DO i=1,4
       SELECT CASE(GetType(this(i)))
       CASE(NO_GRADIENTS)
          CALL CenterBoundary_nogradients(this(i),Mesh,Physics,rvar)
       CASE(PERIODIC)
          CALL CenterBoundary_periodic(this(i),Mesh,Physics,rvar)
       CASE(REFLECTING)
          CALL CenterBoundary_reflecting(this(i),Mesh,Physics,rvar)
       CASE(AXIS)
          CALL CenterBoundary_axis(this(i),Mesh,Physics,rvar)
       CASE(FOLDED)
          CALL CenterBoundary_folded(this(i),Mesh,Physics,rvar)
       CASE(FIXED)
          CALL CenterBoundary_fixed(this(i),Mesh,Physics,rvar)
       CASE(EXTRAPOLATION)
          CALL CenterBoundary_extrapolation(this(i),Mesh,Physics,rvar)
        END SELECT
    END DO
  END SUBROUTINE CenterBoundary


  PURE SUBROUTINE FaceBoundary(this,Mesh,Physics,we,ea,so,no,rstates)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), DIMENSION(4) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: we,ea,so,no
    REAL :: rstates(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,Physics%vnum)
    !------------------------------------------------------------------------!
    INTEGER       :: i
    !------------------------------------------------------------------------!
    INTENT(IN)    :: this,Mesh,Physics,we,ea,so,no
    INTENT(INOUT) :: rstates
    !------------------------------------------------------------------------!
    DO i=1,4
       SELECT CASE(GetType(this(i)))
       CASE(NO_GRADIENTS)
          CALL FaceBoundary_nogradients(this(i),Mesh,Physics,we,ea,so,no,rstates)
       CASE(PERIODIC)
          CALL FaceBoundary_periodic(this(i),Mesh,Physics,we,ea,so,no,rstates)
       CASE(REFLECTING)
          CALL FaceBoundary_reflecting(this(i),Mesh,Physics,we,ea,so,no,rstates)
       CASE(AXIS)
          CALL FaceBoundary_axis(this(i),Mesh,Physics,we,ea,so,no,rstates)
       CASE(FOLDED)
          CALL FaceBoundary_folded(this(i),Mesh,Physics,we,ea,so,no,rstates)
       CASE(FIXED)
          CALL FaceBoundary_fixed(this(i),Mesh,Physics,we,ea,so,no,rstates)
       CASE(EXTRAPOLATION)
          CALL FaceBoundary_extrapolation(this(i),Mesh,Physics,we,ea,so,no,rstates)
       END SELECT
    END DO
  END SUBROUTINE FaceBoundary

  SUBROUTINE CloseBoundary_one(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(NO_GRADIENTS,PERIODIC,EXTRAPOLATION)
       ! do nothing
    CASE(REFLECTING)
       CALL CloseBoundary_reflecting(this)
    CASE(AXIS)
       CALL CloseBoundary_axis(this)
    CASE(FOLDED)
       CALL CloseBoundary_folded(this)
    CASE(FIXED)
       CALL CloseBoundary_fixed(this)
    END SELECT
  END SUBROUTINE CloseBoundary_one


  SUBROUTINE CloseBoundary_all(this,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), DIMENSION(4) :: this
    INTEGER       :: dir
    !------------------------------------------------------------------------!
    INTENT(IN)    :: dir
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    SELECT CASE(dir)
    CASE(WEST,EAST,SOUTH,NORTH)
       CALL CloseBoundary_one(this(dir))
    END SELECT
  END SUBROUTINE CloseBoundary_all


END MODULE boundary_generic
