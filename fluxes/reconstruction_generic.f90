!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: reconstruction_generic.f90                                        #
!#                                                                           #
!# Copyright (C) 2007 Tobias Illenseer <tillense@ita.uni-heidelberg.de>      #
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
! generic module for reconstruction process
!----------------------------------------------------------------------------!
MODULE reconstruction_generic
  USE reconstruction_constant
  USE reconstruction_linear
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: CONSTANT = 1
  INTEGER, PARAMETER :: LINEAR   = 2
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Reconstruction_TYP, &
       ! constants
       CONSTANT, LINEAR, &
       PRIMITIVE, CONSERVATIVE, &
       MINMOD, MONOCENT, SWEBY, SUPERBEE, OSPRE, &
       ! methods
       InitReconstruction, &
       MallocReconstruction, &
       PrimRecon, &
       GetType, &
       CalculateSlopes, &
       CalculateStates, &
       CloseReconstruction
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitReconstruction(this,order,variables,limiter,theta)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP) :: this
    INTEGER                  :: order
    LOGICAL                  :: variables
    INTEGER, OPTIONAL        :: limiter
    REAL, OPTIONAL           :: theta
    !------------------------------------------------------------------------!
    INTEGER                  :: limiter_default
    REAL                     :: theta_default
    !------------------------------------------------------------------------!
    INTENT(IN)               :: order,limiter,variables,theta
    INTENT(INOUT)            :: this
    !------------------------------------------------------------------------!
    ! set defaults for the limiter
    IF (PRESENT(limiter)) THEN
       limiter_default = limiter
    ELSE
       limiter_default = MINMOD
    END IF
    IF (PRESENT(theta)) THEN
       theta_default = theta
    ELSE
       theta_default = 1.0
    END IF

    SELECT CASE(order)
    CASE(CONSTANT)
       CALL InitReconstruction_constant(this,order,variables)
    CASE(LINEAR)
       SELECT CASE(limiter)
       CASE(MINMOD,MONOCENT,SWEBY,SUPERBEE,OSPRE)
          CALL InitReconstruction_linear(this,order,variables,limiter_default,theta_default)
       CASE DEFAULT
          PRINT *, "ERROR in InitReconstruction: unknown limiter"
          STOP
       END SELECT
    CASE DEFAULT
       PRINT *, "ERROR in InitReconstruction: order not supported"
       STOP
    END SELECT

    ! print some information
    PRINT "(A,A)", " RECONSTR-> order:             ", TRIM(GetName(this))
    WRITE (*,FMT='(A)',ADVANCE='NO'), "            primitive:         "
    IF (PrimRecon(this)) THEN
       PRINT "(A)", "true"
    ELSE
       PRINT "(A)", "false"
    END IF
    IF (order.EQ.LINEAR) THEN
       PRINT "(A,A)", "            limiter:           ", TRIM(GetLimiterName(this))
    END IF
  END SUBROUTINE InitReconstruction


  SUBROUTINE MallocReconstruction(this,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP) :: this
    TYPE(Mesh_TYP)           :: Mesh
    TYPE(Physics_TYP)        :: Physics
    !------------------------------------------------------------------------!
    INTEGER       :: err
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    ! allocate memory for reconstruction modules
    SELECT CASE(GetType(this))
    CASE(CONSTANT)
       ! do nothing
    CASE(LINEAR)
       CALL MallocReconstruction_linear(this,Mesh,Physics)
    CASE DEFAULT
       PRINT *, "ERROR in MallocReconstruction: unknown reconstruction type"
       STOP
    END SELECT
  END SUBROUTINE MallocReconstruction


  PURE SUBROUTINE CalculateSlopes(this,Mesh,Physics,rvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP) :: this
    TYPE(Mesh_TYP)           :: Mesh
    TYPE(Physics_TYP)        :: Physics
    REAL :: rvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,rvar
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    
    SELECT CASE(GetType(this))
    CASE(CONSTANT)
       ! do nothing
    CASE(LINEAR)
       CALL CalculateSlopes_linear(this,Mesh,Physics,rvar)
    END SELECT

  END SUBROUTINE CalculateSlopes


  PURE SUBROUTINE CalculateStates(this,Mesh,Physics,npos,pos0,pos,rvar,rstates)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP) :: this
    TYPE(Mesh_TYP)           :: Mesh
    TYPE(Physics_TYP)        :: Physics
    INTEGER                  :: npos
    REAL :: pos0(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2)
    REAL :: pos(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,npos,2)
    REAL :: rvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)
    REAL :: rstates(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,npos,Physics%vnum)
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,npos,pos0,pos,rvar
    INTENT(INOUT) :: this
    INTENT(OUT)   :: rstates
    !------------------------------------------------------------------------!
    
    SELECT CASE(GetType(this))
    CASE(CONSTANT)
       CALL CalculateStates_constant(this,Mesh,Physics,npos,rvar,rstates)
    CASE(LINEAR)
       CALL CalculateStates_linear(this,Mesh,Physics,npos,pos0,pos,rvar,rstates)
    END SELECT
  END SUBROUTINE CalculateStates


  SUBROUTINE CloseReconstruction(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(IN)               :: this
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(CONSTANT)
       ! do nothing
    CASE(LINEAR)
       CALL CloseReconstruction_linear(this)
    CASE DEFAULT
       PRINT *, "ERROR in CloseReconstruction: unknown reconstruction type"
       STOP
    END SELECT
  END SUBROUTINE CloseReconstruction

END MODULE Reconstruction_generic
