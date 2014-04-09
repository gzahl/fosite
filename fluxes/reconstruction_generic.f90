!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: reconstruction_generic.f90                                        #
!#                                                                           #
!# Copyright (C) 2007-2008                                                   #
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
! generic module for reconstruction process
!----------------------------------------------------------------------------!
MODULE reconstruction_generic
  USE reconstruction_constant, InitReconstruction_common => InitReconstruction
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
       MINMOD, MONOCENT, SWEBY, SUPERBEE, OSPRE, PP, &
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
    CHARACTER(LEN=32)        :: infostr
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
       CASE(MINMOD,MONOCENT,SWEBY,SUPERBEE,OSPRE,PP)
          CALL InitReconstruction_linear(this,order,variables,limiter_default,theta_default)
       CASE DEFAULT
          CALL Error(this, "InitReconstruction", "Unknown limiter")
       END SELECT
    CASE DEFAULT
       CALL Error(this,"MallocReconstruction",  "Unknown reconstruction type.")
    END SELECT

    ! print some information
    CALL Info(this, " RECONSTR-> order:             " // TRIM(GetName(this)))
    IF (PrimRecon(this)) THEN
       WRITE (infostr,'(A)') "primitive"
    ELSE
       WRITE (infostr,'(A)') "conservative"
    END IF
    CALL Info(this, "            variables:         " // TRIM(infostr))
    IF (order.EQ.LINEAR) THEN
       CALL Info(this, "            limiter:           " // TRIM(GetLimiterName(this)))
    END IF
  END SUBROUTINE InitReconstruction


  SUBROUTINE MallocReconstruction(this,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP) :: this
    TYPE(Mesh_TYP)           :: Mesh
    TYPE(Physics_TYP)        :: Physics
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
       CALL Error(this,"MallocReconstruction",  "Unknown reconstruction type.")
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
    INTENT(INOUT)            :: this
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(CONSTANT)
       ! do nothing
    CASE(LINEAR)
       CALL CloseReconstruction_linear(this)
    CASE DEFAULT
       CALL Error(this, "CloseReconstruction",  "Unknown reconstruction type.")
    END SELECT
  END SUBROUTINE CloseReconstruction

END MODULE Reconstruction_generic
