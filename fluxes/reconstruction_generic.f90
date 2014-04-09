!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: reconstruction_generic.f90                                        #
!#                                                                           #
!# Copyright (C) 2007-2012                                                   #
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
  USE reconstruction_constant, InitReconstruction_common => InitReconstruction, &
       CloseReconstruction_common => CloseReconstruction
  USE reconstruction_linear
  USE mesh_common, ONLY : Mesh_TYP, Initialized
  USE physics_common, ONLY : Physics_TYP, Initialized
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
       CloseReconstruction, &
       CalculateSlopes, &
       CalculateStates, &
       PrimRecon, &
       GetType, &
       GetName, &
       GetRank, &
       GetNumProcs, &
       Initialized, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitReconstruction(this,Mesh,Physics,order,variables,limiter,theta)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP) :: this
    TYPE(Mesh_TYP)           :: Mesh
    TYPE(Physics_TYP)        :: Physics
    INTEGER, OPTIONAL        :: order
    LOGICAL, OPTIONAL        :: variables
    INTEGER, OPTIONAL        :: limiter
    REAL, OPTIONAL           :: theta
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32)        :: infostr
    INTEGER                  :: order_default,limiter_default
    LOGICAL                  :: variables_default
    REAL                     :: theta_default
    !------------------------------------------------------------------------!
    INTENT(IN)               :: Mesh,Physics,order,limiter,variables,theta
    INTENT(INOUT)            :: this
    !------------------------------------------------------------------------!
    ! check initialization of Mesh and Physics
    IF (.NOT.Initialized(Mesh).OR..NOT.Initialized(Physics)) &
         CALL Error(this,"InitFluxes","mesh and/or physics module uninitialized")

    ! set general reconstruction defaults
    IF (PRESENT(order)) THEN
       order_default = order
    ELSE
       order_default = 2
    END IF
    IF (PRESENT(variables)) THEN
       variables_default = variables
    ELSE
       variables_default = CONSERVATIVE
    END IF
   
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

    SELECT CASE(order_default)
    CASE(CONSTANT)
       CALL InitReconstruction_constant(this,order_default,variables_default)
    CASE(LINEAR)
       CALL InitReconstruction_linear(this,Mesh,Physics,order_default, &
                                      variables_default,limiter_default,theta_default)
    CASE DEFAULT
       CALL Error(this,"InitReconstruction",  "unsupported reconstruction order")
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
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(CONSTANT)
       ! do nothing
    CASE(LINEAR)
!CDIR IEXPAND
       CALL CalculateSlopes_linear(this,Mesh,Physics,rvar)
    END SELECT

  END SUBROUTINE CalculateSlopes


  PURE SUBROUTINE CalculateStates(this,Mesh,Physics,npos,dx,dy,rvar,rstates)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP) :: this
    TYPE(Mesh_TYP)           :: Mesh
    TYPE(Physics_TYP)        :: Physics
    INTEGER                  :: npos
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,npos) &
                             :: dx,dy
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                             :: rvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,npos,Physics%VNUM) &
                             :: rstates
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,npos,dx,dy,rvar
    INTENT(INOUT) :: this
    INTENT(OUT)   :: rstates
    !------------------------------------------------------------------------!
!CDIR IEXPAND    
    SELECT CASE(GetType(this))
    CASE(CONSTANT)
!CDIR IEXPAND
       CALL CalculateStates_constant(this,Mesh,Physics,npos,rvar,rstates)
    CASE(LINEAR)
!CDIR IEXPAND
       CALL CalculateStates_linear(this,Mesh,Physics,npos,dx,dy,rvar,rstates)
    END SELECT
  END SUBROUTINE CalculateStates


  SUBROUTINE CloseReconstruction(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)            :: this
    !------------------------------------------------------------------------!
    IF (.NOT.Initialized(this)) &
         CALL Error(this,"CloseReconstruction","not initialized")
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(CONSTANT)
       ! do nothing
    CASE(LINEAR)
       CALL CloseReconstruction_linear(this)
    CASE DEFAULT
       CALL Error(this, "CloseReconstruction",  "Unknown reconstruction type.")
    END SELECT
    CALL CloseReconstruction_common(this)
  END SUBROUTINE CloseReconstruction

END MODULE Reconstruction_generic
