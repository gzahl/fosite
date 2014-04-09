!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: reconstruction_constant.f90                                       #
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
! basic module for constant (zero order) reconstruction
!----------------------------------------------------------------------------!
MODULE reconstruction_constant
  USE reconstruction_common
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER  :: recontype_name = "constant"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Reconstruction_TYP, &
       ! constants
       PRIMITIVE, CONSERVATIVE, &
       ! methods
       InitReconstruction_constant, &
       CalculateStates_constant, &
       GetType, &
       GetName, &
       PrimRecon
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitReconstruction_constant(this,rtype,pc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP) :: this
    INTEGER                  :: rtype
    LOGICAL                  :: pc
    !------------------------------------------------------------------------!
    INTENT(IN)               :: rtype,pc
    INTENT(INOUT)            :: this
    !------------------------------------------------------------------------!
    CALL InitReconstruction(this,rtype,recontype_name,pc)
  END SUBROUTINE InitReconstruction_constant


  PURE SUBROUTINE CalculateStates_constant(this,Mesh,Physics,npos,rvar,rstates)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP) :: this
    TYPE(Mesh_TYP)           :: Mesh
    TYPE(Physics_TYP)        :: Physics
    INTEGER                  :: npos
    REAL :: rvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)
    REAL :: rstates(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,npos,Physics%vnum)
    !------------------------------------------------------------------------!
    INTEGER :: n
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,Physics,npos,rvar
    INTENT(OUT) :: rstates
    !------------------------------------------------------------------------!

    ! reconstruct boundary states
    FORALL (n=1:npos)
       rstates(:,:,n,:) = rvar(:,:,:)
    END FORALL
  END SUBROUTINE CalculateStates_constant
  
END MODULE reconstruction_constant
