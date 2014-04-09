!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fluxes_common.f90                                                 #
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
! basic fluxes module
!----------------------------------------------------------------------------!
MODULE fluxes_common
  USE common_types, GetType_common => GetType, GetName_common => GetName
  USE reconstruction_common, ONLY : Reconstruction_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE GetType
     MODULE PROCEDURE GetQuadrule, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetQuadruleName, GetName_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  ! fluxes data structure
  TYPE Fluxes_TYP
     TYPE(Common_TYP)           :: quadrule        ! midpoint, trapezoidal.. !
     TYPE(Reconstruction_TYP)   :: Reconstruction  ! store recon. settings   !
     ! various data fields
     REAL, DIMENSION(:,:,:,:), &
          POINTER               :: prim,cons       ! states                  !
     REAL, DIMENSION(:,:,:,:), &
          POINTER               :: rstates         ! prim or cons            !
     REAL, DIMENSION(:,:,:,:), &
          POINTER               :: pfluxes,qfluxes ! physical fluxes         !
     REAL, DIMENSION(:,:,:), &
          POINTER               :: temp1,temp2,temp3,temp4
  END TYPE Fluxes_TYP
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Fluxes_TYP, &
       ! methods
       InitFluxes, &
       GetType, &
       GetName
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitFluxes(this,qrule,qname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP)  :: this
    INTEGER           :: qrule
    CHARACTER(LEN=32) :: qname
    !------------------------------------------------------------------------!
    INTENT(IN)        :: qrule,qname
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%quadrule,qrule,qname)
  END SUBROUTINE InitFluxes


  PURE FUNCTION GetQuadrule(this) RESULT(qr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP), INTENT(IN) :: this
    INTEGER :: qr
    !------------------------------------------------------------------------!
    qr = GetType_common(this%quadrule)
  END FUNCTION GetQuadrule


  PURE FUNCTION GetQuadruleName(this) RESULT(qn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: qn
    !------------------------------------------------------------------------!
    qn = GetName_common(this%quadrule)
  END FUNCTION GetQuadruleName


END MODULE fluxes_common
