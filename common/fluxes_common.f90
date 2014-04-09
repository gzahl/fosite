!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fluxes_common.f90                                                 #
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
!> \defgroup fluxes fluxes
!! \{
!! \brief Family of fluxes modules
!! \}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!!
!! \brief basic fluxes module
!!
!! \extends common_types
!! \ingroup fluxes
!----------------------------------------------------------------------------!
MODULE fluxes_common
  USE common_types, &
       GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Initialized_common => Initialized, Info_common => Info, &
       Warning_common => Warning, Error_common => Error
  USE reconstruction_common, ONLY : Reconstruction_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE GetType
     MODULE PROCEDURE GetQuadrule, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetQuadruleName, GetName_common
  END INTERFACE
  INTERFACE GetRank
     MODULE PROCEDURE GetFluxesRank, GetRank_common
  END INTERFACE
  INTERFACE GetNumProcs
     MODULE PROCEDURE GetFluxesNumProcs, GetNumProcs_common
  END INTERFACE
  INTERFACE Initialized
     MODULE PROCEDURE FluxesInitialized, Initialized_common
  END INTERFACE
  INTERFACE Info
     MODULE PROCEDURE FluxesInfo, Info_common
  END INTERFACE
  INTERFACE Warning
     MODULE PROCEDURE FluxesWarning, Warning_common
  END INTERFACE
  INTERFACE Error
     MODULE PROCEDURE FluxesError, Error_common
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  !> fluxes data structure
  TYPE Fluxes_TYP
     !> \name Variables
     TYPE(Common_TYP)           :: quadrule        ! midpoint, trapezoidal.. !
     TYPE(Reconstruction_TYP)   :: Reconstruction  ! store recon. settings   !
     ! various data fields
     REAL, DIMENSION(:,:,:,:), POINTER &           ! primitive and conserva. !
                                :: prim,cons       !   data on cell faces    !
     REAL, DIMENSION(:,:,:,:), POINTER &           ! reconstructed data:     !
                                :: rstates         !   prim or cons          !
     REAL, DIMENSION(:,:,:,:), POINTER :: pfluxes  ! physical fluxes         !
     REAL, DIMENSION(:,:,:), POINTER &             ! coordinate differences  !
                                :: dx,dy           ! centers to recon. pos.  !
     REAL, DIMENSION(:,:,:), POINTER &             ! boundary fluxes         !
                                :: bxflux,byflux,bxfold,byfold
     REAL                       :: viscosity       ! numerical viscosity     !
                                                   ! parameter               ! 
  END TYPE Fluxes_TYP
  !> \}
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Fluxes_TYP, &
       ! methods
       InitFluxes, &
       CloseFluxes, &
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

  !> \public
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


  !> \public
  SUBROUTINE CloseFluxes(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL CloseCommon(this%quadrule)
  END SUBROUTINE CloseFluxes


  !> \public
  PURE FUNCTION GetQuadrule(this) RESULT(qr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP), INTENT(IN) :: this
    INTEGER :: qr
    !------------------------------------------------------------------------!
    qr = GetType_common(this%quadrule)
  END FUNCTION GetQuadrule


  !> \public
  PURE FUNCTION GetQuadruleName(this) RESULT(qn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: qn
    !------------------------------------------------------------------------!
    qn = GetName_common(this%quadrule)
  END FUNCTION GetQuadruleName


  !> \public
  PURE FUNCTION GetFluxesRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP), INTENT(IN) :: this
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%quadrule)
  END FUNCTION GetFluxesRank


  !> \public
  PURE FUNCTION GetFluxesNumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP), INTENT(IN) :: this
    INTEGER :: p
    !------------------------------------------------------------------------!
    p = GetNumProcs_common(this%quadrule)
  END FUNCTION GetFluxesNumProcs


  !> \public
  PURE FUNCTION FluxesInitialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP), INTENT(IN) :: this
    LOGICAL :: i
    !------------------------------------------------------------------------!
    i = Initialized_common(this%quadrule)
  END FUNCTION FluxesInitialized


  !> \public
  SUBROUTINE FluxesInfo(this,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: msg
    !------------------------------------------------------------------------!
    CALL Info_common(this%quadrule,msg)
  END SUBROUTINE FluxesInfo


  !> \public
  SUBROUTINE FluxesWarning(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Warning_common(this%quadrule,modproc,msg)
  END SUBROUTINE FluxesWarning


  !> \public
  SUBROUTINE FluxesError(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Error_common(this%quadrule,modproc,msg)
  END SUBROUTINE FluxesError


END MODULE fluxes_common
