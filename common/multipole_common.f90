!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: multipole_common.f90                                              #
!#                                                                           #
!# Copyright (C) 2006-2011                                                   #
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
!> basic module for multipole expansion
!----------------------------------------------------------------------------!
MODULE multipole_common
  USE common_types, &
       GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Initialized_common => Initialized, Info_common => Info, &
       Warning_common => Warning, Error_common => Error
  USE mesh_common, ONLY : Selection_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE GetType
     MODULE PROCEDURE GetMultipoleType, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetMultipoleName, GetName_common
  END INTERFACE
  INTERFACE GetRank
     MODULE PROCEDURE GetMultipoleRank, GetRank_common
  END INTERFACE
  INTERFACE GetNumProcs
     MODULE PROCEDURE GetMultipoleNumProcs, GetNumProcs_common
  END INTERFACE
  INTERFACE Initialized
     MODULE PROCEDURE MultipoleInitialized, Initialized_common
  END INTERFACE
  INTERFACE Info
     MODULE PROCEDURE MultipoleInfo, Info_common
  END INTERFACE
  INTERFACE Warning
     MODULE PROCEDURE MultipoleWarning, Warning_common
  END INTERFACE
  INTERFACE Error
     MODULE PROCEDURE MultipoleError, Error_common
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  TYPE GFactors_TYP                             ! storage type for geometrical
     TYPE(Selection_TYP) :: range(2)            !   factors used in mult. exp.
     REAL, DIMENSION(:,:,:,:), POINTER :: data
  END TYPE GFactors_TYP
  TYPE Multipole_TYP
     TYPE(Common_TYP) :: exptype                ! spherical, cylindrical etc.
     TYPE(Selection_TYP) :: iregion             ! selection field for expansion
     TYPE(Selection_TYP), DIMENSION(:), POINTER &
                      :: oregion                ! selection for output regions
     TYPE(GFactors_TYP), DIMENSION(:), POINTER &
                      :: gfactors               ! geometrical factors in expansion
     INTEGER          :: ORDER                  ! highest multipole moment
     INTEGER          :: IMIN,IMAX,JMIN,JMAX    ! density and potential array
                                                !   dimensions
     LOGICAL          :: SPHERICAL              ! = .TRUE. if radius is 
                                                !   independent of polar angle
     REAL, DIMENSION(:,:,:), POINTER &
                      :: coords,&               ! curvilinear coords
                         Pl,&                   ! Legendre Polynomomials
                         PldV,&                 ! Legendre Polynom. * volume
                         wmass,&                ! weighted cell masses
                         temp                   ! temporary storage
     REAL, DIMENSION(:,:), POINTER &
                      :: radius,&               ! radial coordinate
                         invr, &                ! inverse radius
                         invsqrtr, &            ! inverse square root of radius
                         volume                 ! cell volume field
  END TYPE Multipole_TYP
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Multipole_TYP, &
       GFactors_TYP, &
       ! methods
       InitMultipole, &
       CloseMultipole, &
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

  SUBROUTINE InitMultipole(this,etype,ename,imin,imax,jmin,jmax,ireg,oreg,order)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP) :: this
    TYPE(Selection_TYP) :: ireg,oreg(:)
    INTEGER             :: etype,imin,imax,jmin,jmax,order
    CHARACTER(LEN=32)   :: ename
    !------------------------------------------------------------------------!
    INTEGER             :: err
    !------------------------------------------------------------------------!
    INTENT(IN)          :: etype,ename,imin,imax,jmin,jmax,ireg,oreg,order
    INTENT(INOUT)       :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%exptype,etype,ename)
    ! number of multipole moments
    this%ORDER = order
    ! array dimensions of the density and potential array
    this%IMIN = imin
    this%IMAX = imax
    this%JMIN = jmin
    this%JMAX = jmax
    this%iregion = ireg
    ALLOCATE(this%oregion(SIZE(oreg)),STAT=err)
    IF (err.NE.0) CALL Error(this,"InitMultipole","memory allocation failed")
    this%oregion(:) = oreg(:)
  END SUBROUTINE InitMultipole


  SUBROUTINE CloseMultipole(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL CloseCommon(this%exptype)
  END SUBROUTINE CloseMultipole


  PURE FUNCTION GetMultipoleType(this) RESULT(etype)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP), INTENT(IN) :: this
    INTEGER :: etype
    !------------------------------------------------------------------------!
    etype = GetType_common(this%exptype)
  END FUNCTION GetMultipoleType


  PURE FUNCTION GetMultipoleName(this) RESULT(ename)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: ename
    !------------------------------------------------------------------------!
    ename = GetName_common(this%exptype)
  END FUNCTION GetMultipoleName


  PURE FUNCTION GetMultipoleRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP), INTENT(IN) :: this
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%exptype)
  END FUNCTION GetMultipoleRank


  PURE FUNCTION GetMultipoleNumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP), INTENT(IN) :: this
    INTEGER :: p
    !------------------------------------------------------------------------!
    p = GetNumProcs_common(this%exptype)
  END FUNCTION GetMultipoleNumProcs


  PURE FUNCTION MultipoleInitialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP), INTENT(IN) :: this
    LOGICAL :: i
    !------------------------------------------------------------------------!
    i = Initialized_common(this%exptype)
  END FUNCTION MultipoleInitialized

 
  SUBROUTINE MultipoleInfo(this,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: msg
    !------------------------------------------------------------------------!
    CALL Info_common(this%exptype,msg)
  END SUBROUTINE MultipoleInfo


  SUBROUTINE MultipoleWarning(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Warning_common(this%exptype,modproc,msg)
  END SUBROUTINE MultipoleWarning


  SUBROUTINE MultipoleError(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Error_common(this%exptype,modproc,msg)
  END SUBROUTINE MultipoleError

END MODULE multipole_common
