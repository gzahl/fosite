!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: poisson_common.f90                                                #
!#                                                                           #
!# Copyright (C) 2011                                                        #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
! basic Poisson module
!----------------------------------------------------------------------------!
MODULE poisson_common
  USE common_types, &
       GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Initialized_common => Initialized, Info_common => Info, &
       Warning_common => Warning, Error_common => Error
  USE multipole_common, ONLY : Multipole_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE GetType
     MODULE PROCEDURE GetPoissonType, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetPoissonTypeName, GetName_common
  END INTERFACE
  INTERFACE GetRank
     MODULE PROCEDURE GetPoissonRank, GetRank_common
  END INTERFACE
  INTERFACE GetNumProcs
     MODULE PROCEDURE GetPoissonNumProcs, GetNumProcs_common
  END INTERFACE
  INTERFACE Initialized
     MODULE PROCEDURE PoissonInitialized, Initialized_common
  END INTERFACE
  INTERFACE Info
     MODULE PROCEDURE PoissonInfo, Info_common
  END INTERFACE
  INTERFACE Warning
     MODULE PROCEDURE PoissonWarning, Warning_common
  END INTERFACE
  INTERFACE Error
     MODULE PROCEDURE PoissonError_rank0, PoissonError_rankX, Error_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  TYPE Grid_TYP                     ! data type for multigrid poisson solver !
     REAL, DIMENSION(:,:), POINTER   :: u,rho         ! DELTA u = rho        !
     REAL, DIMENSION(:,:), POINTER   :: a,da,b,db,c   !    some geometric    !
     REAL, DIMENSION(:,:), POINTER   :: invc,d        !      values          !
     REAL, DIMENSION(:,:), POINTER   :: vol,bhx,bhy   ! volume, scale factors!
     REAL, DIMENSION(:,:), POINTER   :: tmp           ! temp
     REAL, DIMENSION(:,:,:), POINTER :: bccart,curv   ! cart. and curv. coord.!
     REAL, DIMENSION(:,:,:), POINTER :: tri           ! coefficients for block-iteration !
     INTEGER                         :: ni,nj         ! resolution            !
     INTEGER, DIMENSION(:,:), POINTER :: ij2k, k2ij   ! index field for bound.!
     REAL                            :: hi,hj         ! width of cells; square!  
     REAL                            :: invhi2,invhj2 !     of inverse        !
  END TYPE Grid_TYP
  TYPE Poisson_TYP
     TYPE(Common_TYP)                :: poissontype   ! type of source term   !
     TYPE(Multipole_TYP)             :: multipole     ! multipole expansion   !
     TYPE(Grid_TYP), POINTER         :: grid(:)       ! coarse grids          !
     REAL                            :: MAXRESIDNORM  ! max error of residuum !
     INTEGER                         :: RELAXTYPE     ! type of relaxation method !
     INTEGER                         :: NGRID         ! number of grids       !
     INTEGER                         :: NPRE,NPOST    ! pre and post smoothing!
     INTEGER                         :: NMAXCYCLE     ! max of iterations     !
     INTEGER                         :: MINRES        ! min resolution        !
     REAL, DIMENSION(:,:,:), POINTER :: accel         ! acceleration          !
     REAL, DIMENSION(:,:), POINTER   :: phi           ! potential             !
     INTEGER, DIMENSION(4)           :: Boundary      ! boundary condition    !
     LOGICAL                         :: DIRICHLET     ! true if min ONE bound.!
                                                      !   is set as dirichlet !
!FIXME
INTEGER, DIMENSION(:),POINTER :: relaxcount !test!!!
  END TYPE Poisson_TYP
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Poisson_TYP, &
       Grid_TYP, &
       ! methods
       InitPoisson, &
       ClosePoisson, &
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

  SUBROUTINE InitPoisson(this,stype,sname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Poisson_TYP) :: this
    INTEGER           :: stype
    CHARACTER(LEN=32) :: sname
    !------------------------------------------------------------------------!
    INTENT(IN)        :: stype,sname
    INTENT(OUT)       :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%poissontype,stype,sname)

  END SUBROUTINE InitPoisson


  SUBROUTINE ClosePoisson(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Poisson_TYP), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL CloseCommon(this%poissontype)
  END SUBROUTINE ClosePoisson


 PURE FUNCTION GetPoissonRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Poisson_TYP), INTENT(IN) :: this
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%poissontype)
  END FUNCTION GetPoissonRank


  PURE FUNCTION GetPoissonNumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Poisson_TYP), INTENT(IN) :: this
    INTEGER :: p
    !------------------------------------------------------------------------!
    p = GetNumProcs_common(this%poissontype)
  END FUNCTION GetPoissonNumProcs


  PURE FUNCTION GetPoissonType(this) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Poisson_TYP), INTENT(IN) :: this
    INTEGER :: st
    !------------------------------------------------------------------------!
    st = GetType_common(this%poissontype)
  END FUNCTION GetPoissonType


  PURE FUNCTION GetPoissonTypeName(this) RESULT(sn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Poisson_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: sn
    !------------------------------------------------------------------------!
    sn = GetName_common(this%poissontype)
  END FUNCTION GetPoissonTypeName


  PURE FUNCTION PoissonInitialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Poisson_TYP), INTENT(IN) :: this
    LOGICAL :: i
    !------------------------------------------------------------------------!
    i = Initialized_common(this%poissontype)
  END FUNCTION PoissonInitialized


  SUBROUTINE PoissonInfo(this,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Poisson_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: msg
    !------------------------------------------------------------------------!
    CALL Info_common(this%poissontype,msg)
  END SUBROUTINE PoissonInfo


  SUBROUTINE PoissonWarning(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Poisson_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Warning_common(this%poissontype,modproc,msg)
  END SUBROUTINE PoissonWarning


  SUBROUTINE PoissonError_rank0(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Poisson_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Error_common(this%poissontype,modproc,msg)
  END SUBROUTINE PoissonError_rank0


  SUBROUTINE PoissonError_rankX(this,modproc,msg,rank)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Poisson_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    INTEGER, INTENT(IN)           :: rank
    !------------------------------------------------------------------------!
    CALL Error_common(this%poissontype,modproc,msg,rank)
  END SUBROUTINE PoissonError_rankX

END MODULE poisson_common
