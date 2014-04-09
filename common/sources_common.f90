!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_common.f90                                                #
!#                                                                           #
!# Copyright (C) 2006-2010                                                   #
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
! basic sources module
!----------------------------------------------------------------------------!
MODULE sources_common
  USE common_types, GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Info_common => Info, Warning_common => Warning, Error_common => Error
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE GetType
     MODULE PROCEDURE GetSourceType, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetSourceTypeName, GetName_common
  END INTERFACE
  INTERFACE GetRank
     MODULE PROCEDURE GetSourcesRank, GetRank_common
  END INTERFACE
  INTERFACE GetNumProcs
     MODULE PROCEDURE GetSourcesNumProcs, GetNumProcs_common
  END INTERFACE
  INTERFACE Info
     MODULE PROCEDURE SourcesInfo, Info_common
  END INTERFACE
  INTERFACE Warning
     MODULE PROCEDURE SourcesWarning, Warning_common
  END INTERFACE
  INTERFACE Error
     MODULE PROCEDURE SourcesError_rank0, SourcesError_rankX, Error_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  TYPE Grid_TYP
     REAL, DIMENSION(:,:), POINTER   :: u,rho,a,da,b,db,c,invc,d,&
                                        vol,bhx,bhy,&! coarser grids for multigrid !
                                        mlint,mlext
     REAL, DIMENSION(:,:,:), POINTER :: bccart,curv,QdivsqrtR
     INTEGER                         :: ni,nj
     INTEGER, DIMENSION(:,:), POINTER :: ij2k, k2ij
     REAL                            :: hi,hj,invhi2,invhj2
  END TYPE Grid_TYP
  TYPE Sources_TYP
     TYPE(Common_TYP)                :: sourcetype   ! type of source term   !
     TYPE(Sources_TYP), POINTER      :: next => null() ! next source in list !
     TYPE(Common_TYP)                :: potential    ! newton or wiita       !
     TYPE(Common_TYP)                :: viscosity    ! molecular,alpha,beta  !
     TYPE(Grid_TYP), POINTER         :: grid(:)      ! coarser grids for multigrid !
     REAL                            :: mass         ! mass of point source  !
     REAL                            :: mdot         ! disk accretion rate   !
     REAL                            :: dynconst,bulkconst ! viscosity const.!
     REAL                            :: cvis         ! viscous Courant no.   !
     REAL                            :: MAXRESIDNORM ! max error of residium (multigrid)!
     REAL                            :: MAXAGMNORM   ! max error of arith-geom-mean     !
     INTEGER                         :: outbound     ! outflow boundary      !
     INTEGER                         :: MAXMULT      ! max number of multipol moments   !
     INTEGER                         :: BOUNDARYTYPE ! type of boundary calc (mult.expansion) !
     INTEGER                         :: MGminlevel   ! min multigrid level   !
     INTEGER                         :: ngrid        ! number of grids       !
     INTEGER, DIMENSION(4)           :: Boundary     ! boundary for poisson  !
                                      ! problem (Dirichlet,Neumann,Periodic) !
     REAL, DIMENSION(:,:,:), POINTER :: accel,accart ! acceleration          !
     REAL, DIMENSION(:,:), POINTER   :: radius       ! distance to origin    !
     REAL, DIMENSION(:,:,:), POINTER :: gxr3         ! = GN*x/radius**3      !
     REAL, DIMENSION(:,:), POINTER   :: cellmass     ! rho*dV                !
     REAL, DIMENSION(:,:), POINTER   :: dynvis, &    ! dynamic, kinematic &  !
                                      kinvis,bulkvis !    bulk viscosity     !
     REAL, DIMENSION(:,:), POINTER   :: btxx,btyy,&  ! components of the     !
          btzz,btxy,btxz,btyz                        !    stress tensor      !
     REAL, DIMENSION(:,:,:), POINTER :: ftxx,ftyy,&
          ftzz,ftxy,ftxz,ftyz
  END TYPE Sources_TYP
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       Grid_TYP, &
       ! methods
       InitSources, &
       GetSourcesPointer, &
       GetType, &
       GetName, &
       GetRank, &
       GetNumProcs, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources(this,stype,sname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    INTEGER           :: stype
    CHARACTER(LEN=32) :: sname
    !------------------------------------------------------------------------!
    INTENT(IN)        :: stype,sname
    INTENT(OUT)       :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%sourcetype,stype,sname)
    NULLIFY(this%next)
  END SUBROUTINE InitSources


  FUNCTION GetSourcesPointer(list,stype) RESULT(sp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: list,sp
    INTEGER, INTENT(IN) :: stype
    !------------------------------------------------------------------------!
    sp => list
    DO
       IF (ASSOCIATED(sp).EQV..FALSE.) EXIT
       IF (GetType(sp).EQ.stype) RETURN
       sp => sp%next
    END DO
  END FUNCTION GetSourcesPointer


  PURE FUNCTION GetSourcesRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%sourcetype)
  END FUNCTION GetSourcesRank


  PURE FUNCTION GetSourcesNumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    INTEGER :: p
    !------------------------------------------------------------------------!
    p = GetNumProcs_common(this%sourcetype)
  END FUNCTION GetSourcesNumProcs


  PURE FUNCTION GetSourceType(this) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    INTEGER :: st
    !------------------------------------------------------------------------!
    st = GetType_common(this%sourcetype)
  END FUNCTION GetSourceType


  PURE FUNCTION GetSourceTypeName(this) RESULT(sn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: sn
    !------------------------------------------------------------------------!
    sn = GetName_common(this%sourcetype)
  END FUNCTION GetSourceTypeName


  SUBROUTINE SourcesInfo(this,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: msg
    !------------------------------------------------------------------------!
    CALL Info_common(this%sourcetype,msg)
  END SUBROUTINE SourcesInfo


  SUBROUTINE SourcesWarning(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Warning_common(this%sourcetype,modproc,msg)
  END SUBROUTINE SourcesWarning


  SUBROUTINE SourcesError_rank0(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Error_common(this%sourcetype,modproc,msg)
  END SUBROUTINE SourcesError_rank0


  SUBROUTINE SourcesError_rankX(this,modproc,msg,rank)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    INTEGER, INTENT(IN)           :: rank
    !------------------------------------------------------------------------!
    CALL Error_common(this%sourcetype,modproc,msg,rank)
  END SUBROUTINE SourcesError_rankX

END MODULE sources_common
