!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fosite_common.f 90                                                #
!#                                                                           #
!# Copyright (C) 2011                                                        #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de<                               #
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
! basic fosite module
!----------------------------------------------------------------------------!
MODULE fosite_common
  USE common_types, &
       GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Initialized_common => Initialized, Info_common => Info, &
       Warning_common => Warning, Error_common => Error
  USE mesh_common, ONLY : Mesh_TYP
  USE fluxes_common, ONLY : Fluxes_TYP
  USE physics_common, ONLY : Physics_TYP
  USE fileio_common, ONLY : FileIO_TYP
  USE timedisc_common, ONLY : Timedisc_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE GetType
     MODULE PROCEDURE GetSimType, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetSimName, GetName_common
  END INTERFACE
  INTERFACE GetRank
     MODULE PROCEDURE GetSimRank, GetRank_common
  END INTERFACE
  INTERFACE GetNumProcs
     MODULE PROCEDURE GetSimNumProcs, GetNumProcs_common
  END INTERFACE
  INTERFACE Initialized
     MODULE PROCEDURE SimInitialized, Initialized_common
  END INTERFACE
  INTERFACE Info
     MODULE PROCEDURE SimInfo, Info_common
  END INTERFACE
  INTERFACE Warning
     MODULE PROCEDURE SimWarning, Warning_common
  END INTERFACE
  INTERFACE Error
     MODULE PROCEDURE SimError, Error_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  TYPE Fosite_TYP
     TYPE(Common_TYP)       :: sim                 ! simulation              !
     TYPE(Mesh_TYP)         :: Mesh
     TYPE(Fluxes_TYP)       :: Fluxes
     TYPE(Physics_TYP)      :: Physics
     TYPE(FileIO_TYP)       :: Datafile
     TYPE(Timedisc_TYP)     :: Timedisc
     TYPE(FileIO_TYP)       :: Logfile
     INTEGER                :: iter
     LOGICAL                :: break
     DOUBLE PRECISION       :: wall_time           ! wall clock elapsed time !
     DOUBLE PRECISION       :: log_time            ! time for next log output!
     DOUBLE PRECISION       :: start_time          ! system clock start time !
     DOUBLE PRECISION       :: end_time            ! system clock end time   !
     DOUBLE PRECISION       :: run_time            ! = end_time - start_time !
#ifdef PARALLEL
     INTEGER                :: ierror
     REAL                   :: dt_all              ! min timestep of all     !
                                                   ! processes               ! 
#endif

  END TYPE Fosite_TYP
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Fosite_TYP, &
       ! methods
       InitSim, &
       CloseSim, &
       GetType, &
       GetName, &
       GetRank, &
       GetNumProcs, &
       !GetErrorMap, &
       Initialized, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSim(this,simtype,simname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)  :: this
    INTEGER           :: simtype
    CHARACTER(LEN=32) :: simname
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: simtype,simname
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%sim,simtype,simname)
  END SUBROUTINE InitSim


  PURE FUNCTION GetSimType(this) RESULT(ap)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP), INTENT(IN) :: this
    INTEGER :: ap
    !------------------------------------------------------------------------!
    ap = GetType_common(this%sim)
  END FUNCTION GetSimType


  PURE FUNCTION GetSimName(this) RESULT(an)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: an
    !------------------------------------------------------------------------!
    an = GetName_common(this%sim)
  END FUNCTION GetSimName


  PURE FUNCTION GetSimRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP), INTENT(IN) :: this
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%sim)
  END FUNCTION GetSimRank

  PURE FUNCTION GetSimNumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP), INTENT(IN) :: this
    INTEGER :: p
    !------------------------------------------------------------------------!
    p = GetNumProcs_common(this%sim)
  END FUNCTION GetSimNumProcs


!  PURE FUNCTION GetErrorMap(this, error) RESULT(c)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Fosite_TYP), INTENT(IN) :: this
!    INTEGER, INTENT(IN):: error
!    CHARACTER(LEN=1) :: c
!    !------------------------------------------------------------------------!
!    c = this%errormap(error)
!  END FUNCTION GetErrorMap


  PURE FUNCTION SimInitialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP), INTENT(IN) :: this
    LOGICAL :: i
    !------------------------------------------------------------------------!
    i = Initialized_common(this%sim)
  END FUNCTION SimInitialized

 
  SUBROUTINE SimInfo(this,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP), INTENT(IN)  :: this
    CHARACTER(LEN=*),  INTENT(IN) :: msg
    !------------------------------------------------------------------------!
    CALL Info_common(this%sim,msg)
  END SUBROUTINE SimInfo


  SUBROUTINE SimWarning(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Warning_common(this%sim,modproc,msg)
  END SUBROUTINE SimWarning


  SUBROUTINE SimError(this,modproc,msg,rank)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP), INTENT(IN)  :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    INTEGER, OPTIONAL, INTENT(IN) :: rank
    !------------------------------------------------------------------------!
    IF (PRESENT(rank)) THEN
       CALL Error_common(this%sim,modproc,msg,rank)
    ELSE
       CALL Error_common(this%sim,modproc,msg)
    END IF
  END SUBROUTINE SimError


  SUBROUTINE CloseSim(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP) :: this
    !------------------------------------------------------------------------!
    CALL CloseCommon(this%sim)
  END SUBROUTINE CloseSim


END MODULE fosite_common
