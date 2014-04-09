!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: timedisc_common.f90                                               #
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
! basic module for time dicretization
!----------------------------------------------------------------------------!
MODULE timedisc_common
  USE common_types, &
       GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Initialized_common => Initialized, Info_common => Info, &
       Warning_common => Warning, Error_common => Error
  USE boundary_common, ONLY : Boundary_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE GetType
     MODULE PROCEDURE GetODESolver, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetODESolverName, GetName_common
  END INTERFACE
  INTERFACE GetRank
     MODULE PROCEDURE GetTimediscRank, GetRank_common
  END INTERFACE
  INTERFACE GetNumProcs
     MODULE PROCEDURE GetTimediscNumProcs, GetNumProcs_common
  END INTERFACE
  INTERFACE Initialized
     MODULE PROCEDURE TimediscInitialized, Initialized_common
  END INTERFACE
  INTERFACE Info
     MODULE PROCEDURE TimediscInfo, Info_common
  END INTERFACE
  INTERFACE Warning
     MODULE PROCEDURE TimediscWarning, Warning_common
  END INTERFACE
  INTERFACE Error
     MODULE PROCEDURE TimediscError, Error_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  TYPE Timedisc_TYP
     TYPE(Common_TYP) :: odesolver                     ! Runge-Kutta, etc.   !
     TYPE(Boundary_TYP), DIMENSION(4) :: Boundary  ! one for each boundary   !
     REAL             :: order                         ! time order          !
     REAL             :: cfl                           ! Courant number      !
     REAL             :: dt                            ! actual time step    !
     REAL             :: dtold                         ! last time step      !
     REAL             :: dtmin                         ! min dt of act. calc !
     REAL             :: time                          ! actual time         !
     REAL             :: stoptime                      ! end of simulation   !
     REAL             :: dtlimit                       ! lower limit for dt  !
     REAL             :: tol_rel                       ! rel. error tolerance!
     INTEGER          :: maxiter                       ! maximal iterations  !
     INTEGER          :: n_adj                         ! num. of adjustments !
     INTEGER          :: m                             ! cal steps in emb. RK!
     REAL, DIMENSION(:), POINTER      :: tol_abs       ! abs. error tolerance!
     REAL, DIMENSION(:,:,:), POINTER  :: pvar, cvar    ! prim/cons vars      !
     REAL, DIMENSION(:,:,:), POINTER  :: pold, cold    ! old prim/cons vars  !
     REAL, DIMENSION(:,:,:), POINTER  :: ptmp,ctmp     ! temporary cvars     !
     REAL, DIMENSION(:,:,:,:), POINTER:: coeff         ! coefficents         !
     REAL, DIMENSION(:), POINTER      :: A1,A2,a       ! needed by           !
     REAL, DIMENSION(:,:), POINTER    :: b             ! embedded RK         !
     REAL, DIMENSION(:,:,:), POINTER  :: src,geo_src   ! source terms        !
     REAL, DIMENSION(:,:,:), POINTER  :: rhs,bxrhs,byrhs!ODE right hand side !
     REAL, DIMENSION(:,:,:), POINTER  :: xflux, yflux  ! num. flux func.     !
     REAL, DIMENSION(:,:,:), POINTER  :: amax          ! max. wave speeds    !
  END TYPE Timedisc_TYP
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Timedisc_TYP, &
       ! methods
       InitTimedisc, &
       CloseTimedisc, &
       GetType, &
       GetName, &
       GetOrder, &
       GetCFL, &
       GetRank, &
       GetNumProcs, &
       Initialized, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitTimedisc(this,os,on,order,stoptime,cfl,dtlimit,maxiter)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    INTEGER            :: os
    CHARACTER(LEN=32)  :: on
    INTEGER            :: order,maxiter
    REAL               :: stoptime,cfl,dtlimit
    !------------------------------------------------------------------------!
    INTENT(IN)         :: os,on,order,stoptime,cfl,dtlimit,maxiter
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%odesolver,os,on)

    this%order    = order
    this%cfl      = cfl
    this%stoptime = stoptime
    this%dtlimit  = dtlimit
    this%maxiter  = maxiter

    this%dt       = this%stoptime
    this%dtold    = this%dt
    this%dtmin    = this%dt
    this%time     = 0.0
    this%n_adj    = 0
  END SUBROUTINE InitTimedisc


  SUBROUTINE CloseTimedisc(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL CloseCommon(this%odesolver)
  END SUBROUTINE CloseTimedisc


  PURE FUNCTION GetODESolver(this) RESULT(os)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP), INTENT(IN) :: this
    INTEGER :: os
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    os = GetType_common(this%odesolver)
  END FUNCTION GetODESolver


  PURE FUNCTION GetODESolverName(this) RESULT(on)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: on
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    on = GetName_common(this%odesolver)
  END FUNCTION GetODEsolverName


  PURE FUNCTION GetOrder(this) RESULT(odr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP), INTENT(IN) :: this
    INTEGER :: odr
    !------------------------------------------------------------------------!
    odr = this%order
  END FUNCTION GetOrder


  PURE FUNCTION GetCFL(this) RESULT(cfl)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP), INTENT(IN) :: this
    REAL :: cfl
    !------------------------------------------------------------------------!
    cfl = this%CFL
  END FUNCTION GetCFL

  PURE FUNCTION GetTimediscRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP), INTENT(IN) :: this
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%odesolver)
  END FUNCTION GetTimediscRank


  PURE FUNCTION GetTimediscNumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP), INTENT(IN) :: this
    INTEGER :: p
    !------------------------------------------------------------------------!
    p = GetNumProcs_common(this%odesolver)
  END FUNCTION GetTimediscNumProcs


  PURE FUNCTION TimediscInitialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP), INTENT(IN) :: this
    LOGICAL :: i
    !------------------------------------------------------------------------!
    i = Initialized_common(this%odesolver)
  END FUNCTION TimediscInitialized


  SUBROUTINE TimediscInfo(this,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: msg
    !------------------------------------------------------------------------!
    CALL Info_common(this%odesolver,msg)
  END SUBROUTINE TimediscInfo


  SUBROUTINE TimediscWarning(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Warning_common(this%odesolver,modproc,msg)
  END SUBROUTINE TimediscWarning


  SUBROUTINE TimediscError(this,modproc,msg,rank)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP), INTENT(IN):: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    INTEGER, OPTIONAL, INTENT(IN) :: rank
    !------------------------------------------------------------------------!
    IF (PRESENT(rank)) THEN
       CALL Error_common(this%odesolver,modproc,msg,rank)
    ELSE
       CALL Error_common(this%odesolver,modproc,msg)
    END IF
  END SUBROUTINE TimediscError

END MODULE timedisc_common
