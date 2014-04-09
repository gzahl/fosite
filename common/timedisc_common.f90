!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: timedisc_common.f90                                               #
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
!> \defgroup timedisc timedisc
!! \{
!! \brief Family of timedisc modules
!! \}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!!
!! \brief basic module for time dicretization
!!
!! \extends common_types
!! \ingroup timedisc
!----------------------------------------------------------------------------!
MODULE timedisc_common
  USE common_types, &
       GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Initialized_common => Initialized, Info_common => Info, &
       Warning_common => Warning, Error_common => Error
  USE boundary_common, ONLY : Boundary_TYP
  USE common_dict, ONLY : Dict_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
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
  !> \endcond
  !--------------------------------------------------------------------------!
  TYPE Timedisc_TYP
     !> \name Variables
     TYPE(Common_TYP) :: odesolver                     !< Runge-Kutta, etc.
     TYPE(Boundary_TYP), DIMENSION(4) :: Boundary      !< one for each boundary
     INTEGER          :: order                         !< time order
     REAL             :: cfl                           !< Courant number
     REAL             :: dt                            !< actual time step
     REAL             :: dtold                         !< last time step
     REAL             :: dtmin                         !< min dt of act. calc
     INTEGER          :: dtcause,dtmincause            !< cause of dt/dtmin
     REAL             :: dtmax                         !< max dt of act. calc
     REAL,POINTER     :: dtmean, dtstddev              !< mean and stddev of timestep
     INTEGER          :: dtaccept                      !< no of accepted ts
     REAL             :: time                          !< actual time
     REAL             :: stoptime                      !< end of simulation
     REAL             :: dtlimit                       !< lower limit for dt
     REAL             :: tol_rel                       !< rel. error tolerance
     REAL             :: maxerrold                     !< old maxerr
     LOGICAL          :: break                         !< stop fosite run?
     INTEGER          :: maxiter                       !< maximal iterations
     INTEGER          :: n_adj                         !< num. of adjustments
     INTEGER          :: m                             !< cal steps in emb. RK
     INTEGER          :: degree                        !< polynom degree index in dumka
     INTEGER          :: fargo                         !< = 1 fargo enabled
     !> old error and used in the dumka method
     REAL             :: ERR_N, H_N
     REAL, DIMENSION(:), POINTER      :: tol_abs       !< abs. error tolerance
     REAL                             :: beta          !< time step friction
     REAL, DIMENSION(:,:,:), POINTER  :: pvar, cvar    !< prim/cons vars
     REAL, DIMENSION(:,:,:), POINTER  :: pold, cold    !< old prim/cons vars
     REAL, DIMENSION(:,:,:), POINTER  :: ptmp,ctmp     !< temporary cvars
     REAL, DIMENSION(:,:,:,:), POINTER:: coeff         !< coefficents
     REAL, DIMENSION(:), POINTER      :: A1,A2,a       !<    needed by
     REAL, DIMENSION(:,:), POINTER    :: b             !<    embedded RK
     !> multistep vars
     REAL, DIMENSION(:,:,:,:), POINTER:: phi,oldphi_s,&
                                         newphi_s
     REAL, DIMENSION(:), POINTER      :: gamma
     REAL, DIMENSION(:,:), POINTER    :: c
     INTEGER                          :: pc            !< = 1 predictor-corrector
     REAL, DIMENSION(:,:,:), POINTER  :: src,geo_src   !< source terms
     REAL, DIMENSION(:,:,:), POINTER  :: rhs,rhsold    !< ODE right hand side
     !> numerical fluxes divided by dy or dx
     REAL, DIMENSION(:,:,:), POINTER  :: xfluxdy,yfluxdx
     REAL, DIMENSION(:,:,:), POINTER  :: amax          !< max. wave speeds
     REAL, DIMENSION(:,:), POINTER    :: bflux         !< boundary fluxes for output
     REAL, DIMENSION(:,:,:), POINTER  :: error         !< max. wave speeds
     LOGICAL                          :: write_error   !< enable err writing
     INTEGER, DIMENSION(:), POINTER   :: shift         !< fargo annulus shift
     REAL, DIMENSION(:,:), POINTER    :: buf           !< fargo MPI buffer
     TYPE(elem_TYP),POINTER           :: ts=>Null()    !< timesteps (multistep)
  END TYPE Timedisc_TYP
  TYPE elem_TYP
    TYPE(elem_TYP), POINTER           :: next,prev
    REAL                              :: t             !< time
  END TYPE elem_TYP
  !> \}
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Timedisc_TYP, &
       elem_TYP, &
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
       Error, &
       ! ring structure
       InitRing, Get, Get_t, Set, AddFirst, AddLast, Rotate, RemoveLast
       
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public
  SUBROUTINE InitTimedisc(this,os,on)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    INTEGER            :: os
    CHARACTER(LEN=32)  :: on
    !------------------------------------------------------------------------!
    INTENT(IN)         :: os,on
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%odesolver,os,on)

    this%maxerrold= 0.
    this%dtmean   = 0.
    this%dtstddev = 0.
    this%dtaccept = 0

    this%dt       = this%stoptime
    this%dtold    = this%dt
    this%dtmin    = this%dt
    this%time     = 0.0
    this%n_adj    = 0
  END SUBROUTINE InitTimedisc


  !> \public
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

  SUBROUTINE InitRing(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(elem_TYP), POINTER  :: this
    !------------------------------------------------------------------------!
    INTEGER                  :: err
    !------------------------------------------------------------------------!
    ALLOCATE(this, STAT = err)
    this%next => this
    this%prev => this
    this%t= 0.0
  END SUBROUTINE InitRing

  RECURSIVE FUNCTION Get(this,i) RESULT(elem)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(elem_TYP), POINTER   :: this
    INTEGER, INTENT(IN)       :: i
    TYPE(elem_TYP), POINTER   :: elem
    !------------------------------------------------------------------------!
    IF (i .GT. 0) THEN
      elem=>Get(this%next,i-1)
    ELSE
      elem=>this
    END IF
  END FUNCTION Get

  FUNCTION Get_t(this,i) RESULT(t)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(elem_TYP), POINTER   :: this
    INTEGER, INTENT(IN)       :: i
    TYPE(elem_TYP), POINTER   :: elem
    REAL                      :: t
    !------------------------------------------------------------------------!
    elem=>Get(this,i)
    t = elem%t
  END FUNCTION Get_t

  SUBROUTINE Set(this,i,t) 
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(elem_TYP), POINTER   :: this
    INTEGER, INTENT(IN)       :: i
    REAL, INTENT(IN)          :: t
    TYPE(elem_TYP), POINTER   :: elem
    !------------------------------------------------------------------------!
    elem=>Get(this,i)
    elem%t = t
  END SUBROUTINE Set

  SUBROUTINE AddLast(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(elem_TYP), POINTER :: this
    INTEGER                 :: err
    TYPE(elem_TYP), POINTER :: elem
    !------------------------------------------------------------------------!
    ALLOCATE(elem,STAT = err)
    elem%t = 0.0
    this%prev%next=>elem
    elem%next=>this    
    elem%prev=>this%prev
    this%prev=>elem
  END SUBROUTINE AddLast

  SUBROUTINE AddFirst(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(elem_TYP), POINTER :: this
    TYPE(elem_TYP), POINTER :: elem
    !------------------------------------------------------------------------!
    CALL AddLast(this)
    CALL Rotate(this)
  END SUBROUTINE AddFirst

  SUBROUTINE RemoveLast(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(elem_TYP), POINTER :: this
    INTEGER                 :: err
    TYPE(elem_TYP), POINTER :: elem
    !------------------------------------------------------------------------!
    elem=>this%prev
    this%prev=>this%prev%prev
    this%prev%next=>this
    DEALLOCATE(elem,STAT = err)
  END SUBROUTINE RemoveLast

  SUBROUTINE Rotate(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(elem_TYP), POINTER :: this
    INTEGER                 :: err
    TYPE(elem_TYP), POINTER :: elem
    !------------------------------------------------------------------------!
    this=>this%prev
  END SUBROUTINE Rotate

END MODULE timedisc_common
