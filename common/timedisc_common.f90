!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: timedisc_common.f90                                               #
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
! basic module for time dicretization
!----------------------------------------------------------------------------!
MODULE timedisc_common
  USE common_types, GetType_common => GetType, GetName_common => GetName
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE GetType
     MODULE PROCEDURE GetODESolver, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetODESolverName, GetName_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  TYPE Timedisc_TYP
     TYPE(Common_TYP) :: odesolver                     ! Runge-Kutta, etc.   !
     REAL             :: order                         ! time order          !
     REAL             :: cfl                           ! Courant number      !
     REAL             :: dt                            ! actual time step    !
     REAL             :: dtmin                         ! min dt of act. calc !
     REAL             :: time                          ! actual time         !
     REAL             :: stoptime                      ! end of simulation   !
     REAL             :: dtlimit                       ! lower limit for dt  !
     INTEGER          :: maxiter                       ! maximal iterations  !
     INTEGER          :: n_adj                         ! num. of adjustments !
     REAL, DIMENSION(:,:,:), POINTER :: pvar, cvar     ! prim/cons vars      !
     REAL, DIMENSION(:,:,:), POINTER :: pold, cold     ! old prim/cons vars  !
     REAL, DIMENSION(:,:,:), POINTER :: src, geo_src   ! source terms        !
     REAL, DIMENSION(:,:,:), POINTER :: xflux, yflux   ! num. flux func.     !
     REAL, DIMENSION(:,:,:), POINTER :: amax           ! max. wave speeds    !
     REAL, POINTER    :: eta(:,:)                      ! parameter for rk    !
  END TYPE Timedisc_TYP
  SAVE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Timedisc_TYP, &
       ! methods
       InitTimedisc, &
       GetType, &
       GetName, &
       GetOrder, &
       GetCFL
  !--------------------------------------------------------------------------!

CONTAINS

  PURE SUBROUTINE InitTimedisc(this,os,on,order,cfl,stoptime,dtlimit,maxiter)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    INTEGER            :: os
    CHARACTER(LEN=32)  :: on
    INTEGER            :: order,maxiter
    REAL               :: cfl,stoptime,dtlimit
    !------------------------------------------------------------------------!
    INTENT(IN)         :: os,on,order,cfl,stoptime,dtlimit,maxiter
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%odesolver,os,on)

    this%order    = order
    this%cfl      = cfl
    this%stoptime = stoptime
    this%dtlimit  = dtlimit
    this%maxiter  = maxiter

    this%dt       = this%stoptime
    this%dtmin    = this%dt
    this%time     = 0.0
    this%n_adj    = 0
  END SUBROUTINE InitTimedisc


  PURE FUNCTION GetODESolver(this) RESULT(os)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP), INTENT(IN) :: this
    INTEGER :: os
    !------------------------------------------------------------------------!
    os = GetType_common(this%odesolver)
  END FUNCTION GetODESolver


  PURE FUNCTION GetODESolverName(this) RESULT(on)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: on
    !------------------------------------------------------------------------!
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


END MODULE timedisc_common
