!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: integration.f90                                                   #
!#                                                                           #
!# Copyright (C) 2006-2008                                                   #
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
! Gauss and Romberg integration subroutines
!----------------------------------------------------------------------------!
MODULE Integration
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: GaussN  = 15       ! number of abscissas and weights !
  INTEGER, PARAMETER :: GaussN2 = GaussN / 2
  INTEGER, PARAMETER :: MAXREC = 100              ! max. depth for recursion !
  REAL, DIMENSION(:), ALLOCATABLE :: GaussXk, GaussAk ! abscissas and weights!
  REAL, DIMENSION(:,:), ALLOCATABLE :: stack      ! stack for iterative alg. !
  !--------------------------------------------------------------------------!
  PUBLIC InitIntegration, CloseIntegration, integrate
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitIntegration
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    INTEGER      :: err
    !----------------------------------------------------------------------!
    
    ALLOCATE(GaussXk(GaussN), GaussAk(GaussN), stack(3,MAXREC), &
         STAT = err)
    
    IF (err.NE.0) THEN
       PRINT *, "ERROR in InitIntegration: Can't allocate memory!"
       STOP
    END IF
      
    SELECT CASE(GaussN)
    CASE(2)
       GaussXk = (/ -0.577350269189626, 0.577350269189626 /)
       GaussAk = (/  1.000000000000000, 1.000000000000000 /)
    CASE(4)
       GaussXk = (/ -0.861136311594053, -0.339981043584856, &
                     0.339981043584856,  0.861136311594053 /)
       GaussAk = (/  0.347854845137454,  0.652145154862546, &
                     0.652145154862546,  0.347854845137454 /)
    CASE(8)
       GaussXk = (/ -0.960289856497536, -0.796666477413627, &
                    -0.525532409916329, -0.183434642495650, &
                     0.183434642495650,  0.525532409916329, &
                     0.796666477413627,  0.960289856497536 /)
       GaussAk = (/  0.101228536290376,  0.222381034453374, &
                     0.313706645877887,  0.362683783378362, &
                     0.362683783378362,  0.313706645877887, &
                     0.222381034453374,  0.101228536290376 /)
    CASE(12)
       GaussXk = (/ -0.9815606342467193, -0.9041172563704749, &
                    -0.7699026741943047, -0.5873179542866174, &
                    -0.3678314989981802, -0.1252334085114689, &
                     0.1252334085114689,  0.3678314989981802, &
                     0.5873179542866174,  0.7699026741943047, &
                     0.9041172563704749,  0.9815606342467193 /)
       GaussAk = (/  0.047175336386512,   0.106939325995318, &
                     0.160078328543346,   0.203167426723066, &
                     0.233492536538355,   0.249147045813403, &
                     0.249147045813403,   0.233492536538355, &
                     0.203167426723066,   0.160078328543346, &
                     0.106939325995318,   0.04717533638651 /)
    CASE(15)
       GaussXk = (/ -0.9879925180204854, -0.9372733924007059, &
                    -0.8482065834104272, -0.7244177313601700, &
                    -0.5709721726085388, -0.3941513470775634, &
                    -0.2011940939974345, &
                     0.0, &
                     0.2011940939974345, &
                     0.3941513470775634,  0.5709721726085388, &
                     0.7244177313601700,  0.8482065834104272, &
                     0.9372733924007059,  0.9879925180204854 /)
       GaussAk = (/  0.030753241996117,   0.070366047488108, &
                     0.107159220467172,   0.139570677926154, &
                     0.166269205816994,   0.186161000015562, &
                     0.198431485327112, &
                     0.202578241925561, &
                     0.198431485327112, &
                     0.186161000015562,   0.166269205816994, &
                     0.139570677926154,   0.107159220467172, &
                     0.070366047488108,   0.03075324199612 /)
    CASE DEFAULT
       PRINT *, "ERROR in InitIntegration: select one of 2,4,8,12,15 for GaussN"
       STOP
    END SELECT
     
  END SUBROUTINE InitIntegration


  FUNCTION integrate(fkt,xl,xr,eps,plist,method) RESULT(integral)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    INTERFACE
       FUNCTION fkt(x,plist) RESULT(fx)
         IMPLICIT NONE
         REAL, INTENT(IN) :: x
         REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
         REAL :: fx
       END FUNCTION fkt
    END INTERFACE
    !----------------------------------------------------------------------!
    INTEGER, OPTIONAL :: method                      ! Romberg or Gauss    !
    REAL :: xl, xr                                   ! integration limits  !
    REAL :: eps                                      ! precision           !
    REAL, DIMENSION(:), OPTIONAL :: plist            ! optional parameters !
    REAL :: integral
    !----------------------------------------------------------------------!
    INTEGER          :: meth
    REAL             :: tmpS,err
    !----------------------------------------------------------------------!
    INTENT(IN)       :: plist, xl, xr
    !----------------------------------------------------------------------!
      
    IF (xl.EQ.xr) THEN
       integral = 0.
       RETURN
    END IF

    ! select the method
    ! 1: adapt. Gauss recursive
    ! 2: adapt. Gauss iterative
    ! 3: Romberg
    IF (PRESENT(method).AND.((method.GT.0).AND.method.LE.3)) THEN
       meth = method
    ELSE
       ! default is recursive Gauss
       meth = 1
    END IF
    
    ! tolerance for integration
    err  = ABS(xr-xl) * eps

    SELECT CASE(meth)
    CASE(1)
       ! recursive Gauss integration
       tmpS = 0.
       IF (PRESENT(plist)) THEN
          tmpS = qgauss1D(fkt,MIN(xl,xr),MAX(xl,xr),plist)
          integral = SIGN(1.0,xr-xl)*qadaptive1D_recursive(fkt,tmpS,&
               MIN(xl,xr),MAX(xl,xr),err,plist)
       ELSE
          tmpS = qgauss1D(fkt,MIN(xl,xr),MAX(xl,xr))
          integral = SIGN(1.0,xr-xl)*qadaptive1D_recursive(fkt,tmpS,&
               MIN(xl,xr),MAX(xl,xr),err)
       END IF

    CASE(2)
       ! iterative Gauss integration
       tmpS = 0.
       IF (PRESENT(plist)) THEN
          tmpS = qgauss1D(fkt,MIN(xl,xr),MAX(xl,xr),plist)
          integral = SIGN(1.0,xr-xl)*qadaptive1D_iterative(fkt,tmpS,&
               MIN(xl,xr),MAX(xl,xr),err,plist)
       ELSE
          tmpS = qgauss1D(fkt,MIN(xl,xr),MAX(xl,xr))
          integral = SIGN(1.0,xr-xl)*qadaptive1D_iterative(fkt,tmpS,&
               MIN(xl,xr),MAX(xl,xr),err)
       END IF

    CASE(3)
       ! Romberg integration
       IF (PRESENT(plist)) THEN
          integral = SIGN(1.0,xr-xl)*qromberg(fkt,MIN(xl,xr),MAX(xl,xr),err,plist)
       ELSE
          integral = SIGN(1.0,xr-xl)*qromberg(fkt,MIN(xl,xr),MAX(xl,xr),err)
       END IF

    CASE DEFAULT
       PRINT *, "ERROR in integrate: select one of 1,2 for method"
       STOP
    END SELECT

  END FUNCTION integrate


  !------------------------------------------------------------------------!
  !
  ! 1D Integration with adaptive Gauss
  !
  !------------------------------------------------------------------------!
  RECURSIVE FUNCTION qadaptive1D_recursive(fkt, oldS, xl, xr, tol, plist) RESULT(S)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    INTERFACE
       FUNCTION fkt(x,plist) RESULT(fx)
         IMPLICIT NONE
         REAL, INTENT(IN) :: x
         REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
         REAL :: fx
       END FUNCTION fkt
    END INTERFACE
    !----------------------------------------------------------------------!
    REAL      :: oldS, xl, xr, tol
    REAL, DIMENSION(:), OPTIONAL :: plist
    REAL      :: S
    !----------------------------------------------------------------------!
    REAL :: m, res, resL, resR, qerr
    !----------------------------------------------------------------------!
    INTENT(IN)        :: oldS, xl, xr, tol, plist
    !----------------------------------------------------------------------!
      
    m = 0.5 * (xr+xl)
    resL = qgauss1D(fkt, xl, m, plist)       ! left solution
    resR = qgauss1D(fkt, m, xr, plist)       ! right solution
    res = resL + resR

    ! error estimate
    qerr = (oldS - res) / (4**GaussN - 1.0)

    IF (ABS(qerr).GE.tol) THEN
       ! recursion
       S = qadaptive1D_recursive(fkt, resL, xl, m, 0.5*tol, plist) &
            + qadaptive1D_recursive(fkt, resR, m, xr, 0.5*tol, plist)
    ELSE
       S = res + qerr
    END IF
  END FUNCTION qadaptive1D_recursive


  FUNCTION qadaptive1D_iterative(fkt, oldS, xl, xr, tol, plist) RESULT(S)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    INTERFACE
       FUNCTION fkt(x,plist) RESULT(fx)
         IMPLICIT NONE
         REAL, INTENT(IN) :: x
         REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
         REAL :: fx
       END FUNCTION fkt
    END INTERFACE
    !----------------------------------------------------------------------!
    REAL      :: oldS, xl, xr, tol
    REAL, DIMENSION(:), OPTIONAL :: plist
    REAL      :: S
    !----------------------------------------------------------------------!
    INTEGER   :: i,sptr
    REAL      :: resL, resR, resO, res, qerr
    REAL      :: l,r,m,err
    !----------------------------------------------------------------------!
    INTENT(IN)        :: oldS, tol, xl, xr, plist
    !----------------------------------------------------------------------!

    ! copy basic variables
    l = xl
    r = xr
    resO = oldS
    err = tol

    ! initialize result and stack pointer
    S = 0.
    sptr = 1

    ! main loop
    DO i=1,MAXREC
       m = 0.5 * (r+l)                         ! devide interval
       resL = qgauss1D(fkt, l, m, plist)       ! integrate over [xl..m]
       resR = qgauss1D(fkt, m, r, plist)       ! integrate over [m..xr]
       res = resL + resR                       ! add two results
       qerr = (resO - res) / (4**GaussN - 1.0) ! error estimate
       ! check error
       IF (ABS(qerr) < err) THEN
          ! store result
          S = S + res + qerr
          IF (sptr.GT.1) THEN
             ! pop values from stack
             sptr = sptr-1
             l    = stack(1,sptr)
             r    = stack(2,sptr)
             resO = stack(3,sptr)
          ELSE
             ! finished
             EXIT
          END IF
       ELSE
          ! push right values on stack
          stack(1,sptr) = m
          stack(2,sptr) = r
          stack(3,sptr) = resR
          sptr = sptr+1
          ! continue with left
          r = m
          err = 0.5*err
       END IF
    END DO

    IF (i.GE.MAXREC) THEN
       PRINT *, "ERROR in qadaptive1D_iterative: max recursion reached"
       STOP
    END IF

  END FUNCTION qadaptive1D_iterative


  !------------------------------------------------------------------------!
  !
  ! basic 1D Gauss quadrature function
  !
  !------------------------------------------------------------------------!
  FUNCTION qgauss1D(fkt,xl,xr,plist) RESULT(Sgauss)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    INTERFACE
       FUNCTION fkt(x,plist) RESULT(fx)
         IMPLICIT NONE
         REAL, INTENT(IN) :: x
         REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
         REAL :: fx
       END FUNCTION fkt
    END INTERFACE
    !----------------------------------------------------------------------!
    REAL, INTENT(IN)  :: xl, xr
    REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist ! additional parameters
    REAL :: Sgauss
    REAL :: m, h, t
    INTEGER :: j
    !----------------------------------------------------------------------!
    
    m = 0.5 * (xr+xl)
    h = 0.5 * (xr-xl)
    Sgauss = 0.0
    
    DO j=1,GaussN2
       t  = h * GaussXk(j)
       Sgauss = Sgauss + GaussAk(j) * (fkt(m-t,plist)+fkt(m+t,plist))
    END DO

    IF (MOD(GaussN,2) /= 0) THEN
       Sgauss = Sgauss + GaussAk(GaussN2+1) * fkt(m,plist)
    END IF
    Sgauss = h * Sgauss
  END FUNCTION qgauss1D


  !------------------------------------------------------------------------!
  !
  ! 1D Integration with Romberg algorithm
  !
  !------------------------------------------------------------------------!
  FUNCTION qromberg(fkt, xl, xr, tol, plist) RESULT (Srom)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    INTERFACE
       FUNCTION fkt(x, plist) RESULT(fx)
         IMPLICIT NONE
         REAL, INTENT(IN) :: x
         REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
         REAL :: fx
       END FUNCTION fkt
    END INTERFACE
    !----------------------------------------------------------------------!
    REAL, INTENT(IN) :: xl, xr
    REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist ! additional parameters 
    REAL :: Srom
    INTEGER, PARAMETER :: jmax = 10  ! max refinement
    REAL, DIMENSION(jmax,jmax) :: res
    REAL    :: h
    REAL    :: tol
    INTEGER :: i, j, n
    !----------------------------------------------------------------------!
    
    n = 2
    ! step size
    h = (xr - xl) / n
    res(:,:) = 0.0
    ! sum of trapezoidal rule at the boundaries
    res(1,1) = 0.5 * (fkt(xl,plist) + fkt(xr,plist))
    ! and in between
    DO i=1,n-1
       res(1,1) = res(1,1) + fkt(i*h,plist)
    END DO
    res(1,1) = res(1,1) * h

    DO j=2,jmax
       ! half step size
       h = 0.5 * h
       ! double grid points
       n = 2 * n
       ! sum up contributions due to the new points
       DO i=1,n-1,2
          res(j,1) = res(j,1) + fkt(i*h,plist)
       END DO
       res(j,1) = 0.5 * res(j-1,1) +  h*res(j,1)
       ! extrapolation
       DO i=1,j-1
          res(j,i+1) = res(j,i) + (res(j,i) - res(j-1,i)) / (4**i - 1)
       END DO
       ! check convergence
       IF (ABS(res(j,j)-res(j,j-1))<tol) EXIT
    END DO
    
    Srom = res(j-1,j-1)
  END FUNCTION qromberg


  SUBROUTINE CloseIntegration
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    
    DEALLOCATE(GaussXk,GaussAk,stack)
  END SUBROUTINE CloseIntegration


END MODULE Integration
