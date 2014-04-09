!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: roots.f90                                                         #
!#                                                                           #
!# Copyright (C) 2006-2008,2011                                              #
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
! root finding subroutines:
! 1. Newton's method combined with bisection, for ill-posed problems
! 2. Regula falsi (default)
!----------------------------------------------------------------------------!
MODULE Roots
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  LOGICAL, PARAMETER :: DEBUG_ROOTS = .FALSE.
  INTEGER, PARAMETER :: MAX_ITERATIONS = 1000
  REAL, PARAMETER ::    EPS = 4*EPSILON(EPS)
  !--------------------------------------------------------------------------!
  INTERFACE GetRoot
     MODULE PROCEDURE  GetRoot_regfalsi
  END INTERFACE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       GetRoot, &
       GetRoot_newton, &
       GetRoot_regfalsi
  !--------------------------------------------------------------------------!

CONTAINS

  FUNCTION GetRoot_newton(funcd,x1,x2) RESULT(root)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x1,x2
    REAL :: root
    !------------------------------------------------------------------------!
    INTERFACE
       SUBROUTINE funcd(x,fx,dfx)
         IMPLICIT NONE
         REAL, INTENT(IN)  :: x
         REAL, INTENT(OUT) :: fx,dfx
       END SUBROUTINE funcd
    END INTERFACE
    !------------------------------------------------------------------------!
    REAL    :: fm,dfm,fl,dfl,fr,dfr
    REAL    :: xm,xl,xr,dx
    INTEGER :: i
    !------------------------------------------------------------------------!
    ! compute left and right function values
    xl = MIN(x1,x2)
    xr = MAX(x1,x2)
    CALL funcd(xl,fl,dfl)
    CALL funcd(xr,fr,dfr)
    ! check if root is within the interval [x1,x2]
    IF (fl*fr.GT.0.0) THEN
       WRITE (0,'(A,2(ES13.6,A4),2(ES13.6,A3))') &
            "ERROR in GetRoot_newton: no root " // ACHAR(10) // &
            "      f(",xl,")*f(",xr,") = ", fl, " * ",fr, " >0"
       STOP
    END IF
   
    xm = 0.5*(xl+xr)
    ! main loop
    DO i=1,MAX_ITERATIONS
       ! Newton iteration step
       CALL funcd(xm,fm,dfm)
       dx = fm / (dfm+TINY(dx))  ! avoid division by 0
       root = xm - dx
       IF (ABS(fm).LE.EPS) EXIT
       ! check if we are out of bounds
       IF ((root.LT.xl).OR.(root.GT.xr)) THEN
          ! bisection
          IF ( fl*fm.GT.0.0 ) THEN
             xl = xm
             fl = fm
          ELSE
             xr = xm
             fr = fm
          END IF
          IF (ABS(xr-xl).LT.EPS) THEN
             root = xm
             EXIT
          END IF
          xm = 0.5*(xl+xr)
       ELSE
          xm = root
       END IF
    END DO
    IF (DEBUG_ROOTS) &
         PRINT '(A,I10,4(ES14.6))',"NEWTON:  ",i,xl,xr,xm,fm
    ! check convergence
    IF (i.GE.MAX_ITERATIONS) &
       PRINT *,"WARNING in GetRoot_newton: no convergence, final dx ", dx

  END FUNCTION GetRoot_newton


  FUNCTION GetRoot_regfalsi(func,x1,x2) RESULT(root)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x1,x2
    REAL :: root
    !------------------------------------------------------------------------!
    INTERFACE
       SUBROUTINE func(x,fx)
         IMPLICIT NONE
         REAL, INTENT(IN)  :: x
         REAL, INTENT(OUT) :: fx
       END SUBROUTINE func
    END INTERFACE
    !------------------------------------------------------------------------!
    REAL    :: fm,fl,fr
    REAL    :: xm,xl,xr,dx
    INTEGER :: i
    !------------------------------------------------------------------------!
    ! compute left and right function values
    xl = MIN(x1,x2)
    xr = MAX(x1,x2)
    CALL func(xl,fl)
    CALL func(xr,fr)
    ! check if root is within the interval [x1,x2]
    IF (fl*fr.GT.0.0) THEN
       WRITE (0,'(A,2(ES13.6,A4),2(ES13.6,A3))') &
            "ERROR in GetRoot_regfalsi: no root " // ACHAR(10) // &
            "      f(",xl,")*f(",xr,") = ", fl, " * ",fr, " >0"
       STOP
    END IF
    ! main loop
    DO i=1,MAX_ITERATIONS
       ! regula falsi
       dx = (xl-xr)*fl/(fl - fr + TINY(fl))  ! avoid division by 0
       xm = xl - dx
       root = xm
       CALL func(xm,fm)
       ! check abort criteron
       IF (ABS(fm).LE.EPS) EXIT
       IF (fm*fl.GT.0.0) THEN
          xl=xm
          fl=fm
       ELSE
          xr=xm
          fr=fm
       END IF
    END DO
    IF (DEBUG_ROOTS) &
         PRINT '(A,I10,4(ES14.6))',"REGFAL:  ",i,xl,xr,xm,fm
    ! check convergence
    IF (i.GE.MAX_ITERATIONS) &
       PRINT *,"WARNING in GetRoot_regfalsi: no convergence, final dx ", dx
  END FUNCTION GetRoot_regfalsi

END MODULE Roots
