!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: roots.f90                                                         #
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
! root finding via Newton iteration and bisection
!----------------------------------------------------------------------------!
MODULE Roots
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: MAXIT=100000
  !--------------------------------------------------------------------------!
  PUBLIC :: RtNewtBisec, GetRoot
  !--------------------------------------------------------------------------!

CONTAINS

  FUNCTION RtNewtBisec(funcd,x1,x2,p1,p2,xacc) RESULT(rt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x1,x2,p1,p2,xacc
    REAL :: rt
    !------------------------------------------------------------------------!
    INTERFACE
       SUBROUTINE funcd(x,p1,p2,fx,dfx)
         IMPLICIT NONE
         REAL, INTENT(IN)  :: x,p1,p2
         REAL, INTENT(OUT) :: fx,dfx
       END SUBROUTINE funcd
    END INTERFACE
    !------------------------------------------------------------------------!
    REAL    :: fm,dfm,fl,dfl,fr,dfr
    REAL    :: xm,xl,xr,dx
    INTEGER :: i
    !------------------------------------------------------------------------!

    xm = 0.5*(x1+x2)
    xl = MIN(x1,x2)
    xr = MAX(x1,x2)
    CALL funcd(xl,p1,p2,fl,dfl)          
    CALL funcd(xr,p1,p2,fr,dfr)

    IF ( (fl.GT.0.0 .AND. fr.GT.0.0) .OR. &
         (fl.LT.0.0 .AND. fr.LT.0.0) ) THEN
       PRINT *, "ERROR in RtNewtBisec: root must be bracketed between x1 and x2"
       STOP
    END IF
 
    ! main loop
    DO i=1,MAXIT
       ! Newton iteration step
       CALL funcd(xm,p1,p2,fm,dfm)
       dx = fm/dfm
       rt = xm - dx
       IF (ABS(dx).LT.xacc) EXIT
       ! check if we are out of bounds
       IF ((rt.LT.xl).OR.(rt.GT.xr)) THEN
          ! bisection
          IF ( ((fm.LT.0.0).AND.(fl.LT.0.0)).OR. &
               ((fm.GT.0.0).AND.(fl.GT.0.0)) ) THEN
             xl = xm
             fl = fm
          ELSE
             xr = xm
             fr = fm
          END IF
          IF (ABS(xr-xl).LT.xacc) THEN
             rt = xm
             EXIT
          END IF
          xm = 0.5*(xl+xr)
       ELSE
          xm = rt
       END IF
    END DO

    IF (i.EQ.MAXIT) THEN
       PRINT *, "ERROR in RtNewtBisec: too many iterations, aborting"
       STOP
    END IF

  END FUNCTION RtNewtBisec


  FUNCTION GetRoot(funcd,x1,x2,xacc) RESULT(root)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x1,x2,xacc
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
       WRITE (*,*) "GetRoot Error: f(x1)*f(x2) should be < 0, aborting!"
       STOP
    END IF
    ! main loop
    DO i=1,MAXIT
       ! regular falsi
       dx = fl*(xl-xr)/(fl-fr)
       xm = xl - dx
       root = xm
       CALL funcd(xm,fm,dfm)
       ! check abort criteron
       IF (ABS(fm).LT.xacc) THEN
          EXIT
       END IF
       IF (fm*fl.GT.0.0) THEN
          xl=xm
          fl=fm
       ELSE
          xr=xm
          fr=fm
       END IF
    END DO
    IF (i.GT.MAXIT) THEN
       WRITE (*,*) "WARNING: limit of iterations exceeded!"
    END IF
  END FUNCTION GetRoot

END MODULE Roots
