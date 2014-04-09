!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: functions.f90                                                     #
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
! some special math functions:
! - complete elliptic integrals of first and second kind
!----------------------------------------------------------------------------!
MODULE functions
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  REAL, PARAMETER    :: PI = 3.14159265358979323846
  REAL, PARAMETER    :: SQRT_TWO = 1.41421356237309504880
  REAL, PARAMETER    :: EPS_AGM = EPSILON(EPS_AGM)   ! precision of the AGM
  INTEGER, PARAMETER :: MAX_AGM_ITERATIONS = 20      ! limit iteration steps
  ! store the first 20 coefficients used in computing the Legendre polynomials
  INTEGER :: iii
  INTEGER, PARAMETER :: MAX_PL_COEFF = 20
  REAL, DIMENSION(MAX_PL_COEFF), PARAMETER &
       :: PL_COEFF = (/ (iii/(iii+1.0), iii=1,MAX_PL_COEFF) /)
  !--------------------------------------------------------------------------!
  INTERFACE LegendrePolynomial
     MODULE PROCEDURE LegendrePolynomial_one, LegendrePolynomial_all
  END INTERFACE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       Heaviside, &
       Barrier, &
       EllipticIntegrals, &
       EllipticIntegral_K, &
       EllipticIntegral_E, &
       LegendrePolynomials, &
       LegendrePolynomial, &
       LegendreFunction_QminHalf
  !--------------------------------------------------------------------------!

CONTAINS
  
  ! step function
  !      returns:   0   for x<a
  !                 1/2 for x=a
  !                 1   for x>a
  ELEMENTAL FUNCTION Heaviside(x,a) RESULT(fx)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x,a
    REAL :: fx
    fx = 0.5*(SIGN(1.0,x-a)+1.0)
  END FUNCTION Heaviside

  ! barrier function:
  !         returns    0   for x<a and x>b
  !                    1/2 for x=a and x=b
  !                    1   for a<x<b
  ELEMENTAL FUNCTION Barrier(x,a,b) RESULT(fx)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x,a,b
    REAL :: fx
    fx = Heaviside(x,a)-Heaviside(x,b)
  END FUNCTION Barrier


  ! compute the complete elliptic integrals of first and second kind
  ! using the AGM method (arithmetic-geometric-mean);
  ! 
  ! Function definitions and algorithm are given in
  ! [1] M. Abramowitz and I. A. Stegun: Handbook of Mathematical Functions,
  !     Applied Mathematics Series. National Bureau of Standards, Vol. 55, 1964
  !     online resource: http://people.math.sfu.ca/~cbm/aands/
  ! 
  ! check the value of K(k):
  ! (a) K((SQRT(6)-SQRT(2))/4) = 2**(-7/3) * 3**(1/4) * Gamma(1/3)**3 / PI
  ! (b) K(1/SQRT(2)) = SQRT(PI)/4 * Gamma(1/4)**2
  ! (c) K((SQRT(6)+SQRT(2))/4) = 2**(-7/3) * 3**(3/4) * Gamma(1/3)**3 / PI
  !
  ! check value of E(k):
  ! (a) E((SQRT(6)-SQRT(2))/4) = SQRT(PI)*3**(-1/4) * (
  !        2**(1/3) / SQRT(3) * (SQRT(PI)/Gamma(1/3))**3 
  !      + 0.125*(SQRT(3)+1) / (2**(1/3)) * (Gamma(1/3)/SQRT(PI))**3)
  ! (b) E(1/SQRT(2)) = SQRT(PI)*(PI/(Gamma(1/4)**2)+0.125*(Gamma(1/4)**2)/PI)
  ! (a) E((SQRR(6)+SQRT(2))/4) = SQRT(PI)*3**(1/4) * (
  !        2**(1/3) / SQRT(3) * (SQRT(PI)/Gamma(1/3))**3 
  !      + 0.125*(SQRT(3)-1) / (2**(1/3)) * (Gamma(1/3)/SQRT(PI))**3)
  !
  ! with Gamma(1/3) = 2.6789385347 ... and Gamma(1/4) = 3.6256099082 ...
  !
  ! returns NaN on error, i.e. for |k| > 1
  ELEMENTAL SUBROUTINE EllipticIntegrals(k,Kell,Eell)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: k ! elliptic modulus
    REAL, INTENT(OUT) :: Eell,Kell
    !------------------------------------------------------------------------!
    INTEGER :: i,n
    REAL    :: an,bn,cn,tmp
    !------------------------------------------------------------------------!
    IF (ABS(k).LT.1.0) THEN
       i   = 1
       tmp = 0.0
       an  = 1.0
       bn  = SQRT(1.0-k*k)
       cn  = k
!CDIR UNROLL=20
       DO n=1,MAX_AGM_ITERATIONS ! do not loop forever
          tmp = tmp + cn*cn*i
          IF (ABS(cn).GT.0.5*EPS_AGM*ABS(an)) THEN
             cn = 0.5*(an-bn)
             bn = SQRT(an*bn)
             an = an-cn     ! = 0.5*(an+bn_old)
             i = ISHFT(i,1) ! = 2*i
#if !defined(NECSX8) && !defined(NECSX9)
          ELSE
             EXIT  ! exit prohibits vectorization
#endif
          END IF
       END DO
       ! Kell = 0.5*PI / an better: Kell = 0.5*PI/0.5*(an+bn)
       Kell = PI / (an + bn + TINY(Kell)) ! avoid division by 0
    ELSE
       Kell = SQRT(-1.0*ABS(k)) ! return NaN
    END IF
    Eell = (1.0-0.5*tmp)*Kell
  END SUBROUTINE EllipticIntegrals
  

  ! compute the complete elliptic integral of the first kind
  ! using the AGM method (arithmetic-geometric-mean)
  !
  ! returns NaN on error, i.e. for |k| > 1
  ! compute the complete elliptic integral of the first kind
  ! using the AGM method (arithmetic-geometric-mean)
  !
  ! returns NaN on error, i.e. for |k| > 1
  ELEMENTAL FUNCTION EllipticIntegral_K(k) RESULT(Kell)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: k ! elliptic modulus
    REAL :: Kell
    !------------------------------------------------------------------------!
    REAL :: an,bn,tmp
    INTEGER :: i
    !------------------------------------------------------------------------!
    IF (ABS(k).LT.1.0) THEN
       an = 1.0
       bn = SQRT(1.0-k*k) ! = NaN if |k| > 1
!CDIR UNROLL=20
       DO i=1,MAX_AGM_ITERATIONS ! do not loop forever
          IF (ABS(an-bn).GT.EPS_AGM*ABS(an)) THEN
             tmp = 0.5*(an + bn)
             bn  = SQRT(an*bn)
             an  = tmp
#if !defined(NECSX8) && !defined(NECSX9)
          ELSE
             EXIT  ! exit prohibits vectorization
#endif
          END IF
       END DO
       ! Kell = 0.5*PI / an -> next iteration would compute an = 0.5*(an+bn)
       ! hence return an even better approximation:
       Kell = PI / ( an + bn + TINY(Kell) ) ! avoid division by 0
    ELSE
       Kell = SQRT(-1.0*ABS(k)) ! return NaN       
    END IF
  END FUNCTION EllipticIntegral_K


  ! compute the complete elliptic integral of the second kind
  ! using the AGM method (arithmetic-geometric-mean)
  ! 
  ! returns NaN on error, i.e. for |k| > 1
  ELEMENTAL FUNCTION EllipticIntegral_E(k) RESULT(Eell)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: k ! elliptic modulus 
    REAL :: Eell
    !------------------------------------------------------------------------!
    REAL :: dummy
    !------------------------------------------------------------------------!
    CALL EllipticIntegrals(k,dummy,Eell)
  END FUNCTION EllipticIntegral_E


  PURE SUBROUTINE LegendrePolynomials(l,x,P)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER              :: l
    REAL                 :: x
    REAL, DIMENSION(0:l) :: P
    !------------------------------------------------------------------------!
    INTEGER              :: i
    !------------------------------------------------------------------------!
    INTENT(IN)           :: l,x
    INTENT(OUT)          :: P
    !------------------------------------------------------------------------!
    IF(l.LT.0) THEN
       P(:) = SQRT(1.0*l) ! return NaN
    ELSE
       P(0) = 1.0
       IF (l.GE.1) THEN
          P(1) = x
          ! use the constant coefficients for the first 
          ! MAX_PL_COEFF Legendre Polynomials
          DO i=2,MIN(l,MAX_PL_COEFF)
             P(i) = (1.0+PL_COEFF(i))*x*P(i-1) - PL_COEFF(i)*P(i-2)
          END DO
          ! from MAX_PL_COEFF+1 to l compute the coefficients i/(i+1);
          ! this is probably slower, because of the division
          DO i=MAX_PL_COEFF+1,l
             ! recurrence formula for the Legendre polynomials:
             ! P(i) = ( (2*i+1)*x*P(i-1) - i*P(i-2) ) / (i+1)
             !      = (1+i/(i+1))*x*P(i-1) - i/(i+1)*P(i-2)
             P(i) = i/(i+1.0) ! temporary (is a number close to 1)
             P(i) = (1.0+P(i))*x*P(i-1) - P(i)*P(i-2)
          END DO
       END IF
    END IF
  END SUBROUTINE LegendrePolynomials


  ELEMENTAL FUNCTION LegendrePolynomial_one(l,x,Plminus1,Plminus2) RESULT (Pl)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER    :: l
    REAL       :: x,Plminus1,Plminus2,Pl
    !------------------------------------------------------------------------!
    REAL       :: c
    !------------------------------------------------------------------------!
    INTENT(IN) :: l,x,Plminus1,Plminus2
    !------------------------------------------------------------------------!
    IF(l.LT.0) THEN
       Pl = SQRT(1.0*l) ! return NaN
    ELSE
       SELECT CASE(l)
       CASE(0)
          Pl = 1.0
       CASE(1)
          Pl = x
       CASE DEFAULT
          ! speed up things a little with predefined coefficients
          IF (l.LE.MAX_PL_COEFF) THEN
             c = PL_COEFF(l) 
          ELSE
             c = l/(l+1.0) ! temporary storage
          END IF
          ! do recursion step
          Pl = (1.0+c)*x*Plminus1 - c*Plminus2
       END SELECT
    END IF
  END FUNCTION LegendrePolynomial_one


  ELEMENTAL FUNCTION LegendrePolynomial_all(l,x) RESULT (Pl)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER        :: l
    REAL           :: x,Pl
    !------------------------------------------------------------------------!
    INTEGER        :: i
    REAL           :: Plminus1,Plminus2
    !------------------------------------------------------------------------!
    INTENT(IN)     :: l,x
    !------------------------------------------------------------------------!
!CDIR UNROLL=20
    DO i=0,l
!CDIR IEXPAND
       Pl = LegendrePolynomial_one(i,x,Plminus1,Plminus2)
       Plminus2 = Plminus1
       Plminus1 = Pl
    END DO
  END FUNCTION LegendrePolynomial_all

  ! Computation of the order 0 and -1/2 degree 
  ! Legendre function of the second kind Q_{-1/2}
  ! using the AGM method (arithmetic-geometric-mean)
  ! 
  ! returns NaN on error, i.e. for x <= 1
  ELEMENTAL FUNCTION LegendreFunction_QminHalf(x) RESULT (q)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL       :: x,q
    !------------------------------------------------------------------------!
    INTEGER    :: i
    REAL       :: an,bn,tmp
    !------------------------------------------------------------------------!
    INTENT(IN) :: x
    !------------------------------------------------------------------------!
    IF (x.GE.1.0) THEN
       an = SQRT(x-1.0)
       bn = SQRT(x+1.0) ! if x < 1 bn = NaN, hence q = NaN
!CDIR UNROLL=20
       DO i=1,MAX_AGM_ITERATIONS ! do not loop forever
          IF (ABS(an-bn).GT.EPS_AGM*ABS(an)) THEN
             tmp = 0.5*(an + bn)
             bn = SQRT(an*bn)
             an = tmp
#if !defined(NECSX8) && !defined(NECSX9)
          ELSE
             EXIT  ! exit prohibits vectorization
#endif
          END IF
       END DO
!       q = 0.5*PI*SQRT_TWO / an
       q = PI*SQRT_TWO / (an+bn+TINY(q))
    ELSE
       q = SQRT(-1.0*ABS(x)) ! return NaN
    END IF
  END FUNCTION LegendreFunction_QminHalf

END MODULE functions

