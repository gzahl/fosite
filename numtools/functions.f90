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
       LegendreFunction_QminHalf, &
       Bessel_I0, &
       Bessel_I1, &
       Bessel_K0, &
       Bessel_K0e, &
       Bessel_K1
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
    tmp = 0.0
    IF (ABS(k).LT.1.0) THEN
       i   = 1
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
    IF(l.LT.0) THEN
        Pl = SQRT(1.0*l) ! return NaN
    ELSE
!CDIR UNROLL=20
        DO i=0,l
!CDIR IEXPAND
           Pl = LegendrePolynomial_one(i,x,Plminus1,Plminus2)
           Plminus2 = Plminus1
           Plminus1 = Pl
        END DO
    END IF
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



! Compute the modified Bessel function of the first kind as polynomial
! expansion, which has been proposed by Numerical Recipes in fortran, Second
! Edition on page 229ff. Nontheless the implementation is different and only the
! idea is used.
  ELEMENTAL FUNCTION Bessel_I0(x) RESULT(I0)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL                :: x
    REAL                :: I0
    !------------------------------------------------------------------------!
    INTENT(IN)          :: x
    !------------------------------------------------------------------------!
    REAL, DIMENSION(7), PARAMETER &
                        :: p = (/ 1.0d0, 3.5156229d0, 3.0899424d0, &
                                  1.2067492d0, 0.2659732d0, 0.360768d-1, &
                                  0.45813d-2 /)
    REAL, DIMENSION(9), PARAMETER &
                        :: q = (/ 0.39894228d0, 0.1328592d-1, 0.225319d-2, &
                                  -0.157565d-2, 0.916281d-2, -0.2057706d-1, &
                                  0.2635537d-1, -0.1647633d-1, 0.392377d-2 /)
    REAL                :: t, absx
    !------------------------------------------------------------------------!

    IF(ABS(x).LT.3.75) THEN
        t = (x/3.75)**2
        I0 = p(1)+t*(p(2)+t*(p(3)+t*(p(4)+t*(p(5)+t*(p(6)+t*p(7))))))
    ELSE
        absx = ABS(x)
        t = 3.75/absx
        I0 = (EXP(absx)/sqrt(absx)) &
             * (q(1)+t*(q(2)+t*(q(3)+t*(q(4) &
                + t*(q(5)+t*(q(6)+t*(q(7)+t*(q(8)+t*q(9)))))))))
    ENDIF

    END FUNCTION Bessel_I0



! Compute the modified Bessel function of the first kind as polynomial
! expansion, which has been proposed by Numerical Recipes in fortran, Second
! Edition on page 229ff. Nontheless the implementation is different and only the
! idea is used.
  ELEMENTAL FUNCTION Bessel_I1(x) RESULT(I1)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL                :: x
    REAL                :: I1
    !------------------------------------------------------------------------!
    INTENT(IN)          :: x
    !------------------------------------------------------------------------!
    REAL, DIMENSION(7), PARAMETER &
                        :: p = (/ 0.5d0, 0.87890594d0, 0.51498869d0, &
                                  0.15084934d0, 0.2658733d-1, 0.301532d-2, &
                                  0.32411d-3 /)
    REAL, DIMENSION(9), PARAMETER &
                        :: q = (/ 0.39894228d0, -0.3988024d-1, -0.362018d-2, &
                                  0.163801d-2, -0.1031555d-1, 0.2282967d-1, &
                                  -0.2895312d-1, 0.1787654d-1, -0.420059d-2 /)
    REAL                :: t, absx
    !------------------------------------------------------------------------!

    IF(ABS(x).LT.3.75) THEN
        t = (x/3.75)**2
        I1 = x*(p(1)+t*(p(2)+t*(p(3)+t*(p(4)+t*(p(5)+t*(p(6)+t*p(7)))))))
    ELSE
        absx = ABS(x)
        t = 3.75/absx
        I1 = (EXP(absx)/sqrt(absx)) &
             * (q(1)+t*(q(2)+t*(q(3)+t*(q(4) &
                + t*(q(5)+t*(q(6)+t*(q(7)+t*(q(8)+t*q(9)))))))))
        IF(x.LT.0.) THEN
            I1 = -I1
        END IF

    ENDIF

    END FUNCTION Bessel_I1



! Compute the modified Bessel function of the second kind as polynomial
! expansion, which has been proposed by Numerical Recipes in fortran, Second
! Edition on page 229ff. Nontheless the implementation is different and only
! the idea is used.
  ELEMENTAL FUNCTION Bessel_K0(x) RESULT(K0)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL                :: x
    REAL                :: K0
    !------------------------------------------------------------------------!
    INTENT(IN)          :: x
    !------------------------------------------------------------------------!
    REAL, DIMENSION(7), PARAMETER &
                        :: p = (/ -0.57721566d0, 0.42278420d0, 0.23069756d0, &
                                  0.3488590d-1, 0.262698d-2, 0.10750d-3, &
                                  0.74d-5 /)
    REAL, DIMENSION(7), PARAMETER &
                        :: q = (/ 1.25331414d0, -0.7832358d-1, 0.2189568d-1, &
                                  -0.1062446d-1, 0.587872d-2, -0.251540d-2, &
                                  0.53208d-3 /)
    REAL                :: t
    !------------------------------------------------------------------------!

    IF(x.LE.2.0) THEN
            t = x*x / 4.0
            K0 = (-LOG(x/2.0)*Bessel_I0(x)) &
                 + (p(1)+t*(p(2)+t*(p(3)+t*(p(4)+t*(p(5)+t*(p(6)+t*p(7)))))))
    ELSE
            t = (2.0/x)
            K0 = (EXP(-x)/SQRT(x)) &
                 * (q(1)+t*(q(2)+t*(q(3)+t*(q(4)+t*(q(5)+t*(q(6)+t*q(7)))))))
    ENDIF

    END FUNCTION Bessel_K0

! Compute the exponential scaled modified Bessel function of the second kind 
! e.g. K0e = EXP(x) * K0(x) as polynomial expansion, using coefficients from
! Abramowitz p.379 (http://people.math.sfu.ca/~cbm/aands/page_379.htm) for
! x >= 2, and exp(x)*K0(x) directly for 0<x<2.
  ELEMENTAL FUNCTION Bessel_K0e(x) RESULT(K0e)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL                :: x
    REAL                :: K0e
    !------------------------------------------------------------------------!
    INTENT(IN)          :: x
    !------------------------------------------------------------------------!
    REAL, DIMENSION(7), PARAMETER &
                        :: q = (/ 1.25331414, -0.07832358, 0.02189568,&
                                  -0.01062446, 0.00587872, -0.00251540,&
                                  0.00053208 /)
    REAL                :: t
    !------------------------------------------------------------------------!

    IF(x.LT.2.0) THEN
            K0e = EXP(x) * Bessel_K0(x)
    ELSE
            t = (2.0/x)
            K0e = (1.0/SQRT(x)) &
                 * (q(1)+t*(q(2)+t*(q(3)+t*(q(4)+t*(q(5)+t*(q(6)+t*q(7)))))))
    ENDIF

    END FUNCTION Bessel_K0e

! Compute the modified Bessel function of the second kind as polynomial
! expansion, which has been proposed by Numerical Recipes in fortran, Second
! Edition on page 229ff. Nontheless the implementation is different and only
! the idea is used.
  ELEMENTAL FUNCTION Bessel_K1(x) RESULT(K1)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL                :: x
    REAL                :: K1
    !------------------------------------------------------------------------!
    INTENT(IN)          :: x
    !------------------------------------------------------------------------!
    REAL, DIMENSION(7), PARAMETER &
                        :: p = (/ 1.0d0, 0.15443144d0, -0.67278579d0, &
                                  -0.18156897d0, -0.1919402d-1, -0.110404d-2, &
                                  -0.4686d-4 /)
    REAL, DIMENSION(7), PARAMETER &
                        :: q = (/ 1.25331414d0, 0.23498619d0, -0.3655620d-1, &
                                  0.1504268d-1, -0.780353d-2, 0.325614d-2, &
                                  -0.68245d-3 /)
    REAL                :: t
    !------------------------------------------------------------------------!

    IF(x.LE.2.0) THEN
            t = x*x / 4.0
            K1 = (LOG(x/2.0)*Bessel_I1(x)) &
                 + (1.0/x) &
                    * (p(1)+t*(p(2)+t*(p(3)+t*(p(4)+t*(p(5)+t*(p(6)+t*p(7)))))))
    ELSE
            t = (2.0/x)
            K1 = (EXP(-x)/SQRT(x)) &
                 * (q(1)+t*(q(2)+t*(q(3)+t*(q(4)+t*(q(5)+t*(q(6)+t*q(7)))))))
    ENDIF

    END FUNCTION Bessel_K1

END MODULE functions

