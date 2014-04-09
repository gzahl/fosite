!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_logcylindrical.f90                                       #
!#                                                                           #
!# Copyright (C) 2011,2014                                                   # 
!# Björn Sperling <sperling@astrophysik.uni-kiel.de>                         #
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
!> \author Björn Sperling
!! \author Tobias Illenseer
!!
!! \brief define properties of a 2.5D ln(cosh) cylindrical mesh
!!
!! \warning This geometry lacks rigorous testing.
!!
!! The geometry distinguishes 3 regions
!!   1. region: constant grid spacing (= \f$dx_\mathrm{min}\f$)
!!   2. region: grid spacing increases linear
!!   3. region: constant grid spacing (= \f$dx_\mathrm{max}\f$)
!!
!! The three scaling parameters are:
!! \f{eqnarray*}{
!!    q &=& \frac{dx_\mathrm{max}}{dx_\mathrm{min}} \\
!!  r_0 &=& \text{center of region 2} \\
!!   dl &=& \text{size of region 2}
!! \f}
!! The parameters are transformed to
!! \f[
!!    b = \frac{dl}{4} \qquad a = b\frac{q-1}{q+1} \qquad v = -\frac{r_0}{b}
!! \f]
!! for internal use.
!!
!! Coordinate transformation:
!! \f{eqnarray*}{
!!    x &=& r(\eta) \sin(\varphi) \\
!!    y &=& r(\eta) \cos(\varphi) \\
!!    z &=& z \\
!!    \text{with}\quad r(\eta) &=& a\ln(\cosh(\xi))+b(\xi-v)
!! \f}
!!
!! Scale factors:
!! \f{eqnarray*}{
!!    h_z   &=& 1 \\
!!    h_\xi &=& a \tanh(\xi)+b \\
!!    h_\varphi &=& r(\xi)
!! \f}
!!
!! To calculate the mesh size use asymptotic expansion of the
!! coordinate transformation:
!! \f{eqnarray*}{
!!    r(\xi) &=& (b+a)\xi - (a\ln(2)+b v) \quad\mathrm{for}\quad \xi \gg 0 \\
!!    r(\xi) &=& (b-a)\xi - (a\ln(2)+b v) \quad\mathrm{for}\quad \xi \ll 0 \\
!! \f}
!! Hence:
!! \f{eqnarray*}{
!!    \xi_\mathrm{min} &=& \frac{r_\mathrm{min} + a\ln(2) + b v}{b-a} \\
!!    \xi_\mathrm{max} &=& \frac{r_\mathrm{max} + a\ln(2) + b v}{b+a}
!! \f}
!!
!! \extends geometry_cartesian
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_lncoshcylindrical
  USE geometry_cartesian
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE Convert2Cartesian_lncoshcyl
     MODULE PROCEDURE lncoshcyl2Cartesian_coords, lncoshcyl2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear_lncoshcyl
     MODULE PROCEDURE Cartesian2lncoshcyl_coords, Cartesian2lncoshcyl_vectors
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "lncoshcylindrical"
  INTEGER, PARAMETER :: MAX_ITERATIONS = 1000
  REAL, PARAMETER ::    EPS = 4*EPSILON(EPS)
  !--------------------------------------------------------------------------!
  PUBLIC :: &
      InitGeometry_lncoshcyl, &
      ScaleFactors_lncoshcyl, &
      Radius_lncoshcyl, &
      PositionVector_lncoshcyl, &
      Convert2Cartesian_lncoshcyl, &
      Convert2Curvilinear_lncoshcyl
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_lncoshcyl(this,gt,gp,gp2,gp3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    REAL, INTENT(IN)    :: gp,gp2,gp3 !q,r0,dl
    !------------------------------------------------------------------------!
    REAL                :: a,b,v
    !------------------------------------------------------------------------!
    CALL InitGeometry(this,gt,geometry_name)
    ! b = dl/4.0
    b = gp3/4.0
    ! a = b*(q-1)/(q+1)
    a = b*(gp-1.0)/(gp+1.0)
    ! v = -r0/b 
    v = -gp2/b
    CALL SetScale(this,a,b,v)
    CALL Warning(this,"InitGeometry_lncoshcyl","This geometry has not been tested rigorous.")
  END SUBROUTINE InitGeometry_lncoshcyl
    

  ELEMENTAL SUBROUTINE ScaleFactors_lncoshcyl(a,b,v,xi,hz,hxi,hphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: a,b,v,xi
    REAL, INTENT(OUT) :: hz,hxi,hphi
    !------------------------------------------------------------------------!
    hz   = 1.
    hxi  = a*TANH(xi)+b
    hphi = Radius_lncoshcyl(a,b,v,xi)
  END SUBROUTINE ScaleFactors_lncoshcyl

  ELEMENTAL FUNCTION Radius_lncoshcyl(a,b,v,xi) RESULT(radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: a,b,v,xi
    REAL              :: radius
    !------------------------------------------------------------------------!
    radius = a*LOG(COSH(xi))+b*(xi-v)
  END FUNCTION Radius_lncoshcyl


  ELEMENTAL SUBROUTINE PositionVector_lncoshcyl(a,b,v,xi,r,rx,ry)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: a,b,v,xi,r
    REAL, INTENT(OUT) :: rx,ry
    !------------------------------------------------------------------------!
    rx = Radius_lncoshcyl(a,b,v,xi)
    ry = r
  END SUBROUTINE PositionVector_lncoshcyl


  ! coordinate transformations
  ELEMENTAL SUBROUTINE lncoshcyl2Cartesian_coords(a,b,v,z,xi,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: a,b,v,z,xi
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    x = Radius_lncoshcyl(a,b,v,xi)
    y = z
  END SUBROUTINE lncoshcyl2Cartesian_coords

  ELEMENTAL SUBROUTINE Cartesian2lncoshcyl_coords(a,b,v,x,y,z,xi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: a,b,v,x,y
    REAL, INTENT(OUT) :: z,xi
    !------------------------------------------------------------------------!
    z = y
    xi = invlncosh(a,b,v,x)
  END SUBROUTINE Cartesian2lncoshcyl_coords

  ! vector transformations
  ELEMENTAL SUBROUTINE lncoshcyl2Cartesian_vectors(vz,vr,vx,vy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: vz,vr
    REAL, INTENT(OUT) :: vx,vy
    !------------------------------------------------------------------------!
    vx = vr
    vy = vz
  END SUBROUTINE lncoshcyl2Cartesian_vectors

  ! cartesian -> lncoshcylindrical
  ELEMENTAL SUBROUTINE Cartesian2lncoshcyl_vectors(vx,vy,vz,vr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: vx,vy
    REAL, INTENT(OUT) :: vz,vr
    !------------------------------------------------------------------------!
    vz = vy
    vr = vx
  END SUBROUTINE Cartesian2lncoshcyl_vectors

  !> \private Determine \f$\xi(r)\f$
  !!
  !! Solves the transcendental equation \f$r = a\ln(\cosh(\xi))+b(\xi-v)\f$ for \f$\xi\f$
  !! using regula falsi. We don't use the \link roots::GetRoot_regfalsi \endlink function
  !! from the numtools package because it is not PURE.
  PURE FUNCTION invlncosh(a,b,v,r) RESULT(xi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: a,b,v,r
    REAL :: xi
    !------------------------------------------------------------------------!
    REAL    :: fm,fl,fr
    REAL    :: xm,xl,xr,dx
    REAL    :: ximin,ximax
    INTEGER :: i
    !------------------------------------------------------------------------!
    !> We first estimate the limits for \f$\xi\f$ using asymptotic expansion of
    !! \f$r(\xi)\f$ (see formulas in module description).
    !! Hence the root must be in the interval 
    !! \f[
    !!    \max((b-a)x, (b+a)x)+0.1 < r(\xi) < \max((b-a)x-1.1 (a \ln(2)+bv),(b+a)x-1.1(a\ln(2)+bv))
    !! \f]
    IF (r.GE.0.0) THEN
       xl = 0.0
       xr = (r+1.1*(a*LOG(2.)+b*v)) / (b+a)
    ELSE
       xl = MIN(r+1.1*(a*log(2.)+b*v),r-0.1) / (b-a)
       xr = 0.0
    END IF
    ! compute left and right function values
    fl = r - Radius_lncoshcyl(a,b,v,xl)
    fr = r - Radius_lncoshcyl(a,b,v,xr)
    ! iteration loop to compute the root
    DO i=1,MAX_ITERATIONS
       ! regula falsi
       dx = (xl-xr)*fl/(fl - fr + TINY(fl))  ! avoid division by 0
       xm = xl - dx
       xi = xm
       fm = r - Radius_lncoshcyl(a,b,v,xm)
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
  END FUNCTION invlncosh 

END MODULE geometry_lncoshcylindrical

