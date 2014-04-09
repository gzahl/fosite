!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: multipole_cylindrical.f90                                           #
!#                                                                           #
!# Copyright (C) 2009-2011                                                   #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
!> module for cylindrical multipole expansion
!----------------------------------------------------------------------------!
MODULE multipole_cylindrical
  USE physics_common, ONLY : Physics_TYP
  USE mesh_common, ONLY : Selection_TYP
  USE multipole_spherical
  USE geometry_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: EXP_NAME = "cylindrical"
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE CalculatePotential_cylindrical
     MODULE PROCEDURE CalculatePotential_1_cyl, &
                      CalculatePotential_2_cyl
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! methods
       InitMultipole_cylindrical, &
       CalculatePotential_cylindrical, &
       CalculatePotential_1_cyl, &
       CalculatePotential_2_cyl, &
       CloseMultipole_cylindrical
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitMultipole_cylindrical(this,etype,cart_coords,volume,imin,imax,&
       jmin,jmax,ireg,oreg,order)
    USE functions
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP) :: this
    TYPE(Selection_TYP) :: ireg,oreg(:)
    INTEGER             :: etype,imin,imax,jmin,jmax,order
    REAL, DIMENSION(imin:imax,jmin:jmax,2) :: cart_coords
    REAL, DIMENSION(imin:imax,jmin:jmax)   :: volume
    !------------------------------------------------------------------------!
    INTEGER             :: i,i0,j,j0,k,err
    !------------------------------------------------------------------------!
    INTENT(IN)          :: etype,cart_coords,volume,imin,imax,jmin,jmax,ireg,&
                           oreg,order
    INTENT(INOUT)       :: this
    !------------------------------------------------------------------------!
    CALL InitMultipole(this,etype,EXP_NAME,CYLINDRICAL,cart_coords,volume, &
         imin,imax,jmin,jmax,ireg,oreg,order)
    ! allocate memory
    ALLOCATE(this%invr(this%IMIN:this%IMAX,this%JMIN:this%JMAX), &
         this%invsqrtr(this%IMIN:this%IMAX,this%JMIN:this%JMAX), &
         this%iregion%mask(this%IMIN:this%IMAX,this%JMIN:this%JMAX), &
         this%wmass(this%IMIN:this%IMAX,this%JMIN:this%JMAX,1), &
         this%temp(this%IMIN:this%IMAX,this%JMIN:this%JMAX,1), &
         this%gfactors(SIZE(this%oregion)), &
         STAT=err)
    DO k=1,SIZE(this%oregion)
       IF (err.NE.0) EXIT
       this%gfactors(k)%range(1)%IMIN = this%IMIN
       this%gfactors(k)%range(1)%IMAX = this%IMAX
       this%gfactors(k)%range(1)%JMIN = this%JMIN
       this%gfactors(k)%range(1)%JMAX = this%JMAX
       this%gfactors(k)%range(2) = this%oregion(k)
       ALLOCATE(this%gfactors(k)%data(&
                this%gfactors(k)%range(1)%IMIN:this%gfactors(k)%range(1)%IMAX,&
                this%gfactors(k)%range(1)%JMIN:this%gfactors(k)%range(1)%JMAX,&
                this%gfactors(k)%range(2)%IMIN:this%gfactors(k)%range(2)%IMAX,&
                this%gfactors(k)%range(2)%JMIN:this%gfactors(k)%range(2)%JMAX),&
                STAT=err)
    END DO
    IF (err.NE.0) CALL Error(this,"InitMultipole_cylindrical","unable to allocate memory")
    ! for cylindrical grid the second curvilinear coordinate
    ! is the radial distance to the origin
    this%coords(:,:,2) = ABS( this%coords(:,:,2)) ! = |r|
    this%invr(:,:) = 1./(this%coords(:,:,2) + TINY(1.0)) ! = 1./|r|
    this%invsqrtr(:,:) = SQRT(this%invr(:,:)) ! = 1./sqrt(|r|)
    ! set mask array for input region
    this%iregion%mask(:,:) = .FALSE.
    this%iregion%mask(this%iregion%IMIN:this%iregion%IMAX,&
                      this%iregion%JMIN:this%iregion%JMAX) = .TRUE.
    ! compute geometrical factors
    WHERE (this%iregion%mask(:,:))
       ! using wmass for temporary storage
       this%wmass(:,:,1) = -this%volume(:,:) * this%invsqrtr(:,:) / PI 
    ELSEWHERE
       ! set wmass to 0 outside iregion
       this%wmass(:,:,1) = 0.0
    END WHERE
    DO k=1,SIZE(this%gfactors)
       DO j0=this%gfactors(k)%range(2)%JMIN,this%gfactors(k)%range(2)%JMAX
          DO i0=this%gfactors(k)%range(2)%IMIN,this%gfactors(k)%range(2)%IMAX
!CDIR COLLAPSE
             DO j=this%JMIN,this%JMAX
                DO i=this%IMIN,this%IMAX ! loop over all i to allow collapse
                   ! gfactors are 0 outside iregion, because of wmass (see above)
!CDIR IEXPAND
                   this%gfactors(k)%data(i,j,i0,j0) = &
                          LegendreFunction_QminHalf( 1.0 + &
                        ( (this%coords(i0,j0,2)-this%coords(i,j,2))**2 &
                        + (this%coords(i0,j0,1)-this%coords(i,j,1))**2 ) &
                        * 0.5*this%invr(i0,j0)*this%invr(i,j)) &
                        * this%wmass(i,j,1) * this%invsqrtr(i0,j0)
                END DO
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE InitMultipole_cylindrical


  PURE SUBROUTINE CalculatePotential_1_cyl(this,Physics,Selection,rho,Phi)
    USE functions
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP)  :: this
    TYPE(Physics_TYP)    :: Physics
    TYPE(Selection_TYP)  :: Selection(:)
    REAL, DIMENSION(this%IMIN:this%IMAX,this%JMIN:this%JMAX) :: rho, Phi
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k,i0,j0
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Physics,Selection,rho
    INTENT(INOUT) :: this,Phi
    !------------------------------------------------------------------------!
    ! compute weighted mass
    this%wmass(:,:,1) = this%volume(:,:) * rho(:,:) * this%invsqrtr(:,:)
!CDIR UNROLL=4
    DO k=1,SIZE(Selection)
       DO j0=Selection(k)%JMIN,Selection(k)%JMAX
          DO i0=Selection(k)%IMIN,Selection(k)%IMAX
             Phi(i0,j0) = 0.0
!CDIR COLLAPSE
             DO j=this%iregion%JMIN,this%iregion%JMAX
                DO i=this%IMIN,this%IMAX ! loop over all i to allow collapse
                   Phi(i0,j0) = Phi(i0,j0) &
                        - this%invsqrtr(i0,j0) / PI * this%wmass(i,j,1) &
!CDIR IEXPAND
                        * LegendreFunction_QminHalf( 1.0 + &
                        ( (this%coords(i0,j0,2)-this%coords(i,j,2))**2 &
                        + (this%coords(i0,j0,1)-this%coords(i,j,1))**2) &
                        * 0.5*this%invr(i0,j0)*this%invr(i,j))
                END DO
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE CalculatePotential_1_cyl


  PURE SUBROUTINE CalculatePotential_2_cyl(this,Physics,rho,Phi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP)  :: this
    TYPE(Physics_TYP)    :: Physics
    REAL, DIMENSION(this%IMIN:this%IMAX,this%JMIN:this%JMAX) :: rho, Phi
    !------------------------------------------------------------------------!
    INTEGER :: k,i0,j0,i,j
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Physics,rho
    INTENT(INOUT) :: this,Phi
    !------------------------------------------------------------------------!
    ! compute the potential in all output selections
!CDIR UNROLL=4
    DO k=1,SIZE(this%gfactors)
       DO j0=this%gfactors(k)%range(2)%JMIN,this%gfactors(k)%range(2)%JMAX
          DO i0=this%gfactors(k)%range(2)%IMIN,this%gfactors(k)%range(2)%IMAX
             Phi(i0,j0) = 0.0
!CDIR COLLAPSE
             DO j=this%gfactors(k)%range(1)%JMIN,this%gfactors(k)%range(1)%JMAX
                DO i=this%IMIN,this%IMAX ! loop over all i, i.e. allow collapse
                   ! remark: gfactors are 0 outside iregion
                   Phi(i0,j0) = Phi(i0,j0) + rho(i,j)*this%gfactors(k)%data(i,j,i0,j0)
                END DO
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE CalculatePotential_2_cyl


  SUBROUTINE CloseMultipole_cylindrical(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP) :: this
    !------------------------------------------------------------------------!
    INTEGER :: k
    !------------------------------------------------------------------------!
    INTENT(INOUT)  :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%coords,this%volume,this%invr,this%invsqrtr,this%iregion%mask,&
         this%wmass,this%temp)
    DO k=1,SIZE(this%oregion)
       DEALLOCATE(this%gfactors(k)%data)
    END DO
    DEALLOCATE(this%gfactors)
    CALL CloseMultipole(this)
  END SUBROUTINE CloseMultipole_cylindrical


  ! compute the half degree Legendre function of the second kind Q_{-1/2}
  ! as a function of cylindrical coordinates (z,r) with respect
  ! to some base point (z0,r0)
  ELEMENTAL FUNCTION QminHalf(z0,r0,z,r) RESULT(q)
    USE functions
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: z0,r0,z,r
    REAL :: q
    !------------------------------------------------------------------------!
    REAL :: k
    !------------------------------------------------------------------------!
    k = SQRT( (4.0*r0*r) / ((r0+r)**2 + (z0-z)**2 + TINY(r)) )
!CDIR IEXPAND
    q = k*EllipticIntegral_K(k)
  END FUNCTION QminHalf

END MODULE multipole_cylindrical
