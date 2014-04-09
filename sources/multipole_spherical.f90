!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: multipole_spherical.f90                                           #
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
! module for spherical multipole expansion
!----------------------------------------------------------------------------!
MODULE multipole_spherical
  USE physics_common, ONLY : Physics_TYP
  USE mesh_common, ONLY : Selection_TYP
  USE multipole_common, InitMultipole_common => InitMultipole
  USE geometry_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE InitMultipole
     MODULE PROCEDURE InitMultipole_basic, InitMultipole_common
  END INTERFACE
  CHARACTER(LEN=32), PARAMETER :: EXP_NAME = "spherical"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Multipole_TYP, &
       ! methods
       InitMultipole, &
       CloseMultipole, &
       InitMultipole_spherical, &
       CalculatePotential_spherical, &
       CloseMultipole_spherical, &
       GetType, &
       GetName, &
       GetRank, &
       GetNumProcs, &
       Initialized, &
       Info, &
       Warning, &
       Error       
  !--------------------------------------------------------------------------!

Contains

  SUBROUTINE InitMultipole_basic(this,etype,ename,gtype,cart_coords,volume, &
       imin,imax,jmin,jmax,ireg,oreg,order)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP) :: this
    TYPE(Selection_TYP) :: ireg,oreg(:)
    INTEGER             :: gtype,etype,imin,imax,jmin,jmax,order
    CHARACTER(LEN=*)    :: ename
    REAL, DIMENSION(imin:imax,jmin:jmax,2) :: cart_coords
    REAL, DIMENSION(imin:imax,jmin:jmax)   :: volume
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP)  :: geometry
    INTEGER             :: err
    !------------------------------------------------------------------------!
    INTENT(IN)          :: cart_coords,etype,ename,gtype,imin,imax,jmin,jmax, &
                           ireg,oreg,order
    INTENT(INOUT)       :: this
    !------------------------------------------------------------------------!
    CALL InitMultipole_common(this,etype,ename,imin,imax,jmin,jmax,ireg,oreg,order)
    ! allocate memory for curvilinear coordinate field
    ALLOCATE(this%coords(this%IMIN:this%IMAX,this%JMIN:this%JMAX,2), &
         this%volume(this%IMIN:this%IMAX,this%JMIN:this%JMAX), &
         STAT=err)
    IF (err.NE.0) CALL Error(this,"InitMultipole_basic","unable to allocate memory")
    ! generate temporary geometry object to compute
    ! some geometrical quantities
    CALL InitGeometry(geometry,gtype)
    ! obtain curvilinear coordinates for the given geometry
    CALL Convert2Curvilinear(geometry,cart_coords,this%coords)
    ! assign cell volumes
    this%volume(:,:) = volume(:,:)
  END SUBROUTINE InitMultipole_basic


  SUBROUTINE InitMultipole_spherical(this,etype,cart_coords,volume,imin,imax,&
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
    INTEGER             :: i,j,l
    INTEGER             :: err
    !------------------------------------------------------------------------!
    INTENT(IN)          :: etype,cart_coords,volume,imin,imax,jmin,jmax,ireg,&
                           oreg,order
    INTENT(INOUT)       :: this
    !------------------------------------------------------------------------!
    CALL InitMultipole(this,etype,EXP_NAME,SPHERICAL,cart_coords,volume, &
         imin,imax,jmin,jmax,ireg,oreg,order)
    ! check whether radius is independent of polar angle, i.e. second index
    this%SPHERICAL=.TRUE.
    DO j=this%JMIN+1,this%JMAX
       IF (MAXVAL(ABS(1.0-this%coords(:,j-1,1)/this%coords(:,j,1))) &
            .GT.5*EPSILON(1.0)) THEN
          this%SPHERICAL = .FALSE.
          EXIT
       END IF
    END DO
!this%SPHERICAL = .FALSE.
    ! allocate memory
    IF (this%SPHERICAL) THEN
                                              ! only one entry in j at this%JMIN
       ALLOCATE(this%wmass(this%IMIN:this%IMAX,this%JMIN:this%JMIN,0:this%ORDER), &
            STAT=err)
    ELSE
       ALLOCATE(this%wmass(this%IMIN:this%IMAX,this%JMIN:this%JMAX,0:this%ORDER), &
            STAT=err)
    END IF
    IF (err.EQ.0) &
         ALLOCATE(this%radius(this%IMIN:this%IMAX,this%JMIN:this%JMAX), &
                  this%invr(this%IMIN:this%IMAX,this%JMIN:this%JMAX), &
                  this%Pl(this%IMIN:this%IMAX,this%JMIN:this%JMAX,0:this%ORDER), &
                  this%PldV(this%IMIN:this%IMAX,this%JMIN:this%JMAX,0:this%ORDER), &
                  this%temp(this%IMIN:this%IMAX,this%JMIN:this%JMAX,0:this%ORDER), &
                  this%iregion%mask(this%IMIN:this%IMAX,this%JMIN:this%JMAX), &
                  STAT=err)
    IF (err.NE.0) CALL Error(this,"InitMultipole_spherical","unable to allocate memory")
    ! zero weighted cell masses
    this%wmass(:,:,:) = 0.0
    ! for spherical grid the first curvilinear coordinate
    ! is the radial distance to the origin
    this%radius(:,:) = ABS(this%coords(:,:,1))     ! should allways be > 0
    this%invr(:,:) = 1./(this%radius(:,:)+TINY(1.0))
    ! set mask array for input region
    this%iregion%mask(:,:) = .FALSE.
    this%iregion%mask(this%iregion%IMIN:this%iregion%IMAX,this%iregion%JMIN:this%iregion%JMAX) &
         = .TRUE.
    ! compute the Legendre Polynomials up to order this%ORDER on the whole grid
!CDIR UNROLL=2
    DO l=0,1  ! zero and first order polynomials
!CDIR COLLAPSE
       DO j=this%JMIN,this%JMAX
          DO i=this%IMIN,this%IMAX
!CDIR IEXPAND
             this%Pl(i,j,l) = LegendrePolynomial(l,COS(this%coords(i,j,2)))
          END DO
       END DO
       WHERE (this%iregion%mask(:,:))
          ! multiply polynomials by volume elements
          this%PldV(:,:,l) = this%Pl(:,:,l) * this%volume(:,:) &
                  * (0.25/PI)  ! FIXME: additional factor required for selfgravity module
       ELSEWHERE
          ! zero outside input region
          this%PldV(:,:,l) = 0.0
       END WHERE
    END DO
!CDIR UNROLL=8
    DO l=2,this%ORDER ! second and higher order polynomials
!CDIR COLLAPSE
       DO j=this%JMIN,this%JMAX
!CDIR NODEP
          DO i=this%IMIN,this%IMAX
!CDIR IEXPAND
             this%Pl(i,j,l) = LegendrePolynomial(l,COS(this%coords(i,j,2)), &
                  this%Pl(i,j,l-1),this%Pl(i,j,l-2))
             ! multiply polynomials by volume elements
             this%PldV(i,j,l) = this%Pl(i,j,l) * this%volume(i,j) &
                  * (0.25/PI)  ! FIXME: additional factor required for selfgravity module
          END DO
       END DO
       WHERE (this%iregion%mask(:,:))
          ! multiply polynomials by volume elements
          this%PldV(:,:,l) = this%Pl(:,:,l) * this%volume(:,:) &
                  * (0.25/PI)  ! FIXME: additional factor required for selfgravity module
       ELSEWHERE
          ! zero outside input region
          this%PldV(:,:,l) = 0.0
       END WHERE
    END DO
    DEALLOCATE(this%coords,this%volume,this%iregion%mask)
  END SUBROUTINE InitMultipole_spherical


  PURE SUBROUTINE CalculatePotential_spherical(this,Physics,Selection,rho,Phi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP)  :: this
    TYPE(Physics_TYP)    :: Physics
    TYPE(Selection_TYP)  :: Selection(:)
    REAL, DIMENSION(this%IMIN:this%IMAX,this%JMIN:this%JMAX) :: rho, Phi
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k,i0,j0,l,jmin,jmax
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Physics,Selection,rho
    INTENT(INOUT) :: this,Phi
    !------------------------------------------------------------------------!
    IF (this%SPHERICAL) THEN
       ! if radius is independent of j, i.e. spherical mesh
       ! skip the loop over j
       jmin = this%JMIN
       jmax = this%JMIN
       ! compute mass contained in each cell and multiply by Pl*dV
       ! in case of a spherical grid, one can sum up ring masses
!CDIR UNROLL=8
       DO l=0,this%ORDER
!CDIR NODEP
          DO i=this%iregion%IMIN,this%iregion%IMAX
             ! sum up all masses with the same radius (i.e. ring section)
             ! remark: this%PldV = 0 outside input region
             this%wmass(i,jmin,l) = SUM(rho(i,:)*this%PldV(i,:,l))
          END DO
       END DO
    ELSE
       ! normal input selection indices if the grid is not spherical
       jmin = this%iregion%JMIN
       jmax = this%iregion%JMAX
       ! compute mass contained in each cell and multiply by Pl*dV       
!CDIR UNROLL=8
       DO l=0,this%ORDER
!CDIR COLLAPSE
          DO j=jmin,jmax
             DO i=this%IMIN,this%IMAX
                ! remark: this%PldV = 0 outside input region
                this%wmass(i,j,l) = rho(i,j) * this%PldV(i,j,l)
             END DO
          END DO
       END DO
    END IF

    ! compute the potential in all output selections
!CDIR UNROLL=4
    DO k=1,SIZE(Selection)
       DO j0=Selection(k)%JMIN,Selection(k)%JMAX
          DO i0=Selection(k)%IMIN,Selection(k)%IMAX
             Phi(i0,j0) = 0.0
!CDIR COLLAPSE
             ! compute zero order moment
             DO j=jmin,jmax
                DO i=this%IMIN,this%IMAX
                   this%temp(i,j,0) = MIN(this%invr(i0,j0),this%invr(i,j))
                   ! remark: this%wmass = 0 outside input region
                   Phi(i0,j0) = Phi(i0,j0) - this%Pl(i0,j0,0) &
                        * this%temp(i,j,0)*this%wmass(i,j,0)
                END DO
             END DO
             ! add higher orders
!CDIR UNROLL=8
             DO l=1,this%ORDER
!CDIR COLLAPSE
                DO j=jmin,jmax
                   DO i=this%IMIN,this%IMAX
                      ! remark: this%wmass = 0 outside input region
                      this%temp(i,j,l) = this%temp(i,j,l-1) * this%temp(i,j,0) &
                           * MIN(this%radius(i0,j0),this%radius(i,j))
                      Phi(i0,j0) = Phi(i0,j0) - this%Pl(i0,j0,l) &
                           * this%temp(i,j,l)*this%wmass(i,j,l)
                   END DO
                END DO
             END DO
           END DO
       END DO
    END DO
  END SUBROUTINE CalculatePotential_spherical


  SUBROUTINE CloseMultipole_spherical(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)  :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%radius,this%invr,this%Pl,this%PldV,this%wmass)
    CALL CloseMultipole(this)
  END SUBROUTINE CloseMultipole_spherical

END MODULE multipole_spherical
