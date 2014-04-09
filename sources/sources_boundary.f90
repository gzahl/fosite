!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_boundary.f90                                              #
!#                                                                           #
!# Copyright (C) 2009-2010                                                   #
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
! source terms module for boundary calculation 
!----------------------------------------------------------------------------!
MODULE sources_boundary
  USE sources_common
  USE mesh_common
  USE geometry_common, ONLY : LikeSpherical
  USE geometry_generic
  USE constants_generic
  USE boundary_generic
 IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
       INTEGER, PARAMETER :: SPHERMULTEXPAN       = 1 !very slow!!!!
       INTEGER, PARAMETER :: SPHERMULTEXPANFAST   = 2 !very fast, but only spherical coor.
       INTEGER, PARAMETER :: CYLINMULTEXPAN       = 3 !fast
       INTEGER, PARAMETER :: CYLINMULTEXPANFAST   = 4 !faster then 3, but needs more memory
       INTEGER, PARAMETER :: DIRICHLET            = 0 ! Dirichlet boundary condition 
       INTEGER, PARAMETER :: NEUMANN              = 1 ! Neumann boundary condition
 !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       SPHERMULTEXPAN, SPHERMULTEXPANFAST, &
       CYLINMULTEXPAN, CYLINMULTEXPANFAST, &
       DIRICHLET,NEUMANN,&
       ! methods
       InitSources_boundary, &
       CalcBoundary, &
       CloseSources_boundary
  !--------------------------------------------------------------------------!

Contains

  SUBROUTINE InitSources_boundary(this,Mesh,Boundary,maxmult,bndrytype)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)             :: Mesh
    TYPE(Boundary_TYP), DIMENSION(4) :: Boundary
    INTEGER                    :: maxmult,bndrytype
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP):: geometry         ! geometrical properties (mg) !
    REAL                       :: mu,chi,R1,R2,z1,z2
    INTEGER                    :: ngrid,i,j,k,ni,nj,err
    !------------------------------------------------------------------------!
    INTENT(IN)                 :: Mesh,Boundary,maxmult,bndrytype
    !------------------------------------------------------------------------!
    SELECT CASE(bndrytype)
    CASE(SPHERMULTEXPAN)
       ! initialize geometry
       CALL InitGeometry(geometry,SPHERICAL)
       CALL Info(this, " SOURCES--> boundary:          spherical multipole expansion")
    CASE(SPHERMULTEXPANFAST)
       ! initialize geometry
       CALL InitGeometry(geometry,SPHERICAL)
       IF (.NOT.(LikeSpherical(Mesh%geometry)))&
          CALL Error(this, "InitSources_boundary", "Wrong geometry in case of spher. mult. expansion ")
       CALL Info(this, " SOURCES--> boundary:          spherical multipole expansion (fast)")
    CASE(CYLINMULTEXPAN)
       ! initialize geometry
       CALL InitGeometry(geometry,CYLINDRICAL)
       CALL Info(this, " SOURCES--> boundary:          cylindrical multipole expansion")
    CASE(CYLINMULTEXPANFAST)
       ! initialize geometry
       CALL InitGeometry(geometry,CYLINDRICAL)
       CALL Info(this, " SOURCES--> boundary:          cylindrical multipole expansion (fast)")
    CASE DEFAULT
       CALL Error(this, "InitSources_boundary", "Wrong BOUNDARYTYPE") 
    END SELECT

    this%BOUNDARYTYPE = bndrytype
    this%MAXMULT = maxmult

    DO i=1,4
       SELECT CASE(GetType(Boundary(i)))
       CASE(REFLECTING,AXIS)
          this%Boundary(i) = NEUMANN
          CALL Info(this, " SOURCES--> boundary:          "//TRIM(GetDirectionName(Boundary(i)))//" NEUMANN boundary condition")
       CASE(PERIODIC)
          this%Boundary(i) = PERIODIC
          CALL Info(this, " SOURCES--> boundary:          "//TRIM(GetDirectionName(Boundary(i)))//" PERIODIC boundary condition")
       CASE DEFAULT
          this%Boundary(i) = DIRICHLET
          CALL Info(this, " SOURCES--> boundary:          "//TRIM(GetDirectionName(Boundary(i)))//" DIRICHLET boundary condition")
       END SELECT
    END DO

    ngrid = this%ngrid
       ni = this%grid(ngrid)%ni
       nj = this%grid(ngrid)%nj

       ALLOCATE(this%grid(ngrid)%curv(ni,nj,2),&
                STAT = err)
       IF (err.NE.0) CALL Error(this, "InitSources_boundary", "Unable to allocate memory.")

       CALL Convert2Curvilinear(geometry,this%grid(ngrid)%bccart,&
                                         this%grid(ngrid)%curv)
IF (GetType(geometry) .EQ. CYLINDRICAL) this%grid(ngrid)%curv(:,:,2) = ABS(this%grid(ngrid)%curv(:,:,2))
       
       IF (this%BOUNDARYTYPE == SPHERMULTEXPANFAST) then
          ALLOCATE(this%grid(ngrid)%mlint(ni,0:this%MAXMULT), &
                   this%grid(ngrid)%mlext(ni,0:this%MAXMULT), &
                   STAT = err)
          IF (err.NE.0) CALL Error(this, "InitSources_boundary", "Unable to allocate memory.")
       END IF

       IF (this%BOUNDARYTYPE == CYLINMULTEXPANFAST) THEN
          ALLOCATE(this%grid(ngrid)%QdivsqrtR(ni,nj,2*ni+2*nj-4), &
                   this%grid(ngrid)%ij2k(ni,nj), &
                   this%grid(ngrid)%k2ij(2*ni+2*nj-4,2), &
                   STAT = err)
          IF (err.NE.0) CALL Error(this, "InitSources_boundary", "Unable to allocate memory.")

          this%grid(ngrid)%ij2k(:,:) = -1
          this%grid(ngrid)%k2ij(:,:) = -1
          DO k = 1, nj
             this%grid(ngrid)%ij2k(1,k)  = k
             this%grid(ngrid)%ij2k(ni,k) = nj+k
             this%grid(ngrid)%k2ij(k,1) = 1
             this%grid(ngrid)%k2ij(k,2) = k
             this%grid(ngrid)%k2ij(nj+k,1) = ni
             this%grid(ngrid)%k2ij(nj+k,2) = k
          END DO
          DO k = 2, ni-1
             this%grid(ngrid)%ij2k(k,1)  = 2*nj+k-1
             this%grid(ngrid)%ij2k(k,nj) = 2*nj+ni+k-3
             this%grid(ngrid)%k2ij(2*nj+k-1,1) = k
             this%grid(ngrid)%k2ij(2*nj+k-1,2) = 1
             this%grid(ngrid)%k2ij(2*nj+ni+k-3,1) = k
             this%grid(ngrid)%k2ij(2*nj+ni+k-3,2) = nj
          END DO

!CDIR UNROLL=8
          DO k = 1, 2*ni+2*nj-4
             R1 = this%grid(ngrid)%curv(this%grid(ngrid)%k2ij(k,1),this%grid(ngrid)%k2ij(k,2),2)
             z1 = this%grid(ngrid)%curv(this%grid(ngrid)%k2ij(k,1),this%grid(ngrid)%k2ij(k,2),1)
!CDIR COLLAPSE
             DO j = 2, nj-1
              DO i = 2, ni-1
                R2 = this%grid(ngrid)%curv(i,j,2)
                z2 = this%grid(ngrid)%curv(i,j,1)
                chi = (R1**2+R2**2+(z1-z2)**2)/(2.0*R1*R2)
                !mu = sqrt(4.0*R1*R2/(((R1+R2)**2+(z1-z2)**2)+Tiny(1.0)))
                mu = sqrt(2.0/(1.0+chi))
                this%grid(ngrid)%QdivsqrtR(i,j,k) = &
                           !SumQ(this,0,mu,chi)& !SumQ(..,0,..) = mu kell_agm
                           mu*kell_agm(this,mu)& 
                           / sqrt(R2)&
                           *this%grid(ngrid)%vol(i,j)
             END DO
           END DO
         END DO
      END IF
  END SUBROUTINE InitSources_boundary


  PURE SUBROUTINE CalcBoundary(this)
  IMPLICIT NONE
  !------------------------------------------------------------------------!
  TYPE(Sources_TYP)                :: this
  !------------------------------------------------------------------------!
  INTEGER                          :: bnd
  INTEGER, DIMENSION(4)            :: rg
  !------------------------------------------------------------------------!
  INTENT(INOUT)                    :: this
  !------------------------------------------------------------------------!

  IF (this%BOUNDARYTYPE == SPHERMULTEXPANFAST) &
     CALL calcSPHERMULTEXPANFAST_INIT(this)

  DO bnd=1,4
     IF (this%Boundary(bnd) .EQ. DIRICHLET) THEN
        SELECT CASE(bnd)
        CASE(WEST)
           rg(1) = 1
           rg(2) = 1
           rg(3) = 1
           rg(4) = this%grid(this%ngrid)%nj
        CASE(EAST)
           rg(1) = this%grid(this%ngrid)%ni
           rg(2) = this%grid(this%ngrid)%ni
           rg(3) = 1
           rg(4) = this%grid(this%ngrid)%nj
        CASE(SOUTH)
           rg(1) = 1
           rg(2) = this%grid(this%ngrid)%ni
           rg(3) = 1
           rg(4) = 1
        CASE(NORTH)
           rg(1) = 1
           rg(2) = this%grid(this%ngrid)%ni
           rg(3) = this%grid(this%ngrid)%nj
           rg(4) = this%grid(this%ngrid)%nj
        END SELECT

        SELECT CASE (this%BOUNDARYTYPE)
        CASE(SPHERMULTEXPAN)
           CALL calcSPHERMULTEXPAN(this,rg)
        CASE(CYLINMULTEXPAN)
           CALL calcCYLINMULTEXPAN(this,rg)
        CASE(SPHERMULTEXPANFAST)
           CALL calcSPHERMULTEXPANFAST(this,rg)  
        CASE(CYLINMULTEXPANFAST)
           CALL calcCYLINMULTEXPANFAST(this,rg)   
        END SELECT
     END IF
  END DO

  END SUBROUTINE CalcBoundary


SUBROUTINE CloseSources_boundary(this)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Sources_TYP) :: this
   !------------------------------------------------------------------------!
   INTEGER           :: j
   !------------------------------------------------------------------------!
   INTENT(INOUT)     :: this
   !------------------------------------------------------------------------!
       j = this%ngrid
       DEALLOCATE(this%grid(j)%curv)
       IF (this%BOUNDARYTYPE == SPHERMULTEXPANFAST) THEN
          DEALLOCATE(this%grid(j)%mlext)
          DEALLOCATE(this%grid(j)%mlint)
       ELSE IF (this%BOUNDARYTYPE == CYLINMULTEXPANFAST) THEN
          DEALLOCATE(this%grid(j)%QdivsqrtR)
          DEALLOCATE(this%grid(j)%ij2k)
          DEALLOCATE(this%grid(j)%k2ij)
       END IF
  END SUBROUTINE CloseSources_boundary


PURE SUBROUTINE calcSPHERMULTEXPANFAST_INIT(this)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Sources_TYP)       :: this
   !------------------------------------------------------------------------!
   INTENT(INOUT)           :: this  
   !------------------------------------------------------------------------!
   TYPE(Grid_TYP), POINTER :: pgrid
   INTEGER                 :: ii,jj,l
   !------------------------------------------------------------------------!

      pgrid=>this%grid(this%ngrid)
      pgrid%mlint(:,:) = 0.0
      pgrid%mlext(:,:) = 0.0
      DO l = 0, this%MAXMULT
!CDIR COLLAPSE
         DO ii = 2, pgrid%ni-1
            DO jj = 2, pgrid%nj-1
               pgrid%mlext(ii,l) = pgrid%mlext(ii,l) + 2.0*PI&
                  *pgrid%rho(ii,jj) * pgrid%curv(ii,jj,1)**l&
                  *plgndr_s(l,cos(pgrid%curv(ii,jj,2)))&
                  *pgrid%vol(ii,jj)

               pgrid%mlint(ii,l) = pgrid%mlint(ii,l) + 2.0*PI&
                      *pgrid%rho(ii,jj) / (pgrid%curv(ii,jj,1)**(l+1)+Tiny(1.0))&
                      *plgndr_s(l,cos(pgrid%curv(ii,jj,2)))&
                      *pgrid%vol(ii,jj)
            END DO
         END DO
       END DO
END SUBROUTINE calcSPHERMULTEXPANFAST_INIT

PURE SUBROUTINE calcSPHERMULTEXPANFAST(this,rg)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Sources_TYP)               :: this
   INTEGER,DIMENSION(4)            :: rg
   !------------------------------------------------------------------------!
   TYPE(Grid_TYP), POINTER         :: pgrid
   REAL, PARAMETER                 :: invnorm = 1.0/(4.0*PI)
   INTEGER                         :: i,j,l
   !------------------------------------------------------------------------!
   INTENT(INOUT)                   :: this
   INTENT(IN)                      :: rg
   !------------------------------------------------------------------------!
   pgrid => this%grid(this%ngrid)
!CDIR COLLAPSE 
   DO i = rg(1),rg(2)
!CDIR NODEP
      DO j = rg(3),rg(4)
         pgrid%u(i,j) = 0.0
         DO l = 0, this%MAXMULT
            pgrid%u(i,j) = pgrid%u(i,j) - plgndr_s(l,cos(pgrid%curv(i,j,2)))&
                           *(1.0/(pgrid%curv(i,j,1)**(l+1)+Tiny(1.0))*&
                      SUM(pgrid%mlext(:,l),&
                          pgrid%curv(:,j,1) < pgrid%curv(i,j,1))&
                         +pgrid%curv(i,j,1)**l*&
                      SUM(pgrid%mlint(:,l),&
                          pgrid%curv(:,j,1) >= pgrid%curv(i,j,1))&
                          )*invnorm
         END DO
      END DO
   END DO
 END SUBROUTINE calcSPHERMULTEXPANFAST

 PURE SUBROUTINE calcCYLINMULTEXPAN(this,rg)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Sources_TYP)               :: this
   INTEGER,DIMENSION(4)            :: rg
   !------------------------------------------------------------------------!
   INTENT(INOUT)                   :: this
   INTENT(IN)                      :: rg
   !------------------------------------------------------------------------!
   TYPE(Grid_TYP), POINTER         :: pgrid
   REAL                            :: R1,R2,z1,z2,mu,chi
   INTEGER                         :: i,j,ii,jj,ni,nj
   !------------------------------------------------------------------------!
   pgrid=>this%grid(this%ngrid)
   ni = pgrid%ni
   nj = pgrid%nj
!CDIR COLLAPSE   
   DO i = rg(1),rg(2)
!CDIR NODEP
      DO j = rg(3),rg(4)
         pgrid%u(i,j) = 0.0  
         R1 = pgrid%curv(i,j,2)
         z1 = pgrid%curv(i,j,1)
         !CDIR COLLAPSE
         DO ii = 2,ni-1
            DO jj = 2,nj-1
               R2 = pgrid%curv(ii,jj,2)
               z2 = pgrid%curv(ii,jj,1)
               chi = (R1**2+R2**2+(z1-z2)**2)/(2.0*R1*R2)
               !mu = sqrt(4.0*R1*R2/((R1+R2)**2+(z1-z2)**2+Tiny(1.0)))
               mu = sqrt(2.0/(1.0+chi))
               pgrid%u(i,j) = pgrid%u(i,j)-pgrid%rho(ii,jj)*pgrid%vol(ii,jj)&
                       /sqrt(R2) &
                       !* SumQ(this,0,mu,chi)
                       *mu*kell_agm(this,mu)
            END DO
         END DO
         pgrid%u(i,j) = pgrid%u(i,j)/(sqrt(R1)*2.0*PI)
      END DO
   END DO
 END SUBROUTINE calcCYLINMULTEXPAN

PURE SUBROUTINE calcCYLINMULTEXPANFAST(this,rg)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Sources_TYP)               :: this
   INTEGER,DIMENSION(4)            :: rg
   !------------------------------------------------------------------------!
   INTENT(INOUT)                   :: this
   INTENT(IN)                      :: rg
   !------------------------------------------------------------------------!
   TYPE(Grid_TYP), POINTER         :: pgrid
   INTEGER                         :: i,j,ni,nj
   REAL                            :: R1
   !------------------------------------------------------------------------!
   pgrid=>this%grid(this%ngrid)
   ni = pgrid%ni
   nj = pgrid%nj

!CDIR COLLAPSE
   DO i = rg(1),rg(2)
!CDIR NODEP
      DO j = rg(3),rg(4)
         R1 = pgrid%curv(i,j,2)
         pgrid%u(i,j) = -SUM(pgrid%rho(2:ni-1,2:nj-1)&
          *pgrid%QdivsqrtR(2:ni-1,2:nj-1,pgrid%ij2k(i,j)))&
          /(sqrt(R1)*2.0*PI)
      END DO
   END DO

END SUBROUTINE calcCYLINMULTEXPANFAST

PURE SUBROUTINE calcSPHERMULTEXPAN(this,rg)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Sources_TYP)               :: this
   INTEGER,DIMENSION(4)            :: rg
   !------------------------------------------------------------------------!
   TYPE(Grid_TYP), POINTER         :: pgrid
   REAL, PARAMETER                 :: invnorm = 1.0/(4.0*PI)
   INTEGER                         :: i,j
   !------------------------------------------------------------------------!
   INTENT(INOUT)                   :: this
   INTENT(IN)                      :: rg
   !------------------------------------------------------------------------!
   pgrid=>this%grid(this%ngrid)
!CDIR COLLAPSE   
   DO i = rg(1),rg(2)
!CDIR NODEP
      DO j = rg(3),rg(4)
        pgrid%u(i,j) = (boundarymultext(this,i,j)+boundarymultint(this,i,j))*invnorm
      END DO
   END DO
END SUBROUTINE calcSPHERMULTEXPAN

PURE FUNCTION boundarymultext(this,i,j)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Sources_TYP)     :: this
   INTEGER               :: i,j
   !------------------------------------------------------------------------!
   INTEGER               :: ni,nj,ii,jj,l
   REAL                  :: boundarymultext,Ml
   !------------------------------------------------------------------------!
   INTENT(IN)            :: this,i,j
   !------------------------------------------------------------------------!
   ni = this%grid(this%ngrid)%ni
   nj = this%grid(this%ngrid)%nj

   boundarymultext = 0.0
  
   DO l = 0, this%MAXMULT
      Ml = 0.0
!CDIR COLLAPSE
      DO ii = 2, ni-1
         DO jj = 2, nj-1
            IF ( this%grid(this%ngrid)%curv(ii,jj,1) < (this%grid(this%ngrid)%curv(i,j,1)) ) THEN
               Ml = Ml + 2.0*PI&
                   * this%grid(this%ngrid)%rho(ii,jj) * this%grid(this%ngrid)%curv(ii,jj,1)**l&
                   *plgndr_s(l,cos(this%grid(this%ngrid)%curv(ii,jj,2)))&
                   *this%grid(this%ngrid)%vol(ii,jj)
            END IF
         END DO
      END DO
      boundarymultext = boundarymultext - plgndr_s(l,cos(this%grid(this%ngrid)%curv(i,j,2)))&
                        *Ml/(this%grid(this%ngrid)%curv(i,j,1)+Tiny(1.0))**(l+1)
   END DO
END FUNCTION boundarymultext

PURE FUNCTION boundarymultint(this,i,j)
   IMPLICIT NONE
   !------------------------------------------------------------------------!   
   TYPE(Sources_TYP)   :: this
   INTEGER             :: i,j
   !------------------------------------------------------------------------!
   INTEGER             :: ni,nj, ii, jj, l
   REAL                :: boundarymultint, Ml
   !------------------------------------------------------------------------!
   INTENT(IN)          :: this,i,j 
   !------------------------------------------------------------------------!
   ni = this%grid(this%ngrid)%ni
   nj = this%grid(this%ngrid)%nj

   boundarymultint = 0.0

   DO l = 0, this%MAXMULT
      Ml = 0.0
!CDIR COLLAPSE
      DO ii = 2, ni-1
         DO jj = 2, nj-1
            IF ( this%grid(this%ngrid)%curv(ii,jj,1) >= (this%grid(this%ngrid)%curv(i,j,1)) ) THEN
               Ml = Ml + 2.0*PI &
                    *this%grid(this%ngrid)%rho(ii,jj) / (this%grid(this%ngrid)%curv(ii,jj,1)**(l+1)+Tiny(1.0))&
                    *plgndr_s(l,cos(this%grid(this%ngrid)%curv(ii,jj,2)))&
                    *this%grid(this%ngrid)%vol(ii,jj)
            END IF
         END DO
      END DO
      boundarymultint = boundarymultint - plgndr_s(l,cos(this%grid(this%ngrid)%curv(i,j,2)))&
                       *Ml*this%grid(this%ngrid)%curv(i,j,1)**l
   END DO

 END FUNCTION boundarymultint

 PURE FUNCTION SumQ(this,Mmax,mu,chi)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    INTEGER           :: Mmax
    REAL              :: mu, chi
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mmax,mu,chi
    !------------------------------------------------------------------------!
    REAL              :: SumQ,Kell,Qm_2,Qm_1,Qm
    INTEGER           :: m
    !------------------------------------------------------------------------!

    IF (abs(mu) >= 1.0) THEN
      SumQ = 0.0
      RETURN
    END IF

    Kell = kell_agm(this,mu)
    SumQ = mu * Kell
    IF (Mmax == 0) RETURN
    Qm_2 = SumQ
    Qm_1 = chi * Qm_2 - (1.0 + chi)*mu*Eell_agm(this,Kell,mu)
    SumQ = SumQ + 2.0*Qm_1
    DO m = 2, Mmax
       Qm = 4.0 * (m-1.0)/(2.0*m-1.0) * chi * Qm_1 - (2.0*m-3.0)/(2.0*m-1.0)*Qm_2
       SumQ = SumQ + 2.0*Qm
       Qm_2 = Qm_1
       Qm_1 = Qm
    END DO    
 END FUNCTION SumQ

PURE FUNCTION kell_agm(this,z)
  IMPLICIT NONE
  !------------------------------------------------------------------------!
  TYPE(Sources_TYP)  :: this
  REAL               :: z
  !------------------------------------------------------------------------!
  INTENT(IN)         :: this,z
  !------------------------------------------------------------------------!
  REAL               :: kell_agm,an,bn
  !------------------------------------------------------------------------!
   IF (abs(z) >= 1.0) THEN
      kell_agm = 0.0
      RETURN
   END IF
   an = 1.0
   bn = sqrt(1-z**2)
   DO WHILE (abs(an-bn) > this%MAXAGMNORM)
      kell_agm = 0.5*(an + bn)
      bn = sqrt(an*bn)
      an = kell_agm
   END DO    
   kell_agm = 0.5*PI / an
END FUNCTION kell_agm

PURE FUNCTION Eell_agm(this,Kz,z)
IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Sources_TYP)    :: this
   REAL                 :: Kz,z
   !------------------------------------------------------------------------!
   INTENT(IN)           :: this,Kz,z
   !------------------------------------------------------------------------!
   REAL                 :: Eell_agm,an,bn,cn,tmp
   INTEGER              :: i2
   !------------------------------------------------------------------------!
   IF (z == 1.0) THEN
      Eell_agm = 0.0
      RETURN
   END IF
   an = 1.0
   bn = sqrt(1-z**2)
   cn = z
   Eell_agm = 0.5 * cn**2
   i2 = 1
   DO WHILE ( abs(cn) > this%MAXAGMNORM )
      cn = 0.5*(an-bn)
      tmp = 0.5*(an+bn)
      bn = sqrt(an*bn)
      an = tmp
      i2 = i2*2
      Eell_agm = Eell_agm + 0.5*i2*cn**2
   END DO   
   Eell_agm = Kz*(1.0-Eell_agm)
END FUNCTION Eell_agm

!Test (from Numerical Recipes in Fortran 90)
PURE FUNCTION plgndr_s(l,x)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   INTEGER       :: l
   REAL          :: x
   !------------------------------------------------------------------------!
   INTENT(IN)    :: l,x
   !------------------------------------------------------------------------!
   REAL :: plgndr_s
   INTEGER :: ll
   REAL :: pll,pmm,pmmp1
   !------------------------------------------------------------------------!
   !Computes the ordinary Legendre polynomial Pl (x). Here l is a integer, 
   !while x lies in the range −1 ≤ x ≤ 1.
   if (abs(x) > 1.0) then
   !TODO: call error   
   end if

   pmm=1.0
    if (l == 0) then
       plgndr_s=pmm
    else
        pmmp1=x*pmm
       if (l == 1) then
          plgndr_s=pmmp1
       else
          do ll=2,l
             pll=(x*(2*ll-1)*pmmp1-(ll-1)*pmm)/ll
             pmm=pmmp1
             pmmp1=pll
          end do
          plgndr_s=pll
       end if
    end if
END FUNCTION plgndr_s

END MODULE sources_boundary
