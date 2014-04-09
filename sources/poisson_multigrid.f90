!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: poisson_multigrid.f90                                             #
!#                                                                           #
!# Copyright (C) 2011                                                        #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
! poisson solver via multigrid
!----------------------------------------------------------------------------!
MODULE poisson_multigrid
  USE poisson_common
  USE multipole_generic
  USE mesh_generic
  USE physics_generic
  USE boundary_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: solver_name  = "multigrid"
  INTEGER, PARAMETER :: DIRICHLET              = 1001
  INTEGER, PARAMETER :: NEUMANN                = 1002
  INTEGER, PARAMETER :: RED_BLACK_GAUSS_SEIDEL = 1
  INTEGER, PARAMETER :: BLOCK_GAUSS_SEIDEL     = 2
  INTEGER, PARAMETER :: GAUSS_SEIDEL           = 3
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Poisson_TYP, &
       Grid_TYP, &
       ! constants
       SPHERMULTEXPAN, CYLINMULTEXPAN, &
       RED_BLACK_GAUSS_SEIDEL,BLOCK_GAUSS_SEIDEL,GAUSS_SEIDEL,&
       ! methods
       InitPoisson_multigrid, &
       CalcPotential_multigrid, &
       ClosePoisson_multigrid,&
       GetType, &
       GetName, &
       GetRank, &
       GetNumProcs, &
       Initialized, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!
  
  CONTAINS

  SUBROUTINE InitPoisson_multigrid(this,Mesh,Physics,Boundary,solver,&
                       maxmult,maxresidnorm,bndrytype,relaxtype,&
                       npre,npost,minres,nmaxcycle)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Poisson_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Boundary_TYP), DIMENSION(4) :: Boundary
    INTEGER           :: solver
    !------------------------------------------------------------------------!
    REAL              :: maxresidnorm
    INTEGER           :: maxmult,bndrytype,relaxtype,npre,npost,minres,nmaxcycle
    !------------------------------------------------------------------------!
    TYPE(Selection_TYP) :: iregion
    TYPE(Selection_TYP), ALLOCATABLE :: bnd_region(:)
    INTEGER           :: err,i,j,k,ni,nj,ni0,nj0,jgrid
    CHARACTER(LEN=48)   :: xres_str,yres_str,grid_str,relaxtype_str
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Boundary,solver,maxmult,maxresidnorm,bndrytype,&
                         npre,npost,minres,nmaxcycle
    !------------------------------------------------------------------------!
     CALL InitPoisson(this,solver,solver_name)
#ifdef PARALLEL
    CALL Error(this,"InitPoisson_multigrid", &
         "parallel version has not been implemented for multigrid")
#endif  
    this%MAXRESIDNORM = maxresidnorm
    this%RELAXTYPE    = relaxtype
    IF (npre+npost .LT. 1) CALL Error(this,"InitPoisson_multigrid","NPRE+NPOST >= 1")
    this%NPRE         = npre
    this%NPOST        = npost
    this%minres       = minres
    this%nmaxcycle    = nmaxcycle

    ni = Mesh%INUM
    nj = Mesh%JNUM

    !estimate multigrid steps
    this%NGRID=0
    DO 
       this%NGRID=this%NGRID+1
       ni0=ni
       nj0=nj
       CALL coarsergrid(this%MINRES,ni,nj)
       IF (ni0 .EQ. ni .AND. nj0 .EQ. nj) EXIT
    END DO
    ALLOCATE(this%grid(this%NGRID))

    ! print some information
    SELECT CASE(this%RELAXTYPE)
    CASE(RED_BLACK_GAUSS_SEIDEL)
       WRITE (relaxtype_str, '(A)') "with red-black Gauss-Seidel iteration"
    CASE(BLOCK_GAUSS_SEIDEL)
       WRITE (relaxtype_str, '(A)') "with Gauss-Seidel block iteration (best)"
    CASE(GAUSS_SEIDEL)
       WRITE (relaxtype_str, '(A)') "with Gauss-Seidel iteration (slowest)"
    END SELECT
    WRITE (xres_str, '(I0)') ni0
    WRITE (yres_str, '(I0)') nj0
    WRITE (grid_str, '(I0)') this%NGRID
    CALL Info(this, " POISSON--> multigrid method   " // TRIM(relaxtype_str))
    CALL Info(this, "            grid levels:       " // TRIM(grid_str) )
    CALL Info(this, "            coarsest grid res. " // TRIM(xres_str) // " x " // TRIM(yres_str))
        
IF (ni0 .GT. 19 .OR. nj0 .GT. 19) CALL Warning(this, "InitPoisson_multigrid",&
                "bad multigrid convergence: choose a better resolution e.g. m*2^n+1 with small natural number m")

!  !FIXME convergence test
!  ALLOCATE(this%relaxcount(this%ngrid),STAT = err)
!  IF (err.NE.0) CALL Error(this, "InitPoisson_multigrid", "Unable to allocate memory.")
!  this%relaxcount(:) = 0
!  this%safedmultigrid = 0

    ni = Mesh%INUM
    nj = Mesh%JNUM
    DO jgrid=1, this%NGRID
       ALLOCATE(this%grid(jgrid)%u(0:ni+1,0:nj+1), &
                this%grid(jgrid)%rho(0:ni+1,0:nj+1), &
                this%grid(jgrid)%a(0:ni+1,0:nj+1), &
                this%grid(jgrid)%da(1:ni,1:nj), &
                this%grid(jgrid)%b(0:ni+1,0:nj+1), &
                this%grid(jgrid)%db(1:ni,1:nj), &
                this%grid(jgrid)%c(0:ni+1,0:nj+1), &
                this%grid(jgrid)%invc(0:ni+1,0:nj+1), &
                this%grid(jgrid)%d(0:ni+1,0:nj+1), &
                this%grid(jgrid)%vol(0:ni+1,0:nj+1), &
                this%grid(jgrid)%bccart(0:ni+1,0:nj+1,2), &
                this%grid(jgrid)%bhx(0:ni+1,0:nj+1), &
                this%grid(jgrid)%bhy(0:ni+1,0:nj+1), &
                this%grid(jgrid)%tmp(0:ni+1,0:nj+1), &
                STAT = err)
       IF (err.NE.0) CALL Error(this, "InitPoisson_multigrid", "Unable to allocate memory.")

       this%grid(jgrid)%u(:,:) = 0.0

       this%grid(jgrid)%hi=(Mesh%xmax-Mesh%xmin)/ni
       this%grid(jgrid)%hj=(Mesh%ymax-Mesh%ymin)/nj
       this%grid(jgrid)%ni=ni
       this%grid(jgrid)%nj=nj
       this%grid(jgrid)%invhi2 = 1.0 / this%grid(jgrid)%hi**2
       this%grid(jgrid)%invhj2 = 1.0 / this%grid(jgrid)%hj**2
       IF (jgrid == 1) THEN
          !only at first grid
          this%grid(jgrid)%bccart(0:ni+1,0:nj+1,:)=&
                Mesh%bccart(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1,:)
          this%grid(jgrid)%bhx(0:ni+1,0:nj+1)=&
                     Mesh%bhx(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)
          this%grid(jgrid)%bhy(0:ni+1,0:nj+1)=&
                     Mesh%bhy(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)
!FIXME: reasonable? or bhz and than direct calculation of a,b,c
          this%grid(jgrid)%a(0:ni+1,0:nj+1)=&
                     Mesh%bhy(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)&
                    *Mesh%bhz(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)&
                    /Mesh%bhx(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)
          this%grid(jgrid)%b(0:ni+1,0:nj+1)=&
                     Mesh%bhx(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)&
                    *Mesh%bhz(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)&
                    /Mesh%bhy(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)
          this%grid(jgrid)%c(0:ni+1,0:nj+1)=&
                     Mesh%bhx(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)&
                    *Mesh%bhy(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)&
                    *Mesh%bhz(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)
          !important: first restrict step needs vol
          this%grid(jgrid)%vol(0:ni+1,0:nj+1)= &
                     Mesh%volume(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)
       ELSE
!********************************************************************************
!FIXME: What is better? Restriction or copy?
          !Restriction
          CALL restrict(this%grid(jgrid-1)%bhx,&
                        this%grid(jgrid-1)%ni,this%grid(jgrid-1)%nj,ni,nj,&
                        this%grid(jgrid)%bhx)
          CALL restrict(this%grid(jgrid-1)%bhy,&
                        this%grid(jgrid-1)%ni,this%grid(jgrid-1)%nj,ni,nj,&
                        this%grid(jgrid)%bhy)
          CALL restrict(this%grid(jgrid-1)%a,&
                        this%grid(jgrid-1)%ni,this%grid(jgrid-1)%nj,ni,nj,&
                        this%grid(jgrid)%a)
          CALL restrict(this%grid(jgrid-1)%b,&
                        this%grid(jgrid-1)%ni,this%grid(jgrid-1)%nj,ni,nj,&
                        this%grid(jgrid)%b)
          CALL restrict(this%grid(jgrid-1)%c,&
                        this%grid(jgrid-1)%ni,this%grid(jgrid-1)%nj,ni,nj,&
                        this%grid(jgrid)%c)
           this%grid(jgrid)%vol(:,:) = this%grid(jgrid)%c(:,:)&
                           * this%grid(jgrid)%hi *this%grid(jgrid)%hj


!          ni0 = this%grid(jgrid-1)%ni
!          nj0 = this%grid(jgrid-1)%nj
!          !COPY
!          this%grid(jgrid)%bhx(0:ni+1,0:nj+1)= this%grid(jgrid-1)%bhx(0:ni0+1:2,0:nj0+1:2)
!          this%grid(jgrid)%bhy(0:ni+1,0:nj+1)= this%grid(jgrid-1)%bhy(0:ni0+1:2,0:nj0+1:2)
!          this%grid(jgrid)%a(0:ni+1,0:nj+1)  = this%grid(jgrid-1)%a(0:ni0+1:2,0:nj0+1:2)
!          this%grid(jgrid)%b(0:ni+1,0:nj+1)  = this%grid(jgrid-1)%b(0:ni0+1:2,0:nj0+1:2)
!          this%grid(jgrid)%c(0:ni+1,0:nj+1)  = this%grid(jgrid-1)%c(0:ni0+1:2,0:nj0+1:2)
!          this%grid(jgrid)%vol(:,:)          = this%grid(jgrid-1)%vol(0:ni0+1:2,0:nj0+1:2)
!********************************************************************************


          !boundary volumes of boundary cells
          this%grid(jgrid)%vol(0,:) = this%grid(jgrid)%c(0,:)&
                            * this%grid(1)%hi *this%grid(jgrid)%hj
          this%grid(jgrid)%vol(ni+1,:) = this%grid(jgrid)%c(ni+1,:)&
                            * this%grid(1)%hi *this%grid(jgrid)%hj
          this%grid(jgrid)%vol(:,0) = this%grid(jgrid)%c(:,0)&
                            * this%grid(jgrid)%hi *this%grid(1)%hj
          this%grid(jgrid)%vol(:,nj+1) = this%grid(jgrid)%c(:,nj+1)&
                            * this%grid(jgrid)%hi *this%grid(1)%hj
          !boundary volumes of corner cells
          this%grid(jgrid)%vol(0,0) = this%grid(jgrid)%c(0,0)&
                            * this%grid(1)%hi *this%grid(1)%hj
          this%grid(jgrid)%vol(ni+1,0) = this%grid(jgrid)%c(ni+1,0)&
                            * this%grid(1)%hi *this%grid(1)%hj
          this%grid(jgrid)%vol(ni+1,nj+1) = this%grid(jgrid)%c(ni+1,nj+1)&
                            * this%grid(1)%hi *this%grid(1)%hj
          this%grid(jgrid)%vol(0,nj+1) = this%grid(jgrid)%c(0,nj+1)&
                            * this%grid(1)%hi *this%grid(1)%hj

       END IF

       !FIXME*** 
       this%grid(jgrid)%vol(:,:) = abs(this%grid(jgrid)%vol(:,:))
       !***
       this%grid(jgrid)%tmp(:,:) = 0.0
       this%grid(jgrid)%invc(:,:) = 1.0/this%grid(jgrid)%c(:,:)
       this%grid(jgrid)%d(:,:) = &
                     1.0/( this%grid(jgrid)%a(:,:)*this%grid(jgrid)%invhi2 &
                          +this%grid(jgrid)%b(:,:)*this%grid(jgrid)%invhj2 )
       this%grid(jgrid)%da(1:ni,1:nj) = &
                     this%grid(jgrid)%a(2:ni+1,1:nj)-this%grid(jgrid)%a(0:ni-1,1:nj)
       this%grid(jgrid)%db(1:ni,1:nj) = &
                     this%grid(jgrid)%b(1:ni,2:nj+1)-this%grid(jgrid)%b(1:ni,0:nj-1)

       IF (this%relaxtype .EQ. BLOCK_GAUSS_SEIDEL) THEN
         ALLOCATE(this%grid(jgrid)%tri(0:ni+1,0:nj+1,9),&
                  STAT = err)
         IF (err.NE.0) CALL Error(this, "InitPoisson_multigrid", "Unable to allocate memory.")
         !1. - 5. matrix elements (diag,tri diag, elements); 6. rhs; 7.-9. temp variables (hb,hv,m) 

         this%grid(jgrid)%tri(1:ni,1:nj,1) = -2.0*(this%grid(jgrid)%invhi2*this%grid(jgrid)%a(1:ni,1:nj)&
                                             +this%grid(jgrid)%invhj2*this%grid(jgrid)%b(1:ni,1:nj))
         this%grid(jgrid)%tri(1:ni,1:nj,2) = this%grid(jgrid)%invhj2*(this%grid(jgrid)%b(1:ni,1:nj)&
                                             -this%grid(jgrid)%db(1:ni,1:nj)*0.25)
         this%grid(jgrid)%tri(1:ni,1:nj,3) = this%grid(jgrid)%invhj2*(this%grid(jgrid)%b(1:ni,1:nj)&
                                             +this%grid(jgrid)%db(1:ni,1:nj)*0.25)     
         this%grid(jgrid)%tri(1:ni,1:nj,4) = this%grid(jgrid)%invhi2*(this%grid(jgrid)%a(1:ni,1:nj)&
                                             -this%grid(jgrid)%da(1:ni,1:nj)*0.25)
         this%grid(jgrid)%tri(1:ni,1:nj,5) = this%grid(jgrid)%invhi2*(this%grid(jgrid)%a(1:ni,1:nj)&
                                             +this%grid(jgrid)%da(1:ni,1:nj)*0.25)
       END IF
       CALL coarsergrid(this%MINRES,ni,nj)
    END DO

    ! boundary conditions for poisson solver
    k = 0
    DO i=1,4
       SELECT CASE(GetType(Boundary(i)))
       CASE(REFLECTING,AXIS,NOSLIP)
          this%Boundary(i) = NEUMANN
       CASE(PERIODIC)
          this%Boundary(i) = PERIODIC
       CASE DEFAULT
          this%Boundary(i) = DIRICHLET
          k = k + 1  ! count number of DIRICHLET boundary conditions
       END SELECT
    END DO

    ! initialize multipole expansion for DIRICHLET boundary conditions
    IF (k.GT.0) THEN
       ! set input region indices for multipole expansion
       iregion%IMIN = 1
       iregion%IMAX = this%grid(1)%ni 
       iregion%JMIN = 1
       iregion%JMAX = this%grid(1)%nj

       ! allocate array for output selection indices
       ALLOCATE(bnd_region(k),STAT = err)
       IF (err.NE.0) CALL Error(this, "InitPoisson_multigrid", &
            "Unable to allocate memory.")

       ! set output selection indices for DIRICHLET boundary conditions
       j = 0
       IF (this%Boundary(WEST).EQ.DIRICHLET) THEN
          j = j + 1
          bnd_region(j)%IMIN = 0
          bnd_region(j)%IMAX = 0
          bnd_region(j)%JMIN = 0
          bnd_region(j)%JMAX = this%grid(1)%nj+1
       END IF
       IF (this%Boundary(EAST).EQ.DIRICHLET) THEN
          j = j + 1
          bnd_region(j)%IMIN = this%grid(1)%ni+1
          bnd_region(j)%IMAX = this%grid(1)%ni+1
          bnd_region(j)%JMIN = 0
          bnd_region(j)%JMAX = this%grid(1)%nj+1
       END IF
       IF (this%Boundary(SOUTH).EQ.DIRICHLET) THEN
          j = j + 1
          bnd_region(j)%IMIN = 0
          bnd_region(j)%IMAX = this%grid(1)%ni+1
          bnd_region(j)%JMIN = 0
          bnd_region(j)%JMAX = 0
       END IF
       IF (this%Boundary(NORTH).EQ.DIRICHLET) THEN
          j = j + 1
          bnd_region(j)%IMIN = 0
          bnd_region(j)%IMAX = this%grid(1)%ni+1
          bnd_region(j)%JMIN = this%grid(1)%nj+1
          bnd_region(j)%JMAX = this%grid(1)%nj+1
       END IF

       ! initialize the multipole expansion module
       CALL InitMultipole(this%multipole,bndrytype,this%grid(1)%bccart, &
            2*PI*this%grid(1)%vol,0,this%grid(1)%ni+1, &
            0,this%grid(1)%nj+1,iregion,bnd_region,maxmult)
    END IF
  END SUBROUTINE InitPoisson_multigrid

  PURE SUBROUTINE coarsergrid(MINRES,ni,nj)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER, INTENT(IN)     :: MINRES
    INTEGER, INTENT(INOUT)  :: ni,nj
    !------------------------------------------------------------------------!
    !coarsing in every step if ni is odd
    IF (mod(ni,2) .EQ. 1 .AND. (ni-1)/2 .GE. MINRES) ni=(ni-1)/2
    !coarsing in every step if nj is odd
    IF (mod(nj,2) .EQ. 1 .AND. (nj-1)/2 .GE. MINRES) nj=(nj-1)/2
  END SUBROUTINE coarsergrid

  SUBROUTINE CalcPotential_multigrid(this,Mesh,Physics,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Poisson_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,jcycle
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,pvar
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    !copy rho to local grid
    this%grid(1)%rho(0:this%grid(1)%ni+1,0:this%grid(1)%nj+1) =&
         pvar(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1,Physics%DENSITY)

!only Test: ***********
!multigrid step only if error...
! IF (MAXVAL(ABS(resid(this,this%grid(1)%u,this%grid(1)%rho,1)) .LE. this%MAXRESIDNORM) THEN
!    this%safedmultigrid = this%safedmultigrid + 1
! ELSE
    IF (Initialized(this%multipole)) &
         CALL CalculatePotential(this%multipole,Physics, &
         this%grid(1)%rho,this%grid(1)%u)

    DO j=1,this%NGRID-1
       CALL restrict(this%grid(j)%rho,&
                     this%grid(j)%ni,this%grid(j)%nj,this%grid(j+1)%ni,this%grid(j+1)%nj,&
                     this%grid(j+1)%rho)
       !only boundary values are relevant
       CALL restrict(this%grid(j)%u,&
                     this%grid(j)%ni,this%grid(j)%nj,this%grid(j+1)%ni,this%grid(j+1)%nj,&
                     this%grid(j+1)%u)
    END DO

     DO j=this%NGRID,1,-1
       DO jcycle=1,this%NMAXCYCLE
          CALL multigridREC(this,j,this%grid(j)%u,this%grid(j)%rho)
          IF (MAXVAL(ABS(resid(this,this%grid(j)%u,this%grid(j)%rho,j))&
                         -this%MAXRESIDNORM*ABS(this%grid(j)%u(:,:))) .LE. 0.0  ) EXIT
          IF (jcycle .EQ. this%NMAXCYCLE) &
             CALL Error(this,"CalcPotential_multigrid", "no convergence! (max(resid) > MAXRESIDNORM)")
       END DO
       IF (j .GT. 1)&
             CALL prolong(this%grid(j)%u, this%grid(j)%ni,this%grid(j)%nj,&
                          this%grid(j-1)%ni,this%grid(j-1)%nj,this%grid(j-1)%u)
      END DO
! END IF
!***********************
    !copy u to phi
    this%phi(:,:) = 4.0*PI*Physics%constants%GN*this%grid(1)%u(:,:)
     
  CONTAINS
    RECURSIVE SUBROUTINE multigridREC(this,j,u,rhs)
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      TYPE(Poisson_TYP)                :: this
      INTEGER                          :: j
      REAL, DIMENSION(0:this%grid(j)%ni+1,0:this%grid(j)%nj+1) :: u,rhs
      !------------------------------------------------------------------------!
      INTENT(IN)                       :: rhs,j
      INTENT(INOUT)                    :: this,u
      !------------------------------------------------------------------------!
      INTEGER                          :: jpost,jpre
!CAUTION: slow!
      REAL, DIMENSION(0:this%grid(MIN(j+1,this%ngrid))%ni+1, &
                      0:this%grid(MIN(j+1,this%ngrid))%nj+1) :: res,v
      !------------------------------------------------------------------------!
      DO jpre=1,this%NPRE
         SELECT CASE(this%RELAXTYPE)
         CASE(RED_BLACK_GAUSS_SEIDEL)
!CDIR IEXPAND
             CALL solverbgs(this,u,rhs,j) !red-black Gauss-Seidel relaxation
         CASE(BLOCK_GAUSS_SEIDEL)
!CDIR IEXPAND
            CALL solveazgs(this,u,rhs,j) !zebra line Gauss-Seidel relaxation
         CASE(GAUSS_SEIDEL)
!CDIR IEXPAND
            CALL solvegs(this,u,rhs,j) !Gauss-Seidel relaxation
         END SELECT
      END DO

      IF (j .LT. this%ngrid) THEN
         CALL restrict(resid(this,u,rhs,j),this%grid(j)%ni,this%grid(j)%nj,&
                             this%grid(j+1)%ni,this%grid(j+1)%nj,res(:,:))
         v(:,:)=0.0
         CALL multigridREC(this,j+1,v,res)
 this%grid(j)%tmp(:,:)=0.0
         CALL prolong(v(:,:),this%grid(j+1)%ni,this%grid(j+1)%nj,&
                         this%grid(j)%ni,this%grid(j)%nj,this%grid(j)%tmp(:,:))
         u(:,:) = u(:,:) + this%grid(j)%tmp(:,:)
      END IF

      DO jpost=1,this%NPOST
        SELECT CASE(this%RELAXTYPE)
         CASE(RED_BLACK_GAUSS_SEIDEL)
!CDIR IEXPAND
             CALL solverbgs(this,u,rhs,j) !red-black Gauss-Seidel relaxation
         CASE(BLOCK_GAUSS_SEIDEL)
!CDIR IEXPAND
            CALL solveazgs(this,u,rhs,j) !zebra line Gauss-Seidel relaxation
         CASE(GAUSS_SEIDEL)
!CDIR IEXPAND
            CALL solvegs(this,u,rhs,j) !Gauss-Seidel relaxation
         END SELECT
      END DO

    END SUBROUTINE multigridREC
  END SUBROUTINE CalcPotential_multigrid

PURE SUBROUTINE prolong(uc,nic,njc,nif,njf,uf)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER,INTENT(IN)                 :: nic,njc,nif,njf
    REAL, DIMENSION(0:nic+1,0:njc+1)   :: uc
    REAL, DIMENSION(0:nif+1,0:njf+1)   :: uf
    !------------------------------------------------------------------------!
    INTENT(IN)                         :: uc
    INTENT(INOUT)                      :: uf
    !------------------------------------------------------------------------!

    IF (nif .GT. nic .AND. njf .GT. njc) THEN
      !interpol in both directions
!CDIR NODEP
      uf(2:nif-1:2,2:njf-1:2)= uc(1:nic,1:njc)
!CDIR NODEP
      uf(2:nif-1:2,1:njf:2)  = 0.5*(uc(1:nic,0:njc)+uc(1:nic,1:njc+1))
!CDIR NODEP
      uf(1:nif:2,2:njf-1:2)  = 0.5*(uc(0:nic,1:njc)+uc(1:nic+1,1:njc))
!CDIR NODEP
      uf(1:nif:2,1:njf:2)    = 0.25*(uc(0:nic,0:njc)+uc(1:nic+1,0:njc)&
                                    +uc(0:nic,1:njc+1)+uc(1:nic+1,1:njc+1))
    ELSE IF (nif .GT. nic) THEN
      !interpol in i-direction
!CDIR NODEP
      uf(2:nif-1:2,1:njf)    = uc(1:nic,1:njc) 
!CDIR NODEP
      uf(1:nif:2,1:njf)      = 0.5*(uc(0:nic,1:njc)+uc(1:nic+1,1:njc))
    ELSE
      !interpol in j-direction
!CDIR NODEP
       uf(1:nif,2:njf-1:2)   = uc(1:nic,1:njc) 
!CDIR NODEP
       uf(1:nif,1:njf:2)     = 0.5*(uc(1:nic,0:njc)+uc(1:nic,1:njc+1))
    END IF

  END SUBROUTINE prolong

PURE SUBROUTINE restrict(uf,nif,njf,nic,njc,uc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER,INTENT(IN)                 :: nic,njc,nif,njf
    REAL, DIMENSION(0:nif+1,0:njf+1)   :: uf
    REAL, DIMENSION(0:nic+1,0:njc+1)   :: uc
    !------------------------------------------------------------------------!
    INTENT(IN)             :: uf
    INTENT(INOUT)          :: uc
    !------------------------------------------------------------------------!
    IF (nif .GT. nic .AND. njf .GT. njc) THEN
      !restrict in both directions
!CDIR NODEP
      uc(1:nic,1:njc) = restrict2D(uf(2:nif-1:2,2:njf-1:2),&
                                   uf(1:nif-2:2,2:njf-1:2),&
                                   uf(3:nif:2,2:njf-1:2),&
                                   uf(2:nif-1:2,1:njf-2:2),&
                                   uf(2:nif-1:2,3:njf:2),&
                                   uf(1:nif-2:2,1:njf-2:2),&
                                   uf(3:nif:2,1:njf-2:2),&
                                   uf(3:nif:2,3:njf:2),&
                                   uf(1:nif-2:2,3:njf:2))

      !ghost cells
!CDIR NODEP
      uc(0,1:njc)     = restrict1D(uf(0,2:njf-1:2),&
                                   uf(0,1:njf-2:2),&
                                   uf(0,3:njf:2))

!CDIR NODEP
      uc(nic+1,1:njc) =  restrict1D(uf(nif+1,2:njf-1:2),&
                                    uf(nif+1,1:njf-2:2),&
                                    uf(nif+1,3:njf:2))

!CDIR NODEP
      uc(1:nic,0)     =  restrict1D(uf(2:nif-1:2,0),&
                                    uf(1:nif-2:2,0),&
                                    uf(3:nif:2,0))

!CDIR NODEP
      uc(1:nic,njc+1) =  restrict1D(uf(2:nif-1:2,njf+1),&
                                    uf(1:nif-2:2,njf+1),&
                                    uf(3:nif:2,njf+1))
      !corner ghost cells
      uc(0    ,0)     = uf(0,0)
      uc(nic+1,0)     = uf(nif+1,0)
      uc(0    ,njc+1) = uf(0    ,njf+1)
      uc(nic+1,njc+1) = uf(nif+1,njf+1)
    ELSE IF (nif .GT. nic) THEN
      !restrict in i-direction
      !+boundary
      ! => njf == njc
!CDIR NODEP
      uc(1:nic,0:njf+1) = restrict1D(uf(2:nif-1:2,0:njc+1),&
                                     uf(1:nif-2:2,0:njc+1),&
                                     uf(3:nif  :2,0:njc+1))

      !copy boundary
      uc(0    ,0:njc+1)     = uf(0,0:njf+1)
      uc(nic+1,0:njc+1)     = uf(nif+1,0:njf+1)
    ELSE
      !restrict in j-direction
      !+boundary
      ! => nif == nic
!CDIR NODEP
      uc(0:nic+1,1:njc) = restrict1D(uf(0:nif+1,2:njf-1:2),&
                                     uf(0:nif+1,1:njf-2:2),&
                                     uf(0:nif+1,3:njf:2))

     
      !copy boundary
      uc(0:nic+1,0)     = uf(0:nif+1,0)
      uc(0:nic+1,njc+1) = uf(0:nif+1,njf+1)
    END IF

 END SUBROUTINE restrict

ELEMENTAL FUNCTION restrict2D(u,uE,uW,uS,uN,uSE,uSW,uNE,uNW) RESULT(restrict)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: u,uE,uW,uS,uN,uSE,uSW,uNE,uNW
    REAL              :: restrict
    !------------------------------------------------------------------------!
    restrict = (4.0*u+2.0*(uE+uW+uS+uN)+uSE+uSW+uNE+uNW )*0.0625
 END FUNCTION restrict2D

 ELEMENTAL FUNCTION restrict1D(u,uE,uW) RESULT(restrict)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: u,uE,uW
    REAL              :: restrict
    !------------------------------------------------------------------------!
       restrict = (2.0*u+uE+uW)*0.25
 END FUNCTION restrict1D

! red-black Gauss-Seidel relaxation
PURE SUBROUTINE solverbgs(this,u,rhs,jgrid)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Poisson_TYP), INTENT(INOUT) :: this
   INTEGER, INTENT(IN)              :: jgrid
   REAL, DIMENSION(0:this%grid(jgrid)%ni+1,0:this%grid(jgrid)%nj+1)&
                                    :: u,rhs
   !------------------------------------------------------------------------!
   INTEGER                          :: ni,nj,                      i
   !------------------------------------------------------------------------!
   INTENT(INOUT)                    :: u
   INTENT(IN)                       :: rhs
   !------------------------------------------------------------------------!
!test!
!  this%relaxcount(jgrid) = this%relaxcount(jgrid)+1

   CALL bndryrelax(this,u,jgrid)

   ni = this%grid(jgrid)%ni
   nj = this%grid(jgrid)%nj
 
 !red/black relaxation
    !"black" fields: part 1
!CDIR NODEP
    u(1:ni:2,1:nj:2) = relaxe(u(2:ni+1:2,1:nj  :2),&
                              u(0:ni-1:2,1:nj  :2),&
                              u(1:ni  :2,2:nj+1:2),&
                              u(1:ni  :2,0:nj-1:2),&
                            rhs(1:ni  :2,1:nj  :2),&
             this%grid(jgrid)%a(1:ni  :2,1:nj  :2),&
            this%grid(jgrid)%da(1:ni  :2,1:nj  :2),&
             this%grid(jgrid)%b(1:ni  :2,1:nj  :2),&
            this%grid(jgrid)%db(1:ni  :2,1:nj  :2),&
             this%grid(jgrid)%c(1:ni  :2,1:nj  :2),&
             this%grid(jgrid)%d(1:ni  :2,1:nj  :2),&
             this%grid(jgrid)%invhi2,this%grid(jgrid)%invhj2)  

    !"black" fields: part 2
!CDIR NODEP
    u(2:ni:2,2:nj:2) =     relaxe(u(3:ni+1:2,2:nj:2),&
                                  u(1:ni-1:2,2:nj:2),&
                                  u(2:ni  :2,3:nj+1  :2),&
                                  u(2:ni  :2,1:nj-1:2),&
                                rhs(2:ni  :2,2:nj  :2),&
                 this%grid(jgrid)%a(2:ni  :2,2:nj  :2),&
                this%grid(jgrid)%da(2:ni  :2,2:nj  :2),&
                 this%grid(jgrid)%b(2:ni  :2,2:nj  :2),&
                this%grid(jgrid)%db(2:ni  :2,2:nj  :2),&
                 this%grid(jgrid)%c(2:ni  :2,2:nj  :2),&
                 this%grid(jgrid)%d(2:ni  :2,2:nj  :2),&
                 this%grid(jgrid)%invhi2,this%grid(jgrid)%invhj2)  

    !"red" fields: part 1
!CDIR NODEP
    u(2:ni:2,1:nj:2) =   relaxe(u(3:ni+1:2,1:nj  :2),&
                                u(1:ni-1:2,1:nj  :2),&
                                u(2:ni  :2,2:nj+1:2),&
                                u(2:ni  :2,0:nj-1:2),&
                              rhs(2:ni  :2,1:nj  :2),&
               this%grid(jgrid)%a(2:ni  :2,1:nj  :2),&
              this%grid(jgrid)%da(2:ni  :2,1:nj  :2),&
               this%grid(jgrid)%b(2:ni  :2,1:nj  :2),&
              this%grid(jgrid)%db(2:ni  :2,1:nj  :2),&
               this%grid(jgrid)%c(2:ni  :2,1:nj  :2),&
               this%grid(jgrid)%d(2:ni  :2,1:nj  :2),&
               this%grid(jgrid)%invhi2,this%grid(jgrid)%invhj2)

    !"red" fields: part 2
!CDIR NODEP
    u(1:ni:2,2:nj:2) =   relaxe(u(2:ni+1:2,2:nj  :2),&
                                u(0:ni-1:2,2:nj  :2),&
                                u(1:ni  :2,3:nj+1:2),&
                                u(1:ni  :2,1:nj-1:2),&
                              rhs(1:ni  :2,2:nj  :2),&
               this%grid(jgrid)%a(1:ni  :2,2:nj  :2),&
              this%grid(jgrid)%da(1:ni  :2,2:nj  :2),&
               this%grid(jgrid)%b(1:ni  :2,2:nj  :2),&
              this%grid(jgrid)%db(1:ni  :2,2:nj  :2),&
               this%grid(jgrid)%c(1:ni  :2,2:nj  :2),&
               this%grid(jgrid)%d(1:ni  :2,2:nj  :2),&
               this%grid(jgrid)%invhi2,this%grid(jgrid)%invhj2) 
   
  END SUBROUTINE solverbgs

  ! Gauss-Seidel-type alternating zebra line relaxation (aZGS)
  PURE SUBROUTINE solveazgs(this,u,rhs,jgrid)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Poisson_TYP), INTENT(INOUT) :: this
   INTEGER, INTENT(IN)              :: jgrid
   REAL, DIMENSION(0:this%grid(jgrid)%ni+1,0:this%grid(jgrid)%nj+1)&
                                    :: u,rhs
   REAL                             :: m
   !------------------------------------------------------------------------!
   INTEGER                          :: ni,nj,i,j,k
   !------------------------------------------------------------------------!
   INTENT(INOUT)                    :: u
   INTENT(IN)                       :: rhs
   !------------------------------------------------------------------------!
 !FIXME***
!  this%relaxcount(jgrid) = this%relaxcount(jgrid)+2
 !****

    CALL bndryrelax(this,u,jgrid)

    ni = this%grid(jgrid)%ni
    nj = this%grid(jgrid)%nj

    !odd x-lines 
!CDIR NODEP
    DO i=1,ni,2
      this%grid(jgrid)%tri(i,1:nj,6) = this%grid(jgrid)%c(i,1:nj)*rhs(i,1:nj)&
                             -this%grid(jgrid)%tri(i,1:nj,4)*u(i-1,1:nj)&
                             -this%grid(jgrid)%tri(i,1:nj,5)*u(i+1,1:nj)
      this%grid(jgrid)%tri(i,0,7) = 1.0
      this%grid(jgrid)%tri(i,0,8) = u(i,0)
    END DO
    
    !The first pass (setting coefficients):
    DO k = 1,nj+1
!CDIR NODEP
      DO i=1,ni,2
        this%grid(jgrid)%tri(i,k,9) = this%grid(jgrid)%tri(i,k,2)/this%grid(jgrid)%tri(i,k-1,7)
        this%grid(jgrid)%tri(i,k,7) = this%grid(jgrid)%tri(i,k,1) - this%grid(jgrid)%tri(i,k,9)*this%grid(jgrid)%tri(i,k-1,3)
        this%grid(jgrid)%tri(i,k,8) = this%grid(jgrid)%tri(i,k,6) - this%grid(jgrid)%tri(i,k,9)*this%grid(jgrid)%tri(i,k-1,8)
      END DO
    END DO
 
      !The second pass (back-substition)
    DO k = nj, 1, -1
!CDIR NODEP
      DO i=1,ni,2
        u(i,k) = (this%grid(jgrid)%tri(i,k,8) - this%grid(jgrid)%tri(i,k,3)*u(i,k+1))/this%grid(jgrid)%tri(i,k,7)
      END DO
    END DO

    !even x-lines 
!CDIR NODEP
    DO i=2,ni,2
      this%grid(jgrid)%tri(i,1:nj,6) = this%grid(jgrid)%c(i,1:nj)*rhs(i,1:nj)&
                             -this%grid(jgrid)%tri(i,1:nj,4)*u(i-1,1:nj)&
                             -this%grid(jgrid)%tri(i,1:nj,5)*u(i+1,1:nj)
      this%grid(jgrid)%tri(i,0,7) = 1.0
      this%grid(jgrid)%tri(i,0,8) = u(i,0)
    END DO
 
      !The first pass (setting coefficients):
    DO k = 1,nj+1
!CDIR NODEP
      DO i=2,ni,2
        this%grid(jgrid)%tri(i,k,9) = this%grid(jgrid)%tri(i,k,2)/this%grid(jgrid)%tri(i,k-1,7)
        this%grid(jgrid)%tri(i,k,7) = this%grid(jgrid)%tri(i,k,1) - this%grid(jgrid)%tri(i,k,9)*this%grid(jgrid)%tri(i,k-1,3)
        this%grid(jgrid)%tri(i,k,8) = this%grid(jgrid)%tri(i,k,6) - this%grid(jgrid)%tri(i,k,9)*this%grid(jgrid)%tri(i,k-1,8)
      END DO
    END DO
 
      !The second pass (back-substition)
    DO k = nj, 1, -1
!CDIR NODEP
      DO i=2,ni,2
        u(i,k) = (this%grid(jgrid)%tri(i,k,8) - this%grid(jgrid)%tri(i,k,3)*u(i,k+1))/this%grid(jgrid)%tri(i,k,7)
      END DO
    END DO
  
    !even y-lines 
!CDIR NODEP
    DO j=2,nj,2
      this%grid(jgrid)%tri(1:ni,j,6) = this%grid(jgrid)%c(1:ni,j)*rhs(1:ni,j)&
                                   -this%grid(jgrid)%tri(1:ni,j,2)*u(1:ni,j-1)&
                                   -this%grid(jgrid)%tri(1:ni,j,3)*u(1:ni,j+1)
      this%grid(jgrid)%tri(0,j,7) = 1.0
      this%grid(jgrid)%tri(0,j,8) = u(0,j)
    END DO
 
      !The first pass (setting coefficients):
    DO k = 1,ni+1
!CDIR NODEP
      DO j=2,nj,2
        this%grid(jgrid)%tri(k,j,9) = this%grid(jgrid)%tri(k,j,4)/this%grid(jgrid)%tri(k-1,j,7)
        this%grid(jgrid)%tri(k,j,7) = this%grid(jgrid)%tri(k,j,1) - this%grid(jgrid)%tri(k,j,9)*this%grid(jgrid)%tri(k-1,j,5)
        this%grid(jgrid)%tri(k,j,8) = this%grid(jgrid)%tri(k,j,6) - this%grid(jgrid)%tri(k,j,9)*this%grid(jgrid)%tri(k-1,j,8)
      END DO
    END DO
 
      !The second pass (back-substition)
    DO k = ni, 1, -1
!CDIR NODEP
      DO j=2,nj,2
        u(k,j) = (this%grid(jgrid)%tri(k,j,8) - this%grid(jgrid)%tri(k,j,5)*u(k+1,j))/this%grid(jgrid)%tri(k,j,7)
      END DO
    END DO

    !odd y-lines 
!CDIR NODEP
    DO j=1,nj,2
      this%grid(jgrid)%tri(1:ni,j,6) = this%grid(jgrid)%c(1:ni,j)*rhs(1:ni,j)&
                                   -this%grid(jgrid)%tri(1:ni,j,2)*u(1:ni,j-1)&
                                   -this%grid(jgrid)%tri(1:ni,j,3)*u(1:ni,j+1)
      this%grid(jgrid)%tri(0,j,7) = 1.0
      this%grid(jgrid)%tri(0,j,8) = u(0,j)
    END DO
 
      !The first pass (setting coefficients):
    DO k = 1,ni+1
!CDIR NODEP
      DO j=1,nj,2
        this%grid(jgrid)%tri(k,j,9) = this%grid(jgrid)%tri(k,j,4)/this%grid(jgrid)%tri(k-1,j,7)
        this%grid(jgrid)%tri(k,j,7) = this%grid(jgrid)%tri(k,j,1) - this%grid(jgrid)%tri(k,j,9)*this%grid(jgrid)%tri(k-1,j,5)
        this%grid(jgrid)%tri(k,j,8) = this%grid(jgrid)%tri(k,j,6) - this%grid(jgrid)%tri(k,j,9)*this%grid(jgrid)%tri(k-1,j,8)
      END DO
    END DO
 
      !The second pass (back-substition)
    DO k = ni, 1, -1
!CDIR NODEP
      DO j=1,nj,2
        u(k,j) = (this%grid(jgrid)%tri(k,j,8) - this%grid(jgrid)%tri(k,j,5)*u(k+1,j))/this%grid(jgrid)%tri(k,j,7)
      END DO
    END DO
  END SUBROUTINE solveazgs

  ! Gauss-Seidel relaxation
  PURE SUBROUTINE solvegs(this,u,rhs,jgrid)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Poisson_TYP), INTENT(INOUT) :: this
   INTEGER, INTENT(IN)              :: jgrid
   REAL, DIMENSION(0:this%grid(jgrid)%ni+1,0:this%grid(jgrid)%nj+1)&
                                    :: u,rhs
   !------------------------------------------------------------------------!
   INTEGER                          :: ni,nj,                      i
   !------------------------------------------------------------------------!
   INTENT(INOUT)                    :: u
   INTENT(IN)                       :: rhs
   !------------------------------------------------------------------------!
!test!
!  this%relaxcount(jgrid) = this%relaxcount(jgrid)+1

   CALL bndryrelax(this,u,jgrid)

   ni = this%grid(jgrid)%ni
   nj = this%grid(jgrid)%nj

     u(1:ni,1:nj) =    relaxe(u(2:ni+1,1:nj),&
                              u(0:ni-1,1:nj),&
                              u(1:ni,2:nj+1),&
                              u(1:ni,0:nj-1),&
                            rhs(1:ni,1:nj),&
             this%grid(jgrid)%a(1:ni,1:nj),&
            this%grid(jgrid)%da(1:ni,1:nj),&
             this%grid(jgrid)%b(1:ni,1:nj),&
            this%grid(jgrid)%db(1:ni,1:nj),&
             this%grid(jgrid)%c(1:ni,1:nj),&
             this%grid(jgrid)%d(1:ni,1:nj),&
             this%grid(jgrid)%invhi2,this%grid(jgrid)%invhj2)  
  END SUBROUTINE solvegs

  PURE SUBROUTINE bndryrelax(this,u,jgrid)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Poisson_TYP), INTENT(INOUT) :: this
   INTEGER, INTENT(IN)              :: jgrid
   REAL, DIMENSION(0:this%grid(jgrid)%ni+1,0:this%grid(jgrid)%nj+1)&
                                    :: u
   !------------------------------------------------------------------------!
   INTEGER                          :: ni,nj
   !------------------------------------------------------------------------!
   INTENT(INOUT)                    :: u
   !------------------------------------------------------------------------!
   ni = this%grid(jgrid)%ni
   nj = this%grid(jgrid)%nj

   IF (this%Boundary(WEST) .EQ. NEUMANN) THEN
      ! Neumann
!CDIR NODEP
      u(0,1:nj) = u(1,1:nj)
   ELSE IF (this%Boundary(WEST) .EQ. PERIODIC) THEN
      ! Periodic
!CDIR NODEP
      u(0,1:nj) = u(ni,1:nj)
   END IF

   IF (this%Boundary(EAST) .EQ. NEUMANN) THEN
      ! Neumann
!CDIR NODEP
      u(ni+1,1:nj) = u(ni,1:nj)
   ELSE IF (this%Boundary(EAST) .EQ. PERIODIC) THEN
      ! Periodic
!CDIR NODEP
      u(ni+1,1:nj) = u(1,1:nj)
   END IF

   IF (this%Boundary(SOUTH) .EQ. NEUMANN) THEN
      ! Neumann
!CDIR NODEP
      u(1:ni,0) = u(1:ni,1)
   ELSE IF (this%Boundary(SOUTH) .EQ. PERIODIC) THEN
      ! Periodic
!CDIR NODEP
      u(1:ni,0) = u(1:ni,nj)
   END IF

   IF (this%Boundary(NORTH) .EQ. NEUMANN) THEN
      ! Neumann
!CDIR NODEP
      u(1:ni,nj+1) = u(1:ni,nj)
   ELSE IF (this%Boundary(NORTH) .EQ. PERIODIC) THEN
      ! Periodic
!CDIR NODEP
      u(1:ni,nj+1) = u(1:ni,1)
   END IF

  END SUBROUTINE bndryrelax

  ELEMENTAL FUNCTION relaxe(uE,uW,uN,uS,rhs,a,da,b,db,c,d,invhi2,invhj2) RESULT(u)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: uE,uW,uS,uN,rhs,a,da,b,db,c,d,invhi2,invhj2
    REAL              :: u
    !------------------------------------------------------------------------!
    u = d*(0.125*(invhi2*(da*(uE-uW)+4.0*a*(uE+uW))&
                 +invhj2*(db*(uN-uS)+4.0*b*(uN+uS)))&
          -0.5*c*rhs)
  END FUNCTION relaxe


  PURE FUNCTION resid(this,u,rhs,jgrid)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Poisson_TYP), INTENT(IN)    :: this
   INTEGER, INTENT(IN)              :: jgrid
   REAL, DIMENSION(0:this%grid(jgrid)%ni+1,0:this%grid(jgrid)%nj+1), INTENT(IN) &
                                    :: u,rhs
   !------------------------------------------------------------------------!
   INTEGER                         :: ni,nj
   REAL, DIMENSION(0:this%grid(jgrid)%ni+1,0:this%grid(jgrid)%nj+1) :: resid
   !------------------------------------------------------------------------!
   ni = this%grid(jgrid)%ni
   nj = this%grid(jgrid)%nj
   !rhs - Delta u
   resid(1:ni,1:nj) = rhs(1:ni,1:nj)-&
                     this%grid(jgrid)%invc(1:ni,1:nj)*&
                      (    this%grid(jgrid)%invhi2*( &
                         + 0.25*this%grid(jgrid)%da(1:ni,1:nj) *&
                          (u(2:ni+1,1:nj)-u(0:ni-1,1:nj))&
                         + this%grid(jgrid)%a(1:ni,1:nj)*&
                           (u(2:ni+1,1:nj)+u(0:ni-1,1:nj)-2.0*u(1:ni,1:nj)))&
                         + this%grid(jgrid)%invhj2*( &
                         + 0.25*this%grid(jgrid)%db(1:ni,1:nj) *&
                           (u(1:ni,2:nj+1)-u(1:ni,0:nj-1))&
                         + this%grid(jgrid)%b(1:ni,1:nj)*&
                           (u(1:ni,2:nj+1)+u(1:ni,0:nj-1)-2.0*u(1:ni,1:nj))))
                         

   resid(0     ,0:nj+1)=0.0
   resid(  ni+1,0:nj+1)=0.0
   resid(0:ni+1,0     )=0.0
   resid(0:ni+1,  nj+1)=0.0
END FUNCTION resid



SUBROUTINE ClosePoisson_multigrid(this)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Poisson_TYP) :: this
   !------------------------------------------------------------------------!
   INTEGER           :: j
   !------------------------------------------------------------------------!
   INTENT(INOUT)     :: this
   !------------------------------------------------------------------------!
!********TEST***************
!   DO j=this%ngrid,1,-1
!     print *, j, "Relax steps", this%relaxcount(j), "numerical effort", &
!              this%relaxcount(j)*(this%grid(j)%ni+2)*(this%grid(j)%nj+2)
!   END DO
!   print *, "entire numerical effort    ", SUM(this%relaxcount(:)*(this%grid(:)%ni+2)*(this%grid(:)%nj+2))
!   print *, "effort per cell",  SUM(this%relaxcount(:)*this%grid(:)%ni*this%grid(:)%nj)&
!                                /(this%grid(1)%ni*this%grid(1)%nj)
!   print *, "safed multigrid steps " , this%safedmultigrid
               
!******END TEST**************

   IF (Initialized(this%multipole)) CALL CloseMultipole(this%multipole)
!CDIR NODEP
    DO j=1,this%NGRID
       DEALLOCATE(this%grid(j)%u)
       DEALLOCATE(this%grid(j)%rho)
       DEALLOCATE(this%grid(j)%a)
       DEALLOCATE(this%grid(j)%da)
       DEALLOCATE(this%grid(j)%b)
       DEALLOCATE(this%grid(j)%db)
       DEALLOCATE(this%grid(j)%c)
       DEALLOCATE(this%grid(j)%invc)
       DEALLOCATE(this%grid(j)%d)
       DEALLOCATE(this%grid(j)%vol)
       DEALLOCATE(this%grid(j)%bccart)
       DEALLOCATE(this%grid(j)%bhx)
       DEALLOCATE(this%grid(j)%bhy)
       IF (this%relaxtype .EQ. BLOCK_GAUSS_SEIDEL) DEALLOCATE(this%grid(j)%tri)
    END DO
    DEALLOCATE(this%grid)
  END SUBROUTINE ClosePoisson_multigrid

END MODULE poisson_multigrid
