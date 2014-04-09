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
!red/black order
  INTEGER, PARAMETER, DIMENSION(4)  :: rb_i = (/ 1,2,1,2 /)
  INTEGER, PARAMETER, DIMENSION(4)  :: rb_j = (/ 1,2,2,1 /)  
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
    INTEGER           :: err,i,j,k,ni,nj,jgrid
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
      IF (coarsergrid(this%MINRES,ni) .EQ. ni .AND. coarsergrid(this%MINRES,nj) .EQ. nj) EXIT
      ni = coarsergrid(this%MINRES,ni)
      nj = coarsergrid(this%MINRES,nj)
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
    WRITE (xres_str, '(I0)') ni
    WRITE (yres_str, '(I0)') nj
    WRITE (grid_str, '(I0)') this%NGRID
    CALL Info(this, " POISSON--> multigrid method   " // TRIM(relaxtype_str))
    CALL Info(this, "            grid levels:       " // TRIM(grid_str) )
    CALL Info(this, "            coarsest grid res. " // TRIM(xres_str) // " x " // TRIM(yres_str))
        
    IF (ni .GT. 10 .OR. nj .GT. 10) CALL Warning(this, "InitPoisson_multigrid",&
                "bad multigrid convergence: choose a better resolution e.g. m*2^n+1 with small natural number m")

!  !FIXME convergence test***************************
!  ALLOCATE(this%relaxcount(this%ngrid),STAT = err)
!  IF (err.NE.0) CALL Error(this, "InitPoisson_multigrid", "Unable to allocate memory.")
!  this%relaxcount(:) = 0
!  !*************************************************

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

       CALL INITGRID(this,Mesh,jgrid,ni,nj)

       ni = coarsergrid(this%MINRES,ni)
       nj = coarsergrid(this%MINRES,nj)
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

       !initialize the multipole expansion module
       CALL InitMultipole(this%multipole,bndrytype,this%grid(1)%bccart, &
            2*PI*this%grid(1)%vol,0,this%grid(1)%ni+1, &
            0,this%grid(1)%nj+1,iregion,bnd_region,maxmult)

    END IF
  END SUBROUTINE InitPoisson_multigrid


  SUBROUTINE InitGrid(this,Mesh,jgrid,ni,nj)
  IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Poisson_TYP)      :: this
    TYPE(Mesh_TYP)         :: Mesh
    INTEGER                :: jgrid,ni,nj
    !------------------------------------------------------------------------!
    TYPE(Grid_TYP),POINTER :: grid,pregrid
    INTEGER                :: i,j,err,sp_i, sp_j
    !------------------------------------------------------------------------!
    INTENT(IN)             :: Mesh,jgrid,ni,nj
    INTENT(INOUT)          :: this
    !------------------------------------------------------------------------!
    grid=>this%grid(jgrid)
    
    grid%u(:,:) = 0.0
    !Length incl. 1st ghost cells
    grid%hi=(Mesh%xmax-Mesh%xmin+1.0*Mesh%dx)/(ni+1)
    grid%hj=(Mesh%ymax-Mesh%ymin+1.0*Mesh%dy)/(nj+1)

    grid%ni=ni
    grid%nj=nj
    grid%invhi2 = 1.0 / grid%hi**2
    grid%invhj2 = 1.0 / grid%hj**2
    IF (jgrid == 1) THEN
!CDIR NODEP
          DO i=0,ni+1
!CDIR NODEP
            DO j=0,nj+1
              !only at first grid
              grid%bccart(i,j,:) = Mesh%bccart(i,j,:)
              grid%bhx(i,j)      = Mesh%bhx(i,j)
              grid%bhy(i,j)      = Mesh%bhy(i,j)
              grid%a(i,j)        = Mesh%bhy(i,j)*Mesh%bhz(i,j)/Mesh%bhx(i,j)
              grid%b(i,j)        = Mesh%bhx(i,j)*Mesh%bhz(i,j)/Mesh%bhy(i,j)
              grid%c(i,j)        = Mesh%bhx(i,j)*Mesh%bhy(i,j)*Mesh%bhz(i,j)
              grid%vol(i,j)      = Mesh%volume(i,j)
            END DO
          END DO
       ELSE
          pregrid=>this%grid(jgrid-1)
!           CALL restrict(pregrid%bhx,pregrid%ni,pregrid%nj,ni,nj,grid%bhx)
!           CALL restrict(pregrid%bhy,pregrid%ni,pregrid%nj,ni,nj,grid%bhy)
!           CALL restrict(pregrid%a,pregrid%ni,pregrid%nj,ni,nj,grid%a)
!           CALL restrict(pregrid%b,pregrid%ni,pregrid%nj,ni,nj,grid%b)
!           CALL restrict(pregrid%c,pregrid%ni,pregrid%nj,ni,nj,grid%c)

          sp_i = (pregrid%ni+1)/(ni+1)
          sp_j = (pregrid%nj+1)/(nj+1)
!CDIR NODEP
          DO i=0,ni+1
!CDIR NODEP
            DO j=0,nj+1
              grid%bccart(i,j,:) = pregrid%bccart(sp_i*i,sp_j*j,:) 
              grid%bhx(i,j) = pregrid%bhx(sp_i*i,sp_j*j)
              grid%bhy(i,j) = pregrid%bhy(sp_i*i,sp_j*j)
              grid%a(i,j) = pregrid%a(sp_i*i,sp_j*j)
              grid%b(i,j) = pregrid%b(sp_i*i,sp_j*j)
              grid%c(i,j) = pregrid%c(sp_i*i,sp_j*j)
            END DO
          END DO

          grid%vol(:,:) = grid%c(:,:)*grid%hi*grid%hj

          pregrid=>this%grid(1)
          !boundary volumes of boundary cells
          grid%vol(0,:)       = grid%c(0,:)    * pregrid%hi * grid%hj
          grid%vol(ni+1,:)    = grid%c(ni+1,:) * pregrid%hi * grid%hj
          grid%vol(:,0)       = grid%c(:,0)    * grid%hi * pregrid%hj
          grid%vol(:,nj+1)    = grid%c(:,nj+1) * grid%hi * pregrid%hj
          !boundary volumes of corner cells
          grid%vol(0,0)       = grid%c(0,0)       * pregrid%hi * pregrid%hj
          grid%vol(ni+1,0)    = grid%c(ni+1,0)    * pregrid%hi * pregrid%hj
          grid%vol(ni+1,nj+1) = grid%c(ni+1,nj+1) * pregrid%hi * pregrid%hj
          grid%vol(0,nj+1)    = grid%c(0,nj+1)    * pregrid%hi * pregrid%hj

       END IF

       !FIXME*** 
       grid%vol(:,:) = abs(grid%vol(:,:))
       !***
       grid%tmp(:,:)  = 0.0
       grid%invc(:,:) = 1.0/  grid%c(:,:)
       grid%d(:,:)    = 1.0/( grid%a(:,:)*grid%invhi2 +grid%b(:,:)*grid%invhj2 )

       DO i=1,ni
         DO j=1,nj
           grid%da(i,j) = grid%a(i+1,j)-grid%a(i-1,j)
           grid%db(i,j) = grid%b(i,j+1)-grid%b(i,j-1)
         END DO
       END DO

       
       IF (this%RELAXTYPE .EQ. BLOCK_GAUSS_SEIDEL) THEN
         ALLOCATE(grid%tri(0:ni+1,0:nj+1,9),&
                  STAT = err)
         IF (err.NE.0) CALL Error(this, "InitPoisson_multigrid", "Unable to allocate memory.")

         DO i=1,ni
           DO j=1,nj
             !1. - 5. matrix elements (diag,tri diag, elements); 6. rhs; 7.-9. temp variables (hb,hv,m) 
             grid%tri(i,j,1) = -2.0*(grid%invhi2*grid%a(i,j) + grid%invhj2*grid%b(i,j))
             grid%tri(i,j,2) = grid%invhj2*(grid%b(i,j)-grid%db(i,j)*0.25)
             grid%tri(i,j,3) = grid%invhj2*(grid%b(i,j)+grid%db(i,j)*0.25)     
             grid%tri(i,j,4) = grid%invhi2*(grid%a(i,j)-grid%da(i,j)*0.25)
             grid%tri(i,j,5) = grid%invhi2*(grid%a(i,j)+grid%da(i,j)*0.25)
           END DO
         END DO
       END IF

  END SUBROUTINE InitGrid

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
    REAL              :: mass,sumres
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,pvar
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    !copy rho to local grid
    this%grid(1)%rho(0:this%grid(1)%ni+1,0:this%grid(1)%nj+1) =&
         pvar(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1,Physics%DENSITY)

    mass = SUM(this%grid(1)%rho*this%grid(1)%vol)

    IF (Initialized(this%multipole)) THEN
      CALL CalculatePotential(this%multipole,Physics, &
         this%grid(1)%rho,this%grid(1)%u)
    ELSE
       !TODO: only dirichlet boundary => special cond is required 
       CALL  Error(this,"CalcPotential_multigrid","no multipole")
    END IF

    DO j=1,this%NGRID-1
! *** see below (no need of FMG)
!        CALL restrict(this%grid(j)%rho,&
!                      this%grid(j)%ni,this%grid(j)%nj,this%grid(j+1)%ni,this%grid(j+1)%nj,&
!                      this%grid(j+1)%rho)

       !only boundary values are relevant
       CALL restrict(this%grid(j)%u,&
                     this%grid(j)%ni,this%grid(j)%nj,this%grid(j+1)%ni,this%grid(j+1)%nj,&
                     this%grid(j+1)%u)
    END DO

    !***
    ! no FMG implemented reason: old u as start solution is the best approx (better then nested iteration)
    ! => do not use FMG!
    j = 1
    DO jcycle=1,this%NMAXCYCLE
      CALL multigridREC(this,j,this%grid(j)%u,this%grid(j)%rho)
      sumres = ABS(SUM(resid(this,this%grid(j)%u,this%grid(j)%rho,j)*this%grid(j)%vol))
      IF (sumres/mass .LE. this%MAXRESIDNORM) EXIT
!           IF (MAXVAL(ABS(resid(this,this%grid(j)%u,this%grid(j)%rho,j))) .lt. this%MAXRESIDNORM) EXIT
    END DO
!  print *,jcycle, sumres, MAXVAL(ABS(resid(this,this%grid(j)%u,this%grid(j)%rho,j)/this%grid(j)%rho)), &
!           sumres/mass
    IF (jcycle .GE. this%NMAXCYCLE) &
      CALL Error(this,"CalcPotential_multigrid", "no convergence! a greater NMAXCYCLE or MAXRESIDNORM could perhaps help")

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
    INTEGER                            :: nic,njc,nif,njf
    REAL, DIMENSION(0:nic+1,0:njc+1)   :: uc
    REAL, DIMENSION(0:nif+1,0:njf+1)   :: uf
    !------------------------------------------------------------------------!
    INTENT(IN)                         :: uc,nic,njc,nif,njf
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
    INTENT(IN)                         :: uf
    INTENT(INOUT)                      :: uc
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
   TYPE(Poisson_TYP)                :: this
   INTEGER                          :: jgrid
   REAL, DIMENSION(0:this%grid(jgrid)%ni+1,0:this%grid(jgrid)%nj+1)&
                                    :: u,rhs
   !------------------------------------------------------------------------!
   INTEGER                          :: i,j,k
   !------------------------------------------------------------------------!
   INTENT(INOUT)                    :: this,u
   INTENT(IN)                       :: rhs,jgrid
   !------------------------------------------------------------------------!
!******TEST*****!
!  this%relaxcount(jgrid) = this%relaxcount(jgrid)+1

   CALL bndryrelax(this,u,jgrid)

   DO k=1,4
!CDIR NODEP
     DO i=rb_i(k),this%grid(jgrid)%ni,2
!CDIR NODEP
       DO j=rb_j(k),this%grid(jgrid)%nj,2
         u(i,j) = relaxe(u(i+1,j),&
                         u(i-1,j),&
                         u(i,j+1),&
                         u(i,j-1),&
                         rhs(i,j),&
                         this%grid(jgrid)%a(i,j),&
                         this%grid(jgrid)%da(i,j),&
                         this%grid(jgrid)%b(i,j),&
                         this%grid(jgrid)%db(i,j),&
                         this%grid(jgrid)%c(i,j),&
                         this%grid(jgrid)%d(i,j),&
                         this%grid(jgrid)%invhi2,&
                         this%grid(jgrid)%invhj2)  
       END DO
     END DO
   END DO
   
  END SUBROUTINE solverbgs

  ! Gauss-Seidel-type alternating zebra line relaxation (aZGS)
  PURE SUBROUTINE solveazgs(this,u,rhs,jgrid)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Poisson_TYP)                :: this
   INTEGER                          :: jgrid
   REAL, DIMENSION(0:this%grid(jgrid)%ni+1,0:this%grid(jgrid)%nj+1)&
                                    :: u,rhs
   REAL                             :: m
   !------------------------------------------------------------------------!
   INTEGER                          :: ni,nj,i,j,k,q
   !------------------------------------------------------------------------!
   INTENT(INOUT)                    :: this,u
   INTENT(IN)                       :: rhs,jgrid
   !------------------------------------------------------------------------!
   !******TEST*****!
!  this%relaxcount(jgrid) = this%relaxcount(jgrid)+2


    CALL bndryrelax(this,u,jgrid)

    ni = this%grid(jgrid)%ni
    nj = this%grid(jgrid)%nj

    !odd/even x-lines 
  DO q=1,2
!CDIR NODEP
    DO i=q,ni,2
      this%grid(jgrid)%tri(i,1:nj,6) = this%grid(jgrid)%c(i,1:nj)*rhs(i,1:nj)&
                             -this%grid(jgrid)%tri(i,1:nj,4)*u(i-1,1:nj)&
                             -this%grid(jgrid)%tri(i,1:nj,5)*u(i+1,1:nj)
      this%grid(jgrid)%tri(i,0,7) = 1.0
      this%grid(jgrid)%tri(i,0,8) = u(i,0)
    END DO
    
    !The first pass (setting coefficients):
    DO k = 1,nj+1
!CDIR NODEP
      DO i=q,ni,2
        this%grid(jgrid)%tri(i,k,9) = this%grid(jgrid)%tri(i,k,2)/this%grid(jgrid)%tri(i,k-1,7)
        this%grid(jgrid)%tri(i,k,7) = this%grid(jgrid)%tri(i,k,1) - this%grid(jgrid)%tri(i,k,9)*this%grid(jgrid)%tri(i,k-1,3)
        this%grid(jgrid)%tri(i,k,8) = this%grid(jgrid)%tri(i,k,6) - this%grid(jgrid)%tri(i,k,9)*this%grid(jgrid)%tri(i,k-1,8)
      END DO
    END DO
 
      !The second pass (back-substition)
    DO k = nj, 1, -1
!CDIR NODEP
      DO i=q,ni,2
        u(i,k) = (this%grid(jgrid)%tri(i,k,8) - this%grid(jgrid)%tri(i,k,3)*u(i,k+1))/this%grid(jgrid)%tri(i,k,7)
      END DO
    END DO
  END DO
  
    !even/odd y-lines 
  DO q=2,1,-1
!CDIR NODEP
    DO j=q,nj,2
      this%grid(jgrid)%tri(1:ni,j,6) = this%grid(jgrid)%c(1:ni,j)*rhs(1:ni,j)&
                                   -this%grid(jgrid)%tri(1:ni,j,2)*u(1:ni,j-1)&
                                   -this%grid(jgrid)%tri(1:ni,j,3)*u(1:ni,j+1)
      this%grid(jgrid)%tri(0,j,7) = 1.0
      this%grid(jgrid)%tri(0,j,8) = u(0,j)
    END DO
 
      !The first pass (setting coefficients):
    DO k = 1,ni+1
!CDIR NODEP
      DO j=q,nj,2
        this%grid(jgrid)%tri(k,j,9) = this%grid(jgrid)%tri(k,j,4)/this%grid(jgrid)%tri(k-1,j,7)
        this%grid(jgrid)%tri(k,j,7) = this%grid(jgrid)%tri(k,j,1) - this%grid(jgrid)%tri(k,j,9)*this%grid(jgrid)%tri(k-1,j,5)
        this%grid(jgrid)%tri(k,j,8) = this%grid(jgrid)%tri(k,j,6) - this%grid(jgrid)%tri(k,j,9)*this%grid(jgrid)%tri(k-1,j,8)
      END DO
    END DO
 
      !The second pass (back-substition)
    DO k = ni, 1, -1
!CDIR NODEP
      DO j=q,nj,2
        u(k,j) = (this%grid(jgrid)%tri(k,j,8) - this%grid(jgrid)%tri(k,j,5)*u(k+1,j))/this%grid(jgrid)%tri(k,j,7)
      END DO
    END DO
  END DO

  END SUBROUTINE solveazgs

  ! Gauss-Seidel relaxation
  PURE SUBROUTINE solvegs(this,u,rhs,jgrid)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Poisson_TYP)                :: this
   INTEGER                          :: jgrid
   REAL, DIMENSION(0:this%grid(jgrid)%ni+1,0:this%grid(jgrid)%nj+1)&
                                    :: u,rhs
   !------------------------------------------------------------------------!
   INTEGER                          :: i,j
   !------------------------------------------------------------------------!
   INTENT(INOUT)                    :: this,u
   INTENT(IN)                       :: rhs,jgrid
   !------------------------------------------------------------------------!
!******TEST*****!
!  this%relaxcount(jgrid) = this%relaxcount(jgrid)+1


   CALL bndryrelax(this,u,jgrid)
!CDIR NODEP
   DO i=1,this%grid(jgrid)%ni
!CDIR NODEP
     DO j=1,this%grid(jgrid)%nj
        u(i,j) =  relaxe(u(i+1,j),&
                        u(i-1,j),&
                        u(i,j+1),&
                        u(i,j-1),&
                        rhs(i,j),&
                        this%grid(jgrid)%a(i,j),&
                        this%grid(jgrid)%da(i,j),&
                        this%grid(jgrid)%b(i,j),&
                        this%grid(jgrid)%db(i,j),&
                        this%grid(jgrid)%c(i,j),&
                        this%grid(jgrid)%d(i,j),&
                        this%grid(jgrid)%invhi2,&
                        this%grid(jgrid)%invhj2)  
     END DO
   END DO
  END SUBROUTINE solvegs

  PURE SUBROUTINE bndryrelax(this,u,jgrid)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Poisson_TYP)                :: this
   INTEGER                          :: jgrid
   REAL, DIMENSION(0:this%grid(jgrid)%ni+1,0:this%grid(jgrid)%nj+1)&
                                    :: u
   !------------------------------------------------------------------------!
   INTEGER                          :: ni,nj
   !------------------------------------------------------------------------!
   INTENT(IN)                       :: jgrid
   INTENT(INOUT)                    :: this,u
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
   TYPE(Poisson_TYP)             :: this
   INTEGER                       :: jgrid
   REAL, DIMENSION(0:this%grid(jgrid)%ni+1,0:this%grid(jgrid)%nj+1) &
                                 :: u,rhs
   !------------------------------------------------------------------------!
   INTEGER                       :: i,j,ni,nj
   REAL                          :: du_i,ddu_i,du_j,ddu_j
   REAL, DIMENSION(0:this%grid(jgrid)%ni+1,0:this%grid(jgrid)%nj+1) :: resid
   !------------------------------------------------------------------------!
   INTENT(IN)                    :: this,jgrid,u,rhs
   !------------------------------------------------------------------------!
   ni = this%grid(jgrid)%ni
   nj = this%grid(jgrid)%nj
   !rhs - Delta u
!CDIR NODEP
   DO j=1,nj
!CDIR NODEP
     DO i=1,ni
        du_i  = delta(u(i+1,j),u(i-1,j))
        ddu_i = ddelta(u(i+1,j),u(i,j),u(i-1,j))
        du_j  = delta(u(i,j+1),u(i,j-1))
        ddu_j = ddelta(u(i,j+1),u(i,j),u(i,j-1))

        resid(i,j) = rhs(i,j)-&
                     this%grid(jgrid)%invc(i,j)*&
                      (    this%grid(jgrid)%invhi2*( &
                         + 0.25*this%grid(jgrid)%da(i,j) * du_i &
                         + this%grid(jgrid)%a(i,j) * ddu_i)&
                         + this%grid(jgrid)%invhj2*( &
                         + 0.25*this%grid(jgrid)%db(i,j) * du_j &
                         + this%grid(jgrid)%b(i,j) * ddu_j) )
     END DO
   END DO
 
   resid(0     ,0:nj+1)=0.0
   resid(  ni+1,0:nj+1)=0.0
   resid(0:ni+1,0     )=0.0
   resid(0:ni+1,  nj+1)=0.0
END FUNCTION resid

  ELEMENTAL FUNCTION delta(u1,u2) RESULT(du)
  IMPLICIT NONE
  !------------------------------------------------------------------------!
   REAL, INTENT(IN)           :: u1,u2
  !------------------------------------------------------------------------!
   REAL                       :: du
  !------------------------------------------------------------------------!
    du = u1-u2
  END FUNCTION delta

  ELEMENTAL FUNCTION ddelta(u1,u2,u3) RESULT(du)
  IMPLICIT NONE
  !------------------------------------------------------------------------!
   REAL, INTENT(IN)           :: u1,u2,u3
  !------------------------------------------------------------------------!
   REAL                       :: du
  !------------------------------------------------------------------------!
    du = u1+u3-2.0*u2
  END FUNCTION ddelta

  PURE FUNCTION coarsergrid(nmin,nold) RESULT(nnew)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER, INTENT(IN)     :: nmin
    INTEGER, INTENT(IN)     :: nold
    INTEGER                 :: nnew
    !------------------------------------------------------------------------!
    !coarsing in every step if nold is odd
    IF (mod(nold,2) .EQ. 1 .AND. (nold-1)/2 .GE. nmin) THEN
       nnew=(nold-1)/2
    ELSE 
       nnew = nold
    END IF
  END FUNCTION coarsergrid

SUBROUTINE ClosePoisson_multigrid(this)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Poisson_TYP) :: this
   !------------------------------------------------------------------------!
   INTEGER           :: j
   !------------------------------------------------------------------------!
   INTENT(INOUT)     :: this
   !------------------------------------------------------------------------!
!   !********TEST***************
!   DO j=this%ngrid,1,-1
!     print *, j, "Relax steps", this%relaxcount(j), "numerical effort", &
!              this%relaxcount(j)*(this%grid(j)%ni+2)*(this%grid(j)%nj+2)
!   END DO
!   print *, "entire numerical effort    ", SUM(this%relaxcount(:)*(this%grid(:)%ni+2)*(this%grid(:)%nj+2))
!   print *, "effort per cell",  SUM(this%relaxcount(:)*this%grid(:)%ni*this%grid(:)%nj)&
!                                /(this%grid(1)%ni*this%grid(1)%nj)
!                 
!   !******END TEST**************

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
