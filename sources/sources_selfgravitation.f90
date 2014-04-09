!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_selfgravitation.f90                                       #
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
! source terms module for selfgravitation
!----------------------------------------------------------------------------!
MODULE sources_selfgravitation
  USE sources_pointmass
  USE sources_common, ONLY : Grid_TYP
  USE mesh_common
  USE sources_boundary
  USE physics_generic
  USE boundary_generic
 IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: source_name = "selfgravitation"
  INTEGER, PARAMETER :: NPRE = 3, NPOST = 2, NMAXCYCLE = 1000000
  REAL, PARAMETER    :: NUMNULL = 1.0E-8
 !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! methods
       InitSources_selfgravitation, &
       CalcTimestep_selfgravitation, &
       CalcAccel, &
       ExternalSources_selfgravitation, &
       CloseSources_selfgravitation
  !--------------------------------------------------------------------------!

CONTAINS


  SUBROUTINE InitSources_selfgravitation(this,Mesh,Physics,Boundary,stype,&
                 maxmult,maxresidnorm,maxagmnorm,bndrytype,MGminlevel)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Boundary_TYP), DIMENSION(4) :: Boundary
    INTEGER           :: stype
    !------------------------------------------------------------------------!
    REAL              :: maxresidnorm,maxagmnorm
    INTEGER           :: maxmult,bndrytype,MGminlevel
    !------------------------------------------------------------------------!
    INTEGER           :: err,i,j,k,ni,nj,nig,njg,ngrid,si,sj
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Boundary,stype,maxmult,&
                         maxresidnorm,maxagmnorm,bndrytype,MGminlevel
    !------------------------------------------------------------------------!
    CALL InitSources(this,stype,source_name)
#ifdef PARALLEL
    CALL Error(this, "Selfgravitation", "No parallel version available")
#endif  
    this%MAXRESIDNORM = maxresidnorm
    this%MAXAGMNORM = maxagmnorm

    IF (MGminlevel >= 1) THEN
       this%MGminlevel = MGminlevel
    ELSE
    ! coarsest grid
       this%MGminlevel = 1
    END IF
 
    ALLOCATE(this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
         STAT = err)
    IF (err.NE.0) CALL Error(this, "InitSources_selfgravitation", "Unable to allocate memory.")
!CDIR NODEP
    ! set gravitational acceleration
    this%accel(:,:,:) = 0.

    ni = Mesh%INUM+2
    nj = Mesh%JNUM+2
    
    nig=nint(log(ni-1.0)/log(2.0))
    IF (ni /= 2**nig+1) CALL Error(this, "InitSources_selfgravitation", "INUM not power of 2")
    njg=nint(log(nj-1.0)/log(2.0))
    IF (nj /= 2**njg+1) CALL Error(this, "InitSources_selfgravitation", "JNUM not power of 2")

    this%ngrid=max(nig,njg)
    !to avoid problems
    IF (this%MGminlevel > this%ngrid) this%MGminlevel = this%ngrid

    ALLOCATE(this%grid(this%ngrid))
    DO ngrid=this%ngrid, this%MGminlevel, -1
       ALLOCATE(this%grid(ngrid)%u(ni,nj), &
              this%grid(ngrid)%rho(ni,nj), &
              this%grid(ngrid)%a(ni,nj), &
              this%grid(ngrid)%da(ni,nj), &
              this%grid(ngrid)%b(ni,nj), &
              this%grid(ngrid)%db(ni,nj), &
              this%grid(ngrid)%c(ni,nj), &
              this%grid(ngrid)%invc(ni,nj), &
              this%grid(ngrid)%d(ni,nj), &
              this%grid(ngrid)%vol(ni,nj), &
              this%grid(ngrid)%bccart(ni,nj,2), &
              this%grid(ngrid)%bhx(ni,nj), &
              this%grid(ngrid)%bhy(ni,nj), &
         STAT = err)
       IF (err.NE.0) CALL Error(this, "InitSources_selfgravitation", "Unable to allocate memory.")

       si = (Mesh%INUM+2-1)/(ni-1)
       sj = (Mesh%JNUM+2-1)/(nj-1)

       this%grid(ngrid)%hi=(Mesh%xmax-Mesh%xmin)/(ni-2)
       this%grid(ngrid)%hj=(Mesh%ymax-Mesh%ymin)/(nj-2)
       this%grid(ngrid)%ni=ni
       this%grid(ngrid)%nj=nj
       this%grid(ngrid)%invhi2 = 1.0 / this%grid(ngrid)%hi**2
       this%grid(ngrid)%invhj2 = 1.0 / this%grid(ngrid)%hj**2

       this%grid(ngrid)%bccart(1:ni,1:nj,:)=&
                     Mesh%bccart(Mesh%IMIN-1:Mesh%IMAX+1:si,Mesh%JMIN-1:Mesh%JMAX+1:sj,:)

       IF (ngrid == this%ngrid) THEN
          this%grid(ngrid)%bhx(1:ni,1:nj)=&
                     Mesh%bhx(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)
          this%grid(ngrid)%bhy(1:ni,1:nj)=&
                     Mesh%bhy(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)
          this%grid(ngrid)%a(1:ni,1:nj)=&
                     Mesh%bhy(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)&
                    *Mesh%bhz(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)&
                    /Mesh%bhx(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)
          this%grid(ngrid)%b(1:ni,1:nj)=&
                     Mesh%bhx(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)&
                    *Mesh%bhz(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)&
                    /Mesh%bhy(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)
          this%grid(ngrid)%c(1:ni,1:nj)=&
                     Mesh%bhx(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)&
                    *Mesh%bhy(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)&
                    *Mesh%bhz(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1)
       ELSE
          this%grid(ngrid)%bhx = restrict(this%grid(ngrid+1)%bhx,ni,nj)
          this%grid(ngrid)%bhy = restrict(this%grid(ngrid+1)%bhy,ni,nj)
          this%grid(ngrid)%a = restrict(this%grid(ngrid+1)%a,ni,nj)
          this%grid(ngrid)%b = restrict(this%grid(ngrid+1)%b,ni,nj)
          this%grid(ngrid)%c = restrict(this%grid(ngrid+1)%c,ni,nj)
       END IF

       this%grid(ngrid)%invc(:,:) = 1.0/this%grid(ngrid)%c
       this%grid(ngrid)%d(:,:) = &
                     1.0/( this%grid(ngrid)%a(:,:)*this%grid(ngrid)%invhi2 &
                     + this%grid(ngrid)%b(:,:)*this%grid(ngrid)%invhj2 )
       this%grid(ngrid)%da(2:ni-1,2:nj-1) = &
                     this%grid(ngrid)%a(3:ni,2:nj-1)-this%grid(ngrid)%a(1:ni-2,2:nj-1)
       this%grid(ngrid)%db(2:ni-1,2:nj-1) = &
                     this%grid(ngrid)%b(2:ni-1,3:nj)-this%grid(ngrid)%b(2:ni-1,1:nj-2)
       this%grid(ngrid)%vol(:,:) = abs(this%grid(ngrid)%c(:,:)&
                                * this%grid(ngrid)%hi *this%grid(ngrid)%hj)

       !ghost cells are not resized => this%ngrid hi,hj!!!!!
       this%grid(ngrid)%vol(1:ni,1)  = abs(this%grid(ngrid)%c(1:ni,1)&
                                * this%grid(this%ngrid)%hi *this%grid(this%ngrid)%hj)
       this%grid(ngrid)%vol(1:ni,nj) = abs(this%grid(ngrid)%c(1:ni,nj)&
                                * this%grid(this%ngrid)%hi *this%grid(this%ngrid)%hj)
       this%grid(ngrid)%vol(1,1:nj)  = abs(this%grid(ngrid)%c(1,1:nj)&
                                * this%grid(this%ngrid)%hi *this%grid(this%ngrid)%hj)
       this%grid(ngrid)%vol(ni,1:nj) = abs(this%grid(ngrid)%c(ni,1:nj)&
                                * this%grid(this%ngrid)%hi *this%grid(this%ngrid)%hj)

!two possible options....
!       IF (ni > 2**this%MGminlevel+1) ni=ni/2+1
!       IF (nj > 2**this%MGminlevel+1) nj=nj/2+1
        IF (ni == nj) THEN
           ni=ni/2+1
           nj=nj/2+1
        ELSE IF (ni > nj) THEN
           ni=ni/2+1
        ELSE
           nj=nj/2+1
        END IF
    END DO
    CALL InitSources_boundary(this,Mesh,Boundary,maxmult,bndrytype)

  END SUBROUTINE InitSources_selfgravitation


 SUBROUTINE CalcTimestep_selfgravitation(this,Mesh,Physics,pvar,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    REAL              :: dt,invdt
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: dt
    !------------------------------------------------------------------------!
      invdt = sqrt(MAX(MAXVAL(abs(this%accel(:,:,1)/Mesh%dlx(:,:))),&
                       MAXVAL(abs(this%accel(:,:,2)/Mesh%dly(:,:)))))

      if (invdt > 20.0*Tiny(1.0)) then
         dt = 0.5/invdt
      else
         dt = Huge(1.0)
      endif
   END SUBROUTINE CalcTimestep_selfgravitation

SUBROUTINE CalcAccel(this,Mesh,Physics,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
     REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    TYPE(Grid_TYP), POINTER :: pgrid
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,off
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,pvar
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    off = 1
    pgrid=>this%grid(this%ngrid)
    pgrid%rho(1:Mesh%INUM+2,1:Mesh%JNUM+2) =&
         pvar(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1,Physics%Density)

!Attention: PERIODIC CASE!!!! this%rho = this%rho - mean(this%rho)
    CALL multigrid(this)

    FORALL (i = Mesh%IMIN:Mesh%IMAX, j = Mesh%JMIN:Mesh%JMAX)
       ! g(x) = - grad(phi(x)) 
       this%accel(i,j,2) = 4.0*PI*Physics%constants%GN*(pgrid%u(i+off,j-1+off)&
                                                       -pgrid%u(i+off,j+1+off))/ &
                           (2.0 * Mesh%dy * Mesh%bhy(i,j)) 
       this%accel(i,j,1) = 4.0*PI*Physics%constants%GN*(pgrid%u(i-1+off,j+off)&
                                                       -pgrid%u(i+1+off,j+off))/ &
                           (2.0 * Mesh%dx * Mesh%bhx(i,j))
    END FORALL

    this%accel(Mesh%IGMIN:Mesh%IMIN-1,:,:) = 0.0
    this%accel(Mesh%IMAX+1:Mesh%IGMAX,:,:) = 0.0
    this%accel(:,Mesh%JGMIN:Mesh%JMIN-1,:) = 0.0
    this%accel(:,Mesh%JMAX+1:Mesh%JGMAX,:) = 0.0
  END SUBROUTINE CalcAccel


  SUBROUTINE ExternalSources_selfgravitation(this,Mesh,Physics,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: cvar,pvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,cvar,pvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    CALL CalcAccel(this,Mesh,Physics,pvar)
    ! gravitational source terms
    CALL ExternalSources(Physics,Mesh,this%accel,pvar,cvar,sterm)
  END SUBROUTINE ExternalSources_selfgravitation

SUBROUTINE multigrid(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    !------------------------------------------------------------------------!
    INTEGER           :: j,jcycle,i
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    !ONLY a SERIAL version
    CALL CalcBoundary(this)
    DO j=this%ngrid-1,this%MGminlevel, -1  
       this%grid(j)%rho=restrict(this%grid(j+1)%rho,&
                              this%grid(j)%ni,this%grid(j)%nj)
       !only boundary values are relevant
       this%grid(j)%u=restrict(this%grid(j+1)%u,&
                              this%grid(j)%ni,this%grid(j)%nj)
    END DO
    DO j=this%MGminlevel,this%ngrid
       DO jcycle=1,NMAXCYCLE
          CALL multigridREC(this,j,this%grid(j)%u,this%grid(j)%rho)
          IF (MAXVAL(abs(resid(this,this%grid(j)%u,this%grid(j)%rho,j)&
              /this%grid(j)%u),abs(this%grid(j)%u(:,:)) > NUMNULL) < this%MAXRESIDNORM ) EXIT
          IF (jcycle == NMAXCYCLE)&
             CALL Error(this,"multigrid", "no convergence! (resid > maxresidnorm)")
       END DO
       IF (j < this%ngrid)&
          this%grid(j+1)%u(2:this%grid(j+1)%ni-1,2:this%grid(j+1)%nj-1)=&
              prolong2(this%grid(j)%u,this%grid(j+1)%ni,this%grid(j+1)%nj)
    END DO
 CONTAINS
 RECURSIVE SUBROUTINE multigridREC(this,j,u,rhs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP)                :: this
    INTEGER                          :: j
    REAL, DIMENSION(:,:)             :: u,rhs
    !------------------------------------------------------------------------!
    INTENT(IN)                       :: rhs,j
    INTENT(INOUT)                    :: this,u
    !------------------------------------------------------------------------!
    INTEGER                          :: jpost,jpre
    REAL, DIMENSION( this%grid(MAX(j-1,this%MGminlevel))%ni , this%grid(MAX(j-1,this%MGminlevel))%nj) :: res,v
    !------------------------------------------------------------------------!
    IF (j == this%MGminlevel) THEN
       DO jpre=1,NPRE
!CDIR IEXPAND
          CALL relax(this,u,rhs,j)
       END DO
    ELSE
       DO jpre=1,NPRE
!CDIR IEXPAND
          CALL relax(this,u,rhs,j)
       END DO
       res=restrict(resid(this,u,rhs,j),size(res,1),size(res,2))
       v=0.0
       CALL multigridREC(this,j-1,v,res)
       u(2:size(u,1)-1,2:size(u,2)-1)=u(2:size(u,1)-1,2:size(u,2)-1)+prolong2(v,size(u,1),size(u,2))
       DO jpost=1,NPOST
!CDIR IEXPAND
         CALL relax(this,u,rhs,j)
       END DO
    END IF
 END SUBROUTINE multigridREC
END SUBROUTINE multigrid

 PURE FUNCTION restrict(uf,nic,njc)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   REAL, DIMENSION(:,:)         :: uf
   INTEGER                      :: nic, njc
   !------------------------------------------------------------------------!
   INTENT(IN)                   :: uf, nic, njc
   !------------------------------------------------------------------------!
   INTEGER                      :: nif, njf
   REAL, DIMENSION(nic,njc)     :: restrict
   !------------------------------------------------------------------------!
   nif = size(uf,1)
   njf = size(uf,2)

   IF ( nic < nif .AND. njc < njf ) THEN
   !restrict in both dir
      restrict(2:nic-1,2:njc-1)=0.0625*(&
          4.0*uf(3:nif-2:2,3:njf-2:2)+2.0*(uf(4:nif-1:2,3:njf-2:2)&
               +uf(2:nif-3:2,3:njf-2:2)&
               +uf(3:nif-2:2,4:njf-1:2)&
               +uf(3:nif-2:2,2:njf-3:2))&
            +uf(2:nif-3:2,2:njf-3:2)&
            +uf(2:nif-3:2,4:njf-1:2)&
            +uf(4:nif-1:2,2:njf-3:2)&
            +uf(4:nif-1:2,4:njf-1:2))

      restrict(2:nic-1,1)  =0.5*uf(3:nif-2:2,1)+0.25*(uf(2:nif-3:2,1)&
                          +uf(4:nif-1:2,1))
      restrict(2:nic-1,njc)=0.5*uf(3:nif-2:2,njf)+0.25*(uf(2:nif-3:2,njf)&
                          +uf(4:nif-1:2,njf))
      restrict(1,2:njc-1)  =0.5*uf(1,3:njf-2:2)+0.25*(uf(1,2:njf-3:2)&
                          +uf(1,4:njf-1:2))
      restrict(nic,2:njc-1)=0.5*uf(nif,3:njf-2:2)+0.25*(uf(nif,2:njf-3:2)&
                          +uf(nif,4:njf-1:2))
   ELSE IF ( njc < njf ) THEN
   ! ==> nic == nif
   !restrict only in j-dir
      restrict(2:nic-1,2:njc-1)=0.25*(&
          2.0*uf(2:nif-1,3:njf-2:2)&
            +uf(2:nif-1,4:njf-1:2)+uf(2:nif-1,2:njf-3:2))
      restrict(2:nic-1,1)    =uf(2:nif-1,1)       !only copy
      restrict(2:nic-1,njc)  =uf(2:nif-1,njf)     !only copy
      restrict(1,2:njc-1)  =0.5*uf(1,3:njf-2:2)  +0.25*(uf(1,2:njf-3:2)  +uf(1,4:njf-1:2))
      restrict(nic,2:njc-1)=0.5*uf(nif,3:njf-2:2)+0.25*(uf(nif,2:njf-3:2)+uf(nif,4:njf-1:2))
   ELSE 
   ! ==> njc == njf
   !restrict only in i-dir 
      restrict(2:nic-1,2:njc-1)=0.25*(&
         2.0*uf(3:nif-2:2,2:njf-1)&
            +uf(4:nif-1:2,2:njf-1)+uf(2:nif-3:2,2:njf-1))
      restrict(2:nic-1,1)  =0.5*uf(3:nif-2:2,1)  +0.25*(uf(2:nif-3:2,1)  +uf(4:nif-1:2,1))
      restrict(2:nic-1,njc)=0.5*uf(3:nif-2:2,njf)+0.25*(uf(2:nif-3:2,njf)+uf(4:nif-1:2,njf))
      restrict(1,2:njc-1)=uf(1,2:njf-1)         !only copy
      restrict(nic,2:njc-1)=uf(nif,2:njf-1)     !only copy
   END IF 

   !corner ghost cells
   restrict(1,1)     = uf(1,1)
   restrict(1,njc)   = uf(1,njf)
   restrict(nic,1)   = uf(nif,1)
   restrict(nic,njc) = uf(nif,njf)
 END FUNCTION restrict

PURE FUNCTION restrictBOUND(uf,nic,njc)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   REAL, DIMENSION(:,:)         :: uf
   INTEGER                      :: nic,njc
   !------------------------------------------------------------------------!
   INTENT(IN)                   :: uf,nic,njc
   !------------------------------------------------------------------------!
   INTEGER                      :: nif, njf
   REAL, DIMENSION( nic , njc ) :: restrictBOUND
   !------------------------------------------------------------------------!
   nif = size(uf,1)
   njf = size(uf,2)

   IF ( nic < nif .AND. njc < njf ) THEN
   !restrict in both dir
      restrictBOUND(2:nic-1,1)  =0.5*uf(3:nif-2:2,1)+0.25*(uf(2:nif-3:2,1)&
                          +uf(4:nif-1:2,1))
      restrictBOUND(2:nic-1,njc)=0.5*uf(3:nif-2:2,njf)+0.25*(uf(2:nif-3:2,njf)&
                          +uf(4:nif-1:2,njf))
      restrictBOUND(1,2:njc-1)  =0.5*uf(1,3:njf-2:2)+0.25*(uf(1,2:njf-3:2)&
                          +uf(1,4:njf-1:2))
      restrictBOUND(nic,2:njc-1)=0.5*uf(nif,3:njf-2:2)+0.25*(uf(nif,2:njf-3:2)&
                          +uf(nif,4:njf-1:2))
   ELSE IF ( njc < njf ) THEN
      restrictBOUND(2:nic-1,1)    =uf(2:nif-1,1)       !only copy
      restrictBOUND(2:nic-1,njc)  =uf(2:nif-1,njf)     !only copy
      restrictBOUND(1,2:njc-1)  =0.5*uf(1,3:njf-2:2)  +0.25*(uf(1,2:njf-3:2)  +uf(1,4:njf-1:2))
      restrictBOUND(nic,2:njc-1)=0.5*uf(nif,3:njf-2:2)+0.25*(uf(nif,2:njf-3:2)+uf(nif,4:njf-1:2))
   ELSE
      restrictBOUND(2:nic-1,1)  =0.5*uf(3:nif-2:2,1)  +0.25*(uf(2:nif-3:2,1)  +uf(4:nif-1:2,1))
      restrictBOUND(2:nic-1,njc)=0.5*uf(3:nif-2:2,njf)+0.25*(uf(2:nif-3:2,njf)+uf(4:nif-1:2,njf))
      restrictBOUND(1,2:njc-1)=uf(1,2:njf-1)         !only copy
      restrictBOUND(nic,2:njc-1)=uf(nif,2:njf-1)     !only copy
   END IF

   !corner ghost cells
   restrictBOUND(1,1)     = uf(1,1)
   restrictBOUND(1,njc)   = uf(1,njf)
   restrictBOUND(nic,1)   = uf(nif,1)
   restrictBOUND(nic,njc) = uf(nif,njf)
END FUNCTION restrictBOUND


 PURE FUNCTION prolong(uc,nif,njf)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:)  :: uc
    INTEGER               :: nif,njf 
    !------------------------------------------------------------------------!
    INTENT(IN)            :: uc, nif, njf
    !------------------------------------------------------------------------!
    INTEGER :: nic,njc
    REAL, DIMENSION(nif,njf) :: prolong
    !------------------------------------------------------------------------!
    nic = size(uc,1)
    njc = size(uc,2)
    IF (nif > nic .AND. njf > njc) THEN
       prolong(1:nif:2,1:njf:2)=uc(1:nic,1:njc)
       prolong(2:nif-1:2,1:njf:2)=0.5*(prolong(3:nif:2,1:njf:2)+prolong(1:nif-2:2,1:njf:2))
       prolong(1:nif,2:njf-1:2)=0.5*(prolong(1:nif,3:njf:2)+prolong(1:nif,1:njf-2:2))
    ELSE IF (nif > nic) THEN
       prolong(1:nif:2,1:njf)=uc(1:nic,1:njc) 
       prolong(2:nif-1:2,1:njf)=0.5*(prolong(3:nif:2,1:njf)+prolong(1:nif-2:2,1:njf))
    ELSE
       prolong(1:nif,1:njf:2)=uc(1:nic,1:njc) 
       prolong(1:nif,2:njf-1:2)=0.5*(prolong(1:nif,3:njf:2)+prolong(1:nif,1:njf-2:2))
    END IF
 END FUNCTION prolong

 PURE FUNCTION prolong2(uc,nif,njf)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:)   :: uc
    INTEGER                :: nif,njf
    !------------------------------------------------------------------------!
    INTEGER                :: nic,njc
    !------------------------------------------------------------------------!
    INTENT(IN)             :: uc,nif,njf
    !------------------------------------------------------------------------!
    REAL,DIMENSION(2:nif-1,2:njf-1):: prolong2
    !------------------------------------------------------------------------!
    nic = size(uc,1)
    njc = size(uc,2)
    IF (nif > nic .AND. njf > njc) THEN
       prolong2(3:nif-2:2,3:njf-2:2)=      uc(2:nic-1,2:njc-1)
       prolong2(2:nif-1:2,3:njf-2:2)= 0.5*(uc(1:nic-1,2:njc-1)+uc(2:nic,2:njc-1))
       prolong2(3:nif-2:2,2:njf-1:2)= 0.5*(uc(2:nic-1,1:njc-1)+uc(2:nic-1,2:njc))
       prolong2(2:nif-1:2,2:njf-1:2)=0.25*(uc(1:nic-1,1:njc-1)+uc(2:nic,1:njc-1)&
                                          +uc(1:nic-1,2:njc)  +uc(2:nic,2:njc))
    ELSE IF (nif > nic) THEN
       prolong2(3:nif-2:2,2:njf-1)  = uc(2:nic-1,2:njc-1)
       prolong2(2:nif-1:2,2:njf-1)  = 0.5*(uc(1:nic-1,2:njc-1)+uc(2:nic,2:njc-1))
    ELSE
       prolong2(2:nif-1,3:njf-2:2)  = uc(2:nic-1,2:njc-1)
       prolong2(2:nif-1,2:njf-1:2)  = 0.5*(uc(2:nic-1,1:njc-1)+uc(2:nic-1,2:njc))
    END IF
 END FUNCTION prolong2


 PURE SUBROUTINE relax(this,u,rhs,jgrid)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Sources_TYP)         :: this
   REAL, DIMENSION(:,:)      :: u,rhs
   INTEGER                   :: jgrid
   !------------------------------------------------------------------------!
   TYPE(Grid_TYP), POINTER   :: pgrid
   INTEGER                   :: ni,nj
   !------------------------------------------------------------------------!
   INTENT(INOUT)             :: this,u
   INTENT(IN)                :: rhs,jgrid   
   !------------------------------------------------------------------------!
   pgrid=>this%grid(jgrid)
   ni = pgrid%ni
   nj = pgrid%nj

   IF (this%Boundary(WEST) .EQ. NEUMANN) THEN
      ! Neumann
      u(1,2:nj-1) = u(2,2:nj-1)
   ELSE IF (this%Boundary(WEST) .EQ. PERIODIC) THEN
      ! Periodic
      u(1,2:nj-1) = u(ni-1,2:nj-1)
   END IF

   IF (this%Boundary(EAST) .EQ. NEUMANN) THEN
      ! Neumann
      u(ni,2:nj-1) = u(ni-1,2:nj-1)
   ELSE IF (this%Boundary(EAST) .EQ. PERIODIC) THEN
      ! Periodic
      u(ni,2:nj-1) = u(2,2:nj-1)
   END IF

   IF (this%Boundary(SOUTH) .EQ. NEUMANN) THEN
      ! Neumann
      u(2:ni-1,1) = u(2:ni-1,2)
   ELSE IF (this%Boundary(SOUTH) .EQ. PERIODIC) THEN
      ! Periodic
      u(2:ni-1,1) = u(2:ni-1,nj-1)
   END IF

   IF (this%Boundary(NORTH) .EQ. NEUMANN) THEN
      ! Neumann
      u(2:ni-1,nj) = u(2:ni-1,nj-1)
   ELSE IF (this%Boundary(NORTH) .EQ. PERIODIC) THEN
      ! Periodic
      u(2:ni-1,nj) = u(2:ni-1,2)
   END IF
   u(2:ni-1:2,2:nj-1:2) = pgrid%d(2:ni-1:2,2:nj-1:2) *( &
               + 0.125*pgrid%invhi2*( &
               +pgrid%da(2:ni-1:2,2:nj-1:2) *&
                      (u(3:ni:2,2:nj-1:2)-u(1:ni-2:2,2:nj-1:2))&
               + 4.0*pgrid%a(2:ni-1:2,2:nj-1:2)*&
                      (u(3:ni:2,2:nj-1:2)+u(1:ni-2:2,2:nj-1:2)))&
               + 0.125*pgrid%invhj2*( &
               +pgrid%db(2:ni-1:2,2:nj-1:2) *&
                      (u(2:ni-1:2,3:nj:2)-u(2:ni-1:2,1:nj-2:2))&
               + 4.0*pgrid%b(2:ni-1:2,2:nj-1:2)*&
                      (u(2:ni-1:2,3:nj:2)+u(2:ni-1:2,1:nj-2:2)))&
               - 0.5*pgrid%c(2:ni-1:2,2:nj-1:2)*rhs(2:ni-1:2,2:nj-1:2))
    u(3:ni-2:2,3:nj-2:2) = pgrid%d(3:ni-2:2,3:nj-2:2) *( &
               + 0.125*pgrid%invhi2*( &
               +pgrid%da(3:ni-2:2,3:nj-2:2) *&
                      (u(4:ni-1:2,3:nj-2:2)-u(2:ni-3:2,3:nj-2:2))&
               + 4.0*this%grid(jgrid)%a(3:ni-2:2,3:nj-2:2)*&
                      (u(4:ni-1:2,3:nj-2:2)+u(2:ni-3:2,3:nj-2:2)))&
               + 0.125*pgrid%invhj2*( &
               +pgrid%db(3:ni-2:2,3:nj-2:2) *&
                      (u(3:ni-2:2,4:nj-1:2)-u(3:ni-2:2,2:nj-3:2))&
               + 4.0*pgrid%b(3:ni-2:2,3:nj-2:2)*&
                      (u(3:ni-2:2,4:nj-1:2)+u(3:ni-2:2,2:nj-3:2)))&
               - 0.5*pgrid%c(3:ni-2:2,3:nj-2:2)*rhs(3:ni-2:2,3:nj-2:2))
    u(3:ni-2:2,2:nj-1:2) = pgrid%d(3:ni-2:2,2:nj-1:2) *( &
               + 0.125*pgrid%invhi2*( &
               +pgrid%da(3:ni-2:2,2:nj-1:2) *&
                      (u(4:ni-1:2,2:nj-1:2)-u(2:ni-3:2,2:nj-1:2))&
               + 4.0*pgrid%a(3:ni-2:2,2:nj-1:2)*&
                      (u(4:ni-1:2,2:nj-1:2)+u(2:ni-3:2,2:nj-1:2)))&
               + 0.125*pgrid%invhj2*( &
               +pgrid%db(3:ni-2:2,2:nj-1:2) *&
                      (u(3:ni-2:2,3:nj:2)-u(3:ni-2:2,1:nj-2:2))&
               + 4.0*pgrid%b(3:ni-2:2,2:nj-1:2)*&
                      (u(3:ni-2:2,3:nj:2)+u(3:ni-2:2,1:nj-2:2)))&
               - 0.5*pgrid%c(3:ni-2:2,2:nj-1:2)*rhs(3:ni-2:2,2:nj-1:2))
    u(2:ni-1:2,3:nj-2:2) = pgrid%d(2:ni-1:2,3:nj-2:2) *( &
               + 0.125*pgrid%invhi2*( &
               +pgrid%da(2:ni-1:2,3:nj-2:2) *&
                      (u(3:ni:2,3:nj-2:2)-u(1:ni-2:2,3:nj-2:2))&
               + 4.0*pgrid%a(2:ni-1:2,3:nj-2:2)*&
                      (u(3:ni:2,3:nj-2:2)+u(1:ni-2:2,3:nj-2:2)))&
               + 0.125*pgrid%invhj2*( &
               +pgrid%db(2:ni-1:2,3:nj-2:2) *&
                      (u(2:ni-1:2,4:nj-1:2)-u(2:ni-1:2,2:nj-3:2))&
               + 4.0*pgrid%b(2:ni-1:2,3:nj-2:2)*&
                      (u(2:ni-1:2,4:nj-1:2)+u(2:ni-1:2,2:nj-3:2)))&
              - 0.5*pgrid%c(2:ni-1:2,3:nj-2:2)*rhs(2:ni-1:2,3:nj-2:2))
END SUBROUTINE relax

PURE FUNCTION resid(this,u,rhs,jgrid)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Sources_TYP)               :: this
   REAL, DIMENSION(:,:)            :: u,rhs
   INTEGER                         :: jgrid
   !------------------------------------------------------------------------!
   INTENT(IN)                      :: this,u,rhs,jgrid
   !------------------------------------------------------------------------!
   INTEGER                         :: ni,nj
   REAL, DIMENSION(size(u,1),size(u,2)) :: resid
   !------------------------------------------------------------------------!
   ni = this%grid(jgrid)%ni
   nj = this%grid(jgrid)%nj
   resid(2:ni-1,2:nj-1) = -(this%grid(jgrid)%invc(2:ni-1,2:nj-1)*&
                           ( this%grid(jgrid)%invhi2*( &
                         + 0.25*this%grid(jgrid)%da(2:ni-1,2:nj-1) *&
                          (u(3:ni,2:nj-1)-u(1:ni-2,2:nj-1))&
                         + this%grid(jgrid)%a(2:ni-1,2:nj-1)*&
                           (u(3:ni,2:nj-1)+u(1:ni-2,2:nj-1)-2.0*u(2:ni-1,2:nj-1)))&
                         + this%grid(jgrid)%invhj2*( &
                         + 0.25*this%grid(jgrid)%db(2:ni-1,2:nj-1) *&
                           (u(2:ni-1,3:nj)-u(2:ni-1,1:nj-2))&
                         + this%grid(jgrid)%b(2:ni-1,2:nj-1)*&
                           (u(2:ni-1,3:nj)+u(2:ni-1,1:nj-2)-2.0*u(2:ni-1,2:nj-1))))&
                         - rhs(2:ni-1,2:nj-1))
   resid(1,1:nj) =0.0
   resid(ni,1:nj)=0.0
   resid(1:ni,1) =0.0
   resid(1:ni,nj)=0.0
END FUNCTION resid

SUBROUTINE CloseSources_selfgravitation(this)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Sources_TYP) :: this
   !------------------------------------------------------------------------!
   INTEGER           :: j
   !------------------------------------------------------------------------!
   INTENT(INOUT)     :: this
   !------------------------------------------------------------------------!
   DEALLOCATE(this%accel)
   CALL CloseSources_boundary(this)
!CDIR NODEP
    DO j=this%MGminlevel,this%ngrid
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
    END DO
    DEALLOCATE(this%grid)

  END SUBROUTINE CloseSources_selfgravitation

END MODULE sources_selfgravitation
