!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler2Disothm.f90                                         #
!#                                                                           #
!# Copyright (C) 2007-2012                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
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
!> \addtogroup physics
!! - isothermal gas dynamics
!!   \key{cs,REAL,isothermal sound speed,0.0}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Björn Sperling
!!
!! \brief basic module for 2D isothermal Euler equations
!!
!! \extends physics_common
!! \ingroup physics
!----------------------------------------------------------------------------!
MODULE physics_euler2Disothm
  USE physics_common
  USE sources_common, ONLY : Sources_TYP
  USE mesh_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE CalcWaveSpeeds_euler2Dit
     MODULE PROCEDURE CalcWaveSpeeds_center
     MODULE PROCEDURE CalcWaveSpeeds_faces
  END INTERFACE
     INTERFACE GeometricalSources_euler2Dit
     MODULE PROCEDURE GeometricalSources_center
     MODULE PROCEDURE GeometricalSources_faces
  END INTERFACE
  INTERFACE Convert2Primitive_euler2Dit
     MODULE PROCEDURE Convert2Primitive_center
     MODULE PROCEDURE Convert2Primitive_faces
  END INTERFACE
  INTERFACE Convert2Conservative_euler2Dit
     MODULE PROCEDURE Convert2Conservative_center
     MODULE PROCEDURE Convert2Conservative_faces
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: num_var = 3              ! number of variables       !
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 2D isotherm"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Physics_TYP, &
       ! methods
       InitPhysics_euler2Dit, &
       ClosePhysics_euler2Dit, &
       CalcWaveSpeeds_euler2Dit, &
       CalcFluxesX_euler2Dit, &
       CalcFluxesY_euler2Dit, &
       CalcCharSystemX_euler2Dit, &
       CalcCharSystemY_euler2Dit, &
       CalcBoundaryDataX_euler2Dit, &
       CalcBoundaryDataY_euler2Dit, &
       CalcStresses_euler2Dit, &
       GeometricalSources_euler2Dit, &
       ViscositySources_euler2Dit, &
       ExternalSources_euler2Dit, &
       Convert2Primitive_euler2Dit, &
       Convert2Conservative_euler2Dit, &
       ReflectionMasks_euler2Dit, &
       AxisMasks_euler2Dit, &
       SetWaveSpeeds_euler2Dit, &
       SetEigenValues_euler2Dit, &
       SetBoundaryData_euler2Dit, &
       SetCharVars_euler2Dit, &
       CalcFlux_euler2Dit, &
       MomentumSourcesX_euler2Dit, &
       MomentumSourcesY_euler2Dit, &
       Cons2Prim_euler2Dit, &
       Prim2Cons_euler2Dit, &
       ClosePhysics, &
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

  SUBROUTINE InitPhysics_euler2Dit(this,Mesh,problem,pname,nvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: problem
    CHARACTER(LEN=32), OPTIONAL :: pname
    INTEGER,OPTIONAL  :: nvar
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,problem,pname
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    IF (PRESENT(pname).AND.PRESENT(nvar)) THEN
       CALL InitPhysics(this,problem,pname,nvar)
    ELSE IF (PRESENT(pname).OR.PRESENT(nvar)) THEN
       CALL Error(this, "InitPhysics_euler2Dit", "Both or no optional " &
        // "arguments at all have to be defined.")
    ELSE
       CALL InitPhysics(this,problem,problem_name,num_var)
    END IF
    ! set array indices
    this%DENSITY   = 1                                 ! mass density        !
    this%XVELOCITY = 2                                 ! x-velocity          !
    this%XMOMENTUM = 2                                 ! x-momentum          !
    this%YVELOCITY = 3                                 ! y-velocity          !
    this%YMOMENTUM = 3                                 ! y-momentum          !
    this%ZVELOCITY = 0                                 ! no z-velocity       !
    this%ZMOMENTUM = 0                                 ! no z-momentum       !
    this%PRESSURE  = 0                                 ! no pressure         !
    this%ENERGY    = 0                                 ! no total energy     !
    ! set names for primitive and conservative variables
    this%pvarname(this%DENSITY)   = "density"
    this%pvarname(this%XVELOCITY) = "xvelocity"
    this%pvarname(this%YVELOCITY) = "yvelocity"
    this%cvarname(this%DENSITY)   = "density"
    this%cvarname(this%XMOMENTUM) = "xmomentum"
    this%cvarname(this%YMOMENTUM) = "ymomentum"
    this%DIM = 2

  END SUBROUTINE InitPhysics_euler2Dit

  PURE SUBROUTINE CalcWaveSpeeds_center(this,Mesh,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
         :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! compute minimal and maximal wave speeds at cell centers
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! x-direction
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2Dit(this%bccsound(i,j),pvar(i,j,this%XVELOCITY),&
               this%amin(i,j),this%amax(i,j))
          ! y-direction
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2Dit(this%bccsound(i,j),pvar(i,j,this%YVELOCITY),&
               this%bmin(i,j),this%bmax(i,j))
       END DO
    END DO
  END SUBROUTINE CalcWaveSpeeds_center


  PURE SUBROUTINE CalcWaveSpeeds_faces(this,Mesh,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%VNUM) &
         :: prim
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,prim
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! compute minimal and maximal wave speeds at cell interfaces
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! western
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2Dit(&
               this%fcsound(i,j,1),&
               prim(i,j,1,this%XVELOCITY), &
               this%tmin(i,j,1),this%tmax(i,j,1))
          ! eastern
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2Dit(&
               this%fcsound(i,j,2),&
               prim(i,j,2,this%XVELOCITY), &
               this%amin(i,j),this%amax(i,j))
          ! southern
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2Dit(&
               this%fcsound(i,j,3), &
               prim(i,j,3,this%YVELOCITY), &
               this%tmin(i,j,2),this%tmax(i,j,2))
          ! northern
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2Dit(&
               this%fcsound(i,j,4), &
               prim(i,j,4,this%YVELOCITY), &
               this%bmin(i,j),this%bmax(i,j))
       END DO
    END DO
    ! set minimal and maximal wave speeds at cell interfaces of neighboring cells
    DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
       DO i=Mesh%IMIN-1,Mesh%IMAX
          ! western interfaces
          this%amin(i,j) = MIN(this%tmin(i+1,j,1),this%amin(i,j))
          ! eastern interfaces
          this%amax(i,j) = MAX(this%tmax(i+1,j,1),this%amax(i,j))
       END DO
    END DO
!CDIR COLLAPSE 
    DO j=Mesh%JMIN-1,Mesh%JMAX
!CDIR NODEP
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! southern interfaces
          this%bmin(i,j) = MIN(this%tmin(i,j+1,2),this%bmin(i,j))
          ! northern interfaces
          this%bmax(i,j) = MAX(this%tmax(i,j+1,2),this%bmax(i,j))
       END DO
    END DO
  END SUBROUTINE CalcWaveSpeeds_faces


  PURE SUBROUTINE CalcFluxesX_euler2Dit(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: nmin,nmax
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%VNUM) &
         :: prim,cons,xfluxes
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,nmin,nmax,prim,cons
    INTENT(OUT)       :: xfluxes
    !------------------------------------------------------------------------!
    CALL CalcFlux_euler2Dit(this%fcsound(:,:,nmin:nmax),&
         prim(:,:,nmin:nmax,this%DENSITY), &
         prim(:,:,nmin:nmax,this%XVELOCITY),cons(:,:,nmin:nmax,this%XMOMENTUM),&
         cons(:,:,nmin:nmax,this%YMOMENTUM),xfluxes(:,:,nmin:nmax,this%DENSITY),&
         xfluxes(:,:,nmin:nmax,this%XMOMENTUM),xfluxes(:,:,nmin:nmax,this%YMOMENTUM))
   END SUBROUTINE CalcFluxesX_euler2Dit


  PURE SUBROUTINE CalcFluxesY_euler2Dit(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: nmin,nmax
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%VNUM) &
         :: prim,cons,yfluxes
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,nmin,nmax,prim,cons
    INTENT(OUT)       :: yfluxes
    !------------------------------------------------------------------------!
    CALL CalcFlux_euler2Dit(this%fcsound(:,:,nmin:nmax),&
         prim(:,:,nmin:nmax,this%DENSITY), &
         prim(:,:,nmin:nmax,this%YVELOCITY),cons(:,:,nmin:nmax,this%YMOMENTUM), &
         cons(:,:,nmin:nmax,this%XMOMENTUM),yfluxes(:,:,nmin:nmax,this%DENSITY), &
         yfluxes(:,:,nmin:nmax,this%YMOMENTUM),yfluxes(:,:,nmin:nmax,this%XMOMENTUM))
  END SUBROUTINE CalcFluxesY_euler2Dit


  PURE SUBROUTINE CalcCharSystemX_euler2Dit(this,Mesh,i,dir,pvar,lambda,xvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: i,dir
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: lambda,xvar
    !------------------------------------------------------------------------!
    INTEGER           :: i1,i2
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,i,dir,pvar
    INTENT(INOUT)     :: lambda
    INTENT(OUT)       :: xvar
    !------------------------------------------------------------------------!
    ! compute eigenvalues at i
    CALL SetEigenValues_euler2Dit(this%csiso,pvar(i,:,this%XVELOCITY), &
         lambda(:,1),lambda(:,2),lambda(:,3))
    ! compute characteristic variables
    i1 = i + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
    i2 = MAX(i,i1)
    i1 = MIN(i,i1)
    CALL SetCharVars_euler2Dit(this%csiso,pvar(i1,:,this%DENSITY), &
         pvar(i2,:,this%DENSITY),pvar(i1,:,this%XVELOCITY), &
         pvar(i2,:,this%XVELOCITY),pvar(i1,:,this%YVELOCITY), &
         pvar(i2,:,this%YVELOCITY),xvar(:,1),xvar(:,2),xvar(:,3))
  END SUBROUTINE CalcCharSystemX_euler2Dit


  PURE SUBROUTINE CalcCharSystemY_euler2Dit(this,Mesh,j,dir,pvar,lambda,xvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: j,dir
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,this%VNUM) :: lambda,xvar
    !------------------------------------------------------------------------!
    INTEGER           :: j1,j2
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,j,dir,pvar
    INTENT(INOUT)     :: lambda
    INTENT(OUT)       :: xvar
    !------------------------------------------------------------------------!
    ! compute eigenvalues at j
    CALL SetEigenValues_euler2Dit(this%csiso,pvar(:,j,this%YVELOCITY), &
         lambda(:,1),lambda(:,2),lambda(:,3))
    ! compute characteristic variables
    j1 = j + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
    j2 = MAX(j,j1)
    j1 = MIN(j,j1)
    CALL SetCharVars_euler2Dit(this%csiso,pvar(:,j1,this%DENSITY), &
         pvar(:,j2,this%DENSITY),pvar(:,j1,this%YVELOCITY), &
         pvar(:,j2,this%YVELOCITY),pvar(:,j1,this%XVELOCITY), &
         pvar(:,j2,this%XVELOCITY),xvar(:,1),xvar(:,2),xvar(:,3))
  END SUBROUTINE CalcCharSystemY_euler2Dit


  PURE SUBROUTINE CalcBoundaryDataX_euler2Dit(this,Mesh,i1,dir,xvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: i1,dir
    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: xvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: i2
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,i1,dir,xvar
    INTENT(INOUT)     :: pvar
    !------------------------------------------------------------------------!
    i2 = i1 + SIGN(1,dir)  ! i +/- 1 depending on the sign of dir
    CALL SetBoundaryData_euler2Dit(this%csiso,1.0*SIGN(1,dir), &
         pvar(i1,:,this%DENSITY),pvar(i1,:,this%XVELOCITY), &
         pvar(i1,:,this%YVELOCITY),xvar(:,1),xvar(:,2),xvar(:,3), &
         pvar(i2,:,this%DENSITY),pvar(i2,:,this%XVELOCITY),pvar(i2,:,this%YVELOCITY))
  END SUBROUTINE CalcBoundaryDataX_euler2Dit


  PURE SUBROUTINE CalcBoundaryDataY_euler2Dit(this,Mesh,j1,dir,xvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: j1,dir
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,this%VNUM) :: xvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: j2
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,j1,dir,xvar
    INTENT(INOUT)     :: pvar
    !------------------------------------------------------------------------!
    j2 = j1 + SIGN(1,dir)  ! j +/- 1 depending on the sign of dir
    CALL SetBoundaryData_euler2Dit(this%csiso,1.0*SIGN(1,dir), &
         pvar(:,j1,this%DENSITY),pvar(:,j1,this%YVELOCITY), &
         pvar(:,j1,this%XVELOCITY),xvar(:,1),xvar(:,2),xvar(:,3), &
         pvar(:,j2,this%DENSITY),pvar(:,j2,this%YVELOCITY),pvar(:,j2,this%XVELOCITY))
  END SUBROUTINE CalcBoundaryDataY_euler2Dit


  PURE SUBROUTINE CalcStresses_euler2Dit(this,Mesh,pvar,dynvis,bulkvis, &
       btxx,btxy,btyy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
         dynvis,bulkvis,btxx,btxy,btyy
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,dynvis,bulkvis
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: btxx,btxy,btyy
    !------------------------------------------------------------------------!
    ! compute components of the stress tensor at cell bary centers
    ! inside the computational domain including one slice of ghost cells

    ! compute bulk viscosity first and store the result in this%tmp
!CDIR IEXPAND
    CALL Divergence(Mesh,pvar(:,:,this%XVELOCITY),pvar(:,:,this%YVELOCITY),this%tmp(:,:))
    this%tmp(:,:) = bulkvis(:,:)*this%tmp(:,:)

!CDIR OUTERUNROLL=8
    DO j=Mesh%JMIN-1,Mesh%JMAX+1
!CDIR NODEP
       DO i=Mesh%IMIN-1,Mesh%IMAX+1
          ! compute the diagonal elements of the stress tensor
          btxx(i,j) = dynvis(i,j) * &
               ( (pvar(i+1,j,this%XVELOCITY) - pvar(i-1,j,this%XVELOCITY)) / Mesh%dlx(i,j) &
               + 2.0 * Mesh%cxyx(i,j,1) * pvar(i,j,this%YVELOCITY) ) &
               + this%tmp(i,j)
               
          btyy(i,j) = dynvis(i,j) * &
               ( (pvar(i,j+1,this%YVELOCITY) - pvar(i,j-1,this%YVELOCITY)) / Mesh%dly(i,j) &
               + 2.0 * Mesh%cyxy(i,j,1) * pvar(i,j,this%XVELOCITY) ) &
               + this%tmp(i,j)

          ! compute the off-diagonal elements (no bulk viscosity)
          btxy(i,j) = dynvis(i,j) * ( 0.5 * &
               ( (pvar(i+1,j,this%YVELOCITY) - pvar(i-1,j,this%YVELOCITY)) / Mesh%dlx(i,j) &
               + (pvar(i,j+1,this%XVELOCITY) - pvar(i,j-1,this%XVELOCITY)) / Mesh%dly(i,j) ) &
               - Mesh%cxyx(i,j,1) * pvar(i,j,this%XVELOCITY) &
               - Mesh%cyxy(i,j,1) * pvar(i,j,this%YVELOCITY) )
       END DO
    END DO
  END SUBROUTINE CalcStresses_euler2Dit


  PURE SUBROUTINE GeometricalSources_center(this,Mesh,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
         :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,pvar,cvar
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! no geometrical density sources
          sterm(i,j,this%DENSITY) = 0.
          ! geometrical source terms in momentum equationes
          sterm(i,j,this%XMOMENTUM) = MomentumSourcesX_euler2Dit(&
              cvar(i,j,this%YMOMENTUM),pvar(i,j,this%XVELOCITY),&
              pvar(i,j,this%YVELOCITY),pvar(i,j,this%DENSITY)*this%bccsound(i,j)**2, &
              Mesh%cxyx(i,j,1),Mesh%cyxy(i,j,1),Mesh%czxz(i,j,1))
          sterm(i,j,this%YMOMENTUM) = MomentumSourcesY_euler2Dit(&
              cvar(i,j,this%XMOMENTUM),pvar(i,j,this%XVELOCITY),&
              pvar(i,j,this%YVELOCITY),pvar(i,j,this%DENSITY)*this%bccsound(i,j)**2, &
              Mesh%cxyx(i,j,1),Mesh%cyxy(i,j,1),Mesh%czxz(i,j,1))
       END DO
    END DO
  END SUBROUTINE GeometricalSources_center


  PURE SUBROUTINE GeometricalSources_faces(this,Mesh,prim,cons,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%VNUM) &
         :: prim,cons
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
         :: sterm
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,prim,cons
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! no geometrical density sources
          sterm(i,j,this%DENSITY)   = 0.
          ! momentum sources (sum up corner values, don't use SUM function,
          ! because it prevents COLLAPSING and causes poor vectorization
          sterm(i,j,this%XMOMENTUM) = MomentumSourcesX_euler2Dit(&
              cons(i,j,1,this%YMOMENTUM),prim(i,j,1,this%XVELOCITY),&
              prim(i,j,1,this%YVELOCITY),prim(i,j,1,this%DENSITY)*this%fcsound(i,j,1)**2, &
              Mesh%cxyx(i,j,1),Mesh%cyxy(i,j,1),Mesh%czxz(i,j,1)) &
            + MomentumSourcesX_euler2Dit(&
              cons(i,j,2,this%YMOMENTUM),prim(i,j,2,this%XVELOCITY),&
              prim(i,j,2,this%YVELOCITY),prim(i,j,2,this%DENSITY)*this%fcsound(i,j,2)**2, &
              Mesh%cxyx(i,j,2),Mesh%cyxy(i,j,2),Mesh%czxz(i,j,2)) &
            + MomentumSourcesX_euler2Dit(&
              cons(i,j,3,this%YMOMENTUM),prim(i,j,3,this%XVELOCITY),&
              prim(i,j,3,this%YVELOCITY),prim(i,j,3,this%DENSITY)*this%fcsound(i,j,3)**2, &
              Mesh%cxyx(i,j,3),Mesh%cyxy(i,j,3),Mesh%czxz(i,j,3)) &
            + MomentumSourcesX_euler2Dit(&
              cons(i,j,4,this%YMOMENTUM),prim(i,j,4,this%XVELOCITY),&
              prim(i,j,4,this%YVELOCITY),prim(i,j,4,this%DENSITY)*this%fcsound(i,j,4)**2, &
              Mesh%cxyx(i,j,4),Mesh%cyxy(i,j,4),Mesh%czxz(i,j,4))

          sterm(i,j,this%YMOMENTUM) = MomentumSourcesY_euler2Dit(&
              cons(i,j,1,this%XMOMENTUM),prim(i,j,1,this%XVELOCITY),&
              prim(i,j,1,this%YVELOCITY),prim(i,j,1,this%DENSITY)*this%fcsound(i,j,1)**2, &
              Mesh%cxyx(i,j,1),Mesh%cyxy(i,j,1),Mesh%czyz(i,j,1)) &
            + MomentumSourcesY_euler2Dit(&
              cons(i,j,2,this%XMOMENTUM),prim(i,j,2,this%XVELOCITY),&
              prim(i,j,2,this%YVELOCITY),prim(i,j,2,this%DENSITY)*this%fcsound(i,j,2)**2, &
              Mesh%cxyx(i,j,2),Mesh%cyxy(i,j,2),Mesh%czyz(i,j,2)) &
            + MomentumSourcesY_euler2Dit(&
              cons(i,j,3,this%XMOMENTUM),prim(i,j,3,this%XVELOCITY),&
              prim(i,j,3,this%YVELOCITY),prim(i,j,3,this%DENSITY)*this%fcsound(i,j,3)**2, &
              Mesh%cxyx(i,j,3),Mesh%cyxy(i,j,3),Mesh%czyz(i,j,3)) &
            + MomentumSourcesY_euler2Dit(&
              cons(i,j,4,this%XMOMENTUM),prim(i,j,4,this%XVELOCITY),&
              prim(i,j,4,this%YVELOCITY),prim(i,j,4,this%DENSITY)*this%fcsound(i,j,4)**2, &
              Mesh%cxyx(i,j,4),Mesh%cyxy(i,j,4),Mesh%czyz(i,j,4))
       END DO
    END DO
  END SUBROUTINE GeometricalSources_faces


  ! momentum and energy sources due to external force
  PURE SUBROUTINE ExternalSources_euler2Dit(this,Mesh,accel,pvar,cvar,sterm)
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
         :: accel
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
         :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,accel,pvar,cvar
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          sterm(i,j,this%DENSITY)   = 0.
          sterm(i,j,this%XMOMENTUM) = pvar(i,j,this%DENSITY) * accel(i,j,1)
          sterm(i,j,this%YMOMENTUM) = pvar(i,j,this%DENSITY) * accel(i,j,2)
       END DO
    END DO
  END SUBROUTINE ExternalSources_euler2Dit


  PURE SUBROUTINE ViscositySources_euler2Dit(this,Mesh,pvar,btxx,btxy,btyy,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: &
         pvar,sterm
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
         btxx,btxy,btyy
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,pvar,btxx,btxy,btyy
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! mean values of stress tensor components across the cell interfaces

    ! viscosity source terms
    sterm(:,:,this%DENSITY) = 0.0 

    ! compute viscous momentum sources
    ! divergence of stress tensor with symmetry btyx=btxy
!CDIR IEXPAND
    CALL Divergence(Mesh,btxx,btxy,btxy,btyy,sterm(:,:,this%XMOMENTUM), &
                    sterm(:,:,this%YMOMENTUM))
  END SUBROUTINE ViscositySources_euler2Dit


  PURE SUBROUTINE Convert2Primitive_center(this,Mesh,i1,i2,j1,j2,cvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    INTEGER              :: i1,i2,j1,j2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
                         :: cvar,pvar
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,i1,i2,j1,j2,cvar
    INTENT(OUT) :: pvar
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL Cons2Prim_euler2Dit(cvar(i1:i2,j1:j2,this%DENSITY), &
         cvar(i1:i2,j1:j2,this%XMOMENTUM),cvar(i1:i2,j1:j2,this%YMOMENTUM), &
         pvar(i1:i2,j1:j2,this%DENSITY),pvar(i1:i2,j1:j2,this%XVELOCITY), &
         pvar(i1:i2,j1:j2,this%YVELOCITY))
  END SUBROUTINE Convert2Primitive_center


  PURE SUBROUTINE Convert2Primitive_faces(this,Mesh,i1,i2,j1,j2,cons,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    INTEGER              :: i1,i2,j1,j2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%VNUM) &
                         :: cons,prim
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,i1,i2,j1,j2,cons
    INTENT(OUT) :: prim
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL Cons2Prim_euler2Dit(cons(i1:i2,j1:j2,:,this%DENSITY), &
         cons(i1:i2,j1:j2,:,this%XMOMENTUM),cons(i1:i2,j1:j2,:,this%YMOMENTUM), &
         prim(i1:i2,j1:j2,:,this%DENSITY),prim(i1:i2,j1:j2,:,this%XVELOCITY), &
         prim(i1:i2,j1:j2,:,this%YVELOCITY))
  END SUBROUTINE Convert2Primitive_faces


  PURE SUBROUTINE Convert2Conservative_center(this,Mesh,i1,i2,j1,j2,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    INTEGER              :: i1,i2,j1,j2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
                         :: cvar,pvar
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,i1,i2,j1,j2,pvar
    INTENT(OUT) :: cvar
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL Prim2Cons_euler2Dit(pvar(i1:i2,j1:j2,this%DENSITY), &
         pvar(i1:i2,j1:j2,this%XVELOCITY), pvar(i1:i2,j1:j2,this%YVELOCITY), &
         cvar(i1:i2,j1:j2,this%DENSITY),cvar(i1:i2,j1:j2,this%XMOMENTUM), &
         cvar(i1:i2,j1:j2,this%YMOMENTUM))
  END SUBROUTINE Convert2Conservative_center


  PURE SUBROUTINE Convert2Conservative_faces(this,Mesh,i1,i2,j1,j2,prim,cons)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    INTEGER              :: i1,i2,j1,j2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%VNUM) &
                         :: cons,prim
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,i1,i2,j1,j2,prim
    INTENT(OUT) :: cons
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL Prim2Cons_euler2Dit(prim(i1:i2,j1:j2,:,this%DENSITY), &
         prim(i1:i2,j1:j2,:,this%XVELOCITY),prim(i1:i2,j1:j2,:,this%YVELOCITY), &
         cons(i1:i2,j1:j2,:,this%DENSITY),cons(i1:i2,j1:j2,:,this%XMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%YMOMENTUM))
  END SUBROUTINE Convert2Conservative_faces


  PURE SUBROUTINE ReflectionMasks_euler2Dit(this,reflX,reflY)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    LOGICAL, DIMENSION(this%VNUM) :: reflX,reflY
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this
    INTENT(OUT)       :: reflX,reflY
    !------------------------------------------------------------------------!
    ! western / eastern boundary
    reflX(this%DENSITY)   = .FALSE.
    reflX(this%XVELOCITY) = .TRUE.
    reflX(this%YVELOCITY) = .FALSE.
    ! southern / northern boundary
    reflY(this%DENSITY)   = .FALSE.
    reflY(this%XVELOCITY) = .FALSE.
    reflY(this%YVELOCITY) = .TRUE.
  END SUBROUTINE ReflectionMasks_euler2Dit


  PURE SUBROUTINE AxisMasks_euler2Dit(this,reflX,reflY)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    LOGICAL, DIMENSION(this%VNUM) :: reflX,reflY
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this
    INTENT(OUT)       :: reflX,reflY
    !------------------------------------------------------------------------!
    ! western / eastern boundary
    reflX(this%DENSITY)   = .FALSE.
    reflX(this%XVELOCITY) = .TRUE.
    reflX(this%YVELOCITY) = .TRUE.
    ! southern / northern boundary
    reflY(this%DENSITY)   = .FALSE.
    reflY(this%XVELOCITY) = .TRUE.
    reflY(this%YVELOCITY) = .TRUE.
  END SUBROUTINE AxisMasks_euler2Dit


  ELEMENTAL SUBROUTINE SetWaveSpeeds_euler2Dit(cs,v,amin,amax)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,v
    REAL, INTENT(OUT) :: amin,amax
    !------------------------------------------------------------------------!
    ! minimal and maximal wave speeds
    amin = MIN(0.,v-cs)
    amax = MAX(0.,v+cs)
  END SUBROUTINE SetWaveSpeeds_euler2Dit


  ELEMENTAL SUBROUTINE SetEigenValues_euler2Dit(cs,v,l1,l2,l3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,v
    REAL, INTENT(OUT) :: l1,l2,l3
    !------------------------------------------------------------------------!
    ! all eigenvalues of the isothermal euler problem
    l1 = v - cs
    l2 = v
    l3 = v + cs
  END SUBROUTINE SetEigenValues_euler2Dit


  ELEMENTAL SUBROUTINE SetCharVars_euler2Dit(cs,rho1,rho2,u1,u2,v1,v2,&
       xvar1,xvar2,xvar3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,rho1,rho2,u1,u2,v1,v2
    REAL, INTENT(OUT) :: xvar1,xvar2,xvar3
    !------------------------------------------------------------------------!
    REAL :: dlnrho,du
    !------------------------------------------------------------------------!
    dlnrho = LOG(rho2/rho1)
    du = u2-u1
    ! characteristic variables
    xvar1 = cs*dlnrho - du
    xvar2 = v2-v1
    xvar3 = cs*dlnrho + du
  END SUBROUTINE SetCharVars_euler2Dit


  ELEMENTAL SUBROUTINE SetBoundaryData_euler2Dit(cs,dir,rho1,u1,v1,xvar1,xvar2, &
       xvar3,rho2,u2,v2)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,dir,rho1,u1,v1,xvar1,xvar2,xvar3
    REAL, INTENT(OUT) :: rho2,u2,v2
    !------------------------------------------------------------------------!
    ! extrapolate boundary values using characteristic variables
    rho2 = rho1 * EXP(dir*0.5*(xvar3+xvar1)/cs)
    u2   = u1 + dir*0.5*(xvar3-xvar1)
    v2   = v1 + dir*xvar2
  END SUBROUTINE SetBoundaryData_euler2Dit


  ELEMENTAL SUBROUTINE CalcFlux_euler2Dit(cs,rho,v,m1,m2,f1,f2,f3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,rho,v,m1,m2
    REAL, INTENT(OUT) :: f1, f2, f3
    !------------------------------------------------------------------------!
    f1 = rho*v
    f2 = m1*v + rho*cs*cs
    f3 = m2*v
  END SUBROUTINE CalcFlux_euler2Dit


  ! momentum source terms due to inertial forces
  ! P is the isothermal pressure rho*cs*cs
  ELEMENTAL FUNCTION MomentumSourcesX_euler2Dit(my,vx,vy,P,cxyx,cyxy,czxz) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: my,vx,vy,P,cxyx,cyxy,czxz
    REAL :: st
    !------------------------------------------------------------------------!
    st = -my * (cxyx * vx - cyxy * vy) + (cyxy + czxz) * P
  END FUNCTION MomentumSourcesX_euler2Dit

  ELEMENTAL FUNCTION MomentumSourcesY_euler2Dit(mx,vx,vy,P,cxyx,cyxy,czyz) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: mx,vx,vy,P,cxyx,cyxy,czyz
    REAL :: st
    !------------------------------------------------------------------------!
    st = mx * (cxyx * vx - cyxy * vy) + (cxyx + czyz) * P
  END FUNCTION MomentumSourcesY_euler2Dit


  ELEMENTAL SUBROUTINE Cons2Prim_euler2Dit(rho_in,mu,mv,rho_out,u,v)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho_in,mu,mv
    REAL, INTENT(OUT) :: rho_out,u,v
    !------------------------------------------------------------------------!
    REAL :: inv_rho
    !------------------------------------------------------------------------!
    inv_rho = 1./rho_in
    rho_out = rho_in
    u = mu * inv_rho
    v = mv * inv_rho
  END SUBROUTINE Cons2Prim_euler2Dit

  
  ELEMENTAL SUBROUTINE Prim2Cons_euler2Dit(rho_in,u,v,rho_out,mu,mv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho_in,u,v
    REAL, INTENT(OUT) :: rho_out,mu,mv
    !------------------------------------------------------------------------!
    rho_out = rho_in
    mu = rho_in * u
    mv = rho_in * v
  END SUBROUTINE Prim2Cons_euler2Dit
  

  SUBROUTINE ClosePhysics_euler2Dit(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL ClosePhysics(this)
  END SUBROUTINE ClosePhysics_euler2Dit

END MODULE physics_euler2Disothm
