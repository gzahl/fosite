!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler3Drotsym.f90                                         #
!#                                                                           #
!# Copyright (C) 2007-2013                                                   #
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

!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Björn Sperling
!!
!! \brief basic module for 3D Euler equations with rotational symmetry
!!
!! \extends physics_common
!! \ingroup physics
!----------------------------------------------------------------------------!
MODULE physics_euler3Drotsym
  USE physics_common
  USE mesh_common, ONLY : Mesh_TYP
  USE sources_common, ONLY : Sources_TYP
  USE physics_euler2D, &
       CalcWaveSpeeds_euler3Drs => CalcWaveSpeeds_euler2D, &
       CalcSoundSpeeds_euler3Drs => CalcSoundSpeeds_euler2D
  USE mesh_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE GeometricalSources_euler3Drs
     MODULE PROCEDURE GeometricalSources_center
     MODULE PROCEDURE GeometricalSources_faces
  END INTERFACE
  INTERFACE Convert2Primitive_euler3Drs
     MODULE PROCEDURE Convert2Primitive_center
     MODULE PROCEDURE Convert2Primitive_faces
  END INTERFACE
  INTERFACE Convert2Conservative_euler3Drs
     MODULE PROCEDURE Convert2Conservative_center
     MODULE PROCEDURE Convert2Conservative_faces
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: num_var = 5              ! number of variables       !
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 3D w/ rot. symmetry"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitPhysics_euler3Drs, &
       CalcWaveSpeeds_euler3Drs, &
       CalcSoundSpeeds_euler3Drs, &
       CalcFluxesX_euler3Drs, &
       CalcFluxesY_euler3Drs, &
       CalcCharSystemX_euler3Drs, &
       CalcCharSystemY_euler3Drs, &
       CalcBoundaryDataX_euler3Drs, &
       CalcBoundaryDataY_euler3Drs, &
       CalcPrim2RiemannX_euler3Drs, &
       CalcPrim2RiemannY_euler3Drs, &
       CalcRiemann2PrimX_euler3Drs, &
       CalcRiemann2PrimY_euler3Drs, &
       CalcStresses_euler3drs, &
       GeometricalSources_euler3Drs, &
       ExternalSources_euler3Drs, &
       ViscositySources_euler3Drs, &
       Convert2Primitive_euler3Drs, &
       Convert2Conservative_euler3Drs, &
       ReflectionMasks_euler3Drs, &
       AxisMasks_euler3Drs, &
       SetEigenValues_euler3Drs, &
       ClosePhysics_euler3Drs
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitPhysics_euler3Drs(this,Mesh,problem)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: problem
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,problem
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL InitPhysics(this,problem,problem_name,num_var)
    ! set array indices
    this%DENSITY   = 1                                 ! mass density        !
    this%PRESSURE  = num_var                           ! pressure            !
    this%ENERGY    = num_var                           ! total energy        !
    this%XVELOCITY = 2                                 ! x-velocity          !
    this%XMOMENTUM = 2                                 ! x-momentum          !
    this%YVELOCITY = 3                                 ! y-velocity          !
    this%YMOMENTUM = 3                                 ! y-momentum          !
    this%ZVELOCITY = 4                                 ! rotational velocity !
    this%ZMOMENTUM = 4                                 ! rotational momentum !
    ! set names for primitive and conservative variables
    this%pvarname(this%DENSITY)   = "density"
    this%pvarname(this%XVELOCITY) = "xvelocity"
    this%pvarname(this%YVELOCITY) = "yvelocity"
    this%pvarname(this%ZVELOCITY) = "zvelocity"
    this%pvarname(this%PRESSURE)  = "pressure"
    this%cvarname(this%DENSITY)   = "density"
    this%cvarname(this%XMOMENTUM) = "xmomentum"
    this%cvarname(this%YMOMENTUM) = "ymomentum"
    this%cvarname(this%ZMOMENTUM) = "zmomentum"
    this%cvarname(this%ENERGY)    = "energy"
    this%DIM = 3

    ALLOCATE(this%fcent(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,2), &
         STAT = err)
    ! abort if allocation fails
    IF (err.NE.0) &
         CALL Error(this, "InitPhysics_euler3Drssgs", "Unable to allocate memory.")

  END SUBROUTINE InitPhysics_euler3Drs


  PURE SUBROUTINE CalcFluxesX_euler3Drs(this,Mesh,nmin,nmax,prim,cons,xfluxes)
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
    CALL CalcFlux_euler3Drs(prim(:,:,nmin:nmax,this%DENSITY), &
         prim(:,:,nmin:nmax,this%XVELOCITY),prim(:,:,nmin:nmax,this%PRESSURE),&
         cons(:,:,nmin:nmax,this%XMOMENTUM),cons(:,:,nmin:nmax,this%YMOMENTUM), &
         cons(:,:,nmin:nmax,this%ZMOMENTUM),cons(:,:,nmin:nmax,this%ENERGY), &
         xfluxes(:,:,nmin:nmax,this%DENSITY),xfluxes(:,:,nmin:nmax,this%XMOMENTUM), &
         xfluxes(:,:,nmin:nmax,this%YMOMENTUM),xfluxes(:,:,nmin:nmax,this%ZMOMENTUM), &
         xfluxes(:,:,nmin:nmax,this%ENERGY))
   END SUBROUTINE CalcFluxesX_euler3Drs


  PURE SUBROUTINE CalcFluxesY_euler3Drs(this,Mesh,nmin,nmax,prim,cons,yfluxes)
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
    CALL CalcFlux_euler3Drs(prim(:,:,nmin:nmax,this%DENSITY), &
         prim(:,:,nmin:nmax,this%YVELOCITY),prim(:,:,nmin:nmax,this%PRESSURE),&
         cons(:,:,nmin:nmax,this%YMOMENTUM),cons(:,:,nmin:nmax,this%XMOMENTUM), &
         cons(:,:,nmin:nmax,this%ZMOMENTUM),cons(:,:,nmin:nmax,this%ENERGY), &
         yfluxes(:,:,nmin:nmax,this%DENSITY),yfluxes(:,:,nmin:nmax,this%YMOMENTUM), &
         yfluxes(:,:,nmin:nmax,this%XMOMENTUM),yfluxes(:,:,nmin:nmax,this%ZMOMENTUM), &
         yfluxes(:,:,nmin:nmax,this%ENERGY))
  END SUBROUTINE CalcFluxesY_euler3Drs


  PURE SUBROUTINE CalcCharSystemX_euler3Drs(this,Mesh,i,dir,pvar,lambda,xvar)
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
!CDIR IEXPAND
    CALL SetEigenValues_euler3Drs(this%gamma,pvar(i,:,this%DENSITY), &
         pvar(i,:,this%XVELOCITY),pvar(i,:,this%PRESSURE), &
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4),lambda(:,5))
    ! compute characteristic variables
    i1 = i + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
    i2 = MAX(i,i1)
    i1 = MIN(i,i1)
!CDIR IEXPAND
    CALL SetCharVars_euler3Drs(this%gamma,pvar(i1,:,this%DENSITY), &
         pvar(i2,:,this%DENSITY),pvar(i1,:,this%XVELOCITY), &
         pvar(i2,:,this%XVELOCITY),pvar(i1,:,this%YVELOCITY), &
         pvar(i2,:,this%YVELOCITY),pvar(i1,:,this%ZVELOCITY), &
         pvar(i2,:,this%ZVELOCITY),pvar(i1,:,this%PRESSURE), &
         pvar(i2,:,this%PRESSURE),lambda(:,1),lambda(:,2), &
         lambda(:,3),lambda(:,4),lambda(:,5),xvar(:,1),xvar(:,2), &
         xvar(:,3),xvar(:,4),xvar(:,5))
  END SUBROUTINE CalcCharSystemX_euler3Drs
  

  PURE SUBROUTINE CalcCharSystemY_euler3Drs(this,Mesh,j,dir,pvar,lambda,xvar)
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
!CDIR IEXPAND
    CALL SetEigenValues_euler3Drs(this%gamma,pvar(:,j,this%DENSITY), &
         pvar(:,j,this%YVELOCITY),pvar(:,j,this%PRESSURE), &
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4),lambda(:,5))
    ! compute characteristic variables
    j1 = j + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
    j2 = MAX(j,j1)
    j1 = MIN(j,j1)
!CDIR IEXPAND
    CALL SetCharVars_euler3Drs(this%gamma,pvar(:,j1,this%DENSITY), &
         pvar(:,j2,this%DENSITY),pvar(:,j1,this%YVELOCITY), &
         pvar(:,j2,this%YVELOCITY),pvar(:,j1,this%XVELOCITY), &
         pvar(:,j2,this%XVELOCITY),pvar(:,j1,this%ZVELOCITY), &
         pvar(:,j2,this%ZVELOCITY),pvar(:,j1,this%PRESSURE), &
         pvar(:,j2,this%PRESSURE),lambda(:,1),lambda(:,2), &
         lambda(:,3),lambda(:,4),lambda(:,5),xvar(:,1),xvar(:,2), &
         xvar(:,3),xvar(:,4),xvar(:,5))
  END SUBROUTINE CalcCharSystemY_euler3Drs


  PURE SUBROUTINE CalcBoundaryDataX_euler3Drs(this,Mesh,i1,dir,xvar,pvar)
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
!CDIR IEXPAND
    CALL SetBoundaryData_euler3Drs(this%gamma,1.0*SIGN(1,dir),pvar(i1,:,this%DENSITY), &
         pvar(i1,:,this%XVELOCITY),pvar(i1,:,this%YVELOCITY),pvar(i1,:,this%ZVELOCITY), &
         pvar(i1,:,this%PRESSURE),xvar(:,1),xvar(:,2),xvar(:,3),xvar(:,4),xvar(:,5), &
         pvar(i2,:,this%DENSITY),pvar(i2,:,this%XVELOCITY),pvar(i2,:,this%YVELOCITY), &
         pvar(i2,:,this%ZVELOCITY),pvar(i2,:,this%PRESSURE))
  END SUBROUTINE CalcBoundaryDataX_euler3Drs


  PURE SUBROUTINE CalcBoundaryDataY_euler3Drs(this,Mesh,j1,dir,xvar,pvar)
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
!CDIR IEXPAND
    CALL SetBoundaryData_euler3Drs(this%gamma,1.0*SIGN(1,dir),pvar(:,j1,this%DENSITY), &
         pvar(:,j1,this%YVELOCITY),pvar(:,j1,this%XVELOCITY),pvar(:,j1,this%ZVELOCITY), &
         pvar(:,j1,this%PRESSURE),xvar(:,1),xvar(:,2),xvar(:,3),xvar(:,4),xvar(:,5), &
         pvar(:,j2,this%DENSITY),pvar(:,j2,this%YVELOCITY),pvar(:,j2,this%XVELOCITY), &
         pvar(:,j2,this%ZVELOCITY),pvar(:,j2,this%PRESSURE))
  END SUBROUTINE CalcBoundaryDataY_euler3Drs

 PURE SUBROUTINE CalcPrim2RiemannX_euler3Drs(this,Mesh,i,pvar,lambda,Rinv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: i
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: lambda,Rinv
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,i,pvar
    INTENT(INOUT)     :: lambda
    INTENT(OUT)       :: Rinv
    !------------------------------------------------------------------------!
    ! compute eigenvalues at i
!CDIR IEXPAND
    CALL SetEigenValues_euler3Drs(this%gamma,pvar(i,:,this%DENSITY), &
         pvar(i,:,this%XVELOCITY),pvar(i,:,this%PRESSURE),&
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4),lambda(:,5))
    ! compute Riemann invariants
!CDIR IEXPAND
    CALL Prim2Riemann_euler3Drs(this%gamma,pvar(i,:,this%DENSITY), &
         pvar(i,:,this%XVELOCITY),pvar(i,:,this%YVELOCITY), &
         pvar(i,:,this%ZVELOCITY),pvar(i,:,this%PRESSURE), &
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4),lambda(:,5), &
         Rinv(:,1),Rinv(:,2),Rinv(:,3),Rinv(:,4),Rinv(:,5))
  END SUBROUTINE CalcPrim2RiemannX_euler3Drs


  PURE SUBROUTINE CalcPrim2RiemannY_euler3Drs(this,Mesh,j,pvar,lambda,Rinv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: j
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,this%VNUM) :: lambda,Rinv
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,j,pvar
    INTENT(INOUT)     :: lambda
    INTENT(OUT)       :: Rinv
    !------------------------------------------------------------------------!
    ! compute eigenvalues at j
!CDIR IEXPAND
    CALL SetEigenValues_euler3Drs(this%gamma,pvar(:,j,this%DENSITY), &
         pvar(:,j,this%YVELOCITY),pvar(:,j,this%PRESSURE),&
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4),lambda(:,5))
    ! compute Riemann invariants
!CDIR IEXPAND
    CALL Prim2Riemann_euler3Drs(this%gamma,pvar(:,j,this%DENSITY), &
         pvar(:,j,this%YVELOCITY),pvar(:,j,this%XVELOCITY), &
         pvar(:,j,this%ZVELOCITY),pvar(:,j,this%PRESSURE), &
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4),lambda(:,5), &
         Rinv(:,1),Rinv(:,2),Rinv(:,3),Rinv(:,4),Rinv(:,5))
  END SUBROUTINE CalcPrim2RiemannY_euler3Drs

  PURE SUBROUTINE CalcRiemann2PrimX_euler3Drs(this,Mesh,i,Rinv,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: i
    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: Rinv
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    !------------------------------------------------------------------------!

    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,i,Rinv
    INTENT(INOUT)     :: pvar
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL Riemann2Prim_euler3Drs(this%gamma,Rinv(:,1),Rinv(:,2),&
       Rinv(:,3),Rinv(:,4),Rinv(:,5),pvar(i,:,this%Density), pvar(i,:,this%XVELOCITY), &
       pvar(i,:,this%YVELOCITY),pvar(i,:,this%ZVELOCITY),pvar(i,:,this%PRESSURE))
  END SUBROUTINE CalcRiemann2PrimX_euler3Drs

  PURE SUBROUTINE CalcRiemann2PrimY_euler3Drs(this,Mesh,j,Rinv,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: j
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,this%VNUM) :: Rinv
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,j,Rinv
    INTENT(INOUT)     :: pvar
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL Riemann2Prim_euler3Drs(this%gamma,Rinv(:,1),Rinv(:,2),&
       Rinv(:,3),Rinv(:,4),Rinv(:,5),pvar(:,j,this%Density),pvar(:,j,this%YVELOCITY), &
       pvar(:,j,this%XVELOCITY),pvar(:,j,this%ZVELOCITY),pvar(:,j,this%PRESSURE))
  END SUBROUTINE CalcRiemann2PrimY_euler3Drs


  PURE SUBROUTINE CalcStresses_euler3Drs(this,Mesh,pvar,dynvis,bulkvis, &
       btxx,btxy,btxz,btyy,btyz,btzz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
         dynvis,bulkvis,btxx,btxy,btxz,btyy,btyz,btzz
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,dynvis,bulkvis
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: btxx,btxy,btxz,btyy,btyz,btzz
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
               + this%tmp(i,j) ! bulk viscosity contribution
               
          btyy(i,j) = dynvis(i,j) * &
               ( (pvar(i,j+1,this%YVELOCITY) - pvar(i,j-1,this%YVELOCITY)) / Mesh%dly(i,j) &
               + 2.0 * Mesh%cyxy(i,j,1) * pvar(i,j,this%XVELOCITY) ) &
               + this%tmp(i,j) ! bulk viscosity contribution

          btzz(i,j) = dynvis(i,j) * &
               ( 2.0 * ( Mesh%czxz(i,j,1) * pvar(i,j,this%XVELOCITY) ) &
               + Mesh%czyz(i,j,1) * pvar(i,j,this%YVELOCITY) ) &
               + this%tmp(i,j) ! bulk viscosity contribution

          ! compute the off-diagonal elements (no bulk viscosity)
          btxy(i,j) = dynvis(i,j) * ( 0.5 * &
               ( (pvar(i+1,j,this%YVELOCITY) - pvar(i-1,j,this%YVELOCITY)) / Mesh%dlx(i,j) &
               + (pvar(i,j+1,this%XVELOCITY) - pvar(i,j-1,this%XVELOCITY)) / Mesh%dly(i,j) ) &
               - Mesh%cxyx(i,j,1) * pvar(i,j,this%XVELOCITY) &
               - Mesh%cyxy(i,j,1) * pvar(i,j,this%YVELOCITY) )

          btxz(i,j) = dynvis(i,j) * ( 0.5 * &
               ( (pvar(i+1,j,this%ZVELOCITY) - pvar(i-1,j,this%ZVELOCITY)) / Mesh%dlx(i,j) ) &
               - Mesh%czxz(i,j,1) * pvar(i,j,this%ZVELOCITY) )

          btyz(i,j) = dynvis(i,j) * ( 0.5 * &
               ( (pvar(i,j+1,this%ZVELOCITY) - pvar(i,j-1,this%ZVELOCITY)) / Mesh%dly(i,j) ) &
               - Mesh%czyz(i,j,1) * pvar(i,j,this%ZVELOCITY) )
       END DO
    END DO
  END SUBROUTINE CalcStresses_euler3Drs


  PURE SUBROUTINE GeometricalSources_center(this,Mesh,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
         :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,cvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! calculate centrifugal forces
          CALL CentrifugalForces_euler3Drs(pvar(i,j,this%DENSITY),pvar(i,j,this%ZVELOCITY), &
               Mesh%czxz(i,j,1),Mesh%czyz(i,j,1),this%fcent(i,j,1,1),this%fcent(i,j,2,1))
          
          ! no geometrical sources in continuity and energy equations
          sterm(i,j,this%DENSITY) = 0.0
          sterm(i,j,this%ENERGY)  = 0.0

          ! geometrical source terms in momentum equationes
          ! with centrifugal forces
          sterm(i,j,this%XMOMENTUM) = MomentumSourcesX_euler2D(cvar(i,j,this%YMOMENTUM), &
               pvar(i,j,this%XVELOCITY),pvar(i,j,this%YVELOCITY),pvar(i,j,this%PRESSURE), &
               Mesh%cxyx(i,j,1),Mesh%cyxy(i,j,1),Mesh%czxz(i,j,1)) + this%fcent(i,j,1,1)
          sterm(i,j,this%YMOMENTUM) = MomentumSourcesY_euler2D(cvar(i,j,this%XMOMENTUM), &
               pvar(i,j,this%XVELOCITY),pvar(i,j,this%YVELOCITY),pvar(i,j,this%PRESSURE), &
               Mesh%cxyx(i,j,1),Mesh%cyxy(i,j,1),Mesh%czyz(i,j,1)) + this%fcent(i,j,2,1)
          sterm(i,j,this%ZMOMENTUM) = -pvar(i,j,this%ZVELOCITY)*(Mesh%czxz(i,j,1) * &
               cvar(i,j,this%XMOMENTUM) + Mesh%czyz(i,j,1) * cvar(i,j,this%YMOMENTUM))
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
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,prim,cons
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! no geometrical sources in continuity and energy equations
          sterm(i,j,this%DENSITY) = 0.0
          sterm(i,j,this%ENERGY)  = 0.0
!CDIR UNROLL=4
          DO k=1,4
            ! calculate centrifugal forces
            CALL CentrifugalForces_euler3Drs(prim(i,j,k,this%DENSITY), &
                  prim(i,j,k,this%ZVELOCITY),Mesh%czxz(i,j,k),Mesh%czyz(i,j,k), &
                  this%fcent(i,j,k,1),this%fcent(i,j,k,2))
            ! add momentum sources to centrifugal forces
            this%fcent(i,j,k,1) = this%fcent(i,j,k,1) + MomentumSourcesX_euler2D( &
                  cons(i,j,k,this%YMOMENTUM),prim(i,j,k,this%XVELOCITY), &
                  prim(i,j,k,this%YVELOCITY),prim(i,j,k,this%PRESSURE), &
                  Mesh%cxyx(i,j,k),Mesh%cyxy(i,j,k),Mesh%czxz(i,j,k))

            this%fcent(i,j,k,2) = this%fcent(i,j,k,2) + MomentumSourcesY_euler2D( &
                  cons(i,j,k,this%XMOMENTUM),prim(i,j,k,this%XVELOCITY), &
                  prim(i,j,k,this%YVELOCITY),prim(i,j,k,this%PRESSURE), &
                  Mesh%cxyx(i,j,k),Mesh%cyxy(i,j,k),Mesh%czyz(i,j,k))


          END DO

          ! momentum sources (sum up corner values, don't use SUM function,
          ! because it prevents COLLAPSING and causes poor vectorization
          sterm(i,j,this%XMOMENTUM) = this%fcent(i,j,1,1)+this%fcent(i,j,2,1) &
                  +this%fcent(i,j,3,1)+this%fcent(i,j,4,1)
          sterm(i,j,this%YMOMENTUM) = this%fcent(i,j,1,2)+this%fcent(i,j,2,2) &
                  +this%fcent(i,j,3,2)+this%fcent(i,j,4,2)
          sterm(i,j,this%ZMOMENTUM) = -prim(i,j,1,this%ZVELOCITY) &
                  * (Mesh%czxz(i,j,1)*cons(i,j,1,this%XMOMENTUM) &
                    +Mesh%czyz(i,j,1)*cons(i,j,1,this%YMOMENTUM)) &
                - prim(i,j,2,this%ZVELOCITY) &
                  * (Mesh%czxz(i,j,2)*cons(i,j,2,this%XMOMENTUM) &
                    +Mesh%czyz(i,j,2)*cons(i,j,2,this%YMOMENTUM)) &
                - prim(i,j,3,this%ZVELOCITY) &
                  * (Mesh%czxz(i,j,3)*cons(i,j,3,this%XMOMENTUM) &
                    +Mesh%czyz(i,j,3)*cons(i,j,3,this%YMOMENTUM)) &
                - prim(i,j,4,this%ZVELOCITY) &
                  * (Mesh%czxz(i,j,4)*cons(i,j,4,this%XMOMENTUM) &
                    +Mesh%czyz(i,j,4)*cons(i,j,4,this%YMOMENTUM))
       END DO
    END DO
  END SUBROUTINE GeometricalSources_faces


  PURE SUBROUTINE ExternalSources_euler3Drs(this,Mesh,accel,pvar,cvar,sterm)
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,3) &
         :: accel
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
         :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,accel,pvar,cvar
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    sterm(:,:,this%DENSITY)   = 0.
    sterm(:,:,this%XMOMENTUM) = pvar(:,:,this%DENSITY) * accel(:,:,1)
    sterm(:,:,this%YMOMENTUM) = pvar(:,:,this%DENSITY) * accel(:,:,2)
    sterm(:,:,this%ZMOMENTUM) = pvar(:,:,this%DENSITY) * accel(:,:,3)
    sterm(:,:,this%ENERGY)    = cvar(:,:,this%XMOMENTUM) * accel(:,:,1) &
         +cvar(:,:,this%YMOMENTUM) * accel(:,:,2) &
         +cvar(:,:,this%ZMOMENTUM) * accel(:,:,3)
  END SUBROUTINE ExternalSources_euler3Drs


  PURE SUBROUTINE ViscositySources_euler3Drs(this,Mesh,pvar,btxx,btxy,btxz, &
       btyy,btyz,btzz,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: &
         pvar,sterm
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
         btxx,btxy,btxz,btyy,btyz,btzz
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar, btxx,btxy,btxz,btyy,btyz,btzz
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! viscosity source terms
    sterm(:,:,this%DENSITY) = 0.0 
 
    ! compute 3D tensor divergence for symmetric tensor btxy=btyx, btxz=btzx, etc.
!CDIR IEXPAND
    CALL Divergence(Mesh,btxx,btxy,btxz,btxy,btyy,btyz,btxz,btyz,btzz, &
                    sterm(:,:,this%XMOMENTUM),sterm(:,:,this%YMOMENTUM), &
                    sterm(:,:,this%ZMOMENTUM))
 
    ! compute scalar product of v and btx (x-component)
    ! use this%tmin for temporary storage
!CDIR NODEP
    this%tmin(:,:,1)  = pvar(:,:,this%XVELOCITY)*btxx(:,:) &
                      + pvar(:,:,this%YVELOCITY)*btxy(:,:) &
                      + pvar(:,:,this%ZVELOCITY)*btxz(:,:)

   ! compute scalar product of v and bty (y-component)
   ! use amax for temporary storage
!CDIR NODEP 
    this%tmin(:,:,2) = pvar(:,:,this%XVELOCITY)*btxy(:,:) &
                     + pvar(:,:,this%YVELOCITY)*btyy(:,:) &
                     + pvar(:,:,this%ZVELOCITY)*btyz(:,:) 
 
    ! compute vector divergence of scalar product v_i * bt_ij
!CDIR IEXPAND
    CALL Divergence(Mesh,this%tmin(:,:,1),this%tmin(:,:,2),sterm(:,:,this%ENERGY))
  END SUBROUTINE ViscositySources_euler3Drs


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
    CALL Cons2Prim_euler3Drs(this%gamma,cvar(i1:i2,j1:j2,this%DENSITY), &
         cvar(i1:i2,j1:j2,this%XMOMENTUM),cvar(i1:i2,j1:j2,this%YMOMENTUM), &
         cvar(i1:i2,j1:j2,this%ZMOMENTUM),cvar(i1:i2,j1:j2,this%ENERGY), &
         pvar(i1:i2,j1:j2,this%DENSITY),pvar(i1:i2,j1:j2,this%XVELOCITY), &
         pvar(i1:i2,j1:j2,this%YVELOCITY),pvar(i1:i2,j1:j2,this%ZVELOCITY), &
         pvar(i1:i2,j1:j2,this%PRESSURE))
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
    CALL Cons2Prim_euler3Drs(this%gamma,cons(i1:i2,j1:j2,:,this%DENSITY), &
         cons(i1:i2,j1:j2,:,this%XMOMENTUM),cons(i1:i2,j1:j2,:,this%YMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%ZMOMENTUM),cons(i1:i2,j1:j2,:,this%ENERGY), &
         prim(i1:i2,j1:j2,:,this%DENSITY),prim(i1:i2,j1:j2,:,this%XVELOCITY), &
         prim(i1:i2,j1:j2,:,this%YVELOCITY),prim(i1:i2,j1:j2,:,this%ZVELOCITY), &
         prim(i1:i2,j1:j2,:,this%PRESSURE))
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
    CALL Prim2Cons_euler3Drs(this%gamma,pvar(i1:i2,j1:j2,this%DENSITY), &
         pvar(i1:i2,j1:j2,this%XVELOCITY),pvar(i1:i2,j1:j2,this%YVELOCITY), &
         pvar(i1:i2,j1:j2,this%ZVELOCITY),pvar(i1:i2,j1:j2,this%PRESSURE), &
         cvar(i1:i2,j1:j2,this%DENSITY),cvar(i1:i2,j1:j2,this%XMOMENTUM), &
         cvar(i1:i2,j1:j2,this%YMOMENTUM),cvar(i1:i2,j1:j2,this%ZMOMENTUM), &
         cvar(i1:i2,j1:j2,this%ENERGY))
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
    CALL Prim2Cons_euler3Drs(this%gamma,prim(i1:i2,j1:j2,:,this%DENSITY), &
         prim(i1:i2,j1:j2,:,this%XVELOCITY),prim(i1:i2,j1:j2,:,this%YVELOCITY), &
         prim(i1:i2,j1:j2,:,this%ZVELOCITY),prim(i1:i2,j1:j2,:,this%PRESSURE), &
         cons(i1:i2,j1:j2,:,this%DENSITY),cons(i1:i2,j1:j2,:,this%XMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%YMOMENTUM),cons(i1:i2,j1:j2,:,this%ZMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%ENERGY))
  END SUBROUTINE Convert2Conservative_faces


  PURE SUBROUTINE ReflectionMasks_euler3Drs(this,reflX,reflY)
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
    reflX(this%ZVELOCITY) = .FALSE.
    reflX(this%PRESSURE)  = .FALSE.
    ! southern / northern boundary
    reflY(this%DENSITY)   = .FALSE.
    reflY(this%XVELOCITY) = .FALSE.
    reflY(this%YVELOCITY) = .TRUE.
    reflY(this%ZVELOCITY) = .FALSE.
    reflY(this%PRESSURE)  = .FALSE.
  END SUBROUTINE ReflectionMasks_euler3Drs


  PURE SUBROUTINE AxisMasks_euler3Drs(this,reflX,reflY)
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
    reflX(this%ZVELOCITY) = .TRUE.
    reflX(this%PRESSURE)  = .FALSE.
    ! southern / northern boundary
    reflY(this%DENSITY)   = .FALSE.
    reflY(this%XVELOCITY) = .FALSE.
    reflY(this%YVELOCITY) = .TRUE.
    reflY(this%ZVELOCITY) = .TRUE.
    reflY(this%PRESSURE)  = .FALSE.
  END SUBROUTINE AxisMasks_euler3Drs


  ELEMENTAL SUBROUTINE SetEigenValues_euler3Drs(gamma,rho,v,P,l1,l2,l3,l4,l5)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho,v,P
    REAL, INTENT(OUT) :: l1,l2,l3,l4,l5
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL SetEigenValues_euler2D(gamma,rho,v,P,l1,l2,l3,l5)
    ! set the missing eigenvalue l4
    l4 = v
  END SUBROUTINE SetEigenValues_euler3Drs


  ELEMENTAL SUBROUTINE SetCharVars_euler3Drs(gamma,rho1,rho2,u1,u2,v1,v2,w1,w2, &
       P1,P2,l1,l2,l3,l4,l5,xvar1,xvar2,xvar3,xvar4,xvar5)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho1,rho2,u1,u2,v1,v2,w1,w2,P1,P2,l1,l2,l3,l4,l5
    REAL, INTENT(OUT) :: xvar1,xvar2,xvar3,xvar4,xvar5
    !------------------------------------------------------------------------!
    REAL :: gamcs,dlnP,du
    !------------------------------------------------------------------------!
    gamcs= 2.*gamma / (l5-l1) ! = gamma/cs
    dlnP = LOG(P2/P1)         ! = LOG(P2)-LOG(P1)
    du   = u2-u1
    ! characteristic variables
    xvar1 = dlnP - gamcs * du
    xvar2 = -dlnP + gamma * LOG(rho2/rho1)
    xvar3 = gamcs * (v2-v1)
    xvar4 = gamcs * (w2-w1)
    xvar5 = dlnP + gamcs * du 
  END SUBROUTINE SetCharVars_euler3Drs

  ELEMENTAL SUBROUTINE SetBoundaryData_euler3Drs(gamma,dir,rho1,u1,v1,w1,P1, &
       xvar1,xvar2,xvar3,xvar4,xvar5,rho2,u2,v2,w2,P2)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,dir,rho1,u1,v1,w1,P1,xvar1,xvar2,xvar3,xvar4,xvar5
    REAL, INTENT(OUT) :: rho2,u2,v2,w2,P2
    !------------------------------------------------------------------------!
    REAL :: dlnP,csgam
    !------------------------------------------------------------------------!
    dlnP = 0.5 * (xvar5+xvar1)
    ! extrapolate boundary values using characteristic variables
    rho2 = rho1 * EXP(dir*(xvar2+dlnP)/gamma)
    P2   = P1 * EXP(dir*dlnP)
!CDIR IEXPAND
    csgam= GetSoundSpeed_euler2D(gamma,rho1+rho2,P1+P2) / gamma
    u2   = u1 + dir*csgam * 0.5*(xvar5-xvar1)
    v2   = v1 + dir*csgam * xvar3
    w2   = w1 + dir*csgam * xvar4
  END SUBROUTINE SetBoundaryData_euler3Drs

 ELEMENTAL SUBROUTINE Prim2Riemann_euler3Drs(gamma,rho,vx,vy,vz,p,&
                                       l1,l2,l3,l4,l5,Rminus,Rs,Rvt,Rvt2,Rplus)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho,vx,vy,vz,p,l1,l2,l3,l4,l5
    REAL, INTENT(OUT) :: Rminus,Rs,Rvt,Rvt2,Rplus
    !------------------------------------------------------------------------!
    REAL :: cs
    !------------------------------------------------------------------------!
    cs = l5-l2 ! l2 = v, l6 = v+cs
    ! compute 1st Riemann invariant (R+)
    Rplus = vx + 2./(gamma-1.0) * cs     
    ! compute 2st Riemann invariant (R-) 
    Rminus = vx - 2./(gamma-1.0) * cs
    ! compute entropy
    Rs = p/rho**gamma
    ! tangential velocities
    Rvt = vy   
    Rvt2 = vz     
  END SUBROUTINE Prim2Riemann_euler3Drs

  ELEMENTAL SUBROUTINE Riemann2Prim_euler3Drs(gamma,Rminus,Rs,Rvt,Rvt2,Rplus,&
       rho,vx,vy,vz,p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,Rminus,Rs,Rvt,Rvt2,Rplus
    REAL, INTENT(OUT) :: rho,vx,vy,vz,p
    !------------------------------------------------------------------------!
    REAL :: cs2gam
    !------------------------------------------------------------------------!
    ! tangential velocity    
    vy = Rvt  
    vz = Rvt2
    ! normal velocity
    vx = 0.5*(Rplus+Rminus)
    ! cs**2 / gamma
    cs2gam = (0.25*(gamma-1.0)*(Rplus-Rminus))**2 / gamma
    ! density
    rho = (cs2gam/Rs)**(1./(gamma-1.0))
    ! pressure
    p = cs2gam * rho
  END SUBROUTINE Riemann2Prim_euler3Drs


  ELEMENTAL SUBROUTINE CalcFlux_euler3Drs(rho,v,P,m1,m2,m3,E,f1,f2,f3,f4,f5)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho,v,P,m1,m2,m3,E
    REAL, INTENT(OUT) :: f1, f2, f3, f4, f5
    !------------------------------------------------------------------------!
    f1 = rho*v
    f2 = m1*v + P
    f3 = m2*v
    f4 = m3*v
    f5 = (E+P)*v
  END SUBROUTINE CalcFlux_euler3Drs


  ELEMENTAL SUBROUTINE CentrifugalForces_euler3Drs(rho,vz,czxz,czyz,fcx,fcy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho,vz,czxz,czyz
    REAL, INTENT(OUT) :: fcx,fcy
    !------------------------------------------------------------------------!
    fcx  = rho*vz*vz ! temporary 
    fcy  = fcx * czyz
    fcx  = fcx * czxz
  END SUBROUTINE CentrifugalForces_euler3Drs


  ELEMENTAL SUBROUTINE Cons2Prim_euler3Drs(gamma,rho_in,mu,mv,mw,E,rho_out,u,v,w,P)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho_in,mu,mv,mw,E
    REAL, INTENT(OUT) :: rho_out,u,v,w,P
    !------------------------------------------------------------------------!
    REAL :: inv_rho
    !------------------------------------------------------------------------!
    inv_rho = 1./rho_in
    rho_out = rho_in
    u = mu * inv_rho
    v = mv * inv_rho
    w = mw * inv_rho
    P = (gamma-1.)*(E - 0.5 * inv_rho * (mu*mu+mv*mv+mw*mw))
  END SUBROUTINE Cons2Prim_euler3Drs

  
  ELEMENTAL SUBROUTINE Prim2Cons_euler3Drs(gamma,rho_in,u,v,w,P,rho_out,mu,mv,mw,E)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho_in,u,v,w,P
    REAL, INTENT(OUT) :: rho_out,mu,mv,mw,E
    !------------------------------------------------------------------------!
    rho_out = rho_in
    mu = rho_in * u
    mv = rho_in * v
    mw = rho_in * w
    E  = P/(gamma-1.) + 0.5 * rho_in * (u*u+v*v+w*w)
  END SUBROUTINE Prim2Cons_euler3Drs


  SUBROUTINE ClosePhysics_euler3Drs(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL ClosePhysics_euler2D(this)
    DEALLOCATE(this%fcent)
  END SUBROUTINE ClosePhysics_euler3Drs

END MODULE physics_euler3Drotsym
