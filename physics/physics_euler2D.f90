!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler2D.f90                                               #
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
!! - non-isothermal gas dynamics
!!   \key{gamma,REAL,ratio of specific heats (default is for diatomic
!!      molecular gas),1.4}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Björn Sperling
!!
!! \brief basic module for 2D Euler equations
!!
!! \extends physics_common
!! \ingroup physics
!----------------------------------------------------------------------------!
MODULE physics_euler2D
  USE physics_common
  USE mesh_common, ONLY : Mesh_TYP
  USE sources_common, ONLY : Sources_TYP
  USE physics_euler2Disothm, ONLY : &
       CalcStresses_euler2D => CalcStresses_euler2Dit, &
       MomentumSourcesX_euler2D => MomentumSourcesX_euler2Dit, &
       MomentumSourcesY_euler2D => MomentumSourcesY_euler2Dit, &
       CalcWaveSpeeds_euler2D => CalcWaveSpeeds_euler2Dit, &
       SetEigenValues_euler2Dit, &
       ViscositySources_euler2Dit, &
       ClosePhysics_euler2Dit

  USE mesh_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE CalcSoundSpeeds_euler2D
     MODULE PROCEDURE CalcSoundSpeeds_center
     MODULE PROCEDURE CalcSoundSpeeds_faces
  END INTERFACE
  INTERFACE GeometricalSources_euler2D
     MODULE PROCEDURE GeometricalSources_center
     MODULE PROCEDURE GeometricalSources_faces
  END INTERFACE
  INTERFACE Convert2Primitive_euler2D
     MODULE PROCEDURE Convert2Primitive_center
     MODULE PROCEDURE Convert2Primitive_faces
  END INTERFACE
  INTERFACE Convert2Conservative_euler2D
     MODULE PROCEDURE Convert2Conservative_center
     MODULE PROCEDURE Convert2Conservative_faces
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: num_var = 4              ! number of variables       !
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 2D"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Physics_TYP, &
       ! methods
       InitPhysics_euler2D, &
       ClosePhysics_euler2D, &
       CalcWaveSpeeds_euler2D, &
       CalcSoundSpeeds_euler2D, &
       CalcFluxesX_euler2D, &
       CalcFluxesY_euler2D, &
       CalcCharSystemX_euler2D, &
       CalcCharSystemY_euler2D, &
       CalcBoundaryDataX_euler2D, &
       CalcBoundaryDataY_euler2D, &
       CalcPrim2RiemannX_euler2D, &
       CalcPrim2RiemannY_euler2D, &
       CalcRiemann2PrimX_euler2D, &
       CalcRiemann2PrimY_euler2D, &
       CalcStresses_euler2D, &
       GeometricalSources_euler2D, &
       ViscositySources_euler2D, &
       ExternalSources_euler2D, &
       Convert2Primitive_euler2D, &
       Convert2Conservative_euler2D, &
       ReflectionMasks_euler2D, &
       AxisMasks_euler2D, &
       GetSoundSpeed_euler2D, &
       SetEigenValues_euler2D, &
       SetWaveSpeeds_euler2D, &
       MomentumSourcesX_euler2D, &
       MomentumSourcesY_euler2D
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitPhysics_euler2D(this,Mesh,problem,pname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: problem
    CHARACTER(LEN=32), OPTIONAL :: pname
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,problem,pname
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    IF (PRESENT(pname)) THEN
       CALL InitPhysics(this,problem,pname,num_var)
    ELSE
       CALL InitPhysics(this,problem,problem_name,num_var)
    END IF
    ! set array indices
    this%DENSITY   = 1                                 ! mass density        !
    this%PRESSURE  = num_var                           ! pressure            !
    this%ENERGY    = num_var                           ! total energy        !
    this%XVELOCITY = 2                                 ! x-velocity          !
    this%XMOMENTUM = 2                                 ! x-momentum          !
    this%YVELOCITY = 3                                 ! y-velocity          !
    this%YMOMENTUM = 3                                 ! y-momentum          !
    this%ZVELOCITY = 0                                 ! no z-velocity       !
    this%ZMOMENTUM = 0                                 ! no z-momentum       !
    ! set names for primitive and conservative variables
    this%pvarname(this%DENSITY)   = "density"
    this%pvarname(this%XVELOCITY) = "xvelocity"
    this%pvarname(this%YVELOCITY) = "yvelocity"
    this%pvarname(this%PRESSURE)  = "pressure"
    this%cvarname(this%DENSITY)   = "density"
    this%cvarname(this%XMOMENTUM) = "xmomentum"
    this%cvarname(this%YMOMENTUM) = "ymomentum"
    this%cvarname(this%ENERGY)    = "energy"
    this%DIM = 2

  END SUBROUTINE InitPhysics_euler2D


  PURE SUBROUTINE CalcFluxesX_euler2D(this,Mesh,nmin,nmax,prim,cons,xfluxes)
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
    CALL CalcFlux_euler2D(prim(:,:,nmin:nmax,this%DENSITY), &
         prim(:,:,nmin:nmax,this%XVELOCITY),prim(:,:,nmin:nmax,this%PRESSURE), &
         cons(:,:,nmin:nmax,this%XMOMENTUM),cons(:,:,nmin:nmax,this%YMOMENTUM), &
         cons(:,:,nmin:nmax,this%ENERGY),xfluxes(:,:,nmin:nmax,this%DENSITY),&
         xfluxes(:,:,nmin:nmax,this%XMOMENTUM),xfluxes(:,:,nmin:nmax,this%YMOMENTUM), &
         xfluxes(:,:,nmin:nmax,this%ENERGY))
  END SUBROUTINE CalcFluxesX_euler2D


  PURE SUBROUTINE CalcFluxesY_euler2D(this,Mesh,nmin,nmax,prim,cons,yfluxes)
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
    CALL CalcFlux_euler2D(prim(:,:,nmin:nmax,this%DENSITY), &
         prim(:,:,nmin:nmax,this%YVELOCITY),prim(:,:,nmin:nmax,this%PRESSURE), &
         cons(:,:,nmin:nmax,this%YMOMENTUM),cons(:,:,nmin:nmax,this%XMOMENTUM), &
         cons(:,:,nmin:nmax,this%ENERGY),yfluxes(:,:,nmin:nmax,this%DENSITY), &
         yfluxes(:,:,nmin:nmax,this%YMOMENTUM),yfluxes(:,:,nmin:nmax,this%XMOMENTUM), &
         yfluxes(:,:,nmin:nmax,this%ENERGY))
  END SUBROUTINE CalcFluxesY_euler2D


  PURE SUBROUTINE CalcCharSystemX_euler2D(this,Mesh,i,dir,pvar,lambda,xvar)
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
    CALL SetEigenValues_euler2D(this%gamma,pvar(i,:,this%DENSITY), &
         pvar(i,:,this%XVELOCITY), pvar(i,:,this%PRESSURE),&
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4))
    ! compute characteristic variables
    i1 = i + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
    i2 = MAX(i,i1)
    i1 = MIN(i,i1)
!CDIR IEXPAND
    CALL SetCharVars_euler2D(this%gamma,pvar(i1,:,this%DENSITY), &
         pvar(i2,:,this%DENSITY),pvar(i1,:,this%XVELOCITY), &
         pvar(i2,:,this%XVELOCITY),pvar(i1,:,this%YVELOCITY), &
         pvar(i2,:,this%YVELOCITY),pvar(i1,:,this%PRESSURE), &
         pvar(i2,:,this%PRESSURE),lambda(:,1),lambda(:,2),lambda(:,3), &
         lambda(:,4),xvar(:,1),xvar(:,2),xvar(:,3),xvar(:,4))
  END SUBROUTINE CalcCharSystemX_euler2D


  PURE SUBROUTINE CalcCharSystemY_euler2D(this,Mesh,j,dir,pvar,lambda,xvar)
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
    CALL SetEigenValues_euler2D(this%gamma,pvar(:,j,this%DENSITY), &
         pvar(:,j,this%YVELOCITY), pvar(:,j,this%PRESSURE), &
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4))
    ! compute characteristic variables
    j1 = j + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
    j2 = MAX(j,j1)
    j1 = MIN(j,j1)
!CDIR IEXPAND
   CALL SetCharVars_euler2D(this%gamma,pvar(:,j1,this%DENSITY), &
         pvar(:,j2,this%DENSITY),pvar(:,j1,this%YVELOCITY), &
         pvar(:,j2,this%YVELOCITY),pvar(:,j1,this%XVELOCITY), &
         pvar(:,j2,this%XVELOCITY),pvar(:,j1,this%PRESSURE), &
         pvar(:,j2,this%PRESSURE),lambda(:,1),lambda(:,2),lambda(:,3), &
         lambda(:,4),xvar(:,1),xvar(:,2),xvar(:,3),xvar(:,4))
  END SUBROUTINE CalcCharSystemY_euler2D


  PURE SUBROUTINE CalcBoundaryDataX_euler2D(this,Mesh,i1,dir,xvar,pvar)
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
    CALL SetBoundaryData_euler2D(this%gamma,1.0*SIGN(1,dir),pvar(i1,:,this%DENSITY), &
         pvar(i1,:,this%XVELOCITY),pvar(i1,:,this%YVELOCITY),pvar(i1,:,this%PRESSURE), &
         xvar(:,1),xvar(:,2),xvar(:,3),xvar(:,4),pvar(i2,:,this%DENSITY), &
         pvar(i2,:,this%XVELOCITY),pvar(i2,:,this%YVELOCITY),pvar(i2,:,this%PRESSURE))
  END SUBROUTINE CalcBoundaryDataX_euler2D


  PURE SUBROUTINE CalcBoundaryDataY_euler2D(this,Mesh,j1,dir,xvar,pvar)
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
    CALL SetBoundaryData_euler2D(this%gamma,1.0*SIGN(1,dir),pvar(:,j1,this%DENSITY), &
         pvar(:,j1,this%YVELOCITY),pvar(:,j1,this%XVELOCITY),pvar(:,j1,this%PRESSURE), &
         xvar(:,1),xvar(:,2),xvar(:,3),xvar(:,4),pvar(:,j2,this%DENSITY), &
         pvar(:,j2,this%YVELOCITY),pvar(:,j2,this%XVELOCITY),pvar(:,j2,this%PRESSURE))
  END SUBROUTINE CalcBoundaryDataY_euler2D

  PURE SUBROUTINE CalcPrim2RiemannX_euler2D(this,Mesh,i,pvar,lambda,Rinv)
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
    CALL SetEigenValues_euler2D(this%gamma,pvar(i,:,this%DENSITY), &
         pvar(i,:,this%XVELOCITY),pvar(i,:,this%PRESSURE),&
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4))
    ! compute Riemann invariants
!CDIR IEXPAND
    CALL Prim2Riemann_euler2D(this%gamma,pvar(i,:,this%DENSITY), &
         pvar(i,:,this%XVELOCITY),pvar(i,:,this%YVELOCITY), &
         pvar(i,:,this%PRESSURE),lambda(:,1),lambda(:,2),lambda(:,3),&
         lambda(:,4),Rinv(:,1),Rinv(:,2),Rinv(:,3),Rinv(:,4))
  END SUBROUTINE CalcPrim2RiemannX_euler2D


  PURE SUBROUTINE CalcPrim2RiemannY_euler2D(this,Mesh,j,pvar,lambda,Rinv)
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
    CALL SetEigenValues_euler2D(this%gamma,pvar(:,j,this%DENSITY), &
         pvar(:,j,this%YVELOCITY),pvar(:,j,this%PRESSURE),&
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4))
    ! compute Riemann invariants
!CDIR IEXPAND
    CALL Prim2Riemann_euler2D(this%gamma,pvar(:,j,this%DENSITY), &
         pvar(:,j,this%YVELOCITY),pvar(:,j,this%XVELOCITY), &
         pvar(:,j,this%PRESSURE),lambda(:,1),lambda(:,2),lambda(:,3),&
         lambda(:,4),Rinv(:,1),Rinv(:,2),Rinv(:,3),Rinv(:,4))
  END SUBROUTINE CalcPrim2RiemannY_euler2D

  PURE SUBROUTINE CalcRiemann2PrimX_euler2D(this,Mesh,i,Rinv,pvar)
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
    CALL Riemann2Prim_euler2D(this%gamma,Rinv(:,1),Rinv(:,2),&
       Rinv(:,3),Rinv(:,4),pvar(i,:,this%Density), pvar(i,:,this%XVELOCITY), &
       pvar(i,:,this%YVELOCITY), pvar(i,:,this%PRESSURE))
  END SUBROUTINE CalcRiemann2PrimX_euler2D

  PURE SUBROUTINE CalcRiemann2PrimY_euler2D(this,Mesh,j,Rinv,pvar)
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
    CALL Riemann2Prim_euler2D(this%gamma,Rinv(:,1),Rinv(:,2),&
       Rinv(:,3),Rinv(:,4),pvar(:,j,this%Density), pvar(:,j,this%YVELOCITY), &
       pvar(:,j,this%XVELOCITY), pvar(:,j,this%PRESSURE))
  END SUBROUTINE CalcRiemann2PrimY_euler2D

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
          ! no geometrical density or energy sources
          sterm(i,j,this%DENSITY) = 0.
          sterm(i,j,this%ENERGY) = 0.
          ! geometrical source terms in momentum equationes
          sterm(i,j,this%XMOMENTUM) = MomentumSourcesX_euler2D(&
              cvar(i,j,this%YMOMENTUM),pvar(i,j,this%XVELOCITY),&
              pvar(i,j,this%YVELOCITY),pvar(i,j,this%PRESSURE), &
              Mesh%cxyx(i,j,1),Mesh%cyxy(i,j,1),Mesh%czxz(i,j,1))
          sterm(i,j,this%YMOMENTUM) = MomentumSourcesY_euler2D(&
              cvar(i,j,this%XMOMENTUM),pvar(i,j,this%XVELOCITY),&
              pvar(i,j,this%YVELOCITY),pvar(i,j,this%PRESSURE), &
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
          ! no geometrical density or energy sources
          sterm(i,j,this%DENSITY)   = 0.
          sterm(i,j,this%ENERGY)    = 0.
          ! momentum sources (sum up corner values, don't use SUM function,
          ! because it prevents COLLAPSING and causes poor vectorization
          sterm(i,j,this%XMOMENTUM) = MomentumSourcesX_euler2D(&
              cons(i,j,1,this%YMOMENTUM),prim(i,j,1,this%XVELOCITY),&
              prim(i,j,1,this%YVELOCITY),prim(i,j,1,this%PRESSURE), &
              Mesh%cxyx(i,j,1),Mesh%cyxy(i,j,1),Mesh%czxz(i,j,1)) &
            + MomentumSourcesX_euler2D(&
              cons(i,j,2,this%YMOMENTUM),prim(i,j,2,this%XVELOCITY),&
              prim(i,j,2,this%YVELOCITY),prim(i,j,2,this%PRESSURE), &
              Mesh%cxyx(i,j,2),Mesh%cyxy(i,j,2),Mesh%czxz(i,j,2)) &
            + MomentumSourcesX_euler2D(&
              cons(i,j,3,this%YMOMENTUM),prim(i,j,3,this%XVELOCITY),&
              prim(i,j,3,this%YVELOCITY),prim(i,j,3,this%PRESSURE), &
              Mesh%cxyx(i,j,3),Mesh%cyxy(i,j,3),Mesh%czxz(i,j,3)) &
            + MomentumSourcesX_euler2D(&
              cons(i,j,4,this%YMOMENTUM),prim(i,j,4,this%XVELOCITY),&
              prim(i,j,4,this%YVELOCITY),prim(i,j,4,this%PRESSURE), &
              Mesh%cxyx(i,j,4),Mesh%cyxy(i,j,4),Mesh%czxz(i,j,4)) 

          sterm(i,j,this%YMOMENTUM) = MomentumSourcesY_euler2D(&
              cons(i,j,1,this%XMOMENTUM),prim(i,j,1,this%XVELOCITY),&
              prim(i,j,1,this%YVELOCITY),prim(i,j,1,this%PRESSURE), &
              Mesh%cxyx(i,j,1),Mesh%cyxy(i,j,1),Mesh%czyz(i,j,1)) &
            + MomentumSourcesY_euler2D(&
              cons(i,j,2,this%XMOMENTUM),prim(i,j,2,this%XVELOCITY),&
              prim(i,j,2,this%YVELOCITY),prim(i,j,2,this%PRESSURE), &
              Mesh%cxyx(i,j,2),Mesh%cyxy(i,j,2),Mesh%czyz(i,j,2)) &
            + MomentumSourcesY_euler2D(&
              cons(i,j,3,this%XMOMENTUM),prim(i,j,3,this%XVELOCITY),&
              prim(i,j,3,this%YVELOCITY),prim(i,j,3,this%PRESSURE), &
              Mesh%cxyx(i,j,3),Mesh%cyxy(i,j,3),Mesh%czyz(i,j,3)) &
            + MomentumSourcesY_euler2D(&
              cons(i,j,4,this%XMOMENTUM),prim(i,j,4,this%XVELOCITY),&
              prim(i,j,4,this%YVELOCITY),prim(i,j,4,this%PRESSURE), &
              Mesh%cxyx(i,j,4),Mesh%cyxy(i,j,4),Mesh%czyz(i,j,4))
       END DO
    END DO 
  END SUBROUTINE GeometricalSources_faces


  ! momentum and energy sources due to external force
  PURE SUBROUTINE ExternalSources_euler2D(this,Mesh,accel,pvar,cvar,sterm)
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
          sterm(i,j,this%ENERGY)    = cvar(i,j,this%XMOMENTUM) * accel(i,j,1) + &
               cvar(i,j,this%YMOMENTUM) * accel(i,j,2)
       END DO
    END DO
  END SUBROUTINE ExternalSources_euler2D


  PURE SUBROUTINE ViscositySources_euler2D(this,Mesh,pvar,btxx,btxy,btyy,sterm)         
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: &
         pvar,sterm
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
         btxx,btxy,btyy
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,btxx,btxy,btyy
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL ViscositySources_euler2Dit(this,Mesh,pvar,btxx,btxy,btyy,sterm)
 
    !compute scalar product of v and tau (x-component)
    this%amin(:,:)  = pvar(:,:,this%XVELOCITY)*btxx(:,:) &
                    + pvar(:,:,this%YVELOCITY)*btxy(:,:) 

    !compute scalar product of v and tau (y-component)
    this%amax(:,:) = pvar(:,:,this%XVELOCITY)*btxy(:,:) &
                   + pvar(:,:,this%YVELOCITY)*btyy(:,:)

    ! compute vector divergence of scalar product v and tau
!CDIR IEXPAND
    CALL Divergence(Mesh,this%amin(:,:),this%amax(:,:),sterm(:,:,this%ENERGY))
  END SUBROUTINE ViscositySources_euler2D


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
    CALL Cons2Prim_euler2D(this%gamma,cvar(i1:i2,j1:j2,this%DENSITY), &
         cvar(i1:i2,j1:j2,this%XMOMENTUM),cvar(i1:i2,j1:j2,this%YMOMENTUM), &
         cvar(i1:i2,j1:j2,this%ENERGY),pvar(i1:i2,j1:j2,this%DENSITY), &
         pvar(i1:i2,j1:j2,this%XVELOCITY),pvar(i1:i2,j1:j2,this%YVELOCITY), &
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
    CALL Cons2Prim_euler2D(this%gamma,cons(i1:i2,j1:j2,:,this%DENSITY), &
         cons(i1:i2,j1:j2,:,this%XMOMENTUM),cons(i1:i2,j1:j2,:,this%YMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%ENERGY),prim(i1:i2,j1:j2,:,this%DENSITY), &
         prim(i1:i2,j1:j2,:,this%XVELOCITY),prim(i1:i2,j1:j2,:,this%YVELOCITY), &
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
    CALL Prim2Cons_euler2D(this%gamma,pvar(i1:i2,j1:j2,this%DENSITY), &
         pvar(i1:i2,j1:j2,this%XVELOCITY),pvar(i1:i2,j1:j2,this%YVELOCITY), &
         pvar(i1:i2,j1:j2,this%PRESSURE),cvar(i1:i2,j1:j2,this%DENSITY), &
         cvar(i1:i2,j1:j2,this%XMOMENTUM),cvar(i1:i2,j1:j2,this%YMOMENTUM), &
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
    CALL Prim2Cons_euler2D(this%gamma,prim(i1:i2,j1:j2,:,this%DENSITY), &
         prim(i1:i2,j1:j2,:,this%XVELOCITY),prim(i1:i2,j1:j2,:,this%YVELOCITY), &
         prim(i1:i2,j1:j2,:,this%PRESSURE),cons(i1:i2,j1:j2,:,this%DENSITY), &
         cons(i1:i2,j1:j2,:,this%XMOMENTUM),cons(i1:i2,j1:j2,:,this%YMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%ENERGY))
  END SUBROUTINE Convert2Conservative_faces


  PURE SUBROUTINE ReflectionMasks_euler2D(this,reflX,reflY)
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
    reflX(this%PRESSURE)  = .FALSE.
    ! southern / northern boundary
    reflY(this%DENSITY)   = .FALSE.
    reflY(this%XVELOCITY) = .FALSE.
    reflY(this%YVELOCITY) = .TRUE.
    reflY(this%PRESSURE)  = .FALSE.
  END SUBROUTINE ReflectionMasks_euler2D


  PURE SUBROUTINE AxisMasks_euler2D(this,reflX,reflY)
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
    reflX(this%PRESSURE)  = .FALSE.
    ! southern / northern boundary
    reflY(this%DENSITY)   = .FALSE.
    reflY(this%XVELOCITY) = .TRUE.
    reflY(this%YVELOCITY) = .TRUE.
    reflY(this%PRESSURE)  = .FALSE.
  END SUBROUTINE AxisMasks_euler2D


  PURE SUBROUTINE CalcSoundSpeeds_center(this,Mesh,pvar)
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
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
      DO i=Mesh%IGMIN,Mesh%IGMAX
        this%bccsound(i,j) = GetSoundSpeed_euler2D(&
          this%gamma,&
          pvar(i,j,this%DENSITY),&
          pvar(i,j,this%PRESSURE))
      END DO
    END DO
  END SUBROUTINE CalcSoundSpeeds_center


  PURE SUBROUTINE CalcSoundSpeeds_faces(this,Mesh,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%VNUM) &
         :: prim
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,prim
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
!CDIR COLLAPSE
    DO k=1,4
      DO j=Mesh%JGMIN,Mesh%JGMAX
        DO i=Mesh%IGMIN,Mesh%IGMAX
          this%fcsound(i,j,k) = GetSoundSpeed_euler2D(&
            this%gamma,&
            prim(i,j,k,this%DENSITY),&
            prim(i,j,k,this%PRESSURE))
        END DO
      END DO
    END DO
  END SUBROUTINE CalcSoundSpeeds_faces


  ELEMENTAL FUNCTION GetSoundSpeed_euler2D(gamma,density,pressure) RESULT(cs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,density,pressure
    REAL :: cs
    !------------------------------------------------------------------------!
    cs = SQRT(MAX(2.0*TINY(cs),gamma*pressure/density))
  END FUNCTION GetSoundSpeed_euler2D

  ELEMENTAL SUBROUTINE SetWaveSpeeds_euler2D(cs,v,amin,amax)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,v
    REAL, INTENT(OUT) :: amin,amax
    !------------------------------------------------------------------------!
    ! minimal and maximal wave speeds
    amin = MIN(0.,v-cs)
    amax = MAX(0.,v+cs)
  END SUBROUTINE SetWaveSpeeds_euler2D

  ELEMENTAL SUBROUTINE SetEigenValues_euler2D(gamma,rho,v,P,l1,l2,l3,l4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho,v,P
    REAL, INTENT(OUT) :: l1,l2,l3,l4
    !------------------------------------------------------------------------!
    REAL :: cs
    !------------------------------------------------------------------------!
    ! adiabatic sound speed
!CDIR IEXPAND
    cs = GetSoundSpeed_euler2D(gamma,rho,P)
    ! call subroutine for isothermal case with the adiabatic sound speed
!CDIR IEXPAND
    CALL SetEigenValues_euler2Dit(cs,v,l1,l2,l4)
    ! set the missing eigenvalue l3
    l3 = v
  END SUBROUTINE SetEigenValues_euler2D


  ELEMENTAL SUBROUTINE SetCharVars_euler2D(gamma,rho1,rho2,u1,u2,v1,v2,P1,P2, &
       l1,l2,l3,l4,xvar1,xvar2,xvar3,xvar4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho1,rho2,u1,u2,v1,v2,P1,P2,l1,l2,l3,l4
    REAL, INTENT(OUT) :: xvar1,xvar2,xvar3,xvar4
    !------------------------------------------------------------------------!
    REAL :: gamcs,dlnP,du
    !------------------------------------------------------------------------!
    gamcs= gamma / (l4-l1) ! = 2*gamma/cs
    dlnP = LOG(P2/P1)         ! = LOG(P2)-LOG(P1)
    du   = u2-u1
    ! characteristic variables
    xvar1 = dlnP - gamcs * du
    xvar2 = -dlnP + gamma * LOG(rho2/rho1)
    xvar3 = (v2-v1)
    xvar4 = dlnP + gamcs * du 
  END SUBROUTINE SetCharVars_euler2D


  ELEMENTAL SUBROUTINE SetBoundaryData_euler2D(gamma,dir,rho1,u1,v1,P1,xvar1, &
       xvar2,xvar3,xvar4,rho2,u2,v2,P2)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,dir,rho1,u1,v1,P1,xvar1,xvar2,xvar3,xvar4
    REAL, INTENT(OUT) :: rho2,u2,v2,P2
    !------------------------------------------------------------------------!
    REAL :: dlnP,csgam
    !------------------------------------------------------------------------!
    dlnP = 0.5 * (xvar4+xvar1)
    ! extrapolate boundary values using characteristic variables
    rho2 = rho1 * EXP(dir*(dlnP-xvar2)/gamma)
    P2   = P1 * EXP(dir*dlnP)
!CDIR IEXPAND
    csgam= GetSoundSpeed_euler2D(gamma,rho1+rho2,P1+P2) / gamma
    u2   = u1 + dir*csgam * 0.5*(xvar4-xvar1)
    v2   = v1 + dir*xvar3
  END SUBROUTINE SetBoundaryData_euler2D

 ELEMENTAL SUBROUTINE Prim2Riemann_euler2D(gamma,rho,vx,vy,p,&
                                       l1,l2,l3,l4,Rminus,Rs,Rvt,Rplus)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho,vx,vy,p,l1,l2,l3,l4
    REAL, INTENT(OUT) :: Rminus,Rs,Rvt,Rplus
    !------------------------------------------------------------------------!
    REAL :: cs
    !------------------------------------------------------------------------!
    cs = l4-l2 ! l2 = v, l4 = v+cs
    ! compute 1st Riemann invariant (R+)
    Rplus = vx + 2./(gamma-1.0) * cs     
    ! compute 2st Riemann invariant (R-) 
    Rminus = vx - 2./(gamma-1.0) * cs
    ! compute entropy
    Rs = p/rho**gamma
    ! tangential velocities
    Rvt = vy   
  END SUBROUTINE Prim2Riemann_euler2D

  ELEMENTAL SUBROUTINE Riemann2Prim_euler2D(gamma,Rminus,Rs,Rvt,Rplus,&
       rho,vx,vy,p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,Rminus,Rs,Rvt,Rplus
    REAL, INTENT(OUT) :: rho,vx,vy,p
    !------------------------------------------------------------------------!
    REAL :: cs2gam
    !------------------------------------------------------------------------!
    ! tangential velocity    
    vy = Rvt  
    ! normal velocity
    vx = 0.5*(Rplus+Rminus)
    ! cs**2 / gamma
    cs2gam = (0.25*(gamma-1.0)*(Rplus-Rminus))**2 / gamma
    ! density
    rho = (cs2gam/Rs)**(1./(gamma-1.0))
    ! pressure
    p = cs2gam * rho
  END SUBROUTINE Riemann2Prim_euler2D

  ELEMENTAL SUBROUTINE CalcFlux_euler2D(rho,v,P,m1,m2,E,f1,f2,f3,f4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho,v,P,m1,m2,E
    REAL, INTENT(OUT) :: f1, f2, f3, f4
    !------------------------------------------------------------------------!
    f1 = rho*v
    f2 = m1*v + P
    f3 = m2*v
    f4 = (E+P)*v
  END SUBROUTINE CalcFlux_euler2D


  ELEMENTAL SUBROUTINE Cons2Prim_euler2D(gamma,rho_in,mu,mv,E,rho_out,u,v,P)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho_in,mu,mv,E
    REAL, INTENT(OUT) :: rho_out,u,v,P
    !------------------------------------------------------------------------!
    REAL :: inv_rho
    !------------------------------------------------------------------------!
    inv_rho = 1./rho_in
    rho_out = rho_in
    u = mu * inv_rho
    v = mv * inv_rho
    P = (gamma-1.)*(E - 0.5 * inv_rho * (mu*mu+mv*mv))
  END SUBROUTINE Cons2Prim_euler2D

  
  ELEMENTAL SUBROUTINE Prim2Cons_euler2D(gamma,rho_in,u,v,P,rho_out,mu,mv,E)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho_in,u,v,P
    REAL, INTENT(OUT) :: rho_out,mu,mv,E
    !------------------------------------------------------------------------!
    rho_out = rho_in
    mu = rho_in * u
    mv = rho_in * v
    E = P/(gamma-1.) + 0.5 * rho_in * (u*u+v*v)
  END SUBROUTINE Prim2Cons_euler2D
  

  SUBROUTINE ClosePhysics_euler2D(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL ClosePhysics_euler2Dit(this)
  END SUBROUTINE ClosePhysics_euler2D

END MODULE physics_euler2D
