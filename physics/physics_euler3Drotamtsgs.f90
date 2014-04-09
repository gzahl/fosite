!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler3Drotamtsgs.f90                                      #
!#                                                                           #
!# Copyright (C) 2012                                                        #
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
!> \author Björn Sperling
!! \author Tobias Illenseer
!!
!! \brief basic module for 3D Euler equations with rotational symmetry and
!! angular momentum transport and SGS support - test version
!!
!! \extends physics_common
!! \ingroup physics
!----------------------------------------------------------------------------!
MODULE physics_euler3Drotamtsgs
  USE physics_common
  USE mesh_common, ONLY : Mesh_TYP
  USE sources_common, ONLY : Sources_TYP
  USE physics_euler2D
  USE physics_euler3Drotamt, &
       CalcStresses_euler3Drasgs => CalcStresses_euler3Dra
  USE physics_euler3Drotsymsgs, &
       CalcWaveSpeeds_euler3Drasgs => CalcWaveSpeeds_euler3Drssgs, &
       CalcSoundSpeeds_euler3Drasgs => CalcSoundSpeeds_euler3Drssgs, &
       CalcFluxesX_euler3Drasgs => CalcFluxesX_euler3Drssgs, &
       CalcFluxesY_euler3Drasgs => CalcFluxesY_euler3Drssgs, &
       SetEigenValues_euler3Drasgs => SetEigenValues_euler3Drssgs, &
       ReflectionMasks_euler3Drasgs => ReflectionMasks_euler3Drssgs, &
       AxisMasks_euler3Drasgs => AxisMasks_euler3Drssgs, &
       ClosePhysics_euler3Drasgs => ClosePhysics_euler3Drssgs
  USE mesh_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE GeometricalSources_euler3Drasgs
     MODULE PROCEDURE GeometricalSources_center
     MODULE PROCEDURE GeometricalSources_faces
  END INTERFACE
  INTERFACE Convert2Primitive_euler3Drasgs
     MODULE PROCEDURE Convert2Primitive_center
     MODULE PROCEDURE Convert2Primitive_faces
  END INTERFACE
  INTERFACE Convert2Cons_euler3Drasgs
     MODULE PROCEDURE Convert2Conservative_center
     MODULE PROCEDURE Convert2Conservative_faces
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: num_var = 6
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 3D w/ ang. momentum and SGS model"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitPhysics_euler3Drasgs, &
       CalcWaveSpeeds_euler3Drasgs, &
       CalcSoundSpeeds_euler3Drasgs, &
       CalcFluxesX_euler3Drasgs, &
       CalcFluxesY_euler3Drasgs, &
       CalcCharSystemX_euler3Drasgs, &
       CalcCharSystemY_euler3Drasgs, &
       CalcBoundaryDataX_euler3Drasgs, &
       CalcBoundaryDataY_euler3Drasgs, &
       CalcRiemann2PrimX_euler3Drasgs, &
       CalcRiemann2PrimY_euler3Drasgs, &
       CalcPrim2RiemannX_euler3Drasgs, &
       CalcPrim2RiemannY_euler3Drasgs, &
       CalcStresses_euler3drasgs, &
       GeometricalSources_euler3Drasgs, &
       ExternalSources_euler3Drasgs, &
       ViscositySources_euler3Drasgs, &
       SGSSources_euler3Drasgs, &
       CalcSGSTensor_euler3Drasgs, &
       Convert2Primitive_euler3Drasgs, &
       Convert2Cons_euler3Drasgs, &
       ReflectionMasks_euler3Drasgs, &
       AxisMasks_euler3Drasgs, &
       ClosePhysics_euler3Drasgs
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitPhysics_euler3Drasgs(this,Mesh,problem)
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
    this%DENSITY     = 1                               ! mass density        !
    this%PRESSURE    = num_var                         ! pressure            !
    this%ENERGY      = num_var                         ! total energy        !
    this%SGSPRESSURE = 5                               ! turbulent pressure  !
    this%SGSENERGY   = 5                               ! turbulent energy    !
    this%XVELOCITY   = 2                               ! x-velocity          !
    this%XMOMENTUM   = 2                               ! x-momentum          !
    this%YVELOCITY   = 3                               ! y-velocity          !
    this%YMOMENTUM   = 3                               ! y-momentum          !
    this%ZVELOCITY   = 4                         ! specific angular momentum !
    this%ZMOMENTUM   = 4                         ! angular momentum          !
    ! set names for primitive and conservative variables
    this%pvarname(this%DENSITY)   = "density"
    this%pvarname(this%XVELOCITY) = "xvelocity"
    this%pvarname(this%YVELOCITY) = "yvelocity"
    this%pvarname(this%ZVELOCITY) = "specangmomentum"
    this%pvarname(this%SGSPRESSURE) = "sgspressure"
    this%pvarname(this%PRESSURE)  = "pressure"
    this%cvarname(this%DENSITY)   = "density"
    this%cvarname(this%XMOMENTUM) = "xmomentum"
    this%cvarname(this%YMOMENTUM) = "ymomentum"
    this%cvarname(this%ZMOMENTUM) = "angmomentum"
    this%cvarname(this%SGSENERGY)   = "sgsenergy"
    this%cvarname(this%ENERGY)    = "energy"
    this%DIM = 3

    ALLOCATE(this%fcent(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,2), &
         STAT = err)
    ! abort if allocation fails
    IF (err.NE.0) &
         CALL Error(this, "InitPhysics_euler3Drasgs", "Unable to allocate memory.")

  END SUBROUTINE InitPhysics_euler3Drasgs

   PURE SUBROUTINE CalcCharSystemX_euler3Drasgs(this,Mesh,i,dir,pvar,lambda,xvar)
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
     CALL SetEigenValues_euler3Drasgs(this%gamma,pvar(i,:,this%DENSITY), &
          pvar(i,:,this%XVELOCITY),pvar(i,:,this%SGSPRESSURE),pvar(i,:,this%PRESSURE), &
          lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4),lambda(:,5),lambda(:,6))
     ! compute characteristic variables
     i1 = i + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
     i2 = MAX(i,i1)
     i1 = MIN(i,i1)
 !CDIR IEXPAND
     CALL SetCharVars_euler3DraSGS(this%gamma,pvar(i1,:,this%DENSITY), &
          pvar(i2,:,this%DENSITY),pvar(i1,:,this%XVELOCITY), &
          pvar(i2,:,this%XVELOCITY),pvar(i1,:,this%YVELOCITY), &
          pvar(i2,:,this%YVELOCITY),pvar(i1,:,this%ZVELOCITY), &
          pvar(i2,:,this%ZVELOCITY),pvar(i1,:,this%SGSPRESSURE), &
          pvar(i2,:,this%SGSPRESSURE),pvar(i1,:,this%PRESSURE), &
          pvar(i2,:,this%PRESSURE),lambda(:,1),lambda(:,2), &
          lambda(:,3),lambda(:,4),lambda(:,5),lambda(:,6),xvar(:,1),xvar(:,2), &
          xvar(:,3),xvar(:,4),xvar(:,5),xvar(:,6))
   END SUBROUTINE CalcCharSystemX_euler3Drasgs
  

   PURE SUBROUTINE CalcCharSystemY_euler3Drasgs(this,Mesh,j,dir,pvar,lambda,xvar)
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
     CALL SetEigenValues_euler3Drasgs(this%gamma,pvar(:,j,this%DENSITY), &
          pvar(:,j,this%YVELOCITY),pvar(:,j,this%SGSPRESSURE),pvar(:,j,this%PRESSURE), &
          lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4),lambda(:,5),lambda(:,6))
     ! compute characteristic variables
     j1 = j + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
     j2 = MAX(j,j1)
     j1 = MIN(j,j1)
 !CDIR IEXPAND
     CALL SetCharVars_euler3Drasgs(this%gamma,pvar(:,j1,this%DENSITY), &
          pvar(:,j2,this%DENSITY),pvar(:,j1,this%YVELOCITY), &
          pvar(:,j2,this%YVELOCITY),pvar(:,j1,this%XVELOCITY), &
          pvar(:,j2,this%XVELOCITY),pvar(:,j1,this%ZVELOCITY), &
          pvar(:,j2,this%ZVELOCITY),pvar(:,j1,this%SGSPRESSURE), &
          pvar(:,j2,this%SGSPRESSURE),pvar(:,j1,this%PRESSURE), &
          pvar(:,j2,this%PRESSURE),lambda(:,1),lambda(:,2), &
          lambda(:,3),lambda(:,4),lambda(:,5),lambda(:,6),xvar(:,1),xvar(:,2), &
          xvar(:,3),xvar(:,4),xvar(:,5),xvar(:,6))
   END SUBROUTINE CalcCharSystemY_euler3Drasgs


  PURE SUBROUTINE CalcBoundaryDataX_euler3Drasgs(this,Mesh,i1,dir,xvar,pvar)
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
     CALL SetBoundaryData_euler3Drasgs(this%gamma,1.0*SIGN(1,dir),pvar(i1,:,this%DENSITY), &
          pvar(i1,:,this%XVELOCITY),pvar(i1,:,this%YVELOCITY),pvar(i1,:,this%ZVELOCITY), &
          pvar(i1,:,this%SGSPRESSURE),pvar(i1,:,this%PRESSURE),xvar(:,1),xvar(:,2),xvar(:,3), &
          xvar(:,4),xvar(:,5),xvar(:,6),pvar(i2,:,this%DENSITY),pvar(i2,:,this%XVELOCITY), &
          pvar(i2,:,this%YVELOCITY),pvar(i2,:,this%ZVELOCITY),pvar(i2,:,this%SGSPRESSURE), &
          pvar(i2,:,this%PRESSURE))
   END SUBROUTINE CalcBoundaryDataX_euler3Drasgs
 
   PURE SUBROUTINE CalcBoundaryDataY_euler3Drasgs(this,Mesh,j1,dir,xvar,pvar)
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
     CALL SetBoundaryData_euler3Drasgs(this%gamma,1.0*SIGN(1,dir),pvar(:,j1,this%DENSITY), &
          pvar(:,j1,this%YVELOCITY),pvar(:,j1,this%XVELOCITY),pvar(:,j1,this%ZVELOCITY), &
          pvar(:,j1,this%SGSPRESSURE),pvar(:,j1,this%PRESSURE),xvar(:,1),xvar(:,2),xvar(:,3), &
          xvar(:,4),xvar(:,5),xvar(:,6),pvar(:,j2,this%DENSITY),pvar(:,j2,this%YVELOCITY), &
          pvar(:,j2,this%XVELOCITY), pvar(:,j2,this%ZVELOCITY),pvar(:,j2,this%SGSPRESSURE), &
          pvar(:,j2,this%PRESSURE))
   END SUBROUTINE CalcBoundaryDataY_euler3Drasgs

PURE SUBROUTINE CalcPrim2RiemannX_euler3Drasgs(this,Mesh,i,pvar,lambda,Rinv)
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
    CALL SetEigenValues_euler3Drasgs(this%gamma,pvar(i,:,this%DENSITY), &
         pvar(i,:,this%XVELOCITY),pvar(i,:,this%SGSPRESSURE),pvar(i,:,this%PRESSURE),&
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4),lambda(:,5),lambda(:,6))
    ! compute Riemann invariants
!CDIR IEXPAND
    CALL Prim2Riemann_euler3Drasgs(this%gamma,pvar(i,:,this%DENSITY), &
         pvar(i,:,this%XVELOCITY),pvar(i,:,this%YVELOCITY), &
         pvar(i,:,this%ZVELOCITY),pvar(i,:,this%SGSPRESSURE),pvar(i,:,this%PRESSURE), &
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4),lambda(:,5),lambda(:,6), &
         Rinv(:,1),Rinv(:,2),Rinv(:,3),Rinv(:,4),Rinv(:,5),Rinv(:,6))
  END SUBROUTINE CalcPrim2RiemannX_euler3Drasgs


  PURE SUBROUTINE CalcPrim2RiemannY_euler3Drasgs(this,Mesh,j,pvar,lambda,Rinv)
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
    CALL SetEigenValues_euler3Drasgs(this%gamma,pvar(:,j,this%DENSITY), &
         pvar(:,j,this%YVELOCITY),pvar(:,j,this%SGSPRESSURE),pvar(:,j,this%PRESSURE),&
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4),lambda(:,5),lambda(:,6))
    ! compute Riemann invariants
!CDIR IEXPAND
    CALL Prim2Riemann_euler3Drasgs(this%gamma,pvar(:,j,this%DENSITY), &
         pvar(:,j,this%YVELOCITY),pvar(:,j,this%XVELOCITY), &
         pvar(:,j,this%ZVELOCITY),pvar(:,j,this%SGSPRESSURE),pvar(:,j,this%PRESSURE), &
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4),lambda(:,5),lambda(:,6), &
         Rinv(:,1),Rinv(:,2),Rinv(:,3),Rinv(:,4),Rinv(:,5),Rinv(:,6))
  END SUBROUTINE CalcPrim2RiemannY_euler3Drasgs

  PURE SUBROUTINE CalcRiemann2PrimX_euler3Drasgs(this,Mesh,i,Rinv,pvar)
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
    CALL Riemann2Prim_euler3Drasgs(this%gamma,Rinv(:,1),Rinv(:,2),&
       Rinv(:,3),Rinv(:,4),Rinv(:,5),Rinv(:,6),pvar(i,:,this%Density), &
       pvar(i,:,this%XVELOCITY),pvar(i,:,this%YVELOCITY),pvar(i,:,this%ZVELOCITY),&
       pvar(i,:,this%SGSPRESSURE),pvar(i,:,this%PRESSURE))
  END SUBROUTINE CalcRiemann2PrimX_euler3Drasgs

  PURE SUBROUTINE CalcRiemann2PrimY_euler3Drasgs(this,Mesh,j,Rinv,pvar)
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
    CALL Riemann2Prim_euler3Drasgs(this%gamma,Rinv(:,1),Rinv(:,2),&
       Rinv(:,3),Rinv(:,4),Rinv(:,5),Rinv(:,6),pvar(:,j,this%Density),&
       pvar(:,j,this%YVELOCITY),pvar(:,j,this%XVELOCITY),pvar(:,j,this%ZVELOCITY),&
       pvar(:,j,this%SGSPRESSURE),pvar(:,j,this%PRESSURE))
  END SUBROUTINE CalcRiemann2PrimY_euler3Drasgs

  PURE SUBROUTINE GeometricalSources_center(this,Mesh,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
         :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,cvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! calculate centrifugal forces
!CDIR IEXPAND
    CALL CentrifugalForces(pvar(:,:,this%DENSITY),pvar(:,:,this%ZVELOCITY), &
         Mesh%bhz(:,:),Mesh%czxz(:,:,1),Mesh%czyz(:,:,1),this%fcent(:,:,1,1), &
         this%fcent(:,:,2,1))

    ! no geometrical sources in continuity and angular momentum equation
    sterm(:,:,this%DENSITY)   = 0.0
    sterm(:,:,this%SGSENERGY) = 0.0
    sterm(:,:,this%ZMOMENTUM) = 0.0

    ! geometrical source terms in momentum equationes
    ! with centrifugal forces
!CDIR IEXPAND
    sterm(:,:,this%XMOMENTUM) = MomentumSourcesX_euler2D(cvar(:,:,this%YMOMENTUM), &
         pvar(:,:,this%XVELOCITY),pvar(:,:,this%YVELOCITY),pvar(:,:,this%PRESSURE), &
         Mesh%cxyx(:,:,1),Mesh%cyxy(:,:,1),Mesh%czxz(:,:,1)) + this%fcent(:,:,1,1)
!CDIR IEXPAND
    sterm(:,:,this%YMOMENTUM) = MomentumSourcesY_euler2D(cvar(:,:,this%XMOMENTUM), &
         pvar(:,:,this%XVELOCITY),pvar(:,:,this%YVELOCITY),pvar(:,:,this%PRESSURE), &
         Mesh%cxyx(:,:,1),Mesh%cyxy(:,:,1),Mesh%czyz(:,:,1)) + this%fcent(:,:,2,1)

    ! centrifugal force source terms in energy equation
    sterm(:,:,this%ENERGY) = this%fcent(:,:,1,1) * pvar(:,:,this%XVELOCITY) &
         + this%fcent(:,:,2,1) * pvar(:,:,this%YVELOCITY)
  END SUBROUTINE GeometricalSources_center


  PURE SUBROUTINE GeometricalSources_faces(this,Mesh,prim,cons,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%vnum) &
         :: prim,cons
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
         :: sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,prim,cons
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! calculate centrifugal forces
!CDIR IEXPAND
    CALL CentrifugalForces(prim(:,:,:,this%DENSITY),prim(:,:,:,this%ZVELOCITY), &
         Mesh%chz(:,:,:),Mesh%czxz(:,:,:),Mesh%czyz(:,:,:),this%fcent(:,:,:,1), &
         this%fcent(:,:,:,2))

    ! no geometrical sources in continuity and angular momentum equation
    sterm(:,:,this%DENSITY)   = 0.0
    sterm(:,:,this%SGSENERGY) = 0.0
    sterm(:,:,this%ZMOMENTUM) = 0.0

    ! geometrical source terms in momentum equationes
    ! sum up all four corner values
!CDIR IEXPAND
    sterm(:,:,this%XMOMENTUM) = SUM(MomentumSourcesX_euler2D(cons(:,:,:,this%YMOMENTUM), &
         prim(:,:,:,this%XVELOCITY),prim(:,:,:,this%YVELOCITY),prim(:,:,:,this%PRESSURE), &
         Mesh%cxyx(:,:,:),Mesh%cyxy(:,:,:),Mesh%czxz(:,:,:)) + this%fcent(:,:,:,1),DIM=3)
!CDIR IEXPAND
    sterm(:,:,this%YMOMENTUM) = SUM(MomentumSourcesY_euler2D(cons(:,:,:,this%XMOMENTUM), &
         prim(:,:,:,this%XVELOCITY),prim(:,:,:,this%YVELOCITY),prim(:,:,:,this%PRESSURE), &
         Mesh%cxyx(:,:,:),Mesh%cyxy(:,:,:),Mesh%czyz(:,:,:)) + this%fcent(:,:,:,2),DIM=3)

    ! centrifugal force source terms in energy equation
    sterm(:,:,this%ENERGY) = SUM(this%fcent(:,:,:,1) * prim(:,:,:,this%XVELOCITY) &
         + this%fcent(:,:,:,2) * prim(:,:,:,this%YVELOCITY), DIM=3)
  END SUBROUTINE GeometricalSources_faces

  
 PURE SUBROUTINE ExternalSources_euler3Drasgs(this,Mesh,accel,pvar,cvar,sterm)
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
         :: accel
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
         :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,accel,pvar,cvar
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    CALL ExternalSources_euler3Dra(this,Mesh,accel,pvar,cvar,sterm) 
    sterm(:,:,this%SGSENERGY) = 0.
  END SUBROUTINE ExternalSources_euler3Drasgs

  PURE SUBROUTINE ViscositySources_euler3Drasgs(this,Mesh,pvar,btxx,btxy,btxz, &
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
    CALL ViscositySources_euler3Dra(this,Mesh,pvar,btxx,btxy,btxz, &
       btyy,btyz,btzz,sterm)
    sterm(:,:,this%SGSENERGY) = 0.0
    
  END SUBROUTINE ViscositySources_euler3Drasgs


 PURE SUBROUTINE Convert2Primitive_center(this,Mesh,i1,i2,j1,j2,cvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    INTEGER              :: i1,i2,j1,j2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
                         :: cvar,pvar
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,i1,i2,j1,j2,cvar
    INTENT(OUT) :: pvar
    !------------------------------------------------------------------------!
    CALL Cons2Prim_euler3Drasgs(this%gamma,cvar(i1:i2,j1:j2,this%DENSITY), &
         cvar(i1:i2,j1:j2,this%XMOMENTUM),cvar(i1:i2,j1:j2,this%YMOMENTUM), &
         cvar(i1:i2,j1:j2,this%ZMOMENTUM),cvar(i1:i2,j1:j2,this%SGSENERGY), &
         cvar(i1:i2,j1:j2,this%ENERGY), &
         pvar(i1:i2,j1:j2,this%DENSITY),pvar(i1:i2,j1:j2,this%XVELOCITY), &
         pvar(i1:i2,j1:j2,this%YVELOCITY),pvar(i1:i2,j1:j2,this%ZVELOCITY), &
         pvar(i1:i2,j1:j2,this%SGSPRESSURE),pvar(i1:i2,j1:j2,this%PRESSURE))
  END SUBROUTINE Convert2Primitive_center


  PURE SUBROUTINE Convert2Primitive_faces(this,Mesh,i1,i2,j1,j2,cons,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    INTEGER              :: i1,i2,j1,j2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%vnum) &
                         :: cons,prim
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,i1,i2,j1,j2,cons
    INTENT(OUT) :: prim
    !------------------------------------------------------------------------!
    CALL Cons2Prim_euler3Drasgs(this%gamma,cons(i1:i2,j1:j2,:,this%DENSITY), &
         cons(i1:i2,j1:j2,:,this%XMOMENTUM),cons(i1:i2,j1:j2,:,this%YMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%ZMOMENTUM),cons(i1:i2,j1:j2,:,this%SGSENERGY), &
         cons(i1:i2,j1:j2,:,this%ENERGY), &
         prim(i1:i2,j1:j2,:,this%DENSITY),prim(i1:i2,j1:j2,:,this%XVELOCITY), &
         prim(i1:i2,j1:j2,:,this%YVELOCITY),prim(i1:i2,j1:j2,:,this%ZVELOCITY), &
         prim(i1:i2,j1:j2,:,this%SGSPRESSURE),prim(i1:i2,j1:j2,:,this%PRESSURE))
  END SUBROUTINE Convert2Primitive_faces


  PURE SUBROUTINE Convert2Conservative_center(this,Mesh,i1,i2,j1,j2,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    INTEGER              :: i1,i2,j1,j2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
                         :: cvar,pvar
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,i1,i2,j1,j2,pvar
    INTENT(OUT) :: cvar
    !------------------------------------------------------------------------!
    CALL Prim2Cons_euler3Drasgs(this%gamma,pvar(i1:i2,j1:j2,this%DENSITY), &
         pvar(i1:i2,j1:j2,this%XVELOCITY),pvar(i1:i2,j1:j2,this%YVELOCITY), &
         pvar(i1:i2,j1:j2,this%ZVELOCITY),pvar(i1:i2,j1:j2,this%SGSPRESSURE), &
         pvar(i1:i2,j1:j2,this%PRESSURE), &
         cvar(i1:i2,j1:j2,this%DENSITY),cvar(i1:i2,j1:j2,this%XMOMENTUM), &
         cvar(i1:i2,j1:j2,this%YMOMENTUM),cvar(i1:i2,j1:j2,this%ZMOMENTUM), &
         cvar(i1:i2,j1:j2,this%SGSENERGY),cvar(i1:i2,j1:j2,this%ENERGY))
  END SUBROUTINE Convert2Conservative_center


  PURE SUBROUTINE Convert2Conservative_faces(this,Mesh,i1,i2,j1,j2,prim,cons)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    INTEGER              :: i1,i2,j1,j2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%vnum) &
                         :: cons,prim
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,i1,i2,j1,j2,prim
    INTENT(OUT) :: cons
    !------------------------------------------------------------------------!
    CALL Prim2Cons_euler3Drasgs(this%gamma,prim(i1:i2,j1:j2,:,this%DENSITY), &
         prim(i1:i2,j1:j2,:,this%XVELOCITY),prim(i1:i2,j1:j2,:,this%YVELOCITY), &
         prim(i1:i2,j1:j2,:,this%ZVELOCITY),prim(i1:i2,j1:j2,:,this%SGSPRESSURE), &
         prim(i1:i2,j1:j2,:,this%PRESSURE), &
         cons(i1:i2,j1:j2,:,this%DENSITY),cons(i1:i2,j1:j2,:,this%XMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%YMOMENTUM),cons(i1:i2,j1:j2,:,this%ZMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%SGSENERGY),cons(i1:i2,j1:j2,:,this%ENERGY))
  END SUBROUTINE Convert2Conservative_faces

  PURE SUBROUTINE SGSSources_euler3Drasgs(this,Mesh,Sources,pvar,cvar,sterm)         
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Sources_TYP) :: Sources
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: &
         pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,cvar
    INTENT(INOUT)     :: this,Sources
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    sterm(:,:,this%DENSITY)   = 0.0

    CALL Divergence(Mesh,Sources%btxx,Sources%btxy,Sources%btxz, &
         Sources%btxy,Sources%btyy,Sources%btyz, &
         Sources%btxz,Sources%btyz,Sources%btzz, &
         sterm(:,:,this%XMOMENTUM), &
         sterm(:,:,this%YMOMENTUM), &
         sterm(:,:,this%ZMOMENTUM)) 

    sterm(:,:,this%ZMOMENTUM) = sterm(:,:,this%ZMOMENTUM) * Mesh%bhz(:,:)

    !(v*tau)_x; amin as temporal storage
    this%amin(:,:) = pvar(:,:,this%XVELOCITY)*Sources%btxx(:,:) &
                    +pvar(:,:,this%YVELOCITY)*Sources%btxy(:,:) &
                    +pvar(:,:,this%ZVELOCITY)*Sources%btxz(:,:) / Mesh%bhz(:,:)

    !(v*tau)_y; amax as temporal storage
    this%amax(:,:) = pvar(:,:,this%XVELOCITY)*Sources%btxy(:,:) &
                    +pvar(:,:,this%YVELOCITY)*Sources%btyy(:,:) &
                    +pvar(:,:,this%ZVELOCITY)*Sources%btyz(:,:)  / Mesh%bhz(:,:)
    ! div(v*t)
    CALL Divergence(Mesh,this%amin,this%amax,sterm(:,:,this%ENERGY))

    !add div(v*t) + rhoepsilon - Sigma 
    sterm(:,:,this%ENERGY) = sterm(:,:,this%ENERGY)+Sources%rhoeps(:,:)-Sources%sigma(:,:)
   
    !add D - rhoepsilon + Sigma 
    sterm(:,:,this%SGSENERGY) = Sources%diff(:,:)-Sources%rhoeps(:,:)+Sources%sigma(:,:) 
   
  END SUBROUTINE SGSSources_euler3Drasgs

  PURE SUBROUTINE CalcSGSTensor_euler3Drasgs(this,Mesh,Sources,C,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Sources_TYP) :: Sources
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
                      :: pvar 
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: C 
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,C,pvar
    INTENT(INOUT)     :: this,Sources
    !------------------------------------------------------------------------!
!CDIR OUTERUNROLL=8
    DO j=Mesh%JMIN-1,Mesh%JMAX+1
!CDIR NODEP
       DO i=Mesh%IMIN-1,Mesh%IMAX+1
          ! compute the diagonal elements of the Rate of Strain tensor
          Sources%Sxx(i,j) = &
             (pvar(i+1,j,this%XVELOCITY)-pvar(i-1,j,this%XVELOCITY))*0.5/Mesh%dlx(i,j)&
             +Mesh%cxyx(i,j,1)*pvar(i,j,this%YVELOCITY)

          Sources%Syy(i,j) = &
             (pvar(i,j+1,this%YVELOCITY)-pvar(i,j-1,this%YVELOCITY))*0.5/Mesh%dly(i,j)&
              +Mesh%cyxy(i,j,1)*pvar(i,j,this%XVELOCITY)

          Sources%Szz(i,j) = Mesh%czxz(i,j,1)*pvar(i,j,this%XVELOCITY) &
                            +Mesh%czyz(i,j,1)*pvar(i,j,this%YVELOCITY)

          ! compute the off-diagonal elements the Rate of Strain tensor
          Sources%Sxy(i,j) = &
              (pvar(i+1,j,this%YVELOCITY)-pvar(i-1,j,this%YVELOCITY))*0.25/Mesh%dlx(i,j)&
             +(pvar(i,j+1,this%XVELOCITY)-pvar(i,j-1,this%XVELOCITY))*0.25/Mesh%dly(i,j)&
             - 0.5*Mesh%cxyx(i,j,1)*pvar(i,j,this%XVELOCITY)&
             - 0.5*Mesh%cyxy(i,j,1)*pvar(i,j,this%YVELOCITY)

          Sources%Sxz(i,j) = &
              (pvar(i+1,j,this%ZVELOCITY)/Mesh%bhz(i+1,j) &
              -pvar(i-1,j,this%ZVELOCITY)/Mesh%bhz(i-1,j))*0.25/Mesh%dlx(i,j)&
              - 0.5*Mesh%czxz(i,j,1)*pvar(i,j,this%ZVELOCITY)/Mesh%bhz(i,j)

          Sources%Syz(i,j) = &
              (pvar(i,j+1,this%ZVELOCITY)/Mesh%bhz(i,j+1) &
              -pvar(i,j-1,this%ZVELOCITY)/Mesh%bhz(i,j-1))*0.25/Mesh%dly(i,j) &
              - 0.5*Mesh%czyz(i,j,1)*pvar(i,j,this%ZVELOCITY)/Mesh%bhz(i,j)
       END DO
    END DO

    ! compute divergence first and store the result in this%amin
!CDIR IEXPAND
    CALL Divergence(Mesh,pvar(:,:,this%XVELOCITY),pvar(:,:,this%YVELOCITY),this%amin(:,:))

    Sources%btxx(:,:) = C(:,:)*2.0*(Sources%Sxx(:,:) - 1.0/3.0*this%amin(:,:))
    Sources%btyy(:,:) = C(:,:)*2.0*(Sources%Syy(:,:) - 1.0/3.0*this%amin(:,:))   
    Sources%btzz(:,:) = C(:,:)*2.0*(Sources%Szz(:,:) - 1.0/3.0*this%amin(:,:))   
    Sources%btxy(:,:) = C(:,:)*2.0*Sources%Sxy(:,:)
    Sources%btxz(:,:) = C(:,:)*2.0*Sources%Sxz(:,:)
    Sources%btyz(:,:) = C(:,:)*2.0*Sources%Syz(:,:)

  END SUBROUTINE CalcSGSTensor_euler3Drasgs

  ELEMENTAL SUBROUTINE SetCharVars_euler3Drasgs(gamma,rho1,rho2,u1,u2,v1,v2,w1,w2, &
        sgsP1,sgsP2,P1,P2,l1,l2,l3,l4,l5,l6,xvar1,xvar2,xvar3,xvar4,xvar5,xvar6)
     IMPLICIT NONE
     !------------------------------------------------------------------------!
     REAL, INTENT(IN)  :: gamma,rho1,rho2,u1,u2,v1,v2,w1,w2,sgsP1,sgsP2,P1,P2, &
                          l1,l2,l3,l4,l5,l6
     REAL, INTENT(OUT) :: xvar1,xvar2,xvar3,xvar4,xvar5,xvar6
     !------------------------------------------------------------------------!
     REAL :: gamcs,dlnP,dlnRho,du
     !------------------------------------------------------------------------!
     gamcs= 2.*gamma / (l6-l1) ! = 2*gamma/cs
     dlnP = LOG((P2+sgsP2)/(P1+sgsP1))    ! = LOG(Peff2)-LOG(Peff1)
     dlnRho = LOG(rho2/rho1)
     du   = u2-u1
     ! characteristic variables
     xvar1 = dlnP - gamcs * du
     xvar2 = -dlnP + gamma * dlnRho
     xvar3 = gamcs * (v2-v1)
     xvar4 = -dlnP + gamma * ( LOG(w2/w1) + dlnRho )
     xvar5 = -dlnP + gamma * LOG(sgsP2/sgsP1)
     xvar6 = dlnP + gamcs * du 
   END SUBROUTINE SetCharVars_euler3Drasgs

  ELEMENTAL SUBROUTINE SetBoundaryData_euler3Drasgs(gamma,dir,rho1,u1,v1,w1,sgsP1,P1, &
        xvar1,xvar2,xvar3,xvar4,xvar5,xvar6,rho2,u2,v2,w2,sgsP2,P2)
     IMPLICIT NONE
     !------------------------------------------------------------------------!
     REAL, INTENT(IN)  :: gamma,dir,rho1,u1,v1,w1,sgsP1,P1,xvar1,xvar2,xvar3,xvar4,xvar5,xvar6
     REAL, INTENT(OUT) :: rho2,u2,v2,w2,sgsP2,P2
     !------------------------------------------------------------------------!
     REAL :: dlnP,gamdlnRho,csgam
     !------------------------------------------------------------------------!
     dlnP = 0.5 * (xvar6+xvar1)
     ! extrapolate boundary values using characteristic variables
     gamdlnRho = (xvar2+dlnP)
     rho2 = rho1 * EXP(dir*gamdlnRho/gamma)
     sgsP2 = sgsP1 * EXP(dir*(xvar5+dlnP)/gamma)
     P2   = (sgsP1+P1) * EXP(dir*dlnP) - sgsP2
!CDIR IEXPAND
     csgam= GetSoundSpeed_euler2D(gamma,rho1+rho2,sgsP1+P1+sgsP2+P2) / gamma
     u2   = u1 + dir*csgam * 0.5*(xvar6-xvar1)
     v2   = v1 + dir*csgam * xvar3
     w2   = w1 * EXP(dir*(xvar4+dlnP-gamdlnRho)/gamma)
   END SUBROUTINE SetBoundaryData_euler3Drasgs

 ELEMENTAL SUBROUTINE Prim2Riemann_euler3Drasgs(gamma,rho,vx,vy,lz,q,p,&
                                       l1,l2,l3,l4,l5,l6,Rminus,Rs,Rq,Rvt,Rl,Rplus)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho,vx,vy,lz,p,q,l1,l2,l3,l4,l5,l6
    REAL, INTENT(OUT) :: Rminus,Rs,Rq,Rvt,Rl,Rplus
    !------------------------------------------------------------------------!
    REAL :: cs
    !------------------------------------------------------------------------!
    cs = l6-l2 ! l2 = v, l6 = v+cs
    ! compute 1st Riemann invariant (R+)
    Rplus = vx + 2./(gamma-1.0) * cs     
    ! compute 2st Riemann invariant (R-) 
    Rminus = vx - 2./(gamma-1.0) * cs
    ! compute entropy
    Rs = (p+q)/rho**gamma
    ! compute invariant corresponding to psgs
    Rq = (p+q)/q**gamma
    ! compute invariant corresponding to l 
    Rl = (p+q)/(rho*lz)**gamma
    ! tangential velocities
    Rvt = vy   
  END SUBROUTINE Prim2Riemann_euler3Drasgs

  ELEMENTAL SUBROUTINE Riemann2Prim_euler3Drasgs(gamma,Rminus,Rs,Rq,Rvt,Rl,Rplus,&
       rho,vx,vy,lz,q,p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,Rminus,Rs,Rq,Rvt,Rl,Rplus
    REAL, INTENT(OUT) :: rho,vx,vy,lz,q,p
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
    ! turb. pressure
    q = (Rs/Rq)**(1.0/gamma)*rho
    ! pressure
    p = cs2gam * rho - q
    ! spec. ang. momen.
    lz = ((p+q)/Rl)**(1.0/gamma)/rho
  END SUBROUTINE Riemann2Prim_euler3Drasgs

  ELEMENTAL SUBROUTINE CentrifugalForces(rho,l,r,czxz,czyz,fcx,fcy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho,l,r,czxz,czyz
    REAL, INTENT(OUT) :: fcx,fcy
    !------------------------------------------------------------------------!
    REAL :: tmp2
    !------------------------------------------------------------------------!
    tmp2 = rho*(l/(r + TINY(1.0)))**2
    fcx  = tmp2 * czxz
    fcy  = tmp2 * czyz
  END SUBROUTINE CentrifugalForces


  ELEMENTAL SUBROUTINE Cons2Prim_euler3Drasgs(gamma,rho_in,mu,mv,L,sgsE,E,rho_out,u,v,K,sgsP,P)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho_in,mu,mv,L,sgsE,E
    REAL, INTENT(OUT) :: rho_out,u,v,K,sgsP,P
    !------------------------------------------------------------------------!
    REAL :: inv_rho
    !------------------------------------------------------------------------!

    inv_rho = 1./rho_in
    rho_out = rho_in
    u = mu * inv_rho
    v = mv * inv_rho
    K = L * inv_rho
    sgsP = sgsE/1.5
    P = (gamma-1.)*(E - 0.5 * inv_rho * (mu*mu+mv*mv))
  END SUBROUTINE Cons2Prim_euler3Drasgs

  
  ELEMENTAL SUBROUTINE Prim2Cons_euler3Drasgs(gamma,rho_in,u,v,K,sgsP,P,rho_out,mu,mv,L,sgsE,E)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho_in,u,v,K,sgsP,P
    REAL, INTENT(OUT) :: rho_out,mu,mv,L,sgsE,E
    !------------------------------------------------------------------------!

    rho_out = rho_in
    mu = rho_in * u
    mv = rho_in * v
    L  = rho_in * K
    sgsE = 1.5*sgsP
    E  = P/(gamma-1.) + 0.5 * rho_in * (u*u+v*v)
  END SUBROUTINE Prim2Cons_euler3Drasgs

END MODULE physics_euler3Drotamtsgs
