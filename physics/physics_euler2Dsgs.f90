!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler2Dsgs.f90                                            #
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
!! \brief basic module for 2D Euler equations and SGS model
!!
!! \extends physics_common
!! \ingroup physics
!----------------------------------------------------------------------------!
MODULE physics_euler2Dsgs
  USE physics_common
  USE mesh_common, ONLY : Mesh_TYP
  USE sources_common, ONLY : Sources_TYP
  USE physics_euler2Disothm, &
       CalcWaveSpeeds_euler2Dsgs => CalcWaveSpeeds_euler2Dit, &
       CalcStresses_euler2Dsgs => CalcStresses_euler2Dit, &
       SetWaveSpeeds_euler2Dsgs => SetWaveSpeeds_euler2Dit, &
       MomentumSourcesX_euler2Dsgs => MomentumSourcesX_euler2Dit, &
       MomentumSourcesY_euler2Dsgs => MomentumSourcesY_euler2Dit
  USE physics_euler2D
  USE mesh_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE CalcSoundSpeeds_euler2Dsgs
     MODULE PROCEDURE CalcSoundSpeeds_center
     MODULE PROCEDURE CalcSoundSpeeds_faces
  END INTERFACE
     INTERFACE GeometricalSources_euler2Dsgs
     MODULE PROCEDURE GeometricalSources_center
     MODULE PROCEDURE GeometricalSources_faces
  END INTERFACE
  INTERFACE Convert2Primitive_euler2Dsgs
     MODULE PROCEDURE Convert2Primitive_center
     MODULE PROCEDURE Convert2Primitive_faces
  END INTERFACE
  INTERFACE Convert2Conservative_euler2Dsgs
     MODULE PROCEDURE Convert2Conservative_center
     MODULE PROCEDURE Convert2Conservative_faces
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: num_var = 5              ! number of variables       !
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 2D and SGS model"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Physics_TYP, &
       ! methods
       InitPhysics_euler2Dsgs, &
       ClosePhysics_euler2Dsgs, &
       CalcWaveSpeeds_euler2Dsgs, &
       CalcFluxesX_euler2Dsgs, &
       CalcFluxesY_euler2Dsgs, &
       CalcCharSystemX_euler2Dsgs, &
       CalcCharSystemY_euler2Dsgs, &
       CalcBoundaryDataX_euler2Dsgs, &
       CalcBoundaryDataY_euler2Dsgs, &
       CalcRiemann2PrimX_euler2Dsgs, &
       CalcRiemann2PrimY_euler2Dsgs, &
       CalcPrim2RiemannX_euler2Dsgs, &
       CalcPrim2RiemannY_euler2Dsgs, &
       CalcStresses_euler2Dsgs, &
       CalcSoundSpeeds_euler2Dsgs, &
       GeometricalSources_euler2Dsgs, &
       ViscositySources_euler2Dsgs, &
       ExternalSources_euler2Dsgs, &
       SGSSources_euler2Dsgs, &
       CalcSGSTensor_euler2Dsgs, &
       Convert2Primitive_euler2Dsgs, &
       Convert2Conservative_euler2Dsgs, &
       ReflectionMasks_euler2Dsgs, &
       AxisMasks_euler2Dsgs, &
       SetEigenValues_euler2Dsgs, &
       SetWaveSpeeds_euler2Dsgs, &
       MomentumSourcesX_euler2Dsgs, &
       MomentumSourcesY_euler2Dsgs
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitPhysics_euler2Dsgs(this,Mesh,problem)
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
    this%SGSPRESSURE = 4                               ! turbulent pressure  !
    this%SGSENERGY   = 4                               ! turbulent energy    !
    this%XVELOCITY = 2                                 ! x-velocity          !
    this%XMOMENTUM = 2                                 ! x-momentum          !
    this%YVELOCITY = 3                                 ! y-velocity          !
    this%YMOMENTUM = 3                                 ! y-momentum          !
    this%ZVELOCITY = 0                                 ! no z-velocity       !
    this%ZMOMENTUM = 0                                 ! no z-momentum       !
    ! set names for primitive and conservative variables
    this%pvarname(this%DENSITY)     = "density"
    this%pvarname(this%XVELOCITY)   = "xvelocity"
    this%pvarname(this%YVELOCITY)   = "yvelocity"
    this%pvarname(this%SGSPRESSURE) = "sgspressure"
    this%pvarname(this%PRESSURE)    = "pressure"
    this%cvarname(this%DENSITY)     = "density"
    this%cvarname(this%XMOMENTUM)   = "xmomentum"
    this%cvarname(this%YMOMENTUM)   = "ymomentum"
    this%cvarname(this%SGSENERGY)   = "sgsenergy"
    this%cvarname(this%ENERGY)      = "energy"
    this%DIM = 2

CALL Warning(this, "InitPhysics_euler2Dsgs", "is a experimental module!")

  END SUBROUTINE InitPhysics_euler2Dsgs


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
          pvar(i,j,this%PRESSURE)+pvar(i,j,this%SGSPRESSURE))
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
            prim(i,j,k,this%PRESSURE)+prim(i,j,k,this%SGSPRESSURE))
        END DO
      END DO
    END DO
  END SUBROUTINE CalcSoundSpeeds_faces


  PURE SUBROUTINE CalcFluxesX_euler2Dsgs(this,Mesh,nmin,nmax,prim,cons,xfluxes)
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
    CALL CalcFlux_euler2Dsgs(prim(:,:,nmin:nmax,this%DENSITY), &
         prim(:,:,nmin:nmax,this%XVELOCITY),prim(:,:,nmin:nmax,this%SGSPRESSURE), &
         prim(:,:,nmin:nmax,this%PRESSURE), &
         cons(:,:,nmin:nmax,this%XMOMENTUM),cons(:,:,nmin:nmax,this%YMOMENTUM), &
         cons(:,:,nmin:nmax,this%SGSENERGY),cons(:,:,nmin:nmax,this%ENERGY), &
         xfluxes(:,:,nmin:nmax,this%DENSITY),xfluxes(:,:,nmin:nmax,this%XMOMENTUM), &
         xfluxes(:,:,nmin:nmax,this%YMOMENTUM),xfluxes(:,:,nmin:nmax,this%SGSENERGY), &
         xfluxes(:,:,nmin:nmax,this%ENERGY))
  END SUBROUTINE CalcFluxesX_euler2Dsgs


  PURE SUBROUTINE CalcFluxesY_euler2Dsgs(this,Mesh,nmin,nmax,prim,cons,yfluxes)
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
    CALL CalcFlux_euler2Dsgs(prim(:,:,nmin:nmax,this%DENSITY), &
         prim(:,:,nmin:nmax,this%YVELOCITY),prim(:,:,nmin:nmax,this%SGSPRESSURE), &
         prim(:,:,nmin:nmax,this%PRESSURE),cons(:,:,nmin:nmax,this%YMOMENTUM), &
         cons(:,:,nmin:nmax,this%XMOMENTUM),cons(:,:,nmin:nmax,this%SGSENERGY), &
         cons(:,:,nmin:nmax,this%ENERGY),yfluxes(:,:,nmin:nmax,this%DENSITY), &
         yfluxes(:,:,nmin:nmax,this%YMOMENTUM),yfluxes(:,:,nmin:nmax,this%XMOMENTUM), &
         yfluxes(:,:,nmin:nmax,this%SGSENERGY),yfluxes(:,:,nmin:nmax,this%ENERGY))
  END SUBROUTINE CalcFluxesY_euler2Dsgs


  PURE SUBROUTINE CalcCharSystemX_euler2Dsgs(this,Mesh,i,dir,pvar,lambda,xvar)
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
    CALL SetEigenValues_euler2Dsgs(this%gamma,pvar(i,:,this%DENSITY), &
         pvar(i,:,this%XVELOCITY), pvar(i,:,this%SGSPRESSURE), pvar(i,:,this%PRESSURE),&
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4),lambda(:,5))
    ! compute characteristic variables
    i1 = i + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
    i2 = MAX(i,i1)
    i1 = MIN(i,i1)
!CDIR IEXPAND
    CALL SetCharVars_euler2Dsgs(this%gamma,pvar(i1,:,this%DENSITY), &
         pvar(i2,:,this%DENSITY),pvar(i1,:,this%XVELOCITY), &
         pvar(i2,:,this%XVELOCITY),pvar(i1,:,this%YVELOCITY), &
         pvar(i2,:,this%YVELOCITY),pvar(i1,:,this%SGSPRESSURE), &
         pvar(i2,:,this%SGSPRESSURE),pvar(i1,:,this%PRESSURE), &
         pvar(i2,:,this%PRESSURE),lambda(:,1),lambda(:,2),lambda(:,3), &
         lambda(:,4),lambda(:,5),xvar(:,1),xvar(:,2),xvar(:,3),xvar(:,4),xvar(:,5))
  END SUBROUTINE CalcCharSystemX_euler2Dsgs


  PURE SUBROUTINE CalcCharSystemY_euler2Dsgs(this,Mesh,j,dir,pvar,lambda,xvar)
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
    CALL SetEigenValues_euler2Dsgs(this%gamma,pvar(:,j,this%DENSITY), &
         pvar(:,j,this%YVELOCITY), pvar(:,j,this%SGSPRESSURE), pvar(:,j,this%PRESSURE), &
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4),lambda(:,5))
    ! compute characteristic variables
    j1 = j + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
    j2 = MAX(j,j1)
    j1 = MIN(j,j1)
!CDIR IEXPAND
   CALL SetCharVars_euler2Dsgs(this%gamma,pvar(:,j1,this%DENSITY), &
         pvar(:,j2,this%DENSITY),pvar(:,j1,this%YVELOCITY), &
         pvar(:,j2,this%YVELOCITY),pvar(:,j1,this%XVELOCITY), &
         pvar(:,j2,this%XVELOCITY),pvar(:,j1,this%SGSPRESSURE), &
         pvar(:,j2,this%SGSPRESSURE),pvar(:,j1,this%PRESSURE), &
         pvar(:,j2,this%PRESSURE),lambda(:,1),lambda(:,2),lambda(:,3), &
         lambda(:,4),lambda(:,5),xvar(:,1),xvar(:,2),xvar(:,3),xvar(:,4),xvar(:,5))
  END SUBROUTINE CalcCharSystemY_euler2Dsgs


  PURE SUBROUTINE CalcBoundaryDataX_euler2Dsgs(this,Mesh,i1,dir,xvar,pvar)
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
    CALL SetBoundaryData_euler2Dsgs(this%gamma,1.0*SIGN(1,dir),pvar(i1,:,this%DENSITY), &
         pvar(i1,:,this%XVELOCITY),pvar(i1,:,this%YVELOCITY), &
         pvar(i1,:,this%SGSPRESSURE),pvar(i1,:,this%PRESSURE), &
         xvar(:,1),xvar(:,2),xvar(:,3),xvar(:,4),xvar(:,5),pvar(i2,:,this%DENSITY), &
         pvar(i2,:,this%XVELOCITY),pvar(i2,:,this%YVELOCITY),&
         pvar(i2,:,this%SGSPRESSURE),pvar(i2,:,this%PRESSURE))
  END SUBROUTINE CalcBoundaryDataX_euler2Dsgs


  PURE SUBROUTINE CalcBoundaryDataY_euler2Dsgs(this,Mesh,j1,dir,xvar,pvar)
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
    CALL SetBoundaryData_euler2Dsgs(this%gamma,1.0*SIGN(1,dir),pvar(:,j1,this%DENSITY), &
         pvar(:,j1,this%YVELOCITY),pvar(:,j1,this%XVELOCITY), &
         pvar(:,j1,this%SGSPRESSURE),pvar(:,j1,this%PRESSURE), &
         xvar(:,1),xvar(:,2),xvar(:,3),xvar(:,4),xvar(:,5),pvar(:,j2,this%DENSITY), &
         pvar(:,j2,this%YVELOCITY),pvar(:,j2,this%XVELOCITY), &
         pvar(:,j2,this%SGSPRESSURE),pvar(:,j2,this%PRESSURE))
  END SUBROUTINE CalcBoundaryDataY_euler2Dsgs

 PURE SUBROUTINE CalcPrim2RiemannX_euler2Dsgs(this,Mesh,i,pvar,lambda,Rinv)
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
    CALL SetEigenValues_euler2Dsgs(this%gamma,pvar(i,:,this%DENSITY), &
         pvar(i,:,this%XVELOCITY),pvar(i,:,this%SGSPRESSURE),pvar(i,:,this%PRESSURE),&
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4),lambda(:,5))
    ! compute Riemann invariants
!CDIR IEXPAND
    CALL Prim2Riemann_euler2Dsgs(this%gamma,pvar(i,:,this%DENSITY), &
         pvar(i,:,this%XVELOCITY),pvar(i,:,this%YVELOCITY), &
         pvar(i,:,this%SGSPRESSURE),pvar(i,:,this%PRESSURE), &
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4),lambda(:,5), &
         Rinv(:,1),Rinv(:,2),Rinv(:,3),Rinv(:,4),Rinv(:,5))
  END SUBROUTINE CalcPrim2RiemannX_euler2Dsgs


  PURE SUBROUTINE CalcPrim2RiemannY_euler2Dsgs(this,Mesh,j,pvar,lambda,Rinv)
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
    CALL SetEigenValues_euler2Dsgs(this%gamma,pvar(:,j,this%DENSITY), &
         pvar(:,j,this%YVELOCITY),pvar(:,j,this%SGSPRESSURE),pvar(:,j,this%PRESSURE),&
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4),lambda(:,5))
    ! compute Riemann invariants
!CDIR IEXPAND
    CALL Prim2Riemann_euler2Dsgs(this%gamma,pvar(:,j,this%DENSITY), &
         pvar(:,j,this%YVELOCITY),pvar(:,j,this%XVELOCITY), &
         pvar(:,j,this%SGSPRESSURE),pvar(:,j,this%PRESSURE), &
         lambda(:,1),lambda(:,2),lambda(:,3),lambda(:,4),lambda(:,5), &
         Rinv(:,1),Rinv(:,2),Rinv(:,3),Rinv(:,4),Rinv(:,5))
  END SUBROUTINE CalcPrim2RiemannY_euler2Dsgs

  PURE SUBROUTINE CalcRiemann2PrimX_euler2Dsgs(this,Mesh,i,Rinv,pvar)
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
    CALL Riemann2Prim_euler2Dsgs(this%gamma,Rinv(:,1),Rinv(:,2),&
       Rinv(:,3),Rinv(:,4),Rinv(:,5),pvar(i,:,this%Density), &
       pvar(i,:,this%XVELOCITY),pvar(i,:,this%YVELOCITY),&
       pvar(i,:,this%SGSPRESSURE),pvar(i,:,this%PRESSURE))
  END SUBROUTINE CalcRiemann2PrimX_euler2Dsgs

  PURE SUBROUTINE CalcRiemann2PrimY_euler2Dsgs(this,Mesh,j,Rinv,pvar)
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
    CALL Riemann2Prim_euler2Dsgs(this%gamma,Rinv(:,1),Rinv(:,2),&
       Rinv(:,3),Rinv(:,4),Rinv(:,5),pvar(:,j,this%Density),&
       pvar(:,j,this%YVELOCITY),pvar(:,j,this%XVELOCITY),&
       pvar(:,j,this%SGSPRESSURE),pvar(:,j,this%PRESSURE))
  END SUBROUTINE CalcRiemann2PrimY_euler2Dsgs

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
          sterm(i,j,this%DENSITY)   = 0.
          sterm(i,j,this%SGSENERGY) = 0.
          sterm(i,j,this%ENERGY)    = 0.
          ! geometrical source terms in momentum equationes
          sterm(i,j,this%XMOMENTUM) = MomentumSourcesX_euler2Dsgs(&
              cvar(i,j,this%YMOMENTUM),pvar(i,j,this%XVELOCITY),&
              pvar(i,j,this%YVELOCITY),pvar(i,j,this%PRESSURE)+pvar(i,j,this%SGSPRESSURE), &
              Mesh%cxyx(i,j,1),Mesh%cyxy(i,j,1),Mesh%czxz(i,j,1))
          sterm(i,j,this%YMOMENTUM) = MomentumSourcesY_euler2Dsgs(&
              cvar(i,j,this%XMOMENTUM),pvar(i,j,this%XVELOCITY),&
              pvar(i,j,this%YVELOCITY),pvar(i,j,this%PRESSURE)+pvar(i,j,this%SGSPRESSURE), &
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
          sterm(i,j,this%SGSENERGY) = 0.
          sterm(i,j,this%ENERGY)    = 0.
          ! momentum sources (sum up corner values, don't use SUM function,
          ! because it prevents COLLAPSING and causes poor vectorization
          sterm(i,j,this%XMOMENTUM) = MomentumSourcesX_euler2Dsgs(&
              cons(i,j,1,this%YMOMENTUM),prim(i,j,1,this%XVELOCITY),&
              prim(i,j,1,this%YVELOCITY),prim(i,j,1,this%PRESSURE)+prim(i,j,1,this%SGSPRESSURE), &
              Mesh%cxyx(i,j,1),Mesh%cyxy(i,j,1),Mesh%czxz(i,j,1)) &
            + MomentumSourcesX_euler2Dsgs(&
              cons(i,j,2,this%YMOMENTUM),prim(i,j,2,this%XVELOCITY),&
              prim(i,j,2,this%YVELOCITY),prim(i,j,2,this%PRESSURE)+prim(i,j,2,this%SGSPRESSURE), &
              Mesh%cxyx(i,j,2),Mesh%cyxy(i,j,2),Mesh%czxz(i,j,2)) &
            + MomentumSourcesX_euler2Dsgs(&
              cons(i,j,3,this%YMOMENTUM),prim(i,j,3,this%XVELOCITY),&
              prim(i,j,3,this%YVELOCITY),prim(i,j,3,this%PRESSURE)+prim(i,j,3,this%SGSPRESSURE), &
              Mesh%cxyx(i,j,3),Mesh%cyxy(i,j,3),Mesh%czxz(i,j,3)) &
            + MomentumSourcesX_euler2Dsgs(&
              cons(i,j,4,this%YMOMENTUM),prim(i,j,4,this%XVELOCITY),&
              prim(i,j,4,this%YVELOCITY),prim(i,j,4,this%PRESSURE)+prim(i,j,4,this%SGSPRESSURE), &
              Mesh%cxyx(i,j,4),Mesh%cyxy(i,j,4),Mesh%czxz(i,j,4)) 
          sterm(i,j,this%YMOMENTUM) = MomentumSourcesY_euler2Dsgs(&
              cons(i,j,1,this%XMOMENTUM),prim(i,j,1,this%XVELOCITY),&
              prim(i,j,1,this%YVELOCITY),prim(i,j,1,this%PRESSURE)+prim(i,j,1,this%SGSPRESSURE), &
              Mesh%cxyx(i,j,1),Mesh%cyxy(i,j,1),Mesh%czyz(i,j,1)) &
            + MomentumSourcesY_euler2Dsgs(&
              cons(i,j,2,this%XMOMENTUM),prim(i,j,2,this%XVELOCITY),&
              prim(i,j,2,this%YVELOCITY),prim(i,j,2,this%PRESSURE)+prim(i,j,2,this%SGSPRESSURE), &
              Mesh%cxyx(i,j,2),Mesh%cyxy(i,j,2),Mesh%czyz(i,j,2)) &
            + MomentumSourcesY_euler2Dsgs(&
              cons(i,j,3,this%XMOMENTUM),prim(i,j,3,this%XVELOCITY),&
              prim(i,j,3,this%YVELOCITY),prim(i,j,3,this%PRESSURE)+prim(i,j,3,this%SGSPRESSURE), &
              Mesh%cxyx(i,j,3),Mesh%cyxy(i,j,3),Mesh%czyz(i,j,3)) &
            + MomentumSourcesY_euler2Dsgs(&
              cons(i,j,4,this%XMOMENTUM),prim(i,j,4,this%XVELOCITY),&
              prim(i,j,4,this%YVELOCITY),prim(i,j,4,this%PRESSURE)+prim(i,j,4,this%SGSPRESSURE), &
              Mesh%cxyx(i,j,4),Mesh%cyxy(i,j,4),Mesh%czyz(i,j,4))
       END DO
    END DO 
  END SUBROUTINE GeometricalSources_faces


  ! momentum and energy sources due to external force
  PURE SUBROUTINE ExternalSources_euler2Dsgs(this,Mesh,accel,pvar,cvar,sterm)
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
          sterm(i,j,this%SGSENERGY) = 0.
          sterm(i,j,this%ENERGY)    = cvar(i,j,this%XMOMENTUM) * accel(i,j,1) + &
               cvar(i,j,this%YMOMENTUM) * accel(i,j,2)
       END DO
    END DO
  END SUBROUTINE ExternalSources_euler2Dsgs

  PURE SUBROUTINE  ViscositySources_euler2Dsgs(this,Mesh,pvar,btxx,btxy,btyy,sterm)         
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
    CALL ViscositySources_euler2D(this,Mesh,pvar,btxx,btxy,btyy,sterm)
    sterm(:,:,this%SGSENERGY) = 0.0
   
  END SUBROUTINE ViscositySources_euler2Dsgs

  PURE SUBROUTINE  SGSSources_euler2Dsgs(this,Mesh,Sources,pvar,cvar,sterm)         
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

    CALL Divergence(Mesh,Sources%btxx,Sources%btxy,Sources%btxy,Sources%btyy, &
                     sterm(:,:,this%XMOMENTUM),sterm(:,:,this%YMOMENTUM))  

    !(v*tau)_x; amin as temporal storage
    this%amin(:,:) = pvar(:,:,this%XVELOCITY)*Sources%btxx(:,:) &
                    +pvar(:,:,this%YVELOCITY)*Sources%btxy(:,:) 

    !(v*tau)_y; amax as temporal storage
    this%amax(:,:) = pvar(:,:,this%XVELOCITY)*Sources%btxy(:,:) &
                    +pvar(:,:,this%YVELOCITY)*Sources%btyy(:,:) 

    ! div(v*t)
    CALL Divergence(Mesh,this%amin,this%amax,sterm(:,:,this%ENERGY))

    
    !add div(v*t) + rhoepsilon - Sigma 
    sterm(:,:,this%ENERGY) = sterm(:,:,this%ENERGY)+Sources%rhoeps(:,:)-Sources%sigma(:,:)
    
    !add D - rhoepsilon + Sigma 
    sterm(:,:,this%SGSENERGY) = Sources%diff(:,:)-Sources%rhoeps(:,:)+Sources%sigma(:,:) 
    
  END SUBROUTINE SGSSources_euler2Dsgs

  PURE SUBROUTINE CalcSGSTensor_euler2Dsgs(this,Mesh,Sources,C,pvar)
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
    INTENT(INOUT)     :: Sources,this
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

          ! compute the off-diagonal elements the Rate of Strain tensor
          Sources%Sxy(i,j) = &
              (pvar(i+1,j,this%YVELOCITY)-pvar(i-1,j,this%YVELOCITY))*0.25/Mesh%dlx(i,j)&
             +(pvar(i,j+1,this%XVELOCITY)-pvar(i,j-1,this%XVELOCITY))*0.25/Mesh%dly(i,j)&
             - 0.5*Mesh%cxyx(i,j,1)*pvar(i,j,this%XVELOCITY)&
             - 0.5*Mesh%cyxy(i,j,1)*pvar(i,j,this%YVELOCITY)

       END DO
    END DO

    ! compute divergence first and store the result in this%amin
!CDIR IEXPAND
    CALL Divergence(Mesh,pvar(:,:,this%XVELOCITY),pvar(:,:,this%YVELOCITY),this%amin(:,:))

    Sources%btxx(:,:) = C(:,:)*2.0*(Sources%Sxx(:,:) - 1.0/3.0*this%amin(:,:))
    Sources%btyy(:,:) = C(:,:)*2.0*(Sources%Syy(:,:) - 1.0/3.0*this%amin(:,:))   
    Sources%btxy(:,:) = C(:,:)*2.0*Sources%Sxy(:,:)
   
  END SUBROUTINE CalcSGSTensor_euler2Dsgs

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
    CALL Cons2Prim_euler2Dsgs(this%gamma,cvar(i1:i2,j1:j2,this%DENSITY), &
         cvar(i1:i2,j1:j2,this%XMOMENTUM),cvar(i1:i2,j1:j2,this%YMOMENTUM), &
         cvar(i1:i2,j1:j2,this%SGSENERGY),cvar(i1:i2,j1:j2,this%ENERGY), &
         pvar(i1:i2,j1:j2,this%DENSITY),pvar(i1:i2,j1:j2,this%XVELOCITY), &
         pvar(i1:i2,j1:j2,this%YVELOCITY),pvar(i1:i2,j1:j2,this%SGSPRESSURE), &
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
    CALL Cons2Prim_euler2Dsgs(this%gamma,cons(i1:i2,j1:j2,:,this%DENSITY), &
         cons(i1:i2,j1:j2,:,this%XMOMENTUM),cons(i1:i2,j1:j2,:,this%YMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%SGSENERGY),cons(i1:i2,j1:j2,:,this%ENERGY), &
         prim(i1:i2,j1:j2,:,this%DENSITY),prim(i1:i2,j1:j2,:,this%XVELOCITY), &
         prim(i1:i2,j1:j2,:,this%YVELOCITY),prim(i1:i2,j1:j2,:,this%SGSPRESSURE), &
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
    CALL Prim2Cons_euler2Dsgs(this%gamma,pvar(i1:i2,j1:j2,this%DENSITY), &
         pvar(i1:i2,j1:j2,this%XVELOCITY),pvar(i1:i2,j1:j2,this%YVELOCITY), &
         pvar(i1:i2,j1:j2,this%SGSPRESSURE),pvar(i1:i2,j1:j2,this%PRESSURE), &
         cvar(i1:i2,j1:j2,this%DENSITY),cvar(i1:i2,j1:j2,this%XMOMENTUM), &
         cvar(i1:i2,j1:j2,this%YMOMENTUM),cvar(i1:i2,j1:j2,this%SGSENERGY), &
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
    CALL Prim2Cons_euler2Dsgs(this%gamma,prim(i1:i2,j1:j2,:,this%DENSITY), &
         prim(i1:i2,j1:j2,:,this%XVELOCITY),prim(i1:i2,j1:j2,:,this%YVELOCITY), &
         prim(i1:i2,j1:j2,:,this%SGSPRESSURE),prim(i1:i2,j1:j2,:,this%PRESSURE), &
         cons(i1:i2,j1:j2,:,this%DENSITY),cons(i1:i2,j1:j2,:,this%XMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%YMOMENTUM),cons(i1:i2,j1:j2,:,this%SGSENERGY), &
         cons(i1:i2,j1:j2,:,this%ENERGY))
  END SUBROUTINE Convert2Conservative_faces


  PURE SUBROUTINE ReflectionMasks_euler2Dsgs(this,reflX,reflY)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    LOGICAL, DIMENSION(this%VNUM) :: reflX,reflY
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this
    INTENT(OUT)       :: reflX,reflY
    !------------------------------------------------------------------------!
    CALL ReflectionMasks_euler2D(this,reflX,reflY)
    ! western / eastern boundary
    reflX(this%SGSPRESSURE) = .FALSE.
    ! southern / northern boundary
    reflY(this%SGSPRESSURE) = .FALSE.
  END SUBROUTINE ReflectionMasks_euler2Dsgs


  PURE SUBROUTINE AxisMasks_euler2Dsgs(this,reflX,reflY)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    LOGICAL, DIMENSION(this%VNUM) :: reflX,reflY
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this
    INTENT(OUT)       :: reflX,reflY
    !------------------------------------------------------------------------!
    CALL AxisMasks_euler2D(this,reflX,reflY)
    ! western / eastern boundary
    reflX(this%SGSPRESSURE) = .FALSE.
    ! southern / northern boundary
    reflY(this%SGSPRESSURE) = .FALSE.
  END SUBROUTINE AxisMasks_euler2Dsgs

  ELEMENTAL SUBROUTINE SetEigenValues_euler2Dsgs(gamma,rho,v,sgsP,P,l1,l2,l3,l4,l5)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho,v,sgsP,P
    REAL, INTENT(OUT) :: l1,l2,l3,l4,l5
    !------------------------------------------------------------------------!
    REAL :: cs
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL SetEigenValues_euler2D(gamma,rho,v,sgsP+P,l1,l2,l3,l5)
    ! set the missing eigenvalue l4
    l4 = v
  END SUBROUTINE SetEigenValues_euler2Dsgs


  ELEMENTAL SUBROUTINE SetCharVars_euler2Dsgs(gamma,rho1,rho2,u1,u2,v1,v2,sgsP1,sgsP2,&
       P1,P2,l1,l2,l3,l4,l5,xvar1,xvar2,xvar3,xvar4,xvar5)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho1,rho2,u1,u2,v1,v2,sgsP1,sgsP2,P1,P2,l1,l2,l3,l4,l5
    REAL, INTENT(OUT) :: xvar1,xvar2,xvar3,xvar4,xvar5
    !------------------------------------------------------------------------!
    REAL :: gamcs,dlnP,du
    !------------------------------------------------------------------------!
    gamcs= 2.*gamma / (l5-l1) ! = 2*gamma/cs
    dlnP = LOG((P2+sgsP2)/(P1+sgsP1))    ! = LOG(Peff2)-LOG(Peff1)
    du   = u2-u1
    ! characteristic variables
    xvar1 = dlnP - gamcs * du
    xvar2 = -dlnP + gamma * LOG(rho2/rho1)
    xvar3 = gamcs * (v2-v1)
    xvar4 = -dlnP + gamma * LOG(sgsP2/sgsP1)
    xvar5 = dlnP + gamcs * du 
 END SUBROUTINE SetCharVars_euler2Dsgs


  ELEMENTAL SUBROUTINE SetBoundaryData_euler2Dsgs(gamma,dir,rho1,u1,v1,sgsP1,P1, &
       xvar1,xvar2,xvar3,xvar4,xvar5,rho2,u2,v2,sgsP2,P2)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,dir,rho1,u1,v1,sgsP1,P1,xvar1,xvar2,xvar3,xvar4,xvar5
    REAL, INTENT(OUT) :: rho2,u2,v2,sgsP2,P2
    !------------------------------------------------------------------------!
    REAL :: dlnP,csgam
    !------------------------------------------------------------------------!
    dlnP = 0.5 * (xvar5+xvar1)
    ! extrapolate boundary values using characteristic variables
    rho2 = rho1 * EXP(dir*(xvar2+dlnP)/gamma)
    sgsP2 = sgsP1 * EXP(dir*(xvar4+dlnP)/gamma)
    P2   = (P1+sgsP1) * EXP(dir*dlnP) - sgsP2
!CDIR IEXPAND
    csgam= GetSoundSpeed_euler2D(gamma,rho1+rho2,sgsP1+sgsP2+P1+P2) / gamma
    u2   = u1 + dir*csgam * 0.5*(xvar5-xvar1)
    v2   = v1 + dir*csgam * xvar3
  END SUBROUTINE SetBoundaryData_euler2Dsgs

  ELEMENTAL SUBROUTINE Prim2Riemann_euler2Dsgs(gamma,rho,vx,vy,q,p,&
                                       l1,l2,l3,l4,l5,Rminus,Rs,Rq,Rvt,Rplus)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho,vx,vy,p,q,l1,l2,l3,l4,l5
    REAL, INTENT(OUT) :: Rminus,Rs,Rq,Rvt,Rplus
    !------------------------------------------------------------------------!
    REAL :: cs
    !------------------------------------------------------------------------!
    cs = l5-l2 ! l2 = v, l5 = v+cs
    ! compute 1st Riemann invariant (R+)
    Rplus = vx + 2./(gamma-1.0) * cs     
    ! compute 2st Riemann invariant (R-) 
    Rminus = vx - 2./(gamma-1.0) * cs
    ! compute entropy
    Rs = (p + q)/rho**gamma
    ! compute invariant corresponding to psgs
    Rq = (p + q)/q**gamma
    ! tangential velocities
    Rvt = vy   
  END SUBROUTINE Prim2Riemann_euler2Dsgs

  ELEMENTAL SUBROUTINE Riemann2Prim_euler2Dsgs(gamma,Rminus,Rs,Rq,Rvt,Rplus,&
       rho,vx,vy,q,p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,Rminus,Rs,Rq,Rvt,Rplus
    REAL, INTENT(OUT) :: rho,vx,vy,q,p
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
  END SUBROUTINE Riemann2Prim_euler2Dsgs

  ELEMENTAL SUBROUTINE CalcFlux_euler2Dsgs(rho,v,sgsP,P,m1,m2,sgsE,E,f1,f2,f3,f4,f5)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho,v,sgsP,P,m1,m2,sgsE,E
    REAL, INTENT(OUT) :: f1,f2,f3,f4,f5
    !------------------------------------------------------------------------!
    f1 = rho*v
    f2 = m1*v + P + sgsP
    f3 = m2*v
    f4 = sgsE*v
    f5 = (E+P+sgsP)*v
  END SUBROUTINE CalcFlux_euler2Dsgs


  ELEMENTAL SUBROUTINE Cons2Prim_euler2Dsgs(gamma,rho_in,mu,mv,sgsE,E,rho_out,u,v,sgsP,P)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho_in,mu,mv,sgsE,E
    REAL, INTENT(OUT) :: rho_out,u,v,sgsP,P
    !------------------------------------------------------------------------!
    REAL :: inv_rho
    !------------------------------------------------------------------------!
    inv_rho = 1./rho_in
    rho_out = rho_in
    u = mu * inv_rho
    v = mv * inv_rho
    sgsP = sgsE/1.5
    P = (gamma-1.)*(E - 0.5 * inv_rho * (mu*mu+mv*mv))
  END SUBROUTINE Cons2Prim_euler2Dsgs

  
  ELEMENTAL SUBROUTINE Prim2Cons_euler2Dsgs(gamma,rho_in,u,v,sgsP,P,rho_out,mu,mv,sgsE,E)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho_in,u,v,sgsP,P
    REAL, INTENT(OUT) :: rho_out,mu,mv,sgsE,E
    !------------------------------------------------------------------------!
    rho_out = rho_in
    mu = rho_in * u
    mv = rho_in * v
    sgsE = 1.5*sgsP
    E = P/(gamma-1.) + 0.5 * rho_in * (u*u+v*v)
  END SUBROUTINE Prim2Cons_euler2Dsgs
  

  SUBROUTINE ClosePhysics_euler2Dsgs(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL ClosePhysics_euler2D(this)
  END SUBROUTINE ClosePhysics_euler2Dsgs

END MODULE physics_euler2Dsgs
