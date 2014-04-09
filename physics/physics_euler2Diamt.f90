!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler2Diamt.f90                                           #
!#                                                                           #
!# Copyright (C) 2012                                                        #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
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
!> \author Manuel Jung
!!
!! \brief basic module for 2D Euler equations with inertial angular momentum
!! transport
!!
!! \extends physics_common
!! \ingroup physics
!----------------------------------------------------------------------------!
MODULE physics_euler2Diamt
  USE physics_common
  USE mesh_common, ONLY : Mesh_TYP
  USE sources_common, ONLY : Sources_TYP
  USE physics_euler2Disoiamt, ONLY : &
       CalcStresses_euler2Dia => CalcStresses_euler2Ditia, &
       SetWaveSpeeds_euler2Dia => SetWaveSpeeds_euler2Ditia, &
       MomentumSourcesX_euler2Dia => MomentumSourcesX_euler2Ditia, &
       CalcWaveSpeeds_euler2Dia => CalcWaveSpeeds_euler2Ditia, &
       SetEigenValues_euler2Ditia, &
       ViscositySources_euler2Ditia, &
       InitPhysics_euler2Ditia, &
       ClosePhysics_euler2Ditia
  USE physics_euler2D, ONLY : &
       CalcCharSystemX_euler2Dia => CalcCharSystemX_euler2D, &
       CalcBoundaryDataX_euler2Dia => CalcBoundaryDataX_euler2D, &
       CalcSoundSpeeds_euler2Dia => CalcSoundSpeeds_euler2D, &
       GetSoundSpeed_euler2Dia => GetSoundSpeed_euler2D
  USE mesh_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE GeometricalSources_euler2Dia
     MODULE PROCEDURE GeometricalSources_center
     MODULE PROCEDURE GeometricalSources_faces
  END INTERFACE
  INTERFACE Convert2Primitive_euler2Dia
     MODULE PROCEDURE Convert2Primitive_center
     MODULE PROCEDURE Convert2Primitive_faces
  END INTERFACE
  INTERFACE Convert2Conservative_euler2Dia
     MODULE PROCEDURE Convert2Conservative_center
     MODULE PROCEDURE Convert2Conservative_faces
  END INTERFACE
  INTERFACE SetPotential_euler2Dia
     MODULE PROCEDURE SetPotential_center
     MODULE PROCEDURE SetPotential_faces
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: num_var = 4              ! number of variables       !
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 2D /w inertial amt"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Physics_TYP, &
       ! methods
       InitPhysics_euler2Dia, &
       ClosePhysics_euler2Dia, &
       CalcWaveSpeeds_euler2Dia, &
       CalcSoundSpeeds_euler2Dia, &
       CalcFluxesX_euler2Dia, &
       CalcFluxesY_euler2Dia, &
       CalcCharSystemX_euler2Dia, &
       CalcBoundaryDataX_euler2Dia, &
       CalcStresses_euler2Dia, &
       GeometricalSources_euler2Dia, &
       ViscositySources_euler2Dia, &
       ExternalSources_euler2Dia, &
       Convert2Primitive_euler2Dia, &
       Convert2Conservative_euler2Dia, &
       ReflectionMasks_euler2Dia, &
       AxisMasks_euler2Dia, &
       GetSoundSpeed_euler2Dia, &
       SetEigenValues_euler2Dia, &
       SetWaveSpeeds_euler2Dia, &
       SetPotential_euler2Dia, &
       MomentumSourcesX_euler2Dia
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitPhysics_euler2Dia(this,Mesh,problem)
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
    CALL InitPhysics_euler2Ditia(this,Mesh,problem,problem_name,num_var)
    
    ! set array indices
    this%PRESSURE  = num_var                           ! pressure            !
    this%ENERGY    = num_var                           ! total energy        !
    ! set names for primitive and conservative variables
    this%pvarname(this%PRESSURE)  = "pressure"
    this%cvarname(this%ENERGY)    = "energy"

    ! allocate memory
    ALLOCATE(this%bphi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX),&
             this%fphi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4),&
             STAT = err)
    ! abort if allocation fails
    IF (err.NE.0) &
         CALL Error(this, "InitPhysics_euler2Dia", "Unable to allocate memory.")
    this%bphi = 0.
    this%fphi = 0.

  END SUBROUTINE InitPhysics_euler2Dia


  PURE SUBROUTINE CalcFluxesX_euler2Dia(this,Mesh,nmin,nmax,prim,cons,xfluxes)
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
    CALL CalcFluxX_euler2Dia(&
         this%hy(:,:,nmin:nmax),&
         prim(:,:,nmin:nmax,this%XVELOCITY),prim(:,:,nmin:nmax,this%PRESSURE), &
         cons(:,:,nmin:nmax,this%XMOMENTUM),cons(:,:,nmin:nmax,this%YMOMENTUM), &
         cons(:,:,nmin:nmax,this%ENERGY),xfluxes(:,:,nmin:nmax,this%DENSITY),&
         xfluxes(:,:,nmin:nmax,this%XMOMENTUM),xfluxes(:,:,nmin:nmax,this%YMOMENTUM), &
         xfluxes(:,:,nmin:nmax,this%ENERGY))
  END SUBROUTINE CalcFluxesX_euler2Dia


  PURE SUBROUTINE CalcFluxesY_euler2Dia(this,Mesh,nmin,nmax,prim,cons,yfluxes)
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
    CALL CalcFluxY_euler2Dia(&
         this%hy(:,:,nmin:nmax),&
         this%w(:,:,nmin:nmax), &
         prim(:,:,nmin:nmax,this%DENSITY), &
         prim(:,:,nmin:nmax,this%YVELOCITY),prim(:,:,nmin:nmax,this%PRESSURE), &
         cons(:,:,nmin:nmax,this%XMOMENTUM),cons(:,:,nmin:nmax,this%YMOMENTUM), &
         cons(:,:,nmin:nmax,this%ENERGY),yfluxes(:,:,nmin:nmax,this%DENSITY), &
         yfluxes(:,:,nmin:nmax,this%XMOMENTUM),yfluxes(:,:,nmin:nmax,this%YMOMENTUM), &
         yfluxes(:,:,nmin:nmax,this%ENERGY))
  END SUBROUTINE CalcFluxesY_euler2Dia


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
         
          ! geometrical source terms in momentum equations 
          sterm(i,j,this%XMOMENTUM) = MomentumSourcesX_euler2Dia(&
              cvar(i,j,this%DENSITY),&
              cvar(i,j,this%YMOMENTUM),&
              pvar(i,j,this%YVELOCITY), &
              pvar(i,j,this%PRESSURE), &
              Mesh%bhy(i,j), &
              Mesh%cyxy(i,j,1))
          
          ! no geometrical inertial angular momentum sources
          sterm(i,j,this%YMOMENTUM) = 0.0

          sterm(i,j,this%ENERGY) = EnergySources_euler2Dia(&
              this%Omega, &
              pvar(i,j,this%DENSITY), &
              pvar(i,j,this%XVELOCITY), &
              Mesh%bhy(i,j), &
              Mesh%cyxy(i,j,1))
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
          ! momentum sources (sum up corner values, don't use SUM function,
          ! because it prevents COLLAPSING and causes poor vectorization

          ! momentum sources (sum up corner values, don't use SUM function,
          ! because it prevents COLLAPSING and causes poor vectorization
          sterm(i,j,this%XMOMENTUM) = MomentumSourcesX_euler2Dia(&
              cons(i,j,1,this%DENSITY),&
              cons(i,j,1,this%YMOMENTUM),&
              prim(i,j,1,this%YVELOCITY), &
              prim(i,j,1,this%PRESSURE),&
              this%hy(i,j,1), &
              Mesh%cyxy(i,j,1)) &
            + MomentumSourcesX_euler2Dia(&
              cons(i,j,2,this%DENSITY),&
              cons(i,j,2,this%YMOMENTUM),&
              prim(i,j,2,this%YVELOCITY), &
              prim(i,j,2,this%PRESSURE),&
              this%hy(i,j,2), &
              Mesh%cyxy(i,j,2)) &
            + MomentumSourcesX_euler2Dia(&
              cons(i,j,3,this%DENSITY),&
              cons(i,j,3,this%YMOMENTUM),&
              prim(i,j,3,this%YVELOCITY), &
              prim(i,j,3,this%PRESSURE),&
              this%hy(i,j,3), &
              Mesh%cyxy(i,j,3)) &
            + MomentumSourcesX_euler2Dia(&
              cons(i,j,4,this%DENSITY),&
              cons(i,j,4,this%YMOMENTUM),&
              prim(i,j,4,this%YVELOCITY), &
              prim(i,j,4,this%PRESSURE),&
              this%hy(i,j,4), &
              Mesh%cyxy(i,j,4))

          ! no geometrical inertial angular momentum sources
          sterm(i,j,this%YMOMENTUM) = 0.

          sterm(i,j,this%ENERGY) = &
              EnergySources_euler2Dia(&
              this%Omega, &
              prim(i,j,1,this%DENSITY), &
              prim(i,j,1,this%XVELOCITY), &
              this%hy(i,j,1), &
              Mesh%cyxy(i,j,1)) &
            + EnergySources_euler2Dia(&
              this%Omega, &
              prim(i,j,2,this%DENSITY), &
              prim(i,j,2,this%XVELOCITY), &
              this%hy(i,j,2), &
              Mesh%cyxy(i,j,2)) &
            + EnergySources_euler2Dia(&
              this%Omega, &
              prim(i,j,3,this%DENSITY), &
              prim(i,j,3,this%XVELOCITY), &
              this%hy(i,j,3), &
              Mesh%cyxy(i,j,3)) &
            + EnergySources_euler2Dia(&
              this%Omega, &
              prim(i,j,4,this%DENSITY), &
              prim(i,j,4,this%XVELOCITY), &
              this%hy(i,j,4), &
              Mesh%cyxy(i,j,4))
       END DO
    END DO 
  END SUBROUTINE GeometricalSources_faces


  ! momentum and energy sources due to external force
  PURE SUBROUTINE ExternalSources_euler2Dia(this,Mesh,accel,pvar,cvar,sterm)
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
          sterm(i,j,this%YMOMENTUM) = Mesh%bhy(i,j) &
                                      * pvar(i,j,this%DENSITY) * accel(i,j,2)
          sterm(i,j,this%ENERGY)    = cvar(i,j,this%XMOMENTUM) * accel(i,j,1) &
            + pvar(i,j,this%YVELOCITY) * pvar(i,j,this%DENSITY) * accel(i,j,2)
       END DO
    END DO
  END SUBROUTINE ExternalSources_euler2Dia


  PURE SUBROUTINE ViscositySources_euler2Dia(this,Mesh,pvar,btxx,btxy,btyy,sterm)         
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
    CALL ViscositySources_euler2Ditia(this,Mesh,pvar,btxx,btxy,btyy,sterm)
 
    !compute scalar product of v and tau (x-component)
    this%amin(:,:)  = pvar(:,:,this%XVELOCITY)*btxx(:,:) &
                    + pvar(:,:,this%YVELOCITY)*btxy(:,:) 

    !compute scalar product of v and tau (y-component)
    this%amax(:,:) = pvar(:,:,this%XVELOCITY)*btxy(:,:) &
                   + pvar(:,:,this%YVELOCITY)*btyy(:,:)

    ! compute vector divergence of scalar product v and tau
!CDIR IEXPAND
    CALL Divergence(Mesh,this%amin(:,:),this%amax(:,:),sterm(:,:,this%ENERGY))
  END SUBROUTINE ViscositySources_euler2Dia


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
    CALL Cons2Prim_euler2Dia(&
         this%Omega, &
         Mesh%bhy(i1:i2,j1:j2), &
         this%bphi(i1:i2,j1:j2), &
         this%gamma,&
         cvar(i1:i2,j1:j2,this%DENSITY), &
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
    CALL Cons2Prim_euler2Dia(&
         this%Omega, &
         this%hy(i1:i2,j1:j2,:), &
         this%fphi(i1:i2,j1:j2,:), &
         this%gamma,&
         cons(i1:i2,j1:j2,:,this%DENSITY), &
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
    CALL Prim2Cons_euler2Dia(&
         this%Omega, &
         Mesh%bhy(i1:i2,j1:j2), &
         this%bphi(i1:i2,j1:j2), &
         this%gamma,&
         pvar(i1:i2,j1:j2,this%DENSITY), &
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
    CALL Prim2Cons_euler2Dia(&
         this%Omega, &
         this%hy(i1:i2,j1:j2,:), &
         this%fphi(i1:i2,j1:j2,:), &
         this%gamma,&
         prim(i1:i2,j1:j2,:,this%DENSITY), &
         prim(i1:i2,j1:j2,:,this%XVELOCITY),prim(i1:i2,j1:j2,:,this%YVELOCITY), &
         prim(i1:i2,j1:j2,:,this%PRESSURE),cons(i1:i2,j1:j2,:,this%DENSITY), &
         cons(i1:i2,j1:j2,:,this%XMOMENTUM),cons(i1:i2,j1:j2,:,this%YMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%ENERGY))
  END SUBROUTINE Convert2Conservative_faces


  PURE SUBROUTINE ReflectionMasks_euler2Dia(this,reflX,reflY)
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
  END SUBROUTINE ReflectionMasks_euler2Dia


  PURE SUBROUTINE AxisMasks_euler2Dia(this,reflX,reflY)
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
  END SUBROUTINE AxisMasks_euler2Dia


  ELEMENTAL SUBROUTINE SetEigenValues_euler2Dia(gamma,rho,v,P,l1,l2,l3,l4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho,v,P
    REAL, INTENT(OUT) :: l1,l2,l3,l4
    !------------------------------------------------------------------------!
    REAL :: cs
    !------------------------------------------------------------------------!
    ! adiabatic sound speed
!CDIR IEXPAND
    cs = GetSoundSpeed_euler2Dia(gamma,rho,P)
    ! call subroutine for isothermal case with the adiabatic sound speed
!CDIR IEXPAND
    CALL SetEigenValues_euler2Ditia(cs,v,l1,l2,l4)
    ! set the missing eigenvalue l3
    l3 = v
  END SUBROUTINE SetEigenValues_euler2Dia


  ELEMENTAL FUNCTION EnergySources_euler2Dia(Omega,rho,vx,hy,cyxy) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: Omega,rho,vx,hy,cyxy
    REAL :: st
    !------------------------------------------------------------------------!
    st = 0.0
  END FUNCTION EnergySources_euler2Dia


  ELEMENTAL SUBROUTINE CalcFluxX_euler2Dia(bhy,vx,P,mx,my,E,f1,f2,f3,f4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: bhy,vx,P,mx,my,E
    REAL, INTENT(OUT) :: f1, f2, f3, f4
    !------------------------------------------------------------------------!
    f1 = mx
    f2 = mx*vx + P
    f3 = my*vx
    f4 = (E+P)*vx
  END SUBROUTINE CalcFluxX_euler2Dia


  ELEMENTAL SUBROUTINE CalcFluxY_euler2Dia(hy,w,rho,vy,P,mx,my,E,g1,g2,g3,g4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: hy,rho,vy,P,mx,my,E,w
    REAL, INTENT(OUT) :: g1, g2, g3, g4
    !------------------------------------------------------------------------!
    REAL              :: v
    !------------------------------------------------------------------------!
    v = vy - hy*w
    g1 = rho * v
    g2 = mx  * v
    g3 = my  * v + hy*P
    g4 = E   * v + P * vy
  END SUBROUTINE CalcFluxY_euler2Dia


  ELEMENTAL SUBROUTINE Cons2Prim_euler2Dia(Omega,hy,phi,gamma,rho_in,mu,mv,E,rho_out,u,v,P)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: Omega,hy,phi,gamma,rho_in,mu,mv,E
    REAL, INTENT(OUT) :: rho_out,u,v,P
    !------------------------------------------------------------------------!
    REAL :: inv_rho
    !------------------------------------------------------------------------!
    inv_rho = 1./rho_in
    rho_out = rho_in
    u = mu * inv_rho
    v = mv * inv_rho / hy - hy * Omega
    P = (gamma-1.)*(E - 0.5 * rho_in * (u*u+v*v-(hy*Omega)**2) - rho_in *phi )
  END SUBROUTINE Cons2Prim_euler2Dia

  
  ELEMENTAL SUBROUTINE Prim2Cons_euler2Dia(Omega,hy,phi,gamma,rho_in,u,v,P,rho_out,mu,mv,E)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: Omega,hy,phi,gamma,rho_in,u,v,P
    REAL, INTENT(OUT) :: rho_out,mu,mv,E
    !------------------------------------------------------------------------!
    rho_out = rho_in
    mu = rho_in * u
    mv = rho_in * hy * (v + hy*Omega)
    E = P/(gamma-1.) + 0.5 * rho_in * (u*u+v*v-(hy*Omega)**2) + rho_in * phi
  END SUBROUTINE Prim2Cons_euler2Dia


  PURE SUBROUTINE SetPotential_center(this,Mesh,bphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: bphi
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,bphi
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    this%bphi(:,:) = bphi(:,:)
  END SUBROUTINE SetPotential_center


  PURE SUBROUTINE SetPotential_faces(this,Mesh,fphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4) &
                      :: fphi
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,fphi
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    this%fphi(:,:,:) = fphi(:,:,:)
  END SUBROUTINE SetPotential_faces


  SUBROUTINE ClosePhysics_euler2Dia(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%bphi,this%fphi)
    CALL ClosePhysics_euler2Ditia(this)
  END SUBROUTINE ClosePhysics_euler2Dia

END MODULE physics_euler2Diamt
