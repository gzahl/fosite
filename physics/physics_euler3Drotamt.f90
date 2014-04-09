!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler3Drotamt.f90                                         #
!#                                                                           #
!# Copyright (C) 2007-2008                                                   #
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
! basic module for 3D Euler equations with rotational symmetry and
! angular momentum transport
!----------------------------------------------------------------------------!
MODULE physics_euler3Drotamt
  USE physics_common
  USE physics_euler2D, ReflectionMasks_euler3Dra => ReflectionMasks_euler2D
  USE physics_euler3Drotsym, MallocPhysics_euler3Dra => MallocPhysics_euler3Drs, &
       CheckData_euler3Dra => CheckData_euler3Drs, &
       CalculateWaveSpeeds_euler3Dra => CalculateWaveSpeeds_euler3Drs, &
       CalculateFluxesX_euler3Dra => CalculateFluxesX_euler3Drs, &
       CalculateFluxesY_euler3Dra => CalculateFluxesY_euler3Drs, &
       ExternalSources_euler3Dra => ExternalSources_euler3Drs, &
       AxisMasks_euler3Dra => AxisMasks_euler3Drs, &
       ClosePhysics_euler3Dra => ClosePhysics_euler3Drs
  USE mesh_common, ONLY : Mesh_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTERFACE GeometricalSources_euler3Dra
     MODULE PROCEDURE GeometricalSources_center
     MODULE PROCEDURE GeometricalSources_faces
  END INTERFACE
  INTERFACE Convert2Primitive_euler3Dra
     MODULE PROCEDURE Convert2Primitive_center
     MODULE PROCEDURE Convert2Primitive_faces
  END INTERFACE
  INTERFACE Convert2Conservative_euler3Dra
     MODULE PROCEDURE Convert2Conservative_center
     MODULE PROCEDURE Convert2Conservative_faces
  END INTERFACE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: num_var = 5
  REAL, PARAMETER :: TINY = 1.0E-30              ! to avoid division by 0    !
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 3D w/ ang. momentum"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitPhysics_euler3Dra, &
       MallocPhysics_euler3Dra, &
       CheckData_euler3Dra, &
       CalculateWaveSpeeds_euler3Dra, &
       CalculateFluxesX_euler3Dra, &
       CalculateFluxesY_euler3Dra, &
       GeometricalSources_euler3Dra, &
       ExternalSources_euler3Dra, &
       Convert2Primitive_euler3Dra, &
       Convert2Conservative_euler3Dra, &
       ReflectionMasks_euler3Dra, &
       AxisMasks_euler3Dra, &
       ClosePhysics_euler3Dra
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitPhysics_euler3Dra(this,problem)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    INTEGER           :: problem
    !------------------------------------------------------------------------!
    INTENT(IN)        :: problem
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
    this%ZVELOCITY = 4                           ! specific angular momentum !
    this%ZMOMENTUM = 4                           ! angular momentum          !
  END SUBROUTINE InitPhysics_euler3Dra


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
    CALL CentrifugalForces(pvar(:,:,this%DENSITY),pvar(:,:,this%ZVELOCITY), &
         Mesh%bhz(:,:),Mesh%czxz(:,:,1),Mesh%czyz(:,:,1),this%fcent(:,:,1,1), &
         this%fcent(:,:,2,1))

    ! geometrical source terms in momentum equationes
    ! with centrifugal forces
    sterm(:,:,this%XMOMENTUM) = MomentumSourcesX_euler2D(cvar(:,:,this%YMOMENTUM), &
         pvar(:,:,this%XVELOCITY),pvar(:,:,this%YVELOCITY),pvar(:,:,this%PRESSURE), &
         Mesh%cxyx(:,:,1),Mesh%cyxy(:,:,1),Mesh%czxz(:,:,1)) + this%fcent(:,:,1,1)
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
    CALL CentrifugalForces(prim(:,:,:,this%DENSITY),prim(:,:,:,this%ZVELOCITY), &
         Mesh%chz(:,:,:),Mesh%czxz(:,:,:),Mesh%czyz(:,:,:),this%fcent(:,:,:,1), &
         this%fcent(:,:,:,2))

    ! geometrical source terms in momentum equationes
    ! sum up all four corner values
    sterm(:,:,this%XMOMENTUM) = SUM(MomentumSourcesX_euler2D(cons(:,:,:,this%YMOMENTUM), &
         prim(:,:,:,this%XVELOCITY),prim(:,:,:,this%YVELOCITY),prim(:,:,:,this%PRESSURE), &
         Mesh%cxyx(:,:,:),Mesh%cyxy(:,:,:),Mesh%czxz(:,:,:)) + this%fcent(:,:,:,1),DIM=3)
    sterm(:,:,this%YMOMENTUM) = SUM(MomentumSourcesY_euler2D(cons(:,:,:,this%XMOMENTUM), &
         prim(:,:,:,this%XVELOCITY),prim(:,:,:,this%YVELOCITY),prim(:,:,:,this%PRESSURE), &
         Mesh%cxyx(:,:,:),Mesh%cyxy(:,:,:),Mesh%czyz(:,:,:)) + this%fcent(:,:,:,2),DIM=3)

    ! centrifugal force source terms in energy equation
    sterm(:,:,this%ENERGY) = SUM(this%fcent(:,:,:,1) * prim(:,:,:,this%XVELOCITY) &
         + this%fcent(:,:,:,2) * prim(:,:,:,this%YVELOCITY), DIM=3)
  END SUBROUTINE GeometricalSources_faces


  PURE SUBROUTINE Convert2Primitive_center(this,Mesh,cvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
                         :: cvar,pvar
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,cvar
    INTENT(OUT) :: pvar
    !------------------------------------------------------------------------!
    CALL Cons2Prim_euler3Dra(this%gamma,cvar(:,:,this%DENSITY),cvar(:,:,this%XMOMENTUM), &
         cvar(:,:,this%YMOMENTUM),cvar(:,:,this%ZMOMENTUM),cvar(:,:,this%ENERGY), &
         pvar(:,:,this%DENSITY),pvar(:,:,this%XVELOCITY),pvar(:,:,this%YVELOCITY), &
         pvar(:,:,this%ZVELOCITY),pvar(:,:,this%PRESSURE))
  END SUBROUTINE Convert2Primitive_center


  PURE SUBROUTINE Convert2Primitive_faces(this,Mesh,cons,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%vnum) &
                         :: cons,prim
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,cons
    INTENT(OUT) :: prim
    !------------------------------------------------------------------------!
    CALL Cons2Prim_euler3Dra(this%gamma,cons(:,:,:,this%DENSITY),cons(:,:,:,this%XMOMENTUM), &
         cons(:,:,:,this%YMOMENTUM),cons(:,:,:,this%ZMOMENTUM),cons(:,:,:,this%ENERGY), &
         prim(:,:,:,this%DENSITY),prim(:,:,:,this%XVELOCITY),prim(:,:,:,this%YVELOCITY), &
         prim(:,:,:,this%ZVELOCITY),prim(:,:,:,this%PRESSURE))
  END SUBROUTINE Convert2Primitive_faces


  PURE SUBROUTINE Convert2Conservative_center(this,Mesh,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
                         :: cvar,pvar
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,pvar
    INTENT(OUT) :: cvar
    !------------------------------------------------------------------------!
    CALL Prim2Cons_euler3Dra(this%gamma,pvar(:,:,this%DENSITY),pvar(:,:,this%XVELOCITY), &
         pvar(:,:,this%YVELOCITY),pvar(:,:,this%ZVELOCITY),pvar(:,:,this%PRESSURE), &
         cvar(:,:,this%DENSITY),cvar(:,:,this%XMOMENTUM),cvar(:,:,this%YMOMENTUM), &
         cvar(:,:,this%ZMOMENTUM),cvar(:,:,this%ENERGY))
  END SUBROUTINE Convert2Conservative_center


  PURE SUBROUTINE Convert2Conservative_faces(this,Mesh,prim,cons)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%vnum) &
                         :: cons,prim
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,prim
    INTENT(OUT) :: cons
    !------------------------------------------------------------------------!
    CALL Prim2Cons_euler3Dra(this%gamma,prim(:,:,:,this%DENSITY),prim(:,:,:,this%XVELOCITY), &
         prim(:,:,:,this%YVELOCITY),prim(:,:,:,this%ZVELOCITY),prim(:,:,:,this%PRESSURE), &
         cons(:,:,:,this%DENSITY),cons(:,:,:,this%XMOMENTUM),cons(:,:,:,this%YMOMENTUM), &
         cons(:,:,:,this%ZMOMENTUM),cons(:,:,:,this%ENERGY))
  END SUBROUTINE Convert2Conservative_faces


  ELEMENTAL SUBROUTINE CentrifugalForces(rho,l,r,czxz,czyz,fcx,fcy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho,l,r,czxz,czyz
    REAL, INTENT(OUT) :: fcx,fcy
    !------------------------------------------------------------------------!
    REAL :: tmp2
    !------------------------------------------------------------------------!
    tmp2 = rho*(l/(r + TINY))**2
    fcx  = tmp2 * czxz
    fcy  = tmp2 * czyz
  END SUBROUTINE CentrifugalForces


  ELEMENTAL SUBROUTINE Cons2Prim_euler3Dra(gamma,rho_in,mu,mv,L,E,rho_out,u,v,K,P)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho_in,mu,mv,L,E
    REAL, INTENT(OUT) :: rho_out,u,v,K,P
    !------------------------------------------------------------------------!
    REAL :: inv_rho
    !------------------------------------------------------------------------!

    inv_rho = 1./rho_in
    rho_out = rho_in
    u = mu * inv_rho
    v = mv * inv_rho
    K = L * inv_rho
    P = (gamma-1.)*(E - 0.5 * inv_rho * (mu*mu+mv*mv))
  END SUBROUTINE Cons2Prim_euler3Dra

  
  ELEMENTAL SUBROUTINE Prim2Cons_euler3Dra(gamma,rho_in,u,v,K,P,rho_out,mu,mv,L,E)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho_in,u,v,K,P
    REAL, INTENT(OUT) :: rho_out,mu,mv,L,E
    !------------------------------------------------------------------------!

    rho_out = rho_in
    mu = rho_in * u
    mv = rho_in * v
    L  = rho_in * K
    E  = P/(gamma-1.) + 0.5 * rho_in * (u*u+v*v)
  END SUBROUTINE Prim2Cons_euler3Dra

END MODULE physics_euler3Drotamt
