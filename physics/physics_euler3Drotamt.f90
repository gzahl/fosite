!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler3Drotamt.f90                                         #
!#                                                                           #
!# Copyright (C) 2007-2010                                                   #
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
! basic module for 3D Euler equations with rotational symmetry and
! angular momentum transport
!----------------------------------------------------------------------------!
MODULE physics_euler3Drotamt
  USE physics_common
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_euler2D, ONLY : MomentumSourcesX_euler2D, MomentumSourcesY_euler2D
  USE physics_euler3Drotsym, &
       MallocPhysics_euler3Dra => MallocPhysics_euler3Drs, &
       CheckData_euler3Dra => CheckData_euler3Drs, &
       CalculateWaveSpeeds_euler3Dra => CalculateWaveSpeeds_euler3Drs, &
       CalculateFluxesX_euler3Dra => CalculateFluxesX_euler3Drs, &
       CalculateFluxesY_euler3Dra => CalculateFluxesY_euler3Drs, &
       ExternalSources_euler3Dra => ExternalSources_euler3Drs, &
       ReflectionMasks_euler3Dra => ReflectionMasks_euler3Drs, &
       AxisMasks_euler3Dra => AxisMasks_euler3Drs, &
       ClosePhysics_euler3Dra => ClosePhysics_euler3Drs
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
    INTEGER           :: err
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
    ! set names for primitive and conservative variables
    this%pvarname(this%DENSITY)   = "density"
    this%pvarname(this%XVELOCITY) = "x-velocity"
    this%pvarname(this%YVELOCITY) = "y-velocity"
    this%pvarname(this%ZVELOCITY) = "specang-momentum"
    this%pvarname(this%PRESSURE)  = "pressure"
    this%cvarname(this%DENSITY)   = "density"
    this%cvarname(this%XMOMENTUM) = "x-momentum"
    this%cvarname(this%YMOMENTUM) = "y-momentum"
    this%cvarname(this%ZMOMENTUM) = "ang-momentum"
    this%cvarname(this%ENERGY)    = "energy"

    ALLOCATE(this%structure(4),this%errormap(0:3),STAT = err)
    ! abort if allocation fails
    IF (err.NE.0) &
         CALL Error(this, "InitPhysics_euler3Drotamt", "Unable to allocate memory.")
    this%nstruc = 4
    this%structure(1)%name = "coordinates"
    this%structure(1)%pos = -1
    this%structure(1)%rank = 1
    this%structure(1)%dim  = 2
    this%structure(2)%name = "density"
    this%structure(2)%pos  = this%DENSITY
    this%structure(2)%rank = 0
    this%structure(2)%dim  = 1
    this%structure(3)%name = "velocity"
    this%structure(3)%pos  = this%XVELOCITY
    this%structure(3)%rank = 1
    this%structure(3)%dim  = 3
    this%structure(4)%name = "pressure"
    this%structure(4)%pos  = this%PRESSURE
    this%structure(4)%rank = 0
    this%structure(4)%dim  = 1

    !set errormapping (CheckData Symbols)
    this%errormap(0) = 'X'
    this%errormap(1) = 'D'
    this%errormap(2) = 'P'
    this%errormap(3) = 'B'
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

    ! no geometrical sources in continuity and angular momentum equation
    sterm(:,:,this%DENSITY)   = 0.0
    sterm(:,:,this%ZMOMENTUM) = 0.0

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

    ! no geometrical sources in continuity and angular momentum equation
    sterm(:,:,this%DENSITY)   = 0.0
    sterm(:,:,this%ZMOMENTUM) = 0.0

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
    CALL Cons2Prim_euler3Dra(this%gamma,cvar(i1:i2,j1:j2,this%DENSITY), &
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
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%vnum) &
                         :: cons,prim
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,i1,i2,j1,j2,cons
    INTENT(OUT) :: prim
    !------------------------------------------------------------------------!
    CALL Cons2Prim_euler3Dra(this%gamma,cons(i1:i2,j1:j2,:,this%DENSITY), &
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
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
                         :: cvar,pvar
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,i1,i2,j1,j2,pvar
    INTENT(OUT) :: cvar
    !------------------------------------------------------------------------!
    CALL Prim2Cons_euler3Dra(this%gamma,pvar(i1:i2,j1:j2,this%DENSITY), &
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
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%vnum) &
                         :: cons,prim
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,i1,i2,j1,j2,prim
    INTENT(OUT) :: cons
    !------------------------------------------------------------------------!
    CALL Prim2Cons_euler3Dra(this%gamma,prim(i1:i2,j1:j2,:,this%DENSITY), &
         prim(i1:i2,j1:j2,:,this%XVELOCITY),prim(i1:i2,j1:j2,:,this%YVELOCITY), &
         prim(i1:i2,j1:j2,:,this%ZVELOCITY),prim(i1:i2,j1:j2,:,this%PRESSURE), &
         cons(i1:i2,j1:j2,:,this%DENSITY),cons(i1:i2,j1:j2,:,this%XMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%YMOMENTUM),cons(i1:i2,j1:j2,:,this%ZMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%ENERGY))
  END SUBROUTINE Convert2Conservative_faces


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
