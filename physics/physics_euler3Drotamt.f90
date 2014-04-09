!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler3Drotamt.f90                                         #
!#                                                                           #
!# Copyright (C) 2007 Tobias Illenseer <tillense@ita.uni-heidelberg.de>      #
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
  USE physics_euler2D, GetReflectionMasks_euler3Drotamt => GetReflectionMasks_euler2D
  USE physics_euler3Drotsym, MallocPhysics_euler3Drotamt => MallocPhysics_euler3Drotsym, &
       CheckData_euler3Drotamt => CheckData_euler3Drotsym, &
       CalculateWaveSpeeds_euler3Drotamt => CalculateWaveSpeeds_euler3Drotsym, &
       CalculateFluxesX_euler3Drotamt => CalculateFluxesX_euler3Drotsym, &
       CalculateFluxesY_euler3Drotamt => CalculateFluxesY_euler3Drotsym, &
       ExternalSources_euler3Drotamt => ExternalSources_euler3Drotsym, &
       GetAxisMasks_euler3Drotamt => GetAxisMasks_euler3Drotsym, &
       ClosePhysics_euler3Drotamt => ClosePhysics_euler3Drotsym
  USE mesh_common, ONLY : Mesh_TYP
  USE output_common, ONLY : Datastruct_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTERFACE GeometricalSources_euler3Drotamt
     MODULE PROCEDURE GeometricalSources_center
     MODULE PROCEDURE GeometricalSources_faces
  END INTERFACE
  INTERFACE Convert2Primitive_euler3Drotamt
     MODULE PROCEDURE Convert2Primitive_center
     MODULE PROCEDURE Convert2Primitive_faces
  END INTERFACE
  INTERFACE Convert2Conservative_euler3Drotamt
     MODULE PROCEDURE Convert2Conservative_center
     MODULE PROCEDURE Convert2Conservative_faces
  END INTERFACE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: DEN = 1                  ! mass density              !
  INTEGER, PARAMETER :: V_X = 2                  ! velocity components       !
  INTEGER, PARAMETER :: V_Y = 3
  INTEGER, PARAMETER :: l_Z = 4                  ! specific angular momentum !
  INTEGER, PARAMETER :: PRE = 5                  ! pressure                  !
  INTEGER, PARAMETER :: M_X = 2                  ! momentum components       !
  INTEGER, PARAMETER :: M_Y = 3
  INTEGER, PARAMETER :: M_Z = 4                  ! angular momentum          !
  INTEGER, PARAMETER :: ENE = 5                  ! total energy w/o rot. en. !
  INTEGER, PARAMETER :: num_var = 5
  REAL, PARAMETER :: TINY = 1.0E-30              ! to avoid division by 0    !
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 3D w/ ang. momentum"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitPhysics_euler3Drotamt, &
       MallocPhysics_euler3Drotamt, &
       CheckData_euler3Drotamt, &
       CalculateWaveSpeeds_euler3Drotamt, &
       CalculateFluxesX_euler3Drotamt, &
       CalculateFluxesY_euler3Drotamt, &
       GeometricalSources_euler3Drotamt, &
       ExternalSources_euler3Drotamt, &
       Convert2Primitive_euler3Drotamt, &
       Convert2Conservative_euler3Drotamt, &
       GetDataStruct_euler3Drotamt, &
       GetReflectionMasks_euler3Drotamt, &
       GetAxisMasks_euler3Drotamt, &
       ClosePhysics_euler3Drotamt
  !--------------------------------------------------------------------------!

CONTAINS

  PURE SUBROUTINE InitPhysics_euler3Drotamt(this,problem,gamma,dpmax)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    INTEGER           :: problem
    REAL              :: gamma,dpmax
    !------------------------------------------------------------------------!
    INTENT(IN)        :: problem,gamma,dpmax
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL InitPhysics(this,problem,problem_name,num_var)
    this%gamma  = gamma
    this%dpmax  = dpmax
  END SUBROUTINE InitPhysics_euler3Drotamt


  SUBROUTINE GetDataStruct_euler3Drotamt(this,ds)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Datastruct_TYP), DIMENSION(:), POINTER :: ds
    !------------------------------------------------------------------------!
    INTEGER :: err
    !------------------------------------------------------------------------!
    INTENT(IN) :: this
    !------------------------------------------------------------------------!

    ! 4 entries one for density, x-y-velocity, specific
    ! angular momentum and pressure respectively
    ALLOCATE(ds(4),STAT=err)
    IF (err.NE.0) THEN
       PRINT *, "ERROR in GetDataStruct_euler3Drotamt: Unable to allocate memory "
       STOP
    END IF

    ds(1)%name  = 'density'
    ds(1)%rank  = 0
    ds(1)%shape = 1
    ds(2)%name = 'velocity'
    ds(2)%rank  = 1
    ds(2)%shape = 2
    ds(3)%name  = 'specific angular momentum'
    ds(3)%rank  = 0
    ds(3)%shape = 1
    ds(4)%name  = 'pressure'
    ds(4)%rank  = 0
    ds(4)%shape = 1
  END SUBROUTINE GetDataStruct_euler3Drotamt


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
    CALL CentrifugalForces(pvar(:,:,DEN),pvar(:,:,l_Z),Mesh%bhz(:,:), &
         Mesh%czxz(:,:,1),Mesh%czyz(:,:,1),this%fcent(:,:,1,1),this%fcent(:,:,2,1))

    ! geometrical source terms in momentum equationes
    ! with centrifugal forces
    sterm(:,:,M_X) = MomentumSourcesX_euler2D(cvar(:,:,M_Y),pvar(:,:,V_X), &
         pvar(:,:,V_Y),pvar(:,:,PRE),Mesh%cxyx(:,:,1), &
         Mesh%cyxy(:,:,1),Mesh%czxz(:,:,1)) + this%fcent(:,:,1,1)
    sterm(:,:,M_Y) = MomentumSourcesY_euler2D(cvar(:,:,M_X),pvar(:,:,V_X), &
         pvar(:,:,V_Y),pvar(:,:,PRE),Mesh%cxyx(:,:,1), &
         Mesh%cyxy(:,:,1),Mesh%czyz(:,:,1)) + this%fcent(:,:,2,1)

    ! centrifugal force source terms in energy equation
    sterm(:,:,ENE) = this%fcent(:,:,1,1) * pvar(:,:,V_X) &
         + this%fcent(:,:,2,1) * pvar(:,:,V_Y)
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
    CALL CentrifugalForces(prim(:,:,:,DEN),prim(:,:,:,l_Z),Mesh%chz(:,:,:), &
         Mesh%czxz(:,:,:),Mesh%czyz(:,:,:),this%fcent(:,:,:,1),this%fcent(:,:,:,2))

    ! geometrical source terms in momentum equationes
    ! sum up all four corner values
    sterm(:,:,M_X) = SUM(MomentumSourcesX_euler2D(cons(:,:,:,M_Y), &
         prim(:,:,:,V_X),prim(:,:,:,V_Y),prim(:,:,:,PRE),Mesh%cxyx(:,:,:), &
         Mesh%cyxy(:,:,:),Mesh%czxz(:,:,:)) + this%fcent(:,:,:,1),DIM=3)
    sterm(:,:,M_Y) = SUM(MomentumSourcesY_euler2D(cons(:,:,:,M_X), &
         prim(:,:,:,V_X),prim(:,:,:,V_Y),prim(:,:,:,PRE),Mesh%cxyx(:,:,:), &
         Mesh%cyxy(:,:,:),Mesh%czyz(:,:,:)) + this%fcent(:,:,:,2),DIM=3)

    ! centrifugal force source terms in energy equation
    sterm(:,:,ENE) = SUM(this%fcent(:,:,:,1) * prim(:,:,:,V_X) &
         + this%fcent(:,:,:,2) * prim(:,:,:,V_Y), DIM=3)
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
    CALL Cons2Prim_euler3Drotamt(this%gamma,cvar(:,:,DEN),cvar(:,:,M_X), &
         cvar(:,:,M_Y),cvar(:,:,M_Z),cvar(:,:,ENE),pvar(:,:,DEN),pvar(:,:,V_X), &
         pvar(:,:,V_Y),pvar(:,:,l_Z),pvar(:,:,PRE))
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
    CALL Cons2Prim_euler3Drotamt(this%gamma,cons(:,:,:,DEN),cons(:,:,:,M_X), &
         cons(:,:,:,M_Y),cons(:,:,:,M_Z),cons(:,:,:,ENE),prim(:,:,:,DEN), &
         prim(:,:,:,V_X),prim(:,:,:,V_Y),prim(:,:,:,l_Z),prim(:,:,:,PRE))
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
    CALL Prim2Cons_euler3Drotamt(this%gamma,pvar(:,:,DEN),pvar(:,:,V_X), &
         pvar(:,:,V_Y),pvar(:,:,l_Z),pvar(:,:,PRE),cvar(:,:,DEN),cvar(:,:,M_X), &
         cvar(:,:,M_Y),cvar(:,:,M_Z),cvar(:,:,ENE))
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
    CALL Prim2Cons_euler3Drotamt(this%gamma,prim(:,:,:,DEN),prim(:,:,:,V_X), &
         prim(:,:,:,V_Y),prim(:,:,:,l_Z),prim(:,:,:,PRE),cons(:,:,:,DEN), &
         cons(:,:,:,M_X),cons(:,:,:,M_Y),cons(:,:,:,M_Z),cons(:,:,:,ENE))
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


  ELEMENTAL SUBROUTINE Cons2Prim_euler3Drotamt(gamma,rho_in,mu,mv,L,E,rho_out,u,v,K,P)
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
  END SUBROUTINE Cons2Prim_euler3Drotamt

  
  ELEMENTAL SUBROUTINE Prim2Cons_euler3Drotamt(gamma,rho_in,u,v,K,P,rho_out,mu,mv,L,E)
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
  END SUBROUTINE Prim2Cons_euler3Drotamt

END MODULE physics_euler3Drotamt
