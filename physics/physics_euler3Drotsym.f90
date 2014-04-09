!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler3Drotsym.f90                                         #
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
! basic module for 3D Euler equations with rotational symmetry
!----------------------------------------------------------------------------!
MODULE physics_euler3Drotsym
  USE physics_common
  USE physics_euler2D, GetReflectionMasks_euler3Drotsym => GetReflectionMasks_euler2D
  USE mesh_common, ONLY : Mesh_TYP
  USE output_common, ONLY : Datastruct_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTERFACE CalculateWaveSpeeds_euler3Drotsym
     MODULE PROCEDURE CalculateWaveSpeeds_center
     MODULE PROCEDURE CalculateWaveSpeeds_faces
  END INTERFACE
  INTERFACE GeometricalSources_euler3Drotsym
     MODULE PROCEDURE GeometricalSources_center
     MODULE PROCEDURE GeometricalSources_faces
  END INTERFACE
  INTERFACE Convert2Primitive_euler3Drotsym
     MODULE PROCEDURE Convert2Primitive_center
     MODULE PROCEDURE Convert2Primitive_faces
  END INTERFACE
  INTERFACE Convert2Conservative_euler3Drotsym
     MODULE PROCEDURE Convert2Conservative_center
     MODULE PROCEDURE Convert2Conservative_faces
  END INTERFACE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: DEN = 1                        ! mass density        !
  INTEGER, PARAMETER :: V_X = 2                        ! velocity components !
  INTEGER, PARAMETER :: V_Y = 3
  INTEGER, PARAMETER :: V_Z = 4
  INTEGER, PARAMETER :: PRE = 5                        ! pressure            !
  INTEGER, PARAMETER :: M_X = 2                        ! momentum components !
  INTEGER, PARAMETER :: M_Y = 3
  INTEGER, PARAMETER :: M_Z = 4
  INTEGER, PARAMETER :: ENE = 5                        ! total energy        !
  INTEGER, PARAMETER :: num_var = 5
  REAL, PARAMETER :: TINY = 1.0E-30              ! to avoid division by 0    !
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 3D w/ rot. symmetry"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitPhysics_euler3Drotsym, &
       MallocPhysics_euler3Drotsym, &
       CheckData_euler3Drotsym, &
       CalculateWaveSpeeds_euler3Drotsym, &
       CalculateFluxesX_euler3Drotsym, &
       CalculateFluxesY_euler3Drotsym, &
       GeometricalSources_euler3Drotsym, &
       ExternalSources_euler3Drotsym, &
       Convert2Primitive_euler3Drotsym, &
       Convert2Conservative_euler3Drotsym, &
       GetDataStruct_euler3Drotsym, &
       GetReflectionMasks_euler3Drotsym, &
       GetAxisMasks_euler3Drotsym, &
       ClosePhysics_euler3Drotsym
  !--------------------------------------------------------------------------!

CONTAINS

  PURE SUBROUTINE InitPhysics_euler3Drotsym(this,problem,gamma,dpmax)
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
  END SUBROUTINE InitPhysics_euler3Drotsym


  SUBROUTINE MallocPhysics_euler3Drotsym(this,Mesh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! allocate memory for arrays used in Euler3D w/ rotational symmetry
    CALL MallocPhysics_euler2D(this,Mesh)
    ALLOCATE(this%fcent(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,2), &
         STAT = err)
    IF (err.NE.0) THEN
       PRINT *, "ERROR in MallocPhysics_euler3Drotsym: Can't allocate memory!"
       STOP
    END IF
  END SUBROUTINE MallocPhysics_euler3Drotsym


  SUBROUTINE GetDataStruct_euler3Drotsym(this,ds)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Datastruct_TYP), DIMENSION(:), POINTER :: ds
    !------------------------------------------------------------------------!
    INTEGER :: err
    !------------------------------------------------------------------------!
    INTENT(IN) :: this
    !------------------------------------------------------------------------!

    ! 3 entries one for density, velocity and pressure respectively
    ALLOCATE(ds(3),STAT=err)
    IF (err.NE.0) THEN
       PRINT *, "ERROR in GetDataStruct_euler3Drotsym: Unable to allocate memory "
       STOP
    END IF

    ds(1)%name  = 'density'
    ds(1)%rank  = 0
    ds(1)%shape = 1
    ds(2)%name = 'velocity'
    ds(2)%rank  = 1
    ds(2)%shape = 3
    ds(3)%name  = 'pressure'
    ds(3)%rank  = 0
    ds(3)%shape = 1
  END SUBROUTINE GetDataStruct_euler3Drotsym


  PURE FUNCTION CheckData_euler3Drotsym(this,Mesh,pvar,pold) RESULT(bad_data)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
         :: pvar,pold
    LOGICAL           :: bad_data
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,pvar,pold
    !------------------------------------------------------------------------!
    IF (MINVAL(pvar(:,:,PRE)).LE.0.) THEN
       bad_data = .TRUE.
    ELSE IF (MAXVAL(ABS(1.0-pvar(:,:,PRE)/pold(:,:,PRE))).GT.this%dpmax) THEN
       bad_data = .TRUE.
    ELSE
       bad_data = .FALSE.
    END IF
  END FUNCTION CheckData_euler3Drotsym


  PURE SUBROUTINE CalculateWaveSpeeds_center(this,Mesh,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
         :: pvar
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! sound speed
    this%csound(:,:,1) = GetSoundSpeed_adiabatic(this%gamma,pvar(:,:,DEN), &
         pvar(:,:,PRE))
    
    ! minimal and maximal wave speeds at cell interfaces
    ! x-direction
    CALL SetWaveSpeeds_euler(pvar(:,:,V_X),this%csound(:,:,1), &
         this%tmin(:,:,1),this%tmax(:,:,1))
    ! y-direction
    CALL SetWaveSpeeds_euler(pvar(:,:,V_Y),this%csound(:,:,1), &
         this%tmin(:,:,2),this%tmax(:,:,2))
  END SUBROUTINE CalculateWaveSpeeds_center


  PURE SUBROUTINE CalculateWaveSpeeds_faces(this,Mesh,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%vnum) &
         :: prim
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,prim
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! sound speed
    this%csound(:,:,:) = GetSoundSpeed_adiabatic(this%gamma,prim(:,:,:,DEN), &
         prim(:,:,:,PRE))
    
    ! minimal and maximal wave speeds at cell interfaces
    CALL SetWaveSpeeds_euler(prim(:,:,1:2,V_X),this%csound(:,:,1:2), &
         this%tmin(:,:,1:2),this%tmax(:,:,1:2))

    ! x-direction (west and east states)
    this%amin(Mesh%IMIN-1:Mesh%IMAX,:) = MIN(this%tmin(Mesh%IMIN:Mesh%IMAX+1,:,1), &
         this%tmin(Mesh%IMIN-1:Mesh%IMAX,:,2))
    this%amax(Mesh%IMIN-1:Mesh%IMAX,:) = MAX(this%tmax(Mesh%IMIN:Mesh%IMAX+1,:,1), &
         this%tmax(Mesh%IMIN-1:Mesh%IMAX,:,2))

    ! y-direction (south and north states)
    CALL SetWaveSpeeds_euler(prim(:,:,3:4,V_Y),this%csound(:,:,3:4), &
         this%tmin(:,:,1:2),this%tmax(:,:,1:2))

    this%bmin(:,Mesh%JMIN-1:Mesh%JMAX) = MIN(this%tmin(:,Mesh%JMIN:Mesh%JMAX+1,1), &
         this%tmin(:,Mesh%JMIN-1:Mesh%JMAX,2))
    this%bmax(:,Mesh%JMIN-1:Mesh%JMAX) = MAX(this%tmax(:,Mesh%JMIN:Mesh%JMAX+1,1), &
         this%tmax(:,Mesh%JMIN-1:Mesh%JMAX,2))
  END SUBROUTINE CalculateWaveSpeeds_faces


  PURE SUBROUTINE CalculateFluxesX_euler3Drotsym(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: nmin,nmax
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%vnum) &
         :: prim,cons,xfluxes
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,nmin,nmax,prim,cons
    INTENT(OUT)       :: xfluxes
    !------------------------------------------------------------------------!
    CALL CalculateFlux_euler3Drotsym(prim(:,:,nmin:nmax,DEN),prim(:,:,nmin:nmax,V_X), &
         prim(:,:,nmin:nmax,PRE),cons(:,:,nmin:nmax,M_X),cons(:,:,nmin:nmax,M_Y), &
         cons(:,:,nmin:nmax,M_Z),cons(:,:,nmin:nmax,ENE),xfluxes(:,:,nmin:nmax,1), &
         xfluxes(:,:,nmin:nmax,2),xfluxes(:,:,nmin:nmax,3), &
         xfluxes(:,:,nmin:nmax,4),xfluxes(:,:,nmin:nmax,5))
   END SUBROUTINE CalculateFluxesX_euler3Drotsym


  PURE SUBROUTINE CalculateFluxesY_euler3Drotsym(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: nmin,nmax
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%vnum) &
         :: prim,cons,yfluxes
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,nmin,nmax,prim,cons
    INTENT(OUT)       :: yfluxes
    !------------------------------------------------------------------------!
    CALL CalculateFlux_euler3Drotsym(prim(:,:,nmin:nmax,DEN),prim(:,:,nmin:nmax,V_Y), &
         prim(:,:,nmin:nmax,PRE),cons(:,:,nmin:nmax,M_Y),cons(:,:,nmin:nmax,M_X), &
         cons(:,:,nmin:nmax,M_Z),cons(:,:,nmin:nmax,ENE),yfluxes(:,:,nmin:nmax,1), &
         yfluxes(:,:,nmin:nmax,3), yfluxes(:,:,nmin:nmax,2), &
         yfluxes(:,:,nmin:nmax,4),yfluxes(:,:,nmin:nmax,5))
  END SUBROUTINE CalculateFluxesY_euler3Drotsym


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
    CALL CentrifugalForces(pvar(:,:,DEN),pvar(:,:,V_Z),Mesh%bhz(:,:), &
         Mesh%czxz(:,:,1),Mesh%czyz(:,:,1),this%fcent(:,:,1,1),this%fcent(:,:,2,1))

    ! geometrical source terms in momentum equationes
    ! with centrifugal forces
    sterm(:,:,M_X) = MomentumSourcesX_euler2D(cvar(:,:,M_Y),pvar(:,:,V_X), &
         pvar(:,:,V_Y),pvar(:,:,PRE),Mesh%cxyx(:,:,1), &
         Mesh%cyxy(:,:,1),Mesh%czxz(:,:,1)) + this%fcent(:,:,1,1)
    sterm(:,:,M_Y) = MomentumSourcesY_euler2D(cvar(:,:,M_X),pvar(:,:,V_X), &
         pvar(:,:,V_Y),pvar(:,:,PRE),Mesh%cxyx(:,:,1), &
         Mesh%cyxy(:,:,1),Mesh%czyz(:,:,1)) + this%fcent(:,:,2,1)
    sterm(:,:,M_Z) = -pvar(:,:,V_Z)*(Mesh%czxz(:,:,1)*cvar(:,:,M_X) + &
         Mesh%czyz(:,:,1)*cvar(:,:,M_Y))
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
    CALL CentrifugalForces(prim(:,:,:,DEN),prim(:,:,:,V_Z),Mesh%chz(:,:,:), &
         Mesh%czxz(:,:,:),Mesh%czyz(:,:,:),this%fcent(:,:,:,1),this%fcent(:,:,:,2))

    ! geometrical source terms in momentum equationes
    ! sum up all four corner values
    sterm(:,:,M_X) = SUM(MomentumSourcesX_euler2D(cons(:,:,:,M_Y), &
         prim(:,:,:,V_X),prim(:,:,:,V_Y),prim(:,:,:,PRE),Mesh%cxyx(:,:,:), &
         Mesh%cyxy(:,:,:),Mesh%czxz(:,:,:)) + this%fcent(:,:,:,1),DIM=3)
    sterm(:,:,M_Y) = SUM(MomentumSourcesY_euler2D(cons(:,:,:,M_X), &
         prim(:,:,:,V_X),prim(:,:,:,V_Y),prim(:,:,:,PRE),Mesh%cxyx(:,:,:), &
         Mesh%cyxy(:,:,:),Mesh%czyz(:,:,:)) + this%fcent(:,:,:,2),DIM=3)
    sterm(:,:,M_Z) = SUM(-prim(:,:,:,V_Z)*(Mesh%czxz(:,:,:)*cons(:,:,:,M_X) + &
         Mesh%czyz(:,:,:)*cons(:,:,:,M_Y)),DIM=3)
  END SUBROUTINE GeometricalSources_faces


  PURE SUBROUTINE ExternalSources_euler3Drotsym(this,Mesh,accel,pvar,cvar,sterm)
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
         :: accel
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
         :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,accel,pvar,cvar
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    sterm(:,:,DEN) = 0.
    sterm(:,:,M_X) = pvar(:,:,DEN) * accel(:,:,1)
    sterm(:,:,M_Y) = pvar(:,:,DEN) * accel(:,:,2)
    sterm(:,:,M_Z) = 0.
    sterm(:,:,ENE) = cvar(:,:,M_X) * accel(:,:,1) + cvar(:,:,M_Y) * accel(:,:,2)
  END SUBROUTINE ExternalSources_euler3Drotsym


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
    CALL Cons2Prim_euler3Drotsym(this%gamma,cvar(:,:,DEN),cvar(:,:,M_X), &
         cvar(:,:,M_Y),cvar(:,:,M_Z),cvar(:,:,ENE),pvar(:,:,DEN),pvar(:,:,V_X), &
         pvar(:,:,V_Y),pvar(:,:,V_Z),pvar(:,:,PRE))
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
    CALL Cons2Prim_euler3Drotsym(this%gamma,cons(:,:,:,DEN),cons(:,:,:,M_X), &
         cons(:,:,:,M_Y),cons(:,:,:,M_Z),cons(:,:,:,ENE),prim(:,:,:,DEN), &
         prim(:,:,:,V_X),prim(:,:,:,V_Y),prim(:,:,:,V_Z),prim(:,:,:,PRE))
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
    CALL Prim2Cons_euler3Drotsym(this%gamma,pvar(:,:,DEN),pvar(:,:,V_X), &
         pvar(:,:,V_Y),pvar(:,:,V_Z),pvar(:,:,PRE),cvar(:,:,DEN),cvar(:,:,M_X), &
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
    CALL Prim2Cons_euler3Drotsym(this%gamma,prim(:,:,:,DEN),prim(:,:,:,V_X), &
         prim(:,:,:,V_Y),prim(:,:,:,V_Z),prim(:,:,:,PRE),cons(:,:,:,DEN), &
         cons(:,:,:,M_X),cons(:,:,:,M_Y),cons(:,:,:,M_Z),cons(:,:,:,ENE))
  END SUBROUTINE Convert2Conservative_faces


  PURE SUBROUTINE GetAxisMasks_euler3Drotsym(this,reflX,reflY)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    LOGICAL, DIMENSION(this%vnum) :: reflX,reflY
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this
    INTENT(OUT)       :: reflX,reflY
    !------------------------------------------------------------------------!
    reflX(:) = .FALSE.
    reflY(:) = .FALSE.
    reflX(V_X) = .TRUE.
    reflX(V_Z) = .TRUE.
    reflY(V_Y) = .TRUE.
    reflY(V_Z) = .TRUE.
  END SUBROUTINE GetAxisMasks_euler3Drotsym


  ELEMENTAL SUBROUTINE CalculateFlux_euler3Drotsym(rho,v,P,m1,m2,m3,E,f1,f2,f3,f4,f5)
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
  END SUBROUTINE CalculateFlux_euler3Drotsym


  ELEMENTAL SUBROUTINE CentrifugalForces(rho,vz,r,czxz,czyz,fcx,fcy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho,vz,r,czxz,czyz
    REAL, INTENT(OUT) :: fcx,fcy
    !------------------------------------------------------------------------!
    REAL :: tmp2
    !------------------------------------------------------------------------!
    tmp2 = rho*vz*vz
    fcx  = tmp2 * czxz
    fcy  = tmp2 * czyz
  END SUBROUTINE CentrifugalForces


  ELEMENTAL SUBROUTINE Cons2Prim_euler3Drotsym(gamma,rho_in,mu,mv,mw,E,rho_out,u,v,w,P)
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
  END SUBROUTINE Cons2Prim_euler3Drotsym

  
  ELEMENTAL SUBROUTINE Prim2Cons_euler3Drotsym(gamma,rho_in,u,v,w,P,rho_out,mu,mv,mw,E)
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
  END SUBROUTINE Prim2Cons_euler3Drotsym

  SUBROUTINE ClosePhysics_euler3Drotsym(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL ClosePhysics_euler2D(this)
    DEALLOCATE(this%fcent)
  END SUBROUTINE ClosePhysics_euler3Drotsym


END MODULE physics_euler3Drotsym
