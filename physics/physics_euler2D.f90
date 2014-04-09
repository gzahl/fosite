!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler2D.f90                                               #
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
! basic module for 2D Euler equations
!----------------------------------------------------------------------------!
MODULE physics_euler2D
  USE physics_common
  USE mesh_common, ONLY : Mesh_TYP
  USE output_common, ONLY : Datastruct_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTERFACE CalculateWaveSpeeds_euler2D
     MODULE PROCEDURE CalculateWaveSpeeds_center
     MODULE PROCEDURE CalculateWaveSpeeds_faces
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
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: DEN = 1
  INTEGER, PARAMETER :: V_X = 2
  INTEGER, PARAMETER :: V_Y = 3
  INTEGER, PARAMETER :: PRE = 4
  INTEGER, PARAMETER :: M_X = 2
  INTEGER, PARAMETER :: M_Y = 3
  INTEGER, PARAMETER :: ENE = 4
  INTEGER, PARAMETER :: num_var = 4
  REAL, PARAMETER :: TINY = 1.0E-30              ! to avoid division by 0    !
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 2D"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitPhysics_euler2D, &
       MallocPhysics_euler2D, &
       CheckData_euler2D, &
       CalculateWaveSpeeds_euler2D, &
       CalculateFluxesX_euler2D, &
       CalculateFluxesY_euler2D, &
       GeometricalSources_euler2D, &
       ExternalSources_euler2D, &
       Convert2Primitive_euler2D, &
       Convert2Conservative_euler2D, &
       GetDataStruct_euler2D, &
       GetReflectionMasks_euler2D, &
       GetAxisMasks_euler2D, &
       GetSoundSpeed_adiabatic, &
       SetWaveSpeeds_euler, &
       MomentumSourcesX_euler2D, &
       MomentumSourcesY_euler2D, &
       ClosePhysics_euler2D
  !--------------------------------------------------------------------------!

CONTAINS

  PURE SUBROUTINE InitPhysics_euler2D(this,problem,gamma,dpmax)
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
    this%dpmax = 1.0
  END SUBROUTINE InitPhysics_euler2D


  SUBROUTINE MallocPhysics_euler2D(this,Mesh)
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
    ! allocate memory for arrays used in Euler2D
    ALLOCATE(this%csound(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4), &
         STAT = err)
    IF (err.NE.0) THEN
       PRINT *, "ERROR in MallocPhysics_euler2D: Can't allocate memory!"
       STOP
    END IF    
  END SUBROUTINE MallocPhysics_euler2D

  
  SUBROUTINE GetDataStruct_euler2D(this,ds)
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
       PRINT *, "ERROR in GetDataStruct_euler2D: Unable to allocate memory "
       STOP
    END IF

    ds(1)%name  = 'density'
    ds(1)%rank  = 0
    ds(1)%shape = 1
    ds(2)%name = 'velocity'
    ds(2)%rank  = 1
    ds(2)%shape = 2
    ds(3)%name  = 'pressure'
    ds(3)%rank  = 0
    ds(3)%shape = 1
  END SUBROUTINE GetDataStruct_euler2D


  PURE FUNCTION CheckData_euler2D(this,Mesh,pvar,pold) RESULT (bad_data)
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
  END FUNCTION CheckData_euler2D


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


  PURE SUBROUTINE CalculateFluxesX_euler2D(this,Mesh,nmin,nmax,prim,cons,xfluxes)
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
    CALL CalculateFlux_euler2D(prim(:,:,nmin:nmax,DEN),prim(:,:,nmin:nmax,V_X), &
         prim(:,:,nmin:nmax,PRE),cons(:,:,nmin:nmax,M_X),cons(:,:,nmin:nmax,M_Y), &
         cons(:,:,nmin:nmax,ENE),xfluxes(:,:,nmin:nmax,1),xfluxes(:,:,nmin:nmax,2), &
         xfluxes(:,:,nmin:nmax,3),xfluxes(:,:,nmin:nmax,4))
   END SUBROUTINE CalculateFluxesX_euler2D


  PURE SUBROUTINE CalculateFluxesY_euler2D(this,Mesh,nmin,nmax,prim,cons,yfluxes)
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
    CALL CalculateFlux_euler2D(prim(:,:,nmin:nmax,DEN),prim(:,:,nmin:nmax,V_Y), &
         prim(:,:,nmin:nmax,PRE),cons(:,:,nmin:nmax,M_Y),cons(:,:,nmin:nmax,M_X), &
         cons(:,:,nmin:nmax,ENE),yfluxes(:,:,nmin:nmax,1),yfluxes(:,:,nmin:nmax,3), &
         yfluxes(:,:,nmin:nmax,2),yfluxes(:,:,nmin:nmax,4))
  END SUBROUTINE CalculateFluxesY_euler2D


  PURE SUBROUTINE GeometricalSources_center(this,Mesh,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
         :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,pvar,cvar
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! geometrical source terms in momentum equationes
    sterm(:,:,DEN) = 0.
    sterm(:,:,M_X) = MomentumSourcesX_euler2D(cvar(:,:,M_Y),pvar(:,:,V_X), &
         pvar(:,:,V_Y),pvar(:,:,PRE),Mesh%cxyx(:,:,1),Mesh%cyxy(:,:,1),Mesh%czxz(:,:,1))
    sterm(:,:,M_Y) = MomentumSourcesY_euler2D(cvar(:,:,M_X),pvar(:,:,V_X), &
         pvar(:,:,V_Y),pvar(:,:,PRE),Mesh%cxyx(:,:,1),Mesh%cyxy(:,:,1),Mesh%czyz(:,:,1))
    sterm(:,:,ENE) = 0.
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
    INTENT(IN)        :: this,Mesh,prim,cons
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! geometrical source terms in momentum equationes
    ! sum up all four corner values
    sterm(:,:,DEN) = 0.
    sterm(:,:,M_X) = SUM(MomentumSourcesX_euler2D(cons(:,:,:,M_Y), &
         prim(:,:,:,V_X),prim(:,:,:,V_Y),prim(:,:,:,PRE),Mesh%cxyx(:,:,:), &
         Mesh%cyxy(:,:,:),Mesh%czxz(:,:,:)),DIM=3)
    sterm(:,:,M_Y) = SUM(MomentumSourcesY_euler2D(cons(:,:,:,M_X), &
         prim(:,:,:,V_X),prim(:,:,:,V_Y),prim(:,:,:,PRE),Mesh%cxyx(:,:,:), &
         Mesh%cyxy(:,:,:),Mesh%czyz(:,:,:)),DIM=3)
    sterm(:,:,ENE) = 0.
  END SUBROUTINE GeometricalSources_faces


  ! momentum and energy sources due to external force
  PURE SUBROUTINE ExternalSources_euler2D(this,Mesh,accel,pvar,cvar,sterm)
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
    sterm(:,:,ENE) = cvar(:,:,M_X) * accel(:,:,1) + cvar(:,:,M_Y) * accel(:,:,2)
  END SUBROUTINE ExternalSources_euler2D


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
    CALL Cons2Prim_euler2D(this%gamma,cvar(:,:,DEN),cvar(:,:,M_X),cvar(:,:,M_Y), &
            cvar(:,:,ENE),pvar(:,:,DEN),pvar(:,:,V_X),pvar(:,:,V_Y),pvar(:,:,PRE))
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
    CALL Cons2Prim_euler2D(this%gamma,cons(:,:,:,DEN),cons(:,:,:,M_X), &
         cons(:,:,:,M_Y),cons(:,:,:,ENE),prim(:,:,:,DEN),prim(:,:,:,V_X), &
         prim(:,:,:,V_Y),prim(:,:,:,PRE))
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
    CALL Prim2Cons_euler2D(this%gamma,pvar(:,:,DEN),pvar(:,:,V_X),pvar(:,:,V_Y), &
            pvar(:,:,PRE),cvar(:,:,DEN),cvar(:,:,M_X),cvar(:,:,M_Y),cvar(:,:,ENE))
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
    CALL Prim2Cons_euler2D(this%gamma,prim(:,:,:,DEN),prim(:,:,:,V_X), &
         prim(:,:,:,V_Y),prim(:,:,:,PRE),cons(:,:,:,DEN),cons(:,:,:,M_X), &
         cons(:,:,:,M_Y),cons(:,:,:,ENE))
  END SUBROUTINE Convert2Conservative_faces


  PURE SUBROUTINE GetReflectionMasks_euler2D(this,reflX,reflY)
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
    reflY(V_Y) = .TRUE.
  END SUBROUTINE GetReflectionMasks_euler2D


  PURE SUBROUTINE GetAxisMasks_euler2D(this,reflX,reflY)
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
    reflX(V_Y) = .TRUE.
    reflY(V_Y) = .TRUE.
    reflY(V_X) = .TRUE.
  END SUBROUTINE GetAxisMasks_euler2D


  ELEMENTAL FUNCTION GetSoundSpeed_adiabatic(gamma,density,pressure) RESULT(cs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,density,pressure
    REAL :: cs
    !------------------------------------------------------------------------!
    cs = SQRT(MAX(TINY,gamma*pressure/density))
  END FUNCTION GetSoundSpeed_adiabatic

  
  ELEMENTAL SUBROUTINE SetWaveSpeeds_euler(velocity,cs,amin,amax)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: velocity,cs
    REAL, INTENT(OUT) :: amin,amax
    !------------------------------------------------------------------------!
    ! minimal and maximal wave speeds
    amin = MIN(0.,velocity-cs)
    amax = MAX(0.,velocity+cs)
  END SUBROUTINE SetWaveSpeeds_euler

  
  ELEMENTAL SUBROUTINE CalculateFlux_euler2D(rho,v,P,m1,m2,E,f1,f2,f3,f4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho,v,P,m1,m2,E
    REAL, INTENT(OUT) :: f1, f2, f3, f4
    !------------------------------------------------------------------------!
    f1 = rho*v
    f2 = m1*v + P
    f3 = m2*v
    f4 = (E+P)*v
  END SUBROUTINE CalculateFlux_euler2D


  ! momentum source terms due to inertial forces
  ELEMENTAL FUNCTION MomentumSourcesX_euler2D(my,vx,vy,P,cxyx,cyxy,czxz) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: my,vx,vy,P,cxyx,cyxy,czxz
    REAL :: st
    !------------------------------------------------------------------------!
    st = -my * (cxyx * vx - cyxy * vy) + (cyxy + czxz) * P
  END FUNCTION MomentumSourcesX_euler2D


  ELEMENTAL FUNCTION MomentumSourcesY_euler2D(mx,vx,vy,P,cxyx,cyxy,czyz) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: mx,vx,vy,P,cxyx,cyxy,czyz
    REAL :: st
    !------------------------------------------------------------------------!
    st = mx * (cxyx * vx - cyxy * vy) + (cxyx + czyz) * P
  END FUNCTION MomentumSourcesY_euler2D


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
    DEALLOCATE(this%csound)
  END SUBROUTINE ClosePhysics_euler2D


END MODULE physics_euler2D
