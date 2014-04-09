!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_generic.f90                                               #
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
! generic module for the advection problem
!----------------------------------------------------------------------------!
MODULE physics_generic
  USE physics_euler2D
  USE physics_euler3Drotsym
  USE physics_euler3Drotamt
  USE constants_generic
  USE sources_common, ONLY : Sources_TYP
  USE mesh_common, ONLY : Mesh_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE GeometricalSources
     MODULE PROCEDURE GeometricalSources_center
     MODULE PROCEDURE GeometricalSources_faces
  END INTERFACE
  INTERFACE Convert2Primitive
     MODULE PROCEDURE Convert2Primitive_center
     MODULE PROCEDURE Convert2Primitive_faces
  END INTERFACE
  INTERFACE Convert2Conservative
     MODULE PROCEDURE Convert2Conservative_center
     MODULE PROCEDURE Convert2Conservative_faces
  END INTERFACE
  !--------------------------------------------------------------------------!
  ! flags for advection problems
  INTEGER, PARAMETER :: EULER2D          = 1
  INTEGER, PARAMETER :: EULER3D_ROTSYM   = 2
  INTEGER, PARAMETER :: EULER3D_ROTAMT   = 3
  !--------------------------------------------------------------------------!
  ! basic numerical constants
  REAL, PARAMETER :: PI = 3.141592653589792
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Physics_TYP, &
       ! constants
       PI, &
       EULER2D, EULER3D_ROTSYM, EULER3D_ROTAMT, &
       SI, CGS, GEOMETRICAL, &
       ! methods
       InitPhysics, &
       MallocPhysics, &
       CheckData, &
       CalculateWaveSpeeds, &
       MaxWaveSpeeds, &
       CalculateFluxesX, &
       CalculateFluxesY, &
       GeometricalSources, &
       ExternalSources, &
       ViscositySources, &
       Convert2Primitive, &
       Convert2Conservative, &
       ReflectionMasks, &
       AxisMasks, &
       GetType, &
       GetName, &
       GetRank, &
       Info, &
       Warning, &
       Error, &
       ClosePhysics
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitPhysics(this,problem,units,gamma,mu,rhomin,pmin,dpmax)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    INTEGER           :: problem
    INTEGER, OPTIONAL :: units
    REAL, OPTIONAL    :: gamma,mu,rhomin,pmin,dpmax
    !------------------------------------------------------------------------!
    INTENT(IN)        :: problem,units,gamma,mu,rhomin,pmin,dpmax
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! units
    IF (PRESENT(units)) THEN
       CALL InitConstants(this%constants,units)
    ELSE
       CALL InitConstants(this%constants,SI)
    END IF

    ! ratio of specific heats
    IF (PRESENT(gamma)) THEN
       this%gamma = gamma
    ELSE
       this%gamma = 1.4
    END IF

    ! mean molecular weight
    IF (PRESENT(mu)) THEN
       this%mu = mu
    ELSE
       this%mu = 0.029 ! air
    END IF

    ! density minimum, i.e. vacuum
    IF (PRESENT(rhomin)) THEN
       this%rhomin = rhomin
    ELSE
       this%rhomin = 1.0E-30
    END IF

    ! pressure minimum
    IF (PRESENT(pmin)) THEN
       this%pmin = pmin
    ELSE
       this%pmin = 1.0E-30
    END IF

    ! maximal pressure gradient
    IF (PRESENT(dpmax)) THEN
       this%dpmax = dpmax
    ELSE
       this%dpmax = 1.0
    END IF

    SELECT CASE(problem)
    CASE(EULER2D)
       CALL InitPhysics_euler2D(this,problem)
    CASE(EULER3D_ROTSYM)
       CALL InitPhysics_euler3Drs(this,problem)
    CASE(EULER3D_ROTAMT)
       CALL InitPhysics_euler3Dra(this,problem)
    CASE DEFAULT
       CALL Error(this, "InitPhysics", "Unknown advection problem.")
    END SELECT

    NULLIFY(this%sources)

    ! print some information
    CALL Info(this, " PHYSICS--> advection problem: " // TRIM(GetName(this)))
  END SUBROUTINE InitPhysics


  SUBROUTINE MallocPhysics(this,Mesh)
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
    ! allocate memory for arrays common to all physics
    ALLOCATE(this%amin(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%amax(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%bmin(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%bmax(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%tmin(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
         this%tmax(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
         STAT = err)
    IF (err.NE.0) THEN
       CALL Error(this, "MallocPhysics", "Unable to allocate memory.")
    END IF
    ! call specific allocation procedures
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL MallocPhysics_euler2D(this,Mesh)
    CASE(EULER3D_ROTSYM)
       CALL MallocPhysics_euler3Drs(this,Mesh)
    CASE(EULER3D_ROTAMT)
       CALL MallocPhysics_euler3Dra(this,Mesh)
    END SELECT
  END SUBROUTINE MallocPhysics


  PURE SUBROUTINE AxisMasks(this,reflX,reflY)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    LOGICAL, DIMENSION(this%vnum) :: reflX,reflY
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this
    INTENT(OUT)       :: reflX,reflY
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL AxisMasks_euler2D(this,reflX,reflY)
    CASE(EULER3D_ROTSYM)
       CALL AxisMasks_euler3Drs(this,reflX,reflY)
    CASE(EULER3D_ROTAMT)
       CALL AxisMasks_euler3Dra(this,reflX,reflY)
    END SELECT
  END SUBROUTINE AxisMasks


  PURE FUNCTION CheckData(this,Mesh,pvar,pold) RESULT(bad_data)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
         :: pvar,pold
    INTEGER           :: bad_data
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,pvar,pold
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       bad_data = CheckData_euler2D(this,Mesh,pvar,pold)
    CASE(EULER3D_ROTSYM)
       bad_data = CheckData_euler3Drs(this,Mesh,pvar,pold)
    CASE(EULER3D_ROTAMT)
       bad_data = CheckData_euler3Dra(this,Mesh,pvar,pold)
    END SELECT
  END FUNCTION CheckData


  PURE SUBROUTINE CalculateWaveSpeeds(this,Mesh,prim)
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
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalculateWaveSpeeds_euler2D(this,Mesh,prim)
    CASE(EULER3D_ROTSYM)
       CALL CalculateWaveSpeeds_euler3Drs(this,Mesh,prim)
    CASE(EULER3D_ROTAMT)
       CALL CalculateWaveSpeeds_euler3Dra(this,Mesh,prim)
    END SELECT
  END SUBROUTINE CalculateWaveSpeeds


  PURE SUBROUTINE MaxWaveSpeeds(this,Mesh,pvar,amax)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
         :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
         :: amax
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: amax
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalculateWaveSpeeds_euler2D(this,Mesh,pvar)
    CASE(EULER3D_ROTSYM)
       CALL CalculateWaveSpeeds_euler3Drs(this,Mesh,pvar)
    CASE(EULER3D_ROTAMT)
       CALL CalculateWaveSpeeds_euler3Dra(this,Mesh,pvar)
    END SELECT

    amax = MAX(this%tmax,-this%tmin)
  END SUBROUTINE MaxWaveSpeeds


  PURE SUBROUTINE CalculateFluxesX(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: nmin,nmax
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%vnum) :: &
         prim,cons,xfluxes
    !------------------------------------------------------------------------!
    INTENT(IN)       :: this,Mesh,nmin,nmax,prim,cons
    INTENT(OUT)      :: xfluxes
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalculateFluxesX_euler2D(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    CASE(EULER3D_ROTSYM)
       CALL CalculateFluxesX_euler3Drs(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    CASE(EULER3D_ROTAMT)
       CALL CalculateFluxesX_euler3Dra(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    END SELECT
  END SUBROUTINE CalculateFluxesX


  PURE SUBROUTINE CalculateFluxesY(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)      :: Mesh
    INTEGER             :: nmin,nmax
    REAL, DIMENSION(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1,4,this%vnum) :: &
         prim,cons,yfluxes
    !------------------------------------------------------------------------!
    INTENT(IN)       :: this,Mesh,nmin,nmax,prim,cons
    INTENT(OUT)      :: yfluxes
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalculateFluxesY_euler2D(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    CASE(EULER3D_ROTSYM)
       CALL CalculateFluxesY_euler3Drs(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    CASE(EULER3D_ROTAMT)
       CALL CalculateFluxesY_euler3Dra(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    END SELECT
  END SUBROUTINE CalculateFluxesY


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
    ! calculate geometrical sources depending on the advection problem
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL GeometricalSources_euler2D(this,Mesh,pvar,cvar,sterm)
    CASE(EULER3D_ROTSYM)
       CALL GeometricalSources_euler3Drs(this,Mesh,pvar,cvar,sterm)
    CASE(EULER3D_ROTAMT)
       CALL GeometricalSources_euler3Dra(this,Mesh,pvar,cvar,sterm)
    END SELECT
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
    ! calculate geometrical sources depending on the advection problem
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL GeometricalSources_euler2D(this,Mesh,prim,cons,sterm)
    CASE(EULER3D_ROTSYM)
       CALL GeometricalSources_euler3Drs(this,Mesh,prim,cons,sterm)
    CASE(EULER3D_ROTAMT)
       CALL GeometricalSources_euler3Dra(this,Mesh,prim,cons,sterm)
    END SELECT
  END SUBROUTINE GeometricalSources_faces


  PURE SUBROUTINE ExternalSources(this,Mesh,accel,pvar,cvar,sterm)
    IMPLICIT NONE
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
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL ExternalSources_euler2D(this,Mesh,accel,pvar,cvar,sterm)
    CASE(EULER3D_ROTSYM)
       CALL ExternalSources_euler3Drs(this,Mesh,accel,pvar,cvar,sterm)
    CASE(EULER3D_ROTAMT)
       CALL ExternalSources_euler3Dra(this,Mesh,accel,pvar,cvar,sterm)
    END SELECT
  END SUBROUTINE ExternalSources


  PURE SUBROUTINE ViscositySources(this,Mesh,Sources,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Sources_TYP) :: Sources
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
         :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,pvar,cvar
    INTENT(INOUT)     :: Sources
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL ViscositySources_euler2D(this,Mesh,Sources,pvar,cvar,sterm)
    CASE(EULER3D_ROTSYM)
       CALL ViscositySources_euler3Drs(this,Mesh,Sources,pvar,cvar,sterm)
!    CASE(EULER3D_ROTAMT)
!       CALL ExternalSources_euler3Dra(this,Mesh,accel,pvar,cvar,sterm)
    END SELECT
  END SUBROUTINE ViscositySources


  
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
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL Convert2Primitive_euler2D(this,Mesh,cvar,pvar)
    CASE(EULER3D_ROTSYM)
       CALL Convert2Primitive_euler3Drs(this,Mesh,cvar,pvar)
    CASE(EULER3D_ROTAMT)
       CALL Convert2Primitive_euler3Dra(this,Mesh,cvar,pvar)
    END SELECT
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
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL Convert2Primitive_euler2D(this,Mesh,cons,prim)
    CASE(EULER3D_ROTSYM)
       CALL Convert2Primitive_euler3Drs(this,Mesh,cons,prim)
    CASE(EULER3D_ROTAMT)
       CALL Convert2Primitive_euler3Dra(this,Mesh,cons,prim)
    END SELECT
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
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL Convert2Conservative_euler2D(this,Mesh,pvar,cvar)
    CASE(EULER3D_ROTSYM)
       CALL Convert2Conservative_euler3Drs(this,Mesh,pvar,cvar)
    CASE(EULER3D_ROTAMT)
       CALL Convert2Conservative_euler3Dra(this,Mesh,pvar,cvar)
    END SELECT
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
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL Convert2Conservative_euler2D(this,Mesh,prim,cons)
    CASE(EULER3D_ROTSYM)
       CALL Convert2Conservative_euler3Drs(this,Mesh,prim,cons)
    CASE(EULER3D_ROTAMT)
       CALL Convert2Conservative_euler3Dra(this,Mesh,prim,cons)
    END SELECT
  END SUBROUTINE Convert2Conservative_faces
  

  PURE SUBROUTINE ReflectionMasks(this,reflX,reflY)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    LOGICAL, DIMENSION(this%vnum) :: reflX,reflY
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this
    INTENT(OUT)       :: reflX,reflY
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL ReflectionMasks_euler2D(this,reflX,reflY)
    CASE(EULER3D_ROTSYM)
       CALL ReflectionMasks_euler3Drs(this,reflX,reflY)
    CASE(EULER3D_ROTAMT)
       CALL ReflectionMasks_euler3Dra(this,reflX,reflY)
    END SELECT
  END SUBROUTINE ReflectionMasks


  SUBROUTINE ClosePhysics(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! call specific dallocation procedures
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL ClosePhysics_euler2D(this)
    CASE(EULER3D_ROTSYM)
       CALL ClosePhysics_euler3Drs(this)
    CASE(EULER3D_ROTAMT)
       CALL ClosePhysics_euler3Dra(this)
    END SELECT

    ! deallocate memory for all arrays used in physics module
    DEALLOCATE(this%amin,this%amax,this%bmin,this%bmax, &
         this%tmin,this%tmax)
  END SUBROUTINE ClosePhysics

END MODULE physics_generic
