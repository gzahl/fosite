!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_generic.f90                                               #
!#                                                                           #
!# Copyright (C) 2007 - 2014                                                 #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
!# Manuel Jung      <mjung@astrophysik.uni-kiel.de>                          #
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
!! - general physics settings
!! \key{problem,INTEGER,advection problem
!!      (see \link physics_generic \endlink for a list of currently supported
!!       advection problems)}
!! \key{units,INTEGER,unit system
!!      (see \link constants_generic \endlink for a list of currently supported
!!       unit systems)}
!! \key{mu,REAL,mean molecular weight (default is for air
!!      at normal conditions),0.029}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Björn Sperling
!! \author Manuel Jung
!!
!! \brief generic module for the advection problem
!!
!! \ingroup physics
!----------------------------------------------------------------------------!
MODULE physics_generic
  USE physics_euler2Disothm, ClosePhysics_common => ClosePhysics
  USE physics_euler2D
  USE physics_euler2Disoiamt
  USE physics_euler2Diamt
  USE physics_euler2Diamrot
  USE physics_euler2Disoiamrot
  USE physics_euler2Dsgs
  USE physics_euler3Drotsym
  USE physics_euler3Drotamt
  USE physics_euler3Drotsymsgs
  USE physics_euler3Drotamtsgs
  USE constants_generic
  USE sources_common, ONLY : Sources_TYP
  USE mesh_common, ONLY : Mesh_TYP, Initialized
  USE mesh_generic, ONLY : remap_bounds, CARTESIAN
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE GeometricalSources
     MODULE PROCEDURE GeometricalSources_center
     MODULE PROCEDURE GeometricalSources_faces
  END INTERFACE
  INTERFACE Convert2Primitive
     MODULE PROCEDURE Convert2Primitive_center
     MODULE PROCEDURE Convert2Primitive_centsub
     MODULE PROCEDURE Convert2Primitive_faces
     MODULE PROCEDURE Convert2Primitive_facesub
  END INTERFACE
  INTERFACE Convert2Conservative
     MODULE PROCEDURE Convert2Conservative_center
     MODULE PROCEDURE Convert2Conservative_centsub
     MODULE PROCEDURE Convert2Conservative_faces
     MODULE PROCEDURE Convert2Conservative_facesub
  END INTERFACE
  INTERFACE GetSoundSpeed_adiabatic
     MODULE PROCEDURE GetSoundSpeed_euler2D
  END INTERFACE
  INTERFACE SetSoundSpeeds
     MODULE PROCEDURE SetSoundSpeeds_center
     MODULE PROCEDURE SetSoundSpeeds_faces
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  ! flags for advection problems
  INTEGER, PARAMETER :: EULER2D             = 1
  INTEGER, PARAMETER :: EULER2D_ISOTHERM    = 2
  INTEGER, PARAMETER :: EULER3D_ROTSYM      = 3
  INTEGER, PARAMETER :: EULER3D_ROTAMT      = 4
  INTEGER, PARAMETER :: EULER3D_ROTSYMSGS   = 5
  INTEGER, PARAMETER :: EULER2D_SGS         = 7
  INTEGER, PARAMETER :: EULER3D_ROTAMTSGS   = 8
  INTEGER, PARAMETER :: EULER2D_ISOIAMT     = 9
  INTEGER, PARAMETER :: EULER2D_IAMT        = 11
  INTEGER, PARAMETER :: EULER2D_IAMROT      = 12
  INTEGER, PARAMETER :: EULER2D_ISOIAMROT   = 13
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Physics_TYP, &
       ! constants
       EULER2D, EULER2D_ISOTHERM, EULER3D_ROTSYM, EULER3D_ROTAMT, & 
       EULER2D_SGS, EULER3D_ROTSYMSGS, EULER3D_ROTAMTSGS, &
       SI, CGS, GEOMETRICAL, EULER2D_ISOIAMT, &
       EULER2D_IAMT, EULER2D_IAMROT, EULER2D_ISOIAMROT, &
       ! methods
       InitPhysics, &
       CalculateWaveSpeeds, &
       MaxWaveSpeeds, &
       CalculateFluxesX, &
       CalculateFluxesY, &
       CalculateCharSystemX, &
       CalculateCharSystemY, &
       CalculateBoundaryDataX, &
       CalculateBoundaryDataY, &
       CalculatePrim2RiemannX, &
       CalculatePrim2RiemannY, &
       CalculateRiemann2PrimX, &
       CalculateRiemann2PrimY, &
       CalculateStresses, &
       SetSoundSpeeds, &
       GetSoundSpeeds, &
       UpdateSoundSpeed,&
       GeometricalSources, &
       ExternalSources, &
       ViscositySources, &
       SGSSources, &
       CalculateSGSTensor, &
       Convert2Primitive, &
       Convert2Conservative, &
       ReflectionMasks, &
       AxisMasks, &
       GetSoundSpeed_adiabatic, &
       ClosePhysics, &
       GetType, &
       GetName, &
       GetRank, &
       GetNumProcs, &
       Initialized, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitPhysics(this,Mesh,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Dict_TYP),POINTER &
                      :: config, IO
    !------------------------------------------------------------------------!
    INTEGER           :: problem
    INTEGER           :: units
    INTEGER           :: err, i, valwrite
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this, Mesh
    !------------------------------------------------------------------------!
    ! check initialization of Mesh
    IF (.NOT.Initialized(Mesh)) &
         CALL Error(this,"InitPhysics","mesh module uninitialized")

    ! units
    IF (.NOT.Initialized(this%constants)) THEN
       CALL RequireKey(config, "units", SI)
       CALL GetAttr(config, "units", units)
       CALL InitConstants(this%constants,units)
    END IF

    CALL RequireKey(config, "problem")

    ! ratio of specific heats
    CALL RequireKey(config, "gamma", 1.4)

    ! mean molecular weight
    CALL RequireKey(config, "mu", 0.029)    ! air

    ! isothermal sound speed
    CALL RequireKey(config, "cs", 0.)

    ! speed of rotating frame
    CALL RequireKey(config, "omega", 0.0)

    ! center of rotation (cartesian coordinates)
    CALL RequireKey(config, "centrot_x", 0.0)
    CALL RequireKey(config, "centrot_y", 0.0)

    ! softening parameter to smooth out singularity near center of rotation
    ! (only necessary, if it's inside the computational domain)
    ! set to 0.0 to disable
    ! the softening length is the product of this parameter and the
    ! size of the grid cell next to the center of rotation; thus a value larger
    ! than 1.0 leads to larger softening whereas smaller values will
    ! probably cause odd behaviour due to the 1/r singularity;
    ! if the minimal r on the computational domain is larger than
    ! the size of the associated grid cell, softening is disabled, because
    ! the center of rotation lies outside of the computational domain
    CALL RequireKey(config, "softening", 1.0)
 
    CALL GetAttr(config, "problem", problem)
    CALL GetAttr(config, "gamma", this%gamma)
    CALL GetAttr(config, "mu", this%mu)
    CALL GetAttr(config, "cs", this%csiso)
    CALL GetAttr(config, "omega", this%Omega)
    CALL GetAttr(config, "centrot_x", this%centrot(1))
    CALL GetAttr(config, "centrot_y", this%centrot(2))    
    CALL GetAttr(config, "softening", this%eps)

   ! allocate memory for arrays common to all physics modules
    ALLOCATE(this%amin(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%amax(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%bmin(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%bmax(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%tmp(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%tmin(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4), &
         this%tmax(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4), &
         this%bccsound(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%fcsound(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4), &
         STAT = err)
    IF (err.NE.0) &
         CALL Error(this, "InitPhysics", "Unable to allocate memory.")

    this%amax(:,:) = 0.

    IF(this%csiso.GT.0.) THEN
      this%bccsound(:,:) = this%csiso
      this%fcsound(:,:,:) = this%csiso
    END IF

    SELECT CASE(problem)
    CASE(EULER2D)
       CALL InitPhysics_euler2D(this,Mesh,problem)
    CASE(EULER2D_ISOTHERM)
       CALL InitPhysics_euler2Dit(this,Mesh,problem)
    CASE(EULER3D_ROTSYM)
       CALL InitPhysics_euler3Drs(this,Mesh,problem)
    CASE(EULER3D_ROTAMT)
       CALL InitPhysics_euler3Dra(this,Mesh,problem)
    CASE(EULER3D_ROTSYMSGS)
       CALL InitPhysics_euler3DrsSGS(this,Mesh,problem)
    CASE(EULER3D_ROTAMTSGS)
       CALL InitPhysics_euler3DraSGS(this,Mesh,problem)
    CASE(EULER2D_SGS)
       CALL InitPhysics_euler2Dsgs(this,Mesh,problem)
    CASE(EULER2D_IAMT)
       CALL InitPhysics_euler2Dia(this,Mesh,problem)
    CASE(EULER2D_IAMROT)
       CALL InitPhysics_euler2Diar(this,Mesh,problem)
    CASE(EULER2D_ISOIAMT)
       CALL InitPhysics_euler2Ditia(this,Mesh,problem)
    CASE(EULER2D_ISOIAMROT)
       CALL InitPhysics_euler2Ditiar(this,Mesh,problem)
    CASE DEFAULT
       CALL Error(this, "InitPhysics", "Unknown advection problem.")
    END SELECT

    ! enable/disable absorbing and farfield boundary conditions
    SELECT CASE(problem)
    CASE(EULER2D,EULER2D_ISOTHERM,&
         EULER2D_ISOIAMT,EULER2D_IAMT,EULER2D_SGS, &
         EULER3D_ROTSYM,EULER3D_ROTSYMSGS,EULER3D_ROTAMT,EULER3D_ROTAMTSGS)
       this%supports_absorbing = .TRUE.       
    CASE DEFAULT
       this%supports_absorbing = .FALSE.
    END SELECT
    SELECT CASE(problem)
    CASE(EULER2D,EULER2D_SGS, &
         EULER3D_ROTSYM,EULER3D_ROTSYMSGS,EULER3D_ROTAMT,EULER3D_ROTAMTSGS)
       this%supports_farfield = .TRUE.       
    CASE DEFAULT
       this%supports_farfield  = .FALSE.
    END SELECT

    ! reset source term pointer
    NULLIFY(this%sources)

    ! check if output of sound speeds is requested
    valwrite = 0
    IF (HasKey(config, "output/bccsound")) CALL GetAttr(config, "output/bccsound", valwrite)
    IF (valwrite .EQ. 1) THEN
       CALL AddField(IO,&
                     "bccsound",&
                     remap_bounds(Mesh, this%bccsound),&
                     Dict("name" / "bccsound"))
    END IF
    valwrite = 0
    IF (HasKey(config, "output/fcsound")) CALL GetAttr(config, "output/fcsound", valwrite)
    IF (valwrite .EQ. 1) THEN
       CALL AddField(IO,&
                     "fcsound",&
                     remap_bounds(Mesh, this%fcsound),&
                     Dict("name" / "fcsound"))
    END IF

     this%time = -1.
    ! print some information
    CALL Info(this, " PHYSICS--> advection problem: " // TRIM(GetName(this)))

  END SUBROUTINE InitPhysics


  PURE SUBROUTINE SetSoundSpeeds_center(this,Mesh,bccsound)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: bccsound
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,bccsound
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    this%bccsound(:,:) = bccsound(:,:)
  END SUBROUTINE SetSoundSpeeds_center


  PURE SUBROUTINE SetSoundSpeeds_faces(this,Mesh,fcsound)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4) &
                      :: fcsound
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,fcsound
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    this%fcsound(:,:,:) = fcsound(:,:,:)
  END SUBROUTINE SetSoundSpeeds_faces


  FUNCTION GetSoundSpeeds(this) RESULT(bccsound)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    REAL,DIMENSION(:,:),POINTER &
                      :: bccsound
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this
    !------------------------------------------------------------------------!
    bccsound => this%bccsound
  END FUNCTION GetSoundSpeeds


  PURE SUBROUTINE AxisMasks(this,reflX,reflY)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    LOGICAL, DIMENSION(this%vnum) :: reflX,reflY
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this
    INTENT(OUT)       :: reflX,reflY
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL AxisMasks_euler2D(this,reflX,reflY)
    CASE(EULER2D_ISOTHERM)
       CALL AxisMasks_euler2Dit(this,reflX,reflY)
    CASE(EULER3D_ROTSYM)
       CALL AxisMasks_euler3Drs(this,reflX,reflY)
    CASE(EULER3D_ROTAMT)
       CALL AxisMasks_euler3Dra(this,reflX,reflY)
    CASE(EULER3D_ROTSYMSGS)
       CALL AxisMasks_euler3DrsSGS(this,reflX,reflY)
    CASE(EULER3D_ROTAMTSGS)
       CALL AxisMasks_euler3DraSGS(this,reflX,reflY)
    CASE(EULER2D_SGS)
       CALL AxisMasks_euler2Dsgs(this,reflX,reflY)
    CASE(EULER2D_IAMT)
       CALL AxisMasks_euler2Dia(this,reflX,reflY)
    CASE(EULER2D_IAMROT)
       CALL AxisMasks_euler2Diar(this,reflX,reflY)
    CASE(EULER2D_ISOIAMT)
       CALL AxisMasks_euler2Ditia(this,reflX,reflY)
    CASE(EULER2D_ISOIAMROT)
       CALL AxisMasks_euler2Ditiar(this,reflX,reflY)
    END SELECT
  END SUBROUTINE AxisMasks

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
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalcSoundSpeeds_euler2D(this,Mesh,prim)
       CALL CalcWaveSpeeds_euler2D(this,Mesh,prim)
    CASE(EULER2D_ISOTHERM)
       CALL CalcWaveSpeeds_euler2Dit(this,Mesh,prim)
    CASE(EULER3D_ROTSYM)
       CALL CalcSoundSpeeds_euler3Drs(this,Mesh,prim)
       CALL CalcWaveSpeeds_euler3Drs(this,Mesh,prim)
    CASE(EULER3D_ROTAMT)
       CALL CalcSoundSpeeds_euler3Dra(this,Mesh,prim)
       CALL CalcWaveSpeeds_euler3Dra(this,Mesh,prim)
    CASE(EULER3D_ROTSYMSGS)
       CALL CalcSoundSpeeds_euler3DrsSGS(this,Mesh,prim)
       CALL CalcWaveSpeeds_euler3DrsSGS(this,Mesh,prim)
    CASE(EULER3D_ROTAMTSGS)
       CALL CalcSoundSpeeds_euler3DraSGS(this,Mesh,prim)
       CALL CalcWaveSpeeds_euler3DraSGS(this,Mesh,prim)
    CASE(EULER2D_SGS)
       CALL CalcSoundSpeeds_euler2Dsgs(this,Mesh,prim)
       CALL CalcWaveSpeeds_euler2Dsgs(this,Mesh,prim)
    CASE(EULER2D_IAMT)
       CALL CalcSoundSpeeds_euler2Dia(this,Mesh,prim)
       CALL CalcWaveSpeeds_euler2Dia(this,Mesh,prim)
    CASE(EULER2D_IAMROT)
       CALL CalcSoundSpeeds_euler2Diar(this,Mesh,prim)
       CALL CalcWaveSpeeds_euler2Diar(this,Mesh,prim)
    CASE(EULER2D_ISOIAMT)
       CALL CalcWaveSpeeds_euler2Ditia(this,Mesh,prim)
    CASE(EULER2D_ISOIAMROT)
       CALL CalcWaveSpeeds_euler2Ditiar(this,Mesh,prim)
    END SELECT
  END SUBROUTINE CalculateWaveSpeeds


  PURE SUBROUTINE MaxWaveSpeeds(this,Mesh,time,pvar,amax)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
         :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
         :: amax
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,time
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: amax
    !------------------------------------------------------------------------!
    CALL UpdateSoundSpeed(this,Mesh,time,pvar)
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalcWaveSpeeds_euler2D(this,Mesh,pvar)
    CASE(EULER2D_ISOTHERM)
       CALL CalcWaveSpeeds_euler2Dit(this,Mesh,pvar)
    CASE(EULER3D_ROTSYM)
       CALL CalcWaveSpeeds_euler3Drs(this,Mesh,pvar)
    CASE(EULER3D_ROTAMT)
       CALL CalcWaveSpeeds_euler3Dra(this,Mesh,pvar)
    CASE(EULER3D_ROTSYMSGS)
       CALL CalcWaveSpeeds_euler3DrsSGS(this,Mesh,pvar)
    CASE(EULER3D_ROTAMTSGS)
       CALL CalcWaveSpeeds_euler3DraSGS(this,Mesh,pvar)
    CASE(EULER2D_SGS)
       CALL CalcWaveSpeeds_euler2Dsgs(this,Mesh,pvar)
    CASE(EULER2D_IAMT)
       CALL CalcWaveSpeeds_euler2Dia(this,Mesh,pvar)
    CASE(EULER2D_IAMROT)
       CALL CalcWaveSpeeds_euler2Diar(this,Mesh,pvar)
    CASE(EULER2D_ISOIAMT)
       CALL CalcWaveSpeeds_euler2Ditia(this,Mesh,pvar)
    CASE(EULER2D_ISOIAMROT)
       CALL CalcWaveSpeeds_euler2Ditiar(this,Mesh,pvar)
    END SELECT

    amax(:,:,1) = MAX(this%amax,-this%amin)
    amax(:,:,2) = MAX(this%bmax,-this%bmin)
  END SUBROUTINE MaxWaveSpeeds

  PURE SUBROUTINE UpdateSoundSpeed(this,Mesh,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
         :: pvar
    REAL              :: time
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,time
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!

   IF ((time.NE.this%time) .OR. (time.EQ.0.0)) THEN
      this%time = time
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D_ISOTHERM,EULER2D_ISOIAMT,EULER2D_ISOIAMROT)
       !do nothing 
    CASE(EULER2D)
       CALL CalcSoundSpeeds_euler2D(this,Mesh,pvar)
    CASE(EULER3D_ROTSYM)
       CALL CalcSoundSpeeds_euler3Drs(this,Mesh,pvar)
    CASE(EULER3D_ROTAMT)
       CALL CalcSoundSpeeds_euler3Dra(this,Mesh,pvar)
    CASE(EULER3D_ROTSYMSGS)
       CALL CalcSoundSpeeds_euler3DrsSGS(this,Mesh,pvar)
    CASE(EULER3D_ROTAMTSGS)
       CALL CalcSoundSpeeds_euler3Dra(this,Mesh,pvar)
    CASE(EULER2D_SGS)
       CALL CalcSoundSpeeds_euler2Dsgs(this,Mesh,pvar)
    CASE(EULER2D_IAMT)
       CALL CalcSoundSpeeds_euler2Dia(this,Mesh,pvar)
    CASE(EULER2D_IAMROT)
       CALL CalcSoundSpeeds_euler2Diar(this,Mesh,pvar)
    END SELECT

   END IF 
  END SUBROUTINE UpdateSoundSpeed

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
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalcFluxesX_euler2D(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    CASE(EULER2D_ISOTHERM)
       CALL CalcFluxesX_euler2Dit(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    CASE(EULER3D_ROTSYM)
       CALL CalcFluxesX_euler3Drs(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    CASE(EULER3D_ROTAMT)
       CALL CalcFluxesX_euler3Dra(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    CASE(EULER3D_ROTSYMSGS)
       CALL CalcFluxesX_euler3DrsSGS(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    CASE(EULER3D_ROTAMTSGS)
       CALL CalcFluxesX_euler3DraSGS(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    CASE(EULER2D_SGS)
       CALL CalcFluxesX_euler2Dsgs(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    CASE(EULER2D_IAMT)
       CALL CalcFluxesX_euler2Dia(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    CASE(EULER2D_IAMROT)
       CALL CalcFluxesX_euler2Diar(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    CASE(EULER2D_ISOIAMT)
       CALL CalcFluxesX_euler2Ditia(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    CASE(EULER2D_ISOIAMROT)
       CALL CalculateFluxesX_euler2Ditiar(this,Mesh,nmin,nmax,prim,cons,xfluxes)
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
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalcFluxesY_euler2D(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    CASE(EULER2D_ISOTHERM)
       CALL CalcFluxesY_euler2Dit(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    CASE(EULER3D_ROTSYM)
       CALL CalcFluxesY_euler3Drs(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    CASE(EULER3D_ROTAMT)
       CALL CalcFluxesY_euler3Dra(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    CASE(EULER3D_ROTSYMSGS)
       CALL CalcFluxesY_euler3DrsSGS(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    CASE(EULER3D_ROTAMTSGS)
       CALL CalcFluxesY_euler3DraSGS(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    CASE(EULER2D_SGS)
       CALL CalcFluxesY_euler2Dsgs(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    CASE(EULER2D_IAMT)
       CALL CalcFluxesY_euler2Dia(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    CASE(EULER2D_IAMROT)
       CALL CalcFluxesY_euler2Diar(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    CASE(EULER2D_ISOIAMT)
       CALL CalcFluxesY_euler2Ditia(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    CASE(EULER2D_ISOIAMROT)
       CALL CalculateFluxesY_euler2Ditiar(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    END SELECT
  END SUBROUTINE CalculateFluxesY


  PURE SUBROUTINE CalculateCharSystemX(this,Mesh,i,dir,pvar,lambda,xvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: i,dir
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: lambda,xvar
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,i,dir,pvar
    INTENT(INOUT)     :: lambda
    INTENT(OUT)       :: xvar
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalcCharSystemX_euler2D(this,Mesh,i,dir,pvar,lambda,xvar)
    CASE(EULER2D_ISOTHERM)
       CALL CalcCharSystemX_euler2Dit(this,Mesh,i,dir,pvar,lambda,xvar)
    CASE(EULER3D_ROTSYM)
       CALL CalcCharSystemX_euler3Drs(this,Mesh,i,dir,pvar,lambda,xvar)
    CASE(EULER3D_ROTAMT)
       CALL CalcCharSystemX_euler3Dra(this,Mesh,i,dir,pvar,lambda,xvar)
    CASE(EULER3D_ROTSYMSGS)
       CALL CalcCharSystemX_euler3DrsSGS(this,Mesh,i,dir,pvar,lambda,xvar)
    CASE(EULER3D_ROTAMTSGS)
       CALL CalcCharSystemX_euler3DraSGS(this,Mesh,i,dir,pvar,lambda,xvar)
    CASE(EULER2D_SGS)
       CALL CalcCharSystemX_euler2Dsgs(this,Mesh,i,dir,pvar,lambda,xvar)
    CASE(EULER2D_IAMT)
       CALL CalcCharSystemX_euler2Dia(this,Mesh,i,dir,pvar,lambda,xvar)
    CASE(EULER2D_IAMROT)
!!$ FIXME, not implemented
!!$       CALL CalcCharSystemX_euler2Diar(this,Mesh,i,dir,pvar,lambda,xvar)
    CASE(EULER2D_ISOIAMT)
       CALL CalcCharSystemX_euler2Ditia(this,Mesh,i,dir,pvar,lambda,xvar)
    CASE(EULER2D_ISOIAMROT)
!!$ FIXME, not implemented
!       CALL CalcCharSystemX_euler2Ditiar(this,Mesh,i,dir,pvar,lambda,xvar)
    END SELECT
  END SUBROUTINE CalculateCharSystemX


  PURE SUBROUTINE CalculateCharSystemY(this,Mesh,j,dir,pvar,lambda,xvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: j,dir
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,this%VNUM) :: lambda,xvar
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,j,dir,pvar
    INTENT(INOUT)     :: lambda
    INTENT(OUT)       :: xvar
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalcCharSystemY_euler2D(this,Mesh,j,dir,pvar,lambda,xvar)
    CASE(EULER2D_ISOTHERM)
       CALL CalcCharSystemY_euler2Dit(this,Mesh,j,dir,pvar,lambda,xvar)
    CASE(EULER3D_ROTSYM)
       CALL CalcCharSystemY_euler3Drs(this,Mesh,j,dir,pvar,lambda,xvar)
    CASE(EULER3D_ROTAMT)
       CALL CalcCharSystemY_euler3Dra(this,Mesh,j,dir,pvar,lambda,xvar)
    CASE(EULER3D_ROTSYMSGS)
       CALL CalcCharSystemY_euler3DrsSGS(this,Mesh,j,dir,pvar,lambda,xvar)
    CASE(EULER3D_ROTAMTSGS)
       CALL CalcCharSystemY_euler3DraSGS(this,Mesh,j,dir,pvar,lambda,xvar)
    CASE(EULER2D_SGS)
       CALL CalcCharSystemY_euler2Dsgs(this,Mesh,j,dir,pvar,lambda,xvar)
    CASE(EULER2D_IAMT)
    CASE(EULER2D_IAMROT)
!!$ FIXME, not implemented
!!$       CALL CalcCharSystemY_euler2Diar(this,Mesh,j,dir,pvar,lambda,xvar)
    CASE(EULER2D_ISOIAMT)
    CASE(EULER2D_ISOIAMROT)
    END SELECT
  END SUBROUTINE CalculateCharSystemY


  PURE SUBROUTINE CalculateBoundaryDataX(this,Mesh,i,dir,xvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: i,dir
    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: xvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,i,dir,xvar
    INTENT(INOUT)     :: pvar
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalcBoundaryDataX_euler2D(this,Mesh,i,dir,xvar,pvar)
    CASE(EULER2D_ISOTHERM)
       CALL CalcBoundaryDataX_euler2Dit(this,Mesh,i,dir,xvar,pvar)
    CASE(EULER3D_ROTSYM)
       CALL CalcBoundaryDataX_euler3Drs(this,Mesh,i,dir,xvar,pvar)
    CASE(EULER3D_ROTAMT)
       CALL CalcBoundaryDataX_euler3Dra(this,Mesh,i,dir,xvar,pvar)
    CASE(EULER3D_ROTSYMSGS)
       CALL CalcBoundaryDataX_euler3DrsSGS(this,Mesh,i,dir,xvar,pvar)
    CASE(EULER3D_ROTAMTSGS)
       CALL CalcBoundaryDataX_euler3DraSGS(this,Mesh,i,dir,xvar,pvar)
    CASE(EULER2D_SGS)
       CALL CalcBoundaryDataX_euler2Dsgs(this,Mesh,i,dir,xvar,pvar)
    CASE(EULER2D_IAMT)
       CALL CalcBoundaryDataX_euler2Dia(this,Mesh,i,dir,xvar,pvar)
    CASE(EULER2D_IAMROT)
!!$ FIXME, not implemented
!!$       CALL CalcBoundaryDataX_euler2Diar(this,Mesh,i,dir,xvar,pvar)
    CASE(EULER2D_ISOIAMT)
       CALL CalcBoundaryDataX_euler2Ditia(this,Mesh,i,dir,xvar,pvar)
    CASE(EULER2D_ISOIAMROT)
!!$ FIXME, not implemented
!       CALL CalcBoundaryDataX_euler2Ditiar(this,Mesh,i,dir,xvar,pvar)
    END SELECT
  END SUBROUTINE CalculateBoundaryDataX


  PURE SUBROUTINE CalculateBoundaryDataY(this,Mesh,j,dir,xvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: j,dir
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,this%VNUM) :: xvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,j,dir,xvar
    INTENT(INOUT)     :: pvar
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalcBoundaryDataY_euler2D(this,Mesh,j,dir,xvar,pvar)
    CASE(EULER2D_ISOTHERM)
       CALL CalcBoundaryDataY_euler2Dit(this,Mesh,j,dir,xvar,pvar)
    CASE(EULER3D_ROTSYM)
       CALL CalcBoundaryDataY_euler3Drs(this,Mesh,j,dir,xvar,pvar)
    CASE(EULER3D_ROTAMT)
       CALL CalcBoundaryDataY_euler3Dra(this,Mesh,j,dir,xvar,pvar)
    CASE(EULER3D_ROTSYMSGS)
       CALL CalcBoundaryDataY_euler3Drssgs(this,Mesh,j,dir,xvar,pvar)
    CASE(EULER3D_ROTAMTSGS)
       CALL CalcBoundaryDataY_euler3Drasgs(this,Mesh,j,dir,xvar,pvar)
    CASE(EULER2D_SGS)
       CALL CalcBoundaryDataY_euler2Dsgs(this,Mesh,j,dir,xvar,pvar)
    CASE(EULER2D_IAMT)
    CASE(EULER2D_IAMROT)
!!$ FIXME, not implemented
!!$       CALL CalcBoundaryDataY_euler2Diar(this,Mesh,j,dir,xvar,pvar)
    CASE(EULER2D_ISOIAMT)
    CASE(EULER2D_ISOIAMROT)
    END SELECT
  END SUBROUTINE CalculateBoundaryDataY

  PURE SUBROUTINE CalculatePrim2RiemannX(this,Mesh,i,pvar,lambda,Rinv)
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
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalcPrim2RiemannX_euler2D(this,Mesh,i,pvar,lambda,Rinv)
    CASE(EULER2D_IAMROT)
!!$ FIXME, not implemented
!!$       CALL CalcPrim2RiemannX_euler2Diar(this,Mesh,i,pvar,lambda,Rinv)
    CASE(EULER2D_SGS)
       CALL CalcPrim2RiemannX_euler2Dsgs(this,Mesh,i,pvar,lambda,Rinv)
    CASE(EULER3D_ROTSYM)
       CALL CalcPrim2RiemannX_euler3Drs(this,Mesh,i,pvar,lambda,Rinv)
    CASE(EULER3D_ROTSYMSGS)
       CALL CalcPrim2RiemannX_euler3Drssgs(this,Mesh,i,pvar,lambda,Rinv)
    CASE(EULER3D_ROTAMT)
       CALL CalcPrim2RiemannX_euler3Dra(this,Mesh,i,pvar,lambda,Rinv)
    CASE(EULER3D_ROTAMTSGS)
       CALL CalcPrim2RiemannX_euler3Drasgs(this,Mesh,i,pvar,lambda,Rinv)
    END SELECT
  END SUBROUTINE CalculatePrim2RiemannX

  PURE SUBROUTINE CalculatePrim2RiemannY(this,Mesh,j,pvar,lambda,Rinv)
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
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalcPrim2RiemannY_euler2D(this,Mesh,j,pvar,lambda,Rinv)
    CASE(EULER2D_IAMROT)
!!$ FIXME, not implemented
!!$       CALL CalcPrim2RiemannY_euler2Diar(this,Mesh,j,pvar,lambda,Rinv)
    CASE(EULER2D_SGS)
       CALL CalcPrim2RiemannY_euler2Dsgs(this,Mesh,j,pvar,lambda,Rinv)
    CASE(EULER3D_ROTSYM)
       CALL CalcPrim2RiemannY_euler3Drs(this,Mesh,j,pvar,lambda,Rinv)
    CASE(EULER3D_ROTSYMSGS)
       CALL CalcPrim2RiemannY_euler3Drssgs(this,Mesh,j,pvar,lambda,Rinv)
    CASE(EULER3D_ROTAMT)
       CALL CalcPrim2RiemannY_euler3Dra(this,Mesh,j,pvar,lambda,Rinv) 
    CASE(EULER3D_ROTAMTSGS)
       CALL CalcPrim2RiemannY_euler3Drasgs(this,Mesh,j,pvar,lambda,Rinv) 
    END SELECT
  END SUBROUTINE CalculatePrim2RiemannY

 PURE SUBROUTINE CalculateRiemann2PrimX(this,Mesh,i,Rinv,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: i
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    REAL, DIMENSION(Mesh%JMIN:Mesh%JMAX,this%VNUM) :: Rinv
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,i,Rinv
    INTENT(INOUT)     :: pvar
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalcRiemann2PrimX_euler2D(this,Mesh,i,Rinv,pvar)
    CASE(EULER2D_IAMROT)
!!$ FIXME, not implemented
!!$       CALL CalcRiemann2PrimX_euler2Diar(this,Mesh,i,Rinv,pvar)
    CASE(EULER2D_SGS)
       CALL CalcRiemann2PrimX_euler2Dsgs(this,Mesh,i,Rinv,pvar)
    CASE(EULER3D_ROTSYM)
       CALL CalcRiemann2PrimX_euler3Drs(this,Mesh,i,Rinv,pvar)
    CASE(EULER3D_ROTSYMSGS)
       CALL CalcRiemann2PrimX_euler3Drssgs(this,Mesh,i,Rinv,pvar)
    CASE(EULER3D_ROTAMT)
       CALL CalcRiemann2PrimX_euler3Dra(this,Mesh,i,Rinv,pvar)
    CASE(EULER3D_ROTAMTSGS)
       CALL CalcRiemann2PrimX_euler3Drasgs(this,Mesh,i,Rinv,pvar)
    END SELECT
  END SUBROUTINE CalculateRiemann2PrimX

  PURE SUBROUTINE CalculateRiemann2PrimY(this,Mesh,j,Rinv,pvar)
    IMPLICIT NONE
   !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: j
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    REAL, DIMENSION(Mesh%IMIN:Mesh%IMAX,this%VNUM) :: Rinv
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,j,Rinv
    INTENT(INOUT)     :: pvar
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalcRiemann2PrimY_euler2D(this,Mesh,j,Rinv,pvar)
    CASE(EULER2D_IAMROT)
!!$ FIXME, not implemented
!!$       CALL CalcRiemann2PrimY_euler2Diar(this,Mesh,j,Rinv,pvar)
    CASE(EULER2D_SGS)
       CALL CalcRiemann2PrimY_euler2Dsgs(this,Mesh,j,Rinv,pvar)
    CASE(EULER3D_ROTSYM)
       CALL CalcRiemann2PrimY_euler3Drs(this,Mesh,j,Rinv,pvar)
    CASE(EULER3D_ROTSYMSGS)
       CALL CalcRiemann2PrimY_euler3Drssgs(this,Mesh,j,Rinv,pvar)
    CASE(EULER3D_ROTAMT)
       CALL CalcRiemann2PrimY_euler3Dra(this,Mesh,j,Rinv,pvar)
    CASE(EULER3D_ROTAMTSGS)
       CALL CalcRiemann2PrimY_euler3Drasgs(this,Mesh,j,Rinv,pvar)
    END SELECT
  END SUBROUTINE CalculateRiemann2PrimY

  PURE SUBROUTINE CalculateStresses(this,Mesh,pvar,dynvis,bulkvis, &
       btxx,btxy,btxz,btyy,btyz,btzz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
         dynvis,bulkvis,btxx,btxy,btxz,btyy,btyz,btzz
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,dynvis,bulkvis
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: btxx,btxy,btxz,btyy,btyz,btzz
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalcStresses_euler2D(this,Mesh,pvar,dynvis,bulkvis, &
            btxx,btxy,btyy)
    CASE(EULER2D_ISOTHERM)
       CALL CalcStresses_euler2Dit(this,Mesh,pvar,dynvis,bulkvis, &
            btxx,btxy,btyy)
    CASE(EULER3D_ROTSYM)
       CALL CalcStresses_euler3Drs(this,Mesh,pvar,dynvis,bulkvis, &
            btxx,btxy,btxz,btyy,btyz,btzz)
    CASE(EULER3D_ROTAMT)
       CALL CalcStresses_euler3Dra(this,Mesh,pvar,dynvis,bulkvis, &
            btxx,btxy,btxz,btyy,btyz,btzz)
    CASE(EULER3D_ROTSYMSGS)
       CALL CalcStresses_euler3Drssgs(this,Mesh,pvar,dynvis,bulkvis, &
            btxx,btxy,btxz,btyy,btyz,btzz)
    CASE(EULER3D_ROTAMTSGS)
       CALL CalcStresses_euler3Drasgs(this,Mesh,pvar,dynvis,bulkvis, &
            btxx,btxy,btxz,btyy,btyz,btzz)
    CASE(EULER2D_SGS)
       CALL CalcStresses_euler2Dsgs(this,Mesh,pvar,dynvis,bulkvis, &
            btxx,btxy,btyy)
    CASE(EULER2D_IAMT)
       CALL CalcStresses_euler2Dia(this,Mesh,pvar,dynvis,bulkvis, &
            btxx,btxy,btyy)
    CASE(EULER2D_IAMROT)
       CALL CalcStresses_euler2Diar(this,Mesh,pvar,dynvis,bulkvis, &
            btxx,btxy,btyy)
    CASE(EULER2D_ISOIAMT)
       CALL CalcStresses_euler2Ditia(this,Mesh,pvar,dynvis,bulkvis, &
            btxx,btxy,btyy)
    CASE(EULER2D_ISOIAMROT)
       CALL CalcStresses_euler2Ditiar(this,Mesh,pvar,dynvis,bulkvis, &
            btxx,btxy,btyy)
    END SELECT
  END SUBROUTINE CalculateStresses


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
    ! compute geometrical source only for non-cartesian mesh except for the
    ! EULER2D_IAMROT case for which geometrical sources are always necessary.
    IF ((GetType(Mesh%geometry).NE.CARTESIAN).OR. &
        (GetType(this).EQ.EULER2D_IAMROT).OR. &
        (GetType(this).EQ.EULER2D_ISOIAMROT)) THEN
    ! calculate geometrical sources depending on the advection problem
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
! ATTENTION: Don't use inline expansion here. It will yield false results,
!            because GeometricalSources_XXX is a generic interface for
!            GeometricalSources_faces and GeometricalSources_center in each
!            submodule. However, if one uses inline expansion it allways
!            points on the specific routines in the euler2D module.
!CDIR NOIEXPAND
       CALL GeometricalSources_euler2D(this,Mesh,pvar,cvar,sterm)
    CASE(EULER2D_ISOTHERM)
!CDIR NOIEXPAND
       CALL GeometricalSources_euler2Dit(this,Mesh,pvar,cvar,sterm)
    CASE(EULER3D_ROTSYM)
!CDIR NOIEXPAND
       CALL GeometricalSources_euler3Drs(this,Mesh,pvar,cvar,sterm)
    CASE(EULER3D_ROTAMT)
!CDIR NOIEXPAND
       CALL GeometricalSources_euler3Dra(this,Mesh,pvar,cvar,sterm)
    CASE(EULER3D_ROTSYMSGS)
!CDIR NOIEXPAND
       CALL GeometricalSources_euler3DrsSGS(this,Mesh,pvar,cvar,sterm)
   CASE(EULER3D_ROTAMTSGS)
!CDIR NOIEXPAND
       CALL GeometricalSources_euler3DraSGS(this,Mesh,pvar,cvar,sterm)
    CASE(EULER2D_SGS)
!CDIR NOIEXPAND
       CALL GeometricalSources_euler2Dsgs(this,Mesh,pvar,cvar,sterm)
    CASE(EULER2D_IAMT)
!CDIR NOIEXPAND
       CALL GeometricalSources_euler2Dia(this,Mesh,pvar,cvar,sterm)
    CASE(EULER2D_IAMROT)
!CDIR NOIEXPAND
       CALL GeometricalSources_euler2Diar(this,Mesh,pvar,cvar,sterm)
    CASE(EULER2D_ISOIAMT)
!CDIR NOIEXPAND
       CALL GeometricalSources_euler2Ditia(this,Mesh,pvar,cvar,sterm)
    CASE(EULER2D_ISOIAMROT)
!CDIR NOIEXPAND
       CALL GeometricalSrcs_euler2Ditiar(this,Mesh,pvar,cvar,sterm)
    END SELECT
    ! reset ghost cell data
    sterm(Mesh%IGMIN:Mesh%IMIN-1,:,:) = 0.0
    sterm(Mesh%IMAX+1:Mesh%IGMAX,:,:) = 0.0
    sterm(:,Mesh%JGMIN:Mesh%JMIN-1,:) = 0.0
    sterm(:,Mesh%JMAX+1:Mesh%JGMAX,:) = 0.0
    END IF
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
    ! compute geometrical source only for non-cartesian mesh except for the
    ! EULER2D_IAMROT case for which geometrical sources are always necessary.
    IF ((GetType(Mesh%geometry).NE.CARTESIAN).OR. &
        (GetType(this).EQ.EULER2D_IAMROT)) THEN
    ! calculate geometrical sources depending on the advection problem
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
! ATTENTION: Don't use inline expansion here. It will yield false results,
!            because GeometricalSources_XXX is a generic interface for
!            GeometricalSources_faces and GeometricalSources_center in each
!            submodule. However, if one uses inline expansion it allways
!            points on the specific routines in the euler2D module.
!CDIR NOIEXPAND
       CALL GeometricalSources_euler2D(this,Mesh,prim,cons,sterm)
    CASE(EULER2D_ISOTHERM)
!CDIR NOIEXPAND
       CALL GeometricalSources_euler2Dit(this,Mesh,prim,cons,sterm)
    CASE(EULER3D_ROTSYM)
!CDIR NOIEXPAND
       CALL GeometricalSources_euler3Drs(this,Mesh,prim,cons,sterm)
    CASE(EULER3D_ROTAMT)
!CDIR NOIEXPAND
       CALL GeometricalSources_euler3Dra(this,Mesh,prim,cons,sterm)
    CASE(EULER3D_ROTSYMSGS)
!CDIR NOIEXPAND
       CALL GeometricalSources_euler3DrsSGS(this,Mesh,prim,cons,sterm)
    CASE(EULER3D_ROTAMTSGS)
!CDIR NOIEXPAND
       CALL GeometricalSources_euler3DraSGS(this,Mesh,prim,cons,sterm)
    CASE(EULER2D_SGS)
!CDIR NOIEXPAND
       CALL GeometricalSources_euler2Dsgs(this,Mesh,prim,cons,sterm)
    CASE(EULER2D_IAMT)
!CDIR NOIEXPAND
       CALL GeometricalSources_euler2Dia(this,Mesh,prim,cons,sterm)
    CASE(EULER2D_IAMROT)
!CDIR NOIEXPAND
       CALL GeometricalSources_euler2Diar(this,Mesh,prim,cons,sterm)
    CASE(EULER2D_ISOIAMT)
!CDIR NOIEXPAND
       CALL GeometricalSources_euler2Ditia(this,Mesh,prim,cons,sterm)
    CASE(EULER2D_ISOIAMROT)
!CDIR NOIEXPAND
       CALL GeometricalSrcs_euler2Ditiar(this,Mesh,prim,cons,sterm)
    END SELECT
    ! reset ghost cell data
    sterm(Mesh%IGMIN:Mesh%IMIN-1,:,:) = 0.0
    sterm(Mesh%IMAX+1:Mesh%IGMAX,:,:) = 0.0
    sterm(:,Mesh%JGMIN:Mesh%JMIN-1,:) = 0.0
    sterm(:,Mesh%JMAX+1:Mesh%JGMAX,:) = 0.0
    END IF
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
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
!CDIR IEXPAND
       CALL ExternalSources_euler2D(this,Mesh,accel,pvar,cvar,sterm)
    CASE(EULER2D_ISOTHERM)
!CDIR IEXPAND
       CALL ExternalSources_euler2Dit(this,Mesh,accel,pvar,cvar,sterm)
    CASE(EULER3D_ROTSYM)
!CDIR IEXPAND
       CALL ExternalSources_euler3Drs(this,Mesh,accel,pvar,cvar,sterm)
    CASE(EULER3D_ROTAMT)
!CDIR IEXPAND
       CALL ExternalSources_euler3Dra(this,Mesh,accel,pvar,cvar,sterm)
    CASE(EULER3D_ROTSYMSGS)
!CDIR IEXPAND
       CALL ExternalSources_euler3DrsSGS(this,Mesh,accel,pvar,cvar,sterm)
    CASE(EULER3D_ROTAMTSGS)
!CDIR IEXPAND
       CALL ExternalSources_euler3DraSGS(this,Mesh,accel,pvar,cvar,sterm)
    CASE(EULER2D_SGS)
!CDIR IEXPAND
       CALL ExternalSources_euler2Dsgs(this,Mesh,accel,pvar,cvar,sterm)
    CASE(EULER2D_IAMT)
!CDIR IEXPAND
       CALL ExternalSources_euler2Dia(this,Mesh,accel,pvar,cvar,sterm)
    CASE(EULER2D_IAMROT)
!CDIR IEXPAND
       CALL ExternalSources_euler2Diar(this,Mesh,accel,pvar,cvar,sterm)
    CASE(EULER2D_ISOIAMT)
!CDIR IEXPAND
       CALL ExternalSources_euler2Ditia(this,Mesh,accel,pvar,cvar,sterm)
    CASE(EULER2D_ISOIAMROT)
!CDIR IEXPAND
       CALL ExternalSources_euler2Ditiar(this,Mesh,accel,pvar,cvar,sterm)
    END SELECT
  END SUBROUTINE ExternalSources


  PURE SUBROUTINE ViscositySources(this,Mesh,pvar,btxx,btxy,btxz,btyy,btyz,btzz,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: &
         pvar,sterm
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
         btxx,btxy,btxz,btyy,btyz,btzz
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,btxx,btxy,btxz,btyy,btyz,btzz
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
!CDIR IEXPAND
       CALL ViscositySources_euler2D(this,Mesh,pvar,btxx,btxy,btyy,sterm)
    CASE(EULER2D_ISOTHERM)
!CDIR IEXPAND
       CALL ViscositySources_euler2Dit(this,Mesh,pvar,btxx,btxy,btyy,sterm)
    CASE(EULER3D_ROTSYM)
!CDIR IEXPAND
       CALL ViscositySources_euler3Drs(this,Mesh,pvar,btxx,btxy,btxz,btyy, &
            btyz,btzz,sterm)
    CASE(EULER3D_ROTAMT)
!CDIR IEXPAND
       CALL ViscositySources_euler3Dra(this,Mesh,pvar,btxx,btxy,btxz,btyy, &
            btyz,btzz,sterm)
    CASE(EULER3D_ROTSYMSGS)
!CDIR IEXPAND
       CALL ViscositySources_euler3DrsSGS(this,Mesh,pvar,btxx,btxy,btxz,btyy, &
            btyz,btzz,sterm)
    CASE(EULER3D_ROTAMTSGS)
!CDIR IEXPAND
       CALL ViscositySources_euler3DraSGS(this,Mesh,pvar,btxx,btxy,btxz,btyy, &
            btyz,btzz,sterm)
    CASE(EULER2D_SGS)
!CDIR IEXPAND
       CALL ViscositySources_euler2Dsgs(this,Mesh,pvar,btxx,btxy,btyy, &
            sterm)
    CASE(EULER2D_IAMT)
!CDIR IEXPAND
       CALL ViscositySources_euler2Dia(this,Mesh,pvar,btxx,btxy,btyy,sterm)
    CASE(EULER2D_IAMROT)
!CDIR IEXPAND
       CALL ViscositySources_euler2Diar(this,Mesh,pvar,btxx,btxy,btyy,sterm)
    CASE(EULER2D_ISOIAMT)
!CDIR IEXPAND
       CALL ViscositySources_euler2Ditia(this,Mesh,pvar,btxx,btxy,btyy,sterm)
    CASE(EULER2D_ISOIAMROT)
!CDIR IEXPAND
       CALL ViscositySources_euler2Ditiar(this,Mesh,pvar,btxx,btxy,btyy,sterm)
    END SELECT
  END SUBROUTINE ViscositySources

 PURE SUBROUTINE SGSSources(this,Mesh,Sources,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Sources_TYP) :: Sources
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: &
         pvar,cvar,sterm
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
!         btxx,btxy,btxz,btyy,btyz,btzz
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,cvar
    INTENT(INOUT)     :: this,Sources
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER3D_ROTSYMSGS)
!CDIR IEXPAND
       CALL SGSSources_euler3DrsSGS(this,Mesh,Sources,pvar,cvar,sterm)
    CASE(EULER3D_ROTAMTSGS)
!CDIR IEXPAND
       CALL SGSSources_euler3DraSGS(this,Mesh,Sources,pvar,cvar,sterm)
    CASE(EULER2D_SGS)
!CDIR IEXPAND
       CALL SGSSources_euler2Dsgs(this,Mesh,Sources,pvar,cvar,sterm)
    END SELECT
  END SUBROUTINE SGSSources

  PURE SUBROUTINE CalculateSGSTensor(this,Mesh,Sources,C,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Sources_TYP) :: Sources
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: &
         C,pvar
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,C,pvar
    INTENT(INOUT)     :: this,Sources
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER3D_ROTSYMSGS)
!CDIR IEXPAND
       CALL CalcSGSTensor_euler3DrsSGS(this,Mesh,Sources,C,pvar)
    CASE(EULER3D_ROTAMTSGS)
!CDIR IEXPAND
       CALL CalcSGSTensor_euler3DraSGS(this,Mesh,Sources,C,pvar)
    CASE(EULER2D_SGS)
!CDIR IEXPAND
       CALL CalcSGSTensor_euler2Dsgs(this,Mesh,Sources,C,pvar)
    END SELECT
  END SUBROUTINE CalculateSGSTensor

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
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL Convert2Primitive_euler2D(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cvar,pvar)
    CASE(EULER2D_ISOTHERM)
       CALL Convert2Primitive_euler2Dit(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cvar,pvar)
    CASE(EULER3D_ROTSYM)
       CALL Convert2Primitive_euler3Drs(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cvar,pvar)
    CASE(EULER3D_ROTAMT)
       CALL Convert2Primitive_euler3Dra(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cvar,pvar)
    CASE(EULER3D_ROTSYMSGS)
       CALL Convert2Primitive_euler3DrsSGS(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cvar,pvar)
    CASE(EULER3D_ROTAMTSGS)
       CALL Convert2Primitive_euler3DraSGS(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cvar,pvar)
    CASE(EULER2D_SGS)
       CALL Convert2Primitive_euler2Dsgs(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cvar,pvar)
    CASE(EULER2D_IAMT)
       CALL Convert2Primitive_euler2Dia(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cvar,pvar)
    CASE(EULER2D_IAMROT)
       CALL Convert2Primitive_euler2Diar(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cvar,pvar)
    CASE(EULER2D_ISOIAMT)
       CALL Convert2Primitive_euler2Ditia(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cvar,pvar)
    CASE(EULER2D_ISOIAMROT)
       CALL Convert2Primitive_euler2Ditiar(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cvar,pvar)
    END SELECT
  END SUBROUTINE Convert2Primitive_center

  
  PURE SUBROUTINE Convert2Primitive_centsub(this,Mesh,i1,i2,j1,j2,cvar,pvar)
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
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL Convert2Primitive_euler2D(this,Mesh,i1,i2,j1,j2,cvar,pvar)
    CASE(EULER2D_ISOTHERM)
       CALL Convert2Primitive_euler2Dit(this,Mesh,i1,i2,j1,j2,cvar,pvar)
    CASE(EULER3D_ROTSYM)
       CALL Convert2Primitive_euler3Drs(this,Mesh,i1,i2,j1,j2,cvar,pvar)
    CASE(EULER3D_ROTAMT)
       CALL Convert2Primitive_euler3Dra(this,Mesh,i1,i2,j1,j2,cvar,pvar)
    CASE(EULER3D_ROTSYMSGS)
       CALL Convert2Primitive_euler3DrsSGS(this,Mesh,i1,i2,j1,j2,cvar,pvar)
    CASE(EULER3D_ROTAMTSGS)
       CALL Convert2Primitive_euler3DraSGS(this,Mesh,i1,i2,j1,j2,cvar,pvar)
    CASE(EULER2D_SGS)
       CALL Convert2Primitive_euler2Dsgs(this,Mesh,i1,i2,j1,j2,cvar,pvar)
    CASE(EULER2D_IAMT)
       CALL Convert2Primitive_euler2Dia(this,Mesh,i1,i2,j1,j2,cvar,pvar)
    CASE(EULER2D_IAMROT)
       CALL Convert2Primitive_euler2Diar(this,Mesh,i1,i2,j1,j2,cvar,pvar)
    CASE(EULER2D_ISOIAMT)
       CALL Convert2Primitive_euler2Ditia(this,Mesh,i1,i2,j1,j2,cvar,pvar)
    CASE(EULER2D_ISOIAMROT)
       CALL Convert2Primitive_euler2Ditiar(this,Mesh,i1,i2,j1,j2,cvar,pvar)
    END SELECT
  END SUBROUTINE Convert2Primitive_centsub


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
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL Convert2Primitive_euler2D(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cons,prim)
    CASE(EULER2D_ISOTHERM)
       CALL Convert2Primitive_euler2Dit(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cons,prim)
    CASE(EULER3D_ROTSYM)
       CALL Convert2Primitive_euler3Drs(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cons,prim)
    CASE(EULER3D_ROTAMT)
       CALL Convert2Primitive_euler3Dra(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cons,prim)
    CASE(EULER3D_ROTSYMSGS)
       CALL Convert2Primitive_euler3DrsSGS(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cons,prim)
    CASE(EULER3D_ROTAMTSGS)
       CALL Convert2Primitive_euler3DraSGS(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cons,prim)
    CASE(EULER2D_SGS)
       CALL Convert2Primitive_euler2Dsgs(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cons,prim)
    CASE(EULER2D_IAMT)
       CALL Convert2Primitive_euler2Dia(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cons,prim)
    CASE(EULER2D_IAMROT)
       CALL Convert2Primitive_euler2Diar(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cons,prim)
    CASE(EULER2D_ISOIAMT)
       CALL Convert2Primitive_euler2Ditia(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cons,prim)
    CASE(EULER2D_ISOIAMROT)
       CALL Convert2Primitive_euler2Ditiar(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,cons,prim)
    END SELECT
  END SUBROUTINE Convert2Primitive_faces


  PURE SUBROUTINE Convert2Primitive_facesub(this,Mesh,i1,i2,j1,j2,cons,prim)
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
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL Convert2Primitive_euler2D(this,Mesh,i1,i2,j1,j2,cons,prim)
    CASE(EULER2D_ISOTHERM)
       CALL Convert2Primitive_euler2Dit(this,Mesh,i1,i2,j1,j2,cons,prim)
    CASE(EULER3D_ROTSYM)
       CALL Convert2Primitive_euler3Drs(this,Mesh,i1,i2,j1,j2,cons,prim)
    CASE(EULER3D_ROTAMT)
       CALL Convert2Primitive_euler3Dra(this,Mesh,i1,i2,j1,j2,cons,prim)
    CASE(EULER3D_ROTSYMSGS)
       CALL Convert2Primitive_euler3DrsSGS(this,Mesh,i1,i2,j1,j2,cons,prim)
    CASE(EULER3D_ROTAMTSGS)
       CALL Convert2Primitive_euler3DraSGS(this,Mesh,i1,i2,j1,j2,cons,prim)
    CASE(EULER2D_SGS)
       CALL Convert2Primitive_euler2Dsgs(this,Mesh,i1,i2,j1,j2,cons,prim)
    CASE(EULER2D_IAMT)
       CALL Convert2Primitive_euler2Dia(this,Mesh,i1,i2,j1,j2,cons,prim)
    CASE(EULER2D_IAMROT)
       CALL Convert2Primitive_euler2Diar(this,Mesh,i1,i2,j1,j2,cons,prim)
    CASE(EULER2D_ISOIAMT)
       CALL Convert2Primitive_euler2Ditia(this,Mesh,i1,i2,j1,j2,cons,prim)
    CASE(EULER2D_ISOIAMROT)
       CALL Convert2Primitive_euler2Ditiar(this,Mesh,i1,i2,j1,j2,cons,prim)
    END SELECT
  END SUBROUTINE Convert2Primitive_facesub


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
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL Convert2Conservative_euler2D(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,pvar,cvar)
    CASE(EULER2D_ISOTHERM)
       CALL Convert2Conservative_euler2Dit(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,pvar,cvar)
    CASE(EULER3D_ROTSYM)
       CALL Convert2Conservative_euler3Drs(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,pvar,cvar)
    CASE(EULER3D_ROTSYMSGS)
       CALL Convert2Cons_euler3DrsSGS(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,pvar,cvar)
    CASE(EULER3D_ROTAMTSGS)
       CALL Convert2Cons_euler3DraSGS(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,pvar,cvar)
    CASE(EULER3D_ROTAMT)
       CALL Convert2Conservative_euler3Dra(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,pvar,cvar)
    CASE(EULER2D_SGS)
       CALL Convert2Conservative_euler2Dsgs(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,pvar,cvar)
    CASE(EULER2D_IAMT)
       CALL Convert2Conservative_euler2Dia(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,pvar,cvar)
    CASE(EULER2D_IAMROT)
       CALL Convert2Conservative_euler2Diar(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,pvar,cvar)
    CASE(EULER2D_ISOIAMT)
       CALL Convert2Cons_euler2Ditia(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,pvar,cvar)
    CASE(EULER2D_ISOIAMROT)
       CALL Convert2Cons_euler2Ditiar(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,pvar,cvar)
    END SELECT
  END SUBROUTINE Convert2Conservative_center


  PURE SUBROUTINE Convert2Conservative_centsub(this,Mesh,i1,i2,j1,j2,pvar,cvar)
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
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL Convert2Conservative_euler2D(this,Mesh,i1,i2,j1,j2,pvar,cvar)
    CASE(EULER2D_ISOTHERM)
       CALL Convert2Conservative_euler2Dit(this,Mesh,i1,i2,j1,j2,pvar,cvar)
    CASE(EULER3D_ROTSYM)
       CALL Convert2Conservative_euler3Drs(this,Mesh,i1,i2,j1,j2,pvar,cvar)
    CASE(EULER3D_ROTAMT)
       CALL Convert2Conservative_euler3Dra(this,Mesh,i1,i2,j1,j2,pvar,cvar)
    CASE(EULER3D_ROTSYMSGS)
       CALL Convert2Cons_euler3DrsSGS(this,Mesh,i1,i2,j1,j2,pvar,cvar)
    CASE(EULER3D_ROTAMTSGS)
       CALL Convert2Cons_euler3DraSGS(this,Mesh,i1,i2,j1,j2,pvar,cvar)
    CASE(EULER2D_SGS)
       CALL Convert2Conservative_euler2Dsgs(this,Mesh,i1,i2,j1,j2,pvar,cvar)
    CASE(EULER2D_IAMT)
       CALL Convert2Conservative_euler2Dia(this,Mesh,i1,i2,j1,j2,pvar,cvar)
    CASE(EULER2D_IAMROT)
       CALL Convert2Conservative_euler2Diar(this,Mesh,i1,i2,j1,j2,pvar,cvar)
    CASE(EULER2D_ISOIAMT)
       CALL Convert2Cons_euler2Ditia(this,Mesh,i1,i2,j1,j2,pvar,cvar)
    CASE(EULER2D_ISOIAMROT)
       CALL Convert2Cons_euler2Ditiar(this,Mesh,i1,i2,j1,j2,pvar,cvar)
    END SELECT
  END SUBROUTINE Convert2Conservative_centsub


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
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL Convert2Conservative_euler2D(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,prim,cons)
    CASE(EULER2D_ISOTHERM)
       CALL Convert2Conservative_euler2Dit(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,prim,cons)
    CASE(EULER3D_ROTSYM)
       CALL Convert2Conservative_euler3Drs(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,prim,cons)
    CASE(EULER3D_ROTAMT)
       CALL Convert2Conservative_euler3Dra(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,prim,cons)
    CASE(EULER3D_ROTSYMSGS)
       CALL Convert2Cons_euler3DrsSGS(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,prim,cons)
    CASE(EULER3D_ROTAMTSGS)
       CALL Convert2Cons_euler3DraSGS(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,prim,cons)
    CASE(EULER2D_SGS)
       CALL Convert2Conservative_euler2Dsgs(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,prim,cons)
    CASE(EULER2D_IAMT)
       CALL Convert2Conservative_euler2Dia(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,prim,cons)
    CASE(EULER2D_IAMROT)
       CALL Convert2Conservative_euler2Diar(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,prim,cons)
    CASE(EULER2D_ISOIAMT)
       CALL Convert2Cons_euler2Ditia(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,prim,cons)
    CASE(EULER2D_ISOIAMROT)
       CALL Convert2Cons_euler2Ditiar(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
            Mesh%JGMIN,Mesh%JGMAX,prim,cons)
    END SELECT
  END SUBROUTINE Convert2Conservative_faces
  

  PURE SUBROUTINE Convert2Conservative_facesub(this,Mesh,i1,i2,j1,j2,prim,cons)
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
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL Convert2Conservative_euler2D(this,Mesh,i1,i2,j1,j2,prim,cons)
    CASE(EULER2D_ISOTHERM)
       CALL Convert2Conservative_euler2Dit(this,Mesh,i1,i2,j1,j2,prim,cons)
    CASE(EULER3D_ROTSYM)
       CALL Convert2Conservative_euler3Drs(this,Mesh,i1,i2,j1,j2,prim,cons)
    CASE(EULER3D_ROTAMT)
       CALL Convert2Conservative_euler3Dra(this,Mesh,i1,i2,j1,j2,prim,cons)
    CASE(EULER3D_ROTSYMSGS)
       CALL Convert2Cons_euler3DrsSGS(this,Mesh,i1,i2,j1,j2,prim,cons)
    CASE(EULER3D_ROTAMTSGS)
       CALL Convert2Cons_euler3DraSGS(this,Mesh,i1,i2,j1,j2,prim,cons)
    CASE(EULER2D_SGS)
       CALL Convert2Conservative_euler2Dsgs(this,Mesh,i1,i2,j1,j2,prim,cons)
    CASE(EULER2D_IAMT)
       CALL Convert2Conservative_euler2Dia(this,Mesh,i1,i2,j1,j2,prim,cons)
    CASE(EULER2D_IAMROT)
       CALL Convert2Conservative_euler2Diar(this,Mesh,i1,i2,j1,j2,prim,cons)
    CASE(EULER2D_ISOIAMT)
       CALL Convert2Cons_euler2Ditia(this,Mesh,i1,i2,j1,j2,prim,cons)
    CASE(EULER2D_ISOIAMROT)
       CALL Convert2Cons_euler2Ditiar(this,Mesh,i1,i2,j1,j2,prim,cons)
    END SELECT
  END SUBROUTINE Convert2Conservative_facesub
  

  PURE SUBROUTINE ReflectionMasks(this,reflX,reflY)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    LOGICAL, DIMENSION(this%vnum) :: reflX,reflY
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this
    INTENT(OUT)       :: reflX,reflY
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL ReflectionMasks_euler2D(this,reflX,reflY)
    CASE(EULER2D_ISOTHERM)
       CALL ReflectionMasks_euler2Dit(this,reflX,reflY)
    CASE(EULER3D_ROTSYM)
       CALL ReflectionMasks_euler3Drs(this,reflX,reflY)
    CASE(EULER3D_ROTAMT)
       CALL ReflectionMasks_euler3Dra(this,reflX,reflY)
    CASE(EULER3D_ROTSYMSGS)
       CALL ReflectionMasks_euler3DrsSGS(this,reflX,reflY)
    CASE(EULER3D_ROTAMTSGS)
       CALL ReflectionMasks_euler3DraSGS(this,reflX,reflY)
    CASE(EULER2D_SGS)
       CALL ReflectionMasks_euler2Dsgs(this,reflX,reflY)
    CASE(EULER2D_IAMT)
       CALL ReflectionMasks_euler2Dia(this,reflX,reflY)
    CASE(EULER2D_IAMROT)
       CALL ReflectionMasks_euler2Diar(this,reflX,reflY)
    CASE(EULER2D_ISOIAMT)
       CALL ReflectionMasks_euler2Ditia(this,reflX,reflY)
    CASE(EULER2D_ISOIAMROT)
       CALL ReflectionMasks_euler2Ditiar(this,reflX,reflY)
    END SELECT
  END SUBROUTINE ReflectionMasks


  SUBROUTINE ClosePhysics(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    IF (.NOT.Initialized(this)) &
        CALL Error(this,"ClosePhysics","not initialized")
    ! deallocate pointer variables used in all physics modules
    DEALLOCATE(this%amin,this%amax,this%bmin,this%bmax,this%tmp, &
         this%tmin,this%tmax,this%bccsound,this%fcsound)
    ! call specific dallocation procedures
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL ClosePhysics_euler2D(this)
    CASE(EULER2D_ISOTHERM)
       CALL ClosePhysics_euler2Dit(this)
    CASE(EULER3D_ROTSYM)
       CALL ClosePhysics_euler3Drs(this)
    CASE(EULER3D_ROTAMT)
       CALL ClosePhysics_euler3Dra(this)
    CASE(EULER3D_ROTSYMSGS)
       CALL ClosePhysics_euler3DrsSGS(this)
   CASE(EULER3D_ROTAMTSGS)
       CALL ClosePhysics_euler3DraSGS(this)
    CASE(EULER2D_SGS)
       CALL ClosePhysics_euler2Dsgs(this) 
    CASE(EULER2D_IAMT)
       CALL ClosePhysics_euler2Dia(this)
    CASE(EULER2D_IAMROT)
       CALL ClosePhysics_euler2Diar(this)
    CASE(EULER2D_ISOIAMT)
       CALL ClosePhysics_euler2Ditia(this)
    CASE(EULER2D_ISOIAMROT)
       CALL ClosePhysics_euler2Ditiar(this)
    END SELECT
  END SUBROUTINE ClosePhysics

END MODULE physics_generic
