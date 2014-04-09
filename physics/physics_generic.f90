!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_generic.f90                                               #
!#                                                                           #
!# Copyright (C) 2007 - 2012                                                 #
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
! generic module for the advection problem
!----------------------------------------------------------------------------!
MODULE physics_generic
  USE physics_euler2Disothm, ClosePhysics_common => ClosePhysics
  USE physics_euler2D
  USE physics_euler3Drotsym
  USE physics_euler3Drotamt
  USE constants_generic
  USE sources_common, ONLY : Sources_TYP
  USE mesh_common, ONLY : Mesh_TYP, Initialized
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
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
  !--------------------------------------------------------------------------!
  ! flags for advection problems
  INTEGER, PARAMETER :: EULER2D          = 1
  INTEGER, PARAMETER :: EULER2D_ISOTHERM = 2
  INTEGER, PARAMETER :: EULER3D_ROTSYM   = 3
  INTEGER, PARAMETER :: EULER3D_ROTAMT   = 4
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Physics_TYP, &
       ! constants
       EULER2D, EULER2D_ISOTHERM, EULER3D_ROTSYM, EULER3D_ROTAMT, &
       SI, CGS, GEOMETRICAL, &
       ! methods
       InitPhysics, &
       CheckData, &
       CalculateWaveSpeeds, &
       MaxWaveSpeeds, &
       CalculateFluxesX, &
       CalculateFluxesY, &
       CalculateStresses, &
       GeometricalSources, &
       ExternalSources, &
       ViscositySources, &
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

  SUBROUTINE InitPhysics(this,Mesh,problem,units,gamma,mu,cs,rhomin,pmin,dpmax)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: problem
    INTEGER, OPTIONAL :: units
    REAL, OPTIONAL    :: gamma,mu,cs,rhomin,pmin,dpmax
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,problem,units,gamma,mu,rhomin,pmin,dpmax
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! check initialization of Mesh
    IF (.NOT.Initialized(Mesh)) &
         CALL Error(this,"InitPhysics","mesh module uninitialized")

    ! units
    IF (.NOT.Initialized(this%constants)) THEN
       IF (PRESENT(units)) THEN
          CALL InitConstants(this%constants,units)
       ELSE
          CALL InitConstants(this%constants,SI)
       END IF
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

    ! isothermal sound speed
    IF (PRESENT(cs)) THEN
       this%csiso = cs
    ELSE
       this%csiso = 343.0 ! air at 293 K
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

     ! allocate memory for arrays common to all physics modules
    ALLOCATE(this%amin(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%amax(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%bmin(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%bmax(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%tmin(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4), &
         this%tmax(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4), &
         STAT = err)
    IF (err.NE.0) &
         CALL Error(this, "InitPhysics", "Unable to allocate memory.")

   SELECT CASE(problem)
    CASE(EULER2D)
       CALL InitPhysics_euler2D(this,Mesh,problem)
    CASE(EULER2D_ISOTHERM)
       CALL InitPhysics_euler2Dit(this,Mesh,problem)
    CASE(EULER3D_ROTSYM)
       CALL InitPhysics_euler3Drs(this,Mesh,problem)
    CASE(EULER3D_ROTAMT)
       CALL InitPhysics_euler3Dra(this,Mesh,problem)
    CASE DEFAULT
       CALL Error(this, "InitPhysics", "Unknown advection problem.")
    END SELECT

    NULLIFY(this%sources)

    ! print some information
    CALL Info(this, " PHYSICS--> advection problem: " // TRIM(GetName(this)))

  END SUBROUTINE InitPhysics


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
    END SELECT
  END SUBROUTINE AxisMasks


  PURE FUNCTION CheckData(this,Mesh,pvar,pold,meshrange) RESULT(bad_data)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
                      :: pvar,pold
    INTEGER, DIMENSION(4), OPTIONAL :: meshrange
    INTEGER           :: bad_data
    !------------------------------------------------------------------------!
    INTEGER, DIMENSION(4) :: meshrange_def
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,pvar,pold,meshrange
    !------------------------------------------------------------------------!
    IF (PRESENT(meshrange)) THEN
       meshrange_def = meshrange
    ELSE
       meshrange_def = (/Mesh%IGMIN, Mesh%IGMAX, Mesh%JGMIN, Mesh%JGMAX/)
    END IF

!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       bad_data = CheckData_euler2D(this,Mesh,pvar,pold,meshrange_def)
    CASE(EULER2D_ISOTHERM)
       bad_data = CheckData_euler2Dit(this,Mesh,pvar,pold,meshrange_def)
    CASE(EULER3D_ROTSYM)
       bad_data = CheckData_euler3Drs(this,Mesh,pvar,pold,meshrange_def)
    CASE(EULER3D_ROTAMT)
       bad_data = CheckData_euler3Dra(this,Mesh,pvar,pold,meshrange_def)
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
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalculateWaveSpeeds_euler2D(this,Mesh,prim)
    CASE(EULER2D_ISOTHERM)
       CALL CalculateWaveSpeeds_euler2Dit(this,Mesh,prim)
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
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalculateWaveSpeeds_euler2D(this,Mesh,pvar)
    CASE(EULER2D_ISOTHERM)
       CALL CalculateWaveSpeeds_euler2Dit(this,Mesh,pvar)
    CASE(EULER3D_ROTSYM)
       CALL CalculateWaveSpeeds_euler3Drs(this,Mesh,pvar)
    CASE(EULER3D_ROTAMT)
       CALL CalculateWaveSpeeds_euler3Dra(this,Mesh,pvar)
    END SELECT

    amax(:,:,1) = MAX(this%amax,-this%amin)
    amax(:,:,2) = MAX(this%bmax,-this%bmin)
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
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalculateFluxesX_euler2D(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    CASE(EULER2D_ISOTHERM)
       CALL CalculateFluxesX_euler2Dit(this,Mesh,nmin,nmax,prim,cons,xfluxes)
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
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(EULER2D)
       CALL CalculateFluxesY_euler2D(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    CASE(EULER2D_ISOTHERM)
       CALL CalculateFluxesY_euler2Dit(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    CASE(EULER3D_ROTSYM)
       CALL CalculateFluxesY_euler3Drs(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    CASE(EULER3D_ROTAMT)
       CALL CalculateFluxesY_euler3Dra(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    END SELECT
  END SUBROUTINE CalculateFluxesY


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
       CALL CalculateStresses_euler2D(this,Mesh,pvar,dynvis,bulkvis, &
            btxx,btxy,btyy)
    CASE(EULER2D_ISOTHERM)
       CALL CalculateStresses_euler2Dit(this,Mesh,pvar,dynvis,bulkvis, &
            btxx,btxy,btyy)
    CASE(EULER3D_ROTSYM)
       CALL CalculateStresses_euler3Drs(this,Mesh,pvar,dynvis,bulkvis, &
            btxx,btxy,btxz,btyy,btyz,btzz)
    CASE(EULER3D_ROTAMT)
! ***************************************************************************!
! FIXME: not implemented yet
!       CALL CalculateStresses_euler3Dra(this,Mesh,nmin,nmax,prim,cons,yfluxes)
! ***************************************************************************!
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
! ***************************************************************************!
! FIXME: not implemented yet
!       CALL ExternalSources_euler3Dra(this,Mesh,accel,pvar,cvar,sterm)
! ***************************************************************************!
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
    CASE(EULER3D_ROTAMT)
       CALL Convert2Conservative_euler3Dra(this,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
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
    DEALLOCATE(this%amin,this%amax,this%bmin,this%bmax, &
         this%tmin,this%tmax)
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
    END SELECT
  END SUBROUTINE ClosePhysics

END MODULE physics_generic
