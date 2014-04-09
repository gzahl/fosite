!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler3Drotsym.f90                                         #
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
! basic module for 3D Euler equations with rotational symmetry
!----------------------------------------------------------------------------!
MODULE physics_euler3Drotsym
  USE physics_common
  USE physics_euler2D, ReflectionMasks_euler3Drs => ReflectionMasks_euler2D, &
       CheckData_euler3Drs => CheckData_euler2D, &
       CalculateWaveSpeeds_euler3Drs => CalculateWaveSpeeds_euler2D
  USE mesh_common, ONLY : Mesh_TYP
  USE sources_common, ONLY : Sources_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTERFACE GeometricalSources_euler3Drs
     MODULE PROCEDURE GeometricalSources_center
     MODULE PROCEDURE GeometricalSources_faces
  END INTERFACE
  INTERFACE Convert2Primitive_euler3Drs
     MODULE PROCEDURE Convert2Primitive_center
     MODULE PROCEDURE Convert2Primitive_faces
  END INTERFACE
  INTERFACE Convert2Conservative_euler3Drs
     MODULE PROCEDURE Convert2Conservative_center
     MODULE PROCEDURE Convert2Conservative_faces
  END INTERFACE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: num_var = 5              ! number of variables       !
  REAL, PARAMETER :: TINY = 1.0E-30              ! to avoid division by 0    !
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 3D w/ rot. symmetry"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitPhysics_euler3Drs, &
       MallocPhysics_euler3Drs, &
       CheckData_euler3Drs, &
       CalculateWaveSpeeds_euler3Drs, &
       CalculateFluxesX_euler3Drs, &
       CalculateFluxesY_euler3Drs, &
       GeometricalSources_euler3Drs, &
       ExternalSources_euler3Drs, &
       ViscositySources_euler3Drs, &
       Convert2Primitive_euler3Drs, &
       Convert2Conservative_euler3Drs, &
       ReflectionMasks_euler3Drs, &
       AxisMasks_euler3Drs, &
       ClosePhysics_euler3Drs
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitPhysics_euler3Drs(this,problem)
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
    this%ZVELOCITY = 4                                 ! rotational velocity !
    this%ZMOMENTUM = 4                                 ! rotational momentum !
    ! set names for primitive and conservative variables
    this%pvarname(this%DENSITY)   = "density"
    this%pvarname(this%XVELOCITY) = "x-velocity"
    this%pvarname(this%YVELOCITY) = "y-velocity"
    this%pvarname(this%ZVELOCITY) = "z-velocity"
    this%pvarname(this%PRESSURE)  = "pressure"
    this%cvarname(this%DENSITY)   = "density"
    this%cvarname(this%XMOMENTUM) = "x-momentum"
    this%cvarname(this%YMOMENTUM) = "y-momentum"
    this%cvarname(this%ZMOMENTUM) = "z-momentum"
    this%cvarname(this%ENERGY)    = "energy"
  END SUBROUTINE InitPhysics_euler3Drs


  SUBROUTINE MallocPhysics_euler3Drs(this,Mesh)
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
       CALL Error(this, "MallocPhysics_euler3Drs", "Unable to allocate memory.")
    END IF
  END SUBROUTINE MallocPhysics_euler3Drs


  PURE SUBROUTINE CalculateFluxesX_euler3Drs(this,Mesh,nmin,nmax,prim,cons,xfluxes)
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
    CALL CalculateFlux_euler3Drs(prim(:,:,nmin:nmax,this%DENSITY), &
         prim(:,:,nmin:nmax,this%XVELOCITY),prim(:,:,nmin:nmax,this%PRESSURE),&
         cons(:,:,nmin:nmax,this%XMOMENTUM),cons(:,:,nmin:nmax,this%YMOMENTUM), &
         cons(:,:,nmin:nmax,this%ZMOMENTUM),cons(:,:,nmin:nmax,this%ENERGY), &
         xfluxes(:,:,nmin:nmax,this%DENSITY),xfluxes(:,:,nmin:nmax,this%XMOMENTUM), &
         xfluxes(:,:,nmin:nmax,this%YMOMENTUM),xfluxes(:,:,nmin:nmax,this%ZMOMENTUM), &
         xfluxes(:,:,nmin:nmax,this%ENERGY))
   END SUBROUTINE CalculateFluxesX_euler3Drs


  PURE SUBROUTINE CalculateFluxesY_euler3Drs(this,Mesh,nmin,nmax,prim,cons,yfluxes)
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
    CALL CalculateFlux_euler3Drs(prim(:,:,nmin:nmax,this%DENSITY), &
         prim(:,:,nmin:nmax,this%YVELOCITY),prim(:,:,nmin:nmax,this%PRESSURE),&
         cons(:,:,nmin:nmax,this%YMOMENTUM),cons(:,:,nmin:nmax,this%XMOMENTUM), &
         cons(:,:,nmin:nmax,this%ZMOMENTUM),cons(:,:,nmin:nmax,this%ENERGY), &
         yfluxes(:,:,nmin:nmax,this%DENSITY),yfluxes(:,:,nmin:nmax,this%YMOMENTUM), &
         yfluxes(:,:,nmin:nmax,this%XMOMENTUM),yfluxes(:,:,nmin:nmax,this%ZMOMENTUM), &
         yfluxes(:,:,nmin:nmax,this%ENERGY))
  END SUBROUTINE CalculateFluxesY_euler3Drs


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
         Mesh%czxz(:,:,1),Mesh%czyz(:,:,1),this%fcent(:,:,1,1),this%fcent(:,:,2,1))

    ! geometrical source terms in momentum equationes
    ! with centrifugal forces
    sterm(:,:,this%XMOMENTUM) = MomentumSourcesX_euler2D(cvar(:,:,this%YMOMENTUM), &
         pvar(:,:,this%XVELOCITY),pvar(:,:,this%YVELOCITY),pvar(:,:,this%PRESSURE), &
         Mesh%cxyx(:,:,1),Mesh%cyxy(:,:,1),Mesh%czxz(:,:,1)) + this%fcent(:,:,1,1)
    sterm(:,:,this%YMOMENTUM) = MomentumSourcesY_euler2D(cvar(:,:,this%XMOMENTUM), &
         pvar(:,:,this%XVELOCITY),pvar(:,:,this%YVELOCITY),pvar(:,:,this%PRESSURE), &
         Mesh%cxyx(:,:,1),Mesh%cyxy(:,:,1),Mesh%czyz(:,:,1)) + this%fcent(:,:,2,1)
    sterm(:,:,this%ZMOMENTUM) = -pvar(:,:,this%ZVELOCITY)*(Mesh%czxz(:,:,1) * &
         cvar(:,:,this%XMOMENTUM) + Mesh%czyz(:,:,1) * cvar(:,:,this%YMOMENTUM))
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
    INTENT(IN)        :: Mesh,prim,cons
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! calculate centrifugal forces
    CALL CentrifugalForces(prim(:,:,:,this%DENSITY),prim(:,:,:,this%ZVELOCITY), &
         Mesh%czxz(:,:,:),Mesh%czyz(:,:,:),this%fcent(:,:,:,1),this%fcent(:,:,:,2))

    ! geometrical source terms in momentum equationes
    ! sum up all four corner values
    sterm(:,:,this%XMOMENTUM) = SUM(MomentumSourcesX_euler2D(cons(:,:,:,this%YMOMENTUM), &
         prim(:,:,:,this%XVELOCITY),prim(:,:,:,this%YVELOCITY),prim(:,:,:,this%PRESSURE), &
         Mesh%cxyx(:,:,:),Mesh%cyxy(:,:,:),Mesh%czxz(:,:,:)) + this%fcent(:,:,:,1),DIM=3)
    sterm(:,:,this%YMOMENTUM) = SUM(MomentumSourcesY_euler2D(cons(:,:,:,this%XMOMENTUM), &
         prim(:,:,:,this%XVELOCITY),prim(:,:,:,this%YVELOCITY),prim(:,:,:,this%PRESSURE), &
         Mesh%cxyx(:,:,:),Mesh%cyxy(:,:,:),Mesh%czyz(:,:,:)) + this%fcent(:,:,:,2),DIM=3)
    sterm(:,:,this%ZMOMENTUM) = SUM(-prim(:,:,:,this%ZVELOCITY)*(Mesh%czxz(:,:,:) * &
         cons(:,:,:,this%XMOMENTUM) + Mesh%czyz(:,:,:) * cons(:,:,:,this%YMOMENTUM)),DIM=3)
  END SUBROUTINE GeometricalSources_faces


  PURE SUBROUTINE ExternalSources_euler3Drs(this,Mesh,accel,pvar,cvar,sterm)
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
         :: accel
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
         :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,accel,pvar,cvar
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    sterm(:,:,this%DENSITY)   = 0.
    sterm(:,:,this%XMOMENTUM) = pvar(:,:,this%DENSITY) * accel(:,:,1)
    sterm(:,:,this%YMOMENTUM) = pvar(:,:,this%DENSITY) * accel(:,:,2)
    sterm(:,:,this%ZMOMENTUM) = 0.
    sterm(:,:,this%ENERGY)    = cvar(:,:,this%XMOMENTUM) * accel(:,:,1) + &
         cvar(:,:,this%YMOMENTUM) * accel(:,:,2)
  END SUBROUTINE ExternalSources_euler3Drs


 PURE SUBROUTINE ViscositySources_euler3Drs(this,Mesh,Sources,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Sources_TYP) :: Sources
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: &
         pvar,cvar,sterm
   !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,pvar,cvar
    INTENT(INOUT)     :: Sources
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! bulk viscosity contribution to the stress tensor
!CDIR NOUNROLL
    DO j=Mesh%JMIN-1,Mesh%JMAX+1
!CDIR NODEP
       DO i=Mesh%IMIN-1,Mesh%IMAX+1
          ! using btxy as temporary storage
          Sources%btxy(i,j) = Sources%bulkvis(i,j) * 0.5 * (&
                ( Mesh%fhy(i,j,2) * Mesh%fhz(i,j,2) * &
                ( pvar(i+1,j,this%XVELOCITY)+pvar(i,j,this%XVELOCITY) ) &
                - Mesh%fhy(i,j,1) * Mesh%fhz(i,j,1) * &
                ( pvar(i,j,this%XVELOCITY)+pvar(i-1,j,this%XVELOCITY) ) ) &
                * Mesh%dydV(i,j) &
                +(Mesh%fhx(i,j,4) * Mesh%fhz(i,j,4) * &
                ( pvar(i,j+1,this%YVELOCITY)+pvar(i,j,this%YVELOCITY) ) &
                - Mesh%fhx(i,j,3) * Mesh%fhz(i,j,3) * &
                ( pvar(i,j,this%YVELOCITY)+pvar(i,j-1,this%YVELOCITY) ) ) &
                * Mesh%dxdV(i,j) )
       END DO
    END DO

    ! stress tensor at cell bary centers
!CDIR NOUNROLL
    DO j=Mesh%JMIN-1,Mesh%JMAX+1
!CDIR NODEP
       DO i=Mesh%IMIN-1,Mesh%IMAX+1
          Sources%btxx(i,j) = Sources%dynvis(i,j) * &
               ( (pvar(i+1,j,this%XVELOCITY) - pvar(i-1,j,this%XVELOCITY)) / Mesh%dlx(i,j) &
               + 2.0 * Mesh%cxyx(i,j,1) * pvar(i,j,this%YVELOCITY) ) &
               + Sources%btxy(i,j) ! bulk viscosity contribution
               
          Sources%btyy(i,j) = Sources%dynvis(i,j) * &
               ( (pvar(i,j+1,this%YVELOCITY) - pvar(i,j-1,this%YVELOCITY)) / Mesh%dly(i,j) &
               + 2.0 * Mesh%cyxy(i,j,1) * pvar(i,j,this%XVELOCITY) ) &
               + Sources%btxy(i,j) ! bulk viscosity contribution

          Sources%btzz(i,j) = Sources%dynvis(i,j) * &
               ( 2.0 * ( Mesh%czxz(i,j,1) * pvar(i,j,this%XVELOCITY) ) &
               + Mesh%czyz(i,j,1) * pvar(i,j,this%YVELOCITY) ) &
               + Sources%btxy(i,j) ! bulk viscosity contribution

          Sources%btxy(i,j) = Sources%dynvis(i,j) * ( 0.5 * &
               ( (pvar(i+1,j,this%YVELOCITY) - pvar(i-1,j,this%YVELOCITY)) / Mesh%dlx(i,j) &
               + (pvar(i,j+1,this%XVELOCITY) - pvar(i,j-1,this%XVELOCITY)) / Mesh%dly(i,j) ) &
               - Mesh%cxyx(i,j,1) * pvar(i,j,this%XVELOCITY) &
               - Mesh%cyxy(i,j,1) * pvar(i,j,this%YVELOCITY) )

          Sources%btxz(i,j) = Sources%dynvis(i,j) * ( 0.5 * &
               ( (pvar(i+1,j,this%ZVELOCITY) - pvar(i-1,j,this%ZVELOCITY)) / Mesh%dlx(i,j) ) &
               - Mesh%czxz(i,j,1) * pvar(i,j,this%ZVELOCITY) )

          Sources%btyz(i,j) = Sources%dynvis(i,j) * ( 0.5 * &
               ( (pvar(i,j+1,this%ZVELOCITY) - pvar(i,j-1,this%ZVELOCITY)) / Mesh%dly(i,j) ) &
               - Mesh%czyz(i,j,1) * pvar(i,j,this%ZVELOCITY) )
       END DO
    END DO

    ! mean values of stress tensor components across the cell interfaces
!CDIR NOUNROLL
    DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
       DO i=Mesh%IMIN,Mesh%IMAX
          Sources%ftxx(i,j,1) = 0.5 * ( Sources%btxx(i-1,j) + Sources%btxx(i,j) )
          Sources%ftxx(i,j,2) = 0.5 * ( Sources%btxx(i+1,j) + Sources%btxx(i,j) )
          Sources%ftxx(i,j,3) = 0.5 * ( Sources%btxx(i,j-1) + Sources%btxx(i,j) )
          Sources%ftxx(i,j,4) = 0.5 * ( Sources%btxx(i,j+1) + Sources%btxx(i,j) )

          Sources%ftyy(i,j,1) = 0.5 * ( Sources%btyy(i-1,j) + Sources%btyy(i,j) )
          Sources%ftyy(i,j,2) = 0.5 * ( Sources%btyy(i+1,j) + Sources%btyy(i,j) )
          Sources%ftyy(i,j,3) = 0.5 * ( Sources%btyy(i,j-1) + Sources%btyy(i,j) )
          Sources%ftyy(i,j,4) = 0.5 * ( Sources%btyy(i,j+1) + Sources%btyy(i,j) )

          Sources%ftzz(i,j,1) = 0.5 * ( Sources%btzz(i-1,j) + Sources%btzz(i,j) )
          Sources%ftzz(i,j,2) = 0.5 * ( Sources%btzz(i+1,j) + Sources%btzz(i,j) )
          Sources%ftzz(i,j,3) = 0.5 * ( Sources%btzz(i,j-1) + Sources%btzz(i,j) )
          Sources%ftzz(i,j,4) = 0.5 * ( Sources%btzz(i,j+1) + Sources%btzz(i,j) )

          Sources%ftxy(i,j,1) = 0.5 * ( Sources%btxy(i-1,j) + Sources%btxy(i,j) )
          Sources%ftxy(i,j,2) = 0.5 * ( Sources%btxy(i+1,j) + Sources%btxy(i,j) )
          Sources%ftxy(i,j,3) = 0.5 * ( Sources%btxy(i,j-1) + Sources%btxy(i,j) )
          Sources%ftxy(i,j,4) = 0.5 * ( Sources%btxy(i,j+1) + Sources%btxy(i,j) )

          Sources%ftxz(i,j,1) = 0.5 * ( Sources%btxz(i-1,j) + Sources%btxz(i,j) )
          Sources%ftxz(i,j,2) = 0.5 * ( Sources%btxz(i+1,j) + Sources%btxz(i,j) )
          Sources%ftxz(i,j,3) = 0.5 * ( Sources%btxz(i,j-1) + Sources%btxz(i,j) )
          Sources%ftxz(i,j,4) = 0.5 * ( Sources%btxz(i,j+1) + Sources%btxz(i,j) )
          
          Sources%ftyz(i,j,1) = 0.5 * ( Sources%btyz(i-1,j) + Sources%btyz(i,j) )
          Sources%ftyz(i,j,2) = 0.5 * ( Sources%btyz(i+1,j) + Sources%btyz(i,j) )
          Sources%ftyz(i,j,3) = 0.5 * ( Sources%btyz(i,j-1) + Sources%btyz(i,j) )
          Sources%ftyz(i,j,4) = 0.5 * ( Sources%btyz(i,j+1) + Sources%btyz(i,j) )
       END DO
    END DO

    ! viscosity source terms
    ! (a) momentum sources
    sterm(:,:,this%XMOMENTUM) = Mesh%dydV(:,:) * &
         ( Mesh%fhy(:,:,2) * Mesh%fhz(:,:,2) * Sources%ftxx(:,:,2) &
         - Mesh%fhy(:,:,1) * Mesh%fhz(:,:,1) * Sources%ftxx(:,:,1) &
         - Mesh%bhz(:,:) * (Mesh%fhy(:,:,2) - Mesh%fhy(:,:,1)) * Sources%btyy(:,:) &
         - Mesh%bhy(:,:) * (Mesh%fhz(:,:,2) - Mesh%fhz(:,:,1)) * Sources%btzz(:,:) ) &
         + Mesh%dxdV(:,:) * &
         ( Mesh%fhx(:,:,4) * Mesh%fhz(:,:,4) * Sources%ftxy(:,:,4) &
         - Mesh%fhx(:,:,3) * Mesh%fhz(:,:,3) * Sources%ftxy(:,:,3) &
         + Mesh%bhz(:,:) * (Mesh%fhx(:,:,4) - Mesh%fhx(:,:,3)) * Sources%btxy(:,:) )

    sterm(:,:,this%YMOMENTUM) = Mesh%dydV(:,:) * &
         ( Mesh%fhy(:,:,2) * Mesh%fhz(:,:,2) * Sources%ftxy(:,:,2) &
         - Mesh%fhy(:,:,1) * Mesh%fhz(:,:,1) * Sources%ftxy(:,:,1) &
         + Mesh%bhz(:,:) * (Mesh%fhy(:,:,2) - Mesh%fhy(:,:,1)) * Sources%btxy(:,:) ) &
         + Mesh%dxdV(:,:) * &
         ( Mesh%fhx(:,:,4) * Mesh%fhz(:,:,4) * Sources%ftyy(:,:,4) &
         - Mesh%fhx(:,:,3) * Mesh%fhz(:,:,3) * Sources%ftyy(:,:,3) &
         - Mesh%bhz(:,:) * (Mesh%fhx(:,:,4) - Mesh%fhx(:,:,3)) * Sources%btxx(:,:) &
         - Mesh%bhx(:,:) * (Mesh%fhz(:,:,4) - Mesh%fhz(:,:,3)) * Sources%btzz(:,:) )
    
    sterm(:,:,this%ZMOMENTUM) = Mesh%dydV(:,:) * &
         ( Mesh%fhy(:,:,2) * Mesh%fhz(:,:,2) * Sources%ftxz(:,:,2) &
         - Mesh%fhy(:,:,1) * Mesh%fhz(:,:,1) * Sources%ftxz(:,:,1) &
         + Mesh%bhy(:,:) * ( Mesh%fhz(:,:,2) - Mesh%fhz(:,:,1) ) * Sources%btxz(:,:) ) &
         + Mesh%dxdV(:,:) * &
         ( Mesh%fhx(:,:,4) * Mesh%fhz(:,:,4) * Sources%ftyz(:,:,4) &
         - Mesh%fhx(:,:,3) * Mesh%fhz(:,:,3) * Sources%ftyz(:,:,3) &
         + Mesh%bhx(:,:) * (Mesh%fhz(:,:,4) - Mesh%fhz(:,:,3)) * Sources%btyz(:,:) )

    ! (b) energy sources
!CDIR NOUNROLL
    DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
       DO i=Mesh%IMIN,Mesh%IMAX
          sterm(i,j,this%ENERGY) = 0.5 * (&
                 Mesh%dydV(i,j) * ( Mesh%fhy(i,j,2) * Mesh%fhz(i,j,2) * &
               ( (pvar(i+1,j,this%XVELOCITY)+pvar(i,j,this%XVELOCITY)) * Sources%ftxx(i,j,2) &
               + (pvar(i+1,j,this%YVELOCITY)+pvar(i,j,this%YVELOCITY)) * Sources%ftxy(i,j,2) &
               + (pvar(i+1,j,this%ZVELOCITY)+pvar(i,j,this%ZVELOCITY)) * Sources%ftxz(i,j,2) ) &
               - Mesh%fhy(i,j,1) * Mesh%fhz(i,j,1) * &
               ( (pvar(i,j,this%XVELOCITY)+pvar(i-1,j,this%XVELOCITY)) * Sources%ftxx(i,j,1) &
               + (pvar(i,j,this%YVELOCITY)+pvar(i-1,j,this%YVELOCITY)) * Sources%ftxy(i,j,1) &
               + (pvar(i,j,this%ZVELOCITY)+pvar(i-1,j,this%ZVELOCITY)) * Sources%ftxz(i,j,1) ) ) &
               + Mesh%dxdV(i,j) * ( Mesh%fhx(i,j,4) * Mesh%fhz(i,j,4) * &
               ( (pvar(i,j+1,this%XVELOCITY)+pvar(i,j,this%XVELOCITY)) * Sources%ftxy(i,j,4) &
               + (pvar(i,j+1,this%YVELOCITY)+pvar(i,j,this%YVELOCITY)) * Sources%ftyy(i,j,4) &
               + (pvar(i,j+1,this%ZVELOCITY)+pvar(i,j,this%ZVELOCITY)) * Sources%ftyz(i,j,4) ) &
               - Mesh%fhx(i,j,3) * Mesh%fhz(i,j,3) * &
               ( (pvar(i,j,this%XVELOCITY)+pvar(i,j-1,this%XVELOCITY)) * Sources%ftxy(i,j,3) &
               + (pvar(i,j,this%YVELOCITY)+pvar(i,j-1,this%YVELOCITY)) * Sources%ftyy(i,j,3) &
               + (pvar(i,j,this%ZVELOCITY)+pvar(i,j-1,this%ZVELOCITY)) * Sources%ftyz(i,j,3) ) ) )
       END DO
    END DO
  END SUBROUTINE ViscositySources_euler3Drs


  PURE SUBROUTINE Convert2Primitive_center(this,Mesh,cvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
                         :: cvar,pvar
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,cvar
    INTENT(OUT) :: pvar
    !------------------------------------------------------------------------!
    CALL Cons2Prim_euler3Drs(this%gamma,cvar(:,:,this%DENSITY),cvar(:,:,this%XMOMENTUM), &
         cvar(:,:,this%YMOMENTUM),cvar(:,:,this%ZMOMENTUM),cvar(:,:,this%ENERGY), &
         pvar(:,:,this%DENSITY),pvar(:,:,this%XVELOCITY),pvar(:,:,this%YVELOCITY), &
         pvar(:,:,this%ZVELOCITY),pvar(:,:,this%PRESSURE))
  END SUBROUTINE Convert2Primitive_center


  PURE SUBROUTINE Convert2Primitive_faces(this,Mesh,cons,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%VNUM) &
                         :: cons,prim
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,cons
    INTENT(OUT) :: prim
    !------------------------------------------------------------------------!
    CALL Cons2Prim_euler3Drs(this%gamma,cons(:,:,:,this%DENSITY),cons(:,:,:,this%XMOMENTUM), &
         cons(:,:,:,this%YMOMENTUM),cons(:,:,:,this%ZMOMENTUM),cons(:,:,:,this%ENERGY), &
         prim(:,:,:,this%DENSITY),prim(:,:,:,this%XVELOCITY),prim(:,:,:,this%YVELOCITY), &
         prim(:,:,:,this%ZVELOCITY),prim(:,:,:,this%PRESSURE))
  END SUBROUTINE Convert2Primitive_faces


  PURE SUBROUTINE Convert2Conservative_center(this,Mesh,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
                         :: cvar,pvar
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,pvar
    INTENT(OUT) :: cvar
    !------------------------------------------------------------------------!
    CALL Prim2Cons_euler3Drs(this%gamma,pvar(:,:,this%DENSITY),pvar(:,:,this%XVELOCITY), &
         pvar(:,:,this%YVELOCITY),pvar(:,:,this%ZVELOCITY),pvar(:,:,this%PRESSURE), &
         cvar(:,:,this%DENSITY),cvar(:,:,this%XMOMENTUM),cvar(:,:,this%YMOMENTUM), &
         cvar(:,:,this%ZMOMENTUM),cvar(:,:,this%ENERGY))
  END SUBROUTINE Convert2Conservative_center


  PURE SUBROUTINE Convert2Conservative_faces(this,Mesh,prim,cons)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%VNUM) &
                         :: cons,prim
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,prim
    INTENT(OUT) :: cons
    !------------------------------------------------------------------------!
    CALL Prim2Cons_euler3Drs(this%gamma,prim(:,:,:,this%DENSITY),prim(:,:,:,this%XVELOCITY), &
         prim(:,:,:,this%YVELOCITY),prim(:,:,:,this%ZVELOCITY),prim(:,:,:,this%PRESSURE), &
         cons(:,:,:,this%DENSITY),cons(:,:,:,this%XMOMENTUM),cons(:,:,:,this%YMOMENTUM), &
         cons(:,:,:,this%ZMOMENTUM),cons(:,:,:,this%ENERGY))
  END SUBROUTINE Convert2Conservative_faces


  PURE SUBROUTINE AxisMasks_euler3Drs(this,reflX,reflY)
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
    reflX(this%XVELOCITY) = .TRUE.
    reflX(this%ZVELOCITY) = .TRUE.
    reflY(this%YVELOCITY) = .TRUE.
    reflY(this%ZVELOCITY) = .TRUE.
  END SUBROUTINE AxisMasks_euler3Drs


  ELEMENTAL SUBROUTINE CalculateFlux_euler3Drs(rho,v,P,m1,m2,m3,E,f1,f2,f3,f4,f5)
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
  END SUBROUTINE CalculateFlux_euler3Drs


  ELEMENTAL SUBROUTINE CentrifugalForces(rho,vz,czxz,czyz,fcx,fcy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho,vz,czxz,czyz
    REAL, INTENT(OUT) :: fcx,fcy
    !------------------------------------------------------------------------!
    REAL :: tmp2
    !------------------------------------------------------------------------!
    tmp2 = rho*vz*vz
    fcx  = tmp2 * czxz
    fcy  = tmp2 * czyz
  END SUBROUTINE CentrifugalForces


  ELEMENTAL SUBROUTINE Cons2Prim_euler3Drs(gamma,rho_in,mu,mv,mw,E,rho_out,u,v,w,P)
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
  END SUBROUTINE Cons2Prim_euler3Drs

  
  ELEMENTAL SUBROUTINE Prim2Cons_euler3Drs(gamma,rho_in,u,v,w,P,rho_out,mu,mv,mw,E)
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
  END SUBROUTINE Prim2Cons_euler3Drs


  SUBROUTINE ClosePhysics_euler3Drs(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL ClosePhysics_euler2D(this)
    DEALLOCATE(this%fcent)
  END SUBROUTINE ClosePhysics_euler3Drs


END MODULE physics_euler3Drotsym
