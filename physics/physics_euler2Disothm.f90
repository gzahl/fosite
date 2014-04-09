!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler2Disothm.f90                                         #
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
! basic module for 2D isothermal Euler equations
!----------------------------------------------------------------------------!
MODULE physics_euler2Disothm
  USE physics_common
  USE physics_euler2D,  ReflectionMasks_euler2Dit =>  ReflectionMasks_euler2D, &
       AxisMasks_euler2Dit =>  AxisMasks_euler2D
  USE mesh_common, ONLY : Mesh_TYP
  USE sources_common, ONLY : Sources_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTERFACE CalculateWaveSpeeds_euler2Dit
     MODULE PROCEDURE CalculateWaveSpeeds_center
     MODULE PROCEDURE CalculateWaveSpeeds_faces
  END INTERFACE
     INTERFACE GeometricalSources_euler2Dit
     MODULE PROCEDURE GeometricalSources_center
     MODULE PROCEDURE GeometricalSources_faces
  END INTERFACE
  INTERFACE Convert2Primitive_euler2Dit
     MODULE PROCEDURE Convert2Primitive_center
     MODULE PROCEDURE Convert2Primitive_faces
  END INTERFACE
  INTERFACE Convert2Conservative_euler2Dit
     MODULE PROCEDURE Convert2Conservative_center
     MODULE PROCEDURE Convert2Conservative_faces
  END INTERFACE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: num_var = 3              ! number of variables       !
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 2D isotherm"
  CHARACTER(LEN=16), PARAMETER, DIMENSION(num_var) :: varname &
       = (/ "density         ", &
            "x-velocity      ", &
            "y-velocity      "/)
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Physics_TYP, &
       ! methods
       InitPhysics_euler2Dit, &
       ClosePhysics_euler2Dit, &
       CheckData_euler2Dit, &
       CalculateWaveSpeeds_euler2Dit, &
       CalculateFluxesX_euler2Dit, &
       CalculateFluxesY_euler2Dit, &
       GeometricalSources_euler2Dit, &
       ViscositySources_euler2Dit, &
       ExternalSources_euler2Dit, &
       Convert2Primitive_euler2Dit, &
       Convert2Conservative_euler2Dit, &
       ReflectionMasks_euler2Dit, &
       AxisMasks_euler2Dit, &
       SetWaveSpeeds_euler, &
       MomentumSourcesX_euler2Dit, &
       MomentumSourcesY_euler2Dit, &
       GetType, &
       GetName, &
       GetRank, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitPhysics_euler2Dit(this,problem)
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
    this%XVELOCITY = 2                                 ! x-velocity          !
    this%XMOMENTUM = 2                                 ! x-momentum          !
    this%YVELOCITY = num_var                           ! y-velocity          !
    this%YMOMENTUM = num_var                           ! y-momentum          !
    this%ZVELOCITY = 0                                 ! no z-velocity          !
    this%ZMOMENTUM = 0                                 ! no z-momentum          !
    this%PRESSURE  = 0                                 ! no pressure            !
    this%ENERGY    = 0                                 ! no total energy        !
    ! set names for primitive and conservative variables
    this%pvarname(this%DENSITY)   = "density"
    this%pvarname(this%XVELOCITY) = "x-velocity"
    this%pvarname(this%YVELOCITY) = "y-velocity"
    this%cvarname(this%DENSITY)   = "density"
    this%cvarname(this%XMOMENTUM) = "x-momentum"
    this%cvarname(this%YMOMENTUM) = "y-momentum"
  END SUBROUTINE InitPhysics_euler2Dit


  PURE FUNCTION CheckData_euler2Dit(this,Mesh,pvar,pold) RESULT (bad_data)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    ! checks for bad density or pressure on the whole grid including ghost cells
    ! return values:
    !   0 : valid density and pressure data
    !  -1 : density < this%rhomin, i.e. vacuum generated
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
         :: pvar,pold
    INTEGER           :: bad_data
    !------------------------------------------------------------------------!
    REAL              :: rhomin
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,pvar,pold
    !------------------------------------------------------------------------!
    ! density minimum
    rhomin = MINVAL(pvar(:,:,this%DENSITY))
    IF (rhomin.LE.this%rhomin) THEN
       bad_data = -1
       RETURN
    END IF
    ! everything ok
    bad_data = 0
  END FUNCTION CheckData_euler2Dit


  PURE SUBROUTINE CalculateWaveSpeeds_center(this,Mesh,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
         :: pvar
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! minimal and maximal wave speeds at cell interfaces
    ! x-direction
    CALL SetWaveSpeeds_euler(pvar(:,:,this%XVELOCITY),this%csiso, &
         this%tmin(:,:,1),this%tmax(:,:,1))
    ! y-direction
    CALL SetWaveSpeeds_euler(pvar(:,:,this%YVELOCITY),this%csiso, &
         this%tmin(:,:,2),this%tmax(:,:,2))
  END SUBROUTINE CalculateWaveSpeeds_center


  PURE SUBROUTINE CalculateWaveSpeeds_faces(this,Mesh,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%VNUM) &
         :: prim
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,prim
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! minimal and maximal wave speeds at cell interfaces
    CALL SetWaveSpeeds_euler(prim(:,:,1:2,this%XVELOCITY),this%csiso, &
         this%tmin(:,:,1:2),this%tmax(:,:,1:2))

    ! x-direction (west and east states)
    this%amin(Mesh%IMIN-1:Mesh%IMAX,:) = MIN(this%tmin(Mesh%IMIN:Mesh%IMAX+1,:,1), &
         this%tmin(Mesh%IMIN-1:Mesh%IMAX,:,2))
    this%amax(Mesh%IMIN-1:Mesh%IMAX,:) = MAX(this%tmax(Mesh%IMIN:Mesh%IMAX+1,:,1), &
         this%tmax(Mesh%IMIN-1:Mesh%IMAX,:,2))

    ! y-direction (south and north states)
    CALL SetWaveSpeeds_euler(prim(:,:,3:4,this%YVELOCITY),this%csiso, &
         this%tmin(:,:,1:2),this%tmax(:,:,1:2))

    this%bmin(:,Mesh%JMIN-1:Mesh%JMAX) = MIN(this%tmin(:,Mesh%JMIN:Mesh%JMAX+1,1), &
         this%tmin(:,Mesh%JMIN-1:Mesh%JMAX,2))
    this%bmax(:,Mesh%JMIN-1:Mesh%JMAX) = MAX(this%tmax(:,Mesh%JMIN:Mesh%JMAX+1,1), &
         this%tmax(:,Mesh%JMIN-1:Mesh%JMAX,2))
  END SUBROUTINE CalculateWaveSpeeds_faces


  PURE SUBROUTINE CalculateFluxesX_euler2Dit(this,Mesh,nmin,nmax,prim,cons,xfluxes)
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
    CALL CalculateFlux_euler2Dit(this%csiso,prim(:,:,nmin:nmax,this%DENSITY), &
         prim(:,:,nmin:nmax,this%XVELOCITY),cons(:,:,nmin:nmax,this%XMOMENTUM),&
         cons(:,:,nmin:nmax,this%YMOMENTUM),xfluxes(:,:,nmin:nmax,this%DENSITY),&
         xfluxes(:,:,nmin:nmax,this%XMOMENTUM),xfluxes(:,:,nmin:nmax,this%YMOMENTUM))
   END SUBROUTINE CalculateFluxesX_euler2Dit


  PURE SUBROUTINE CalculateFluxesY_euler2Dit(this,Mesh,nmin,nmax,prim,cons,yfluxes)
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
    CALL CalculateFlux_euler2Dit(this%csiso,prim(:,:,nmin:nmax,this%DENSITY), &
         prim(:,:,nmin:nmax,this%YVELOCITY),cons(:,:,nmin:nmax,this%YMOMENTUM), &
         cons(:,:,nmin:nmax,this%XMOMENTUM),yfluxes(:,:,nmin:nmax,this%DENSITY), &
         yfluxes(:,:,nmin:nmax,this%YMOMENTUM),yfluxes(:,:,nmin:nmax,this%XMOMENTUM))
  END SUBROUTINE CalculateFluxesY_euler2Dit


  PURE SUBROUTINE GeometricalSources_center(this,Mesh,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
         :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,pvar,cvar
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! geometrical source terms in momentum equationes
    sterm(:,:,this%DENSITY)   = 0.
    sterm(:,:,this%XMOMENTUM) = MomentumSourcesX_euler2Dit(this%csiso, &
         cvar(:,:,this%DENSITY),cvar(:,:,this%YMOMENTUM), &
         pvar(:,:,this%XVELOCITY),pvar(:,:,this%YVELOCITY), &
         Mesh%cxyx(:,:,1),Mesh%cyxy(:,:,1),Mesh%czxz(:,:,1))
    sterm(:,:,this%YMOMENTUM) = MomentumSourcesY_euler2Dit(this%csiso, &
         cvar(:,:,this%DENSITY),cvar(:,:,this%XMOMENTUM), &
         pvar(:,:,this%XVELOCITY),pvar(:,:,this%YVELOCITY), &
         Mesh%cxyx(:,:,1),Mesh%cyxy(:,:,1),Mesh%czyz(:,:,1))
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
    INTENT(IN)        :: this,Mesh,prim,cons
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! geometrical source terms in momentum equationes
    ! sum up all four corner values
    sterm(:,:,this%DENSITY)   = 0.
    sterm(:,:,this%XMOMENTUM) = SUM(MomentumSourcesX_euler2Dit(this%csiso, &
         cons(:,:,:,this%density),cons(:,:,:,this%YMOMENTUM), &
         prim(:,:,:,this%XVELOCITY),prim(:,:,:,this%YVELOCITY), &
         Mesh%cxyx(:,:,:),Mesh%cyxy(:,:,:),Mesh%czxz(:,:,:)),DIM=3)
    sterm(:,:,this%YMOMENTUM) = SUM(MomentumSourcesY_euler2Dit(this%csiso, &
         cons(:,:,:,this%DENSITY),cons(:,:,:,this%XMOMENTUM), &
         prim(:,:,:,this%XVELOCITY),prim(:,:,:,this%YVELOCITY), &
         Mesh%cxyx(:,:,:),Mesh%cyxy(:,:,:),Mesh%czyz(:,:,:)),DIM=3)
  END SUBROUTINE GeometricalSources_faces


  ! momentum and energy sources due to external force
  PURE SUBROUTINE ExternalSources_euler2Dit(this,Mesh,accel,pvar,cvar,sterm)
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
  END SUBROUTINE ExternalSources_euler2Dit


  PURE SUBROUTINE ViscositySources_euler2Dit(this,Mesh,Sources,pvar,cvar,sterm)
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
                ( Mesh%fhy(i,j,2) * ( pvar(i+1,j,this%XVELOCITY)+pvar(i,j,this%XVELOCITY) ) &
                - Mesh%fhy(i,j,1) * ( pvar(i,j,this%XVELOCITY)+pvar(i-1,j,this%XVELOCITY) ) ) &
                * Mesh%dydV(i,j) &
                +(Mesh%fhx(i,j,4) * ( pvar(i,j+1,this%YVELOCITY)+pvar(i,j,this%YVELOCITY) ) &
                - Mesh%fhx(i,j,3) * ( pvar(i,j,this%YVELOCITY)+pvar(i,j-1,this%YVELOCITY) ) ) &
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
               + Sources%btxy(i,j)
               
          Sources%btyy(i,j) = Sources%dynvis(i,j) * &
               ( (pvar(i,j+1,this%YVELOCITY) - pvar(i,j-1,this%YVELOCITY)) / Mesh%dly(i,j) &
               + 2.0 * Mesh%cyxy(i,j,1) * pvar(i,j,this%XVELOCITY) ) &
               + Sources%btxy(i,j)

          Sources%btxy(i,j) = Sources%dynvis(i,j) * ( 0.5 * &
               ( (pvar(i+1,j,this%YVELOCITY) - pvar(i-1,j,this%YVELOCITY)) / Mesh%dlx(i,j) &
               + (pvar(i,j+1,this%XVELOCITY) - pvar(i,j-1,this%XVELOCITY)) / Mesh%dly(i,j) ) &
               - Mesh%cxyx(i,j,1) * pvar(i,j,this%XVELOCITY) &
               - Mesh%cyxy(i,j,1) * pvar(i,j,this%YVELOCITY) )
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

          Sources%ftxy(i,j,1) = 0.5 * ( Sources%btxy(i-1,j) + Sources%btxy(i,j) )
          Sources%ftxy(i,j,2) = 0.5 * ( Sources%btxy(i+1,j) + Sources%btxy(i,j) )
          Sources%ftxy(i,j,3) = 0.5 * ( Sources%btxy(i,j-1) + Sources%btxy(i,j) )
          Sources%ftxy(i,j,4) = 0.5 * ( Sources%btxy(i,j+1) + Sources%btxy(i,j) )
       END DO
    END DO

    ! viscosity source terms
    sterm(:,:,this%DENSITY) = 0.0 

    ! (a) momentum sources
    sterm(:,:,this%XMOMENTUM) = Mesh%dydV(:,:) * &
         ( Mesh%fhy(:,:,2) * Sources%ftxx(:,:,2) - Mesh%fhy(:,:,1) * Sources%ftxx(:,:,1) &
         - Mesh%bhz(:,:) * (Mesh%fhy(:,:,2) - Mesh%fhy(:,:,1)) * Sources%btyy(:,:) ) &
         + Mesh%dxdV(:,:) * &
         ( Mesh%fhx(:,:,4) * Sources%ftxy(:,:,4) - Mesh%fhx(:,:,3) * Sources%ftxy(:,:,3) &
         + Mesh%bhz(:,:) * (Mesh%fhx(:,:,4) - Mesh%fhx(:,:,3)) * Sources%btxy(:,:) )

    sterm(:,:,this%YMOMENTUM) = Mesh%dydV(:,:) * &
         ( Mesh%fhy(:,:,2) * Sources%ftxy(:,:,2) - Mesh%fhy(:,:,1) * Sources%ftxy(:,:,1) &
         + Mesh%bhz(:,:) * (Mesh%fhy(:,:,2) - Mesh%fhy(:,:,1)) * Sources%btxy(:,:) ) &
         + Mesh%dxdV(:,:) * &
         ( Mesh%fhx(:,:,4) * Sources%ftyy(:,:,4) - Mesh%fhx(:,:,3) * Sources%ftyy(:,:,3) &
         - Mesh%bhz(:,:) * (Mesh%fhx(:,:,4) - Mesh%fhx(:,:,3)) * Sources%btxx(:,:) )
  END SUBROUTINE ViscositySources_euler2Dit


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
    CALL Cons2Prim_euler2Dit(cvar(:,:,this%DENSITY),cvar(:,:,this%XMOMENTUM), &
         cvar(:,:,this%YMOMENTUM),pvar(:,:,this%DENSITY), &
         pvar(:,:,this%XVELOCITY),pvar(:,:,this%YVELOCITY))
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
    CALL Cons2Prim_euler2Dit(cons(:,:,:,this%DENSITY),cons(:,:,:,this%XMOMENTUM), &
         cons(:,:,:,this%YMOMENTUM),prim(:,:,:,this%DENSITY), &
         prim(:,:,:,this%XVELOCITY),prim(:,:,:,this%YVELOCITY))
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
    CALL Prim2Cons_euler2Dit(pvar(:,:,this%DENSITY),pvar(:,:,this%XVELOCITY), &
         pvar(:,:,this%YVELOCITY),cvar(:,:,this%DENSITY), &
         cvar(:,:,this%XMOMENTUM),cvar(:,:,this%YMOMENTUM))
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
    CALL Prim2Cons_euler2Dit(prim(:,:,:,this%DENSITY),prim(:,:,:,this%XVELOCITY), &
         prim(:,:,:,this%YVELOCITY),cons(:,:,:,this%DENSITY), &
         cons(:,:,:,this%XMOMENTUM),cons(:,:,:,this%YMOMENTUM))
  END SUBROUTINE Convert2Conservative_faces


  ELEMENTAL SUBROUTINE CalculateFlux_euler2Dit(cs,rho,v,m1,m2,f1,f2,f3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,rho,v,m1,m2
    REAL, INTENT(OUT) :: f1, f2, f3
    !------------------------------------------------------------------------!
    f1 = rho*v
    f2 = m1*v + rho*cs*cs
    f3 = m2*v
  END SUBROUTINE CalculateFlux_euler2Dit


  ! momentum source terms due to inertial forces
  ELEMENTAL FUNCTION MomentumSourcesX_euler2Dit(cs,rho,my,vx,vy,cxyx,cyxy,czxz) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,rho,my,vx,vy,cxyx,cyxy,czxz
    REAL :: st
    !------------------------------------------------------------------------!
    st = -my * (cxyx * vx - cyxy * vy) + (cyxy + czxz) * rho * cs*cs
  END FUNCTION MomentumSourcesX_euler2Dit


  ELEMENTAL FUNCTION MomentumSourcesY_euler2Dit(cs,rho,mx,vx,vy,cxyx,cyxy,czyz) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,rho,mx,vx,vy,cxyx,cyxy,czyz
    REAL :: st
    !------------------------------------------------------------------------!
    st = mx * (cxyx * vx - cyxy * vy) + (cxyx + czyz) * rho * cs*cs
  END FUNCTION MomentumSourcesY_euler2Dit


  ELEMENTAL SUBROUTINE Cons2Prim_euler2Dit(rho_in,mu,mv,rho_out,u,v)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho_in,mu,mv
    REAL, INTENT(OUT) :: rho_out,u,v
    !------------------------------------------------------------------------!
    REAL :: inv_rho
    !------------------------------------------------------------------------!
    inv_rho = 1./rho_in
    rho_out = rho_in
    u = mu * inv_rho
    v = mv * inv_rho
  END SUBROUTINE Cons2Prim_euler2Dit

  
  ELEMENTAL SUBROUTINE Prim2Cons_euler2Dit(rho_in,u,v,rho_out,mu,mv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho_in,u,v
    REAL, INTENT(OUT) :: rho_out,mu,mv
    !------------------------------------------------------------------------!
    rho_out = rho_in
    mu = rho_in * u
    mv = rho_in * v
  END SUBROUTINE Prim2Cons_euler2Dit
  

  SUBROUTINE ClosePhysics_euler2Dit(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL ClosePhysics(this)
  END SUBROUTINE ClosePhysics_euler2Dit


END MODULE physics_euler2Disothm
