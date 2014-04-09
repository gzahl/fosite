!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler2Disothm.f90                                         #
!#                                                                           #
!# Copyright (C) 2007 - 2010                                                 #
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
! basic module for 2D isothermal Euler equations
!----------------------------------------------------------------------------!
MODULE physics_euler2Disothm
  USE physics_common
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
       CalculateStresses_euler2Dit, &
       GeometricalSources_euler2Dit, &
       ViscositySources_euler2Dit, &
       ExternalSources_euler2Dit, &
       Convert2Primitive_euler2Dit, &
       Convert2Conservative_euler2Dit, &
       ReflectionMasks_euler2Dit, &
       AxisMasks_euler2Dit, &
       SetWaveSpeeds_euler2Dit, &
       MomentumSourcesX_euler2Dit, &
       MomentumSourcesY_euler2Dit, &
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

  SUBROUTINE InitPhysics_euler2Dit(this,problem)
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
    this%XVELOCITY = 2                                 ! x-velocity          !
    this%XMOMENTUM = 2                                 ! x-momentum          !
    this%YVELOCITY = num_var                           ! y-velocity          !
    this%YMOMENTUM = num_var                           ! y-momentum          !
    this%ZVELOCITY = 0                                 ! no z-velocity       !
    this%ZMOMENTUM = 0                                 ! no z-momentum       !
    this%PRESSURE  = 0                                 ! no pressure         !
    this%ENERGY    = 0                                 ! no total energy     !
    ! set names for primitive and conservative variables
    this%pvarname(this%DENSITY)   = "density"
    this%pvarname(this%XVELOCITY) = "x-velocity"
    this%pvarname(this%YVELOCITY) = "y-velocity"
    this%cvarname(this%DENSITY)   = "density"
    this%cvarname(this%XMOMENTUM) = "x-momentum"
    this%cvarname(this%YMOMENTUM) = "y-momentum"

    ALLOCATE(this%structure(3),this%errormap(0:1),STAT = err)
    ! abort if allocation fails
    IF (err.NE.0) &
         CALL Error(this, "InitPhysics_euler2Dit", "Unable to allocate memory.")
    this%nstruc = 3
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
    this%structure(3)%dim  = 2

    ! set errormapping (CheckData Symbols)
    this%errormap(0) = 'X'
    this%errormap(1) = 'D'
  END SUBROUTINE InitPhysics_euler2Dit


  PURE FUNCTION CheckData_euler2Dit(this,Mesh,pvar,pold,mr) RESULT (bad_data)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    ! checks for bad density or pressure on the whole grid including ghost cells
    ! return values:
    !   0 : valid density and pressure data
    !   1 : density < this%rhomin, i.e. vacuum generated
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
         :: pvar,pold
    INTEGER,DIMENSION(4) :: mr
    INTEGER           :: bad_data
    !------------------------------------------------------------------------!
    REAL              :: rhomin
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,pvar,pold,mr
    !------------------------------------------------------------------------!
    ! density minimum
    rhomin = MINVAL(pvar(mr(1):mr(2),mr(3):mr(4),this%DENSITY))
    IF (rhomin.LE.this%rhomin) THEN
       bad_data = 1
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
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! compute minimal and maximal wave speeds at cell centers
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! x-direction
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2Dit(this%csiso,pvar(i,j,this%XVELOCITY),&
               this%amin(i,j),this%amax(i,j))
          ! y-direction
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2Dit(this%csiso,pvar(i,j,this%YVELOCITY),&
               this%bmin(i,j),this%bmax(i,j))
       END DO
    END DO
  END SUBROUTINE CalculateWaveSpeeds_center


  PURE SUBROUTINE CalculateWaveSpeeds_faces(this,Mesh,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%VNUM) &
         :: prim
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,prim
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! compute minimal and maximal wave speeds at cell interfaces
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! western
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2Dit(this%csiso,prim(i,j,1,this%XVELOCITY), &
               this%tmin(i,j,1),this%tmax(i,j,1))
          ! eastern
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2Dit(this%csiso,prim(i,j,2,this%XVELOCITY), &
               this%amin(i,j),this%amax(i,j))
          ! southern
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2Dit(this%csiso,prim(i,j,3,this%YVELOCITY), &
               this%tmin(i,j,2),this%tmax(i,j,2))
          ! northern
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2Dit(this%csiso,prim(i,j,4,this%YVELOCITY), &
               this%bmin(i,j),this%bmax(i,j))
       END DO
    END DO
    ! set minimal and maximal wave speeds at cell interfaces of neighboring cells
    ! western interfaces
    this%amin(Mesh%IMIN-1:Mesh%IMAX,:) = MIN(this%tmin(Mesh%IMIN:Mesh%IMAX+1,:,1), &
         this%amin(Mesh%IMIN-1:Mesh%IMAX,:))
    ! eastern interfaces
    this%amax(Mesh%IMIN-1:Mesh%IMAX,:) = MAX(this%tmax(Mesh%IMIN:Mesh%IMAX+1,:,1), &
         this%amax(Mesh%IMIN-1:Mesh%IMAX,:))
    ! southern interfaces
    this%bmin(:,Mesh%JMIN-1:Mesh%JMAX) = MIN(this%tmin(:,Mesh%JMIN:Mesh%JMAX+1,2), &
         this%bmin(:,Mesh%JMIN-1:Mesh%JMAX))
    ! northern interfaces
    this%bmax(:,Mesh%JMIN-1:Mesh%JMAX) = MAX(this%tmax(:,Mesh%JMIN:Mesh%JMAX+1,2), &
         this%bmax(:,Mesh%JMIN-1:Mesh%JMAX))
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


  PURE SUBROUTINE CalculateStresses_euler2Dit(this,Mesh,pvar,dynvis,bulkvis, &
       btxx,btxy,btyy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
         dynvis,bulkvis,btxx,btxy,btyy
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    REAL              :: bt_bulk
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,pvar,dynvis,bulkvis
    INTENT(INOUT)     :: btxx,btxy,btyy
    !------------------------------------------------------------------------!
    ! compute components of the stress tensor at cell bary centers
    ! inside the computational domain including one slice of ghost cells
!CDIR UNROLL=8
    DO j=Mesh%JMIN-1,Mesh%JMAX+1
!CDIR NODEP
       DO i=Mesh%IMIN-1,Mesh%IMAX+1
          ! compute bulk viscosity first contribution
          bt_bulk = bulkvis(i,j) * 0.5 * (&
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
          
          ! compute the diagonal elements of the stress tensor
          btxx(i,j) = dynvis(i,j) * &
               ( (pvar(i+1,j,this%XVELOCITY) - pvar(i-1,j,this%XVELOCITY)) / Mesh%dlx(i,j) &
               + 2.0 * Mesh%cxyx(i,j,1) * pvar(i,j,this%YVELOCITY) ) &
               + bt_bulk ! bulk viscosity contribution
               
          btyy(i,j) = dynvis(i,j) * &
               ( (pvar(i,j+1,this%YVELOCITY) - pvar(i,j-1,this%YVELOCITY)) / Mesh%dly(i,j) &
               + 2.0 * Mesh%cyxy(i,j,1) * pvar(i,j,this%XVELOCITY) ) &
               + bt_bulk ! bulk viscosity contribution

          ! compute the off-diagonal elements (no bulk viscosity)
          btxy(i,j) = dynvis(i,j) * ( 0.5 * &
               ( (pvar(i+1,j,this%YVELOCITY) - pvar(i-1,j,this%YVELOCITY)) / Mesh%dlx(i,j) &
               + (pvar(i,j+1,this%XVELOCITY) - pvar(i,j-1,this%XVELOCITY)) / Mesh%dly(i,j) ) &
               - Mesh%cxyx(i,j,1) * pvar(i,j,this%XVELOCITY) &
               - Mesh%cyxy(i,j,1) * pvar(i,j,this%YVELOCITY) )
       END DO
    END DO
  END SUBROUTINE CalculateStresses_euler2Dit


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
    sterm(:,:,this%XMOMENTUM) = MomentumSourcesX_euler2Dit(&
         cvar(:,:,this%YMOMENTUM),pvar(:,:,this%XVELOCITY),&
         pvar(:,:,this%YVELOCITY),pvar(:,:,this%DENSITY)*this%csiso**2, &
         Mesh%cxyx(:,:,1),Mesh%cyxy(:,:,1),Mesh%czxz(:,:,1))
    sterm(:,:,this%YMOMENTUM) = MomentumSourcesY_euler2Dit(&
         cvar(:,:,this%XMOMENTUM),pvar(:,:,this%XVELOCITY),&
         pvar(:,:,this%YVELOCITY),pvar(:,:,this%DENSITY)*this%csiso**2, &
         Mesh%cxyx(:,:,1),Mesh%cyxy(:,:,1),Mesh%czxz(:,:,1))
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
    sterm(:,:,this%XMOMENTUM) = SUM(MomentumSourcesX_euler2Dit(&
         cons(:,:,:,this%YMOMENTUM),prim(:,:,:,this%XVELOCITY),&
         prim(:,:,:,this%YVELOCITY),prim(:,:,:,this%DENSITY)*this%csiso**2, &
         Mesh%cxyx(:,:,:),Mesh%cyxy(:,:,:),Mesh%czxz(:,:,:)),DIM=3)
    sterm(:,:,this%YMOMENTUM) = SUM(MomentumSourcesY_euler2Dit(&
         cons(:,:,:,this%XMOMENTUM),prim(:,:,:,this%XVELOCITY),&
         prim(:,:,:,this%YVELOCITY),prim(:,:,:,this%DENSITY)*this%csiso**2, &
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


  PURE SUBROUTINE ViscositySources_euler2Dit(this,Mesh,pvar,btxx,btxy,btyy, &
       ftxx,ftxy,ftyy,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: &
         pvar,sterm
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
         btxx,btxy,btyy
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4) :: &
         ftxx,ftxy,ftyy
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,pvar,btxx,btxy,btyy
    INTENT(INOUT)     :: ftxx,ftxy,ftyy
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! mean values of stress tensor components across the cell interfaces
!CDIR UNROLL=8
    DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
       DO i=Mesh%IMIN,Mesh%IMAX
          ftxx(i,j,1) = 0.5 * ( btxx(i-1,j) + btxx(i,j) )
          ftxx(i,j,2) = 0.5 * ( btxx(i+1,j) + btxx(i,j) )
          ftxx(i,j,3) = 0.5 * ( btxx(i,j-1) + btxx(i,j) )
          ftxx(i,j,4) = 0.5 * ( btxx(i,j+1) + btxx(i,j) )

          ftyy(i,j,1) = 0.5 * ( btyy(i-1,j) + btyy(i,j) )
          ftyy(i,j,2) = 0.5 * ( btyy(i+1,j) + btyy(i,j) )
          ftyy(i,j,3) = 0.5 * ( btyy(i,j-1) + btyy(i,j) )
          ftyy(i,j,4) = 0.5 * ( btyy(i,j+1) + btyy(i,j) )

          ftxy(i,j,1) = 0.5 * ( btxy(i-1,j) + btxy(i,j) )
          ftxy(i,j,2) = 0.5 * ( btxy(i+1,j) + btxy(i,j) )
          ftxy(i,j,3) = 0.5 * ( btxy(i,j-1) + btxy(i,j) )
          ftxy(i,j,4) = 0.5 * ( btxy(i,j+1) + btxy(i,j) )

          ! viscosity source terms
          sterm(i,j,this%DENSITY) = 0.0 

          ! (a) momentum sources
          sterm(i,j,this%XMOMENTUM) = Mesh%dydV(i,j) * &
               ( Mesh%fhy(i,j,2) * ftxx(i,j,2) - Mesh%fhy(i,j,1) * ftxx(i,j,1) &
               - Mesh%bhz(i,j) * (Mesh%fhy(i,j,2) - Mesh%fhy(i,j,1)) * btyy(i,j) ) &
               + Mesh%dxdV(i,j) * &
               ( Mesh%fhx(i,j,4) * ftxy(i,j,4) - Mesh%fhx(i,j,3) * ftxy(i,j,3) &
               + Mesh%bhz(i,j) * (Mesh%fhx(i,j,4) - Mesh%fhx(i,j,3)) * btxy(i,j) )

          sterm(i,j,this%YMOMENTUM) = Mesh%dydV(i,j) * &
               ( Mesh%fhy(i,j,2) * ftxy(i,j,2) - Mesh%fhy(i,j,1) * ftxy(i,j,1) &
               + Mesh%bhz(i,j) * (Mesh%fhy(i,j,2) - Mesh%fhy(i,j,1)) * btxy(i,j) ) &
               + Mesh%dxdV(i,j) * &
               ( Mesh%fhx(i,j,4) * ftyy(i,j,4) - Mesh%fhx(i,j,3) * ftyy(i,j,3) &
               - Mesh%bhz(i,j) * (Mesh%fhx(i,j,4) - Mesh%fhx(i,j,3)) * btxx(i,j) )
       END DO
    END DO

  END SUBROUTINE ViscositySources_euler2Dit


  PURE SUBROUTINE Convert2Primitive_center(this,Mesh,i1,i2,j1,j2,cvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    INTEGER              :: i1,i2,j1,j2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
                         :: cvar,pvar
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,i1,i2,j1,j2,cvar
    INTENT(OUT) :: pvar
    !------------------------------------------------------------------------!
    CALL Cons2Prim_euler2Dit(cvar(i1:i2,j1:j2,this%DENSITY), &
         cvar(i1:i2,j1:j2,this%XMOMENTUM),cvar(i1:i2,j1:j2,this%YMOMENTUM), &
         pvar(i1:i2,j1:j2,this%DENSITY),pvar(i1:i2,j1:j2,this%XVELOCITY), &
         pvar(i1:i2,j1:j2,this%YVELOCITY))
  END SUBROUTINE Convert2Primitive_center


  PURE SUBROUTINE Convert2Primitive_faces(this,Mesh,i1,i2,j1,j2,cons,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    INTEGER              :: i1,i2,j1,j2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%VNUM) &
                         :: cons,prim
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,i1,i2,j1,j2,cons
    INTENT(OUT) :: prim
    !------------------------------------------------------------------------!
    CALL Cons2Prim_euler2Dit(cons(i1:i2,j1:j2,:,this%DENSITY), &
         cons(i1:i2,j1:j2,:,this%XMOMENTUM),cons(i1:i2,j1:j2,:,this%YMOMENTUM), &
         prim(i1:i2,j1:j2,:,this%DENSITY),prim(i1:i2,j1:j2,:,this%XVELOCITY), &
         prim(i1:i2,j1:j2,:,this%YVELOCITY))
  END SUBROUTINE Convert2Primitive_faces


  PURE SUBROUTINE Convert2Conservative_center(this,Mesh,i1,i2,j1,j2,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    INTEGER              :: i1,i2,j1,j2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
                         :: cvar,pvar
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,i1,i2,j1,j2,pvar
    INTENT(OUT) :: cvar
    !------------------------------------------------------------------------!
    CALL Prim2Cons_euler2Dit(pvar(i1:i2,j1:j2,this%DENSITY), &
         pvar(i1:i2,j1:j2,this%XVELOCITY), pvar(i1:i2,j1:j2,this%YVELOCITY), &
         cvar(i1:i2,j1:j2,this%DENSITY),cvar(i1:i2,j1:j2,this%XMOMENTUM), &
         cvar(i1:i2,j1:j2,this%YMOMENTUM))
  END SUBROUTINE Convert2Conservative_center


  PURE SUBROUTINE Convert2Conservative_faces(this,Mesh,i1,i2,j1,j2,prim,cons)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)    :: this
    TYPE(Mesh_TYP)       :: Mesh
    INTEGER              :: i1,i2,j1,j2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%VNUM) &
                         :: cons,prim
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,i1,i2,j1,j2,prim
    INTENT(OUT) :: cons
    !------------------------------------------------------------------------!
    CALL Prim2Cons_euler2Dit(prim(i1:i2,j1:j2,:,this%DENSITY), &
         prim(i1:i2,j1:j2,:,this%XVELOCITY),prim(i1:i2,j1:j2,:,this%YVELOCITY), &
         cons(i1:i2,j1:j2,:,this%DENSITY),cons(i1:i2,j1:j2,:,this%XMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%YMOMENTUM))
  END SUBROUTINE Convert2Conservative_faces


  PURE SUBROUTINE ReflectionMasks_euler2Dit(this,reflX,reflY)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    LOGICAL, DIMENSION(this%VNUM) :: reflX,reflY
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this
    INTENT(OUT)       :: reflX,reflY
    !------------------------------------------------------------------------!
    ! western / eastern boundary
    reflX(this%DENSITY)   = .FALSE.
    reflX(this%XVELOCITY) = .TRUE.
    reflX(this%YVELOCITY) = .FALSE.
    ! southern / northern boundary
    reflY(this%DENSITY)   = .FALSE.
    reflY(this%XVELOCITY) = .FALSE.
    reflY(this%YVELOCITY) = .TRUE.
  END SUBROUTINE ReflectionMasks_euler2Dit


  PURE SUBROUTINE AxisMasks_euler2Dit(this,reflX,reflY)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    LOGICAL, DIMENSION(this%VNUM) :: reflX,reflY
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this
    INTENT(OUT)       :: reflX,reflY
    !------------------------------------------------------------------------!
    ! western / eastern boundary
    reflX(this%DENSITY)   = .FALSE.
    reflX(this%XVELOCITY) = .TRUE.
    reflX(this%YVELOCITY) = .TRUE.
    ! southern / northern boundary
    reflY(this%DENSITY)   = .FALSE.
    reflY(this%XVELOCITY) = .TRUE.
    reflY(this%YVELOCITY) = .TRUE.
  END SUBROUTINE AxisMasks_euler2Dit


  ELEMENTAL SUBROUTINE SetWaveSpeeds_euler2Dit(cs,v,amin,amax)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,v
    REAL, INTENT(OUT) :: amin,amax
    !------------------------------------------------------------------------!
    ! minimal and maximal wave speeds
    amin = MIN(0.,v-cs)
    amax = MAX(0.,v+cs)
  END SUBROUTINE SetWaveSpeeds_euler2Dit


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
  ! P is the isothermal pressure rho*cs*cs
  ELEMENTAL FUNCTION MomentumSourcesX_euler2Dit(my,vx,vy,P,cxyx,cyxy,czxz) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: my,vx,vy,P,cxyx,cyxy,czxz
    REAL :: st
    !------------------------------------------------------------------------!
    st = -my * (cxyx * vx - cyxy * vy) + (cyxy + czxz) * P
  END FUNCTION MomentumSourcesX_euler2Dit

  ELEMENTAL FUNCTION MomentumSourcesY_euler2Dit(mx,vx,vy,P,cxyx,cyxy,czyz) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: mx,vx,vy,P,cxyx,cyxy,czyz
    REAL :: st
    !------------------------------------------------------------------------!
    st = mx * (cxyx * vx - cyxy * vy) + (cxyx + czyz) * P
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
    DEALLOCATE(this%structure,this%errormap)
    CALL ClosePhysics(this)
  END SUBROUTINE ClosePhysics_euler2Dit


END MODULE physics_euler2Disothm
