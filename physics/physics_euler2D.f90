!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler2D.f90                                               #
!#                                                                           #
!# Copyright (C) 2007-2012                                                   #
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
! basic module for 2D Euler equations
!----------------------------------------------------------------------------!
MODULE physics_euler2D
  USE physics_common
  USE mesh_common, ONLY : Mesh_TYP
  USE sources_common, ONLY : Sources_TYP
  USE physics_euler2Disothm, &
       CalculateStresses_euler2D => CalculateStresses_euler2Dit, &
       SetWaveSpeeds_euler2D => SetWaveSpeeds_euler2Dit, &
       MomentumSourcesX_euler2D => MomentumSourcesX_euler2Dit, &
       MomentumSourcesY_euler2D => MomentumSourcesY_euler2Dit
  USE mesh_generic
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
  INTEGER, PARAMETER :: num_var = 4              ! number of variables       !
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 2D"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Physics_TYP, &
       ! methods
       InitPhysics_euler2D, &
       ClosePhysics_euler2D, &
       CheckData_euler2D, &
       CalculateWaveSpeeds_euler2D, &
       CalculateFluxesX_euler2D, &
       CalculateFluxesY_euler2D, &
       CalculateStresses_euler2D, &
       GeometricalSources_euler2D, &
       ViscositySources_euler2D, &
       ExternalSources_euler2D, &
       Convert2Primitive_euler2D, &
       Convert2Conservative_euler2D, &
       ReflectionMasks_euler2D, &
       AxisMasks_euler2D, &
       GetSoundSpeed_euler2D, &
       SetWaveSpeeds_euler2D, &
       MomentumSourcesX_euler2D, &
       MomentumSourcesY_euler2D
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitPhysics_euler2D(this,Mesh,problem)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: problem
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,problem
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
    this%ZVELOCITY = 0                                 ! no z-velocity       !
    this%ZMOMENTUM = 0                                 ! no z-momentum       !
    ! set names for primitive and conservative variables
    this%pvarname(this%DENSITY)   = "density"
    this%pvarname(this%XVELOCITY) = "x-velocity"
    this%pvarname(this%YVELOCITY) = "y-velocity"
    this%pvarname(this%PRESSURE)  = "pressure"
    this%cvarname(this%DENSITY)   = "density"
    this%cvarname(this%XMOMENTUM) = "x-momentum"
    this%cvarname(this%YMOMENTUM) = "y-momentum"
    this%cvarname(this%ENERGY)    = "energy"
    
    ! allocate memory for arrays used in Euler2D
    ALLOCATE(this%csound(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4), &
         this%structure(4),this%errormap(0:3),STAT = err)
    ! abort if allocation fails
    IF (err.NE.0) &
      CALL Error(this,"InitPhysics_euler2D","Unable to allocate memory.")

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
    this%structure(3)%dim  = 2
    this%structure(4)%name = "pressure"
    this%structure(4)%pos  = this%PRESSURE
    this%structure(4)%rank = 0
    this%structure(4)%dim  = 1

    !set errormapping (CheckData Symbols)
    this%errormap(0) = 'X'
    this%errormap(1) = 'D'
    this%errormap(2) = 'P'
    this%errormap(3) = 'B'
  END SUBROUTINE InitPhysics_euler2D


  PURE FUNCTION CheckData_euler2D(this,Mesh,pvar,pold,mr) RESULT (bad_data)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    ! checks for bad density or pressure on the whole grid including ghost cells
    ! return values:
    !   0 : valid density and pressure data
    !   1 : density < this%rhomin, i.e. vacuum generated
    !   2 : pressure < this%pmin
    !   3 : (pold - pvar)/pold < dpmax pressure variations to large
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
         :: pvar,pold
    INTEGER, DIMENSION(4) :: mr
    INTEGER, DIMENSION(4) :: def_mr
    INTEGER           :: bad_data
    !------------------------------------------------------------------------!
    REAL              :: pmin,dpmax,rhomin
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,pvar,pold,mr
    !------------------------------------------------------------------------!
    ! density minimum
    rhomin = MINVAL(pvar(mr(1):mr(2),mr(3):mr(4),this%DENSITY))
    ! check for density underflow or irregular values, i.e. inf/nan
    IF ((rhomin.LE.this%rhomin).OR.(rhomin*0.0.NE.0.0)) THEN
       bad_data = 1
       RETURN
    END IF
    ! pressure minimum
    pmin  = MINVAL(pvar(mr(1):mr(2),mr(3):mr(4),this%PRESSURE))
    ! check for pressure underflow or irregular values, i.e. inf/nan
    IF ((pmin.LE.this%pmin).OR.(pmin*0.0.NE.0.0)) THEN
       bad_data = 2
       RETURN
    END IF
    ! pressure variations (check only within computational domain, otherwise
    ! farfield (and maybe other) boundary conditions will not work)
    !Meshrange: IMIN,IMAX,JMIN,JMAX (no ghost cells!)
    def_mr(1) = max(Mesh%IMIN,mr(1))
    def_mr(2) = min(Mesh%IMAX,mr(2))
    def_mr(3) = max(Mesh%JMIN,mr(3))
    def_mr(4) = min(Mesh%JMAX,mr(4))
    dpmax = MAXVAL(ABS(1.0-pvar(def_mr(1):def_mr(2),def_mr(3):def_mr(4),this%PRESSURE) &
         / pold(def_mr(1):def_mr(2),def_mr(3):def_mr(4),this%PRESSURE)))
    IF (dpmax.GT.this%dpmax) THEN
       bad_data = 3
       RETURN
    END IF
    ! everything ok
    bad_data = 0
  END FUNCTION CheckData_euler2D


  PURE SUBROUTINE CalculateWaveSpeeds_center(this,Mesh,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
         :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    REAL              :: cs
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! compute minimal and maximal wave speeds at cell centers
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! sound speed
!CDIR IEXPAND
          cs = GetSoundSpeed_euler2D(this%gamma, &
               pvar(i,j,this%DENSITY),pvar(i,j,this%PRESSURE))
          ! x-direction
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2D(cs,pvar(i,j,this%XVELOCITY),&
               this%amin(i,j),this%amax(i,j))
          ! y-direction
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2D(cs,pvar(i,j,this%YVELOCITY),&
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
          CALL SetWaveSpeeds_euler2D(GetSoundSpeed_euler2D(this%gamma, &
               prim(i,j,1,this%DENSITY),prim(i,j,1,this%PRESSURE)), &
               prim(i,j,1,this%XVELOCITY), &
               this%tmin(i,j,1),this%tmax(i,j,1))
          ! eastern
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2D(GetSoundSpeed_euler2D(this%gamma, &
               prim(i,j,2,this%DENSITY),prim(i,j,2,this%PRESSURE)), &
               prim(i,j,2,this%XVELOCITY), &
               this%amin(i,j),this%amax(i,j))
          ! southern
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2D(GetSoundSpeed_euler2D(this%gamma, &
               prim(i,j,3,this%DENSITY),prim(i,j,3,this%PRESSURE)), &
               prim(i,j,3,this%YVELOCITY), &
               this%tmin(i,j,2),this%tmax(i,j,2))
          ! northern
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2D(GetSoundSpeed_euler2D(this%gamma, &
               prim(i,j,4,this%DENSITY),prim(i,j,4,this%PRESSURE)), &
               prim(i,j,4,this%YVELOCITY), &
               this%bmin(i,j),this%bmax(i,j))
       END DO
    END DO
    ! set minimal and maximal wave speeds at cell interfaces of neighboring cells
    DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
       DO i=Mesh%IMIN-1,Mesh%IMAX
          ! western interfaces
          this%amin(i,j) = MIN(this%tmin(i+1,j,1),this%amin(i,j))
          ! eastern interfaces
          this%amax(i,j) = MAX(this%tmax(i+1,j,1),this%amax(i,j))
       END DO
    END DO
!CDIR COLLAPSE
    DO j=Mesh%JMIN-1,Mesh%JMAX
!CDIR NODEP
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! southern interfaces
          this%bmin(i,j) = MIN(this%tmin(i,j+1,2),this%bmin(i,j))
          ! northern interfaces
          this%bmax(i,j) = MAX(this%tmax(i,j+1,2),this%bmax(i,j))
       END DO
    END DO
  END SUBROUTINE CalculateWaveSpeeds_faces


  PURE SUBROUTINE CalculateFluxesX_euler2D(this,Mesh,nmin,nmax,prim,cons,xfluxes)
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
    CALL CalculateFlux_euler2D(prim(:,:,nmin:nmax,this%DENSITY), &
         prim(:,:,nmin:nmax,this%XVELOCITY),prim(:,:,nmin:nmax,this%PRESSURE), &
         cons(:,:,nmin:nmax,this%XMOMENTUM),cons(:,:,nmin:nmax,this%YMOMENTUM), &
         cons(:,:,nmin:nmax,this%ENERGY),xfluxes(:,:,nmin:nmax,this%DENSITY),&
         xfluxes(:,:,nmin:nmax,this%XMOMENTUM),xfluxes(:,:,nmin:nmax,this%YMOMENTUM), &
         xfluxes(:,:,nmin:nmax,this%ENERGY))
  END SUBROUTINE CalculateFluxesX_euler2D


  PURE SUBROUTINE CalculateFluxesY_euler2D(this,Mesh,nmin,nmax,prim,cons,yfluxes)
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
    CALL CalculateFlux_euler2D(prim(:,:,nmin:nmax,this%DENSITY), &
         prim(:,:,nmin:nmax,this%YVELOCITY),prim(:,:,nmin:nmax,this%PRESSURE), &
         cons(:,:,nmin:nmax,this%YMOMENTUM),cons(:,:,nmin:nmax,this%XMOMENTUM), &
         cons(:,:,nmin:nmax,this%ENERGY),yfluxes(:,:,nmin:nmax,this%DENSITY), &
         yfluxes(:,:,nmin:nmax,this%YMOMENTUM),yfluxes(:,:,nmin:nmax,this%XMOMENTUM), &
         yfluxes(:,:,nmin:nmax,this%ENERGY))
  END SUBROUTINE CalculateFluxesY_euler2D


  PURE SUBROUTINE GeometricalSources_center(this,Mesh,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
         :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,pvar,cvar
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! no geometrical density or energy sources
          sterm(i,j,this%DENSITY) = 0.
          sterm(i,j,this%ENERGY) = 0.
          ! geometrical source terms in momentum equationes
          sterm(i,j,this%XMOMENTUM) = MomentumSourcesX_euler2D(&
              cvar(i,j,this%YMOMENTUM),pvar(i,j,this%XVELOCITY),&
              pvar(i,j,this%YVELOCITY),pvar(i,j,this%PRESSURE), &
              Mesh%cxyx(i,j,1),Mesh%cyxy(i,j,1),Mesh%czxz(i,j,1))
          sterm(i,j,this%YMOMENTUM) = MomentumSourcesY_euler2D(&
              cvar(i,j,this%XMOMENTUM),pvar(i,j,this%XVELOCITY),&
              pvar(i,j,this%YVELOCITY),pvar(i,j,this%PRESSURE), &
              Mesh%cxyx(i,j,1),Mesh%cyxy(i,j,1),Mesh%czxz(i,j,1))
       END DO
    END DO
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
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,prim,cons
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! no geometrical density or energy sources
          sterm(i,j,this%DENSITY)   = 0.
          sterm(i,j,this%ENERGY)    = 0.
          ! momentum sources (sum up corner values, don't use SUM function,
          ! because it prevents COLLAPSING and causes poor vectorization
          sterm(i,j,this%XMOMENTUM) = MomentumSourcesX_euler2D(&
              cons(i,j,1,this%YMOMENTUM),prim(i,j,1,this%XVELOCITY),&
              prim(i,j,1,this%YVELOCITY),prim(i,j,1,this%PRESSURE), &
              Mesh%cxyx(i,j,1),Mesh%cyxy(i,j,1),Mesh%czxz(i,j,1)) &
            + MomentumSourcesX_euler2D(&
              cons(i,j,2,this%YMOMENTUM),prim(i,j,2,this%XVELOCITY),&
              prim(i,j,2,this%YVELOCITY),prim(i,j,2,this%PRESSURE), &
              Mesh%cxyx(i,j,2),Mesh%cyxy(i,j,2),Mesh%czxz(i,j,2)) &
            + MomentumSourcesX_euler2D(&
              cons(i,j,3,this%YMOMENTUM),prim(i,j,3,this%XVELOCITY),&
              prim(i,j,3,this%YVELOCITY),prim(i,j,3,this%PRESSURE), &
              Mesh%cxyx(i,j,3),Mesh%cyxy(i,j,3),Mesh%czxz(i,j,3)) &
            + MomentumSourcesX_euler2D(&
              cons(i,j,4,this%YMOMENTUM),prim(i,j,4,this%XVELOCITY),&
              prim(i,j,4,this%YVELOCITY),prim(i,j,4,this%PRESSURE), &
              Mesh%cxyx(i,j,4),Mesh%cyxy(i,j,4),Mesh%czxz(i,j,4)) 

          sterm(i,j,this%YMOMENTUM) = MomentumSourcesY_euler2D(&
              cons(i,j,1,this%XMOMENTUM),prim(i,j,1,this%XVELOCITY),&
              prim(i,j,1,this%YVELOCITY),prim(i,j,1,this%PRESSURE), &
              Mesh%cxyx(i,j,1),Mesh%cyxy(i,j,1),Mesh%czyz(i,j,1)) &
            + MomentumSourcesY_euler2D(&
              cons(i,j,2,this%XMOMENTUM),prim(i,j,2,this%XVELOCITY),&
              prim(i,j,2,this%YVELOCITY),prim(i,j,2,this%PRESSURE), &
              Mesh%cxyx(i,j,2),Mesh%cyxy(i,j,2),Mesh%czyz(i,j,2)) &
            + MomentumSourcesY_euler2D(&
              cons(i,j,3,this%XMOMENTUM),prim(i,j,3,this%XVELOCITY),&
              prim(i,j,3,this%YVELOCITY),prim(i,j,3,this%PRESSURE), &
              Mesh%cxyx(i,j,3),Mesh%cyxy(i,j,3),Mesh%czyz(i,j,3)) &
            + MomentumSourcesY_euler2D(&
              cons(i,j,4,this%XMOMENTUM),prim(i,j,4,this%XVELOCITY),&
              prim(i,j,4,this%YVELOCITY),prim(i,j,4,this%PRESSURE), &
              Mesh%cxyx(i,j,4),Mesh%cyxy(i,j,4),Mesh%czyz(i,j,4))
       END DO
    END DO 
  END SUBROUTINE GeometricalSources_faces


  ! momentum and energy sources due to external force
  PURE SUBROUTINE ExternalSources_euler2D(this,Mesh,accel,pvar,cvar,sterm)
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
         :: accel
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
         :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,accel,pvar,cvar
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          sterm(i,j,this%DENSITY)   = 0.
          sterm(i,j,this%XMOMENTUM) = pvar(i,j,this%DENSITY) * accel(i,j,1)
          sterm(i,j,this%YMOMENTUM) = pvar(i,j,this%DENSITY) * accel(i,j,2)
          sterm(i,j,this%ENERGY)    = cvar(i,j,this%XMOMENTUM) * accel(i,j,1) + &
               cvar(i,j,this%YMOMENTUM) * accel(i,j,2)
       END DO
    END DO
  END SUBROUTINE ExternalSources_euler2D


  PURE SUBROUTINE ViscositySources_euler2D(this,Mesh,pvar,btxx,btxy,btyy,sterm)         
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: &
         pvar,sterm
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
         btxx,btxy,btyy
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,btxx,btxy,btyy
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL ViscositySources_euler2Dit(this,Mesh,pvar,btxx,btxy,btyy,sterm)
 
    !compute scalar product of v and tau (x-component)
    this%amin(:,:)  = pvar(:,:,this%XVELOCITY)*btxx(:,:) &
                    + pvar(:,:,this%YVELOCITY)*btxy(:,:) 

   !compute scalar product of v and tau (y-component)
    this%amax(:,:) = pvar(:,:,this%XVELOCITY)*btxy(:,:) &
                   + pvar(:,:,this%YVELOCITY)*btyy(:,:)

    ! compute vector divergence of scalar product v and tau
!CDIR IEXPAND
    CALL Divergence(Mesh,this%amin(:,:),this%amax(:,:),sterm(:,:,this%ENERGY))
  END SUBROUTINE ViscositySources_euler2D


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
    CALL Cons2Prim_euler2D(this%gamma,cvar(i1:i2,j1:j2,this%DENSITY), &
         cvar(i1:i2,j1:j2,this%XMOMENTUM),cvar(i1:i2,j1:j2,this%YMOMENTUM), &
         cvar(i1:i2,j1:j2,this%ENERGY),pvar(i1:i2,j1:j2,this%DENSITY), &
         pvar(i1:i2,j1:j2,this%XVELOCITY),pvar(i1:i2,j1:j2,this%YVELOCITY), &
         pvar(i1:i2,j1:j2,this%PRESSURE))
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
    CALL Cons2Prim_euler2D(this%gamma,cons(i1:i2,j1:j2,:,this%DENSITY), &
         cons(i1:i2,j1:j2,:,this%XMOMENTUM),cons(i1:i2,j1:j2,:,this%YMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%ENERGY),prim(i1:i2,j1:j2,:,this%DENSITY), &
         prim(i1:i2,j1:j2,:,this%XVELOCITY),prim(i1:i2,j1:j2,:,this%YVELOCITY), &
         prim(i1:i2,j1:j2,:,this%PRESSURE))
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
    CALL Prim2Cons_euler2D(this%gamma,pvar(i1:i2,j1:j2,this%DENSITY), &
         pvar(i1:i2,j1:j2,this%XVELOCITY),pvar(i1:i2,j1:j2,this%YVELOCITY), &
         pvar(i1:i2,j1:j2,this%PRESSURE),cvar(i1:i2,j1:j2,this%DENSITY), &
         cvar(i1:i2,j1:j2,this%XMOMENTUM),cvar(i1:i2,j1:j2,this%YMOMENTUM), &
         cvar(i1:i2,j1:j2,this%ENERGY))
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
    CALL Prim2Cons_euler2D(this%gamma,prim(i1:i2,j1:j2,:,this%DENSITY), &
         prim(i1:i2,j1:j2,:,this%XVELOCITY),prim(i1:i2,j1:j2,:,this%YVELOCITY), &
         prim(i1:i2,j1:j2,:,this%PRESSURE),cons(i1:i2,j1:j2,:,this%DENSITY), &
         cons(i1:i2,j1:j2,:,this%XMOMENTUM),cons(i1:i2,j1:j2,:,this%YMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%ENERGY))
  END SUBROUTINE Convert2Conservative_faces


  PURE SUBROUTINE ReflectionMasks_euler2D(this,reflX,reflY)
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
    reflX(this%PRESSURE)  = .FALSE.
    ! southern / northern boundary
    reflY(this%DENSITY)   = .FALSE.
    reflY(this%XVELOCITY) = .FALSE.
    reflY(this%YVELOCITY) = .TRUE.
    reflY(this%PRESSURE)  = .FALSE.
  END SUBROUTINE ReflectionMasks_euler2D


  PURE SUBROUTINE AxisMasks_euler2D(this,reflX,reflY)
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
    reflX(this%PRESSURE)  = .FALSE.
    ! southern / northern boundary
    reflY(this%DENSITY)   = .FALSE.
    reflY(this%XVELOCITY) = .TRUE.
    reflY(this%YVELOCITY) = .TRUE.
    reflY(this%PRESSURE)  = .FALSE.
  END SUBROUTINE AxisMasks_euler2D


  ELEMENTAL FUNCTION GetSoundSpeed_euler2D(gamma,density,pressure) RESULT(cs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,density,pressure
    REAL :: cs
    !------------------------------------------------------------------------!
    cs = SQRT(MAX(2.0*TINY(cs),gamma*pressure/density))
  END FUNCTION GetSoundSpeed_euler2D

  
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
    CALL ClosePhysics_euler2Dit(this)
  END SUBROUTINE ClosePhysics_euler2D

END MODULE physics_euler2D
