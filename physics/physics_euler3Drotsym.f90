!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler3Drotsym.f90                                         #
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
! basic module for 3D Euler equations with rotational symmetry
!----------------------------------------------------------------------------!
MODULE physics_euler3Drotsym
  USE physics_common
  USE mesh_common, ONLY : Mesh_TYP
  USE sources_common, ONLY : Sources_TYP
  USE physics_euler2D, &
       CheckData_euler3Drs => CheckData_euler2D, &
       CalculateWaveSpeeds_euler3Drs => CalculateWaveSpeeds_euler2D
  USE mesh_generic
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
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 3D w/ rot. symmetry"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitPhysics_euler3Drs, &
       CheckData_euler3Drs, &
       CalculateWaveSpeeds_euler3Drs, &
       CalculateFluxesX_euler3Drs, &
       CalculateFluxesY_euler3Drs, &
       CalculateStresses_euler3drs, &
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

  SUBROUTINE InitPhysics_euler3Drs(this,Mesh,problem)
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

    ! allocate memory for arrays used in Euler3Drotsym
    ALLOCATE(this%csound(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4), &
         this%fcent(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,2), &
         this%structure(4),this%errormap(0:3),STAT = err)
    ! abort if allocation fails
    IF (err.NE.0) &
         CALL Error(this, "InitPhysics_euler3Drs", "Unable to allocate memory.")

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
    this%structure(3)%dim  = 3
    this%structure(4)%name = "pressure"
    this%structure(4)%pos  = this%PRESSURE
    this%structure(4)%rank = 0
    this%structure(4)%dim  = 1
                                                                              
    !set errormapping (CheckData Symbols)                                     
    this%errormap(0) = 'X'
    this%errormap(1) = 'D'
    this%errormap(2) = 'P'
    this%errormap(3) = 'B'
  END SUBROUTINE InitPhysics_euler3Drs


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


  PURE SUBROUTINE CalculateStresses_euler3Drs(this,Mesh,pvar,dynvis,bulkvis, &
       btxx,btxy,btxz,btyy,btyz,btzz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
         dynvis,bulkvis,btxx,btxy,btxz,btyy,btyz,btzz
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,dynvis,bulkvis
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: btxx,btxy,btxz,btyy,btyz,btzz
    !------------------------------------------------------------------------!
    ! compute components of the stress tensor at cell bary centers
    ! inside the computational domain including one slice of ghost cells

    ! compute bulk viscosity first and store the result in this%amin
!CDIR IEXPAND
    CALL Divergence(Mesh,pvar(:,:,this%XVELOCITY),pvar(:,:,this%YVELOCITY),this%amin(:,:))
    this%amin(:,:) = bulkvis(:,:)*this%amin(:,:)

!CDIR OUTERUNROLL=8
    DO j=Mesh%JMIN-1,Mesh%JMAX+1
!CDIR NODEP
       DO i=Mesh%IMIN-1,Mesh%IMAX+1
          ! compute the diagonal elements of the stress tensor
          btxx(i,j) = dynvis(i,j) * &
               ( (pvar(i+1,j,this%XVELOCITY) - pvar(i-1,j,this%XVELOCITY)) / Mesh%dlx(i,j) &
               + 2.0 * Mesh%cxyx(i,j,1) * pvar(i,j,this%YVELOCITY) ) &
               + this%amin(i,j) ! bulk viscosity contribution
               
          btyy(i,j) = dynvis(i,j) * &
               ( (pvar(i,j+1,this%YVELOCITY) - pvar(i,j-1,this%YVELOCITY)) / Mesh%dly(i,j) &
               + 2.0 * Mesh%cyxy(i,j,1) * pvar(i,j,this%XVELOCITY) ) &
               + this%amin(i,j) ! bulk viscosity contribution

          btzz(i,j) = dynvis(i,j) * &
               ( 2.0 * ( Mesh%czxz(i,j,1) * pvar(i,j,this%XVELOCITY) ) &
               + Mesh%czyz(i,j,1) * pvar(i,j,this%YVELOCITY) ) &
               + this%amin(i,j) ! bulk viscosity contribution

          ! compute the off-diagonal elements (no bulk viscosity)
          btxy(i,j) = dynvis(i,j) * ( 0.5 * &
               ( (pvar(i+1,j,this%YVELOCITY) - pvar(i-1,j,this%YVELOCITY)) / Mesh%dlx(i,j) &
               + (pvar(i,j+1,this%XVELOCITY) - pvar(i,j-1,this%XVELOCITY)) / Mesh%dly(i,j) ) &
               - Mesh%cxyx(i,j,1) * pvar(i,j,this%XVELOCITY) &
               - Mesh%cyxy(i,j,1) * pvar(i,j,this%YVELOCITY) )

          btxz(i,j) = dynvis(i,j) * ( 0.5 * &
               ( (pvar(i+1,j,this%ZVELOCITY) - pvar(i-1,j,this%ZVELOCITY)) / Mesh%dlx(i,j) ) &
               - Mesh%czxz(i,j,1) * pvar(i,j,this%ZVELOCITY) )

          btyz(i,j) = dynvis(i,j) * ( 0.5 * &
               ( (pvar(i,j+1,this%ZVELOCITY) - pvar(i,j-1,this%ZVELOCITY)) / Mesh%dly(i,j) ) &
               - Mesh%czyz(i,j,1) * pvar(i,j,this%ZVELOCITY) )
       END DO
    END DO
  END SUBROUTINE CalculateStresses_euler3Drs


  PURE SUBROUTINE GeometricalSources_center(this,Mesh,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
         :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,cvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! calculate centrifugal forces
          CALL CentrifugalForces_euler3Drs(pvar(i,j,this%DENSITY),pvar(i,j,this%ZVELOCITY), &
               Mesh%czxz(i,j,1),Mesh%czyz(i,j,1),this%fcent(i,j,1,1),this%fcent(i,j,2,1))
          
          ! no geometrical sources in continuity and energy equations
          sterm(i,j,this%DENSITY) = 0.0
          sterm(i,j,this%ENERGY)  = 0.0

          ! geometrical source terms in momentum equationes
          ! with centrifugal forces
          sterm(i,j,this%XMOMENTUM) = MomentumSourcesX_euler2D(cvar(i,j,this%YMOMENTUM), &
               pvar(i,j,this%XVELOCITY),pvar(i,j,this%YVELOCITY),pvar(i,j,this%PRESSURE), &
               Mesh%cxyx(i,j,1),Mesh%cyxy(i,j,1),Mesh%czxz(i,j,1)) + this%fcent(i,j,1,1)
          sterm(i,j,this%YMOMENTUM) = MomentumSourcesY_euler2D(cvar(i,j,this%XMOMENTUM), &
               pvar(i,j,this%XVELOCITY),pvar(i,j,this%YVELOCITY),pvar(i,j,this%PRESSURE), &
               Mesh%cxyx(i,j,1),Mesh%cyxy(i,j,1),Mesh%czyz(i,j,1)) + this%fcent(i,j,2,1)
          sterm(i,j,this%ZMOMENTUM) = -pvar(i,j,this%ZVELOCITY)*(Mesh%czxz(i,j,1) * &
               cvar(i,j,this%XMOMENTUM) + Mesh%czyz(i,j,1) * cvar(i,j,this%YMOMENTUM))
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
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,prim,cons
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! no geometrical sources in continuity and energy equations
          sterm(i,j,this%DENSITY) = 0.0
          sterm(i,j,this%ENERGY)  = 0.0
!CDIR UNROLL=4
          DO k=1,4
            ! calculate centrifugal forces
            CALL CentrifugalForces_euler3Drs(prim(i,j,k,this%DENSITY), &
                  prim(i,j,k,this%ZVELOCITY),Mesh%czxz(i,j,k),Mesh%czyz(i,j,k), &
                  this%fcent(i,j,k,1),this%fcent(i,j,k,2))
            ! add momentum sources to centrifugal forces
            this%fcent(i,j,k,1) = this%fcent(i,j,k,1) + MomentumSourcesX_euler2D( &
                  cons(i,j,k,this%YMOMENTUM),prim(i,j,k,this%XVELOCITY), &
                  prim(i,j,k,this%YVELOCITY),prim(i,j,k,this%PRESSURE), &
                  Mesh%cxyx(i,j,k),Mesh%cyxy(i,j,k),Mesh%czxz(i,j,k))

            this%fcent(i,j,k,2) = this%fcent(i,j,k,2) + MomentumSourcesY_euler2D( &
                  cons(i,j,k,this%XMOMENTUM),prim(i,j,k,this%XVELOCITY), &
                  prim(i,j,k,this%YVELOCITY),prim(i,j,k,this%PRESSURE), &
                  Mesh%cxyx(i,j,k),Mesh%cyxy(i,j,k),Mesh%czyz(i,j,k))


          END DO

          ! momentum sources (sum up corner values, don't use SUM function,
          ! because it prevents COLLAPSING and causes poor vectorization
          sterm(i,j,this%XMOMENTUM) = this%fcent(i,j,1,1)+this%fcent(i,j,2,1) &
                  +this%fcent(i,j,3,1)+this%fcent(i,j,4,1)
          sterm(i,j,this%YMOMENTUM) = this%fcent(i,j,1,2)+this%fcent(i,j,2,2) &
                  +this%fcent(i,j,3,2)+this%fcent(i,j,4,2)
          sterm(i,j,this%ZMOMENTUM) = -prim(i,j,1,this%ZVELOCITY) &
                  * (Mesh%czxz(i,j,1)*cons(i,j,1,this%XMOMENTUM) &
                    +Mesh%czyz(i,j,1)*cons(i,j,1,this%YMOMENTUM)) &
                - prim(i,j,2,this%ZVELOCITY) &
                  * (Mesh%czxz(i,j,2)*cons(i,j,2,this%XMOMENTUM) &
                    +Mesh%czyz(i,j,2)*cons(i,j,2,this%YMOMENTUM)) &
                - prim(i,j,3,this%ZVELOCITY) &
                  * (Mesh%czxz(i,j,3)*cons(i,j,3,this%XMOMENTUM) &
                    +Mesh%czyz(i,j,3)*cons(i,j,3,this%YMOMENTUM)) &
                - prim(i,j,4,this%ZVELOCITY) &
                  * (Mesh%czxz(i,j,4)*cons(i,j,4,this%XMOMENTUM) &
                    +Mesh%czyz(i,j,4)*cons(i,j,4,this%YMOMENTUM))
       END DO
    END DO
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


  PURE SUBROUTINE ViscositySources_euler3Drs(this,Mesh,pvar,btxx,btxy,btxz, &
       btyy,btyz,btzz,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: &
         pvar,sterm
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
         btxx,btxy,btxz,btyy,btyz,btzz
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar, btxx,btxy,btxz,btyy,btyz,btzz
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! viscosity source terms
    sterm(:,:,this%DENSITY) = 0.0 
 
    ! compute 3D tensor divergence for symmetric tensor btxy=btyx, btxz=btzx, etc.
!CDIR IEXPAND
    CALL Divergence(Mesh,btxx,btxy,btxz,btxy,btyy,btyz,btxz,btyz,btzz, &
                    sterm(:,:,this%XMOMENTUM),sterm(:,:,this%YMOMENTUM), &
                    sterm(:,:,this%ZMOMENTUM))
 
    ! compute scalar product of v and btx (x-component)
    ! use amin for temporary storage
!CDIR NODEP
    this%amin(:,:)  = pvar(:,:,this%XVELOCITY)*btxx(:,:) &
                    + pvar(:,:,this%YVELOCITY)*btxy(:,:) &
                    + pvar(:,:,this%ZVELOCITY)*btxz(:,:)

   ! compute scalar product of v and bty (y-component)
   ! use amax for temporary storage
!CDIR NODEP 
    this%amax(:,:) = pvar(:,:,this%XVELOCITY)*btxy(:,:) &
                   + pvar(:,:,this%YVELOCITY)*btyy(:,:) &
                   + pvar(:,:,this%ZVELOCITY)*btyz(:,:) 
 
    ! compute vector divergence of scalar product v_i * bt_ij
!CDIR IEXPAND
    CALL Divergence(Mesh,this%amin,this%amax,sterm(:,:,this%ENERGY))
  END SUBROUTINE ViscositySources_euler3Drs


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
    CALL Cons2Prim_euler3Drs(this%gamma,cvar(i1:i2,j1:j2,this%DENSITY), &
         cvar(i1:i2,j1:j2,this%XMOMENTUM),cvar(i1:i2,j1:j2,this%YMOMENTUM), &
         cvar(i1:i2,j1:j2,this%ZMOMENTUM),cvar(i1:i2,j1:j2,this%ENERGY), &
         pvar(i1:i2,j1:j2,this%DENSITY),pvar(i1:i2,j1:j2,this%XVELOCITY), &
         pvar(i1:i2,j1:j2,this%YVELOCITY),pvar(i1:i2,j1:j2,this%ZVELOCITY), &
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
    CALL Cons2Prim_euler3Drs(this%gamma,cons(i1:i2,j1:j2,:,this%DENSITY), &
         cons(i1:i2,j1:j2,:,this%XMOMENTUM),cons(i1:i2,j1:j2,:,this%YMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%ZMOMENTUM),cons(i1:i2,j1:j2,:,this%ENERGY), &
         prim(i1:i2,j1:j2,:,this%DENSITY),prim(i1:i2,j1:j2,:,this%XVELOCITY), &
         prim(i1:i2,j1:j2,:,this%YVELOCITY),prim(i1:i2,j1:j2,:,this%ZVELOCITY), &
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
    CALL Prim2Cons_euler3Drs(this%gamma,pvar(i1:i2,j1:j2,this%DENSITY), &
         pvar(i1:i2,j1:j2,this%XVELOCITY),pvar(i1:i2,j1:j2,this%YVELOCITY), &
         pvar(i1:i2,j1:j2,this%ZVELOCITY),pvar(i1:i2,j1:j2,this%PRESSURE), &
         cvar(i1:i2,j1:j2,this%DENSITY),cvar(i1:i2,j1:j2,this%XMOMENTUM), &
         cvar(i1:i2,j1:j2,this%YMOMENTUM),cvar(i1:i2,j1:j2,this%ZMOMENTUM), &
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
    CALL Prim2Cons_euler3Drs(this%gamma,prim(i1:i2,j1:j2,:,this%DENSITY), &
         prim(i1:i2,j1:j2,:,this%XVELOCITY),prim(i1:i2,j1:j2,:,this%YVELOCITY), &
         prim(i1:i2,j1:j2,:,this%ZVELOCITY),prim(i1:i2,j1:j2,:,this%PRESSURE), &
         cons(i1:i2,j1:j2,:,this%DENSITY),cons(i1:i2,j1:j2,:,this%XMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%YMOMENTUM),cons(i1:i2,j1:j2,:,this%ZMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%ENERGY))
  END SUBROUTINE Convert2Conservative_faces


  PURE SUBROUTINE ReflectionMasks_euler3Drs(this,reflX,reflY)
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
    reflX(this%ZVELOCITY) = .FALSE.
    reflX(this%PRESSURE)  = .FALSE.
    ! southern / northern boundary
    reflY(this%DENSITY)   = .FALSE.
    reflY(this%XVELOCITY) = .FALSE.
    reflY(this%YVELOCITY) = .TRUE.
    reflY(this%ZVELOCITY) = .FALSE.
    reflY(this%PRESSURE)  = .FALSE.
  END SUBROUTINE ReflectionMasks_euler3Drs


  PURE SUBROUTINE AxisMasks_euler3Drs(this,reflX,reflY)
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
    reflX(this%ZVELOCITY) = .TRUE.
    reflX(this%PRESSURE)  = .FALSE.
    ! southern / northern boundary
    reflY(this%DENSITY)   = .FALSE.
    reflY(this%XVELOCITY) = .FALSE.
    reflY(this%YVELOCITY) = .TRUE.
    reflY(this%ZVELOCITY) = .TRUE.
    reflY(this%PRESSURE)  = .FALSE.
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


  ELEMENTAL SUBROUTINE CentrifugalForces_euler3Drs(rho,vz,czxz,czyz,fcx,fcy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho,vz,czxz,czyz
    REAL, INTENT(OUT) :: fcx,fcy
    !------------------------------------------------------------------------!
    fcx  = rho*vz*vz ! temporary 
    fcy  = fcx * czyz
    fcx  = fcx * czxz
  END SUBROUTINE CentrifugalForces_euler3Drs


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
