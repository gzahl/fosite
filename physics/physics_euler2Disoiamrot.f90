!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler2Dlocisoiamt2.f90                                   #
!#                                                                           #
!# Copyright (C) 2013                                                        #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
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
!! - gas dynamics with rotating frame of reference
!!   \key{omega,REAL,angular speed of rotating frame of reference,0.0}
!!   \key{centrot_x,REAL,cartesian x-coordiante for center of rotation,0.0}
!!   \key{centrot_y,REAL,cartesian y-coordiante for center of rotation,0.0}
!!   \key{softening,REAL,softening parameter to smooth out singularity
!!        near center of rotation,1.0}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Manuel Jung
!!
!! \brief basic module for 2D isothermal Euler problem with inertial angular momentum transport
!!
!! \extends physics_common
!! \ingroup physics
!----------------------------------------------------------------------------!
MODULE physics_euler2Disoiamrot
  USE physics_common
  USE mesh_common, ONLY : Mesh_TYP
  USE sources_common, ONLY : Sources_TYP
  USE common_dict
  USE physics_euler2Disothm, ONLY : &
       InitPhysics_euler2Dit, &
       ClosePhysics_euler2Dit, &
       CalcWaveSpeeds_euler2Ditiar => CalcWaveSpeeds_euler2Dit, &
       CalcStresses_euler2Ditiar => CalcStresses_euler2Dit, &
       ViscositySources_euler2Dit, &
       ReflectionMasks_euler2Ditiar =>  ReflectionMasks_euler2Dit, &
       AxisMasks_euler2Ditiar => AxisMasks_euler2Dit, &
       SetWaveSpeeds_euler2Ditiar => SetWaveSpeeds_euler2Dit
  USE mesh_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE GeometricalSrcs_euler2Ditiar
     MODULE PROCEDURE GeometricalSources_center
     MODULE PROCEDURE GeometricalSources_faces
  END INTERFACE
  INTERFACE Convert2Primitive_euler2Ditiar
     MODULE PROCEDURE Convert2Primitive_center
     MODULE PROCEDURE Convert2Primitive_faces
  END INTERFACE
  INTERFACE Convert2Cons_euler2Ditiar
     MODULE PROCEDURE Convert2Conservative_center
     MODULE PROCEDURE Convert2Conservative_faces
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: num_var = 3              ! number of variables       !
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 2D lociso w/ ang. mom."
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Physics_TYP, &
       ! methods
       InitPhysics_euler2Ditiar, &
       ClosePhysics_euler2Ditiar, &
       CalcWaveSpeeds_euler2Ditiar, &
       CalculateFluxesX_euler2Ditiar, &
       CalculateFluxesY_euler2Ditiar, &
       CalcStresses_euler2Ditiar, &
       GeometricalSrcs_euler2Ditiar, &
       ViscositySources_euler2Ditiar, &
       ExternalSources_euler2Ditiar, &
       Convert2Primitive_euler2Ditiar, &
       Convert2Cons_euler2Ditiar, &
       ReflectionMasks_euler2Ditiar, &
       AxisMasks_euler2Ditiar, &
       SetWaveSpeeds_euler2Ditiar
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitPhysics_euler2Ditiar(this,Mesh,problem)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: problem
    !------------------------------------------------------------------------!
    INTEGER           :: err,mask(4)
    REAL              :: dr,x,y
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: bccart
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,2) :: fccart
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,problem
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL InitPhysics_euler2Dit(this,Mesh,problem,problem_name,num_var)

    ! set names for conservative variables
    this%cvarname(this%DENSITY)   = "density"
    !this%cvarname(this%XMOMENTUM) = "radialmomentum"
    !this%cvarname(this%YMOMENTUM) = "iangularmomentum"
    this%cvarname(this%XMOMENTUM) = "xmomentum"
    this%cvarname(this%YMOMENTUM) = "ymomentum"

    ! allocate memory
    ALLOCATE(this%bcposvec(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
        this%bcradius(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
        this%fposvec(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,2), &
        this%fradius(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4), &
        this%divposvec(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
        STAT = err)
    ! abort if allocation fails
    IF (err.NE.0) &
         CALL Error(this, "InitPhysics_euler2Ditiar", "Unable to allocate memory.")
  
    ! translate the position vector to the center of rotation
    bccart(:,:,1) = Mesh%bccart(:,:,1) - this%centrot(1)
    bccart(:,:,2) = Mesh%bccart(:,:,2) - this%centrot(2)
    fccart(:,:,:,1) = Mesh%fccart(:,:,:,1) - this%centrot(1)
    fccart(:,:,:,2) = Mesh%fccart(:,:,:,2) - this%centrot(2)  
    ! compute curvilinear components of translated position vectors
    CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,bccart,this%bcposvec)
    ! get curvilinear coordinates of the center of rotation
    CALL Convert2Curvilinear(Mesh%geometry,this%centrot(1),this%centrot(2),x,y)
    ! determine the softening length to avoid devision by zero when r -> 0
    this%bcradius(:,:) = this%bcposvec(:,:,1)**2 + this%bcposvec(:,:,2)**2
    dr = SQRT(MINVAL(this%bcradius(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)))
    ! check if the center of rotation lies within the computational domain
    ! including the innermost ghost cells
    mask(1) = Mesh%IMIN-1
    mask(2) = Mesh%IMAX+1
    mask(3) = Mesh%JMIN-1
    mask(4) = Mesh%JMAX+1
    IF (InternalPoint(Mesh,x,y,mask)) THEN
       dr = this%eps*0.5*dr ! factor 0.5 was determined experimental
    ELSE
       ! otherwise disable softening
       dr = 0.0
    END IF
    ! compute curvilinear components of the translated position vector using its
    ! cartesian values and the distance to the center of rotation using the
    ! softening length to smooth out the 1/r singularity as r -> 0
    this%bcradius(:,:) = SQRT(this%bcradius(:,:)+dr**2)
    this%bcposvec(:,:,1) = this%bcposvec(:,:,1)/this%bcradius(:,:)
    this%bcposvec(:,:,2) = this%bcposvec(:,:,2)/this%bcradius(:,:)
    ! curvilinear components of face centered translated position vector
    CALL Convert2Curvilinear(Mesh%geometry,Mesh%fpos,fccart,this%fposvec)
    this%fradius(:,:,:) = SQRT(this%fposvec(:,:,:,1)**2 + this%fposvec(:,:,:,2)**2+dr**2)
    this%fposvec(:,:,:,1) = this%fposvec(:,:,:,1)/this%fradius(:,:,:)
    this%fposvec(:,:,:,2) = this%fposvec(:,:,:,2)/this%fradius(:,:,:)
    ! compute the divergence of the translated position vector
    CALL Divergence(Mesh,this%fposvec,this%divposvec)
!!$    CALL Divergence(Mesh,this%bcposvec(:,:,1),this%bcposvec(:,:,2),this%divposvec)
!!$    this%divposvec = 2.0
    ! overwrite face centered values of pos. vect. with corner values for
    ! trapezoidal mesh
    SELECT CASE(GetType(Mesh))
    CASE(MIDPOINT)
       ! do nothing
    CASE(TRAPEZOIDAL)
!!$ FIXME: implement source terms for trapezoidal rule
       CALL Error(this,"InitPhysics_euler2Ditiar", "Trapezoidal mesh not supported yet.")
       fccart(:,:,:,1) = Mesh%ccart(:,:,:,1) - this%centrot(1)
       fccart(:,:,:,2) = Mesh%ccart(:,:,:,2) - this%centrot(2)  
       CALL Convert2Curvilinear(Mesh%geometry,Mesh%cpos,fccart,this%fposvec)
       this%fradius(:,:,:) = SQRT(this%fposvec(:,:,:,1)**2 + this%fposvec(:,:,:,2)**2+dr**2)
       this%fposvec(:,:,:,1) = this%fposvec(:,:,:,1)/this%fradius(:,:,:)
       this%fposvec(:,:,:,2) = this%fposvec(:,:,:,2)/this%fradius(:,:,:)
    CASE DEFAULT
       CALL Error(this,"InitPhysics_euler2Ditiar", "Mesh not supported.")
    END SELECT
  END SUBROUTINE InitPhysics_euler2Ditiar


  PURE SUBROUTINE CalculateFluxesX_euler2Ditiar(this,Mesh,nmin,nmax,prim,cons,xfluxes)
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
    CALL CalculateFlux_euler2Ditiar(this%fcsound(:,:,nmin:nmax), &
         this%fposvec(:,:,nmin:nmax,1), &
         this%fposvec(:,:,nmin:nmax,2),this%fradius(:,:,nmin:nmax), &
         prim(:,:,nmin:nmax,this%DENSITY), &
         prim(:,:,nmin:nmax,this%XVELOCITY), &
         cons(:,:,nmin:nmax,this%XMOMENTUM),cons(:,:,nmin:nmax,this%YMOMENTUM), &
         xfluxes(:,:,nmin:nmax,this%DENSITY),&
         xfluxes(:,:,nmin:nmax,this%XMOMENTUM),xfluxes(:,:,nmin:nmax,this%YMOMENTUM))
  END SUBROUTINE CalculateFluxesX_euler2Ditiar


  PURE SUBROUTINE CalculateFluxesY_euler2Ditiar(this,Mesh,nmin,nmax,prim,cons,yfluxes)
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
    ! the x- and y-flux functions are almost identical, hence we use the
    ! same elemental subroutine, but replace 
    !    rx -> ry and ry -> -rx
    !    vx -> vy
    CALL CalculateFlux_euler2Ditiar(this%fcsound(:,:,nmin:nmax),&
         this%fposvec(:,:,nmin:nmax,2), &
         -this%fposvec(:,:,nmin:nmax,1),this%fradius(:,:,nmin:nmax), &
         prim(:,:,nmin:nmax,this%DENSITY), &
         prim(:,:,nmin:nmax,this%YVELOCITY), &
         cons(:,:,nmin:nmax,this%XMOMENTUM),cons(:,:,nmin:nmax,this%YMOMENTUM), &
         yfluxes(:,:,nmin:nmax,this%DENSITY), &
         yfluxes(:,:,nmin:nmax,this%XMOMENTUM),yfluxes(:,:,nmin:nmax,this%YMOMENTUM))
  END SUBROUTINE CalculateFluxesY_euler2Ditiar


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
          ! no geometrical density, angular momentum and rothalpy sources
          sterm(i,j,this%DENSITY)   = 0.
          sterm(i,j,this%YMOMENTUM) = 0.
          ! geometrical source terms in momentum equationes
          sterm(i,j,this%XMOMENTUM) = MomentumSources_euler2Ditiar(&
              this%bccsound(i,j),&
              this%bcposvec(i,j,1),this%bcposvec(i,j,2), &
              this%bcradius(i,j),this%divposvec(i,j),this%Omega, &
              pvar(i,j,this%DENSITY),pvar(i,j,this%XVELOCITY), &
              pvar(i,j,this%YVELOCITY))
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
          sterm(i,j,this%YMOMENTUM) = 0.
!!! FIXME: implement momentum source for trapezoidal mesh using corner values
!!$          sterm(i,j,this%XMOMENTUM) = ...
       END DO
    END DO 
  END SUBROUTINE GeometricalSources_faces


  ! momentum and energy sources due to external force
  PURE SUBROUTINE ExternalSources_euler2Ditiar(this,Mesh,accel,pvar,cvar,sterm)
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
          sterm(i,j,this%XMOMENTUM) = pvar(i,j,this%DENSITY) * ( &
              this%bcposvec(i,j,1)*accel(i,j,1) + this%bcposvec(i,j,2)*accel(i,j,2))
          sterm(i,j,this%YMOMENTUM) = pvar(i,j,this%DENSITY) * this%bcradius(i,j) * ( &
              this%bcposvec(i,j,1)*accel(i,j,2) - this%bcposvec(i,j,2)*accel(i,j,1))
       END DO
    END DO
  END SUBROUTINE ExternalSources_euler2Ditiar


  PURE SUBROUTINE ViscositySources_euler2Ditiar(this,Mesh,pvar,btxx,btxy,btyy,sterm)         
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: &
         pvar,sterm
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
         btxx,btxy,btyy
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    REAL              :: tmp
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,btxx,btxy,btyy
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL ViscositySources_euler2Dit(this,Mesh,pvar,btxx,btxy,btyy,sterm)
 
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! viscous radial momentum sources:
          ! project the 2D source term onto radial direction
          tmp = this%bcposvec(i,j,1)*sterm(i,j,this%XVELOCITY) &
                + this%bcposvec(i,j,2)*sterm(i,j,this%YVELOCITY)
          ! viscous angular momentum sources:
          ! project the cross product of r and 2D source term onto vertical direction
          sterm(i,j,this%YVELOCITY) = this%bcradius(i,j) * &
                ( this%bcposvec(i,j,1)*sterm(i,j,this%YVELOCITY) &
                - this%bcposvec(i,j,2)*sterm(i,j,this%XVELOCITY) )
          ! copy temporary result
          sterm(i,j,this%XVELOCITY) = tmp
       END DO
    END DO
  END SUBROUTINE ViscositySources_euler2Ditiar


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
    CALL Cons2Prim_euler2Ditiar(this%bcposvec(i1:i2,j1:j2,1), &
         this%bcposvec(i1:i2,j1:j2,2),this%bcradius(i1:i2,j1:j2),this%Omega, &
         this%gamma,cvar(i1:i2,j1:j2,this%DENSITY), &
         cvar(i1:i2,j1:j2,this%XMOMENTUM),cvar(i1:i2,j1:j2,this%YMOMENTUM), &
         pvar(i1:i2,j1:j2,this%DENSITY), &
         pvar(i1:i2,j1:j2,this%XVELOCITY),pvar(i1:i2,j1:j2,this%YVELOCITY))
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
    CALL Cons2Prim_euler2Ditiar(this%fposvec(i1:i2,j1:j2,:,1), &
         this%fposvec(i1:i2,j1:j2,:,2),this%fradius(i1:i2,j1:j2,:),this%Omega, &
         this%gamma,cons(i1:i2,j1:j2,:,this%DENSITY), &
         cons(i1:i2,j1:j2,:,this%XMOMENTUM),cons(i1:i2,j1:j2,:,this%YMOMENTUM), &
         prim(i1:i2,j1:j2,:,this%DENSITY), &
         prim(i1:i2,j1:j2,:,this%XVELOCITY),prim(i1:i2,j1:j2,:,this%YVELOCITY))
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
    CALL Prim2Cons_euler2Ditiar(this%bcposvec(i1:i2,j1:j2,1), &
         this%bcposvec(i1:i2,j1:j2,2),this%bcradius(i1:i2,j1:j2),this%Omega, &
         this%gamma,pvar(i1:i2,j1:j2,this%DENSITY), &
         pvar(i1:i2,j1:j2,this%XVELOCITY),pvar(i1:i2,j1:j2,this%YVELOCITY), &
         cvar(i1:i2,j1:j2,this%DENSITY), &
         cvar(i1:i2,j1:j2,this%XMOMENTUM),cvar(i1:i2,j1:j2,this%YMOMENTUM))
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
    CALL Prim2Cons_euler2Ditiar(this%fposvec(i1:i2,j1:j2,:,1), &
         this%fposvec(i1:i2,j1:j2,:,2),this%fradius(i1:i2,j1:j2,:),this%Omega, &
         this%gamma,prim(i1:i2,j1:j2,:,this%DENSITY), &
         prim(i1:i2,j1:j2,:,this%XVELOCITY),prim(i1:i2,j1:j2,:,this%YVELOCITY), &
         cons(i1:i2,j1:j2,:,this%DENSITY), &
         cons(i1:i2,j1:j2,:,this%XMOMENTUM),cons(i1:i2,j1:j2,:,this%YMOMENTUM))
  END SUBROUTINE Convert2Conservative_faces


  ELEMENTAL SUBROUTINE CalculateFlux_euler2Ditiar(cs,erx,ery,r,rho,v,Mr,Lz,f1,f2,f3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL        :: P
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: erx,ery,r,rho,v,Mr,Lz,cs
    REAL, INTENT(OUT) :: f1,f2,f3
    !------------------------------------------------------------------------!
    P = rho * cs * cs
    f1 = rho*v
    f2 = Mr*v + P*erx
    f3 = Lz*v - P*r*ery
  END SUBROUTINE CalculateFlux_euler2Ditiar


  ! radial momentum source term
  ELEMENTAL FUNCTION MomentumSources_euler2Ditiar(cs,erx,ery,r,divr,Omega,rho,u,v) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: erx,ery,r,divr,Omega,rho,u,v,cs
    REAL :: st
    !------------------------------------------------------------------------!
    REAL :: Omega_if
    !------------------------------------------------------------------------!
    ! inertial frame angular frequency
    Omega_if = (-ery*u + erx*v) / r + Omega
    st = divr*rho*cs*cs + rho*r*Omega_if*Omega_if
  END FUNCTION MomentumSources_euler2Ditiar


  ELEMENTAL SUBROUTINE Cons2Prim_euler2Ditiar(erx,ery,r,Omega,gamma,rho_in,Mr,Lz,rho_out,u,v)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)    :: erx,ery,r,Omega,gamma,rho_in,Mr,Lz
    REAL, INTENT(OUT)   :: rho_out,u,v
! REAL, INTENT(INOUT) :: u,v
! REAL, INTENT(OUT)   :: rho_out,P
    !------------------------------------------------------------------------!
    REAL :: inv_rho,vphi,vr
    !------------------------------------------------------------------------!
    ! some temporary variables
    inv_rho = 1./rho_in                    ! inverse of density
    vr   = Mr*inv_rho                      ! radial velocity
    vphi = inv_rho*(Lz/r)-r*Omega          ! rotating frame angular velocity
    ! conversion to density, 2D curvilinear velocity components and pressure
    rho_out = rho_in
    u = erx*vr - ery*vphi
    v = ery*vr + erx*vphi
  END SUBROUTINE Cons2Prim_euler2Ditiar

  
  ELEMENTAL SUBROUTINE Prim2Cons_euler2Ditiar(erx,ery,r,Omega,gamma,rho_in,u,v,rho_out,Mr,Lz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: erx,ery,r,Omega,gamma,rho_in,u,v
    REAL, INTENT(OUT) :: rho_out,Mr,Lz
    !------------------------------------------------------------------------!
    rho_out = rho_in
    Mr = rho_in * (erx*u + ery*v)
    Lz = rho_in * r * (erx*v - ery*u + r*Omega)
  END SUBROUTINE Prim2Cons_euler2Ditiar
  

  SUBROUTINE ClosePhysics_euler2Ditiar(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%bcposvec,this%bcradius,this%fposvec,this%fradius, &
               this%divposvec)
    CALL ClosePhysics_euler2Dit(this)
  END SUBROUTINE ClosePhysics_euler2Ditiar

END MODULE physics_euler2Disoiamrot
