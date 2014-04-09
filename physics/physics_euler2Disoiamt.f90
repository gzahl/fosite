!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler2Disoiamt.f90                                        #
!#                                                                           #
!# Copyright (C) 2012, 2013                                                  #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
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
!> \author Manuel Jung
!! \author Tobias Illenseer
!!
!! \brief basic module for 2D isothermal Euler equations with inertial angular
!! momentum transport
!!
!! \extends physics_common
!! \ingroup physics
!----------------------------------------------------------------------------!
MODULE physics_euler2Disoiamt
  USE physics_common
  USE physics_euler2Disothm, ONLY : &
    InitPhysics_euler2Dit, &
    ViscositySources_euler2Dit, &
    ExternalSources_euler2Dit, &
    CalcFluxesX_euler2Ditia => CalcFluxesX_euler2Dit, &
    CalcStresses_euler2Ditia => CalcStresses_euler2Dit, &
    ReflectionMasks_euler2Ditia => ReflectionMasks_euler2Dit, &
    AxisMasks_euler2Ditia => AxisMasks_euler2Dit, &
    CalcFluxX_euler2Ditia => CalcFlux_euler2Dit, &
    SetWaveSpeeds_euler2Ditia => SetWaveSpeeds_euler2Dit, &
    SetEigenValues_euler2Ditia => SetEigenValues_euler2Dit, &
    SetBoundaryData_euler2Ditia => SetBoundaryData_euler2Dit, &
    SetCharVars_euler2Ditia => SetCharVars_euler2Dit
  USE sources_common, ONLY : Sources_TYP
  USE mesh_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE CalcWaveSpeeds_euler2Ditia
     MODULE PROCEDURE CalcWaveSpeeds_center
     MODULE PROCEDURE CalcWaveSpeeds_faces
  END INTERFACE
  INTERFACE GeometricalSources_euler2Ditia
     MODULE PROCEDURE GeometricalSources_center
     MODULE PROCEDURE GeometricalSources_faces
  END INTERFACE
  INTERFACE Convert2Primitive_euler2Ditia
     MODULE PROCEDURE Convert2Primitive_center
     MODULE PROCEDURE Convert2Primitive_faces
  END INTERFACE
  INTERFACE Convert2Cons_euler2Ditia
     MODULE PROCEDURE Convert2Conservative_center
     MODULE PROCEDURE Convert2Conservative_faces
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: num_var = 3              ! number of variables       !
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 2D iso /w inertial amt"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Physics_TYP, &
       ! methods
       InitPhysics_euler2Ditia, &
       ClosePhysics_euler2Ditia, &
       CalcWaveSpeeds_euler2Ditia, &
       CalcFluxesX_euler2Ditia, &
       CalcFluxesY_euler2Ditia, &
       CalcCharSystemX_euler2Ditia, &
       CalcBoundaryDataX_euler2Ditia, &
       CalcStresses_euler2Ditia, &
       GeometricalSources_euler2Ditia, &
       ViscositySources_euler2Ditia, &
       ExternalSources_euler2Ditia, &
       Convert2Primitive_euler2Ditia, &
       Convert2Cons_euler2Ditia, &
       ReflectionMasks_euler2Ditia, &
       AxisMasks_euler2Ditia, &
       SetWaveSpeeds_euler2Ditia, &
       SetEigenValues_euler2Ditia, &
       CalcFluxX_euler2Ditia, &
       CalcFluxY_euler2Ditia, &
       MomentumSourcesX_euler2Ditia, &
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

  SUBROUTINE InitPhysics_euler2Ditia(this,Mesh,problem,pname,nvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: problem
    CHARACTER(LEN=32), OPTIONAL :: pname
    INTEGER,OPTIONAL  :: nvar
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,problem,pname
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    IF (PRESENT(pname).AND.PRESENT(nvar)) THEN
      CALL InitPhysics_euler2Dit(this,Mesh,problem,pname,nvar)
    ELSE IF (PRESENT(pname).OR.PRESENT(nvar)) THEN
       CALL Error(this, "InitPhysics_euler2Dit", "Both or no optional " &
        // "arguments at all have to be defined.")
    ELSE
      CALL InitPhysics_euler2Dit(this,Mesh,problem,problem_name,num_var)
    END IF

    ! set pointer to either face or corner scale factors
    SELECT CASE(GetType(Mesh))
    CASE(MIDPOINT)
       this%hy => Mesh%fhy
    CASE(TRAPEZOIDAL)
       this%hy => Mesh%chy
    CASE DEFAULT
       CALL Error(this,"InitPhysics_euler2Ditia", "Mesh not supported.")
    END SELECT

    ! allocate memory
    ALLOCATE(this%w(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4),&
             STAT = err)
    ! abort if allocation fails
    IF (err.NE.0) &
         CALL Error(this, "InitPhysics_euler2Ditia", "Unable to allocate memory.")

    this%w = 0.
  END SUBROUTINE InitPhysics_euler2Ditia

  PURE SUBROUTINE CalcWaveSpeeds_center(this,Mesh,pvar)
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
          CALL SetWaveSpeeds_euler2Ditia(&
               this%bccsound(i,j),&
               pvar(i,j,this%XVELOCITY),&
               this%amin(i,j),this%amax(i,j))
          ! y-direction
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2Ditia( &
               this%bccsound(i,j),&
               pvar(i,j,this%YVELOCITY)-Mesh%bhy(i,j)*this%w(i,j,1),&
               this%bmin(i,j),this%bmax(i,j))
       END DO
    END DO
  END SUBROUTINE CalcWaveSpeeds_center

  PURE SUBROUTINE CalcWaveSpeeds_faces(this,Mesh,prim)
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
          CALL SetWaveSpeeds_euler2Ditia(&
               this%fcsound(i,j,1),&
               prim(i,j,1,this%XVELOCITY), &
               this%tmin(i,j,1),this%tmax(i,j,1))
          ! eastern
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2Ditia(&
               this%fcsound(i,j,2),&
               prim(i,j,2,this%XVELOCITY), &
               this%amin(i,j),this%amax(i,j))
          ! southern
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2Ditia(&
               this%fcsound(i,j,3), &
               prim(i,j,3,this%YVELOCITY)-this%hy(i,j,3)*this%w(i,j,3), &
               this%tmin(i,j,2),this%tmax(i,j,2))
          ! northern
!CDIR IEXPAND
          CALL SetWaveSpeeds_euler2Ditia(&
               this%fcsound(i,j,4), &
               prim(i,j,4,this%YVELOCITY)-this%hy(i,j,4)*this%w(i,j,4), &
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
  END SUBROUTINE CalcWaveSpeeds_faces


  PURE SUBROUTINE CalcFluxesY_euler2Ditia(this,Mesh,nmin,nmax,prim,cons,yfluxes)
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
    CALL CalcFluxY_euler2Ditia(this%fcsound(:,:,nmin:nmax),&
         this%hy(:,:,nmin:nmax),&
         this%w(:,:,nmin:nmax), &
         prim(:,:,nmin:nmax,this%DENSITY), &
         prim(:,:,nmin:nmax,this%YVELOCITY),cons(:,:,nmin:nmax,this%XMOMENTUM), &
         cons(:,:,nmin:nmax,this%YMOMENTUM),yfluxes(:,:,nmin:nmax,this%DENSITY), &
         yfluxes(:,:,nmin:nmax,this%XMOMENTUM),yfluxes(:,:,nmin:nmax,this%YMOMENTUM))
  END SUBROUTINE CalcFluxesY_euler2Ditia


  PURE SUBROUTINE CalcCharSystemX_euler2Ditia(this,Mesh,i,dir,pvar,lambda,xvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: i,dir
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: lambda,xvar
    !------------------------------------------------------------------------!
    INTEGER           :: i1,i2
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,i,dir,pvar
    INTENT(INOUT)     :: lambda
    INTENT(OUT)       :: xvar
    !------------------------------------------------------------------------!
    ! compute eigenvalues at i
    CALL SetEigenValues_euler2Ditia(this%csiso,pvar(i,:,this%XVELOCITY), &
         lambda(:,1),lambda(:,2),lambda(:,3))
    ! compute characteristic variables
    i1 = i + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
    i2 = MAX(i,i1)
    i1 = MIN(i,i1)
    CALL SetCharVars_euler2Ditia(this%csiso,pvar(i1,:,this%DENSITY), &
         pvar(i2,:,this%DENSITY),pvar(i1,:,this%XVELOCITY), &
         pvar(i2,:,this%XVELOCITY),pvar(i1,:,this%YVELOCITY), &
         pvar(i2,:,this%YVELOCITY),xvar(:,1),xvar(:,2),xvar(:,3))
  END SUBROUTINE CalcCharSystemX_euler2Ditia


  PURE SUBROUTINE CalcCharSystemY_euler2Ditia(this,Mesh,j,dir,pvar,lambda,xvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: j,dir
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,this%VNUM) :: lambda,xvar
    !------------------------------------------------------------------------!
    INTEGER           :: j1,j2
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,j,dir,pvar
    INTENT(INOUT)     :: lambda
    INTENT(OUT)       :: xvar
    !------------------------------------------------------------------------!
    ! compute eigenvalues at j
    CALL SetEigenValues_euler2Ditia(this%csiso,&
         pvar(:,j,this%YVELOCITY)-Mesh%bhy(:,j)*this%w(:,j,1), &
         lambda(:,1),lambda(:,2),lambda(:,3))
    ! compute characteristic variables
    j1 = j + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
    j2 = MAX(j,j1)
    j1 = MIN(j,j1)
    CALL SetCharVars_euler2Ditia(this%csiso,pvar(:,j1,this%DENSITY), &
         pvar(:,j2,this%DENSITY),pvar(:,j1,this%YVELOCITY), &
         pvar(:,j2,this%YVELOCITY),pvar(:,j1,this%XVELOCITY), &
         pvar(:,j2,this%XVELOCITY),xvar(:,1),xvar(:,2),xvar(:,3))
  END SUBROUTINE CalcCharSystemY_euler2Ditia


  PURE SUBROUTINE CalcBoundaryDataX_euler2Ditia(this,Mesh,i1,dir,xvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: i1,dir
    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: xvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: i2
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,i1,dir,xvar
    INTENT(INOUT)     :: pvar
    !------------------------------------------------------------------------!
    i2 = i1 + SIGN(1,dir)  ! i +/- 1 depending on the sign of dir
    CALL SetBoundaryData_euler2Ditia(this%csiso,1.0*SIGN(1,dir), &
         pvar(i1,:,this%DENSITY),pvar(i1,:,this%XVELOCITY), &
         pvar(i1,:,this%YVELOCITY),xvar(:,1),xvar(:,2),xvar(:,3), &
         pvar(i2,:,this%DENSITY),pvar(i2,:,this%XVELOCITY),pvar(i2,:,this%YVELOCITY))
  END SUBROUTINE CalcBoundaryDataX_euler2Ditia


  PURE SUBROUTINE CalcBoundaryDataY_euler2Ditia(this,Mesh,j1,dir,xvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    INTEGER           :: j1,dir
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,this%VNUM) :: xvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: j2
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,j1,dir,xvar
    INTENT(INOUT)     :: pvar
    !------------------------------------------------------------------------!
    j2 = j1 + SIGN(1,dir)  ! j +/- 1 depending on the sign of dir
    CALL SetBoundaryData_euler2Ditia(this%csiso,1.0*SIGN(1,dir), &
         pvar(:,j1,this%DENSITY),pvar(:,j1,this%YVELOCITY), &
         pvar(:,j1,this%XVELOCITY),xvar(:,1),xvar(:,2),xvar(:,3), &
         pvar(:,j2,this%DENSITY),pvar(:,j2,this%YVELOCITY),pvar(:,j2,this%XVELOCITY))
  END SUBROUTINE CalcBoundaryDataY_euler2Ditia


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
          ! no geometrical density sources
          sterm(i,j,this%DENSITY) = 0.
          ! geometrical source terms in momentum equationes
          sterm(i,j,this%XMOMENTUM) = MomentumSourcesX_euler2Ditia(&
              cvar(i,j,this%DENSITY),&
              cvar(i,j,this%YMOMENTUM),&
              pvar(i,j,this%YVELOCITY), &
              pvar(i,j,this%DENSITY)*this%bccsound(i,j)**2, &
              Mesh%bhy(i,j), &
              Mesh%cyxy(i,j,1))

          ! no geometrical inertial angular momentum sources
          sterm(i,j,this%YMOMENTUM) = 0.0
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
          ! no geometrical density sources
          sterm(i,j,this%DENSITY)   = 0.
          ! momentum sources (sum up corner values, don't use SUM function,
          ! because it prevents COLLAPSING and causes poor vectorization
          sterm(i,j,this%XMOMENTUM) = MomentumSourcesX_euler2Ditia(&
              cons(i,j,1,this%DENSITY),&
              cons(i,j,1,this%YMOMENTUM),&
              prim(i,j,1,this%YVELOCITY), &
              prim(i,j,1,this%DENSITY)*this%fcsound(i,j,1)**2, &
              this%hy(i,j,1), &
              Mesh%cyxy(i,j,1)) &
            + MomentumSourcesX_euler2Ditia(&
              cons(i,j,2,this%DENSITY),&
              cons(i,j,2,this%YMOMENTUM),&
              prim(i,j,2,this%YVELOCITY), &
              prim(i,j,1,this%DENSITY)*this%fcsound(i,j,3)**2, &
              this%hy(i,j,2), &
              Mesh%cyxy(i,j,2)) &
            + MomentumSourcesX_euler2Ditia(&
              cons(i,j,3,this%DENSITY),&
              cons(i,j,3,this%YMOMENTUM),&
              prim(i,j,3,this%YVELOCITY), &
              prim(i,j,1,this%DENSITY)*this%fcsound(i,j,3)**2, &
              this%hy(i,j,3), &
              Mesh%cyxy(i,j,3)) &
            + MomentumSourcesX_euler2Ditia(&
              cons(i,j,4,this%DENSITY),&
              cons(i,j,4,this%YMOMENTUM),&
              prim(i,j,4,this%YVELOCITY), &
              prim(i,j,1,this%DENSITY)*this%fcsound(i,j,4)**2, &
              this%hy(i,j,4), &
              Mesh%cyxy(i,j,4))

          ! no geometrical inertial angular momentum sources
          sterm(i,j,this%YMOMENTUM) = 0.
       END DO
    END DO
  END SUBROUTINE GeometricalSources_faces


  ! momentum and energy sources due to external force
  PURE SUBROUTINE ExternalSources_euler2Ditia(this,Mesh,accel,pvar,cvar,sterm)
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
!CDIR IEXPAND
    CALL ExternalSources_euler2Dit(this,Mesh,accel,pvar,cvar,sterm)
    sterm(:,:,this%YMOMENTUM) = sterm(:,:,this%YMOMENTUM) * Mesh%bhy(:,:)
  END SUBROUTINE ExternalSources_euler2Ditia


  PURE SUBROUTINE ViscositySources_euler2Ditia(this,Mesh,pvar,btxx,btxy,btyy,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: &
         pvar,sterm
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
         btxx,btxy,btyy
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,pvar,btxx,btxy,btyy
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL ViscositySources_euler2Dit(this,Mesh,pvar,btxx,btxy,btyy,sterm)
    sterm(:,:,this%YMOMENTUM) = sterm(:,:,this%YMOMENTUM) * Mesh%bhy(:,:)
  END SUBROUTINE ViscositySources_euler2Ditia


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
!CDIR IEXPAND
    CALL Cons2Prim_euler2Ditia(&
         this%Omega, &
         Mesh%bhy(i1:i2,j1:j2), &
         cvar(i1:i2,j1:j2,this%DENSITY), &
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
!CDIR IEXPAND
    CALL Cons2Prim_euler2Ditia(&
         this%Omega, &
         this%hy(i1:i2,j1:j2,:), &
         cons(i1:i2,j1:j2,:,this%DENSITY), &
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
!CDIR IEXPAND
    CALL Prim2Cons_euler2Ditia(&
         this%Omega, &
         Mesh%bhy(i1:i2,j1:j2), &
         pvar(i1:i2,j1:j2,this%DENSITY), &
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
!CDIR IEXPAND
    CALL Prim2Cons_euler2Ditia(&
         this%Omega, &
         this%hy(i1:i2,j1:j2,:), &
         prim(i1:i2,j1:j2,:,this%DENSITY), &
         prim(i1:i2,j1:j2,:,this%XVELOCITY),prim(i1:i2,j1:j2,:,this%YVELOCITY), &
         cons(i1:i2,j1:j2,:,this%DENSITY),cons(i1:i2,j1:j2,:,this%XMOMENTUM), &
         cons(i1:i2,j1:j2,:,this%YMOMENTUM))
  END SUBROUTINE Convert2Conservative_faces


  ELEMENTAL SUBROUTINE CalcFluxY_euler2Ditia(cs,hy,w,rho,vy,mx,my,g1,g2,g3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: hy,cs,rho,vy,mx,my,w
    REAL, INTENT(OUT) :: g1, g2, g3
    !------------------------------------------------------------------------!
    REAL              :: v
    !------------------------------------------------------------------------!
    v  = vy - hy*w
    g1 = rho * v
    g2 = mx  * v
    g3 = my  * v + hy*rho*cs*cs
  END SUBROUTINE CalcFluxY_euler2Ditia

  ! momentum source terms due to inertial forces
  ! P is the isothermal pressure rho*cs*cs
  ELEMENTAL FUNCTION MomentumSourcesX_euler2Ditia(rho,my,vy,P,hy,cyxy) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: my,vy,P,hy,cyxy,rho
    REAL :: st
    !------------------------------------------------------------------------!
    st = my*my / rho * cyxy**3 + P * cyxy
    !st = cyxy * (my*vy/hy + P)
  END FUNCTION MomentumSourcesX_euler2Ditia

  ELEMENTAL SUBROUTINE Cons2Prim_euler2Ditia(Omega,hy,rho_in,mu,mv,rho_out,u,v)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho_in,mu,mv,hy,Omega
    REAL, INTENT(OUT) :: rho_out,u,v
    !------------------------------------------------------------------------!
    REAL :: inv_rho
    !------------------------------------------------------------------------!
    inv_rho = 1./rho_in
    rho_out = rho_in
    u = mu * inv_rho
    v = mv * inv_rho / hy - hy * Omega
  END SUBROUTINE Cons2Prim_euler2Ditia

  
  ELEMENTAL SUBROUTINE Prim2Cons_euler2Ditia(Omega,hy,rho_in,u,v,rho_out,mu,mv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho_in,u,v,hy,Omega
    REAL, INTENT(OUT) :: rho_out,mu,mv
    !------------------------------------------------------------------------!
    rho_out = rho_in
    mu = rho_in * u
    mv = rho_in * hy * (v + hy*Omega)
  END SUBROUTINE Prim2Cons_euler2Ditia
  

  SUBROUTINE ClosePhysics_euler2Ditia(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%w)
    CALL ClosePhysics(this)
  END SUBROUTINE ClosePhysics_euler2Ditia

END MODULE physics_euler2Disoiamt
