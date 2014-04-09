!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_rotframe.f90                                              #
!#                                                                           #
!# Copyright (C) 2010-2011                                                   #
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
!> \author Tobias Illenseer
!!
!! \brief source terms module for inertial forces caused by a rotating grid
!!
!! \extends sources_c_accel
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_rotframe
  USE sources_c_accel
  USE fluxes_common, ONLY : Fluxes_TYP
  USE physics_generic
  USE mesh_generic
  USE common_dict
#ifdef PARALLEL
#ifdef HAVE_MPI_MOD
  USE mpi
#endif
#endif
  IMPLICIT NONE
#ifdef PARALLEL
#ifdef HAVE_MPIF_H
  include 'mpif.h'
#endif
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: source_name = "inertial forces"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! methods
       InitSources_rotframe, &
       InfoSources_rotframe, &
       ExternalSources_rotframe, &
       Convert2RotatingFrame_rotframe, &
       CloseSources_rotframe
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources_rotframe(this,Mesh,Physics,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Physics_TYP) :: Physics
    TYPE(Dict_TYP),POINTER :: config
    INTEGER           :: stype
    REAL              :: omega,gparam
    !------------------------------------------------------------------------!
    INTEGER           :: err
    INTEGER           :: i,j
    REAL              :: x, y
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "stype", stype)
    CALL InitSources(this,stype,source_name)

    SELECT CASE(GetType(Physics))
    CASE(EULER2D,EULER2D_ISOTHERM,EULER2D_IAMT)
       ! do nothing
    CASE DEFAULT
       CALL Error(this,"ExternalSources_rotframe","physics not supported")
    END SELECT

    ALLOCATE(this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM), &
         this%cent(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
         this%centproj(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
         this%cos1(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
         this%omega(1,1), &
         STAT = err)
    IF (err.NE.0) CALL Error(this,"InitSources_rotframe", "Unable allocate memory!")

    ! set angular velocity of the rotating frame
    CALL RequireKey(config, "omega", 0.0)
    CALL GetAttr(config, "omega", omega)
    this%omega(1,1) = omega
    
    CALL RequireKey(config, "gparam", 1.0)
    CALL GetAttr(config, "gparam", gparam)

    ! shifting values
    CALL RequireKey(config, "x", 0.0)
    CALL RequireKey(config, "y", 0.0)
    CALL GetAttr(config, "x", x)
    CALL GetAttr(config, "y", y)
    ! define position vectors
    ! for bianglespherical no shift of axis possible 
    ! (for a planet reasonable)
    IF (GetType(Mesh%geometry).EQ.BIANGLESPHERICAL) THEN

      this%centproj(:,:,1) = gparam*SIN(Mesh%bcenter(:,:,1))*&
                              COS(Mesh%bcenter(:,:,1))
      ! for better performance
      this%cos1(:,:,1)  = COS(Mesh%bcenter(:,:,1))
      this%cos1(:,:,2)  = COS(Mesh%bcenter(:,:,2))
    ELSE
      IF (ABS(x).LE.TINY(x).AND.ABS(y).LE.TINY(y)) THEN
        ! no shift of point mass: set position vector to Mesh defaults
        this%cent(:,:,:) = Mesh%bposvec(:,:,:)
      ELSE
        ! shifted center of rotation:
        ! compute curvilinear components of shift vector
        this%cent(:,:,1) = x
        this%cent(:,:,2) = y
        CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,this%cent,this%cent)
        ! subtract the result from the position vector:
        ! this gives you the curvilinear components of all vectors pointing 
        ! from the center of rotation to the bary center of any cell on the mesh
        this%cent(:,:,:) = Mesh%bposvec(:,:,:) - this%cent(:,:,:)
      END IF
    END IF

    ! reset acceleration term
    this%accel(:,:,:) = 0.0
  END SUBROUTINE InitSources_rotframe


  SUBROUTINE InfoSources_rotframe(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: omega_str
    !------------------------------------------------------------------------!
    WRITE (omega_str,'(ES8.2)') this%omega(1,1)
    CALL Info(this,"            angular velocity:  " // TRIM(omega_str))
  END SUBROUTINE InfoSources_rotframe


  SUBROUTINE ExternalSources_rotframe(this,Mesh,Physics,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics

    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: cvar,pvar,sterm
    REAL              :: gparam
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,cvar,pvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! Two cases due to different angles between angular velocity to mesh
    ! 1. only a projected part plays role for bianglespherical geometry
    IF (GetType(Mesh%geometry).EQ.BIANGLESPHERICAL) THEN
!CDIR OUTERUNROLL=8
      DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
         DO i=Mesh%IMIN,Mesh%IMAX
            this%accel(i,j,1) = this%omega(1,1)*(this%omega(1,1)*&
                   this%centproj(i,j,1) + 2.0*this%cos1(i,j,1)*&
                   pvar(i,j,Physics%YVELOCITY))
            this%accel(i,j,2) = -this%omega(1,1)*2.0*this%cos1(i,j,1)*&
                   pvar(i,j,Physics%XVELOCITY)
         END DO
      END DO
    ! 2. omega is always perpendicular to other curvilinear coordinates
    ELSE
!CDIR OUTERUNROLL=8
      DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
         DO i=Mesh%IMIN,Mesh%IMAX
            ! components of centrifugal and coriolis acceleration
!!$            this%accel(i,j,1) = this%omega(1,1) * 2*pvar(i,j,Physics%YVELOCITY)
!!$            this%accel(i,j,2) = -this%omega(1,1) * 2*pvar(i,j,Physics%XVELOCITY)
!!$            this%accel(i,j,:) = this%omega(1,1) * this%omega(1,1)*this%cent(i,j,:)
            this%accel(i,j,1) = this%omega(1,1)*(this%omega(1,1)*this%cent(i,j,1) &
                 + 2.0*pvar(i,j,Physics%YVELOCITY)) 
            this%accel(i,j,2) = this%omega(1,1)*(this%omega(1,1)*this%cent(i,j,2) &
                 - 2.0*pvar(i,j,Physics%XVELOCITY)) 
         END DO
      END DO
    END IF
!!$          ! components of centrifugal and coriolis acceleration
!!$          PRINT '(6(ES28.20))',Mesh%bcenter(i,j,1), Mesh%bcenter(i,j,2), &
!!$               this%omega(1,1)**2 * this%cent(i,j,1),this%omega(1,1)**2 &
!!$               * this%cent(i,j,2),this%accel(i,j,1),this%accel(i,j,2)
!!$       END DO
!!$       PRINT *,""
!!$    END DO
!!$    STOP
    ! inertial forces source terms
    CALL ExternalSources(Physics,Mesh,this%accel,pvar,cvar,sterm)
  END SUBROUTINE ExternalSources_rotframe

  SUBROUTINE Convert2RotatingFrame_rotframe(this,Mesh,Physics,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    REAL              :: gparam
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,this
    INTENT(INOUT)     :: pvar
    !------------------------------------------------------------------------!
    ! Convert velocities to the rotating frame
    IF (GetType(Mesh%geometry).EQ.BIANGLESPHERICAL) THEN     
      pvar(:,:,Physics%XVELOCITY) = &
      pvar(:,:,Physics%XVELOCITY)
      pvar(:,:,Physics%YVELOCITY) = &
      pvar(:,:,Physics%YVELOCITY) &
      - this%omega(1,1)*SIN(Mesh%bcenter(:,:,1))*gparam
    ELSE
      pvar(:,:,Physics%XVELOCITY) = &
      pvar(:,:,Physics%XVELOCITY) &
      + this%omega(1,1) * this%cent(:,:,2)
      pvar(:,:,Physics%YVELOCITY) = &
      pvar(:,:,Physics%YVELOCITY) &
      - this%omega(1,1) * this%cent(:,:,1)
    END IF
  END SUBROUTINE Convert2RotatingFrame_rotframe
 
  SUBROUTINE CloseSources_rotframe(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%accel,this%cent,this%omega,this%cos1,this%centproj)
    CALL CloseSources(this)
  END SUBROUTINE CloseSources_rotframe

END MODULE sources_rotframe
