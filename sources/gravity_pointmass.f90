!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: gravity_pointmass.f90                                             #
!#                                                                           #
!# Copyright (C) 2007-2014                                                   #
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
!> \addtogroup gravity
!! - parameters of \link gravity_pointmass \endlink as key-values
!! \key{potential,INTEGER,type of the potential}
!! \key{mass,REAL,mass of the point mass, 1.0}
!! \key{x,REAL,cartesian x-position of the point mass,0.0}
!! \key{y,REAL,cartesian y-position of the point mass,0.0}
!! \key{softening,REAL,Softening (e.g. for planets inside the computational domain),0.0}
!! \key{switchon,REAL,soft switch on,-1.0}
!! \key{outbound,INTEGER,enable mass accretion by setting "outbound" to one of the four boundaries,depends on mesh}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!!
!! \brief source terms module for gravitational acceleration due to a point
!! mass at the center of the coordinate system
!!
!! \extends gravity_common
!! \ingroup gravity
!----------------------------------------------------------------------------!
MODULE gravity_pointmass
  USE gravity_common
  USE common_types, ONLY : Common_TYP, InitCommon
  USE boundary_common, ONLY : WEST,EAST,SOUTH,NORTH
  USE fluxes_generic, ONLY : Fluxes_TYP, GetBoundaryFlux
  USE physics_generic
  USE geometry_generic
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
  CHARACTER(LEN=32), PARAMETER :: gravity_name = "central point mass"
  CHARACTER(LEN=16), PARAMETER :: potential_name(2) = (/ &
                                  "Newton          ", &
                                  "Paczinski-Wiita " /)
  INTEGER, PARAMETER :: NEWTON = 1
  INTEGER, PARAMETER :: WIITA  = 2
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Gravity_TYP, &
       ! constants
       NEWTON, WIITA, &
       ! methods
       InitGravity_pointmass, &
       InfoGravity_pointmass, &
       CalcDiskHeight_pointmass, &
       GetAccelGravity_pointmass, &
       GetDiskHeight_pointmass, &
       GetinvDistanceCO_pointmass, &
       GetDistVector_pointmass,&
       CloseGravity_pointmass, &
       GetGravityPointer, &
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

  SUBROUTINE InitGravity_pointmass(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Dict_TYP),POINTER :: config,IO
    INTEGER           :: gtype
    !------------------------------------------------------------------------!
    INTEGER           :: potential_def, valwrite
    INTEGER           :: err
    INTEGER           :: i,j
    REAL              :: r,a=0.0,eps
    REAL              :: invdt_x, invdt_y
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM) :: accel
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "gtype", gtype)
    CALL InitGravity(this,gtype,gravity_name)

    ! type of the potential
    CALL RequireKey(config, "potential", NEWTON)
    CALL GetAttr(config, "potential", potential_def)

    ! set mass
    CALL RequireKey(config, "mass", 1.0)
    CALL GetAttr(config, "mass", this%mass)

    ! cartesian position of the point mass
    CALL RequireKey(config, "x", 0.0)
    CALL RequireKey(config, "y", 0.0)
    CALL GetAttr(config, "x", this%r0(1))
    CALL GetAttr(config, "y", this%r0(2))

    ! Softening (e.g. for planets inside the computational domain)
    CALL RequireKey(config, "softening", 0.0)
    CALL GetAttr(config, "softening", eps)

    ! soft switch on
    CALL RequireKey(config, "switchon", -1.0)
    CALL GetAttr(config, "switchon", this%switchon)
    
    ! reset mass flux
    this%mdot = 0.0

    this%cs => GetSoundSpeeds(Physics)

    ALLOCATE(this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM), &
             this%omega(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
             this%height(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
             this%invr(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
             this%r_prim(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2),&
             this%scaled(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM),&
         STAT = err)
    IF (err.NE.0) CALL Error(this,"InitGravity_pointmass", "Unable allocate memory!")

    ! initialize potential
    SELECT CASE(potential_def)
    CASE(NEWTON) ! newtonian gravity
       CALL InitCommon(this%potential,NEWTON,potential_name(NEWTON))
       a = 0.0
    CASE(WIITA) ! pseudo-Newton Paczinski-Wiita potential
       CALL InitCommon(this%potential,WIITA,potential_name(WIITA))
       ! Schwarzschild radius
       a = 2.*Physics%constants%GN * this%mass / Physics%constants%C**2
    CASE DEFAULT
       CALL Error(this,"InitGravity_pointmass", "potential must be either NEWTON or WIITA")       
    END SELECT

    ! define position vector and inverse radius
    IF (ABS(this%r0(1)).LE.TINY(this%r0(1)).AND.ABS(this%r0(2)).LE.TINY(this%r0(2))) THEN
       ! no shift of point mass: set position vector and inverse radius to Mesh defaults
       this%r_prim(:,:,:) = Mesh%bposvec(:,:,:)
       IF(eps.GT.0.) THEN
         this%invr(:,:) = 1.0 / SQRT(Mesh%bradius(:,:)**2+eps*eps)
       ELSE
         this%invr(:,:) = 1.0 / Mesh%bradius(:,:)
       END IF
    ELSE
       ! shifted point mass position:
       ! compute curvilinear components of shift vector
       this%r_prim(:,:,1) = this%r0(1)
       this%r_prim(:,:,2) = this%r0(2)
       CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,this%r_prim,this%r_prim)
       ! subtract the result from the position vector:
       ! this gives you the curvilinear components of all vectors pointing 
       ! from the point mass to the bary center of any cell on the mesh
       this%r_prim(:,:,:) = Mesh%bposvec(:,:,:) - this%r_prim(:,:,:)
       ! compute the inverse of its absolute value
       IF(eps.GT.0.) THEN
         this%invr(:,:) = 1.0 / SQRT(this%r_prim(:,:,1)**2+this%r_prim(:,:,2)**2+eps*eps)
       ELSE
         this%invr(:,:) = 1.0 / SQRT(this%r_prim(:,:,1)**2+this%r_prim(:,:,2)**2)
       END IF
    END IF

    ! initialize gravitational acceleration an Keplerian angular velocity
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! compute Keplerian angular velocity
          ! omega = sqrt(GM/r*(r-a)^2) = sqrt(GM/r) / (r-a) = sqrt(GM/r) / r / (1-a/r)
          ! the last term vanishes for newtonian gravity (a=0)
          this%omega(i,j) = SQRT(Physics%Constants%GN*this%mass*this%invr(i,j)) &
                           * this%invr(i,j) * 1.0 / (1.0-a*this%invr(i,j))
          ! curvilinear components of the gravitational acceleration
          ! -d Phi / dr = -r * omega^2 * e_r = -omega^2 * r_prim
          this%accel(i,j,1:2) = -this%omega(i,j)*this%omega(i,j) * this%r_prim(i,j,1:2)
       END DO
    END DO

    ! enable mass accretion by setting "outbound" to one of the four boundaries
    ! of the computational domain (depends on mesh geometry)
    SELECT CASE(GetType(Mesh%geometry))
    CASE(POLAR,LOGPOLAR,TANPOLAR,SINHPOLAR,SPHERICAL,OBLATE_SPHEROIDAL,SINHSPHERICAL)
       CALL RequireKey(config, "outbound", WEST)
    CASE(CYLINDRICAL,TANCYLINDRICAL)
       CALL RequireKey(config, "outbound", SOUTH)
    CASE DEFAULT
       CALL RequireKey(config, "outbound", 0)! disable growth of central point mass
       CALL WARNING(this,"GravitySources_pointmass","geometry does not support accretion")
    END SELECT
    CALL GetAttr(config, "outbound", this%outbound)
  END SUBROUTINE InitGravity_pointmass


  SUBROUTINE InfoGravity_pointmass(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: mass_str
    !------------------------------------------------------------------------!
    WRITE (mass_str,'(ES8.2)') this%mass
    CALL Info(this,"            potential:         " // &
         TRIM(GetName(this%potential)))
    CALL Info(this,"            mass:              " // TRIM(mass_str))
  END SUBROUTINE InfoGravity_pointmass

  SUBROUTINE UpdatePointmass(this,Mesh,Physics,Fluxes,time,pvar) 
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    !------------------------------------------------------------------------!
    REAL, DIMENSION(Physics%VNUM) :: bflux
    REAL              :: oldmass,r
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,time,pvar
    INTENT(INOUT)     :: Physics
    !------------------------------------------------------------------------!
    ! update accel and omega only in case of accretion
    IF (this%outbound.NE.0) THEN
        ! get boundary flux
#ifdef PARALLEL
        bflux(:)  = GetBoundaryFlux(Fluxes,Mesh,Physics,this%outbound,MPI_COMM_WORLD)
#else
!CDIR IEXPAND
        bflux(:)  = GetBoundaryFlux(Fluxes,Mesh,Physics,this%outbound)
#endif
        ! store old mass
        oldmass = this%mass
        ! compute new mass
        this%mass = this%mass + (bflux(Physics%DENSITY) - this%mdot)
        ! scale gravitational acceleration and Keplerian angular velocity
        ! with newmass/oldmass to account for accretion
        this%accel(:,:,:) = this%mass/oldmass * this%accel(:,:,:)
        this%omega(:,:)   = SQRT(this%mass/oldmass) * this%omega(:,:)
        ! store actual total mass flux
        this%mdot = bflux(Physics%DENSITY)
      END IF
  END SUBROUTINE UpdatePointmass

 PURE SUBROUTINE UpdateDiskHeight_pointmass(this,Mesh,Physics,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    REAL              :: cs,sqrtgamma
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,pvar
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! compute disk height
    sqrtgamma = SQRT(Physics%gamma)
    DO j=Mesh%JGMIN,Mesh%JGMAX
      DO i=Mesh%IGMIN,Mesh%IGMAX
        this%height(i,j) &
          = CalcDiskHeight_pointmass(sqrtgamma,this%cs(i,j),this%omega(i,j))
      END DO
    END DO
  END SUBROUTINE UpdateDiskHeight_pointmass

  ELEMENTAL FUNCTION CalcDiskHeight_pointmass(sqrtgamma,cs,omega) RESULT(height)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: sqrtgamma,cs,omega
    REAL             :: height
    !------------------------------------------------------------------------!
     height = cs / (sqrtgamma * omega)
  END FUNCTION CalcDiskHeight_pointmass

  FUNCTION GetAccelGravity_pointmass(this,Mesh,Physics,Fluxes,time,pvar) RESULT(ac)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    REAL, DIMENSION(:,:,:), POINTER &
                      :: ac
    !------------------------------------------------------------------------!
    REAL              :: switch
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,pvar,time
    INTENT(INOUT)     :: this,Physics
    !------------------------------------------------------------------------!
    CALL UpdatePointmass(this,Mesh,Physics,Fluxes,time,pvar)

    IF (time.LE.this%switchon) THEN
      switch = SIN(0.5*PI*time/this%switchon)**2
      this%scaled(:,:,:) = this%accel(:,:,:) * switch
      ac => this%scaled
    ELSE
      ac => this%accel
    END IF

  END FUNCTION GetAccelGravity_pointmass

  FUNCTION GetDiskHeight_pointmass(this,Mesh,Physics,pvar) RESULT(height)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                      :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: height
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar
    INTENT(INOUT)     :: Physics
    !------------------------------------------------------------------------!
    ! Update disk height
    CALL UpdateDiskHeight_pointmass(this,Mesh,Physics,pvar)
!FIXME? Pointer?
    height(:,:) = this%height(:,:)
  END FUNCTION GetDiskHeight_pointmass

 FUNCTION GetinvDistanceCO_pointmass(this,Mesh) RESULT(invdis)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: invdis
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    !------------------------------------------------------------------------!
!FIXME?: Pointer?
    invdis(:,:) = this%invr(:,:)
   END FUNCTION GetinvDistanceCO_pointmass

! returns 2D vector from point mass position to cell bary centers
! in curvlinear coordinates
  FUNCTION GetDistVector_pointmass(this,Mesh) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
                      :: r
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    !------------------------------------------------------------------------!
!FIXME?: Pointer?
    r(:,:,:) = Mesh%bposvec(:,:,:)
  END FUNCTION GetDistVector_pointmass

  SUBROUTINE CloseGravity_pointmass(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP) :: this
    !------------------------------------------------------------------------!
    CHARACTER(LEN=128):: buffer
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    IF (GetRank(this).EQ.0) THEN
      CALL Info(this,"-------------------------------------------------------------------")
       WRITE(buffer, "(A,(ES11.3))") " central mass: ", this%mass
       CALL Info(this,buffer)
    END IF
    DEALLOCATE(this%accel,this%invr,this%r_prim,this%scaled)
    CALL CloseGravity(this)
  END SUBROUTINE CloseGravity_pointmass

END MODULE gravity_pointmass
