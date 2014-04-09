!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: gravity_generic.f90                                               #
!#                                                                           #
!# Copyright (C) 2014                                                        #
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
!! - general parameters of gravity group as key-values
!! \key{gtype,INTEGER,Type of gravity source}
!! \key{energy,INTEGER,Add source terms to energy equation?}
!! \key{output/accel,INTEGER,enable(=1) output of acceleration}
!! \key{output/height,INTEGER,enable(=1) output of disc height}
!----------------------------------------------------------------------------!
!> \author Björn Sperling
!!
!! \brief generic gravity terms module providing functionaly common to all
!! gravity terms
!!
!! \ingroup gravity
!----------------------------------------------------------------------------!
MODULE gravity_generic
  USE mesh_common, ONLY : Mesh_TYP
  USE boundary_common, ONLY : Boundary_TYP
  USE timedisc_common, ONLY : Timedisc_TYP
  USE sources_common, ONLY: Sources_TYP, InitSources, Error
  USE gravity_common, ONLY: Gravity_TYP
  USE gravity_pointmass
  USE gravity_binary
  USE gravity_monopol
  USE gravity_multigrid
  USE gravity_spectral
  USE gravity_potential
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic, ONLY : PI
  USE geometry_generic
  USE common_dict
  USE boundary_generic
  USE mesh_generic, ONLY : PI, Divergence
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  USE boundary_common, ONLY : Boundary_TYP
  USE physics_generic
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: source_name = "gravity"
  ! flags for source terms
  INTEGER, PARAMETER :: POINTMASS        = 1
  INTEGER, PARAMETER :: POINTMASS_BINARY = 2
  INTEGER, PARAMETER :: MONOPOL          = 3
  INTEGER, PARAMETER :: MULTIGRID        = 4
  INTEGER, PARAMETER :: SPECTRAL         = 5
  INTEGER, PARAMETER :: POTENTIAL        = 6
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! constants
       POINTMASS, POINTMASS_BINARY, MONOPOL, &
       NEWTON, WIITA, &
       MULTIGRID, SPECTRAL, POTENTIAL, &
       RED_BLACK_GAUSS_SEIDEL,BLOCK_GAUSS_SEIDEL,GAUSS_SEIDEL, &
       SPHERMULTEXPAN, CYLINMULTEXPAN, &
       ! methods
       InitGravity, &
       CloseGravity, &
       GravitySources, &
       GetGravityPointer, &
       GetDiskHeight, &
       GetDiskinvHeight2, &
       GetAccelGravity, &
       GetinvDistanceCO, &
       GetDistVector,&
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

  SUBROUTINE InitGravity(this,Mesh,Fluxes,Physics,Boundary,stype,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Physics_TYP) :: Physics
    TYPE(Boundary_TYP), DIMENSION(4) :: Boundary
    TYPE(Dict_TYP),POINTER :: config,IO
    INTEGER           :: stype
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: gp
    TYPE(Dict_TYP),POINTER :: dir,src,IOsrc
    INTEGER           :: gtype,err,k,i
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,Physics,stype
    INTENT(INOUT)     :: Boundary
    !------------------------------------------------------------------------!
    ! gravity is a source module
    CALL InitSources(this,stype,source_name)

    ALLOCATE(this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM), &
             this%invheight2(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
             this%height(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
             this%tmp(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
             this%tmp2(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
             this%tmp4(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM), &
             this%bcposvec(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
             this%bccart(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
             this%calc_boundary(4),&
         STAT=err)
    IF (err.NE.0) CALL Error(this, "InitGravity", "Unable allocate memory!")

    this%accel(:,:,:) = 0.

    ! Add source terms to energy equation?
    ! Set this to zero, if a potential is defined in physics_euler2Diamt
    CALL RequireKey(config, "energy", 1)
    CALL GetAttr(config, "energy", i)
    IF(i.EQ.0) THEN
      this%addtoenergy = .FALSE.
    ELSE
      this%addtoenergy = .TRUE.
    END IF

    dir => config
    DO WHILE(ASSOCIATED(dir))
      IF ((GetDataType(dir).eq.DICT_DIR) .AND. (TRIM(GetKey(dir)).NE."output")) THEN
        CALL GetAttr(config, GetKey(dir), src)
        CALL GetAttr(IO, GetKey(dir), IOsrc)

        CALL RequireKey(src, "gtype")
        CALL GetAttr(src, "gtype", gtype)

        SELECT CASE(gtype)
        CASE(POINTMASS)
           ! gravitational acceleration due to point mass
           CALL InitGravity_pointmass(this%glist,Mesh,Physics,src,IOsrc)
        CASE(MONOPOL)
           ! monopol approximation for gravity
           CALL InitGravity_monopol(this%glist,Mesh,Physics,src)
        CASE(POINTMASS_BINARY)
           ! gravitational source term of rotating binary system
           CALL InitGravity_binary(this%glist,Mesh,Physics,src)
        CASE(MULTIGRID)
           ! self-gravitation in geometries with rotational symmetry
           CALL InitGravity_multigrid(this%glist,Mesh,Physics,Boundary,src)
        CASE(SPECTRAL)
           ! self-gravitation in flat geometries
           CALL InitGravity_spectral(this%glist,Mesh,Physics,Boundary,src,IOsrc)
        CASE(POTENTIAL)
           ! gravitation due to one or multiple constant potentials
           CALL InitGravity_potential(this%glist,Mesh,Physics,Boundary,src,IOsrc)
        CASE DEFAULT
           CALL Error(this%glist,"InitGravity", "unknown gravity term")
        END SELECT
        ! print some information
        IF (ASSOCIATED(this%glist)) THEN
           CALL Info(this%glist, " GRAVITY--> Gravity term:      " // GetName(this%glist))
           ! print setup information of the individual Gravity terms
           CALL InfoGravity(this%glist)
        END IF
      END IF
      dir => GetNext(dir)
    END DO

!FIXME: good solution??????
    ! set accel pointer in boundary module
    Boundary(1)%accel => this%accel
    Boundary(2)%accel => this%accel

! FIXME: this only works for pointmass or binary 
! The boundary values must not be calculated if the boundaries are either periodic or 
! belong to mpi halos in the parallel version of the code (boundary=NONE).
! In these 2 cases the height is already calculated from the correct (periodic) boundary values
! of the pvars. 
! Only the case of non-periodic outer boundaries has to be treated seperately as the height 
! calculated for artificial bourdary pvars may be meaningless.
! The calculation does not change if we go from serial to parallel execution. 
    this%calc_boundary(:) = 0
#ifdef PARALLEL
    DO k=1,4
       IF(.not.(GetType(boundary(k)).EQ.PERIODIC .or. GetType(boundary(k)).EQ.NONE)) &
          this%calc_boundary(k) = 1
    END DO    
#else
    DO k=1,4
       IF(.not.(GetType(boundary(k)).EQ.PERIODIC )) &
          this%calc_boundary(k) = 1
    END DO       
#endif
    CALL SetOutput(this,config,IO)
    ! reset start value for time variable
    this%time = -1.0
    this%timeid = 0
  END SUBROUTINE InitGravity
  

  SUBROUTINE SetOutput(this,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP)       :: this
    TYPE(Dict_TYP),POINTER  :: config,IO
    !------------------------------------------------------------------------!
    INTEGER                 :: valwrite
    !------------------------------------------------------------------------! 
    valwrite = 0
    IF (HasKey(config, "output/accel")) CALL GetAttr(config, "output/accel", valwrite)
    IF (valwrite .EQ. 1) THEN
       CALL AddField(IO, &
               "accel", &
               this%accel, &
               Dict("name" / "accel"))
    END IF

    valwrite = 0
    IF (HasKey(config, "output/height")) CALL GetAttr(config, "output/height", valwrite)
    IF (valwrite .EQ. 1) THEN
       CALL AddField(IO, &
               "height", &
               this%height, &
               Dict("name" / "height"))
    END IF

  END SUBROUTINE SetOutput

  SUBROUTINE InfoGravity(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this
    !------------------------------------------------------------------------!
    IF (ASSOCIATED(this)) THEN
!CDIR IEXPAND
       SELECT CASE(GetType(this))
       CASE(MONOPOL)
          ! do nothing
       CASE(POINTMASS)
          CALL InfoGravity_pointmass(this)
       CASE(POINTMASS_BINARY)
          CALL InfoGravity_binary(this)
       END SELECT
    END IF
  END SUBROUTINE InfoGravity

  SUBROUTINE GravitySources(this,Mesh,Physics,Fluxes,time,dt,pvar,cvar,gterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL              :: time,dt
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: cvar,pvar,gterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,time,pvar,cvar
    INTENT(INOUT)     :: Physics
    INTENT(OUT)       :: gterm
    !------------------------------------------------------------------------!
    ! update acceleration of all gravity sources
    CALL CalcAccelGravity(this,Mesh,Physics,Fluxes,time,pvar)

    ! gravitational source terms
    CALL ExternalSources(Physics,Mesh,this%accel,pvar,cvar,gterm)

    ! Set src term in energy equation to zero, if it is handeled in the physics
    ! module
    IF((.NOT.this%addtoenergy).AND.(Physics%ENERGY.GT.0)) THEN
      gterm(:,:,Physics%ENERGY) = 0.
    END IF
  END SUBROUTINE GravitySources

  SUBROUTINE CalcDiskinvHeight2(this,Mesh,Physics,Fluxes,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: gravptr
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,time,pvar
    INTENT(INOUT)     :: Physics
    !------------------------------------------------------------------------!
    ! update acceleration of all gravity sources
    CALL CalcAccelGravity(this,Mesh,Physics,Fluxes,time,pvar)

    IF (Update(this,time,1)) THEN
      ! make sure sound speed is up to date
      CALL UpdateSoundSpeed(Physics,Mesh,time,pvar)
      this%invheight2(:,:) = 0.0
      ! go through all gravity terms in the list
      gravptr => this%glist
      DO WHILE(ASSOCIATED(gravptr))
        ! call specific subroutine
!CDIR IEXPAND
       SELECT CASE(GetType(gravptr))
        CASE(POINTMASS)
          this%tmp(:,:) = GetDiskHeight_pointmass(gravptr,Mesh,Physics,pvar)
        CASE(POINTMASS_BINARY)
          this%tmp(:,:) = GetDiskHeight_binary(gravptr,Mesh,Physics,time,pvar)
        CASE DEFAULT
          CALL Error(gravptr,"CalcDiskinvHeight2", "unknown gravity term")
        END SELECT
        this%invheight2(:,:) = this%invheight2(:,:) + 1./this%tmp(:,:)**2
        ! next source term
        gravptr => gravptr%next
      END DO    
    END IF
  END SUBROUTINE CalcDiskinvHeight2

  SUBROUTINE CalcBoundary(this,Mesh,var)
      IMPLICIT NONE
      !-----------------------------------------------------------------------!
      TYPE(Sources_TYP),POINTER :: this
      TYPE(Mesh_TYP)    :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX):: var
      !-----------------------------------------------------------------------!
      INTEGER ::  i,j
      !-----------------------------------------------------------------------!     
      INTENT(IN)  :: Mesh
      INTENT(INOUT) :: var 
      !-----------------------------------------------------------------------!     
      
      ! UNROLL=Mesh%GNUM would be sufficient, but the compiler does
       ! not know the value of Mesh%GNUM, hence we set UNROLL=4 and
       ! hope that nobody sets Mesh%GNUM to a value greater than 4
      IF(this%calc_boundary(WEST).EQ.1)THEN
!CDIR UNROLL=4
         DO i=1,Mesh%GNUM
            var(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX) = (i+1)*var(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX) &
                 - i*var(Mesh%IMIN+1,Mesh%JMIN:Mesh%JMAX)
         END DO
       END IF

      IF(this%calc_boundary(EAST).EQ.1)THEN
!CDIR UNROLL=4
          DO i=1,Mesh%GNUM
             var(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX) = (i+1)*var(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) &
             - i*var(Mesh%IMAX-1,Mesh%JMIN:Mesh%JMAX)
          END DO
      END IF

      IF(this%calc_boundary(SOUTH).EQ.1)THEN
!CDIR UNROLL=4
         DO j=1,Mesh%GNUM
            var(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j) = (j+1)*var(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN) &
            - j*var(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+1)
         END DO
      END IF

      IF(this%calc_boundary(NORTH).EQ.1)THEN
!CDIR UNROLL=4
         DO j=1,Mesh%GNUM
            var(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j) = (j+1)*var(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX) &
               - j*var(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-1)
         END DO
      END IF

  END SUBROUTINE CalcBoundary

  SUBROUTINE CalcAccelGravity(this,Mesh,Physics,Fluxes,time,pvar) 
  IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: gravptr
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,time,pvar
    INTENT(INOUT)     :: Physics
    !------------------------------------------------------------------------!
    IF (Update(this,time,2)) THEN
      ! reset gterm
      this%accel(:,:,:) = 0.
      ! go through all gravity terms in the list
      gravptr => this%glist
      DO WHILE (ASSOCIATED(gravptr))
        ! call specific subroutine
!CDIR IEXPAND
        SELECT CASE(GetType(gravptr))
        CASE(POINTMASS)
          this%ptr3 => &
              GetAccelGravity_pointmass(gravptr,Mesh,Physics,Fluxes,time,pvar)
        CASE(POINTMASS_BINARY)
          this%ptr3 => &
              GetAccelGravity_binary(gravptr,Mesh,Physics,time,pvar)
        CASE(MONOPOL)
          this%ptr3 => &
              GetAccelGravity_monopol(gravptr,Mesh,Physics,time,pvar)
        CASE(MULTIGRID)
          this%ptr3 => &
              GetAccelGravity_multigrid(gravptr,Mesh,Physics,pvar)
        CASE(SPECTRAL)
          this%ptr3 => &
              GetAccelGravity_spectral(gravptr,Mesh,Physics,pvar)
        CASE(POTENTIAL)
          this%ptr3 => &
              GetAccelGravity_potential(gravptr,Mesh,Physics,time)
        CASE DEFAULT
          CALL Error(gravptr,"GravitySource", "unknown gravity term")
        END SELECT
        ! add to the sources
        this%accel(:,:,:) = this%accel(:,:,:) + this%ptr3(:,:,:)
        ! next source term
        gravptr => gravptr%next
      END DO
    END IF
  END SUBROUTINE CalcAccelGravity


  FUNCTION GetDiskinvHeight2(this,Mesh,Physics,Fluxes,time,pvar) RESULT(invheight2)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: invheight2
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,time,pvar
    INTENT(INOUT)     :: Physics
    !------------------------------------------------------------------------!
    IF (Physics%DIM .NE. 2) &
      CALL Error(this,"GetDiskinvHeight2", "DiskHeight is only supported in 2D")
    ! inverse of disk height 1/h**2
    CALL CalcDiskinvHeight2(this,Mesh,Physics,Fluxes,time,pvar)
!FIXME: Pointer?
    ! Get 1/h**2
    invheight2(:,:) = this%invheight2(:,:)
  END FUNCTION GetDiskinvHeight2

  FUNCTION GetDiskHeight(this,Mesh,Physics,Fluxes,time,pvar) RESULT(height)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: height
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,time,pvar
    INTENT(INOUT)     :: Physics
    !------------------------------------------------------------------------!
    IF (Physics%DIM .NE. 2) &
      CALL Error(this,"GetDiskHeight", "DiskHeight is only supported in 2D")
    IF (Update(this,time,3)) THEN
      ! inverse of disk height 1/h**2
      CALL CalcDiskinvHeight2(this,Mesh,Physics,Fluxes,time,pvar)
      ! disk height
      this%height(:,:) = 1./SQRT(this%invheight2(:,:))
      ! set ghost cell data
      CALL CalcBoundary(this,Mesh,this%height(:,:))
    END IF
!FIXME: Pointer?
    ! Get h
    height(:,:) = this%height(:,:)
  END FUNCTION GetDiskHeight

  FUNCTION GetAccelGravity(this,Mesh,Physics,Fluxes,time,pvar) RESULT(ac)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM) &
                      :: ac
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes,time,pvar
    INTENT(INOUT)     :: Physics
    !------------------------------------------------------------------------!
    ! update acceleration of all gravity sources
    CALL CalcAccelGravity(this,Mesh,Physics,Fluxes,time,pvar)
!FIXME: Pointer?
    ! Get acceleration
    ac(:,:,:) = this%accel(:,:,:)
  END FUNCTION GetAccelGravity

  FUNCTION GetinvDistanceCO(this,Mesh,Physics,time,pvar,component) RESULT(invdis)
  IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: invdis
    INTEGER, OPTIONAL :: component
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: gravptr
    INTEGER           :: found
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,time,component,pvar
    INTENT(INOUT)     :: Physics
    !------------------------------------------------------------------------!
    ! go through all gravity terms in the list
    gravptr => this%glist
    ! how many found?
    found = 0
    DO WHILE (ASSOCIATED(gravptr))
       ! call specific subroutine
!CDIR IEXPAND
       SELECT CASE(GetType(gravptr))
       CASE(POINTMASS)
          IF (present(component) .AND. component .NE. 1) & 
            CALL Error(this, "GetinvDistanceCO", "POINTMASS has only ONE component!")
          invdis(:,:) =&
              GetinvDistanceCO_pointmass(gravptr,Mesh)
          found = found + 1
       CASE(POINTMASS_BINARY)
          invdis(:,:) =&
              GetinvDistanceCO_binary(gravptr,Mesh,Physics,time,pvar,component)
          found = found + 1
       END SELECT
       ! next source term
       gravptr => gravptr%next
    END DO
    IF (found .NE. 1) CALL Error(this, "GetinvDistanceCO", &
                 "Huh! There isn't exactly one pointmass source!")
  END FUNCTION GetinvDistanceCO

  FUNCTION GetDistVector(this,Mesh,Physics,time,pvar,component) RESULT(r)
  IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
                      :: r
    INTEGER, OPTIONAL :: component
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: gravptr
    INTEGER           :: found
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,time,component,pvar
    INTENT(INOUT)     :: Physics
    !------------------------------------------------------------------------!
    ! go through all gravity terms in the list
    gravptr => this%glist
    ! how many found?
    found = 0
    DO WHILE (ASSOCIATED(gravptr))
       ! call specific subroutine
!CDIR IEXPAND
       SELECT CASE(GetType(gravptr))
       CASE(POINTMASS)
          IF (present(component) .AND. component .NE. 1) & 
            CALL Error(this, "GetDistVector", "POINTMASS has only ONE component!")
            r(:,:,:) =&
                GetDistVector_pointmass(gravptr,Mesh)
          found = found + 1
       CASE(POINTMASS_BINARY)
          CALL Error(this, "GetDistVector", "No distvector calculation for binary!")
!           r(:,:,:) =&
!                 GetDistVector_binary(gravptr,Mesh,Physics,time,pvar,component)
!           found = found + 1
       END SELECT
       ! next source term
       gravptr => gravptr%next
    END DO
    IF (found .NE. 1) CALL Error(this, "GetDistanceCO", &
                 "Huh! There isn't exactly one pointmass source!")
  END FUNCTION GetDistVector

  FUNCTION Update(this,time,idpos) RESULT(makeupdate)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    REAL              :: time
    INTEGER           :: idpos    ! use different ids for different routines !
    LOGICAL           :: makeupdate
    !------------------------------------------------------------------------!

    !------------------------------------------------------------------------!
    INTENT(IN)        :: time,idpos
    !------------------------------------------------------------------------!
    IF ((time.NE.this%time) .OR. (time.EQ.0.0)) THEN
      this%time = time
      ! reset id
      this%timeid = 0
      this%timeid = IBSET(this%timeid,idpos)
      makeupdate = .TRUE.
    ELSEIF (.NOT.BTEST(this%timeid,idpos)) THEN
      this%timeid = IBSET(this%timeid,idpos)
      makeupdate = .TRUE.
    ELSE
      makeupdate = .FALSE.
    END IF
  END FUNCTION

  SUBROUTINE CloseGravity(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: gravptr,ptemp
    !------------------------------------------------------------------------!
    ! release temporary/global storage
    DEALLOCATE(this%accel,this%invheight2,this%height,this%tmp,&
               this%tmp2,this%tmp4,this%calc_boundary)
    gravptr => this%glist
    ! call deallocation procedures for all source terms
    DO
       IF (.NOT.ASSOCIATED(gravptr)) EXIT
       IF (.NOT.Initialized(gravptr)) &
            CALL Error(gravptr,"CloseGravity","not initialized")
       ! call specific deconstructor
!CDIR IEXPAND
       SELECT CASE(GetType(gravptr))
       CASE(POINTMASS)
          CALL CloseGravity_pointmass(gravptr)
       CASE(POINTMASS_BINARY)
          CALL CloseGravity_binary(gravptr)
       CASE(MONOPOL)
          CALL CloseGravity_monopol(gravptr)
       CASE(MULTIGRID)
          CALL CloseGravity_multigrid(gravptr)
       CASE(SPECTRAL)
          CALL CloseGravity_spectral(gravptr)
       CASE(POTENTIAL)
         CALL CloseGravity_potential(gravptr)
       END SELECT
       ! deallocate source term structure
       ptemp=>gravptr
       gravptr=>gravptr%next
       DEALLOCATE(ptemp)
    END DO
  END SUBROUTINE CloseGravity

END MODULE gravity_generic
