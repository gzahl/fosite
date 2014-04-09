!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_pointmass.f90                                             #
!#                                                                           #
!# Copyright (C) 2007-2011                                                   #
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
! source terms module for gravitational acceleration due to
! a point mass at the center of the coordinate system
!----------------------------------------------------------------------------!
MODULE sources_pointmass
  USE common_types, ONLY : Common_TYP, InitCommon
  USE sources_common, InitSources_common => InitSources
  USE boundary_common, ONLY : WEST,EAST,SOUTH,NORTH
  USE fluxes_generic, ONLY : Fluxes_TYP, GetBoundaryFlux
  USE physics_generic
  USE geometry_generic
  USE mesh_generic
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
  CHARACTER(LEN=32), PARAMETER :: source_name = "central point mass"
  CHARACTER(LEN=16), PARAMETER :: potential_name(2) = (/ &
                                  "Newton          ", &
                                  "Paczinski-Wiita " /)
  INTEGER, PARAMETER :: NEWTON = 1
  INTEGER, PARAMETER :: WIITA  = 2
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! constants
       NEWTON, WIITA, &
       ! methods
       InitSources, &
       CloseSources, &
       InitSources_pointmass, &
       InfoSources_pointmass, &
       ExternalSources_pointmass, &
       CalcTimestep_pointmass, &
       CloseSources_pointmass, &
       GetSourcesPointer, &
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

  SUBROUTINE InitSources(this,stype,sname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    INTEGER           :: stype
    CHARACTER(LEN=32) :: sname
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: newsrc, tmpsrc
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: stype,sname
    !------------------------------------------------------------------------!
    ! allocate memory for new source term
    ALLOCATE(newsrc,STAT=err)
    IF (err.NE.0) CALL Error(this,"InitSources", "Unable allocate memory!")
    
    ! basic initialization
    CALL InitSources_common(newsrc,stype,sname)

    ! add new source term to beginning of
    ! list of source terms
    IF (.NOT.ASSOCIATED(this)) THEN
       this => newsrc
       NULLIFY(this%next)
    ELSE
       tmpsrc => this
       this => newsrc
       this%next => tmpsrc
    END IF

  END SUBROUTINE InitSources


  SUBROUTINE InitSources_pointmass(this,Mesh,Physics,stype,potential,mass,cvis)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Physics_TYP) :: Physics
    INTEGER           :: stype
    INTEGER,OPTIONAL  :: potential
    REAL,OPTIONAL     :: mass,cvis
    !------------------------------------------------------------------------!
    INTEGER           :: potential_def
    INTEGER           :: err
    INTEGER           :: i,j
    REAL              :: r,a
    REAL              :: invdt_x, invdt_y
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: accel
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,stype,potential,mass,cvis
    !------------------------------------------------------------------------!
    CALL InitSources(this,stype,source_name)

    ! type of the potential
    IF (PRESENT(potential)) THEN
       potential_def = potential
    ELSE
       potential_def = NEWTON
    END IF

    ! set mass
    IF (PRESENT(mass)) THEN
       this%mass = mass
    ELSE
       this%mass = 1.0
    END IF

    ! Courant number for time step estimate
    IF (PRESENT(cvis)) THEN
       this%cvis = cvis
    ELSE
       this%cvis = 0.9
    END IF

    ! reset mass flux
    this%mdot = 0.0

    ALLOCATE(this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
         STAT = err)
    IF (err.NE.0) CALL Error(this,"InitSources_pointmass", "Unable allocate memory!")

    ! initialize potential
    SELECT CASE(potential_def)
    CASE(NEWTON) ! newtonian gravity
       CALL InitCommon(this%potential,NEWTON,potential_name(NEWTON))
       a = 0.0
    CASE(WIITA) ! pseudo-Newton Paczinski-Wiita potential
       CALL InitCommon(this%potential,WIITA,potential_name(WIITA))
       ! Schwarzschild radius
       a = 2*Physics%constants%GN * this%mass / Physics%constants%C**2
    CASE DEFAULT
       CALL Error(this,"InitSources_pointmass", "potential must be either NEWTON or WIITA")       
    END SELECT

    ! initialize gravitational acceleration
!CDIR COLLAPSE
    DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! calculate the distance to the center
          r = SQRT(Mesh%bccart(i,j,1)**2 + Mesh%bccart(i,j,2)**2)
          ! calculate cartesian components
          accel(i,j,:) = -(Physics%constants%GN * this%mass) &
               * Mesh%bccart(i,j,:) / ((r+TINY(1.0))*((r-a)**2+TINY(1.0)))
       END DO
    END DO
    accel(Mesh%IGMIN:Mesh%IMIN-1,:,:) = 0.0
    accel(Mesh%IMAX+1:Mesh%IGMAX,:,:) = 0.0
    accel(:,Mesh%JGMIN:Mesh%JMIN-1,:) = 0.0
    accel(:,Mesh%JMAX+1:Mesh%JGMAX,:) = 0.0    

    ! convert to curvilinear vector components
    CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,accel,this%accel)

    ! timestep control
    ! x-direction
    IF (Mesh%INUM.GT.1) THEN
       invdt_x = MAXVAL(ABS(this%accel(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1) &
            / Mesh%dlx(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)))
    ELSE
       ! set to zero, i.e. no limit in x-direction due to pointmass
       invdt_x = 0.0
    END IF
    ! y-direction
    IF (Mesh%JNUM.GT.1) THEN
       invdt_y = MAXVAL(ABS(this%accel(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2) &
            / Mesh%dly(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)))
    ELSE
       ! set to zero, i.e. no limit in y-direction due to pointmass
       invdt_y = 0.0
    END IF
    invdt_x = MAX(invdt_x, invdt_y)
    IF (invdt_x.LT.TINY(invdt_x)) THEN
       ! set time step to huge value, i.e. no limitation due do the
       ! pointmass module
       this%cvis = HUGE(this%cvis)
    ELSE
       ! this factor is constant throughout the simulation;
       ! devide it by SQRT(this%mass) to get the time step limit
       this%cvis = this%cvis * SQRT(this%mass/MAX(invdt_x, invdt_y))
    END IF

    SELECT CASE(GetType(Mesh%geometry))
    CASE(POLAR,LOGPOLAR,TANPOLAR,SINHPOLAR,SPHERICAL,OBLATE_SPHEROIDAL,SINHSPHERICAL)
       this%outbound = WEST
    CASE(CYLINDRICAL,TANCYLINDRICAL)
       this%outbound = SOUTH
    CASE DEFAULT
       this%outbound = 0 ! disable growth of central point mass
       CALL WARNING(this,"ExternalSources_pointmass","geometry not supported")
    END SELECT
  END SUBROUTINE InitSources_pointmass


  SUBROUTINE InfoSources_pointmass(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: mass_str
    !------------------------------------------------------------------------!
    WRITE (mass_str,'(ES8.2)') this%mass
    CALL Info(this,"            potential:         " // &
         TRIM(GetName(this%potential)) // ACHAR(10) // &
         "            mass:              " // TRIM(mass_str))
  END SUBROUTINE InfoSources_pointmass


  SUBROUTINE ExternalSources_pointmass(this,Mesh,Physics,Fluxes,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: cvar,pvar,sterm
    !------------------------------------------------------------------------!
    REAL, DIMENSION(Physics%VNUM) :: bflux
    REAL              :: oldmass
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Fluxes,cvar,pvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    IF (this%outbound.NE.0) THEN
       ! get boundary flux
#ifdef PARALLEL
       bflux(:)  = GetBoundaryFlux(Fluxes,Mesh,Physics,this%outbound,MPI_COMM_WORLD)
#else
       bflux(:)  = GetBoundaryFlux(Fluxes,Mesh,Physics,this%outbound)
#endif
       ! store old mass
       oldmass = this%mass
       ! compute new mass
       this%mass = this%mass + (bflux(Physics%DENSITY) - this%mdot)
       ! multiply gravitational acceleration with newmass/oldmass
       this%accel(:,:,:) = this%accel(:,:,:) * this%mass/oldmass
       ! store actual total mass flux
       this%mdot = bflux(Physics%DENSITY)
    END IF
    ! gravitational source terms
    CALL ExternalSources(Physics,Mesh,this%accel,pvar,cvar,sterm)
  END SUBROUTINE ExternalSources_pointmass


  PURE SUBROUTINE CalcTimestep_pointmass(this,Mesh,Physics,pvar,cvar,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar,cvar
    REAL              :: dt
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,pvar,cvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: dt
    !------------------------------------------------------------------------!
    ! largest time step due to gravitational time scale
    dt = this%cvis / SQRT(this%mass)
  END SUBROUTINE CalcTimestep_pointmass


  SUBROUTINE CloseSources_pointmass(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    IF (GetRank(this).EQ.0) THEN
       PRINT *, "-------------------------------------------------------------------"
       PRINT "(A,(ES11.3))", " central mass: ", this%mass
    END IF
    DEALLOCATE(this%accel)
    CALL CloseSources(this)
  END SUBROUTINE CloseSources_pointmass

END MODULE sources_pointmass
