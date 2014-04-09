!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: gauss3d.f90                                                       #
!#                                                                           #
!# Copyright (C) 2006-2014                                                   #
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
!> 3D Gaussian pressure or density pulse with and without rotation
!! \author Tobias Illenseer
!!
!! [1] Illenseer, T. F., Duschl, W. J.: Two-dimensional central-upwind schemes
!!     for curvilinear grids and application to gas dynamics with angular momentum,
!!     Comput. Phys. Comm. 180 (2009), 2283-2302
!!     DOI: 10.1016/j.cpc.2009.07.016
!----------------------------------------------------------------------------!
PROGRAM gauss3d
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE timedisc_generic
  USE sources_generic
  USE common_dict
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameter
  REAL, PARAMETER :: TSIM   = 0.6         ! simulation time
  REAL, PARAMETER :: GAMMA  = 1.4         ! ratio of specific heats
  ! initial condition (dimensionless units)
  REAL, PARAMETER :: RHO0   = 1.0         ! ambient density
  REAL, PARAMETER :: RHO1   = 0.0         ! peak density above RHO0
  REAL, PARAMETER :: RWIDTH = 0.06        ! half width of the Gaussian
  REAL, PARAMETER :: P0     = 1.0         ! ambient pressure
  REAL, PARAMETER :: P1     = 1.0         ! peak pressure above P0
  REAL, PARAMETER :: PWIDTH = 0.06        ! half width of the Gaussian
  REAL, PARAMETER :: OMEGA0 = 0.0         ! angular velocity 
  REAL, PARAMETER :: ETA    = 0.0         ! dynamic viscosity (0.0 disables)
  ! location of the pulse in cylindrical coordinates
  REAL, PARAMETER :: R0     = 0.0         ! radial position 
  REAL, PARAMETER :: Z0     = 0.0         ! vertical position
  ! mesh settings
!!$  INTEGER, PARAMETER :: MGEO = SPHERICAL   ! geometry
  INTEGER, PARAMETER :: MGEO = CYLINDRICAL
!!$  INTEGER, PARAMETER :: MGEO = OBLATE_SPHEROIDAL
!!$  INTEGER, PARAMETER :: MGEO = SINHSPHERICAL
  INTEGER, PARAMETER :: XRES  = 100       ! x-resolution
  INTEGER, PARAMETER :: YRES  = 100       ! y-resolution
  REAL, PARAMETER    :: RMAX  = 1.0       ! width of square that fits into
                                          !   computational domain
  REAL, PARAMETER    :: GPAR  = 0.8       ! geometry scaling parameter
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 10         ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &         ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &         ! output data file name
                     :: OFNAME = 'gauss3d' 
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  TAP_PLAN(1)

  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)
  
  CALL RunFosite(Sim)

  CALL CloseFosite(Sim)

  TAP_CHECK(.TRUE.,"Finished simulation")
  TAP_DONE

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    USE functions, ONLY : Asinh, Acosh
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(4)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, sources, &
                               timedisc, fluxes, vis
    REAL              :: x1,x2,y1,y2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(SPHERICAL)
       x1 = 0.0
       x2 = SQRT(2.0)*RMAX
       y1 = 0.0
       y2 = 0.5*PI
       bc(WEST)  = REFLECTING
       bc(EAST)  = ABSORBING
       bc(SOUTH) = AXIS
       bc(NORTH) = REFLECTING
    CASE(CYLINDRICAL)
       x1 = 0.0
       x2 = RMAX
       y1 = 0.0
       y2 = RMAX
       bc(WEST)  = REFLECTING
       bc(EAST)  = ABSORBING
       bc(SOUTH) = AXIS
       bc(NORTH) = ABSORBING
    CASE(OBLATE_SPHEROIDAL)
       x1 = 0.0
       x2 = 0.5*Acosh(2./GPAR**2 * (1.0 + SQRT(1.0+0.25*GPAR**4)))
       y1 = 0.0
       y2 = 0.5*PI
       bc(WEST)  = REFLECTING
       bc(EAST)  = ABSORBING
       bc(SOUTH) = REFLECTING
       bc(NORTH) = AXIS
    CASE(SINHSPHERICAL)
       x1 = 0.0
       x2 = Asinh(SQRT(2.0)*RMAX/GPAR)
       y1 = 0.0
       y2 = 0.5*PI
       bc(WEST)  = REFLECTING
       bc(EAST)  = ABSORBING
       bc(SOUTH) = AXIS
       bc(NORTH) = REFLECTING
    CASE DEFAULT
       CALL Error(Sim%Physics,"InitProgram","geometry not supported for this test")
    END SELECT

    !mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / MGEO, &
           "inum"     / XRES, &
           "jnum"     / YRES, &
           "xmin"     / x1, &
           "xmax"     / x2, &
           "ymin"     / y1, &
           "ymax"     / y2, &
           "gparam"   / GPAR)

    ! boundary conditions
    boundary => Dict("western" / bc(WEST), &
               "eastern" / bc(EAST), &
               "southern" / bc(SOUTH), &
               "northern" / bc(NORTH))

    ! physics settings
    physics => Dict("problem" / EULER3D_ROTSYM, &
              "gamma"   / GAMMA)                ! ratio of specific heats        !

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "variables" / CONSERVATIVE, &        ! vars. to use for reconstruction!
             "limiter"   / MONOCENT, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    NULLIFY(sources)
    ! viscosity source term
    IF (ETA.GT.TINY(ETA)) THEN
       vis => Dict("stype"    / VISCOSITY, &
          "vismodel" / MOLECULAR, &
          "dynconst" / ETA)
       CALL SetAttr(sources, "vis", vis)
    END IF

    ! time discretization settings
    timedisc => Dict("method"    / MODIFIED_EULER, &
               "order"     / 3, &
               "cfl"       / 0.4, &
               "stoptime"  / TSIM, &
               "dtlimit"   / 1.0E-8, &
               "maxiter"   / 10000000)

    ! initialize data input/output
!    datafile => Dict("fileformat" / VTK, &
    datafile => Dict("fileformat" / GNUPLOT, "filecycles" / 0, &
               "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
               "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "timedisc" / timedisc, &
!             "logfile"  / logfile, &
             "datafile" / datafile)

    IF (ASSOCIATED(sources)) &
        CALL SetAttr(config, "sources", sources)

  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Timedisc_TYP):: Timedisc
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: radius
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: posvec
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    IF (ABS(R0).LE.TINY(R0).AND.ABS(Z0).LE.TINY(Z0)) THEN
       ! no shift of point mass set radius and posvec to Mesh defaults
       radius(:,:) = Mesh%bradius(:,:)
       posvec(:,:,:) = Mesh%bposvec(:,:,:)
    ELSE
       ! shifted point mass position:
       ! compute curvilinear components of shift vector
       posvec(:,:,1) = R0
       posvec(:,:,2) = Z0
       CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,posvec,posvec)
       ! subtract the result from the position vector:
       ! this gives you the curvilinear components of all vectors pointing 
       ! from the point mass to the bary center of any cell on the mesh
       posvec(:,:,:) = Mesh%bposvec(:,:,:) - posvec(:,:,:)
       ! compute its absolute value
       radius(:,:) = SQRT(posvec(:,:,1)**2+posvec(:,:,2)**2)
    END IF

    ! initial density and pressure
    FORALL (i=Mesh%IGMIN:Mesh%IGMAX,j=Mesh%JGMIN:Mesh%JGMAX)
       Timedisc%pvar(i,j,Physics%DENSITY) = RHO0 + RHO1*EXP(-LOG(2.0) &
            * (radius(i,j)/RWIDTH)**2)
       Timedisc%pvar(i,j,Physics%PRESSURE) = P0 + P1*EXP(-LOG(2.0) &
            * (radius(i,j)/PWIDTH)**2)
    END FORALL

    ! velocities in the x-y-plane
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.

    ! rotational velocity
    Timedisc%pvar(:,:,Physics%ZVELOCITY) = OMEGA0 * Mesh%bhz(:,:)

    ! for specific angular momentum transport
    IF (GetType(Physics).EQ.EULER3D_ROTAMT) THEN
       ! specific angular momentum
       Timedisc%pvar(:,:,Physics%ZVELOCITY) = Timedisc%pvar(:,:,Physics%ZVELOCITY)*Mesh%bhz(:,:)
    END IF
    
    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)

    CALL Info(Mesh, " DATA-----> initial condition: 3D Gaussian pulse")

  END SUBROUTINE InitData

END PROGRAM gauss3d
