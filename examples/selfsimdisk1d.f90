!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: selfsimdisk1d.f90                                                 #
!#                                                                           #
!# Copyright (C) 2008-2012                                                   #
!# Marc Junker      <maj@astrophysik.uni-kiel.de>                            #
!# Bjoern Sperling  <sperling@astrophysik.uni-kiel.de>                       #
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
!> 1D selfsimilar beta disk solutions
!----------------------------------------------------------------------------!
PROGRAM Init
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE sources_generic
  USE timedisc_generic
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
  ! the equations are solved in non-dimensional units, setting the isothermal
  ! speed of sound to 1 and using the inverse Mach number at (t=0,R=1)
  ! (see below) as velocity scale; we use geometrical units to set the
  ! gravitational constant to unity
  ! simulation parameters
  REAL, PARAMETER :: TSIM        = 1.0E+0        ! simulation time [dyn. time]
  REAL, PARAMETER :: RE          = 1.0E+3        ! Reynolds number (at t=0,R=1)
  REAL, PARAMETER :: CHI0        = -0.5          ! power law index (-1.5<CHI0<0)
  REAL, PARAMETER :: ETA0        = 1.0E+3        ! inital mass ratio M_D / M_BH > 1
  ! mesh settings
  INTEGER, PARAMETER :: MGEO     = LOGPOLAR      ! geometry
!!$  INTEGER, PARAMETER :: MGEO     = TANPOLAR
!!$  INTEGER, PARAMETER :: MGEO     = SINHPOLAR
  INTEGER, PARAMETER :: XRES = 50          ! x-resolution
  INTEGER, PARAMETER :: YRES = 1           ! y-resolution (1D!)
  REAL, PARAMETER    :: RMIN = 1e-3        ! min radius of comp. domain
  REAL, PARAMETER    :: RMAX = 1e+1        ! max radius of comp. domain
  REAL, PARAMETER    :: GPAR = 1.0         ! geometry scaling parameter
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 100         ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'betadisk'
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim%Timedisc, Sim%Mesh, Sim%Physics, Sim%Fluxes)
  
  CALL RunFosite(Sim)

  CALL CloseFosite(Sim)

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(4),sgbc
    DOUBLE PRECISION, PARAMETER :: GN = 6.6742E-11     ! [m^3/kg/s^2]
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, sources, &
                               grav, timedisc, fluxes, vis
    REAL              :: x1,x2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(LOGPOLAR)
       x1 = LOG(RMIN/GPAR)
       x2 = LOG(RMAX/GPAR)
    CASE(TANPOLAR)
       x1 = ATAN(RMIN/GPAR)
       x2 = ATAN(RMAX/GPAR)
    CASE(SINHPOLAR)
       x1 = LOG(RMIN/GPAR+SQRT(1.0+(RMIN/GPAR)**2))  ! = ASINH(RIN))
       x2 = LOG(RMAX/GPAR+SQRT(1.0+(RMAX/GPAR)**2))
    END SELECT
    

    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / MGEO, &
           "inum"     / XRES, &
           "jnum"     / YRES, &
           "xmin"     / x1, &
           "xmax"     / x2, &
           "ymin"     / 0.0, &
           "ymax"     / (2.0*PI), &
           "gparam"   / GPAR)

    ! physics settings
    physics => Dict("problem" / EULER2D_ISOTHERM, &
              "units"   / GEOMETRICAL, &        ! Newton constant = 1 !
              "cs"      / 1.0)                 ! fixed at 1.0

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "variables" / PRIMITIVE, &   ! vars. to use for reconstruction!
             "limiter"   / MINMOD, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    boundary => Dict( &
           "western"  / ABSORBING, &
!!$           "western"  / NO_GRADIENTS, &
!!$           "western"  / EXTRAPOLATION, &
!!$           "eastern"  / REFLECTING, &
!!$           "eastern"  / ABSORBING, &
           "eastern"  / CUSTOM, &
!!$           "eastern"  / NO_GRADIENTS, &
           "southern" / PERIODIC, &
           "northern" / PERIODIC)

    ! source term due to a point mass
    grav => Dict("stype" / GRAVITY, &
                  "cvis" / 0.9, &
            "pmass/gtype"     / POINTMASS, &                ! grav. accel. of a point mass    !
            "pmass/potential" / NEWTON, &                   ! type of gravitational potential !
            "pmass/mass"      / 1., & 
!            "pmass/outbound" / 0, &      ! turn off accretion
            "selfg/gtype" / MONOPOL)  ! selfgravitational source term

    ! viscosity source term
    vis => Dict("stype"     / VISCOSITY, &
          "vismodel"  / BETA, &
          "dynconst"  / (1./RE), &
          "bulkconst" / 0.0, &
          "cvis"      / 0.5)

    sources => Dict("grav"       / grav, &
              "vis"         / vis)

    ! time discretization settings
    timedisc => Dict( &
           "method"   / MODIFIED_EULER, &
           "order"    / 3, &
           "cfl"      / 0.3, &
           "stoptime" / TSIM, &
           "dtlimit"  / (1.0E-8*TSIM), &
           "maxiter"  / 2000000000)

    ! initialize log input/output
!!$    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,, &
!!$         fileformat = BINARY,, &
!!$         filename   => Dict("betadisklog",, &
!!$         dtwall     = 3600,, &
!!$         filecycles = 1)         

    ! initialize data input/output
    datafile => Dict("fileformat" / GNUPLOT, "filecycles" / 0, &
               "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
               "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "sources"  / sources, &
             "timedisc" / timedisc, &
!             "logfile"  / logfile, &
             "datafile" / datafile)
  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Timedisc,Mesh,Physics,Fluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP):: Timedisc
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    !------------------------------------------------------------------------!
    ! Local variable declaration
    CHARACTER(LEN=32) :: info_str
    INTEGER           :: i,j,i0,j0
    REAL              :: r,Sigma0,mass
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: radius
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! compute distance to the origin
    radius(:,:) = SQRT(Mesh%bccart(:,:,1)**2+Mesh%bccart(:,:,2)**2)

    ! radial velocity vanishes
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.0

    ! set dimensionless surface density to power law
    Sigma0 = ETA0/RMAX**(2*CHI0+3.) * (CHI0+1.5)/PI
    DO i=Mesh%IGMIN,Mesh%IGMAX
       DO j=Mesh%JGMIN,Mesh%JGMAX
          Timedisc%pvar(i,j,Physics%DENSITY) = Sigma0 * radius(i,j)**(2*CHI0+1.)
       END DO
    END DO

    ! set rotational velocity to balance gravity
    Timedisc%pvar(Mesh%IGMIN:Mesh%IMIN-1,:,Physics%YVELOCITY) = 0.0
    DO j0=Mesh%JMIN,Mesh%JMAX
       DO i0=Mesh%IMIN,Mesh%IGMAX
          ! compute enclosed mass
          mass = 0.0
          DO j=Mesh%JMIN,Mesh%JMAX
             DO i=Mesh%IMIN,Mesh%IGMAX
                IF (radius(i,j).LT.radius(i0,j0)) THEN
                   mass = mass + Timedisc%pvar(i,j,Physics%DENSITY) &
                        * Mesh%volume(i,j)
                END IF
             END DO
          END DO
          Timedisc%pvar(i0,j0,Physics%YVELOCITY) = SQRT( & 
               (1.0+mass) / (radius(i0,j0)+TINY(radius(i0,j0))))
       END DO
    END DO

    ! custom boundary conditions if requested
    IF (GetType(Timedisc%boundary(EAST)).EQ.CUSTOM) THEN
       Timedisc%boundary(EAST)%cbtype(:,Physics%DENSITY) = CUSTOM_REFLECT
       Timedisc%boundary(EAST)%cbtype(:,Physics%XVELOCITY) = CUSTOM_REFLECT
       Timedisc%boundary(EAST)%cbtype(:,Physics%YVELOCITY) = CUSTOM_FIXED
       Timedisc%Boundary(EAST)%data(:,:,Physics%YVELOCITY) = &
            Timedisc%pvar(Mesh%IMAX+1:Mesh%IGMAX,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY)   
    END IF

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh," DATA-----> initial condition: " // "power law rotation curve")
  END SUBROUTINE InitData
 END PROGRAM Init
