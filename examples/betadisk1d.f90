!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: betadisk1d.f90                                                    #
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
!> 1D accretion disk
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
  USE gravity_generic
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
  ! general constants
  REAL, PARAMETER :: YEAR        = 3.15576E+7    ! Julian year [sec]
  REAL, PARAMETER :: PARSEC      = 3.0857E+16    ! parsec [m]
  REAL, PARAMETER :: AU          = 1.49598E+11   ! astronomical unit [m]
  REAL, PARAMETER :: MSUN        = 1.989E+30     ! solar mass [kg]
  REAL, PARAMETER :: RG          = 8.31          ! gas constant
  ! simulation parameters
  REAL, PARAMETER :: SIMTIME     = 1.0E7*YEAR    ! simulation time
  REAL, PARAMETER :: CENTRALMASS = 1.0E4*MSUN    ! central mass
  REAL, PARAMETER :: MDISK       = 1.0E6*MSUN    ! initial disk mass
  REAL, PARAMETER :: MU          = 0.602E-03     ! mean molecular mass
  REAL, PARAMETER :: TEMP        = 1.0E2         ! disk temperature [K]
  REAL, PARAMETER :: BETA_VIS    = 1.0E-3        ! viscosity parameter
  ! computational domain
  REAL, PARAMETER :: RMIN        = 1.0E-1*PARSEC ! inner radius [m]
  REAL, PARAMETER :: RMAX        = 1.0E+5*PARSEC ! outer radius [m]
  REAL, PARAMETER :: RSCALE      = 1.0E-0*PARSEC ! geometry scale
  ! mesh settings
!  INTEGER, PARAMETER :: MGEO     = LOGPOLAR      ! geometry
!!$  INTEGER, PARAMETER :: MGEO     = TANPOLAR
  INTEGER, PARAMETER :: MGEO     = SINHPOLAR
  ! grid resolution
  INTEGER, PARAMETER :: RES      = 252            ! radial grid points
  ! number of output time steps
  INTEGER, PARAMETER :: ONUM     = 100
  ! output file name (without extension)
  CHARACTER(LEN=*), PARAMETER  :: OFNAME = "betadisk_isotherm-t1e9-m1e3-md1e9-T1e2-beta1e-3"
  ! initial condition
  !   1 : power low disk
  !   2 : exponential disk
  !   3 : Kalnajs disk
  INTEGER, PARAMETER :: IC       = 1
  REAL, PARAMETER :: RCORE       = 1.0E-2*PARSEC ! constant density core radius
  REAL, PARAMETER :: RTRAN       = 1.0E+2*PARSEC ! transition radius
  REAL, PARAMETER :: RDISK       = 1.0E+3*PARSEC ! edge of the disk
  ! slopes for power law disk
  REAL, PARAMETER :: AA          = 1.0           ! between RCORE and RTRANS
  REAL, PARAMETER :: BB          = 6.0           ! between RTRANS and RDISK
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)          :: Sim
  !--------------------------------------------------------------------------!

CALL InitFosite(Sim)

CALL MakeConfig(Sim%config)

!CALL PrintDict(config)

CALL SetupFosite(Sim)

! set initial condition
CALL InitData(Sim%Timedisc,Sim%Mesh,Sim%Physics)

CALL RunFosite(Sim)

CALL CloseFosite(Sim)

CONTAINS
  SUBROUTINE MakeConfig(config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL              :: x1,x2
    TYPE(Dict_TYP),POINTER :: fluxes,mesh,boundary,timedisc,logfile,datafile,&
                              sources,reconstruction,grav,vis,cooling,&
                              selfgravity,physics
    !------------------------------------------------------------------------!
    ! physics settings
!     physics => Dict("problem" / EULER2D_ISOTHERM, "cs" / SQRT(RG*TEMP/MU), &
    physics => Dict("problem" / EULER2D, "gamma" / 1.4, &
              "mu" / MU)

    ! numerical scheme for flux calculation
    fluxes => Dict("order"     / LINEAR, &
                   "variables" / PRIMITIVE, & ! vars. to use for reconstruction!
                   "limiter" / MINMOD, &      ! one of: minmod, monocent,...   !
                   "theta" / 1.2)            ! optional parameter for limiter !

    ! mesh settings and boundary conditions

    SELECT CASE(MGEO)
    CASE(LOGPOLAR)
       x1 = LOG(RMIN/RSCALE)
       x2 = LOG(RMAX/RSCALE)
    CASE(TANPOLAR)
       x1 = ATAN(RMIN/RSCALE)
       x2 = ATAN(RMAX/RSCALE)
    CASE(SINHPOLAR)
       x1 = LOG(RMIN/RSCALE+SQRT(1.0+(RMIN/RSCALE)**2))  ! = ASINH(RIN))
       x2 = LOG(RMAX/RSCALE+SQRT(1.0+(RMAX/RSCALE)**2))
    END SELECT

    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / MGEO, &
           "inum" / RES, &
           "jnum" / 1, &
           "xmin" / x1, &
           "xmax" / x2, &
           "ymin" / 0.0, &
           "ymax" / (2.0*PI), &
           "gparam" / RSCALE)

    boundary => Dict(&
               "western" / ABSORBING, &
!               "western" / NO_GRADIENTS, &
               "eastern" / REFLECTING, &
!               "eastern" / NO_GRADIENTS, &
!               "eastern" / FIXED, &
!               "eastern" / FARFIELD, &
               "southern" / NO_GRADIENTS, &
               "northern" / NO_GRADIENTS)

    ! source term due to gravity
    grav => Dict("stype" / GRAVITY, &
                 "cvis" / 0.9, &
            "pmass/gtype" / POINTMASS, &     ! grav. accel. of a point mass    !
            "pmass/potential" / NEWTON, &    ! type of gravitational potential !
            "pmass/mass" / CENTRALMASS, &    ! mass [kg]
            "gravmono/gtype" / MONOPOL)

!!$    ! turn off accretion
!!$    Physics%sources%outbound = 0

    ! viscosity source term
    vis => Dict("stype" / VISCOSITY, &
          "vismodel" / BETA, &
!          "vismodel" / ALPHA, &
          "dynconst" / BETA_VIS, &
          "cvis" / 0.3)

    ! cooling source term
    cooling => Dict("stype" / DISK_COOLING, &
              "cvis" / 0.4)            ! cooling Courant number !

    sources => Dict("grav" / grav, &
!              "cooling" / cooling, &
              "viscosity" / vis)

    ! time discretization settings
    timedisc => Dict(&
               "method" / DORMAND_PRINCE, &
               "cfl" / 0.4, &
               "tol_rel" / 1.0E-3, &
!               "tol_abs" / (/ 0.0, 1.0E-3, 1.0E-3 /), &        ! EULER2D_ISOTHERM
               "tol_abs" / (/ 0.0, 1.0E-3, 1.0E-3, 1.0E-0 /), & ! EULER2D
               "stoptime" / SIMTIME, &
               "dtlimit" / 1.0E+3, &
               "maxiter" / 2000000000) 

    ! initialize data input/output
    datafile => Dict("fileformat" / GNUPLOT, &
               "filename" / OFNAME, &
               "count" / ONUM, &
               "filecycles" / 0)
    
    config => Dict("physics" / physics, &
             "fluxes" / fluxes, &
             "mesh" / mesh, &
             "boundary" / boundary, &
             "sources" / sources, &
             "timedisc" / timedisc, &
             "datafile" / datafile)

  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Timedisc,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP):: Timedisc
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Sources_TYP),POINTER :: sp
    INTEGER           :: i0,j0,i,j
    REAL              :: a0, a1
    REAL              :: sigma0,sigma_inf
    REAL              :: epsilon
    REAL              :: Rs
    REAL              :: v_r
    REAL              :: r,rho,v,omega,omega0,mass,am
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: sigma,radius
    CHARACTER(LEN=80) :: teststr
    CHARACTER(LEN=32) :: info_str   
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! Schwarzschildradius of the black hole
    Rs = 2.0 * Physics%constants%GN * CENTRALMASS / Physics%constants%C**2
    ! distance to the origin i.e. the axis
    radius(:,:) = SQRT(Mesh%bccart(:,:,1)**2 + Mesh%bccart(:,:,2)**2)

    ! initial condition for surface density
    SELECT CASE(IC)
    CASE(1)
       teststr = "power law disk"
       !                     { 1                                 for RMIN  < r <= RCORE
       ! sigma(r) = sigma0 * { (RCORE/r)**AA                     for RCORE < r <= RTRAN
       !                     { (RCORE/RTRAN)**AA * (RTRAN/r)**BB for RTRAN < r <= RDISK
       a0      = RCORE / RTRAN
       a1      = RTRAN / RDISK
       sigma0 = MDISK / (PI * RCORE**2) / ( -AA/(2.0-AA) + a0**(AA-2.0) &
            * (2.0/(2.0-AA) - 2.0/(2.0-BB) * (1.0 - a1**(BB-2.0)) ) )
       WHERE (radius(:,:).LE.RCORE)
          sigma(:,:) = sigma0
       ELSEWHERE (radius(:,:).LE.RTRAN)
          sigma(:,:) = sigma0 * (RCORE/radius(:,:))**AA
       ELSEWHERE (radius(:,:).LE.RDISK)
          sigma(:,:) = sigma0 * a0**AA * (RTRAN/radius(:,:))**BB
       ELSEWHERE
          sigma(:,:) = sigma0 * a0**AA * a1**BB
       END WHERE
    CASE(2)
       teststr = "exponential disk"
       ! sigma(r) = sigma0 * exp(-r/R0)
       sigma_inf = 1.0E-10           ! minimal sigma
       sigma0    = 1.0               ! is rescaled below
!!$       sigma0 = MDISK / (2 * PI * RCORE**2) / ( 1.0 - (1.0-RDISK/RCORE) * EXP(-RDISK/RCORE)  )
       sigma(:,:) = sigma_inf + sigma0 * radius(:,:)/RCORE
    CASE(3)
       teststr = "Kalnajs disk"
       sigma_inf = 1.0E-10           ! minimal sigma
       sigma0    = MDISK * 3.0 / (2 * PI * RDISK**2)
       WHERE (radius(:,:).LE.RDISK)
          sigma(:,:) = sigma0 * SQRT(1.0 - (radius(:,:)/RDISK)**2)
       ELSEWHERE
          sigma(:,:) = sigma_inf
       END WHERE
    END SELECT

    ! set initial data
    Timedisc%pvar(:,:,Physics%DENSITY) = sigma(:,:)
    ! vanishing radial velocity
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.0
    ! keplerian velocity (with respect to monopol approximation)
    sp => GetSourcesPointer(Physics%sources,GRAVITY)
    ! compute enclosed mass
    DO j0=Mesh%JMIN,Mesh%JMAX
       DO i0=Mesh%IMIN,Mesh%IMAX
          ! compute enclosed mass
          mass = 0.0
          IF (ASSOCIATED(GetGravityPointer(sp%glist,MONOPOL))) THEN
             DO j=Mesh%JMIN,Mesh%JMAX
                DO i=Mesh%IMIN,Mesh%IGMAX
                   IF (radius(i,j).LT.radius(i0,j0)) THEN
                      mass = mass + Timedisc%pvar(i,j,Physics%DENSITY) &
                           * Mesh%volume(i,j)
                   END IF
                END DO
             END DO
          END IF
!!$          Timedisc%pvar(i0,j0,Physics%XVELOCITY) = -1.44E+02 * SQRT(PARSEC/radius(i0,j0))
          Timedisc%pvar(i0,j0,Physics%YVELOCITY) = 0.7*SQRT(Physics%Constants%GN* &
               (mass + CENTRALMASS) / (radius(i0,j0)+TINY(radius(i0,j0))))
!!$          Timedisc%pvar(i0,j0,Physics%YVELOCITY) = SQRT(Physics%Constants%GN* &
!!$               CENTRALMASS / (radius(i0,j0)+TINY(radius(i0,j0))))
       END DO
    END DO

    ! to make sure the total mass is really the disks mass
    DO j=Mesh%JMIN,Mesh%JMAX
       DO i=Mesh%IMIN,Mesh%IMAX
          IF (radius(i,j).LT.RDISK) THEN
             mass = mass + Timedisc%pvar(i,j,Physics%DENSITY)* Mesh%volume(i,j)
          END IF
       END DO
    END DO
    ! rescale surface density
    Timedisc%pvar(:,:,Physics%DENSITY) = MDISK/mass * Timedisc%pvar(:,:,Physics%DENSITY)

    ! pressure for isothermal disk
    IF (GetType(Physics).EQ.EULER2D) &
       Timedisc%pvar(:,:,Physics%PRESSURE) = Physics%Constants%RG * TEMP / Physics%mu &
       *  Timedisc%pvar(:,:,Physics%DENSITY)
    
    ! far-field boundary conditions if requested
    IF (GetType(Timedisc%Boundary(EAST)).EQ.FIXED)  THEN
       DO j=Mesh%JMIN,Mesh%JMAX
          DO i=1,Mesh%GNUM
             Timedisc%Boundary(EAST)%data(i,j,:)   = Timedisc%pvar(Mesh%IMAX+i,j,:)
!!$             Timedisc%Boundary(EAST)%data(i,j,Physics%DENSITY)   = 1.37E+03
!!$             Timedisc%Boundary(EAST)%data(i,j,Physics%XVELOCITY) = -1.00E+02
!!$             Timedisc%Boundary(EAST)%data(i,j,Physics%YVELOCITY) = 1.1E+5
             IF (Physics%PRESSURE.NE.0) &
                  Timedisc%Boundary(EAST)%data(i,j,Physics%PRESSURE)   = 1.0E-10
          END DO
          Timedisc%Boundary(EAST)%fixed(:,:) = .TRUE.
!!$          Timedisc%Boundary(EAST)%fixed(:,Physics%DENSITY) = .FALSE.
          Timedisc%Boundary(EAST)%fixed(:,Physics%XVELOCITY) = .FALSE.
       END DO
    END IF

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh," DATA-----> initial condition: " // TRIM(teststr))
    WRITE (info_str, '(ES8.2)') MDISK/MSUN
    CALL Info(Mesh,"            disk mass:         " // TRIM(info_str) // " M_sun")
    mass = SUM(Timedisc%pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,Physics%DENSITY) &
         * Mesh%volume(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN))
    WRITE (info_str, '(ES8.2)') mass/MSUN
    CALL Info(Mesh,"            total mass:        " // TRIM(info_str) // " M_sun")
    ! compute angular momentum within the disk
    DO j=Mesh%JMIN,Mesh%JMAX
       DO i=Mesh%IMIN,Mesh%IMAX
          IF (radius(i,j).LT.RDISK) THEN
             am = am + Timedisc%pvar(i,j,Physics%DENSITY) * Mesh%bhy(i,j) &
                  * Timedisc%pvar(i,j,Physics%YVELOCITY) * Mesh%volume(i,j)
          END IF
       END DO
    END DO
    WRITE (info_str, '(ES8.2)') am
    CALL Info(Mesh,"            total angular mom. " // TRIM(info_str) // " kg*m^2/s")
  END SUBROUTINE InitData
END PROGRAM Init
