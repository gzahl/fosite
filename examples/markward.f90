!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: markward.f90                                                      #
!#                                                                           #
!# Copyright (C) 2012-2014                                                   #
!# Manuel Jung    <mjung@astrophysik.uni-kiel.de>                            #
!# Björn Sperling <sperling@astrophysik.uni-kiel.de>                         #
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
!> markward standard disk (self gravitating accretion disk around a SMBH)
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
  USE sources_rotframe, ONLY : Convert2RotatingFrame_rotframe
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
  REAL, PARAMETER :: GN = 6.67384E-11
  REAL, PARAMETER :: YEAR = 3.15576E+7
  REAL, PARAMETER :: PARSEC = 3.0857E+16
  REAL, PARAMETER :: AU = 1.49598E+11
  REAL, PARAMETER :: MSUN = 1.989E+30
  REAL, PARAMETER :: RG = 8.31              ! gas constant
  REAL, PARAMETER :: SIMTIME = 3.1E+2*YEAR
  REAL, PARAMETER :: MBH = 1.54E+7*MSUN
  REAL, PARAMETER :: MDISK = 6.0E+5*MSUN
  REAL, PARAMETER :: MINNERDISK = 7.9E+7*MSUN
  REAL, PARAMETER :: CENTRALMASS = 1.0*(MBH+MINNERDISK)
  REAL, PARAMETER :: MU = 6.02E-4
  REAL, PARAMETER :: TEMP = 40.0 !K
  REAL, PARAMETER :: BETA_VIS = 1.0E-5

  REAL, PARAMETER :: RMIN =  1.0E-1 * PARSEC
  REAL, PARAMETER :: RMAX = 1.E+0 * PARSEC
  INTEGER, PARAMETER :: MGEO = LOGPOLAR
  INTEGER, PARAMETER :: XRES = 128
  INTEGER, PARAMETER :: YRES = 256
  INTEGER, PARAMETER :: ONUM = 10

  REAL, PARAMETER :: noise = 0.3
  INTEGER, PARAMETER :: GREEN = 1
!  REAL, PARAMETER :: SIGMA = 3.5E-3 * PARSEC
  
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)          :: Sim
  !--------------------------------------------------------------------------!

CALL InitFosite(Sim)

CALL MakeConfig(Sim%config)

!CALL PrintDict(Sim%config)

CALL SetupFosite(Sim)

! set initial condition
CALL InitData(Sim%Timedisc,Sim%Mesh,Sim%Physics,Sim%Fluxes)

CALL RunFosite(Sim)

CALL CloseFosite(Sim)

CONTAINS
  SUBROUTINE MakeConfig(config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Dict_TYP),POINTER :: mesh,boundary,timedisc,datafile,&
                              sources,fluxes,grav,physics,rotframe,&
                              vis,psgs
    !------------------------------------------------------------------------!
    physics => Dict( &
!        "problem" / EULER2D_SGS, &
        "problem" / EULER2D_IAMT, &
        "mu" / MU, &
        "gamma" / 1.4, &
        "units" / SI)

    fluxes => Dict( &
        "order" / LINEAR, &
        "variables" / PRIMITIVE, &
        "limiter" / MINMOD, &
!        "limiter" / SWEBY, &
!        "limiter" / MONOCENT, &
!        "limiter" / OSPRE, &
        "theta" / 1.2)

    mesh => Dict( &
        "meshtype" / MIDPOINT, &
        "geometry" / MGEO, &
        "inum" / XRES, &
        "jnum" / YRES, &
        "xmin" / LOG(RMIN/PARSEC), &
        "xmax" / LOG(RMAX/PARSEC), &
        "ymin" / (-PI), &
        "ymax" / ( PI), &
        "gparam" / PARSEC, &
        "decomposition" / (/-1,1/), &
        "output/volume" / 1 )

    boundary => Dict( &
        "western" / CUSTOM, &
        "eastern" / NO_GRADIENTS, &
!        "eastern" / FARFIELD, &
        "southern" / PERIODIC, &
        "northern" / PERIODIC)

    grav => Dict("stype" / GRAVITY, &
                 "cvis" / 0.9, &
        "pmass/gtype" / POINTMASS, & 
        "pmass/potential" / NEWTON, &
        "pmass/mass" / CENTRALMASS, &
        "pmass/outbound" / 0 , &
        "self/gtype" / SPECTRAL, &
        "self/green" / GREEN)!, &
!        "self/sigma" / SIGMA)

    vis => Dict( &
        "stype" / VISCOSITY, &
        "vismodel" / BETA, &
        "dynconst" / BETA_VIS, &
!        "bulkconst" / 0.3, &
        "output/stress" / 1, &
        "output/dynvis" / 1, &
        "output/kinvis" / 1, &
        "cvis" / 0.3)

   psgs => Dict(&
         "stype" / SGS,&
         "output/dynvis" / 1, &
         "output/kinvis" / 1)

    sources => Dict( &
!        "sgs" / psgs, &
!        "viscosity" / vis, &
        "grav" / grav)



    timedisc => Dict( &
!        "method" / MODIFIED_EULER, & ! use a high order method
        "method" / DORMAND_PRINCE, &  ! <= factor 2 faster!!!)
        "fargo" / 1, &
        "tol_rel" / 1.0E-2, &
        "cfl" / 0.3, &
        "stoptime" / SIMTIME, &
        "dtlimit" / 1.0E-10, &
        "maxiter" / 2000000000)

    datafile => Dict( &
        "fileformat" / VTK, &
        "filename" / "markward", &
        "count" / ONUM)
    
    config => Dict( &
        "physics" / physics, &
        "fluxes" / fluxes, &
        "mesh" / mesh, &
        "boundary" / boundary, &
        "sources" / sources, &
        "timedisc" / timedisc, &
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
    INTEGER           :: i,j,dir,ig,ierror
    REAL              :: sigma0, mdot, sigma_inf, mass
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: rands
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    INTENT(INOUT)     :: Timedisc, Physics,Fluxes
    !------------------------------------------------------------------------!

    sigma0 = 14.1
    sigma_inf = 1.E-10

    ! Set a constant surface density with a little noise and scale it to the
    ! desired disk mass MDISK
    CALL InitRandSeed(Timedisc)
    CALL RANDOM_NUMBER(rands)
    rands = rands * noise * 2.0 + (1.0 - noise)
    Timedisc%pvar(:,:,Physics%DENSITY) = sigma_inf + sigma0 * rands*(RMIN/Mesh%bradius)

    mass = SUM(Mesh%volume * Timedisc%pvar(:,:,Physics%DENSITY))
#ifdef PARALLEL
       CALL MPI_AllReduce(MPI_IN_PLACE,mass,1,DEFAULT_MPI_REAL,MPI_SUM, &
            Mesh%comm_cart,ierror)
#endif
    Timedisc%pvar(:,:,Physics%DENSITY) = Timedisc%pvar(:,:,Physics%DENSITY) &
                                         * MDISK / mass

    ! Since the initial density is constant everywhere, there is no pressure
    ! gradient force contributing to the initial YVELOCITY
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.0
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.0

    Timedisc%pvar(:,:,Physics%PRESSURE) = TEMP*RG/MU * Timedisc%pvar(:,:,Physics%DENSITY)
    IF(GetType(Physics) .EQ. Euler2d_SGS)&
        Timedisc%pvar(:,:,Physics%SGSPRESSURE) = 1.0E-9*Timedisc%pvar(:,:,Physics%PRESSURE)

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    ! set the velocity due to the centrifugal force
    Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%XVELOCITY+Physics%DIM) = &
          GetCentrifugalVelocity(Timedisc,Mesh,Physics,Fluxes,(/0.,0.,1./))


    IF(GetType(Timedisc%Boundary(WEST)).EQ.CUSTOM) THEN
      Timedisc%boundary(WEST)%cbtype(:,Physics%DENSITY) = CUSTOM_NOGRAD
      Timedisc%boundary(WEST)%cbtype(:,Physics%XVELOCITY) = CUSTOM_OUTFLOW
      Timedisc%boundary(WEST)%cbtype(:,Physics%YVELOCITY) = CUSTOM_KEPLER
      Timedisc%boundary(WEST)%cbtype(:,Physics%Pressure) = CUSTOM_NOGRAD
      IF(GetType(Physics) .EQ. Euler2d_SGS)&
                Timedisc%boundary(WEST)%cbtype(:,Physics%SGSPRESSURE) = CUSTOM_NOGRAD
    END IF

    !v_phi at Boundary (because accel is zero at b)
    DO j=Mesh%JMIN,Mesh%JMAX
      DO i=1,Mesh%GNUM
        Timedisc%pvar(Mesh%IMAX+i,j,Physics%YVELOCITY) = &
           Timedisc%pvar(Mesh%IMAX,j,Physics%YVELOCITY) &
           * SQRT(Mesh%bradius(Mesh%IMAX,j)/Mesh%bradius(Mesh%IMAX+i,j))
      END DO
    END DO

   ! eastern
    IF (GetType(Timedisc%Boundary(EAST)).EQ.FARFIELD) THEN
       DO j=Mesh%JMIN,Mesh%JMAX
          DO i=1,Mesh%GNUM
             Timedisc%Boundary(EAST)%data(i,j,:) = Timedisc%pvar(Mesh%IMAX+i,j,:)
          END DO
       END DO
    END IF

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
  END SUBROUTINE InitData

  SUBROUTINE InitRandSeed(Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP),INTENT(IN) :: Timedisc
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    !------------------------------------------------------------------------!
    ! Initialize random number generator with a seed based on the systems time
    ! source: http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
#ifdef PARALLEL
    seed = seed + GetRank(Timedisc)
#endif
    CALL RANDOM_SEED(PUT = seed)
    DEALLOCATE(seed)
  END SUBROUTINE InitRandSeed
END PROGRAM Init
