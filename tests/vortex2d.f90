!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: vortex2d.f90                                                      #
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
!> 2D isentropic vortex
!! \author Tobias Illenseer
!!
!! References:
!! [1] Yee, H. C. et al.: Low-dissipative high-order shock-capturing methods
!!     using characteristic-based filters, J. Comput. Phys. 150 (1999), 199-238
!!     DOI: 10.1006/jcph.1998.6177
!----------------------------------------------------------------------------!
PROGRAM vortex2d
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE timedisc_generic
  USE sources_generic
  USE sources_rotframe, ONLY : convert2RotatingFrame_rotframe
  USE common_dict
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: TSIM    = 30.0     ! simulation stop time
  REAL, PARAMETER    :: GAMMA   = 1.4      ! ratio of specific heats
  REAL, PARAMETER    :: CSISO   = 0.0      ! if .ne. 0.0 -> isothermal simulation
  ! initial condition (dimensionless units)
  REAL, PARAMETER    :: RHOINF  = 1.       ! ambient density
  REAL, PARAMETER    :: PINF    = 1.       ! ambient pressure
  REAL, PARAMETER    :: VSTR    = 5.0      ! nondimensional vortex strength
  REAL, PARAMETER    :: UINF    = 0.0      ! cartesian components of constant
  REAL, PARAMETER    :: VINF    = 0.0      !   global velocity field
  REAL, PARAMETER    :: X0      = 0.0      ! vortex position (cart. coords.)
  REAL, PARAMETER    :: Y0      = 0.0
  REAL, PARAMETER    :: R0      = 1.0      ! size of vortex
  REAL, PARAMETER    :: OMEGA   = 0.0      ! angular speed of rotational frame
                                           ! around [X0,Y0]
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CARTESIAN   ! geometry
!!$  INTEGER, PARAMETER :: MGEO = POLAR    
!!$  INTEGER, PARAMETER :: MGEO = LOGPOLAR
!!$  INTEGER, PARAMETER :: MGEO = TANPOLAR
!!$  INTEGER, PARAMETER :: MGEO = SINHPOLAR
!!$  INTEGER, PARAMETER :: MGEO = ELLIPTIC
  INTEGER, PARAMETER :: XRES = 100         ! x-resolution
  INTEGER, PARAMETER :: YRES = 100         ! y-resolution
  REAL, PARAMETER    :: RMIN = 1.0E-2      ! inner radius for polar grids
  REAL, PARAMETER    :: RMAX = 5.0         ! outer radius
  REAL, PARAMETER    :: GPAR = 1.0         ! geometry scaling parameter     !
  ! physics settings
!!$  LOGICAL, PARAMETER :: WITH_IAR = .TRUE.  ! use EULER2D_IAMROT
  LOGICAL, PARAMETER :: WITH_IAR = .FALSE.
  REAL, PARAMETER    :: PTB_AMP = 0.0E-02  ! amplitude of velocity perturbations
                                           !   set to zero to disable this
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 1           ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'vortex2d'
  !--------------------------------------------------------------------------!
  REAL               :: sigma
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  TAP_PLAN(3)

  CALL InitFosite(Sim)
  CALL MakeConfig(Sim, Sim%config)
  CALL SetupFosite(Sim)
  sigma = Run(Sim%Mesh, Sim%Physics, Sim%Timedisc)
  TAP_CHECK_SMALL(sigma,3.8E-3,"PP")

  CALL InitFosite(Sim)
  CALL MakeConfig(Sim, Sim%config)
  CALL SetAttr(Sim%config, "/fluxes/limiter", OSPRE)
  CALL SetupFosite(Sim)
  sigma = Run(Sim%Mesh, Sim%Physics, Sim%Timedisc)
  TAP_CHECK_SMALL(sigma,8.9E-3,"OSPRE")

  CALL InitFosite(Sim)
  CALL MakeConfig(Sim, Sim%config)
  CALL SetAttr(Sim%config, "/mesh/inum", 200)
  CALL SetAttr(Sim%config, "/mesh/jnum", 200)
  CALL SetupFosite(Sim)
  sigma = Run(Sim%Mesh, Sim%Physics, Sim%Timedisc)
  TAP_CHECK_SMALL(sigma,5.1E-4,"res=200")

  CALL CloseFosite(Sim)
  TAP_DONE

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    USE functions, ONLY : ASINH
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(4)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, &
                               timedisc, fluxes, sources, rotframe
    REAL              :: x1,x2,y1,y2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(CARTESIAN)
       x1 =-RMAX
       x2 = RMAX
       y1 =-RMAX 
       y2 = RMAX
    CASE(POLAR)
       x1 = RMIN
       x2 = RMAX
       y1 = 0.0 
       y2 = 2.0*PI       
    CASE(LOGPOLAR)
       x1 = LOG(RMIN/GPAR)
       x2 = LOG(RMAX/GPAR)
       y1 = 0.0 
       y2 = 2.0*PI       
    CASE(TANPOLAR)
       x1 = ATAN(RMIN/GPAR)
       x2 = ATAN(RMAX/GPAR)
       y1 = 0.0 
       y2 = 2.0*PI       
    CASE(SINHPOLAR)
       x1 = ASINH(RMIN/GPAR)
       x2 = ASINH(RMAX/GPAR)
       y1 = 0.0 
       y2 = 2.0*PI       
    CASE(ELLIPTIC)
       x1 = RMIN
       x2 = ASINH(RMAX/GPAR)
       y1 = 0.0 
       y2 = 2.0*PI       
    CASE DEFAULT
       CALL Error(Sim%Physics,"InitProgram","mesh geometry not supported for 2D isentropic vortex")
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


    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(CARTESIAN)
       bc(WEST)  = PERIODIC
       bc(EAST)  = PERIODIC
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(POLAR)
       bc(WEST)  = NOSLIP
       bc(EAST)  = REFLECTING
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(LOGPOLAR)
       bc(WEST)  = NOSLIP
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(TANPOLAR)
       bc(WEST)  = NOSLIP
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(SINHPOLAR)
       bc(WEST)  = NOSLIP
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(ELLIPTIC)
       bc(WEST)  = FOLDED
       bc(EAST)  = REFLECTING
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE DEFAULT
       CALL Error(Sim%Physics,"InitProgram","mesh geometry not supported for 2D isentropic vortex")
    END SELECT

    ! boundary conditions
    boundary => Dict("western" / bc(WEST), &
               "eastern" / bc(EAST), &
               "southern" / bc(SOUTH), &
               "northern" / bc(NORTH))

    ! physics settings
    IF (CSISO.GT.TINY(CSISO)) THEN
       physics => Dict("problem" / EULER2D_ISOTHERM, &
                 "cs"      / CSISO)                      ! isothermal sound speed  !
    ELSE
       IF (WITH_IAR) THEN
          ! REMARK: the optimal softening parameter depends on mesh geometry, limiter and
          ! possibly other settings; modify this starting with the default of 1.0, if the
          ! results show odd behaviour near the center of rotation; larger values increase
          ! softening; 0.5 give reasonable results for PP limiter on cartesian mesh
          physics => Dict("problem"   / EULER2D_IAMT, &
                     "centrot_x" / X0,"centrot_y" / Y0, & ! center of rotation      !
                     "softening" / 0.5, &                 ! softening parameter     !
                     "gamma"     / GAMMA)                 ! ratio of specific heats !
       ELSE
          physics => Dict("problem"   / EULER2D, &
                    "gamma"     / GAMMA)
       END IF
    END IF

    ! flux calculation and reconstruction method
    fluxes => Dict( &
!                "order"     / CONSTANT, &
              "order"     / LINEAR, &
!               "variables" / CONSERVATIVE, &
              "variables" / PRIMITIVE, &
              "limiter"   / PP, &                    ! PP gives better results than
!               "limiter"   / OSPRE, &                 ! any of the other limiters   
!               "limiter"   / VANLEER, &               ! on cartesian mesh, but is
!               "limiter"   / MONOCENT, &              ! rather costly
!               "output/slopes" / 1, &
              "theta"     / 1.2)                     ! optional parameter for limiter !

    ! activate inertial forces due to rotating frame if OMEGA > 0
    NULLIFY(sources)
    IF (OMEGA.GT.TINY(OMEGA)) THEN
       rotframe => Dict("stype" / ROTATING_FRAME, &
               "omega" / OMEGA, &
               "x"     / X0, &
               "y"     / Y0)
       sources => Dict("rotframe" / rotframe)
    END IF

    ! time discretization settings
    timedisc => Dict(&
         "method"   / MODIFIED_EULER, &
         "order"    / 3, &
         "cfl"      / 0.4, &
         "stoptime" / TSIM, &
         "dtlimit"  / 1.0E-4, &
         "maxiter"  / 1000000, &
!          "output/fluxes" / 1, &
         "output/iangularmomentum" / 1, &
         "output/rothalpy" / 1)

    ! initialize data input/output
     datafile => Dict(&
          "fileformat" / GNUPLOT, "filecycles" / 0, &
          "decimals" / 15, &
          "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
          "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "timedisc" / timedisc, &
             "datafile" / datafile)

    ! add sources terms
    IF (ASSOCIATED(sources)) &
        CALL SetAttr(config, "sources", sources)


  END SUBROUTINE MakeConfig


  FUNCTION Run(Mesh,Physics,Timedisc) RESULT(sigma)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j,dir,ig
    INTEGER           :: n,clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: radius
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
                      :: posvec,ephi,v0
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                      :: pvar0
    REAL              :: csinf,domega,sigma
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    IF (ABS(X0).LE.TINY(X0).AND.ABS(Y0).LE.TINY(Y0)) THEN
       ! no shift of point mass set radius and posvec to Mesh defaults
       radius(:,:) = Mesh%bradius(:,:)
       posvec(:,:,:) = Mesh%bposvec(:,:,:)
    ELSE
       ! shifted point mass position:
       ! compute curvilinear components of shift vector
       posvec(:,:,1) = X0
       posvec(:,:,2) = Y0
       CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,posvec,posvec)
       ! subtract the result from the position vector:
       ! this gives you the curvilinear components of all vectors pointing 
       ! from the point mass to the bary center of any cell on the mesh
       posvec(:,:,:) = Mesh%bposvec(:,:,:) - posvec(:,:,:)
       ! compute its absolute value
       radius(:,:) = SQRT(posvec(:,:,1)**2+posvec(:,:,2)**2)
    END IF

    ! curvilinear components of azimuthal unit vector
    ! (maybe with respect to shifted origin)
    ! from ephi = ez x er = ez x posvec/radius = ez x (rxi*exi + reta*eeta)/r
    !             = rxi/r*(ez x exi) + reta/r*(ez x eeta) = rxi/r*eeta - reta/r*exi
    ! because (ez,exi,eeta) is right handed orthonormal set of basis vectors
    ephi(:,:,1) = -posvec(:,:,2)/radius(:,:)
    ephi(:,:,2) = posvec(:,:,1)/radius(:,:)

    csinf = SQRT(GAMMA*PINF/RHOINF) ! sound speed at infinity (isentropic vortex)
    ! initial condition depends on physics;
    ! could be either isothermal or isentropic
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! local angular velocity of the vortex
          domega = 0.5*VSTR/PI*EXP(0.5*(1.-(radius(i,j)/R0)**2))
          SELECT CASE(GetType(Physics))
          CASE(EULER2D_ISOTHERM,EULER2D_ISOIAMT)
             ! density
             Timedisc%pvar(i,j,Physics%DENSITY) = RHOINF * EXP(-0.5*(R0*domega/CSISO)**2)
          CASE(EULER2D,EULER2D_IAMROT,EULER2D_IAMT)
             ! density
             ! ATTENTION: there's a factor of 1/PI missing in the density
             ! formula  eq. (3.3) in [1]
             Timedisc%pvar(i,j,Physics%DENSITY) = RHOINF * (1.0 - &
                  0.5*(GAMMA-1.0)*(R0*domega/csinf)**2  )**(1./(GAMMA-1.))
             ! pressure
             Timedisc%pvar(i,j,Physics%PRESSURE) = PINF &
                  * (Timedisc%pvar(i,j,Physics%DENSITY)/RHOINF)**GAMMA
          CASE DEFAULT
             CALL Error(Physics,"InitData","Physics must be either EULER2D or EULER2D_ISOTHERM")
          END SELECT
          Timedisc%pvar(i,j,Physics%XVELOCITY:Physics%YVELOCITY) = &
                  domega*radius(i,j)*ephi(i,j,1:2)
       END DO
    END DO

    ! compute curvilinear components of constant background velocity field
    ! and add to the vortex velocity field
    IF (ABS(UINF).GT.TINY(UINF).OR.ABS(VINF).GT.TINY(VINF)) THEN
       v0(:,:,1) = UINF
       v0(:,:,2) = VINF
       CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,v0,v0)
       Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY) = &
           Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY) + v0(:,:,1:2)
    END IF

    IF (ASSOCIATED(Physics%sources)) &
      CALL Convert2RotatingFrame_rotframe(&
          GetSourcesPointer(Physics%sources, ROTATING_FRAME),&
          Mesh,&
          Physics,&
          Timedisc%pvar)

    ! boundary conditions
    DO dir=WEST,EAST
       DO j=Mesh%JMIN,Mesh%JMAX
          DO ig=1,Mesh%GNUM
             SELECT CASE(dir)
             CASE(WEST)
                i = Mesh%IMIN-ig
             CASE(EAST)
                i = Mesh%IMAX+ig
             END SELECT
             SELECT CASE(GetType(Timedisc%Boundary(dir)))
             CASE(NOSLIP)
                Timedisc%Boundary(dir)%data(ig,j,Physics%YVELOCITY) = &
                    Timedisc%pvar(i,j,Physics%YVELOCITY)
             CASE(CUSTOM)
                Timedisc%boundary(dir)%cbtype(j,Physics%DENSITY) = CUSTOM_REFLECT
                Timedisc%boundary(dir)%cbtype(j,Physics%XVELOCITY) = CUSTOM_REFLNEG
                Timedisc%boundary(dir)%cbtype(j,Physics%YVELOCITY) = CUSTOM_FIXED
                Timedisc%Boundary(dir)%data(ig,j,Physics%YVELOCITY) = &
                    Timedisc%pvar(i,j,Physics%YVELOCITY)
             END SELECT
          END DO
       END DO
    END DO
    DO dir=SOUTH,NORTH
       DO i=Mesh%IMIN,Mesh%IMAX
          DO ig=1,Mesh%GNUM
             SELECT CASE(dir)
             CASE(SOUTH)
                j = Mesh%JMIN-ig
             CASE(NORTH)
                j = Mesh%JMAX+ig
             END SELECT
             SELECT CASE(GetType(Timedisc%Boundary(dir)))
             CASE(NOSLIP)
                Timedisc%Boundary(dir)%data(i,ig,Physics%XVELOCITY) = &
                    Timedisc%pvar(i,j,Physics%XVELOCITY)
             CASE(CUSTOM)
                Timedisc%boundary(dir)%cbtype(i,Physics%DENSITY) = CUSTOM_REFLECT
                Timedisc%boundary(dir)%cbtype(i,Physics%XVELOCITY) = CUSTOM_FIXED
                Timedisc%boundary(dir)%cbtype(i,Physics%YVELOCITY) = CUSTOM_REFLNEG
                Timedisc%Boundary(dir)%data(i,ig,Physics%XVELOCITY) = &
                    Timedisc%pvar(i,j,Physics%XVELOCITY)
             END SELECT
          END DO
       END DO
    END DO

    ! add velocity perturbations if requested
    IF (PTB_AMP.GT.TINY(PTB_AMP)) THEN
       ! Seed the random number generator with a mix from current time and mpi rank
       n = 4       ! has been choosen arbitrary
       CALL RANDOM_SEED(size=n)
       ALLOCATE(seed(n))
       CALL SYSTEM_CLOCK(COUNT=clock)
       seed = clock + (GetRank(Timedisc)+1) * (/(i-1, i=1,n)/)
       CALL RANDOM_SEED(PUT=seed)
       DEALLOCATE(seed)

       CALL RANDOM_NUMBER(v0)
       Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY) = (v0(:,:,1:2)-0.5)*PTB_AMP &
            + Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY)
    END IF

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh," DATA-----> initial condition: 2D vortex")

    pvar0 = Timedisc%pvar
    CALL RunFosite(Sim)
    sigma = SQRT(SUM((Timedisc%pvar(:,:,:)-pvar0(:,:,:))**2)/SIZE(pvar0))
  END Function Run
END PROGRAM vortex2d
