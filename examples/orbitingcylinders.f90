!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: orbitingcylinders.f90                                             #
!#                                                                           #
!# Copyright (C) 2013                                                        #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
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
!> self-gravitating orbiting cylinder test from
!! Chan, Chi-kwan; Psaltis, Dimitrios; Özel, Feryal, 2006
!! Spectral Methods for Time-dependent Studies of Accretion Flows. II. 
!! Two-dimensional Hydrodynamic Disks with Self-Gravity
!! http://adsabs.harvard.edu/abs/2006ApJ...645..506C
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
  USE functions, ONLY : Ei
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! problem setup
  REAL, PARAMETER :: sigma       = 0.1
  REAL, DIMENSION(2), PARAMETER &
!                  :: x1          = (/ COS(1.E-3), SIN(1.E-3) /)
                  :: x1          = (/ 1., 0. /)
  REAL, DIMENSION(2), PARAMETER &
                  :: x2          = -x1

  ! general constants
  REAL, PARAMETER :: GN          = 1.
  REAL, PARAMETER :: P0          = 1.
  REAL, PARAMETER :: P           = 10.
  REAL, PARAMETER :: t           = P*P0

  REAL, PARAMETER :: flaring     = 0.05
  REAL, PARAMETER :: RG          = 8.31
  REAL, PARAMETER :: MU          = 6.02E-04
  REAL, PARAMETER :: GAMMA       = 5./3.

  REAL, PARAMETER :: RMIN        = 0.2
  REAL, PARAMETER :: RMAX        = 1.8
  REAL, PARAMETER :: GPAR        = 1.0

  INTEGER, PARAMETER :: XRES     = 128
  INTEGER, PARAMETER :: YRES     = 3*XRES
  INTEGER, PARAMETER :: ONUM     = P * 1
  REAL, PARAMETER    :: omega    = 0.0

  INTEGER, PARAMETER :: GREEN    = 3
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)          :: Sim
  !--------------------------------------------------------------------------!
  
CALL InitFosite(Sim)

CALL MakeConfig(Sim%config)

!CALL PrintDict(Sim%config)

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
    TYPE(Dict_TYP),POINTER :: mesh,boundary,timedisc,datafile,&
                              fluxes,physics,grav
    REAL    :: a1,a2,a3
    !------------------------------------------------------------------------!
    physics => Dict( &
        "problem" / EULER2D_IAMT, &
        "omega" / omega, &
        "mu" / MU, &
        "gamma" / GAMMA, &
        "units" / GEOMETRICAL)

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
        "geometry" / POLAR, &
        "xmin" / RMIN, &
        "xmax" / RMAX, &
        "inum" / XRES, &
        "jnum" / YRES, &
        "ymin" / (PI*(-1.0+1.0/XRES)), &
        "ymax" / (PI*( 1.0+1.0/XRES)), &
        "decomposition" / (/ -1, 1/), &
        "output/volume" / 1 )

    boundary => Dict( &
!        "western" / REFLECTING, &
!        "eastern" / REFLECTING, &
!        "western" / NO_GRADIENTS, &
!        "eastern" / NO_GRADIENTS, &
        "western" / NOSLIP, &
        "eastern" / NOSLIP, &
        "southern" / PERIODIC, &
        "northern" / PERIODIC)

    grav => Dict( "stype" / GRAVITY, &
        "self/gtype" / SPECTRAL, &
        "self/green" / GREEN)!, &
!        "self/sigma" / SIGMA)

    timedisc => Dict( &
        "method" / MODIFIED_EULER, &
!        "method" / DUMKA, &
!        "dtmax" / 8., &
        "order" / 3, &
        "fargo" / 1, &
!        "method" / RK_FEHLBERG, &
!        "order" / 5, &
!        "tol_rel" / 1.0E-8, &
        "cfl" / 0.75, &
        "stoptime" / t, &
        "dtlimit" / 1.0E-10, &
        "maxiter" / 2000000000, &
        "output/xmomentum" / 1, &
        "output/ymomentum" / 1, &
!        "output/energy" / 1, &
        "output/geometrical_sources" / 1, &
        "output/external_sources" / 1, &
        "output/bflux" / 1)

    datafile => Dict( &
        "fileformat" / XDMF, &
        "filename" / "orbitingcylinders", &
        "count" / ONUM)
    
    config => Dict( &
        "physics" / physics, &
        "fluxes" / fluxes, &
        "mesh" / mesh, &
        "boundary" / boundary, &
        "sources/gravity" / grav, &
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
    INTEGER           :: i,j,k,dir,ig
    REAL              :: cs
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: r, r1, r2
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    INTENT(INOUT)     :: Timedisc,Physics
    !------------------------------------------------------------------------!
    r = SQRT(Mesh%bccart(:,:,1)**2+Mesh%bccart(:,:,2)**2)
    r1 = SQRT((Mesh%bccart(:,:,1)-x1(1))**2+(Mesh%bccart(:,:,2)-x1(2))**2)
    r2 = SQRT((Mesh%bccart(:,:,1)-x2(1))**2+(Mesh%bccart(:,:,2)-x2(2))**2)
    
    Timedisc%pvar(:,:,Physics%DENSITY)  = 0.02/(3.2*PI) &
      + 0.99 * (  EXP(-0.5*(r1/sigma)**2)/(2.*PI*sigma**2) &
                + EXP(-0.5*(r2/sigma)**2)/(2.*PI*sigma**2))

    Timedisc%pvar(:,:,Physics%PRESSURE) = 0.02/(3.2*PI) &
      + GN / (2.*PI*sigma**2) &
         * (  Ei(-(r1/sigma)**2) - Ei(-0.5*(r1/sigma)**2) &
            + Ei(-(r2/sigma)**2) - Ei(-0.5*(r2/sigma)**2))
         
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.0
    Timedisc%pvar(:,:,Physics%YVELOCITY) = r - r*omega

    DO dir=WEST,EAST
        IF(GetType(Timedisc%Boundary(dir)).EQ.NOSLIP) THEN
            DO j=Mesh%JMIN,Mesh%JMAX
                DO ig=1,Mesh%GNUM
                    SELECT CASE (dir)
                    CASE(WEST)
                        i = Mesh%IMIN-ig
                    CASE(EAST)
                        i = Mesh%IMAX+ig
                    END SELECT
                    Timedisc%Boundary(dir)%data(ig,j,Physics%YVELOCITY) &
                      = Timedisc%pvar(i,j,Physics%YVELOCITY)
                END DO
            END DO
        END IF
    END DO

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
  END SUBROUTINE InitData
END PROGRAM Init
