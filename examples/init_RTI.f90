!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_RTI.f90                                                      #
!#                                                                           #
!# Copyright (C) 2008-2010                                                   # 
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
! Program and data initialization for Rayleight-Taylor instability
! References: 
! [1] Rayleigh, Lord (1883) "Investigation of the character of the equilibrium 
!     of an incompressible heavy fluid of variable density"
!     Proceedings of the London Mathematical Society 14: 170–177 
!     DOI: 10.1112/plms/s1-14.1.170
! [2] Taylor, Sir Geoffrey Ingram (1950). "The instability of liquid surfaces
!     when accelerated in a direction perpendicular to their planes"
!     Proceedings of the Royal Society of London. Series A, (1065): 192–196
!     DOI: 10.1098/rspa.1950.0052
! [3] D.H. Sharp "An overview of Rayleigh-Taylor instability" 
!     Physica D: Nonlinear Phenomena, vol. 12, Issues 1-3, July 1984, Pages 3-10
!     DOI: 10.1016/0167-2789(84)90510-4
!----------------------------------------------------------------------------!
MODULE Init
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE sources_generic
  USE timedisc_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! simulation parameters
  REAL,    PARAMETER  :: TSIM    = 10.0     ! simulation time
  REAL,    PARAMETER  :: DYNVIS  = 0.0      ! dynamic viscosity constant
  REAL,    PARAMETER  :: BULKVIS = 0.0      ! bulk viscosity constant
!!$  REAL,    PARAMETER  :: DYNVIS  = 1.0E-4   
!!$  REAL,    PARAMETER  :: BULKVIS = -6.67E-5
  ! initial condition (SI units)
  REAL,    PARAMETER  :: RHO0    = 2.0      ! density: upper region
  REAL,    PARAMETER  :: RHO1    = 1.0      ! density: lower region
  REAL,    PARAMETER  :: YACC    = 0.2      ! grav. acceleration
  REAL,    PARAMETER  :: P0      = 1.2      ! pressure at the top
  REAL,    PARAMETER  :: A0      = 0.02     ! amplitude of initial disturbance
  ! mesh settings
  INTEGER, PARAMETER  :: XRES    = 50       ! resolution in x
  INTEGER, PARAMETER  :: YRES    = 100      ! resolution in y
  REAL, PARAMETER     :: WIDTH   = 1.0      ! width of comp. domain
  REAL, PARAMETER     :: HEIGHT  = 2.0      ! height of comp. domain
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 10           ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &           ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &           ! output data file name
                     :: OFNAME = 'RTI' 
 !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! methods
       InitProgram
  !--------------------------------------------------------------------------!


CONTAINS

  SUBROUTINE InitProgram(Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Timedisc_TYP):: Timedisc
    TYPE(FILEIO_TYP)  :: Datafile
    TYPE(FILEIO_TYP)  :: Logfile
    !------------------------------------------------------------------------!
    ! Local variable declaration
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!

    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER2D, &
         gamma   = 1.4, &           ! ratio of specific heats        !
         dpmax   = 1.0)             ! for advanced time step control !

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)         ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         variables = CONSERVATIVE, &! vars. to use for reconstruction!
         limiter   = MONOCENT, &    ! one of: minmod, monocent,...   !
         theta     = 1.2)           ! optional parameter for limiter !

    ! mesh settings
    CALL InitMesh(Mesh,Fluxes, &
         geometry = CARTESIAN, &
             inum = XRES, &          ! resolution in x and            !
             jnum = YRES, &          !   y direction                  !             
             xmin = 0., &
             xmax = WIDTH, &
             ymin = 0., &
             ymax = HEIGHT)

    ! boundary conditions
    CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = REFLECTING, &
         eastern  = REFLECTING, &
         southern = REFLECTING, &
         northern = REFLECTING)

    CALL InitSources(Physics%sources,Mesh,Fluxes,Physics,Timedisc%Boundary, &
         stype  = C_ACCEL, & 
         yaccel = -YACC )            ! acceleration in y-direction

    ! viscosity source term
    IF ((DYNVIS.GT.TINY(DYNVIS)).OR.(BULKVIS.GT.TINY(BULKVIS))) &
       CALL InitSources(Physics%sources,Mesh,Fluxes,Physics,Timedisc%Boundary, &
            stype    = VISCOSITY, &
            vismodel = MOLECULAR, &
            dynconst = DYNVIS, &
           bulkconst = BULKVIS)

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = TSIM, &
         dtlimit  = 1.0E-4, &
         maxiter  = 100000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc%pvar,Timedisc%cvar)


    ! initialize log input/output
    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,&
         fileformat = BINARY, &
         filename   = TRIM(ODIR) // TRIM(OFNAME) // "log", &
         dtwall     = 1800,&
         filecycles = 1)

    ! initialize data input/output
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc, &
         fileformat = VTK, &
!!$         fileformat = GNUPLOT, &
!!$         filecycles = 0, &
         filename   = TRIM(ODIR) // TRIM(OFNAME), &
         count      = ONUM)
  END SUBROUTINE InitProgram


  SUBROUTINE InitData(Mesh,Physics,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)&
                      :: pvar,cvar
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX):: y0
    REAL, DIMENSION(2) :: accel
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh, Physics
    INTENT(OUT)       :: pvar, cvar
    !------------------------------------------------------------------------!
    ! this marks the line between the two fluids
    y0(:,:) = 0.5*Mesh%ymax + A0*COS(2*PI*Mesh%bcenter(:,:,1)/Mesh%xmax)

    ! initial hydrostatic stratification
    WHERE (Mesh%bcenter(:,:,2).GT.y0(:,:))
       ! upper fluid
       pvar(:,:,Physics%DENSITY)  = RHO0
       pvar(:,:,Physics%PRESSURE) = P0 + YACC * RHO0 * (Mesh%ymax-Mesh%bcenter(:,:,2))
    ELSEWHERE
       ! lower fluid
       pvar(:,:,Physics%DENSITY)  = RHO1
       pvar(:,:,Physics%PRESSURE) = P0 + YACC * (RHO1 * (y0(:,:)-Mesh%bcenter(:,:,2)) &
            + RHO0 * (Mesh%ymax-y0(:,:)))
    END WHERE

    ! velocity vanishes everywhere
    pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY) = 0.

    CALL Convert2Conservative(Physics,Mesh,pvar,cvar)
    CALL Info(Mesh, " DATA-----> initial condition: " // &
         "Rayleigh–Taylor instability")

  END SUBROUTINE InitData
END MODULE Init
