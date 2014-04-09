!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_pringle2d.f90                                                #
!#                                                                           #
!# Copyright (C) 2008                                                        #
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
! Program and data initialization for pringle disk 
!----------------------------------------------------------------------------!

!**************************************!
!* IMPORTANT:                         *!
!* - compile with autodouble          *!
!* - use primitive reconstruction     *!
!**************************************!

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
#ifdef PARALLEL
  include 'mpif.h'
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
  REAL, PARAMETER :: YEAR        = 3.15576E+7    !Julian year [sec]
  REAL, PARAMETER :: PARSEC      = 3.0857E+16    !Parsec [m]
  REAL, PARAMETER :: AU          = 1.49598E+11   !AU [m] 
  REAL, PARAMETER :: CENTRALMASS = 1.98E+30      !Solarmass [kg]
  REAL, PARAMETER :: MU          = 1.0 * AU      !gauss mu
  REAL, PARAMETER :: SIGMA       = 1.0E-2 * AU   !gauss sigma
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
    INTEGER           :: onum
    REAL              :: test_stoptime, radius_s
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!
    test_stoptime = YEAR * 2000
    onum = 20

    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER2D_ISOTHERM, &
         gamma   = 1.4, &           ! ratio of specific heats        !
         cs      = 400.0, &
         dpmax   = 1.0)             ! for advanced time step control !

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)         ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         variables = PRIMITIVE, &   ! vars. to use for reconstruction!
         limiter   = MONOCENT, &    ! one of: minmod, monocent,...   !
         theta     = 1.2)           ! optional parameter for limiter !

    CALL InitMesh(Mesh,Fluxes, &
            geometry = POLAR, &
                inum = 50, &        ! resolution in x and            !
                jnum = 12, &        !   y direction                  !
                xmin = 0.1 * AU, &
                xmax = 2.6 * AU, &
                ymin = 0.0, &
                ymax = 2*PI)

    ! boundary conditions
    CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = NO_GRADIENTS, &
         eastern  = NO_GRADIENTS, &
         southern = PERIODIC, &
         northern = PERIODIC)

    ! source term due to a point mass
    CALL InitSources(Physics%sources,Mesh,Fluxes,Physics, &
         stype  = POINTMASS, &            ! grav. accel. of a point mass     !
         mass = CENTRALMASS)              ! mass of the accreting object[kg] !

    ! source term due to viscosity
    CALL InitSources(Physics%sources,Mesh,Fluxes,Physics, &
           stype    = VISCOSITY, &
           vismodel = PRINGLE, &
           cvis     = 0.5, &
           dynconst = 1.0E+10)

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = test_stoptime, &
         dtlimit  = 1E-8, &
         maxiter  = 100000000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc%pvar,Timedisc%cvar)

    ! initialize log input/output
    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc, &
         fileformat = BINARY, &
#ifdef PARALLEL
         filename   = "/tmp/pringle_disklog", &
#else
         filename   = "pringle_disklog", &
#endif
         filecycles = 1)

    ! initialize data input/output
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc, &
         fileformat = GNUPLOT, &
#ifdef PARALLEL
         filename   = "/tmp/pringle_disk", &
#else
         filename   = "pringle_disk", &
#endif
         count      = onum)
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
    INTEGER           :: ierror
    REAL              :: rho0, P0, V0_local,V0, radius_s
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(OUT)       :: pvar,cvar
    !------------------------------------------------------------------------!
    V0_local = SUM(exp(-.5 *((Mesh%bcenter(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1) &
         - MU) / SIGMA)**2))
#ifdef PARALLEL
    CALL MPI_Allreduce(V0_local,V0,1,DEFAULT_MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
#else
    V0=V0_local
#endif

    ! Initial Density
    rho0 = CENTRALMASS * 1E-1 / V0

    pvar(:,:,Physics%DENSITY) = rho0 * &
        exp(-.5 *((Mesh%bcenter(:,:,1) - MU) / SIGMA)**2) + rho0 * 1E-9
    
    ! velocities in the x-y-plane
    pvar(:,:,2) = 0.
    ! rotational velocity
    pvar(:,:,3) = SQRT(Physics%constants%GN  * CENTRALMASS / Mesh%bcenter(:,:,1))

    CALL Convert2Conservative(Physics,Mesh,pvar,cvar)

    CALL Info(Mesh, " DATA-----> initial condition: Pringle disk")
  END SUBROUTINE InitData
 END MODULE Init
