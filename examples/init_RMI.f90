!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_RMI.f90                                                      #
!#                                                                           #
!# Copyright (C) 2006-2012                                                   #
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
! Program and data initialization for Richtmyer-Meshkov instability
! [1] Richtmyer, R. D.: Taylor instability in a shock acceleration of
!     compressible fluids, Comm. Pure. Appl. Math. 13 (1960), 297-319;
!     DOI: 10.1002/cpa.3160130207
!  [2] E. E. Meshkov (Евгений Евграфович Мешков)
!      "Instability of the Interface of Two Gases Accelerated by a Shock Wave"
!      Soviet Fluid Dynamics 4 ,101-104 (1969)
!----------------------------------------------------------------------------!
PROGRAM Init
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE timedisc_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: TSIM  = 200.0      ! simulation time
  REAL, PARAMETER    :: GAMMA = 1.4        ! ratio of specific heats
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CARTESIAN   ! geometry of the mesh
  INTEGER, PARAMETER :: XRES = 300         ! resolution
  INTEGER, PARAMETER :: YRES = 200
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 100         ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'RMI' 
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

CALL InitFosite(Sim)

CALL InitProgram(Sim%Mesh, Sim%Physics, Sim%Fluxes, Sim%Timedisc, &
                 Sim%Datafile, Sim%Logfile)

! set initial condition
CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)

CALL RunFosite(Sim)

CALL CloseFosite(Sim)

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
    ! mesh settings
    CALL InitMesh(Mesh,&
         meshtype = MIDPOINT, &
         geometry = CARTESIAN, &
             inum = XRES, &                 ! resolution in x and            !
             jnum = YRES, &                 !   y direction                  !             
             xmin = 0.0, &
             xmax = 60.0, &
             ymin = 0.0, &
             ymax = 40.0)

    ! physics settings
    CALL InitPhysics(Physics,Mesh, &
         problem = EULER2D, &
         gamma   = GAMMA, &                 ! ratio of specific heats        !
         dpmax   = 1.0)                     ! for advanced time step control !

    ! flux calculation and reconstruction method
    CALL InitFluxes(Fluxes,Mesh,Physics, &
         order     = LINEAR, &
         variables = CONSERVATIVE, &        ! vars. to use for reconstruction!
         limiter   = MONOCENT, &    ! one of: minmod, monocent,...   !
         theta     = 1.2)           ! optional parameter for limiter !

    ! boundary conditions
    CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = NO_GRADIENTS, &
         eastern  = NO_GRADIENTS, &
         southern = PERIODIC, &
         northern = PERIODIC)
    
    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = TSIM, &
         dtlimit  = 1.0E-4, &
         maxiter  = 1000000)

     ! initialize log input/output
!!$    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,&
!!$         fileformat = BINARY, &
!!$         filename   = TRIM(ODIR) // TRIM(OFNAME) // 'log', &
!!$         filecycles = 1)

    ! initialize data input/output
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc, &
         fileformat = VTK, &
!!$         fileformat = NETCDF, &
!!$         fileformat = GNUPLOT, filecycles = 0, &
         filename   = TRIM(ODIR) // TRIM(OFNAME), &
         count      = ONUM)

  END SUBROUTINE InitProgram


  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL              :: xlength,ylength
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!

    ! Richtmyer-Meshkov instability
    xlength = Mesh%xmax-Mesh%xmin
    ylength = Mesh%ymax-Mesh%ymin

    ! density and pressure
    WHERE ( Mesh%bcenter(:,:,1) >= 0.3*xlength + &
         (1./30.)*xlength*COS(2*PI*3./ylength*Mesh%bcenter(:,:,2)) )
       Timedisc%pvar(:,:,Physics%DENSITY)  = 0.25
       Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
    ELSEWHERE ( (Mesh%bcenter(:,:,1) <= 0.1*xlength).AND. &
         (Mesh%bcenter(:,:,1) >= (1./30.)*xlength)  )
       Timedisc%pvar(:,:,Physics%DENSITY)  = 4.22
       Timedisc%pvar(:,:,Physics%PRESSURE) = 4.9
    ELSEWHERE
       Timedisc%pvar(:,:,Physics%DENSITY)  = 1.
       Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
    END WHERE

    ! velocities vanish
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: Richtmyer-Meshkov instability")

  END SUBROUTINE InitData

END PROGRAM Init
