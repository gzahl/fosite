!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_riemann3d.f90                                                #
!#                                                                           #
!# Copyright (C) 2006-2010                                                   #
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
! Program and data initialization for 3D Riemann problem
! References:
! [1] Langseth, J. O., LeVeque, R. J.: A wave propagation method for three-
!     dimensional hyperbolic conservation laws, J. Comput. Phys., 165(1) (2000)
!     DOI: 10.1006/jcph.2000.6606
!----------------------------------------------------------------------------!
MODULE Init
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE timedisc_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! simulation parameters
  REAL, PARAMETER    :: TSIM  = 0.7        ! simulation time
  REAL, PARAMETER    :: GAMMA = 1.4        ! ratio of specific heats
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CYLINDRICAL ! geometry
  INTEGER, PARAMETER :: XRES  = 100        ! x-resolution
  INTEGER, PARAMETER :: YRES  = 150        ! y-resolution
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 7           ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'riemann3d' 
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
    TYPE(FileIO_TYP)  :: Datafile
    TYPE(FileIO_TYP)  :: Logfile
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL              :: z0,z1,r0,r1
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!
    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER3D_ROTSYM, &
         gamma   = 1.4, &           ! ratio of specific heats        !
         mu      = 0.029, &         ! mean molecular weight          !
         dpmax   = 10.0)            ! for advanced time step control !

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
         geometry = CYLINDRICAL, &
             inum = XRES, &
             jnum = YRES, &
             xmin = 0.0, &
             xmax = 1.0, &
             ymin = 0.0, &
             ymax = 1.5, &
           gparam = 1.0)

    ! boundary conditions
    CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = REFLECTING, &
         eastern  = REFLECTING, &
         southern = AXIS, &
         northern = REFLECTING)

   ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = TSIM, &
         dtlimit  = 1.0E-9, &
         maxiter  = 1000000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc)

    ! initialize log input/output
!!$    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,&
!!$         fileformat = BINARY, &
!!$         filename   = TRIM(ODIR) // TRIM(OFNAME) // 'log', &
!!$         filecycles = 1)

    ! initialize data input/output
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc, &
!!$         fileformat = VTK, &
         fileformat = GNUPLOT, filecycles = 0, &
         filename   = TRIM(ODIR) // TRIM(OFNAME), &
         count      = ONUM)
  END SUBROUTINE InitProgram


  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL              :: rmax,x0,y0
    REAL              :: P_in, P_out
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! radius and position of the bubble
    rmax = 0.2
    x0 = 0.0
    y0 = 0.4
    ! inside and ambient pressure
    P_in  = 5.0
    P_out = 1.0

    ! density
    Timedisc%pvar(:,:,Physics%DENSITY) = 1.0
    ! velocities
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%ZVELOCITY) = 0.
    ! pressure
    WHERE (((Mesh%bccart(:,:,1)-x0)**2+(Mesh%bccart(:,:,2)-y0)**2).LT.rmax**2)
       Timedisc%pvar(:,:,Physics%PRESSURE)  = P_in
    ELSEWHERE
       Timedisc%pvar(:,:,Physics%PRESSURE)  = P_out
    END WHERE

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: Spherical pressure discontinuity between walls")
  END SUBROUTINE InitData

END MODULE Init
