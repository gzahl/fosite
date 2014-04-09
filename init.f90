!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_gauss2d.f90                                                  #
!#                                                                           #
!# Copyright (C) 2006-2008                                                   #
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
! Program and data initialization for 2D Gaussian pressure pulse
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
#ifdef PARALLEL
  include 'mpif.h'
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
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
    INTEGER           :: geometry
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!

    ! set the geometry
    geometry = CARTESIAN
!    geometry = POLAR
!    geometry = LOGPOLAR

    ! physics settings
    CALL InitPhysics(Physics, &
         problem   = EULER2D, &
         gamma     = 1.4, &         ! ratio of specific heats        !
         dpmax     = 1.0)           ! for advanced time step control !

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)         ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         variables = CONSERVATIVE, &! vars. to use for reconstruction!
         limiter   = MONOCENT, &    ! one of: minmod, monocent,...   !
         theta     = 1.3)           ! optional parameter for limiter !

    SELECT CASE(geometry)
    CASE(CARTESIAN)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes, &
            geometry = CARTESIAN, &
                inum = 100, &       ! resolution in x and            !
                jnum = 100, &       !   y direction                  !             
                xmin = -0.5, &
                xmax = 0.5, &
                ymin = -0.5, &
                ymax = 0.5)
       ! boundary conditions
    CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = NO_GRADIENTS, &
         eastern  = NO_GRADIENTS, &
         southern = NO_GRADIENTS, &
         northern = NO_GRADIENTS)
    CASE(POLAR)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes, &
            geometry = POLAR, &
                inum = 50, &        ! resolution in x and            !
                jnum = 30, &        !   y direction                  !             
                xmin = 0.001, &
                xmax = 0.5, &
                ymin = 0.0, &
                ymax = 2*PI)
       ! boundary conditions
       CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = AXIS, &
         eastern  = NO_GRADIENTS, &
         southern = PERIODIC, &
         northern = PERIODIC)
    CASE(LOGPOLAR)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes, &
            geometry = LOGPOLAR, &
                inum = 50, &        ! resolution in x and            !
                jnum = 120, &       !   y direction                  !             
                xmin = LOG(0.1), &
                xmax = LOG(5.0), &
                ymin = 0.0, &
                ymax = 2*PI, &
              gparam = 0.1)
       ! boundary conditions
       CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = AXIS, &
         eastern  = NO_GRADIENTS, &
         southern = PERIODIC, &
         northern = PERIODIC)
    CASE DEFAULT
       CALL Error(Mesh,"InitProgram","geometry should be either cartesian or (log)polar")
    END SELECT

    ! viscosity source term
!!$    CALL InitSources(Physics%sources,Mesh,Fluxes,Physics, &
!!$         stype    = VISCOSITY, &
!!$         vismodel = MOLECULAR, &
!!$         dynconst = 1.0E-1)

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = 0.3, &
         dtlimit  = 1.0E-4, &
         maxiter  = 1000000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc%pvar,Timedisc%cvar)

    ! initialize log input/output
    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,&
         fileformat = BINARY, &
#ifdef PARALLEL
         filename   = "/tmp/gauss2dlog", &
#else
         filename   = "gauss2dlog", &
#endif
         stoptime   = Timedisc%stoptime, &
         dtwall     = 1800,&
         filecycles = 1)

    ! initialize data input/output
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc, &
         fileformat = BINARY, &
#ifdef PARALLEL
         filename   = "/tmp/gauss2d", &
#else
         filename   = "gauss2d", &
#endif
         stoptime   = Timedisc%stoptime, &
         count      = 10, &
         filecycles = 0)
  END SUBROUTINE InitProgram


  SUBROUTINE InitData(Mesh,Physics,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM)&
                      :: pvar,cvar
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j
    REAL              :: amplitude, hwidth, rho0, P0
#ifdef PARALLEL
    REAL              :: hwidth_all
    INTEGER           :: ierror
#endif
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: cart
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(OUT)       :: pvar,cvar
    !------------------------------------------------------------------------!

    CALL Convert2Cartesian(Mesh%geometry,Mesh%bcenter,cart)
    
    ! 2D gaussian pressure pulse
    rho0      = 1.0
    P0        = 1.0
    amplitude = 1.0
    hwidth    = 0.06*MAXVAL(cart(:,:,:))

#ifdef PARALLEL
    CALL MPI_Allreduce(hwidth,hwidth_all,1,DEFAULT_MPI_REAL,MPI_MAX, &
         Mesh%comm_cart,ierror)
    hwidth = hwidth_all
#endif

    pvar(:,:,Physics%DENSITY)   = rho0
    pvar(:,:,Physics%XVELOCITY) = 0.
    pvar(:,:,Physics%YVELOCITY) = 0.

    FORALL (i=Mesh%IGMIN:Mesh%IGMAX,j=Mesh%JGMIN:Mesh%JGMAX)
       pvar(i,j,Physics%PRESSURE) = P0 + amplitude*EXP(-LOG(2.0) * &
            (cart(i,j,1)**2+cart(i,j,2)**2)/hwidth**2)
    END FORALL

    CALL Convert2Conservative(Physics,Mesh,pvar,cvar)
    CALL Info(Mesh, " DATA-----> initial condition: " // &
         "2D gaussian pressure pulse")

  END SUBROUTINE InitData

END MODULE Init
