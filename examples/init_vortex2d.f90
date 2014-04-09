!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_vortex2d.f90                                                 #
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
! Program and data initialization for 2D isentropic vortex
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
    INTEGER           :: geometry
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!

    ! set the geometry
    geometry = CARTESIAN
!    geometry = POLAR

    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER2D, &
         gamma   = 1.4, &           ! ratio of specific heats        !
         mu      = 0.602E-03, &     ! mean molecular weight          !
         dpmax   = 1.0)             ! for advanced time step control !

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
                xmin = -5.0, &
                xmax = 5.0, &
                ymin = -5.0, &
                ymax = 5.0)
       ! boundary conditions
       CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = PERIODIC, &
         eastern  = PERIODIC, &
         southern = PERIODIC, &
         northern = PERIODIC)
    CASE(POLAR)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes, &
            geometry = POLAR, &
                inum = 50, &        ! resolution in x and            !
                jnum = 30, &        !   y direction                  !             
                xmin = 0.0, &
                xmax = 5.0, &
                ymin = 0.0, &
                ymax = 2*PI)
       ! boundary conditions
       CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = NO_GRADIENTS, &
         eastern  = NO_GRADIENTS, &
         southern = PERIODIC, &
         northern = PERIODIC)
    CASE DEFAULT
       CALL Error(Physics,"InitProgram", &
            "geometry should be either cartesian or polar")
    END SELECT

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = 100.0, &
         dtlimit  = 1.0E-4, &
         maxiter  = 1000000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc)

    ! initialize log input/output
    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,&
         fileformat = BINARY, &
#ifdef PARALLEL
         filename   = "/tmp/vortex2dlog", &
#else
         filename   = "vortex2dlog", &
#endif
         filecycles = 1)

    ! initialize data input/output
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc, &
         fileformat = GNUPLOT, &
#ifdef PARALLEL
         filename   = "/tmp/vortex2d", &
#else
         filename   = "vortex2d", &
#endif
         count      = 10)

  END SUBROUTINE InitProgram


  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: cart
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: vxy
    REAL              :: x0,y0,x,y,r2
    REAL              :: rho0,u0,v0,P0,T,T0,T0_nd,du,dv,T_nd,dT_nd,beta
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!

    ! ambient values of primitive variables
    rho0  = 1.
    P0    = 1.
    T0_nd = .2       ! nondimensional temperature
    T0    = Physics%mu/Physics%constants%RG * P0/rho0
    
    ! vortex strength
    beta = 5.

    CALL Convert2Cartesian(Mesh%geometry,Mesh%bcenter,cart)
    
    SELECT CASE(GetType(Mesh%geometry))
    CASE(CARTESIAN)
       ! initial vortex position
       x0 = 0.5*(Mesh%xmax+Mesh%xmin)
       y0 = 0.5*(Mesh%ymax+Mesh%ymin)
       ! global velocity field
       u0 = 0.1
       v0 = 0.1
    CASE(POLAR)
       ! vortex position
       x0 = 0.0
       y0 = 0.0
       ! global velocity field should be zero in polar coords
       u0 = 0.
       v0 = 0.
    CASE DEFAULT
       CALL Error(Mesh,"InitProgram", "geometry should be either cartesian or polar")
    END SELECT

    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          x  = cart(i,j,1) - x0
          y  = cart(i,j,2) - y0
          r2 = x*x + y*y
          du = -y*beta/(2.*PI)*EXP(0.5*(1.-r2))
          dv = -du*x/y
          dT_nd = -(Physics%gamma-1.)*beta / (8.*Physics%gamma*PI*PI) * EXP(1.-r2)
          T_nd = T0_nd + dT_nd
          T = T0/T0_nd * T_nd
          ! set primitive variables
          Timedisc%pvar(i,j,Physics%DENSITY) = rho0*(T/T0)**(1./(Physics%gamma-1.))
          vxy(i,j,1) = u0 + du
          vxy(i,j,2) = v0 + dv
          Timedisc%pvar(i,j,Physics%PRESSURE) = &
               Physics%constants%RG/Physics%mu * T * Timedisc%pvar(i,j,Physics%DENSITY)
       END DO
    END DO

    ! transform velocities
    CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,vxy(:,:,:),&
         Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY))

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh," DATA-----> initial condition: 2D isentropic vortex")

  END SUBROUTINE InitData

END MODULE Init
