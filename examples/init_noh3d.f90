!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_noh3d.f90                                                    #
!#                                                                           #
!# Copyright (C) 2008 Tobias Illenseer <tillense@astrophysik.uni-kiel.de>    #
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
! Program and data initialization for 3D Noh problem
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
    CHARACTER(LEN=8)  :: geo_ext
    INTEGER           :: geometry
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!

    ! set the geometry
    geometry = CYLINDRICAL
!    geometry = SPHERICAL

    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER3D_ROTSYM, &
         gamma   = 5./3., &         ! ratio of specific heats        !
         dpmax   = 1.0E+12)         ! for advanced time step control !

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)         ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         ! IMPORTANT: always use primitive reconstruction for NOH problem
         variables = PRIMITIVE, &   ! vars. to use for reconstruction!
         limiter   = MONOCENT, &    ! one of: minmod, monocent,...   !
         theta     = 1.2)           ! optional parameter for limiter !

    SELECT CASE(geometry)
    CASE(CYLINDRICAL)
       ! mesh settings
       ! file name extension
       geo_ext = "_cyl"
       CALL InitMesh(Mesh,Fluxes, &
            geometry = CYLINDRICAL, &
                inum = 200, &       ! resolution in x and            !
                jnum = 200, &       !   y direction                  !             
                xmin = 0.0, &
                xmax = 1.0, &
                ymin = 0.0, &
                ymax = 1.0)
       ! boundary conditions
       CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = REFLECTING, &
         eastern  = NOH3D, &
         southern = REFLECTING, &
         northern = NOH3D)
    CASE(SPHERICAL)
       ! file name extension
       geo_ext = "_spher"
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes,&
            geometry = SPHERICAL, &
                inum = 200, &       ! resolution in x and            !
                jnum = 12,  &       !   y direction                  !             
                xmin = 0.0, &
                xmax = 1.0, &
                ymin = 0.0, &
                ymax = 0.5*PI)
       ! boundary conditions
       CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = REFLECTING, &
         eastern  = NOH3D, &
         southern = AXIS, &
         northern = REFLECTING)   
    CASE DEFAULT
       CALL Error(Mesh,"InitProgram", "geometry not supported")
    END SELECT

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = 1.8, &
         dtlimit  = 1.0E-6, &
         maxiter  = 100000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc)

    ! initialize log input/output
    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,&
         fileformat = BINARY, &
#ifdef PARALLEL
         filename   = "/tmp/noh3dlog" // TRIM(geo_ext), &
#else
         filename   = "noh3dlog" // TRIM(geo_ext), &
#endif
         filecycles = 1)

    ! initialize data input/output
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc,&
         fileformat = GNUPLOT, &
#ifdef PARALLEL
         filename   = "/tmp/noh3d" // TRIM(geo_ext), &
#else
         filename   = "noh3d" // TRIM(geo_ext), &
#endif
         count      = 6)

  END SUBROUTINE InitProgram


  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)  :: Physics
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Timedisc_TYP) :: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER            :: i,j
    REAL               :: rho0, vr0, P0
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: vcart
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! constant initial data
    rho0 = 1.0
    vr0  = -1.0
    P0   = 1.0E-12

    ! initial condition
    Timedisc%pvar(:,:,Physics%DENSITY)   = rho0
    Timedisc%pvar(:,:,Physics%ZVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%PRESSURE)  = P0
    ! cartesian components of the velocity
    vcart(:,:,1) = vr0 * Mesh%bccart(:,:,1) / &
         SQRT(Mesh%bccart(:,:,1)**2 + Mesh%bccart(:,:,2)**2)
    vcart(:,:,2) = vcart(:,:,1) * Mesh%bccart(:,:,2) / Mesh%bccart(:,:,1)
    CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,vcart,&
         Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY))

    ! supersonic inflow boundary conditions at the outer boundaries
    SELECT CASE(GetType(Mesh%geometry))
    CASE(CYLINDRICAL)
       ! set boundary data equal to initial values
       Timedisc%boundary(EAST)%data(:,:,:) = Timedisc%pvar(Mesh%IMAX+1:Mesh%IMAX+2,:,:)
       Timedisc%boundary(NORTH)%data(:,:,:) = Timedisc%pvar(:,Mesh%JMAX+1:Mesh%JMAX+2,:)
    CASE(SPHERICAL)
       ! set boundary data equal to initial values
       Timedisc%boundary(EAST)%data(:,:,:) = Timedisc%pvar(Mesh%IMAX+1:Mesh%IMAX+2,:,:)
    CASE DEFAULT
       CALL Error(Mesh,"InitProgram", "geometry not supported")
    END SELECT

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: 3D Noh problem")
  END SUBROUTINE InitData

END MODULE Init
