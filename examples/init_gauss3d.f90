!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_gauss3d.f90                                                  #
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
! Program and data initialization for 3D Gaussian pressure or density pulse
!----------------------------------------------------------------------------!
MODULE Init
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE timedisc_generic
  USE sources_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER :: testnum
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
    INTEGER           :: onum
    REAL              :: test_stoptime
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!

    ! 3D tests with continous initial condition
    ! 1. Gaussian pressure pulse
    ! 2. Rotating Gaussian density pulse
    testnum=1

    ! set the geometry
    geometry = CYLINDRICAL
!    geometry = SPHERICAL
!    geometry = OBLATE_SPHEROIDAL

    SELECT CASE(testnum)
    CASE(1) ! Gaussian pressure pulse
       test_stoptime = 0.6
       onum = 6
    CASE(2) ! Rotating Gaussian density pulse
       test_stoptime = 1.0
       onum = 10
     CASE DEFAULT
       PRINT *, "ERROR in InitProgram: testnum should be one of 1, 2"
       STOP
    END SELECT

    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER3D_ROTSYM, &
         gamma   = 1.4, &           ! ratio of specific heats        !
         dpmax   = 1.0)             ! for advanced time step control !

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)         ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         variables = PRIMITIVE, &   ! vars. to use for reconstruction!
         limiter   = MONOCENT, &    ! one of: minmod, monocent,...   !
         theta     = 1.3)           ! optional parameter for limiter !

    SELECT CASE(geometry)
    CASE(CYLINDRICAL)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes, &
            geometry = CYLINDRICAL, &
                inum = 200, &       ! resolution in x and            !
                jnum = 200, &       !   y direction                  !             
                xmin = 0.0, &
                xmax = 1.0, &
                ymin = 0.0, &
                ymax = 1.0)
       ! boundary conditions
       IF (testnum.EQ.1) THEN
          CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
               western  = REFLECTING, &
               eastern  = NO_GRADIENTS, &
               southern = AXIS, &
               northern = NO_GRADIENTS)
       ELSE
          CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
               western  = REFLECTING, &
               eastern  = REFLECTING, &
               southern = AXIS, &
               northern = REFLECTING)
       END IF
    CASE(SPHERICAL)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes,SPHERICAL,60,30,0.0,1.5,0.0,0.5*PI)
       ! boundary conditions
       IF (testnum.EQ.1) THEN
          CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
               western  = REFLECTING, &
               eastern  = REFLECTING, &
               southern = AXIS, &
               northern = REFLECTING)
       ELSE
          CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
               western  = REFLECTING, &
               eastern  = NO_GRADIENTS, &
               southern = AXIS, &
               northern = REFLECTING)
       END IF
    CASE(OBLATE_SPHEROIDAL)
       ! mesh settings
       CALL InitMesh(Mesh,Fluxes,OBLATE_SPHEROIDAL,50,30, &
            0.0,1.4,0.0,0.5*PI, &
            0.75)                    ! optional geometry parameter    !
       ! boundary conditions
       CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = REFLECTING, &
         eastern  = REFLECTING, &
         southern = AXIS, &
         northern = AXIS)
    CASE DEFAULT
       PRINT *, "ERROR in InitProgram: geometry should be either ", ACHAR(13), &
            "cylindrical, spherical or oblate spheroidal"
       STOP
    END SELECT

    ! viscosity source term
!!$    CALL InitSources(Physics%sources,Mesh,Fluxes,Physics, &
!!$         stype    = VISCOSITY, &
!!$         vismodel = MOLECULAR, &
!!$         cvis     = 0.5, &
!!$         dynconst = 1.0E-2)

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = test_stoptime, &
         dtlimit  = 1.0E-8, &
         maxiter  = 10000000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc%pvar,Timedisc%cvar)

    ! initialize log input/output
    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc, &
         fileformat = BINARY, &
#ifdef PARALLEL
         filename   = "/tmp/gauss3dlog", &
#else
         filename   = "gauss3dlog", &
#endif
         dtwall     = 1740, &
         filecycles = 1)

    ! initialize data input/output
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc, &
         fileformat = GNUPLOT, &
#ifdef PARALLEL
         filename   = "/tmp/gauss3d", &
#else
         filename   = "gauss3d", &
#endif
         stoptime   = Timedisc%stoptime, &
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
    REAL              :: x0, y0
    REAL              :: rho0, rho1, P0, P1, hwidth, omega
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: cart
    CHARACTER(LEN=80) :: teststr
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(OUT)       :: pvar,cvar
    !------------------------------------------------------------------------!

    CALL Convert2Cartesian(Mesh%geometry,Mesh%bcenter,cart)
    
    SELECT CASE(testnum)
    CASE(1) ! 3D pressure pulse
       teststr = "3D Gaussian pressure pulse"
       ! position
       x0        = 0.0
       y0        = 0.0
       ! density
       rho0      = 1.0
       rho1      = 0.0
       ! pressure
       P0        = 1.0
       P1        = 1.0
       hwidth    = 0.06*MAXVAL(cart(:,:,:))
       ! angular velocity
       omega     = 0.0
    CASE(2) ! 3D density pulse with angular motion
       teststr = "3D Gaussian density pulse w/ rotation"
       ! position
       x0        = 0.0
       y0        = 0.4
       ! density
       rho0      = 1.0e-2
       rho1      = 10.0
       ! pressure
       P0        = 1.0
       P1        = 0.0
       hwidth    = 0.1
       ! angular velocity
       omega     = 10.0
    END SELECT

    FORALL (i=Mesh%IGMIN:Mesh%IGMAX,j=Mesh%JGMIN:Mesh%JGMAX)
       pvar(i,j,Physics%DENSITY) = rho0 + rho1*EXP(-LOG(2.0) &
            * ((cart(i,j,1)-x0)**2 + (cart(i,j,2)-y0)**2)/hwidth**2)
       pvar(i,j,Physics%PRESSURE) = P0 + P1*EXP(-LOG(2.0) &
            * ((cart(i,j,1)-x0)**2 + (cart(i,j,2)-y0)**2)/hwidth**2)
    END FORALL

    ! velocities in the x-y-plane
    pvar(:,:,Physics%XVELOCITY) = 0.
    pvar(:,:,Physics%YVELOCITY) = 0.

    ! rotational velocity
    pvar(:,:,Physics%ZVELOCITY) = omega * Mesh%bhz(:,:)

    ! for specific angular momentum transport
    IF (GetType(Physics).EQ.EULER3D_ROTAMT) THEN
       ! specific angular momentum
       pvar(:,:,Physics%ZVELOCITY) = pvar(:,:,Physics%ZVELOCITY)*Mesh%bhz(:,:)
    END IF
    
    CALL Convert2Conservative(Physics,Mesh,pvar,cvar)

    CALL Info(Mesh, " DATA-----> initial condition: " // TRIM(teststr))

  END SUBROUTINE InitData

END MODULE Init
