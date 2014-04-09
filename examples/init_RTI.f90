!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_RTI.f90                                                  #
!#                                                                           #
!# Copyright (C) 2006-2008                                                   # 
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
  INTEGER :: testnum
  REAL    :: g
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

    ! 1. no viscosity old version
    ! 2. viscosity
     !air               18.27E-6
     !nitrogen          17.81E-6
     !oxygen            20.18E-6
     !carbon dioxide    14.8E-6
     !water by 10°C     1.308E−3

    g = 0.2
    testnum=2

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
         theta     = 1.3)           ! optional parameter for limiter !

    ! mesh settings
    CALL InitMesh(Mesh,Fluxes, &
         geometry = CARTESIAN, &
             inum = 50, &          ! resolution in x and            !
             jnum = 200, &          !   y direction                  !             
             xmin = 0., &
             xmax = 1.0/6.0, &
             ymin = 0., &
             ymax = 1.0)

    ! boundary conditions
    CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = REFLECTING, &
         eastern  = REFLECTING, &
         southern = REFLECTING, &
         northern = REFLECTING)

    CALL InitSources(Physics%sources,Mesh,Fluxes,Physics, &
         stype  = C_ACCEL, & 
         yaccel = -g )      ! yaccel 

    IF (testnum .eq. 2) then
    ! viscosity source term
    CALL InitSources(Physics%sources,Mesh,Fluxes,Physics, &
         stype    = VISCOSITY, &
         vismodel = MOLECULAR, &
         dynconst = 1E-3, &
         bulkconst = -6.67E-4)
    ENDIF

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = 10.0, &
         dtlimit  = 1.0E-4, &
         maxiter  = 100000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc%pvar,Timedisc%cvar)


  ! initialize log input/output
    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,&
         fileformat = BINARY, &
#ifdef PARALLEL
         filename   = "/tmp/RTIlog", &
#else
         filename   = "RTIlog", &
#endif
         stoptime   = Timedisc%stoptime, &
         dtwall     = 1800,&
         filecycles = 1)

    ! initialize data input/output
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc, &
         fileformat = GNUPLOT, &
#ifdef PARALLEL
         filename   = "/tmp/RTI", &
#else
         filename   = "RTI", &
#endif
         stoptime   = Timedisc%stoptime, &
         count      = 25, &
         filecycles = 0)
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
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX):: d
    REAL              :: rho0, rho1, p0, A
    REAL, DIMENSION(2) :: accel
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh, Physics
    INTENT(OUT)       :: pvar, cvar
    !------------------------------------------------------------------------!

    p0 = 1.2

    ! upper regions
    rho0 = 2.0 

    ! lower region
    rho1 = 1.0

    A=0.01

    do j = Mesh%JGMIN, Mesh%JGMAX, 1
       do i = Mesh%IGMIN, Mesh%IGMAX, 1
          if (Mesh%bcenter(i,j,2) .LT. Mesh%ymax/2.0) then
            pvar(i,j,4) = rho1 * g * (Mesh%ymax/2.0 - Mesh%bcenter(i,j,2))&
                          + rho0 * g * Mesh%ymax/2.0 +p0
          else
            pvar(i,j,4) = rho0 * g * (Mesh%ymax - Mesh%bcenter(i,j,2)) + p0
          end if
       end do
    end do

    do  j = Mesh%JGMIN, Mesh%JGMAX, 1
       do i = Mesh%IGMIN, Mesh%IGMAX, 1
          if (Mesh%bcenter(i,j,2) .LT. A*cos(2*3.1415*Mesh%bcenter(i,j,1)/Mesh%xmax)+Mesh%ymax/2.0) then
            pvar(i,j,Physics%DENSITY) = rho1 
          else
            pvar(i,j,Physics%DENSITY) = rho0 
          end if
       end do
    end do

    ! velocity vanishes everywhere
    pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY) = 0.

    CALL Convert2Conservative(Physics,Mesh,pvar,cvar)
    CALL Info(Mesh, " DATA-----> initial condition: " // &
         "Rayleigh–Taylor instability")

  END SUBROUTINE InitData
END MODULE Init
