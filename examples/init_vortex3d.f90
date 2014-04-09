!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_vortex3d.f90                                                 #
!#                                                                           #
!# Copyright (C) 2006 Tobias Illenseer <tillense@ita.uni-heidelberg.de>      #
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
! Program and data initialization for 3D isentropic vortex
!----------------------------------------------------------------------------!
MODULE Init
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE output_generic
  USE logio_generic
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

  SUBROUTINE InitProgram(Mesh,Physics,Fluxes,Timedisc,Output,Logio)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Timedisc_TYP):: Timedisc
    TYPE(Output_TYP)  :: Output
    TYPE(Logio_TYP)   :: Logio
    !------------------------------------------------------------------------!
    ! Local variable declaration
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Output,Logio
    !------------------------------------------------------------------------!

    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER3D_ROTSYM, &
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

    ! mesh settings
    CALL InitMesh(Mesh,Fluxes, &
         geometry = CYLINDRICAL, &
             inum = 1, &         ! resolution in x and            !
             jnum = 100, &       !   y direction                  !             
             xmin = -0.1, &
             xmax = 0.1, &
             ymin = 0.0, &
             ymax = 5.0)

    ! boundary conditions
    CALL InitBoundary(Mesh%boundary,Mesh,Physics,PERIODIC,WEST)
    CALL InitBoundary(Mesh%boundary,Mesh,Physics,PERIODIC,EAST)
    CALL InitBoundary(Mesh%boundary,Mesh,Physics,AXIS,SOUTH)
    CALL InitBoundary(Mesh%boundary,Mesh,Physics,REFLECTING,NORTH)

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = 100.0, &
         dtlimit  = 1.0E-4, &
         maxiter  = 100000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc%pvar,Timedisc%cvar)

    ! initialize log input/output
    CALL InitLogio(Logio,Mesh,Physics,Timedisc, &
         logformat = NOLOG, &
         filename  = "vortex3d.log", &
         logdt     = 300)

    ! set parameters for data output
    CALL InitOutput(Output,Mesh,Physics,Timedisc,&
         filetype  = GNUPLOT, &
         filename  = "vortex3d.dat", &
         mode      = OVERWRITE, &
         starttime = 0.0, &
         stoptime  = Timedisc%stoptime, &
         count     = 10)

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
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2)&
                      :: cart
    INTEGER           :: i,j
    REAL              :: r,r2
    REAL              :: rho0,P0,T,T0,T0_nd,du,dv,T_nd,dT_nd,beta
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(OUT)       :: pvar,cvar
    !------------------------------------------------------------------------!

    ! ambient values of primitive variables
    rho0  = 1.
    P0    = 1.
    T0_nd = .2       ! nondimensional temperature
    T0    = Physics%mu/Physics%constants%RG * P0/rho0
    
    ! vortex strength
    beta = 5.

    CALL Convert2Cartesian(Mesh%geometry,Mesh%bcenter,cart)

    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          r  = cart(i,j,1)                      ! distance to axis
          r2 = r*r
          dv = r*beta/(2.*PI)*EXP(0.5*(1.-r2))
          dT_nd = -(Physics%gamma-1.)*beta / (8.*Physics%gamma*PI*PI) * EXP(1.-r2)
          T_nd = T0_nd + dT_nd
          T = T0/T0_nd * T_nd
          ! set primitive variables
          pvar(i,j,1) = rho0*(T/T0)**(1./(Physics%gamma-1.))
          pvar(i,j,2) = 0.
          pvar(i,j,3) = 0.
          pvar(i,j,4) = dv
          pvar(i,j,5) = Physics%constants%RG/Physics%mu * T * pvar(i,j,1)
       END DO
    END DO
    
    IF (GetType(Physics).EQ.EULER3D_ROTAMT) THEN
       pvar(:,:,4) = pvar(:,:,4)*Mesh%bhz(:,:)
    END IF
    
    CALL Convert2Conservative(Physics,Mesh,pvar,cvar)
    PRINT "(A,A)", " DATA-----> initial condition: ", &
         "3D isentropic vortex"

  END SUBROUTINE InitData

END MODULE Init
