!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_bondi2d.f90                                                  #
!#                                                                           #
!# Copyright (C) 2006 - 2008                                                 #
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
! Program and data initialization for 2D Bondi accretion
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
  !--------------------------------------------------------------------------!
  PRIVATE
  ! general constants
  REAL, PARAMETER :: MSUN    = 1.989E+30   ! solar mass [kg]                 !
  ! test boundary conditions at infinity
  REAL, PARAMETER :: RHOINF  = 1.0E-20     ! density at infinity [kg/m^3]    !
  REAL, PARAMETER :: CSINF   = 1.0E+04     ! sound speed at infinity [m/s]   !
  REAL, PARAMETER :: ACCMASS = 1.0 * MSUN  ! mass of the accreting object    !
  REAL, PARAMETER :: GAMMA   = 1.4         ! ratio of specific heats         !
  ! some derives quandities
  REAL            :: RB                    ! Bondi radius                    !
  REAL            :: RIN, ROUT             ! inner & outer radius            !
  REAL            :: TS                    ! time scale                      !
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
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!
    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER2D, &
         gamma   = GAMMA, &                 ! ratio of specific heats        !
         dpmax   = 1.0)                     ! for advanced time step control !

    ! derived constants
    RB   = Physics%Constants%GN * ACCMASS / CSINF**2 ! bondi radius [m]      !
    RIN  = 1.0E-1 * RB                      ! inner & outer radius of the    !
    ROUT = 1.0E+1 * RB                      !    computational domain [m]    !
    TS   = ROUT / CSINF                     ! time scale [s]                 !

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)                 ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         variables = PRIMITIVE, &           ! vars. to use for reconstruction!
         limiter   = MONOCENT, &            ! one of: minmod, monocent,...   !
         theta     = 1.3)                   ! optional parameter for limiter !

    ! mesh settings
    CALL InitMesh(Mesh,Fluxes, &
         geometry = POLAR, &
             inum = 64, &                   ! resolution in x and            !
             jnum = 6, &                    !   y direction                  !
             xmin = RIN, &
             xmax = ROUT, &
             ymin = 0.0, &
             ymax = 2*PI)

    ! source term due to a point mass
    CALL InitSources(Physics%sources,Mesh,Fluxes,Physics, &
         stype  = POINTMASS, &            ! grav. accel. of a point mass     !
           mass = ACCMASS)                ! mass of the accreting object[kg] !

    ! boundary conditions
    CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = EXTRAPOLATION, &
         eastern  = FIXED, &
         southern = PERIODIC, &
         northern = PERIODIC)

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = 30 * TS, &
         dtlimit  = 1.0E-6 * TS, &
         maxiter  = 1000000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Fluxes,Timedisc)

  ! initialize log input/output
    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc, &
         fileformat = BINARY, &
#ifdef PARALLEL
         filename   = "/tmp/bondi2dlog", &
#else
         filename   = "bondi2d", &
#endif
         filecycles = 1)

    ! initialize data input/output
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc, &
         fileformat = GNUPLOT, &
#ifdef PARALLEL
         filename   = "/tmp/bondi2d", &
#else
         filename   = "bondi2d", &
#endif
         count      = 30)

  END SUBROUTINE InitProgram


  SUBROUTINE InitData(Mesh,Physics,Fluxes,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL              :: rho,vr,cs2
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Fluxes
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! Bondi solution at the outer boundary
    CALL bondi(Mesh%xmax/RB,GAMMA,RHOINF,CSINF,rho,vr)
    cs2 = CSINF**2 * (rho/RHOINF)**(GAMMA-1.0)
    
    ! initial condition
    Timedisc%pvar(:,:,Physics%DENSITY)   = rho
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%PRESSURE)  = rho * cs2 / GAMMA
    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)

    ! boundary condition: subsonic inflow according to Bondi's solution
    ! calculate Bondi solution for y=ymin..ymax at xmax
#ifdef PARALLEL
    IF (GetType(Timedisc%Boundary(EAST)).EQ.FIXED) THEN
#endif
       DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=1,2
             CALL bondi(Mesh%bcenter(Mesh%IMAX+i,j,1)/RB,GAMMA,RHOINF,CSINF,rho,vr)
             cs2 = CSINF**2 * (rho/RHOINF)**(GAMMA-1.0)
             ! set boundary data to either primitive or conservative values
             ! depending on the reconstruction
             IF (PrimRecon(Fluxes%reconstruction).EQV.PRIMITIVE) THEN
                Timedisc%Boundary(EAST)%data(i,j,Physics%DENSITY)   = rho
                Timedisc%Boundary(EAST)%data(i,j,Physics%XVELOCITY) = vr
                Timedisc%Boundary(EAST)%data(i,j,Physics%YVELOCITY) = 0.
                Timedisc%Boundary(EAST)%data(i,j,Physics%PRESSURE)  = rho * cs2 / GAMMA
             ELSE
                Timedisc%Boundary(EAST)%data(i,j,Physics%DENSITY)   = rho
                Timedisc%Boundary(EAST)%data(i,j,Physics%XMOMENTUM) = rho*vr
                Timedisc%Boundary(EAST)%data(i,j,Physics%YMOMENTUM) = 0.
                Timedisc%Boundary(EAST)%data(i,j,Physics%ENERGY)    = rho * &
                     (cs2 / (GAMMA*(GAMMA-1.0)) + 0.5*vr*vr)
             END IF
          END DO
       END DO
       ! this tells the boundary routine which values to fix (.TRUE.)
       ! and which to extrapolate (.FALSE.)
       Timedisc%Boundary(EAST)%fixed = (/ .TRUE., .FALSE., .TRUE., .TRUE. /)
#ifdef PARALLEL
    END IF
#endif

    CALL Info(Mesh," DATA-----> initial condition: 2D Bondi accretion")
  END SUBROUTINE InitData


  SUBROUTINE bondi(r,gamma,rhoinf,csinf,rho,vr)
    USE roots
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    ! computes the Bondi solution for rotationally symmetric accretion in    !    
    !   planar geometry                                                      !
    ! uses funcd, GetRoot                                                    !
    ! INPUT paramter:                                                        !
    !   r     : radius in units of the Bondi radius r_b = G*M/csinf**2       !
    !   gamma : ratio of specific heats (1 < gamma < 5/3)                    !
    !   rhoinf: density at infinity                                          !
    !   csinf : speed of sound at infinity                                   !
    ! OUTPUT paramter:                                                       !
    !   rho   : density @ r                                                  !
    !   vr    : radial velocity @ r                                          !
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r,gamma,rhoinf,csinf
    REAL, INTENT(OUT) :: rho,vr
    !------------------------------------------------------------------------!
    REAL, PARAMETER :: xacc = 1.0E-6     ! accuracy for root finding
    REAL :: gp1,gm1,rc,chi,lambda,psi,gr
    COMMON /funcd_parameter/ gm1, gr
    !------------------------------------------------------------------------!
    ! for convenience
    gm1 = gamma - 1.0
    gp1 = gamma + 1.0

    ! critical radius
    rc = (3.0-gamma) / 2.0

    ! critical dimensionless accretion rate
    lambda = rc**(-rc/gm1)

    ! Newton-Raphson to solve Bondis equations for psi
    chi = r / lambda
    gr  = chi**(2.*gm1/gp1) * (1./gm1 + 1./r)
    IF (r.LT.rc) THEN
       psi = GetRoot(funcd,1.0,gr,xacc)
    ELSE
       psi = GetRoot(funcd,1.0E-5,1.0,xacc)
    END IF
    
    ! return values
    rho = rhoinf * chi**(-2./gp1) / psi        ! density
    vr  = -csinf * chi**(-gm1/gp1) * psi       ! radial velocity
  END SUBROUTINE bondi


  ! find the root of fy to compute the exact Bondi solution
  SUBROUTINE funcd(y,fy,dfy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: y
    REAL, INTENT(OUT) :: fy,dfy
    !------------------------------------------------------------------------!
    REAL :: gm1,gr
    COMMON /funcd_parameter/ gm1,gr
    !------------------------------------------------------------------------!
    fy  = 0.5*y*y + y**(-gm1) / gm1 - gr
    dfy = y - y**(-gm1-1.)
  END SUBROUTINE funcd

END MODULE Init
