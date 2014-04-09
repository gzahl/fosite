!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_bondi2d.f90                                                  #
!#                                                                           #
!# Copyright (C) 2006 - 2010                                                 #
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
! References:
! [1] Bondi, H.: On spherically symmetrical accretion,
!     Mon. Not. Roy. Astron. Soc., 112 (1951)
!     ADS link: http://adsabs.harvard.edu/abs/1952MNRAS.112..195B
! [2] Padmanabhan, T.:Theoretical Astrophysics, Vol. I: Astrophysical
!     Processes, Cambridge University Press (2000), Chapter 8.9.2
!----------------------------------------------------------------------------!

!**************************************!
!* IMPORTANT:                         *!
!* - compile with autodouble          *!
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
  REAL, PARAMETER    :: MSUN = 1.989E+30   ! solar mass [kg]
  ! simulation parameters
  REAL, PARAMETER    :: TSIM    = 10.0     ! simulation time [TAU] (free fall)
  REAL, PARAMETER    :: ACCMASS = 1.0*MSUN ! mass of the accreting object
  REAL, PARAMETER    :: GAMMA   = 1.4      ! ratio of specific heats 
  ! boundary conditions
  REAL, PARAMETER    :: RHOINF  = 1.0E-20  ! density at infinity [kg/m^3]
  REAL, PARAMETER    :: CSINF   = 1.0E+04  ! sound speed at infinity [m/s]
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = POLAR       ! geometry
!!$  INTEGER, PARAMETER :: MGEO = LOGPOLAR
!!$  INTEGER, PARAMETER :: MGEO = TANPOLAR
!!$  INTEGER, PARAMETER :: MGEO = SINHPOLAR
  INTEGER, PARAMETER :: XRES = 50          ! x-resolution
  INTEGER, PARAMETER :: YRES = 1           ! y-resolution
  REAL, PARAMETER    :: RIN  = 0.1         ! inner/outer radii in terms of
  REAL, PARAMETER    :: ROUT = 2.0         !   the Bondi radius RB, ROUT > 1
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 10          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'bondi2d' 
  ! some derives quandities
  REAL               :: RB                 ! Bondi radius
  REAL               :: TAU                ! free fall time scale
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
    REAL              :: x1,x2,scale
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!
    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER2D, &
         gamma   = GAMMA, &                 ! ratio of specific heats        !
         dpmax   = 1.0)                     ! for advanced time step control !

    ! derived constants
    RB  = Physics%Constants%GN * ACCMASS / CSINF**2  ! bondi radius [m]      !
    TAU = RB / CSINF                        ! free fall time scale [s]       !

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)                 ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         variables = CONSERVATIVE, &        ! vars. to use for reconstruction!
         limiter   = MONOCENT, &            ! one of: minmod, monocent,...   !
         theta     = 1.2)                   ! optional parameter for limiter !

    ! mesh settings
    SELECT CASE(MGEO)
    CASE(POLAR)
       x1 = RIN * RB
       x2 = ROUT * RB
       scale = 1.0
    CASE(LOGPOLAR)
       x1 = LOG(RIN)
       x2 = LOG(ROUT)
       scale = RB
    CASE(TANPOLAR)
       x1 = ATAN(RIN)
       x2 = ATAN(ROUT)
       scale = RB
    CASE(SINHPOLAR)
       x1 = LOG(RIN+SQRT(1.0+RIN*RIN))  ! = ASINH(RIN))
       x2 = LOG(ROUT+SQRT(1.0+ROUT*ROUT))
       scale = RB
    CASE DEFAULT
       CALL Error(Physics,"InitProgram","mesh geometry not supported for 2D Bondi accretion")
    END SELECT
    CALL InitMesh(Mesh,Fluxes, &
         geometry = MGEO, &
             inum = XRES, &
             jnum = YRES, &
             xmin = x1, &
             xmax = x2, &
             ymin = 0.0, &
             ymax = 2*PI, &
           gparam = scale)

    ! source term due to a point mass
    CALL InitSources(Physics%sources,Mesh,Fluxes,Physics,Timedisc%boundary, &
         stype  = POINTMASS, &            ! grav. accel. of a point mass     !
           mass = ACCMASS)                ! mass of the accreting object[kg] !
    Physics%sources%outbound = 0          ! disable accretion

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
         stoptime = TSIM * TAU, &
         dtlimit  = 1.0E-6 * TAU, &
         maxiter  = 1000000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Fluxes,Timedisc)

    ! initialize log input/output
!!$    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,&
!!$         fileformat = BINARY, &
!!$         filename   = TRIM(ODIR) // TRIM(OFNAME) // 'log', &
!!$         filecycles = 1)                  ! just one log file

    ! initialize data input/output
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc, &
         fileformat = GNUPLOT, &
         filecycles = 0, &                ! all time steps in one file
         filename   = TRIM(ODIR) // TRIM(OFNAME), &
         count      = ONUM)
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
    REAL              :: r,rho,vr,cs2
    INTEGER           :: i,j
    CHARACTER(LEN=64) :: info_str
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
    IF (GetType(Timedisc%Boundary(EAST)).EQ.FIXED) THEN
       DO j=Mesh%JMIN,Mesh%JMAX
          DO i=1,2
             r = SQRT(Mesh%bccart(Mesh%IMAX+i,j,1)**2+Mesh%bccart(Mesh%IMAX+i,j,2)**2)
             CALL bondi(r/RB,GAMMA,RHOINF,CSINF,rho,vr)
             cs2 = CSINF**2 * (rho/RHOINF)**(GAMMA-1.0)
             ! set boundary data to either primitive or conservative values
             ! depending on the reconstruction
             Timedisc%Boundary(EAST)%data(i,j,Physics%DENSITY)   = rho
             Timedisc%Boundary(EAST)%data(i,j,Physics%XVELOCITY) = vr
             Timedisc%Boundary(EAST)%data(i,j,Physics%YVELOCITY) = 0.
             Timedisc%Boundary(EAST)%data(i,j,Physics%PRESSURE)  = rho * cs2 / GAMMA
          END DO
          ! this tells the boundary routine which values to fix (.TRUE.)
          ! and which to extrapolate (.FALSE.)
          Timedisc%Boundary(EAST)%fixed(j,:) = (/ .TRUE., .FALSE., .TRUE., .TRUE. /)
      END DO
    END IF

    CALL Info(Mesh," DATA-----> initial condition: 2D Bondi accretion")
    WRITE(info_str,"(ES9.3)") RB
    CALL Info(Mesh, "                               " // "Bondi radius:       " &
         // TRIM(info_str) // " m")
    WRITE(info_str,"(ES9.3)") TAU
    CALL Info(Mesh, "                               " // "Free fall time:     " &
         // TRIM(info_str) // " s")
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
