!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_bondi3d.f90                                                  #
!#                                                                           #
!# Copyright (C) 2006-2011                                                   #
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
! Program and data initialization for 3D Bondi accretion
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
  REAL, PARAMETER    :: MSUN = 1.989E+30   ! solar mass [kg]                 !
  ! simulation parameters
  REAL, PARAMETER    :: TSIM = 10.0        ! simulation time [TAU] (free fall)
  REAL, PARAMETER    :: ACCMASS = 1.0*MSUN ! mass of the accreting object    !
  REAL, PARAMETER    :: GAMMA = 1.4        ! ratio of specific heats         !
  ! boundary conditions
  REAL, PARAMETER    :: RHOINF = 3.351E-17 ! density at infinity [kg/m^3]    !
  REAL, PARAMETER    :: CSINF = 537.0      ! sound speed at infinity [m/s]   !
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = SPHERICAL   ! geometry of the mesh            !
!!$  INTEGER, PARAMETER :: MGEO = OBLATE_SPHEROIDAL
!!$  INTEGER, PARAMETER :: MGEO = CYLINDRICAL
!!$  INTEGER, PARAMETER :: MGEO = TANCYLINDRICAL
!!$  INTEGER, PARAMETER :: MGEO = SINHSPHERICAL
  INTEGER, PARAMETER :: XRES = 100         ! x-resolution                    !
  INTEGER, PARAMETER :: YRES = 1           ! y-resolution                    !
  REAL, PARAMETER    :: RIN  = 0.1         ! inner/outer radii in terms of   !
  REAL, PARAMETER    :: ROUT = 2.0         !   the Bondi radius RB, ROUT > 1 !
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 10          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'bondi3d'
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
    REAL              :: x1,x2,y1,y2
    INTEGER           :: bc(4)
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile
    !------------------------------------------------------------------------!
    ! physics settings
    CALL InitPhysics(Physics, &
         problem = EULER3D_ROTSYM, &
         gamma   = GAMMA, &                 ! ratio of specific heats        !
         dpmax   = 1.0E+3)                  ! for advanced time step control !

    ! derived constants
    RB  = Physics%Constants%GN * ACCMASS / CSINF**2  ! bondi radius [m]      !
    TAU = RB / CSINF                        ! free fall time scale [s]       !

    ! numerical scheme for flux calculation
    CALL InitFluxes(Fluxes, &
         scheme = MIDPOINT)                 ! quadrature rule                !

    ! reconstruction method
    CALL InitReconstruction(Fluxes%reconstruction, &
         order     = LINEAR, &
         variables = PRIMITIVE, &           ! vars. to use for reconstruction!
         limiter   = MONOCENT, &            ! one of: minmod, monocent,...   !
!!$         limiter   = SUPERBEE, &            ! for better entropy conservation !
         theta     = 1.2)                   ! optional parameter for limiter !

    ! geometry dependent setttings
    SELECT CASE(MGEO)
    CASE(SPHERICAL)
       x1 = RIN * RB
       x2 = ROUT * RB
       y1 = 0.0
       y2 = PI
       bc(WEST)  = EXTRAPOLATION
       bc(EAST)  = FARFIELD
       bc(SOUTH) = AXIS
       bc(NORTH) = AXIS
    CASE(SINHSPHERICAL)
       x1 = LOG(RIN+SQRT(RIN**2+1.0))   ! = ASINH(RIN)
       x2 = LOG(ROUT+SQRT(ROUT**2+1.0)) ! = ASINH(ROUT)
       y1 = 0.0
       y2 = PI
       bc(WEST)  = EXTRAPOLATION
       bc(EAST)  = FARFIELD
       bc(SOUTH) = AXIS
       bc(NORTH) = AXIS
    CASE(CYLINDRICAL)
       x1 = -ROUT * RB
       x2 = ROUT * RB
       y1 = RIN * RB
       y2 = ROUT * RB
       bc(WEST)  = FARFIELD
       bc(EAST)  = FARFIELD
       bc(SOUTH) = NO_GRADIENTS
       bc(NORTH) = FARFIELD
    CASE(OBLATE_SPHEROIDAL)
       x1 = LOG(RIN+SQRT(RIN**2+1.0))   ! = ASINH(RIN)
       x2 = LOG(ROUT+SQRT(ROUT**2-1.0)) ! = ACOSH(ROUT)
       y1 = -0.5*PI
       y2 = 0.5*PI
       bc(WEST)  = FARFIELD
       bc(EAST)  = FARFIELD
       bc(SOUTH) = AXIS
       bc(NORTH) = AXIS
    CASE(TANCYLINDRICAL)
       x1 = ATAN(-ROUT)
       x2 = ATAN(ROUT)
       y1 = RIN * RB
       y2 = ROUT * RB
       bc(WEST)  = FARFIELD
       bc(EAST)  = FARFIELD
       bc(SOUTH) = FARFIELD
       bc(NORTH) = FARFIELD
    CASE DEFAULT
       CALL Error(Physics,"InitProgram","mesh geometry not supported for Bondi accretion")
    END SELECT

    ! mesh settings
    CALL InitMesh(Mesh,Fluxes, &
         geometry = MGEO, &
             inum = XRES, &
             jnum = YRES, &
             xmin = x1, &
             xmax = x2, &
             ymin = y1, &
             ymax = y2, &
           gparam = RB)

    ! boundary conditions
    CALL InitBoundary(Timedisc%boundary,Mesh,Physics, &
         western  = bc(WEST), &
         eastern  = bc(EAST), &
         southern = bc(SOUTH), &
         northern = bc(NORTH))

    ! source term due to a point mass
    CALL InitSources(Physics%sources,Mesh,Fluxes,Physics,Timedisc%boundary, &
          stype = POINTMASS, &            ! grav. accel. of a point mass     !
           mass = ACCMASS)                ! mass of the accreting object[kg] !
    Physics%sources%outbound = 0          ! disable accretion

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
!!$         filecycles = 1)

    ! initialize data input/output
    CALL InitFileIO(Datafile,Mesh,Physics,Timedisc, &
         fileformat = GNUPLOT, &
         filecycles = 0, &
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
    REAL              :: r,x,y,rho,vr,cs2
    REAL, DIMENSION(Mesh%GNUM,Mesh%JMIN:Mesh%JMAX,2) :: v1,bc1
    REAL, DIMENSION(Mesh%IMIN:Mesh%IMAX,Mesh%GNUM,2) :: v2,bc2
    INTEGER           :: i,j
    CHARACTER(LEN=64) :: info_str
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Fluxes
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! initial condition: use data at infinity
    Timedisc%pvar(:,:,Physics%DENSITY)   = RHOINF
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%ZVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%PRESSURE)  = RHOINF * CSINF**2 / GAMMA
    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)

    ! boundary conditions
    ! outflow condition, this is only for transition; in the end the flow
    ! becomes supersonic and all variables are extrapolated from the interior
    IF (GetType(Timedisc%Boundary(WEST)).EQ.FARFIELD) THEN
       DO j=Mesh%JMIN,Mesh%JMAX
          DO i=1,Mesh%GNUM
             x = Mesh%bccart(Mesh%IMIN-i,j,1)
             y = Mesh%bccart(Mesh%IMIN-i,j,2)
             r = SQRT(x**2+y**2)
             CALL bondi(r/RB,GAMMA,RHOINF,CSINF,rho,vr)
             cs2 = CSINF**2 * (rho/RHOINF)**(GAMMA-1.0)
             v1(i,j,1) = x/r * vr
             v1(i,j,2) = y/r * vr
             bc1(i,j,:) = Mesh%bcenter(Mesh%IMIN-i,j,:)
             Timedisc%Boundary(WEST)%data(i,j,Physics%DENSITY)   = rho
             Timedisc%Boundary(WEST)%data(i,j,Physics%ZVELOCITY) = 0.
             Timedisc%Boundary(WEST)%data(i,j,Physics%PRESSURE)  = rho * cs2 / GAMMA
          END DO
       END DO
       CALL Convert2Curvilinear(Mesh%geometry,bc1,v1, &
            Timedisc%Boundary(WEST)%data(:,:,Physics%XVELOCITY:Physics%YVELOCITY))
    END IF
    ! subsonic inflow according to Bondi's solution
    ! calculate Bondi solution for y=ymin..ymax at xmax
    IF (GetType(Timedisc%Boundary(EAST)).EQ.FARFIELD) THEN
       DO j=Mesh%JMIN,Mesh%JMAX
          DO i=1,Mesh%GNUM
             x = Mesh%bccart(Mesh%IMAX+i,j,1)
             y = Mesh%bccart(Mesh%IMAX+i,j,2)
             r = SQRT(x**2+y**2)
             CALL bondi(r/RB,GAMMA,RHOINF,CSINF,rho,vr)
             cs2 = CSINF**2 * (rho/RHOINF)**(GAMMA-1.0)
             v1(i,j,1) = x/r * vr
             v1(i,j,2) = y/r * vr
             bc1(i,j,:) = Mesh%bcenter(Mesh%IMAX+i,j,:)
             Timedisc%Boundary(EAST)%data(i,j,Physics%DENSITY)   = rho
             Timedisc%Boundary(EAST)%data(i,j,Physics%ZVELOCITY) = 0.
             Timedisc%Boundary(EAST)%data(i,j,Physics%PRESSURE)  = rho * cs2 / GAMMA
          END DO
       END DO
       CALL Convert2Curvilinear(Mesh%geometry,bc1,v1, &
            Timedisc%Boundary(EAST)%data(:,:,Physics%XVELOCITY:Physics%YVELOCITY))
    END IF
    IF (GetType(Timedisc%Boundary(SOUTH)).EQ.FARFIELD) THEN
       DO j=1,Mesh%GNUM
          DO i=Mesh%IMIN,Mesh%IMAX
             x = Mesh%bccart(i,Mesh%JMIN-j,1)
             y = Mesh%bccart(i,Mesh%JMIN-j,2)
             r = SQRT(x**2+y**2)
             CALL bondi(r/RB,GAMMA,RHOINF,CSINF,rho,vr)
             cs2 = CSINF**2 * (rho/RHOINF)**(GAMMA-1.0)
             v2(i,j,1) = x/r * vr
             v2(i,j,2) = y/r * vr
             bc2(i,j,:) = Mesh%bcenter(i,Mesh%JMIN-j,:)
             Timedisc%Boundary(SOUTH)%data(i,j,Physics%DENSITY)   = rho
             Timedisc%Boundary(SOUTH)%data(i,j,Physics%ZVELOCITY) = 0.
             Timedisc%Boundary(SOUTH)%data(i,j,Physics%PRESSURE)  = rho * cs2 / GAMMA
          END DO
       END DO
       CALL Convert2Curvilinear(Mesh%geometry,bc2,v2, &
            Timedisc%Boundary(SOUTH)%data(:,:,Physics%XVELOCITY:Physics%YVELOCITY))
    END IF
    IF (GetType(Timedisc%Boundary(NORTH)).EQ.FARFIELD) THEN
       DO j=1,Mesh%GNUM
          DO i=Mesh%IMIN,Mesh%IMAX
             x = Mesh%bccart(i,Mesh%JMAX+j,1)
             y = Mesh%bccart(i,Mesh%JMAX+j,2)
             r = SQRT(x**2+y**2)
             CALL bondi(r/RB,GAMMA,RHOINF,CSINF,rho,vr)
             cs2 = CSINF**2 * (rho/RHOINF)**(GAMMA-1.0)
             v2(i,j,1) = x/r * vr
             v2(i,j,2) = y/r * vr
             bc2(i,j,:) = Mesh%bcenter(i,Mesh%JMAX+j,:)
             Timedisc%Boundary(NORTH)%data(i,j,Physics%DENSITY)   = rho
             Timedisc%Boundary(NORTH)%data(i,j,Physics%ZVELOCITY) = 0.
             Timedisc%Boundary(NORTH)%data(i,j,Physics%PRESSURE)  = rho * cs2 / GAMMA
          END DO
       END DO
       CALL Convert2Curvilinear(Mesh%geometry,bc2,v2, &
            Timedisc%Boundary(NORTH)%data(:,:,Physics%XVELOCITY:Physics%YVELOCITY))
    END IF
    CALL Info(Mesh," DATA-----> initial condition: 3D Bondi accretion")
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
    ! computes the Bondi solution for spherically symmetric accretion        !
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
    REAL :: gp1,gm1,g35,rc,chi,lambda,psi,gr
    COMMON /funcd_parameter/ gm1, gr
    !------------------------------------------------------------------------!
    ! for convenience
    gm1 = gamma - 1.0
    gp1 = gamma + 1.0
    g35 = 0.5 * (5.0-3.0*gamma)

    ! critical radius
    rc = 0.5 * g35
    ! critical dimensionless accretion rate
    lambda = 0.25 * g35**(-g35/gm1)

    ! Newton-Raphson to solve Bondis equations for psi
    chi = r**2 / lambda
    gr  = chi**(2.*gm1/gp1) * (1./gm1 + 1./r)
    IF (r.LT.rc) THEN
       psi = GetRoot_newton(funcd,1.0,gr)
    ELSE
       psi = GetRoot_newton(funcd,1.0E-6,1.0)
    END IF
    
    ! return values
    rho = rhoinf * chi**(-2./gp1) / psi        ! density
    vr  = -csinf * chi**(-gm1/gp1) * psi       ! radial velocity
  END SUBROUTINE bondi


  ! for exact Bondi solution at the outer boundary
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
