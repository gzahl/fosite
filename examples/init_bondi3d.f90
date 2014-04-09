!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: init_bondi3d.f90                                                  #
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
! Program and data initialization for 3D Bondi accretion
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
  USE output_generic
  USE logio_generic
  USE sources_generic
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

    ! mesh settings
    CALL InitMesh(Mesh,Fluxes, &
         geometry = SPHERICAL, &
            inum = 50, &
            jnum = 6, &
            xmin = 2.0E+18, &
            xmax = 4.77E+19, &
            ymin = 0.0, &
            ymax = PI)

    ! source term due to a point mass
    CALL InitSources(Physics%sources,Mesh,Fluxes,Physics, &
         stype  = POINTMASS, &      ! grav. accel. of a point mass   !
         sparam = 1.0E+38)          ! mass [kg]                      !

    ! boundary conditions
    CALL InitBoundary(Mesh%boundary,Mesh,Physics,EXTRAPOLATION,WEST)
    CALL InitBoundary(Mesh%boundary,Mesh,Physics,FIXED,EAST)
    CALL InitBoundaryData(Mesh%boundary,Mesh,Fluxes,Physics,EAST)
    CALL InitBoundary(Mesh%boundary,Mesh,Physics,AXIS,SOUTH)
    CALL InitBoundary(Mesh%boundary,Mesh,Physics,AXIS,NORTH)

    ! time discretization settings
    CALL InitTimedisc(Timedisc,Mesh,Physics,&
         method   = MODIFIED_EULER, &
         order    = 3, &
         cfl      = 0.4, &
         stoptime = 5.0E+16, &
         dtlimit  = 1.0E+8, &
         maxiter  = 100000)

    ! set initial condition
    CALL InitData(Mesh,Physics,Timedisc%pvar,Timedisc%cvar)

    ! initialize log input/output
    CALL InitLogio(Logio,Mesh,Physics,Timedisc, &
         logformat = NOLOG, &
         filename  = "bondi3d.log", &
         logdt     = 300)

    ! set parameters for data output
    CALL InitOutput(Output,Mesh,Physics,Timedisc,&
         filetype  = GNUPLOT, &
         filename  = "bondi3d.dat", &
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
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(OUT)       :: pvar,cvar
    !------------------------------------------------------------------------!

    ! ensure that all cells have pressure and density != 0 and vanishing velocities
    pvar(:,:,1)   = 1.0E-20
    pvar(:,:,2:4) = 0.
    pvar(:,:,5)   = 1.0E-12

    CALL Convert2Conservative(Physics,Mesh,pvar,cvar)

    PRINT "(A,A)", " DATA-----> initial condition: ", &
         "3D Bondi accretion"
  END SUBROUTINE InitData


  SUBROUTINE InitBoundaryData(Boundary,Mesh,Fluxes,Physics,direction)
    USE roots, ONLY : RtNewtBisec
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), DIMENSION(4) :: Boundary
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    INTEGER           :: direction
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(SOURCES_TYP), POINTER :: srcptr
    REAL, DIMENSION(2,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) :: pdata
    REAL              :: rhoinf,pinf,cinf
    REAL              :: gamma,gp1,gm1
    REAL              :: rb,rc,rhoc,vc
    REAL              :: lambda,chi,psi
    REAL              :: r,gr
    REAL              :: xacc
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Fluxes,Physics,direction
    INTENT(INOUT)     :: Boundary,Mesh
    !------------------------------------------------------------------------!

    ! calculate exact solution for 2D bondi problem
    ! and initialize the boundary data for subsonic inflow
    rhoinf = 1.0E-20
    pinf   = 1.0E-12

    ! some constants
    gamma  = Physics%gamma
    gp1    = gamma + 1.
    gm1    = gamma - 1.

    ! sound speed @ infitiy
    cinf = SQRT(gamma * pinf / rhoinf)
    ! critical radius
    rc = 0.25 * (5.-3.*gamma)
    ! critical density
    rhoc = (0.5/rc)**(1./gm1)
    ! critical velocity
    vc = SQRT(0.5/rc)
    ! critical dimensionless accretion rate
    lambda = rc**2 * rhoc * vc

    ! bondi radius
    rb = 0.
    srcptr => GetSourcesPointer(Physics%sources,POINTMASS)
    IF (ASSOCIATED(srcptr)) THEN
       rb = Physics%Constants%GN*srcptr%mass / cinf**2
    ELSE
       PRINT *, "ERROR in InitBoundaryData: no gravitational point sources defined"
       STOP
    END IF

    ! accuracy
    xacc = 1.0E-6

    ! calculate boundary solution for y=ymin..ymax at xmax
    DO j=Mesh%JMIN,Mesh%JMAX
       DO i=1,2
          r = Mesh%bcenter(Mesh%IMAX+i,j,1) / rb
          chi = r**2 / lambda
          gr  = chi**(2.*gm1/gp1) * (1./gm1 + 1./r)
          ! Newton-Raphson to solve Bondis equations for psi
          IF (r.LT.rc) THEN
             psi = RtNewtBisec(funcd,1.0,gr,gm1,gr,xacc)
          ELSE
             psi = RtNewtBisec(funcd,1.0E-10,1.0,gm1,gr,xacc)
          END IF
          pdata(i,j,1) = rhoinf * chi**(-2./gp1) / psi       ! density
          pdata(i,j,2) = -cinf * chi**(-gm1/gp1) * psi       ! radial velocity
          pdata(i,j,3) = 0.                                  ! azimuthal velocity
          pdata(i,j,4) = 0.                                  ! spec. angular momentum
          pdata(i,j,5) = pinf * (pdata(i,j,1)/rhoinf)**gamma ! pressure
       END DO
    END DO

    ! set boundary data to either primitive or conservative values
    ! depending on the reconstruction
    IF (PrimRecon(Fluxes%reconstruction).EQV.PRIMITIVE) THEN
       Boundary(direction)%data = pdata
    ELSE
       Boundary(direction)%data(:,:,1) = pdata(:,:,1)              ! density
       Boundary(direction)%data(:,:,2) = pdata(:,:,1)*pdata(:,:,2) ! radial momentum
       Boundary(direction)%data(:,:,3) = pdata(:,:,1)*pdata(:,:,3) ! azimuthal momentum
       Boundary(direction)%data(:,:,4) = pdata(:,:,1)*pdata(:,:,4) ! angular momentum
       Boundary(direction)%data(:,:,5) = pdata(:,:,5)/gm1 &        ! total energy density
            + 0.5*pdata(:,:,1)*pdata(:,:,2)**2
    END IF

    ! this tells the boundary routine which values to fix (.TRUE.)
    ! and which to extrapolate (.FALSE.)
    Boundary(direction)%fixed = (/ .TRUE., .FALSE., .TRUE., .TRUE., .TRUE. /)
  END SUBROUTINE InitBoundaryData


  ! for exact Bondi solution at the outer boundary
  SUBROUTINE funcd(y,gm1,gx,fy,dfy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: y,gm1,gx
    REAL, INTENT(OUT) :: fy,dfy
    !------------------------------------------------------------------------!
    fy  = 0.5*y*y + y**(-gm1) / gm1 - gx
    dfy = y - y**(-gm1-1.)
!    PRINT '(5(ES14.6))',y,fy,dfy,gm1,gx
  END SUBROUTINE funcd

END MODULE Init
