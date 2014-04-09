!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_viscosity.f90                                             #
!#                                                                           #
!# Copyright (C) 2008 Bjoern Sperling <sperling@astrophysik.uni-kiel.de>     #
!#                    Tobias Illenseer <tillense@astrophysik.uni-kiel.de>    #
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
!# along with Mesh program; if not, write to the Free Software               #
!# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
!#                                                                           #
!#############################################################################
!----------------------------------------------------------------------------!
! assumption: Newton Fluid (isotropic)
! 
!----------------------------------------------------------------------------!
MODULE sources_viscosity
  USE common_types, ONLY : Common_TYP, InitCommon
  USE sources_pointmass
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  REAL, PARAMETER :: TINY = 1.0E-30              ! to avoid division by 0    !
  CHARACTER(LEN=32), PARAMETER :: source_name = "viscosity of newton fluid"
  ! flags for viscosity model
  INTEGER, PARAMETER :: MOLECULAR = 1
  INTEGER, PARAMETER :: ALPHA     = 2
  INTEGER, PARAMETER :: BETA      = 3
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! constants
       MOLECULAR, ALPHA, BETA, &
       ! methods
       InitSources_viscosity, &
       ExternalSources_viscosity, &
       CalcTimestep_viscosity, &
       CloseSources_viscosity
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources_viscosity(this,Mesh,Physics,Fluxes,stype,model,dynconst, &
       bulkconst,cvis,pmsrc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this, pmsrc
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    INTEGER           :: stype
    INTEGER           :: model
    REAL              :: dynconst,bulkconst,cvis
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    INTEGER           :: err
    REAL              :: r,a,RS
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: cart
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Fluxes,stype,model,dynconst,bulkconst,cvis
    !------------------------------------------------------------------------!
    CALL InitSources(this,stype,source_name)

    ! check mesh
    IF (GetType(Fluxes).NE.MIDPOINT) &
         CALL Error(this,"InitSources_viscosity","only midpoint rule is currently supported")

    ALLOCATE(this%dynvis(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%bulkvis(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btxx(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btyy(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btzz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btxy(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btxz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%btyz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%ftxx(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4), &
         this%ftyy(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4), &
         this%ftzz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4), &
         this%ftxy(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4), &
         this%ftxz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4), &
         this%ftyz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4), &
         STAT=err)
    IF (err.NE.0) CALL Error(this,"InitSource_viscosity", "Unable to allocate memory.")

    ! set viscosity constants
    this%dynconst  = dynconst
    this%bulkconst = bulkconst

    ! do the model specific initialization
    SELECT CASE(model)
    CASE(MOLECULAR)
       CALL InitCommon(this%viscosity,MOLECULAR,"molecular")
    CASE(ALPHA)
       CALL InitCommon(this%viscosity,ALPHA,"alpha")
       ! check for point mass source term
       IF (ASSOCIATED(pmsrc).EQV..FALSE.) &
            CALL Error(this,"InitSources_viscosity","no point mass source term defined")
       this%mass = pmsrc%mass
       ! Schwarzschild radius
       RS = 2*Physics%constants%GN * this%mass / Physics%constants%C**2
       ! check type of the potential
       IF (GetType(pmsrc%potential).EQ.WIITA) THEN
          ! post-newtonian Paczinski-Wiita potential
          a = RS
       ELSE
          ! newtonian potential
          a = 0.0
       END IF
       ! store keplerian angular velocity
       ALLOCATE(this%omega(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX),STAT=err)
       IF (err.NE.0) CALL Error(this,"InitSource_viscosity", "Unable to allocate memory.")
       CALL Convert2Cartesian(Mesh%geometry,Mesh%bcenter,cart)
       DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
             ! distance to the origin of the coordinate system
             r  = SQRT(cart(i,j,1)**2 + cart(i,j,2)**2)
             this%omega(i,j) = SQRT(0.5*RS/r) * Physics%constants%C / (r-a)
          END DO
       END DO
    CASE(BETA)
       CALL InitCommon(this%viscosity,BETA,"beta")
    END SELECT

    ! set viscous courant number
    this%cvis = cvis
    ! initialize viscosity arrays
    this%dynvis(:,:)  = dynconst
    this%bulkvis(:,:) = bulkconst
    ! set stress tensor arrays to zero
    this%btxx(:,:) = 0.0
    this%btyy(:,:) = 0.0
    this%btzz(:,:) = 0.0
    this%btxy(:,:) = 0.0
    this%btxz(:,:) = 0.0
    this%btyz(:,:) = 0.0
    this%ftxx(:,:,:) = 0.0
    this%ftyy(:,:,:) = 0.0
    this%ftzz(:,:,:) = 0.0
    this%ftxy(:,:,:) = 0.0
    this%ftxz(:,:,:) = 0.0
    this%ftyz(:,:,:) = 0.0    
  END SUBROUTINE InitSources_viscosity


  SUBROUTINE UpdateViscosity(this,Mesh,Physics,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,pvar
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!

    SELECT CASE(GetType(this%viscosity))
    CASE(MOLECULAR)
       ! do nothing, since it has already been initialized in InitSources_viscosity
    CASE(ALPHA)
       this%dynvis(:,:) = etafkt_alpha(this%dynconst,Physics%gamma,this%omega(:,:),&
            pvar(:,:,Physics%PRESSURE))
    END SELECT

  CONTAINS

    ELEMENTAL FUNCTION etafkt(x,y,rho,vx,vy,vz,p) RESULT(eta)
      IMPLICIT NONE
      !----------------------------------------------------------------------!
      REAL, INTENT(IN) :: x, y, rho, vx, vy, vz, p
      REAL :: eta
      !----------------------------------------------------------------------!
      eta = 0.
    END FUNCTION etafkt
    
    ELEMENTAL FUNCTION zetafkt(x,y,rho,vx,vy,vz,p) RESULT(zeta)
      IMPLICIT NONE
      !----------------------------------------------------------------------!
      REAL, INTENT(IN) :: x, y, rho, vx, vy, vz, p
      REAL :: zeta
      !----------------------------------------------------------------------!
      zeta = 0.
    END FUNCTION zetafkt
    
    ! Shakura-Sunyaev type alpha viscosity
    ELEMENTAL FUNCTION etafkt_alpha(alpha,gamma,omega,p) RESULT(eta)
      IMPLICIT NONE
      !----------------------------------------------------------------------!
      REAL, INTENT(IN) :: alpha,gamma,omega,p
      REAL :: eta
      !----------------------------------------------------------------------!
      eta = alpha * SQRT(gamma) * p / omega
    END FUNCTION etafkt_alpha

  END SUBROUTINE UpdateViscosity


  SUBROUTINE ExternalSources_viscosity(this,Mesh,Physics,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: cvar,pvar,sterm
    !------------------------------------------------------------------------!
    INTEGER,SAVE      :: count=0
    INTEGER           :: error
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,cvar,pvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    CALL UpdateViscosity(this,Mesh,Physics,pvar)
    CALL ViscositySources(Physics,Mesh,this,pvar,cvar,sterm)
!!$    count = count+1
!!$    IF (count.GT.1) THEN
!!$       PRINT '(I1,A,2(I4),A,3(ES15.7))', GetRank(this), " ***", Mesh%IMIN,Mesh%JMIN, &
!!$            ACHAR(10), &
!!$           ( pvar(Mesh%IMIN+1,Mesh%JMIN-1,Physics%YVELOCITY) &
!!$           - pvar(Mesh%IMIN-1,Mesh%JMIN-1,Physics%YVELOCITY)) / Mesh%dlx(Mesh%IMIN,Mesh%JMIN-1), &
!!$           ( pvar(Mesh%IMIN,Mesh%JMIN,Physics%XVELOCITY) &
!!$           - pvar(Mesh%IMIN,Mesh%JMIN-2,Physics%XVELOCITY)) / Mesh%dly(Mesh%IMIN,Mesh%JMIN-1)
!!$       PRINT '(3(ES15.7),A,3(ES15.7))', sterm(Mesh%IMIN,Mesh%JMIN,Physics%XMOMENTUM), &
!!$            sterm(Mesh%IMIN,Mesh%JMIN,Physics%YMOMENTUM), &
!!$            sterm(Mesh%IMIN,Mesh%JMIN,Physics%ENERGY), &
!!$            ACHAR(10), &
!!$            sterm(Mesh%IMAX,Mesh%JMAX,Physics%XMOMENTUM), &
!!$            sterm(Mesh%IMAX,Mesh%JMAX,Physics%YMOMENTUM), &
!!$            sterm(Mesh%IMAX,Mesh%JMAX,Physics%ENERGY)
!!$    END IF
!!$    IF (count.EQ.3) THEN
!!$       CALL MPI_Finalize(error)
!!$       STOP
!!$    END IF
  END SUBROUTINE ExternalSources_viscosity

  
  SUBROUTINE CalcTimestep_viscosity(this,Mesh,Physics,pvar,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    REAL              :: dt
    !------------------------------------------------------------------------!
    REAL              :: invdt_x, invdt_y
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: dt
    !------------------------------------------------------------------------!
    ! FIXME: this call necessary, but its inefficient, since we call it
    !        again, when computing the viscous source term
    CALL UpdateViscosity(this,Mesh,Physics,pvar)

    ! using this%btxx temporary for kinematic viscosity
    this%btxx(:,:) = this%dynvis(:,:) / pvar(:,:,Physics%DENSITY)

    ! x-direction
    IF (Mesh%INUM.GT.1) THEN
       invdt_x = MAXVAL(this%btxx(:,:) / Mesh%dlx(:,:)**2)
    ELSE
       ! set to zero, i.e. no limit in x-direction due to diffusion
       invdt_x = 0.0
    END IF
    ! y-direction
    IF (Mesh%JNUM.GT.1) THEN
       invdt_y = MAXVAL(this%btxx(:,:) / Mesh%dly(:,:)**2)
    ELSE
       ! set to zero, i.e. no limit in x-direction due to diffusion
       invdt_y = 0.0
    END IF
    ! largest time step due to diffusion
    dt = this%cvis / MAX(invdt_x, invdt_y)
  END SUBROUTINE CalcTimestep_viscosity


  SUBROUTINE CloseSources_viscosity(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    IF (GetType(this%viscosity).EQ.ALPHA) DEALLOCATE(this%omega)
    DEALLOCATE(this%dynvis,this%bulkvis, &
               this%btxx,this%btyy,this%btzz,this%btxy,this%btxz,this%btyz,&
               this%ftxx,this%ftyy,this%ftzz,this%ftxy,this%ftxz,this%ftyz)
  END SUBROUTINE CloseSources_viscosity

END MODULE sources_viscosity
