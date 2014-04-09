!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_viscosity.f90                                             #
!#                                                                           #
!# Copyright (C) 2008-2009                                                   #
!# Bjoern Sperling <sperling@astrophysik.uni-kiel.de>                        #
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
!# along with Mesh program; if not, write to the Free Software               #
!# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
!#                                                                           #
!#############################################################################
!----------------------------------------------------------------------------!
! viscosity of Newtonian fluid
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
  CHARACTER(LEN=32), PARAMETER :: source_name = "viscosity of Newtonian fluid"
  ! flags for viscosity model
  INTEGER, PARAMETER :: MOLECULAR = 1
  INTEGER, PARAMETER :: ALPHA     = 2
  INTEGER, PARAMETER :: BETA      = 3
  INTEGER, PARAMETER :: PRINGLE   = 4
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! constants
       MOLECULAR, ALPHA, BETA, PRINGLE, &
       ! methods
       InitSources_viscosity, &
       ExternalSources_viscosity, &
       CalcTimestep_viscosity, &
       CloseSources_viscosity
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources_viscosity(this,Mesh,Physics,Fluxes,stype,model,dynconst, &
       bulkconst,cvis)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    INTEGER           :: stype
    INTEGER           :: model
    REAL              :: dynconst,bulkconst,cvis
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Fluxes,stype,model,dynconst,bulkconst,cvis
    !------------------------------------------------------------------------!
    IF (.NOT.Initialized(Fluxes)) &
         CALL Error(this,"InitSources_viscosity","fluxes module uninitialized")

    CALL InitSources(this,stype,source_name)

    ! check mesh
    IF (GetType(Fluxes).NE.MIDPOINT) &
         CALL Error(this,"InitSources_viscosity","only midpoint rule is currently supported")

    ALLOCATE(this%dynvis(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%kinvis(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
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
    IF (err.NE.0) CALL Error(this,"InitSources_viscosity", "Unable to allocate memory.")

    ! set viscosity constants
    this%dynconst  = dynconst
    this%bulkconst = bulkconst

    ! do the model specific initialization
    SELECT CASE(model)
    CASE(MOLECULAR)
       CALL InitCommon(this%viscosity,MOLECULAR,"molecular")
    CASE(ALPHA)
!!$       IF (Physics%PRESSURE.EQ.0) &
!!$            CALL Error(this,"InitSources_viscosity", "Physics not supported for alpha-viscosity")
       SELECT CASE(GetType(Mesh%geometry))
       CASE(POLAR,LOGPOLAR,TANPOLAR,SINHPOLAR,CYLINDRICAL,TANCYLINDRICAL,SPHERICAL,OBLATE_SPHEROIDAL)
          CALL InitCommon(this%viscosity,ALPHA,"alpha")
       CASE DEFAULT
          CALL Error(this,"InitSources_viscosity", "Geometry not supported for alpha-viscosity")
       END SELECT
    CASE(BETA)
       SELECT CASE(GetType(Physics))
       CASE (EULER2D,EULER2D_ISOTHERM,EULER3D_ROTSYM)
          CALL InitCommon(this%viscosity,BETA,"beta")
       CASE DEFAULT
          CALL Error(this,"InitSources_viscosity", "Physics not supported for beta-viscosity")          
       END SELECT
    CASE(PRINGLE)
       CALL InitCommon(this%viscosity,PRINGLE,"const. kinematic viscosity")
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



  PURE SUBROUTINE UpdateViscosity(this,Mesh,Physics,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,pvar,cvar
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!

    SELECT CASE(GetType(this%viscosity))
    CASE(MOLECULAR)
       ! do nothing, since it has already been initialized in InitSources_viscosity
    CASE(ALPHA)
       ! Shakura-Sunyaev type alpha viscosity
       ! standard alpha prescription: nu = alpha*cs*h
       ! or nu = alpha / A * cs**2 / omega with A = -d ln(omega)/d ln(r)
       ! (see Kato, Fukue & Minishige: Black Hole Accretion Disks, 2008; equation (3.46))
       ! 
       ! this is a rough estimation assuming that the logarithmic derivative
       ! of the angular velocity (d ln(omega)/d ln(r)) is of the order of -1 
       SELECT CASE(GetType(Mesh%geometry))
       CASE(POLAR,LOGPOLAR,TANPOLAR,SINHPOLAR)
          SELECT CASE(GetType(Physics))
          CASE(EULER2D)
!CDIR OUTERUNROLL=8
             DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
                DO i=Mesh%IMIN,Mesh%IMAX
                   this%dynvis(i,j) = etafkt_alpha(this%dynconst,pvar(i,j,Physics%PRESSURE), &
                        pvar(i,j,Physics%YVELOCITY)/SQRT(Mesh%bccart(i,j,1)**2+Mesh%bccart(i,j,2)**2))
                END DO
             END DO
          CASE(EULER2D_ISOTHERM)
!CDIR OUTERUNROLL=8
             DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
                DO i=Mesh%IMIN,Mesh%IMAX
                   this%dynvis(i,j) = etafkt_alpha(this%dynconst,Physics%csiso**2*pvar(i,j,Physics%DENSITY), &
                        pvar(i,j,Physics%YVELOCITY)/SQRT(Mesh%bccart(i,j,1)**2+Mesh%bccart(i,j,2)**2))
                END DO
             END DO             
          END SELECT
       CASE (CYLINDRICAL,TANCYLINDRICAL)
!CDIR OUTERUNROLL=8
          DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
             DO i=Mesh%IMIN,Mesh%IMAX
                this%dynvis(i,j) = etafkt_alpha(this%dynconst,pvar(i,j,Physics%PRESSURE), &
                     pvar(i,j,Physics%ZVELOCITY)/Mesh%bccart(i,j,1))
             END DO
          END DO
       CASE(SPHERICAL,OBLATE_SPHEROIDAL,SINHSPHERICAL)
!CDIR OUTERUNROLL=8
          DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
             DO i=Mesh%IMIN,Mesh%IMAX
                this%dynvis(i,j) =  etafkt_alpha(this%dynconst,pvar(i,j,Physics%PRESSURE), &
                     pvar(i,j,Physics%ZVELOCITY)/Mesh%bccart(i,j,1))
             END DO
          END DO
       END SELECT
    CASE(BETA)
       ! Duschl type beta viscosity
       SELECT CASE(GetType(Physics))
       CASE (EULER2D,EULER2D_ISOTHERM)
          this%dynvis(:,:) = etafkt_beta(this%dynconst, &
               Mesh%bhy(:,:)*cvar(:,:,Physics%YMOMENTUM))
       CASE (EULER3D_ROTSYM)
          this%dynvis(:,:) = etafkt_beta(this%dynconst, &
               Mesh%bhz(:,:)*cvar(:,:,Physics%ZMOMENTUM))
       CASE (EULER3D_ROTAMT)
          this%dynvis(:,:) = etafkt_beta(this%dynconst, &
               cvar(:,:,Physics%ZMOMENTUM))
       END SELECT
    CASE(PRINGLE)
       ! constant kinematic viscosity
       this%dynvis(:,:) = etafkt_pringle(pvar(:,:,Physics%DENSITY))
    END SELECT

  CONTAINS
    ! some elemental functions for computation of the viskosity

    ! Alpha viscosity
    ELEMENTAL FUNCTION etafkt_alpha(alpha,P,omega) RESULT(eta)
      IMPLICIT NONE
      !----------------------------------------------------------------------!
      REAL, INTENT(IN) :: alpha,P,omega
      REAL :: eta
      !----------------------------------------------------------------------!
      eta = alpha * P / omega
    END FUNCTION etafkt_alpha
    
    ! Beta viscosity
    ELEMENTAL FUNCTION etafkt_beta(beta,L) RESULT(eta)
      IMPLICIT NONE
      !----------------------------------------------------------------------!
      REAL, INTENT(IN) :: beta,L
      REAL :: eta
      !----------------------------------------------------------------------!
      ! L = r * rho * v_phi is the angular momentum
      eta = beta * abs(L)
    END FUNCTION etafkt_beta
    
    ! Pringle disk (kinematic viscosity is constant)
    ELEMENTAL FUNCTION etafkt_pringle(rho) RESULT(eta)
      IMPLICIT NONE
      !----------------------------------------------------------------------!
      REAL, INTENT(IN) :: rho
      REAL :: eta
      !----------------------------------------------------------------------!
      eta = this%dynconst * rho   !here: dynconst is kinematic...
    END FUNCTION etafkt_pringle
    
  END SUBROUTINE UpdateViscosity

 
  PURE SUBROUTINE ExternalSources_viscosity(this,Mesh,Physics,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: cvar,pvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,cvar,pvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    CALL UpdateViscosity(this,Mesh,Physics,pvar,cvar)
    CALL CalculateStresses(Physics,Mesh,pvar,this%dynvis,this%bulkvis, &
            this%btxx,this%btxy,this%btxz,this%btyy,this%btyz,this%btzz)
    CALL ViscositySources(Physics,Mesh,pvar,this%btxx,this%btxy,this%btxz, &
         this%btyy,this%btyz,this%btzz,this%ftxx,this%ftxy,this%ftxz,this%ftyy, &
         this%ftyz,this%ftzz,sterm)
  END SUBROUTINE ExternalSources_viscosity

  
  PURE SUBROUTINE CalcTimestep_viscosity(this,Mesh,Physics,pvar,cvar,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar,cvar
    REAL              :: dt
    !------------------------------------------------------------------------!
    REAL              :: invdt_x, invdt_y
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,pvar,cvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: dt
    !------------------------------------------------------------------------!
    ! FIXME: this call is necessary, but its inefficient, since we call it
    !        again, when computing the viscous source term
    CALL UpdateViscosity(this,Mesh,Physics,pvar,cvar)

    ! compute kinematic viscosity
    this%kinvis(:,:) = this%dynvis(:,:) / pvar(:,:,Physics%DENSITY)

    ! x-direction
    IF (Mesh%INUM.GT.1) THEN
       invdt_x = MAXVAL(this%kinvis(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) &
            / Mesh%dlx(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)**2)
    ELSE
       ! set to zero, i.e. no limit in x-direction due to diffusion
       invdt_x = 0.0
    END IF
    ! y-direction
    IF (Mesh%JNUM.GT.1) THEN
       invdt_y = MAXVAL(this%kinvis(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) &
            / Mesh%dly(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)**2)
    ELSE
       ! set to zero, i.e. no limit in y-direction due to diffusion
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
    DEALLOCATE(this%dynvis,this%kinvis,this%bulkvis, &
               this%btxx,this%btyy,this%btzz,this%btxy,this%btxz,this%btyz,&
               this%ftxx,this%ftyy,this%ftzz,this%ftxy,this%ftxz,this%ftyz)
  END SUBROUTINE CloseSources_viscosity
 

END MODULE sources_viscosity
