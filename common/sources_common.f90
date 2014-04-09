!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_common.f90                                                #
!#                                                                           #
!# Copyright (C) 2006-2014                                                   #
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
!> \defgroup sources sources
!! \{
!! \brief Family of sources modules
!! \}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!!
!! \brief basic sources module
!!
!! \extends common_types
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_common
  USE common_types, &
       GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Initialized_common => Initialized, Info_common => Info, &
       Warning_common => Warning, Error_common => Error
  USE mesh_common, ONLY : Selection_TYP
  USE common_dict, ONLY : Dict_TYP
  USE gravity_common, ONLY: Gravity_TYP
#ifdef HAVE_FFTW  
  USE fftw, ONLY : C_PTR, C_DOUBLE, C_DOUBLE_COMPLEX
#endif
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE GetType
     MODULE PROCEDURE GetSourceType, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetSourceTypeName, GetName_common
  END INTERFACE
  INTERFACE GetRank
     MODULE PROCEDURE GetSourcesRank, GetRank_common
  END INTERFACE
  INTERFACE GetNumProcs
     MODULE PROCEDURE GetSourcesNumProcs, GetNumProcs_common
  END INTERFACE
  INTERFACE Initialized
     MODULE PROCEDURE SourcesInitialized, Initialized_common
  END INTERFACE
  INTERFACE Info
     MODULE PROCEDURE SourcesInfo, Info_common
  END INTERFACE
  INTERFACE Warning
     MODULE PROCEDURE SourcesWarning, Warning_common
  END INTERFACE
  INTERFACE Error
     MODULE PROCEDURE SourcesError_rank0, SourcesError_rankX, Error_common
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  TYPE Sources_TYP
     !> \name Variables
     TYPE(Common_TYP)                :: sourcetype   !< type of source term
     TYPE(Sources_TYP), POINTER      :: next => null() !< next source in list
     TYPE(Gravity_TYP), POINTER      :: glist => null()!< gravity list
     TYPE(Common_TYP)                :: viscosity    !< molecular,alpha,beta
     REAL                            :: time         !< simulation time
     INTEGER                         :: timeid       !<    update of this id?
     !> 0: no src term in energy equation
     !! 1: src term in energy equation
     LOGICAL                         :: addtoenergy
!FIXME
     REAL                            :: mass         !< mass for diskthomson
     REAL                            :: mdot         !< disk accretion rate
     REAL                            :: dynconst,bulkconst ! viscosity const.
     REAL                            :: cvis         !< viscous Courant no.
     REAL                            :: eps1,eps2    !< softening parameter
     !> angular velocity of the rotating frame
     REAL                            :: rotomega
     !> effective surface of the dust grains
     REAL                            :: a_eff
     !> optical properties of the disk surface
     REAL                            :: kappa
     REAL                            :: T_star       !< temperature star
     REAL                            :: R_star       !< radius star
     REAL                            :: Qabs         !< dust absorption efficiency
     REAL                            :: T_sublim_min !< dust starts to sublimate
     !> dust density = 0.0 due to sublimation
     REAL                            :: T_sublim_max
     !> \name
     !!#### wave_damping

     !> inner and outer ware damping boundaries
     REAL, DIMENSION(2)              :: r
     !> time of a orbital period at the inner and outer boundaries
     REAL, DIMENSION(2)              :: tau
     !> \name
     !!####
     INTEGER                         :: dust_type    !< select mc3d dust catalogue
     !> 1: binary primary component or single star heating
     !! 2: secondary
     INTEGER                         :: star
     !< 1: calc boundary height,
     !! 0: boundaries are already correct
     INTEGER,DIMENSION(:),POINTER    :: calc_boundary
     REAL, DIMENSION(:,:,:), POINTER :: accel,accart !< acceleration
     REAL, DIMENSION(:,:,:), POINTER :: bcposvec,bccart !< position vector
     REAL, DIMENSION(:,:), POINTER   :: radius       !< distance to origin
     REAL, DIMENSION(:,:), POINTER   :: invr         !< 1./radius
     REAL, DIMENSION(:,:), POINTER   :: cs           !< speed of sound
     REAL, DIMENSION(:,:), POINTER   :: omega        !< angular velocity
     REAL, DIMENSION(:,:,:), POINTER :: omega2       !< Omega Kepler squared
     REAL, DIMENSION(:,:,:), POINTER :: gxr3         !< = GN*x/radius**3
     REAL, DIMENSION(:,:), POINTER   :: cellmass     !< rho*dV
     ! --------------------------------------------------------------------- !
     !> \name
     !!#### planet_heating/cooling
     REAL, DIMENSION(:,:), POINTER   :: Energy       !< energy
     REAL, DIMENSION(:,:), POINTER   :: RHO_s        !< surface density
     REAL, DIMENSION(:,:), POINTER   :: P_s          !< surface pressure
     REAL, DIMENSION(:,:), POINTER   :: T_s          !< temperatur
     REAL                            :: tau_inf      !< optical depth
     REAL                            :: T_0          !< equilibrium temp
     REAL, DIMENSION(:,:), POINTER   :: T_init       !< initial temperatur
     REAL                            :: intensity
     REAL                            :: albedo
     REAL                            :: distance     !< distance planet<->star
     REAL                            :: R_planet     !< radius planet
     REAL                            :: mu           !< molar mass
     REAL                            :: theta0       !< shifting theta-angle
     REAL                            :: phi0         !< shifting phi-angle
     REAL                            :: omegasun     !< day-night omega
     REAL                            :: year         !< trop. year of a planet
     REAL, DIMENSION(:,:,:), POINTER :: centproj     !< rot.frame centr.3d->2d
     REAL, DIMENSION(:,:,:), POINTER :: cos1
     REAL, DIMENSION(:,:,:), POINTER :: sin1
     REAL                            :: c_p          !< spec. heat cap.
     REAL                            :: gacc         !< grav. acceleration
     REAL                            :: gamma        !< rat. spec. heats
     ! --------------------------------------------------------------------- !
     REAL, DIMENSION(:,:), POINTER   :: height       !< disk height h
     REAL, DIMENSION(:,:), POINTER   :: Sigma_dust   !< dust surface desnity
     REAL, DIMENSION(:,:), POINTER   :: invheight2   !< 1/h**2
     REAL, DIMENSION(:,:), POINTER   :: Qcool        !< cooling source
     REAL, DIMENSION(:,:), POINTER   :: Qstar        !< stellar heating source
     REAL, DIMENSION(:,:), POINTER   :: H_tau1       !< z(tau=1)
     REAL, DIMENSION(:,:), POINTER   :: rescale      !< convert tau_z to tau_s
     !> angle between disk surface and line of sight to the star
     REAL, DIMENSION(:,:), POINTER   :: flaring_angle
     !> inverse 3d distance to heating central object
     REAL, DIMENSION(:,:), POINTER   :: invdr_3d
     !> vector from heating central object to mesh points, 3dim
     REAL, DIMENSION(:,:,:), POINTER :: Distance_3d
     REAL, DIMENSION(:,:,:), POINTER :: D_3d_norm    !< nornmalized Distance_3d
     !> normal vector of the disk surface in curvilinear coordinates
     REAL, DIMENSION(:,:,:), POINTER :: n
     REAL, DIMENSION(:,:,:), POINTER :: n_cart       !< cartesian n
     REAL, DIMENSION(:,:), POINTER   :: dynvis, &    !< dynamic viscosity
                                        kinvis, &    !< kinematic viscosity
                                        bulkvis      !< bulk viscosity
     !> components of the stress tensor
     REAL, DIMENSION(:,:), POINTER   :: btxx,btyy,&
          btzz,btxy,btxz,btyz,tmp,tmp2,tmp3
     REAL, DIMENSION(:,:,:), POINTER :: tmp4         !<    temp array
     REAL, DIMENSION(:,:), POINTER   :: Sxx,Syy,&
          Szz,Sxy,Sxz,Syz
     !> source terms of sgs module
     REAL, DIMENSION(:,:), POINTER   :: diff,rhoeps,&
                                        sigma
     REAL, DIMENSION(:,:,:), POINTER :: ftxx,ftyy,&
          ftzz,ftxy,ftxz,ftyz
     REAL, DIMENSION(:,:), POINTER   :: delta        !< half-width of the filter
     REAL, DIMENSION(:,:,:), POINTER :: cent         !< rot. frame centrifugal
     REAL, DIMENSION(:,:,:), POINTER :: init_pvar    !< pvar state from init
     REAL, DIMENSION(:,:,:), POINTER  :: ptr3        !< 3d pointer
#ifdef HAVE_FFTW
     REAL, DIMENSION(:,:,:), POINTER :: fk, rand     !< forcing (fourier)
     !> \name
     !!#### forcing
     REAL                            :: L,K0,F0,CHI,&
                                        T,invsqrtN,&
                                        stoptime
     TYPE(C_PTR)                     :: plan_r2r     !< fftw plan (real to real)
     REAL(C_DOUBLE), POINTER         :: temp_c(:,:)  !< fftw temp storage
     REAL(C_DOUBLE), POINTER         :: Ftemp_c(:,:)
#endif
  END TYPE Sources_TYP
  !> \}
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! methods
       InitSources, &
       CloseSources, &
       GetSourcesPointer, &
       GetType, &
       GetName, &
       GetRank, &
       GetNumProcs, &
       Initialized, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public
  SUBROUTINE InitSources(this,stype,sname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    INTEGER           :: stype
    CHARACTER(LEN=32) :: sname
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: newsrc, tmpsrc
    TYPE(Sources_TYP) :: errsrc      ! we need this only for error reporting !
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: stype,sname
    !------------------------------------------------------------------------!
    ! allocate memory for new source term
    ALLOCATE(newsrc,STAT=err)
    IF (err.NE.0) CALL Error(errsrc,"InitSources", "Unable allocate memory!")
    
     ! basic initialization
    CALL InitCommon(newsrc%sourcetype,stype,sname)

    ! add new source term to beginning of
    ! list of source terms
    IF (.NOT.ASSOCIATED(this)) THEN
       this => newsrc
       NULLIFY(this%next)
    ELSE
       tmpsrc => this
       this => newsrc
       this%next => tmpsrc
    END IF
  END SUBROUTINE InitSources


  !> \public
  SUBROUTINE CloseSources(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL CloseCommon(this%sourcetype)
  END SUBROUTINE CloseSources


  !> \public
  FUNCTION GetSourcesPointer(list,stype) RESULT(sp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: list,sp
    INTEGER, INTENT(IN) :: stype
    !------------------------------------------------------------------------!
    sp => list
    DO
       IF (ASSOCIATED(sp).EQV..FALSE.) EXIT
!CDIR IEXPAND
       IF (GetType(sp).EQ.stype) RETURN
       sp => sp%next
    END DO
  END FUNCTION GetSourcesPointer


  PURE FUNCTION GetSourcesRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%sourcetype)
  END FUNCTION GetSourcesRank


  PURE FUNCTION GetSourcesNumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    INTEGER :: p
    !------------------------------------------------------------------------!
    p = GetNumProcs_common(this%sourcetype)
  END FUNCTION GetSourcesNumProcs


  PURE FUNCTION GetSourceType(this) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    INTEGER :: st
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    st = GetType_common(this%sourcetype)
  END FUNCTION GetSourceType


  PURE FUNCTION GetSourceTypeName(this) RESULT(sn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: sn
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    sn = GetName_common(this%sourcetype)
  END FUNCTION GetSourceTypeName

  PURE FUNCTION SourcesInitialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    LOGICAL :: i
    !------------------------------------------------------------------------!
    i = Initialized_common(this%sourcetype)
  END FUNCTION SourcesInitialized


  SUBROUTINE SourcesInfo(this,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: msg
    !------------------------------------------------------------------------!
    CALL Info_common(this%sourcetype,msg)
  END SUBROUTINE SourcesInfo


  SUBROUTINE SourcesWarning(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Warning_common(this%sourcetype,modproc,msg)
  END SUBROUTINE SourcesWarning


  SUBROUTINE SourcesError_rank0(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Error_common(this%sourcetype,modproc,msg)
  END SUBROUTINE SourcesError_rank0


  SUBROUTINE SourcesError_rankX(this,modproc,msg,rank)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    INTEGER, INTENT(IN)           :: rank
    !------------------------------------------------------------------------!
    CALL Error_common(this%sourcetype,modproc,msg,rank)
  END SUBROUTINE SourcesError_rankX

END MODULE sources_common
