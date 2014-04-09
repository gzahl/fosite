!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: gravity_common.f90                                                #
!#                                                                           #
!# Copyright (C) 2014                                                        #
!# Björn Sperling <sperling@astrophysik.uni-kiel.de>                         #
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
!> \defgroup gravity gravity
!! \{
!! \brief Family of gravity modules
!! \}
!----------------------------------------------------------------------------!
!> \author Björn Sperling
!!
!! \brief basic gravity module
!!
!! \extends common_types
!----------------------------------------------------------------------------!
MODULE gravity_common
  USE common_types, &
       GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Initialized_common => Initialized, Info_common => Info, &
       Warning_common => Warning, Error_common => Error
  USE mesh_common, ONLY : Selection_TYP
  USE common_dict, ONLY : Dict_TYP
  USE multipole_common, ONLY : Multipole_TYP
#ifdef HAVE_FFTW  
  USE fftw, ONLY : C_PTR, C_DOUBLE, C_DOUBLE_COMPLEX
#endif
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE GetType
     MODULE PROCEDURE GetGravityType, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetGravityTypeName, GetName_common
  END INTERFACE
  INTERFACE GetRank
     MODULE PROCEDURE GetGravityRank, GetRank_common
  END INTERFACE
  INTERFACE GetNumProcs
     MODULE PROCEDURE GetGravityNumProcs, GetNumProcs_common
  END INTERFACE
  INTERFACE Initialized
     MODULE PROCEDURE GravityInitialized, Initialized_common
  END INTERFACE
  INTERFACE Info
     MODULE PROCEDURE GravityInfo, Info_common
  END INTERFACE
  INTERFACE Warning
     MODULE PROCEDURE GravityWarning, Warning_common
  END INTERFACE
  INTERFACE Error
     MODULE PROCEDURE GravityError_rank0, GravityError_rankX, Error_common
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  !> data type for the multigrid poisson solver
  TYPE Grid_TYP
     !> \name Variables
     REAL, DIMENSION(:,:), POINTER   :: u,rho         !< DELTA u = rho
     !> some geometric values
     REAL, DIMENSION(:,:), POINTER   :: a,da,b,db,c,invc,d
     REAL, DIMENSION(:,:), POINTER   :: vol,bhx,bhy   !< volume, scale factors
     REAL, DIMENSION(:,:), POINTER   :: tmp           !< temp
     REAL, DIMENSION(:,:,:), POINTER :: bccart,curv   !< cart. and curv. coord.
     REAL, DIMENSION(:,:,:), POINTER :: tri           !< coefficients for block-iteration
     INTEGER                         :: ni,nj         !< resolution
     INTEGER, DIMENSION(:,:), POINTER :: ij2k, k2ij   !< index field for bound.
     REAL                            :: hi,hj         !< width of cells
     REAL                            :: invhi2,invhj2 !< square of inverse
  END TYPE Grid_TYP
  !> data type for all gravity modules
  TYPE Gravity_TYP
     !> \name Variables
     TYPE(Common_TYP)                :: gravitytype  !< type of gravity term
     TYPE(Gravity_TYP), POINTER      :: next => null() !< next gravity in list
     TYPE(Common_TYP)                :: potential    !< newton or wiita
     TYPE(Gravity_TYP), POINTER      :: pm => null() !< pointmass term
     REAL                            :: time         !< last update
     REAL                            :: mass         !< mass of pointmass
     REAL                            :: mass2        !< 2nd mass for binaries
     REAL                            :: excent       !< excentricity
     REAL                            :: semaaxis     !< semi major axis
     REAL                            :: period       !< period of binaries
     REAL, DIMENSION(2,2)            :: binpos       !< 2D cart. positions
     REAL                            :: mdot         !< disk accretion rate
     REAL                            :: eps1,eps2    !< softening parameter
     REAL                            :: r0(2)        !< location of pointmass
     !> time when the pointmass is fully switched on
     REAL                            :: switchon
     REAL, DIMENSION(:,:,:),POINTER  :: scaled       !< scaled acceleration field
     !> time of a orbital period at the inner and outer boundaries
     REAL, DIMENSION(2)              :: tau
     INTEGER                         :: outbound     !< outflow boundary
     REAL, DIMENSION(:,:,:), POINTER :: accel        !< acceleration
     REAL, DIMENSION(:,:), POINTER   :: radius,radius3!< distance to origin
     REAL, DIMENSION(:,:), POINTER   :: invr,invr_sec!< 1./radius
     !> x,y components of the vector between mesh points and primary/secondary
     !! star
     REAL, DIMENSION(:,:,:), POINTER :: r_prim,r_sec
     REAL, DIMENSION(:,:), POINTER   :: omega        !< angular velocity
     REAL, DIMENSION(:,:,:), POINTER :: omega2       !< Omega Kepler squared
     REAL, DIMENSION(:,:,:), POINTER :: gposvecr3    !< = GN*x/radius**3
     REAL, DIMENSION(:,:), POINTER   :: cellmass     !< rho*dV
     REAL, DIMENSION(:,:), POINTER   :: enclmass     !< enclosed mass
     REAL, DIMENSION(:,:), POINTER   :: height       !< disk height h
     REAL, DIMENSION(:,:), POINTER   :: tmp,tmp2,tmp3!< temp arrays
     REAL, DIMENSION(:,:), POINTER   :: cs           !< speed of sound pointer

     !> \name
     !!#### Poisson module
     TYPE(Common_TYP)                :: poissontype   !< type of source term
     TYPE(Multipole_TYP)             :: multipole     !< multipole expansion
     TYPE(Grid_TYP), POINTER         :: grid(:)       !< coarse grids
     REAL                            :: MAXRESIDNORM  !< max error of residuum
     !> type of relaxation method
     INTEGER                         :: RELAXTYPE
     INTEGER                         :: NGRID         !< number of grids
     INTEGER                         :: NPRE,NPOST    !< pre and post smoothing
     INTEGER                         :: NMAXCYCLE     !< max of iterations
     INTEGER                         :: MINRES        !< min resolution
     REAL, DIMENSION(:,:), POINTER   :: phi           !< potential
     REAL, DIMENSION(:,:,:,:), POINTER :: mphi        !< multiple potentials
     INTEGER                         :: n             !< number of potentials
     REAL,DIMENSION(:),POINTER       :: s0, sdelta    !< ramp fn: s0 + sdelta*t
     REAL,DIMENSION(:),POINTER       :: lastfac       !< last switchon factors
     INTEGER, DIMENSION(4)           :: Boundary      !< boundary condition
     LOGICAL                         :: DIRICHLET     !< true if min ONE bound.
                                                      !< boundary cond. is set
#ifdef HAVE_FFTW
    !> \name
    !!#### spectral poisson solver
    !> plan for real to complex fourier transforms
    TYPE(C_PTR)                      :: plan_r2c
    !> plan for complex to real fourier transforms
    TYPE(C_PTR)                      :: plan_c2r
    REAL(C_DOUBLE), POINTER          :: row(:)        !< temporary variable
    TYPE(C_PTR)                      :: p_row         !< temporary variable
    COMPLEX(C_DOUBLE_COMPLEX), POINTER &
                                     :: Frow(:)       !< temporary variable
    TYPE(C_PTR)                      :: p_Frow        !< temporary variable
    !> Important precalculated matrix - fourier transformed I
    REAL(C_DOUBLE), DIMENSION(:,:,:), POINTER &
                                     :: FI
    TYPE(C_PTR)                      :: p_FI
#endif
    INTEGER                          :: green
    REAL                             :: sigma
    !> local IMAX, INUM
    INTEGER                          :: IMAX, INUM
    !> \name Variables in Parallel Mode
#ifdef PARALLEL
    INTEGER                          :: error         !< MPI error flag
    REAL,DIMENSION(:),POINTER        :: sbuf1,sbuf2,rbuf1,rbuf2
    !> displacment and length of domain
    INTEGER,DIMENSION(:),POINTER     :: displ, num
#endif
!FIXME
    INTEGER, DIMENSION(:),POINTER :: relaxcount !test!!!
  END TYPE Gravity_TYP
  !> \}
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Gravity_TYP, &
       Grid_TYP, &
       ! methods
       InitGravity, &
       CloseGravity, &
       GetGravityPointer, &
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
  SUBROUTINE InitGravity(this,stype,sname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this
    INTEGER           :: stype
    CHARACTER(LEN=32) :: sname
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: newgrav, tmpgrav
    TYPE(Gravity_TYP) :: errgrav      ! we need this only for error reporting !
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: stype,sname
    !------------------------------------------------------------------------!
    ! allocate memory for new gravity term
    ALLOCATE(newgrav,STAT=err)
    IF (err.NE.0) CALL Error(errgrav,"InitGravity", "Unable allocate memory!")
    
     ! basic initialization
    CALL InitCommon(newgrav%gravitytype,stype,sname)

    ! add new gravity term to beginning of
    ! list of gravity terms
    IF (.NOT.ASSOCIATED(this)) THEN
       this => newgrav
       NULLIFY(this%next)
    ELSE
       tmpgrav => this
       this => newgrav
       this%next => tmpgrav
    END IF
  END SUBROUTINE InitGravity


  !> \public
  SUBROUTINE CloseGravity(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL CloseCommon(this%gravitytype)
  END SUBROUTINE CloseGravity


  !> \public
  FUNCTION GetGravityPointer(list,stype) RESULT(gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: list,gp
    INTEGER, INTENT(IN) :: stype
    !------------------------------------------------------------------------!
    gp => list
    DO
       IF (ASSOCIATED(gp).EQV..FALSE.) EXIT
!CDIR IEXPAND
       IF (GetType(gp).EQ.stype) RETURN
       gp => gp%next
    END DO
  END FUNCTION GetGravityPointer


  PURE FUNCTION GetGravityRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), INTENT(IN) :: this
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%gravitytype)
  END FUNCTION GetGravityRank


  PURE FUNCTION GetGravityNumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), INTENT(IN) :: this
    INTEGER :: p
    !------------------------------------------------------------------------!
    p = GetNumProcs_common(this%gravitytype)
  END FUNCTION GetGravityNumProcs


  PURE FUNCTION GetGravityType(this) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), INTENT(IN) :: this
    INTEGER :: st
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    st = GetType_common(this%gravitytype)
  END FUNCTION GetGravityType


  PURE FUNCTION GetGravityTypeName(this) RESULT(sn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: sn
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    sn = GetName_common(this%gravitytype)
  END FUNCTION GetGravityTypeName

  PURE FUNCTION GravityInitialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), INTENT(IN) :: this
    LOGICAL :: i
    !------------------------------------------------------------------------!
    i = Initialized_common(this%gravitytype)
  END FUNCTION GravityInitialized


  SUBROUTINE GravityInfo(this,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: msg
    !------------------------------------------------------------------------!
    CALL Info_common(this%gravitytype,msg)
  END SUBROUTINE GravityInfo


  SUBROUTINE GravityWarning(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Warning_common(this%gravitytype,modproc,msg)
  END SUBROUTINE GravityWarning


  SUBROUTINE GravityError_rank0(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Error_common(this%gravitytype,modproc,msg)
  END SUBROUTINE GravityError_rank0


  SUBROUTINE GravityError_rankX(this,modproc,msg,rank)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    INTEGER, INTENT(IN)           :: rank
    !------------------------------------------------------------------------!
    CALL Error_common(this%gravitytype,modproc,msg,rank)
  END SUBROUTINE GravityError_rankX

END MODULE gravity_common
