!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_common.f90                                                #
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
!> \defgroup physics physics
!! \{
!! \brief Family of physics modules
!! \}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!!
!! \brief basic physics module
!!
!! \extends common_types
!! \ingroup physics
!----------------------------------------------------------------------------!
MODULE physics_common
  USE common_types, &
       GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Initialized_common => Initialized, Info_common => Info, &
       Warning_common => Warning, Error_common => Error
  USE sources_common, ONLY : Sources_TYP
  USE constants_common, ONLY : Constants_TYP
  USE common_dict, ONLY : Dict_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE GetType
     MODULE PROCEDURE GetAdvProblem, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetAdvProblemName, GetName_common
  END INTERFACE
  INTERFACE GetRank
     MODULE PROCEDURE GetPhysicsRank, GetRank_common
  END INTERFACE
  INTERFACE GetNumProcs
     MODULE PROCEDURE GetPhysicsNumProcs, GetNumProcs_common
  END INTERFACE
  INTERFACE Initialized
     MODULE PROCEDURE PhysicsInitialized, Initialized_common
  END INTERFACE
  INTERFACE Info
     MODULE PROCEDURE PhysicsInfo, Info_common
  END INTERFACE
  INTERFACE Warning
     MODULE PROCEDURE PhysicsWarning, Warning_common
  END INTERFACE
  INTERFACE Error
     MODULE PROCEDURE PhysicsError, Error_common
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  TYPE Physics_TYP
     !> \name Variables
     TYPE(Common_TYP)    :: advproblem            !< advection problem
     TYPE(Constants_TYP) :: constants             !< physical constants
     TYPE(Sources_TYP), POINTER &
                         :: sources => null()     !< list of source terms
     REAL                :: gamma,&               !< ratio of spec. heats
                            time,&                !< simulation time       
                            mu, &                 !< mean molecular weight
                            csiso, &              !< isothermal sound speed
                            Omega, &              !< speed of rotating
                            centrot(2), &         !< center of rotation
                            eps                   !< softening length
     INTEGER             :: VNUM, &               !< number of variables
                            DIM, &                !< Dimension (2 or 3)
                            DENSITY,PRESSURE,ENERGY,SGSPRESSURE,SGSENERGY, &
                            XVELOCITY,XMOMENTUM,YVELOCITY,YMOMENTUM,&
                            ZVELOCITY,ZMOMENTUM   !< array indicies for primitive and conservative variables
     LOGICAL             :: supports_absorbing    !< absorbing boundary conditions supported
                            !! \details .TRUE. if absorbing boundary conditions are supported by the physics module
     LOGICAL             :: supports_farfield     !< farfield boundary conditions supported
                            !! \details .TRUE. if farfield boundary conditions are supported by the physics module
     CHARACTER(LEN=16), DIMENSION(:), POINTER &
                         :: pvarname,cvarname     !< names of variables
     REAL, DIMENSION(:,:), POINTER &
                         :: bccsound, &           !< bary centered speed of sound
                            amin, amax, &
                            bmin, bmax, &         !< wave speeds
                            bcradius, &           !< distance to the origin bary center values
                            divposvec, &          !< divergence of the position vector
                            bphi, &               !< bary centered constant gravitational potential
                            tmp                   !< temporary storage
     REAL, DIMENSION(:,:,:), POINTER &
                         :: fcsound, &            !< speed of sound faces
                            fradius, &            !< distance to the origin face values
                            tmin, tmax, &         !< temporary storage
                            bcposvec, &           !< curvilinear components of the position vector bary center values
                            w, &                  !< fargo bulk velocity
                            fphi, &               !< face centered constant gravitational potential
                            hy                    !< chy or fhy depending on reconstruction
     REAL, DIMENSION(:,:,:,:), POINTER &
                         :: fcent, &              !< centrifugal force
                            fposvec               !< curvilinear components of the position vector face values
  END TYPE Physics_TYP
  !> \}
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Physics_TYP, &
       ! methods
       InitPhysics, &
       ClosePhysics, &
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

  !> \public Constructor of basic physics module
  SUBROUTINE InitPhysics(this,atype,aname,vnum)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    INTEGER           :: atype
    CHARACTER(LEN=32) :: aname
    INTEGER           :: vnum
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: atype,aname,vnum
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%advproblem,atype,aname)
    this%vnum=vnum
    ALLOCATE(this%pvarname(this%vnum),this%cvarname(this%vnum),STAT=err)
    IF (err.NE.0) CALL Error(this,"InitPhysics", &
         "unable to allocate memory")
  END SUBROUTINE InitPhysics


  !> \public
  PURE FUNCTION GetAdvProblem(this) RESULT(ap)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP), INTENT(IN) :: this
    INTEGER :: ap
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    ap = GetType_common(this%advproblem)
  END FUNCTION GetAdvProblem


  !> \public
  PURE FUNCTION GetAdvProblemName(this) RESULT(an)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: an
    !------------------------------------------------------------------------!
    an = GetName_common(this%advproblem)
  END FUNCTION GetAdvProblemName


  !> \public
  PURE FUNCTION GetPhysicsRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP), INTENT(IN) :: this
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%advproblem)
  END FUNCTION GetPhysicsRank

  !> \public
  PURE FUNCTION GetPhysicsNumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP), INTENT(IN) :: this
    INTEGER :: p
    !------------------------------------------------------------------------!
    p = GetNumProcs_common(this%advproblem)
  END FUNCTION GetPhysicsNumProcs

  !> \public
  PURE FUNCTION PhysicsInitialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP), INTENT(IN) :: this
    LOGICAL :: i
    !------------------------------------------------------------------------!
    i = Initialized_common(this%advproblem)
  END FUNCTION PhysicsInitialized

  !> \public
  SUBROUTINE PhysicsInfo(this,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: msg
    !------------------------------------------------------------------------!
    CALL Info_common(this%advproblem,msg)
  END SUBROUTINE PhysicsInfo


  !> \public
  SUBROUTINE PhysicsWarning(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Warning_common(this%advproblem,modproc,msg)
  END SUBROUTINE PhysicsWarning


  !> \public
  SUBROUTINE PhysicsError(this,modproc,msg,rank)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    INTEGER, OPTIONAL, INTENT(IN) :: rank
    !------------------------------------------------------------------------!
    IF (PRESENT(rank)) THEN
       CALL Error_common(this%advproblem,modproc,msg,rank)
    ELSE
       CALL Error_common(this%advproblem,modproc,msg)
    END IF
  END SUBROUTINE PhysicsError


  !> \public Destructor of basic physics module
  SUBROUTINE ClosePhysics(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%pvarname,this%cvarname)
    CALL CloseCommon(this%advproblem)
  END SUBROUTINE ClosePhysics


END MODULE physics_common
