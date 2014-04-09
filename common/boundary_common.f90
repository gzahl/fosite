!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_common.f90                                               #
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
!> \defgroup boundary boundary
!! \{
!! \brief Family of boundary modules
!!
!! This is the family of boundary modules.The generic interface routines are
!! defined in the module \link boundary_generic \endlink. The basic boundary
!! data type and common basic subroutines and functions are defined in
!! \link boundary_common \endlink. Any other module of the family specifies
!! distinct boundary conditions.
!! \}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!!
!! \brief Basic module for boundary conditions.
!!
!! This module defines the basic data type, functions and subroutines
!! for all boundary modules.
!!
!! \extends common_types
!! \ingroup boundary
!----------------------------------------------------------------------------!
MODULE boundary_common
  USE common_types, &
       GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Initialized_common => Initialized, Info_common => Info, &
       Warning_common => Warning, Error_common => Error
  IMPLICIT NONE
#ifdef PARALLEL
#ifdef HAVE_MPIF_H
  include 'mpif.h'
#endif
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE GetType
     MODULE PROCEDURE GetCondition, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetConditionName, GetName_common
  END INTERFACE
  INTERFACE GetRank
     MODULE PROCEDURE GetBoundaryRank, GetRank_common
  END INTERFACE
  INTERFACE GetNumProcs
     MODULE PROCEDURE GetBoundaryNumProcs, GetNumProcs_common
  END INTERFACE
  INTERFACE Initialized
     MODULE PROCEDURE BoundaryInitialized, Initialized_common
  END INTERFACE
  INTERFACE Info
     MODULE PROCEDURE BoundaryInfo, Info_common
  END INTERFACE
  INTERFACE Warning
     MODULE PROCEDURE BoundaryWarning, Warning_common
  END INTERFACE
  INTERFACE Error
     MODULE PROCEDURE BoundaryError, Error_common
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  !> boundary data structure
  !!
  !! This data type stores all information on the boundary condition.
  TYPE Boundary_TYP
     !> \name Variables

     !> boundary type: no_gradients, reflecting, periodic, etc.
     TYPE(Common_TYP)  :: condition
     !> boundary orientation: west, east, south, north
     TYPE(Common_TYP)  :: direction
     INTEGER           :: IMID, &         !< i index of cell in the middle
                          JMID, &         !< j index of cell in the middle
                          nohdim          !< dimension of Noh problem
     LOGICAL           :: first_call      !< used in far-field bc
     LOGICAL,DIMENSION(:),POINTER :: &
                          reflX, &        !< mask array for reflecting bc
                          reflY           !< mask array for reflecting bc
     LOGICAL,DIMENSION(:,:),POINTER :: &
                          fixed           !< mask array for fixed bc
     INTEGER,DIMENSION(:,:),POINTER :: &
                          cbtype          !< custom boundary condition type
     REAL,DIMENSION(:,:),POINTER :: &
                          invr, &         !< inverse distance to center
                          Rscale, &       !< radial scaling constants
                          invRscale, &    !< inverse radial scaling constants
                          xvar, &         !< characteristic variables for absorbing bc
                          lambda, &       !< eigenvalues for absorbing bc
                          Rinv, &         !< Riemann invariants at the boundary
                          RinvInf         !< far field Riemann invariants
     !> boundary data array, e.g. for fixed,
     !! custom and noslip and boundary conditions
     REAL,DIMENSION(:,:,:),POINTER :: &
                          data
     !> \todo Probably obsolte variable. Remove it!
     REAL,DIMENSION(:,:,:),POINTER :: &
                          accel => NULL() !< pointer to gravitational accel
#ifdef PARALLEL
     !> \name Variables in Parallel Mode
     REAL,DIMENSION(:,:,:),POINTER :: &
                          sendbuf, &      !< send buffer for boundary data
                          recvbuf         !< receive buffer for boundary data
#endif
  END TYPE Boundary_TYP
  !--------------------------------------------------------------------------!
  !> \name Public Attributes
  CHARACTER(LEN=32), DIMENSION(4), PARAMETER :: & 
        direction_name = (/' west', ' east', 'south', 'north' /)
  !< string literal for each orientation
  INTEGER, PARAMETER :: &
     WEST  = 1, & !< named constant for western boundary
     EAST  = 2, & !< named constant for eastern boundary
     SOUTH = 3, & !< named constant for southern boundary
     NORTH = 4    !< named constant for northern boundary
  !> \}
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Boundary_TYP, &
       ! constants
#ifdef PARALLEL
       DEFAULT_MPI_REAL, &
#endif
       WEST, EAST, SOUTH, NORTH, &
       ! methods
       InitBoundary, &
       CloseBoundary, &
       GetType, &
       GetName, &
       GetDirection, &
       GetDirectionName, &
       GetRank, &
       GetNumProcs, &
       Initialized, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor of common boundary class.
  !!
  !! Initializes the boundary condition and direction.
  SUBROUTINE InitBoundary(this,bctype,bcname,direction)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this       !< \param [in,out] this boundary type
    INTEGER            :: bctype     !< \param [in] bctype boundary condition number
    INTEGER            :: direction  !< \param [in] direction west, east, south or north
    CHARACTER(LEN=*)   :: bcname     !< \param [in] bcname boundary condition name
    !------------------------------------------------------------------------!
    INTENT(IN)    :: bctype,bcname,direction
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    ! set boundary condition
    CALL InitCommon(this%condition,bctype,bcname)
    ! check for wrong direction
    SELECT CASE(direction)
    CASE(WEST,EAST,NORTH,SOUTH)
       ! ok
    CASE DEFAULT
       CALL Error(this,"InitBoundary_common", "Unknown direction")
    END SELECT
    ! set direction
    CALL InitCommon(this%direction,direction,direction_name(direction))
  END SUBROUTINE InitBoundary


  !> \public Destructor of common boundary class.
  SUBROUTINE CloseBoundary(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(INOUT) :: this !< \param [in] this boundary type
    !------------------------------------------------------------------------!
    CALL CloseCommon(this%condition)
  END SUBROUTINE CloseBoundary


  !> \public Get the boundary condition number; 
  !! overloads \b GetType from \link common_types::GetType \endlink
  !! \return boundary condition number
  PURE FUNCTION GetCondition(this) RESULT(bt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this !< \param [in] this boundary type
    INTEGER :: bt
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    bt = GetType_common(this%condition)
  END FUNCTION GetCondition


  !> \public Get the boundary condition name; 
  !! overloads \b GetName from \link common_types::GetName \endlink
  !! \return boundary condition name
  PURE FUNCTION GetConditionName(this) RESULT(bn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this !< \param [in] this boundary type
    CHARACTER(LEN=32) :: bn
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    bn = GetName_common(this%condition)
  END FUNCTION GetConditionName


  !> \public Get the direction number.
  !! \return direction number
  PURE FUNCTION GetDirection(this) RESULT(dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this !< \param [in] this boundary type
    INTEGER :: dir
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    dir = GetType_common(this%direction)
  END FUNCTION GetDirection


  !> \public Get the direction name.
  !! \return direction name
  PURE FUNCTION GetDirectionName(this) RESULT(dn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this !< \param [in] this boundary type
    CHARACTER(LEN=32) :: dn
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    dn = GetName_common(this%direction)
  END FUNCTION GetDirectionName


  !> \public Query initialization status;
  !! overloads \b Initialized from \link common_types::Initialized \endlink
  PURE FUNCTION BoundaryInitialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this !< \param [in] this boundary type
    LOGICAL :: i
    !------------------------------------------------------------------------!
    i = Initialized_common(this%condition)
  END FUNCTION BoundaryInitialized

 
  !> \public Get the MPI rank;
  !! overloads \b GetRank from \link common_types::GetRank \endlink
  !! \return MPI rank
  PURE FUNCTION GetBoundaryRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this !< \param [in] this boundary type
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%condition)
  END FUNCTION GetBoundaryRank


  !> \public Get the total number of MPI processes;
  !! overloads \b GetNumProcs from \link common_types::GetNumProcs \endlink
  !! \return number of MPI processes
  PURE FUNCTION GetBoundaryNumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this !< \param [in] this boundary type
    INTEGER :: p
    !------------------------------------------------------------------------!
    p = GetNumProcs_common(this%condition)
  END FUNCTION GetBoundaryNumProcs


  !> \public Print information on standard output;
  !! overloads \b Info from \link common_types::Info \endlink
  SUBROUTINE BoundaryInfo(this,msg,rank,node_info,tostderr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this !< \param [in] this boundary type
    CHARACTER(LEN=*),  INTENT(IN) :: msg  !< \param [in] msg info message
    INTEGER, OPTIONAL, INTENT(IN) :: rank !< \param [in] rank MPI rank
    LOGICAL, OPTIONAL, INTENT(IN) :: node_info !< \param [in] node_info enable rank output
    LOGICAL, OPTIONAL, INTENT(IN) :: tostderr  !< \param [in] tostderr enable STDERR output
    !------------------------------------------------------------------------!
    CALL Info_common(this%condition,msg,rank,node_info,tostderr)
  END SUBROUTINE BoundaryInfo


  !> \public Print warning message on standard error;
  !! overloads \b Warning from \link common_types::Warning \endlink
  SUBROUTINE BoundaryWarning(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this !< \param [in] this boundary type
    CHARACTER(LEN=*),  INTENT(IN) :: modproc !< \param [in] modproc name of module procedure
    CHARACTER(LEN=*),  INTENT(IN) :: msg !< \param [in] msg warning message
    !------------------------------------------------------------------------!
    CALL Warning_common(this%condition,modproc,msg)
  END SUBROUTINE BoundaryWarning


  !> \public Print error message on standard error and terminate the program;
  !! overloads \b Error from \link common_types::Error \endlink
  SUBROUTINE BoundaryError(this,modproc,msg,rank,node_info)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), INTENT(IN) :: this !< \param [in] this boundary type
    CHARACTER(LEN=*),  INTENT(IN) :: modproc !< \param [in] modproc name of module procedure
    CHARACTER(LEN=*),  INTENT(IN) :: msg !< \param [in] msg warning message
    INTEGER, OPTIONAL, INTENT(IN) :: rank !< \param [in] rank MPI rank
    LOGICAL, OPTIONAL, INTENT(IN) :: node_info !< \param [in] node_info enable rank output
    !------------------------------------------------------------------------!
    CALL Error_common(this%condition,modproc,msg,rank,node_info)
  END SUBROUTINE BoundaryError

END MODULE boundary_common
