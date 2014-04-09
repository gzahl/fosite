!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: output_common.f90                                                 #
!#                                                                           #
!# Copyright (C) 2007 Tobias Illenseer <tillense@ita.uni-heidelberg.de>      #
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
! basic output module
!----------------------------------------------------------------------------!
MODULE output_common
  USE common_types, GetType_common => GetType, GetName_common => GetName
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE GetType
     MODULE PROCEDURE GetOutputFormat, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetOutputFormatName, GetName_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  TYPE Datastruct_TYP
     CHARACTER(LEN=32)          :: name            ! name of output field    !
     INTEGER                    :: rank            ! scalar, vector, etc.    !
     INTEGER                    :: shape           ! e.g. vector components  !
  END TYPE Datastruct_TYP
  TYPE Output_TYP
     TYPE(Common_TYP)       :: outputformat        ! output data format      !
     TYPE(Datastruct_TYP), DIMENSION(:), POINTER  &
                            :: datastruct       ! structure of output fields !
     INTEGER                :: outunit          ! output file unit           !
     CHARACTER(LEN=256)     :: filename         ! output file name           !
     CHARACTER(LEN=64)      :: formatstring     ! for ASCII data output      !
     LOGICAL                :: openmode         ! append data to ex. file?   !
     REAL                   :: time             ! next output time           !
     INTEGER                :: step             ! time step counter          !
     REAL                   :: start,end        ! output time intervall      !
     INTEGER                :: count            ! number of output datasets  !
  END TYPE Output_TYP
  !--------------------------------------------------------------------------!
  LOGICAL, PARAMETER :: APPEND    = .TRUE.
  LOGICAL, PARAMETER :: OVERWRITE = .FALSE.
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Datastruct_TYP, &
       Output_TYP, &
       ! constants
       APPEND, OVERWRITE, &
       ! methods
       InitOutput, &
       GetType, &
       GetName, &
       GetFilename, &
       GetUnit, &
       SetFormatString, &
       GetFormatString, &
       GetMode, &
       GetStart, &
       GetEnd, &
       IncTime
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitOutput(this,oformat,formatname,ounit,ofname,omode,start,end, &
       count)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP)   :: this
    INTEGER            :: oformat,ounit
    CHARACTER(LEN=32)  :: formatname
    CHARACTER(LEN=*)   :: ofname
    LOGICAL            :: omode
    REAL               :: start,end
    INTEGER            :: count
    !------------------------------------------------------------------------!
    INTEGER            :: ios
    LOGICAL            :: ex
    !------------------------------------------------------------------------!
    INTENT(IN)         :: oformat,formatname,ounit,ofname,omode,start,end,count
    INTENT(OUT)        :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%outputformat,oformat,formatname)
    this%outunit  = ounit
    this%filename = ofname
    this%openmode = omode
    this%start    = start
    this%end      = end
    this%count    = count
    this%time     = this%start
    this%step     = 0

    ! check if file already exists
    INQUIRE(FILE=this%filename, EXIST=ex)

    IF (.NOT. ex) THEN
       ! test if we can create a new file
       OPEN(this%outunit, IOSTAT=ios, FILE=this%filename, STATUS="NEW", ACTION="WRITE")
       CLOSE(this%outunit)
       IF (ios > 0) THEN
          PRINT *, "ERROR in CheckFile: Can't create new file ", &
               ACHAR(34), TRIM(this%filename), ACHAR(34)
          STOP
       END IF
    ELSE
       ! if file exists and we are not in APPEND mode
       ! check if we can overwrite the file
       IF (this%openmode.NEQV.APPEND) THEN
          OPEN(this%outunit, IOSTAT=ios, FILE=this%filename, STATUS="OLD", &
               POSITION="REWIND", ACTION="WRITE")
          CLOSE(this%outunit)
          IF (ios > 0) THEN
             PRINT *, "ERROR in CheckFile: Can't overwrite old file ", &
                  ACHAR(34), TRIM(this%filename), ACHAR(34)
             STOP
          END IF
       END IF
    END IF
  END SUBROUTINE InitOutput


  PURE FUNCTION GetOutputFormat(this) RESULT(of)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP), INTENT(IN) :: this
    INTEGER :: of
    !------------------------------------------------------------------------!
    of = GetType_common(this%outputformat)
  END FUNCTION GetOutputFormat


  PURE FUNCTION GetOutputFormatName(this) RESULT(fn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: fn
    !------------------------------------------------------------------------!
    fn = GetName_common(this%outputformat)
  END FUNCTION GetOutputFormatName


  PURE FUNCTION GetFilename(this) RESULT(fn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP), INTENT(IN) :: this
    CHARACTER(LEN=256) :: fn
    !------------------------------------------------------------------------!
    fn = this%filename
  END FUNCTION GetFilename


  PURE FUNCTION GetUnit(this) RESULT(ou)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP), INTENT(IN) :: this
    INTEGER :: ou
    !------------------------------------------------------------------------!
    ou = this%outunit
  END FUNCTION GetUnit


  PURE SUBROUTINE SetFormatString(this,fs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP), INTENT(OUT)  :: this
    CHARACTER(LEN=64), INTENT(IN)  :: fs
    !------------------------------------------------------------------------!
    this%formatstring = fs
  END SUBROUTINE SetFormatString


  PURE FUNCTION GetFormatString(this) RESULT(fs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP), INTENT(IN) :: this
    CHARACTER(LEN=64) :: fs
    !------------------------------------------------------------------------!
    fs = this%formatstring
  END FUNCTION GetFormatString


  PURE FUNCTION GetMode(this) RESULT(om)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP), INTENT(IN) :: this
    LOGICAL :: om
    !------------------------------------------------------------------------!
    om = this%openmode
  END FUNCTION GetMode


  PURE FUNCTION GetStart(this) RESULT(start)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP), INTENT(IN) :: this
    REAL :: start
    !------------------------------------------------------------------------!
    start = this%start
  END FUNCTION GetStart


  PURE FUNCTION GetEnd(this) RESULT(end)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP), INTENT(IN) :: this
    REAL :: end
    !------------------------------------------------------------------------!
    end = this%end
  END FUNCTION GetEnd

  
  PURE SUBROUTINE IncTime(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    this%time = this%time + ABS(this%end-this%start) / this%count
    this%step = this%step + 1
  END SUBROUTINE IncTime

END MODULE output_common
