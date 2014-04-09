!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: logio_common.f90                                                  #
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
! basic module for I/O of log files
!----------------------------------------------------------------------------!
MODULE logio_common
  USE common_types, GetType_common => GetType, GetName_common => GetName
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE GetType
     MODULE PROCEDURE GetLogFormat, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetLogFormatName, GetName_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  TYPE Logio_TYP
     PRIVATE
     TYPE(Common_TYP)       :: logformat           ! logfile format          !
     CHARACTER(LEN=256)     :: filename            ! logfile name            !
     INTEGER                :: logunit             ! logfile i/o unit        !
     INTEGER                :: logdt               ! timestep for log output !
  END TYPE Logio_TYP
  SAVE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Logio_TYP, &
       ! methods
       InitLogio, &
       GetType, &
       GetName, &
       GetFilename, &
       GetUnit, &
       GetLogstep
  !--------------------------------------------------------------------------!

CONTAINS

  PURE SUBROUTINE InitLogio(this,logfmt,fmtname,logunit,filename,logdt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Logio_TYP)    :: this
    INTEGER            :: logfmt,logunit
    CHARACTER(LEN=32)  :: fmtname
    CHARACTER(LEN=*)   :: filename
    INTEGER            :: logdt
    !------------------------------------------------------------------------!
    INTENT(IN)         :: logfmt,fmtname,logunit,filename,logdt
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%logformat,logfmt,fmtname)
    this%logunit  = logunit
    this%logdt    = logdt
    this%filename = filename
  END SUBROUTINE InitLogio


  PURE FUNCTION GetLogFormat(this) RESULT(lf)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Logio_TYP), INTENT(IN) :: this
    INTEGER :: lf
    !------------------------------------------------------------------------!
    lf = GetType_common(this%logformat)
  END FUNCTION GetLogFormat


  PURE FUNCTION GetLogFormatName(this) RESULT(fn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Logio_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: fn
    !------------------------------------------------------------------------!
    fn = GetName_common(this%logformat)
  END FUNCTION GetLogFormatName


  PURE FUNCTION GetFilename(this) RESULT(fn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Logio_TYP), INTENT(IN) :: this
    CHARACTER(LEN=256) :: fn
    !------------------------------------------------------------------------!
    fn = this%filename
  END FUNCTION GetFilename


  PURE FUNCTION GetUnit(this) RESULT(lu)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Logio_TYP), INTENT(IN) :: this
    INTEGER :: lu
    !------------------------------------------------------------------------!
    lu = this%logunit
  END FUNCTION GetUnit


  PURE FUNCTION GetLogstep(this) RESULT(ls)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Logio_TYP), INTENT(IN) :: this
    INTEGER :: ls
    !------------------------------------------------------------------------!
    ls = this%logdt
  END FUNCTION GetLogstep

END MODULE logio_common
