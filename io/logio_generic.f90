!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: logio_generic.f90                                                 #
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
! Subroutines for logio of simulation data
!----------------------------------------------------------------------------!
MODULE logio_generic
  USE logio_binary, InitLogio_common => InitLogio
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  USE timedisc_common, ONLY : Timedisc_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: logformat_name = "no logio"
  INTEGER, PARAMETER :: logunit  = 13            ! log file unit number      !
  ! log file format flags
  INTEGER, PARAMETER :: NOLOG    = 1
  INTEGER, PARAMETER :: BINARY   = 2
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Logio_TYP, &
       ! constants
       NOLOG, BINARY, &
       ! methods 
       InitLogio, &
       ReadLogdata, &
       WriteLogdata, &
       GetLogstep, &
       CloseLogio
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitLogio(this,Mesh,Physics,Timedisc,logformat,filename,logdt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Logio_TYP)    :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Timedisc_TYP) :: Timedisc
    INTEGER            :: logformat
    CHARACTER(LEN=*)   :: filename
    INTEGER            :: logdt
    !------------------------------------------------------------------------!
    LOGICAL            :: ok
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,Physics,logformat,filename,logdt
    INTENT(INOUT)      :: this,Timedisc
    !------------------------------------------------------------------------!

    SELECT CASE(logformat)
    CASE(NOLOG)
       CALL InitLogio_common(this,logformat,logformat_name,logunit,filename,logdt)
       RETURN
    CASE(BINARY)
       CALL InitLogio_binary(this,Mesh,Physics,logformat,logunit,filename,logdt)
    CASE DEFAULT
       PRINT *,"ERROR in InitLogio: Unknown log file format"
       STOP
    END SELECT

    ! check if log file already exists
    INQUIRE(FILE=GetFilename(this), EXIST=ok)
    
    IF (ok) THEN
       CALL ReadLogdata(this,Mesh,Physics,Timedisc,ok)
       IF (.NOT.ok) THEN
          PRINT *,"ERROR in InitLogio: Unable to read log file data"
          STOP
       END IF
    ELSE
       CALL WriteHeader(this,Mesh,Physics)
       CALL WriteLogdata(this,Mesh,Physics,Timedisc)
    END IF

    ! print some information
    PRINT "(A,A)", " LOGIO----> data format:       ", TRIM(GetName(this))
    PRINT "(A,A)", "            file name:         ", TRIM(GetFilename(this))
  END SUBROUTINE InitLogio


  SUBROUTINE ReadLogdata(this,Mesh,Physics,Timedisc,ok)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Logio_TYP)    :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Timedisc_TYP) :: Timedisc
    LOGICAL            :: ok
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this,Mesh,Physics
    INTENT(OUT)        :: Timedisc
    !------------------------------------------------------------------------!
    
    SELECT CASE(GetType(this))
    CASE(NOLOG)
       ! do nothing
    CASE(BINARY)
       CALL ReadLogdata_binary(this,Mesh,Physics,Timedisc,ok)
    CASE DEFAULT
       PRINT *,"ERROR in ReadLogdata: Unknown log file format"
       STOP
    END SELECT
    
  END SUBROUTINE ReadLogdata


  SUBROUTINE WriteHeader(this,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Logio_TYP)    :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this,Mesh,Physics
    !------------------------------------------------------------------------!

    SELECT CASE(GetType(this))
    CASE(NOLOG)
       ! do nothing
    CASE(BINARY)
       CALL WriteHeader_binary(this,Mesh,Physics)
    CASE DEFAULT
       PRINT *,"ERROR in WriteHeader: Unknown log file format"
       STOP
    END SELECT
  END SUBROUTINE WriteHeader


  SUBROUTINE WriteLogdata(this,Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Logio_TYP)    :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Timedisc_TYP) :: Timedisc
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this,Mesh,Physics,Timedisc
    !------------------------------------------------------------------------!
    
    SELECT CASE(GetType(this))
    CASE(NOLOG)
       ! do nothing
    CASE(BINARY)
       CALL WriteLogdata_binary(this,Mesh,Physics,Timedisc)
    CASE DEFAULT
       PRINT *,"ERROR in ReadLogdata: Unknown log file format"
       STOP
    END SELECT
  END SUBROUTINE WriteLogdata


  SUBROUTINE CloseLogio
    IMPLICIT NONE
    !------------------------------------------------------------------------!

  END SUBROUTINE CloseLogio

END MODULE logio_generic
