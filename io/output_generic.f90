!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: output_generic.f90                                                #
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
! Subroutines for data output
!----------------------------------------------------------------------------!
MODULE output_generic
  USE output_gnuplot
  USE output_opendx
  USE physics_generic, ONLY : Physics_TYP, GetDataStruct
  USE mesh_common, ONLY : Mesh_TYP
  USE timedisc_common, ONLY : Timedisc_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: GNUPLOT = 1
  INTEGER, PARAMETER :: OPENDX  = 2
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Output_TYP, &
       ! constants
       GNUPLOT, OPENDX, &
       APPEND, OVERWRITE, &
       ! methods 
       InitOutput, &
       WriteOutput, &
       GetFilename, &
       GetStart, &
       GetEnd, &
       IncTime, &
       CloseOutput
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitOutput(this,Mesh,Physics,Timedisc, &
       filetype,filename,mode,starttime,stoptime,count)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP)   :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Timedisc_TYP) :: Timedisc
    INTEGER            :: filetype
    CHARACTER(LEN=*)   :: filename
    LOGICAL            :: mode
    REAL               :: starttime,stoptime
    INTEGER            :: count
    !------------------------------------------------------------------------!
    LOGICAL            :: ex
    INTEGER            :: ios
    INTEGER            :: err
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,Physics,Timedisc,filetype,filename,mode,starttime,stoptime,count
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!

    ! get info about the output data structure
    CALL GetDataStruct(Physics,this%datastruct)

    SELECT CASE(filetype)
    CASE(GNUPLOT)
       CALL InitOutput_gnuplot(this,Mesh,Physics,filetype,filename,mode, &
            starttime,stoptime,count)
    CASE(OPENDX)
       CALL InitOutput_opendx(this,Mesh,Physics,filetype,filename,mode, &
            starttime,stoptime,count)
       PRINT "(A)", "WARNING: OpenDX output not fully supported"
    CASE DEFAULT
       PRINT *,"ERROR in InitOutput: Unknown output file format"
       STOP
    END SELECT

    ! for log file input data
    DO WHILE (Timedisc%time.GT.this%time)
       CALL IncTime(this)
    END DO

    ! write initial data
    IF (this%step.EQ.0) THEN
       CALL WriteOutput(this,Mesh,Physics,Timedisc%pvar)
    ELSE
       CALL IncTime(this)
    END IF

    ! print some information
    PRINT "(A,A)", " OUTPUT---> data file type:    ", TRIM(GetName(this))
    PRINT "(A,A)", "            file name:         ", TRIM(GetFilename(this))
  END SUBROUTINE InitOutput


  SUBROUTINE WriteOutput(this,Mesh,Physics,ovar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP)   :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) :: ovar
    !------------------------------------------------------------------------!
    INTEGER            :: ios
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,Physics,ovar
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!

    ! open output file for writing
    OPEN(GetUnit(this), IOSTAT=ios, FILE=GetFilename(this), STATUS="OLD", &
         POSITION="APPEND", ACTION="WRITE")

    ! write output for specific file format
    IF (ios.LE.0) THEN
       SELECT CASE(GetType(this))
       CASE(GNUPLOT)
          CALL WriteOutput_gnuplot(this,Mesh,Physics,ovar)
       CASE(OPENDX)
          CALL WriteOutput_opendx(this,Mesh,Physics,ovar)
       CASE DEFAULT
          PRINT *,"ERROR in WriteOutput: Unknown output file format"
          CLOSE(GetUnit(this))
          STOP
       END SELECT
    ELSE
       PRINT *, "ERROR in WriteOutput:  Unable to write data to file ", &
            ACHAR(34), TRIM(GetFilename(this)), ACHAR(34) 
    END IF

    ! close output file
    CLOSE(GetUnit(this))

    ! increment output time and time step
    CALL IncTime(this)
  END SUBROUTINE WriteOutput


  SUBROUTINE CloseOutput(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP)   :: this
    !------------------------------------------------------------------------!
    INTEGER            :: ios
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this
    !------------------------------------------------------------------------!

    OPEN(GetUnit(this), IOSTAT=ios, FILE=GetFilename(this), STATUS="OLD", &
         POSITION="APPEND", ACTION="WRITE")

    IF (ios.LE.0) THEN
       SELECT CASE(GetType(this))
       CASE(1)
          ! gnuplot style plain ASCII table -> do nothing
       CASE(2)
          ! Open Data Explorer file format (dx)
          CALL CloseOutput_opendx(this)
       CASE DEFAULT
          PRINT *, "ERROR in CloseData: unknown file format"
          CLOSE(GetUnit(this))
          STOP
       END SELECT
    END IF

    CLOSE(GetUnit(this))
    
    ! deallocate memory for datastruct
    DEALLOCATE(this%datastruct)
  END SUBROUTINE CloseOutput

END MODULE output_generic
