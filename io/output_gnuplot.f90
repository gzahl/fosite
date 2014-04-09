!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: output_gnuplot.f90                                                #
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
! Subroutines for gnuplot tabular data output
!----------------------------------------------------------------------------!
MODULE output_gnuplot
  USE output_common
  USE physics_generic
  USE mesh_generic
  USE boundary_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: OUNIT = 12                    ! output unit number   !
  CHARACTER(LEN=32), PARAMETER :: format_name = "gnuplot"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Output_TYP, &
       ! constants
       APPEND, OVERWRITE, &
       ! methods 
       InitOutput_gnuplot, &
       WriteOutput_gnuplot, &
       GetType, &
       GetName, &
       GetFilename, &
       GetStart, &
       GetEnd, &
       GetUnit, &
       IncTime
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitOutput_gnuplot(this,Mesh,Physics,oformat,fname,omode,start,end, &
       count)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP)   :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: oformat
    CHARACTER(LEN=*)   :: fname
    CHARACTER(LEN=64)  :: fmt
    LOGICAL            :: omode
    REAL               :: start,end
    INTEGER            :: count
    !------------------------------------------------------------------------!
    INTEGER            :: ios
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,Physics,oformat,fname,omode,start,end,count
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    CALL InitOutput(this,oformat,format_name,OUNIT,fname,omode,start,end,count)

    ! output format string
    fmt = "(ES14.6)"
    CALL SetFormatString(this,fmt)

    ! write info to file header
    IF (GetMode(this).NEQV.APPEND) THEN
       OPEN(GetUnit(this), IOSTAT=ios, FILE=GetFilename(this), STATUS="OLD", &
            POSITION="REWIND", ACTION="WRITE")
       IF (ios.GT.0) THEN
          PRINT *, "ERROR in InitOutput_gnuplot: Unable to open output file "
          CLOSE(GetUnit(this))
          STOP
       ELSE
          WRITE (GetUnit(this),FMT="(A)",IOSTAT=ios) "# Data output of fosite"
          CLOSE(GetUnit(this))
          IF (ios.GT.0) THEN
             PRINT *, "ERROR in InitOutput_gnuplot: Unable to write file header "
             STOP
          END IF
       END IF
    END IF
  END SUBROUTINE InitOutput_gnuplot

  
  SUBROUTINE WriteOutput_gnuplot(this,Mesh,Physics,ovar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP)   :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) :: ovar
    !------------------------------------------------------------------------!
    INTEGER     :: i,j,k
    INTEGER     :: ios
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,Physics,ovar
    !------------------------------------------------------------------------!

    WRITE (GetUnit(this), "(A)", IOSTAT=ios) ACHAR(10)

    DO i=Mesh%IMIN,Mesh%IMAX
       DO j=Mesh%JMIN,Mesh%JMAX
          ! write coordinates to file
          WRITE (GetUnit(this), FMT='(ES13.6, ES14.6)', ADVANCE='NO', IOSTAT=ios) &
               Mesh%bcenter(i,j,1), Mesh%bcenter(i,j,2)
          ! write the data
          DO k=1,Physics%vnum
             WRITE (GetUnit(this), FMT=GetFormatString(this), ADVANCE='NO', IOSTAT=ios) &
                  ovar(i,j,k)
          END DO
          WRITE (GetUnit(this), FMT='(A)', IOSTAT=ios)
       END DO
       ! close 2D polar (fix for pm3d plots with gnuplot)
       IF ((GetType(Mesh%geometry).EQ.POLAR) &
            .AND.(GetType(Mesh%boundary(SOUTH)).EQ.PERIODIC)) THEN
          ! write coordinates to file
          WRITE (GetUnit(this), FMT='(ES13.6, ES14.6)', ADVANCE='NO', IOSTAT=ios) &
               Mesh%bcenter(i,Mesh%JMIN,1), (Mesh%bcenter(i,Mesh%JMIN,2)+2.*PI)
          ! write the data
          DO k=1,Physics%vnum
             WRITE (GetUnit(this), FMT=GetFormatString(this), ADVANCE='NO', IOSTAT=ios) &
                  ovar(i,Mesh%JMIN,k)
          END DO
          WRITE (GetUnit(this), FMT='(A)', IOSTAT=ios)
       END IF
       ! omit new line for 1-dimensional problems
       IF (Mesh%JMIN.NE.Mesh%JMAX) WRITE (GetUnit(this), "(A)", IOSTAT=ios)
    END DO

  END SUBROUTINE WriteOutput_gnuplot

END MODULE output_gnuplot
