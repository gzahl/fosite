!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: output_opendx.f90                                                 #
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
! Subroutines for OpenDX data output
!----------------------------------------------------------------------------!
MODULE output_opendx
  USE output_common
  USE physics_generic
  USE mesh_common, ONLY : Mesh_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: OUNIT = 12                    ! output unit number   !
  CHARACTER(LEN=32), PARAMETER :: format_name = "OpenDX"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Output_TYP, &
       ! constants
       APPEND, OVERWRITE, &
       ! methods 
       InitOutput_opendx, &
       WriteOutput_opendx, &
       GetType, &
       GetName, &
       GetFilename, &
       GetUnit, &
       CloseOutput_opendx
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitOutput_opendx(this,Mesh,Physics,oformat,fname,omode,start,end, &
       count)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP)   :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: oformat
    CHARACTER(LEN=*)   :: fname
    LOGICAL            :: omode
    REAL               :: start,end
    INTEGER            :: count
    !------------------------------------------------------------------------!
    LOGICAL            :: ex
    INTEGER            :: ios
    INTEGER            :: err
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,Physics,oformat,fname,omode,start,end,count
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    CALL InitOutput(this,oformat,format_name,OUNIT,fname,omode,start,end,count)

    ! output format string
    this%formatstring="(ES14.6)"

    ! write info to file header
    IF (this%openmode.NEQV.APPEND) THEN
       OPEN(GetUnit(this), IOSTAT=ios, FILE=GetFilename(this), STATUS="OLD", &
            POSITION="REWIND", ACTION="WRITE")
       IF (ios.GT.0) THEN
          PRINT *, "ERROR in InitOutput_opendx: Unable to open output file "
          CLOSE(GetUnit(this))
          STOP
       ELSE
          ! write positions information
          WRITE (GetUnit(this), &
               "(7(A),I5,I5,A, A,ES16.7,ES16.7,A, A,ES16.7,ES16.7,A,A,ES16.7,ES16.7)", &
               IOSTAT=ios) "# regular positions", ACHAR(10), &
               "object ", ACHAR(34), "grid", ACHAR(34), " class gridpositions counts ", &
               Mesh%INUM, Mesh%JNUM, ACHAR(10), &
               "origin ", Mesh%center(Mesh%IMIN,Mesh%JMIN,1), &
               Mesh%center(Mesh%IMIN,Mesh%JMIN,2), ACHAR(10), &
               "delta  ", Mesh%dx, 0.0, ACHAR(10), &
               "delta  ", 0.0, Mesh%dy
          IF (ios.GT.0) THEN
             PRINT *, "ERROR in InitOutput_opendx: Unable to write positions information "
             CLOSE(GetUnit(this))
             STOP
          END IF
          ! write connections information
          WRITE (GetUnit(this), "(7(A),I5,I5)", IOSTAT=ios) &
               "# regular connections", ACHAR(10), &
               "object ", ACHAR(34), "grid_con", ACHAR(34), &
               " class gridconnections counts ", Mesh%INUM, Mesh%JNUM
          CLOSE(GetUnit(this))
          IF (ios.GT.0) THEN
             PRINT *, "ERROR in InitOutput_opendx: Unable to write connections information "
             STOP
          END IF
       END IF
    END IF

  END SUBROUTINE InitOutput_opendx
  

  SUBROUTINE WriteOutput_opendx(this,Mesh,Physics,ovar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP)   :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) :: ovar
    !------------------------------------------------------------------------!
    INTEGER            :: ijnum
    INTEGER            :: i,j,k,n,kk
    CHARACTER(LEN=4)   :: cstep
    CHARACTER(LEN=12)  :: fmt
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this,Mesh,Physics,ovar
    !------------------------------------------------------------------------!
    ! time step information
    WRITE (cstep, "(I4.4)") this%step
    WRITE (GetUnit(this), "(A,A)") "# time step No. ", cstep

    ! total number of data points
    ijnum = Mesh%INUM * Mesh%JNUM

    ! loop over all data structures
    k=1
    DO n=1,SIZE(this%datastruct)
       IF (this%datastruct(n)%rank.EQ.0) THEN
          ! for scalar data
          WRITE (GetUnit(this), "(7(A),I2,A,I7,A)") &
               "object ", ACHAR(34), TRIM(this%datastruct(n)%name), &
               "data_", cstep, ACHAR(34), &
               " class array type float rank ", this%datastruct(n)%rank, &
               " items ", ijnum, " data follows"
       ELSE
          ! for vector data
          WRITE (GetUnit(this), "(7(A),I2,A,I2,A,I7,A)") &
               "object ", ACHAR(34), TRIM(this%datastruct(n)%name), &
               "data_", cstep, ACHAR(34), &
               " class array type float rank ", this%datastruct(n)%rank, &
               " shape ", this%datastruct(n)%shape, &
               " items ", ijnum, " data follows"
       END IF
       ! the data itself
       kk = k+this%datastruct(n)%shape-1
       WRITE (fmt, "(A,I2,A)") "(", this%datastruct(n)%shape, "(E15.7))"
       DO i=Mesh%IMIN,Mesh%IMAX
          DO j=Mesh%JMIN,Mesh%JMAX
             WRITE (this%outunit, fmt) ovar(i,j,k:kk)
          END DO
       END DO
       ! dependency of the data
       WRITE (GetUnit(this), "(8(A))") "attribute ", ACHAR(34), "dep", ACHAR(34), &
            " string ", ACHAR(34), "positions", ACHAR(34)
       ! put everything together
       WRITE (GetUnit(this), "(46(A))") &
            "object ", ACHAR(34), TRIM(this%datastruct(n)%name), &
            "_", cstep, ACHAR(34), " class field", ACHAR(10), &
            "component ", ACHAR(34), "positions", ACHAR(34), &
            " value ", ACHAR(34), "grid", ACHAR(34), ACHAR(10), &
            "component ", ACHAR(34), "connections", ACHAR(34), &
            " value ", ACHAR(34), "grid_con", ACHAR(34), ACHAR(10), &
            "component ", ACHAR(34), "data", ACHAR(34), &
            " value ", ACHAR(34), TRIM(this%datastruct(n)%name), "data_", &
            cstep, ACHAR(34), ACHAR(10), &
            "attribute ", ACHAR(34), "name", ACHAR(34), &
            " string ", ACHAR(34), TRIM(this%datastruct(n)%name), ACHAR(34), ACHAR(10)
       k = kk+1
    END DO
  END SUBROUTINE WriteOutput_opendx


  SUBROUTINE CloseOutput_opendx(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Output_TYP)   :: this
    !------------------------------------------------------------------------!
    CHARACTER(LEN=4)   :: cm
    INTEGER            :: n,m
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this
    !------------------------------------------------------------------------!
    ! create the time series
    DO n=1,SIZE(this%datastruct)
       WRITE (GetUnit(this), "(5(A))") &
            "object ", ACHAR(34), TRIM(this%datastruct(n)%name), ACHAR(34), &
            " class series"
       DO m=0, this%step
          WRITE (cm, "(I4.4)") m
          WRITE (GetUnit(this), "(A,I4,6(A))") &
               "member ", m, " value ", ACHAR(34), &
               TRIM(this%datastruct(n)%name), "_", cm, ACHAR(34)
       END DO
       WRITE (GetUnit(this), "(9(A))") &
            "attribute ", ACHAR(34), "name", ACHAR(34), &
            " string ", ACHAR(34), TRIM(this%datastruct(n)%name), ACHAR(34), ACHAR(10)
    END DO
    ! put everything together
    WRITE (GetUnit(this), "(5(A))") &
         "object ", ACHAR(34), "default", ACHAR(34), " class group"
    DO n=1,SIZE(this%datastruct)
       WRITE (GetUnit(this), "(8(A))") &
            "member ", ACHAR(34), TRIM(this%datastruct(n)%name), ACHAR(34), &
            " value ", ACHAR(34), TRIM(this%datastruct(n)%name), ACHAR(34)
    END DO
    WRITE (GetUnit(this), "(A)") "end"
  END SUBROUTINE CloseOutput_opendx

END MODULE output_opendx
