!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: logio_binary.f90                                                  #
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
! Specific module for binary log data
!----------------------------------------------------------------------------!
MODULE logio_binary
  USE logio_common
  USE mesh_generic, ONLY : Mesh_TYP, GetType
  USE physics_generic, ONLY : Physics_TYP, GetType, Convert2Primitive
  USE timedisc_common, ONLY : Timedisc_TYP  
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: logformat_name = "binary"
  INTEGER, PARAMETER :: reclen   = 16            ! data record length        !
  ! program control records
  INTEGER, PARAMETER :: ctrlgeo  = 1             ! geometry                  !
  INTEGER, PARAMETER :: ctrlphys = 2             ! physics                   !
  INTEGER, PARAMETER :: ctrlnvar = 3             ! number of variables       !
  INTEGER, PARAMETER :: ctrlnx   = 4             ! number of cells in x-dir. !
  INTEGER, PARAMETER :: ctrlny   = 5             ! number of cells in y-dir. !
  ! computational domain records
  INTEGER, PARAMETER :: ctrlxmin = 11
  INTEGER, PARAMETER :: ctrlxmax = 12
  INTEGER, PARAMETER :: ctrlymin = 13
  INTEGER, PARAMETER :: ctrlymax = 14
  ! data set records
  INTEGER, PARAMETER :: trec     = 21            ! time stamp                !
  INTEGER, PARAMETER :: datarec  = 22            ! data records              !
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Logio_TYP, &
       ! constants
       ! methods 
       InitLogio, &
       InitLogio_binary, &
       ReadLogdata_binary, &
       WriteLogdata_binary, &
       WriteHeader_binary,&
       GetType, &
       GetName, &
       GetFilename, &
       GetLogstep, &
       CloseLogio_binary
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitLogio_binary(this,Mesh,Physics,logformat,logunit,filename,logdt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Logio_TYP)    :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: logformat,logunit
    INTEGER            :: logdt
    CHARACTER(LEN=*)   :: filename
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,Physics,logformat,logunit,filename,logdt
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    CALL InitLogio(this,logformat,logformat_name,logunit,filename,logdt)

  END SUBROUTINE InitLogio_binary

  
  SUBROUTINE ReadLogdata_binary(this,Mesh,Physics,Timedisc,ok)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Logio_TYP)    :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Timedisc_TYP) :: Timedisc
    LOGICAL            :: ok
    !------------------------------------------------------------------------!
    REAL               :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) :: cvar
    INTEGER            :: ios
    INTEGER            :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this,Mesh,Physics
    INTENT(OUT)        :: Timedisc
    !------------------------------------------------------------------------!
    
    ! open and rewind log file
    OPEN(GetUnit(this), IOSTAT=ios, FILE=GetFilename(this), STATUS="OLD", &
         ACCESS="DIRECT", RECL=reclen, ACTION="READ")

    IF (ios.LE.0) THEN
       ok = .TRUE.
       ! read log file data if possible
       CALL ReadHeader(this,Mesh,Physics,ok)
    
       IF (ok) THEN
          ! read timestamp
          READ (GetUnit(this),REC=trec,IOSTAT=ios) time
          IF (ios.GT.0) THEN 
             ok=.FALSE.
          ELSE
             ! read data
             DO j=0,(Mesh%JGMAX)-(Mesh%JGMIN)
                DO i=0,(Mesh%IGMAX)-(Mesh%IGMIN)
                   READ (GetUnit(this),REC=datarec+i+j*((Mesh%IGMAX)-(Mesh%IGMIN)),IOSTAT=ios) &
                        cvar(i+(Mesh%IGMIN),j+(Mesh%JGMIN),:)
                   IF (ios.GT.0) THEN
                      ok=.FALSE.
                      PRINT *, "WARING: Unable to read data from old log file "
                      EXIT
                   END IF
                END DO
             END DO
             IF (ok) THEN
                Timedisc%time = time
                Timedisc%cvar(:,:,:) = cvar(:,:,:)
                CALL Convert2Primitive(Physics,Mesh,Timedisc%cvar,Timedisc%pvar)
             END IF
          END IF
       END IF
    ELSE
       ok = .FALSE.
       PRINT *, "WARNING: Unable to open old log file ", &
            ACHAR(34), TRIM(GetFilename(this)), ACHAR(34)
    END IF

    ! close log file
    CLOSE(GetUnit(this))

  END SUBROUTINE ReadLogdata_binary


  SUBROUTINE WriteLogdata_binary(this,Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Logio_TYP)    :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Timedisc_TYP) :: Timedisc
    !------------------------------------------------------------------------!
    LOGICAL            :: ok
    INTEGER            :: ios
    INTEGER            :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this,Mesh,Physics,Timedisc
    !------------------------------------------------------------------------!

    OPEN(GetUnit(this), IOSTAT=ios, ACCESS="DIRECT", FILE=GetFilename(this), &
         STATUS="OLD", RECL=reclen, ACTION="WRITE")

    ok = .TRUE.
    IF (ios.LE.0) THEN
       ! write timestamp
       WRITE (GetUnit(this),REC=trec,IOSTAT=ios) Timedisc%time
       IF (ios.GT.0) ok = .FALSE.
       ! write data
       DO j=0,(Mesh%JGMAX)-(Mesh%JGMIN)
          DO i=0,(Mesh%IGMAX)-(Mesh%IGMIN)
             WRITE (GetUnit(this),REC=datarec+i+j*((Mesh%IGMAX)-(Mesh%IGMIN)),IOSTAT=ios) &
                  Timedisc%cvar(i+Mesh%IGMIN,j+Mesh%JGMIN,:)
             IF (ios.GT.0) THEN
                ok=.FALSE.
                EXIT
             END IF
          END DO
       END DO
    END IF

    CLOSE(GetUnit(this))

    IF (.NOT.ok) THEN
       PRINT *, "ERROR in WriteLogdata_binary: Unable to write log data to file"
       STOP
    END IF

  END SUBROUTINE WriteLogdata_binary


  SUBROUTINE WriteHeader_binary(this,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Logio_TYP)    :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    !------------------------------------------------------------------------!
    LOGICAL            :: ok
    INTEGER            :: ios
    INTEGER            :: i
    INTEGER            :: idata(10)
    REAL               :: rdata(10)
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this,Mesh,Physics
    !------------------------------------------------------------------------!

    ! create new file
    OPEN(GetUnit(this), IOSTAT=ios, ACCESS="DIRECT", RECL=reclen, &
         FILE=GetFilename(this), STATUS="NEW", ACTION="WRITE")

    IF (ios.LE.0) THEN
       ok = .TRUE.
       ! write header
       idata(1) = GetType(Mesh%geometry)
       idata(2) = GetType(Physics)
       idata(3) = Physics%vnum
       idata(4) = Mesh%INUM
       idata(5) = Mesh%JNUM
       idata(6:10) = 0   ! empty records
       
       rdata(1) = Mesh%xmin
       rdata(2) = Mesh%xmax
       rdata(3) = Mesh%ymin
       rdata(4) = Mesh%ymax
       rdata(5:10) = 0. ! empty records
       
       DO i=0,9
          WRITE (GetUnit(this),REC=ctrlgeo+i,IOSTAT=ios) idata(i+1)
          WRITE (GetUnit(this),REC=ctrlxmin+i,IOSTAT=ios) rdata(i+1)
          IF (ios.GT.0) THEN
             ok=.FALSE.
             EXIT
          END IF
       END DO
    ELSE
       ok = .FALSE.
    END IF

    CLOSE(GetUnit(this))

    IF (.NOT.ok) THEN
       PRINT *, "ERROR in WriteHeader_binary: Unable to write log file header"
       STOP
    END IF

  END SUBROUTINE WriteHeader_binary


  SUBROUTINE ReadHeader(this,Mesh,Physics,ok)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Logio_TYP)   :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    LOGICAL           :: ok
    !------------------------------------------------------------------------!
    INTEGER           :: ios
    INTEGER           :: i
    INTEGER           :: idata(10)
    REAL              :: rdata(10)
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Mesh,Physics
    INTENT(INOUT)     :: ok
    !------------------------------------------------------------------------!

    ! read header data
    DO i=0,9
       READ (GetUnit(this),REC=ctrlgeo+i,IOSTAT=ios) idata(i+1)
       READ (GetUnit(this),REC=ctrlxmin+i,IOSTAT=ios) rdata(i+1)
       ! abort on read error
       IF (ios.GT.0) THEN 
          ok = .FALSE.
          PRINT *, "WARNING: Unable to read header data from old log file"
       END IF
    END DO

    ! check header data
    IF (GetType(Mesh%geometry).NE.idata(1)) THEN
       ok = .FALSE.
       PRINT "(A,I2)", "WARNING: geometry mismatch in log file, should be ", idata(1)
    END IF
    IF (GetType(Physics).NE.idata(2)) THEN
       ok = .FALSE.
       PRINT "(A,I2)", "WARNING: physics mismatch in log file, should be ", idata(2)
    END IF
    IF (Physics%vnum.NE.idata(3)) THEN
       ok = .FALSE.
       PRINT "(A,I2)", "WARNING: vnum mismatch in log file, should be ", idata(3)
    END IF
    IF ((Mesh%INUM.NE.idata(4)).OR.(Mesh%JNUM.NE.idata(5))) THEN
       ok = .FALSE.
       PRINT "(A,I3,A,I3)", "WARNING: mesh resolution mismatch in log file, should be ", &
            idata(4), " x ", idata(5)
    END IF
    IF ((Mesh%xmin.NE.rdata(1)).OR.(Mesh%xmax.NE.rdata(2)) &
         .OR.(Mesh%ymin.NE.rdata(3)).OR.(Mesh%ymax.NE.rdata(4))) THEN
       ok = .FALSE.
       PRINT "(3(A),ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A)", &
            "WARNING: computational domain mismatch in log file,", ACHAR(10), &
            "         should be ", &
            rdata(1), " ..", rdata(2), "] x [", rdata(3), " ..", rdata(4), "]"
    END IF
  END SUBROUTINE ReadHeader


  SUBROUTINE CloseLogio_binary
    IMPLICIT NONE
    !------------------------------------------------------------------------!

  END SUBROUTINE CloseLogio_binary

END MODULE logio_binary
