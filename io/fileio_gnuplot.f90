!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fileio_gnuplot.f90                                                #
!#                                                                           #
!# Copyright (C) 2008 Tobias Illenseer <tillense@astrophysik.uni-kiel.de>    #
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
! module for GNUPLOT readable tabular file I/O
!----------------------------------------------------------------------------!
MODULE fileio_gnuplot
  USE common_types, ONLY : Error_common => Error
  USE fileio_common, InitFileIO_common => InitFileIO, Error_fileio => Error
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  USE timedisc_common, ONLY : Timedisc_TYP
  IMPLICIT NONE
#ifdef PARALLEL
  include 'mpif.h'
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE Error
     MODULE PROCEDURE Error_gnuplot, Error_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  ! some string lengths
  INTEGER, PARAMETER   :: HLEN = 25            ! header                      !
  INTEGER, PARAMETER   :: FLEN = 12            ! one data field              !
  ! some special strings
  CHARACTER, PARAMETER :: LF = ACHAR(10)       ! line feed                   !
  CHARACTER, PARAMETER :: SP = ACHAR(32)       ! space                       !
  CHARACTER*2, PARAMETER :: RECSEP = SP // SP  ! data record seperator       !
  CHARACTER*2, PARAMETER :: LINSEP = SP // LF  ! line seperator              !
  CHARACTER*2, PARAMETER :: BLKSEP = LF // LF  ! block seperator             !
  CHARACTER(LEN=HLEN), PARAMETER :: &          ! the header string           !
       header_string = "# Data output of fosite" // LF // LF
 !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       FileIO_TYP, &
       ! constants
       READONLY, READEND, REPLACE, APPEND, &
       ASCII, BIN, &
       FILE_EXISTS, &
#ifdef PARALLEL
       DEFAULT_MPI_REAL,&
#endif
       ! methods
       InitFileIO, &
       InitFileIO_gnuplot, &
       CloseFileIO_gnuplot, &
       WriteHeader_gnuplot, &
       ReadHeader_gnuplot, &
       WriteTimestamp_gnuplot,&
       ReadTimestamp_gnuplot,&
       WriteDataset_gnuplot, &
       ReadDataset_gnuplot, &
       OpenFile, &
       OpenFile_gnuplot, &
       CloseFile_gnuplot, &
       RewindFile, &
       AdjustTimestep, &
       IncTime, &
       GetFilename, &
       GetFilestatus, &
       GetType, &
       GetName, &
       GetRank, &
       GetNumProcs, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS
  
  SUBROUTINE InitFileIO(this,Mesh,Physics,fmt,fmtname,filename,extension, &
       stoptime,dtwall,count,fcycles)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    INTEGER           :: fmt
    CHARACTER(LEN=*)  :: fmtname,filename,extension
    REAL              :: stoptime
    INTEGER           :: dtwall
    INTEGER           :: count
    INTEGER           :: fcycles
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,fmt,fmtname,filename,extension,stoptime, &
         dtwall,count,fcycles
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    ! basic FileIO initialization
    CALL InitFileIO_common(this,fmt,fmtname,filename,extension,fcycles)
    this%stoptime = stoptime
    this%dtwall   = dtwall
    this%time     = 0.
    this%count    = count
    this%step     = 0
    ! count the number of output columns, i.e. fields per data point
    IF ((Mesh%INUM.EQ.1).OR.(Mesh%JNUM.EQ.1)) THEN
       ! 1D mesh
       this%cols = 1 + Physics%vnum
    ELSE
       ! 2D mesh
       this%cols = 2 + Physics%vnum
    END IF
#ifdef PARALLEL
    ! check data type extents in files
    ! first try to create a new dummy file
    CALL MPI_File_open(MPI_COMM_WORLD,GetFilename(this),IOR(IOR(MPI_MODE_RDWR,&
         MPI_MODE_CREATE),IOR(MPI_MODE_EXCL,MPI_MODE_DELETE_ON_CLOSE)),&
         MPI_INFO_NULL,this%handle,this%error)
    ! maybe file exists
    IF (this%error.NE.0) CALL MPI_File_open(MPI_COMM_WORLD,GetFilename(this),&
         MPI_MODE_RDONLY,MPI_INFO_NULL,this%handle,this%error)
    ! then check the data type sizes
    ! extent of integer in file
    IF (this%error.EQ.0) CALL MPI_File_get_type_extent(this%handle,MPI_INTEGER,&
         this%intext,this%error)
    ! extent of real in file
    IF (this%error.EQ.0) CALL MPI_File_get_type_extent(this%handle,DEFAULT_MPI_REAL,&
         this%realext,this%error)
    IF (this%error.NE.0) CALL Error(this,"InitFileIO","unable to check file properties")
    CALL MPI_File_close(this%handle,this%error)
#endif
  END SUBROUTINE InitFileIO


  SUBROUTINE InitFileIO_gnuplot(this,Mesh,Physics,fmt,filename,stoptime,dtwall,&
       count,fcycles)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    INTEGER           :: fmt
    CHARACTER(LEN=*)  :: filename
    REAL              :: stoptime
    INTEGER           :: dtwall
    INTEGER           :: count
    INTEGER           :: fcycles
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER           :: i
    INTEGER, DIMENSION(Mesh%IMAX-Mesh%IMIN+1) :: blocklen,indices
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,fmt,filename,stoptime,dtwall,count,fcycles
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL InitFileIO(this,Mesh,Physics,fmt,"GNUPLOT",filename,"dat",stoptime,&
         dtwall,count,fcycles)

    ! length of one output line
    this%linelen = this%cols * FLEN

#ifdef PARALLEL
    ! create new data type handle for one line
    CALL MPI_Type_contiguous(this%linelen,MPI_CHARACTER,this%basictype,this%error)
    CALL MPI_Type_commit(this%basictype,this%error)

    ! number of output blocks
    this%blocknum = Mesh%IMAX - Mesh%IMIN + 1
    ! size of the output buffer
    this%bufsize  = Mesh%JMAX - Mesh%JMIN + 1

    ! allocate memory for output buffer and displacement records
    ALLOCATE(this%outbuf(this%linelen,Mesh%JMIN:Mesh%JMAX),&
         STAT=this%error)
    IF (this%error.NE.0) THEN
       CALL Error(this,"InitFileIO_gnuplot","Unable to allocate memory.")
    END IF

    blocklen(:) = this%bufsize
    DO i=Mesh%IMIN,Mesh%IMAX
        indices(i-Mesh%IMIN+1) = (i-1)*Mesh%JNUM + Mesh%JMIN - 1
    END DO

    ! new file type for the staggered data
    CALL MPI_Type_indexed(this%blocknum,blocklen,indices,&
         this%basictype,this%filetype,this%error)
    CALL MPI_Type_commit(this%filetype,this%error)
#endif
    ! write the format string for one entry in the data file:
    ! FLEN-2 characters for the number and 2 for the separators
    WRITE (this%fmtstr,'(A3,I2,A,I1,A5)') '(ES', FLEN-2, '.', FLEN-9,',A,A)'
    ! write format string for one output line
    WRITE (this%linefmt, '(A,I0,A)') "(A", this%linelen-1, ")"
  END SUBROUTINE InitFileIO_gnuplot


 PURE SUBROUTINE AdjustTimestep(this,time,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    REAL             :: time,dt
    !------------------------------------------------------------------------!
    INTENT(IN)       :: this
    INTENT(INOUT)    :: time,dt
    !------------------------------------------------------------------------!
    IF ((time+dt)/this%time.GT.1.0) THEN
       dt = this%time - time
    ELSE IF((time+1.5*dt)/this%time.GT.1.0) THEN
       dt = 0.5*(this%time - time)
    END IF
  END SUBROUTINE AdjustTimestep


  PURE SUBROUTINE IncTime(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    this%time = this%time + ABS(this%stoptime) / this%count
    this%step = this%step + 1
  END SUBROUTINE IncTime


  SUBROUTINE OpenFile(this,action,fformat)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    INTEGER          :: action
    CHARACTER(LEN=*) :: fformat
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)       :: action,fformat
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    SELECT CASE(action)
    CASE(READONLY)
#ifdef PARALLEL
       CALL MPI_File_open(MPI_COMM_WORLD,GetFilename(this),MPI_MODE_RDONLY, &
            MPI_INFO_NULL,this%handle,this%error)
#else
       OPEN(this%unit,FILE=GetFilename(this),FORM=fformat,STATUS="OLD", &
            ACTION="READ",POSITION="REWIND",IOSTAT=this%error)
#endif
    CASE(READEND)
#ifdef PARALLEL
       CALL MPI_File_open(MPI_COMM_WORLD,GetFilename(this),IOR(MPI_MODE_RDONLY,&
            MPI_MODE_APPEND),MPI_INFO_NULL,this%handle,this%error)
       ! opening in append mode doesn't seem to work for pvfs2, hence ...
       offset = 0
       CALL MPI_File_seek(this%handle,offset,MPI_SEEK_END,this%error)
       CALL MPI_File_sync(this%handle,this%error)
#else
       OPEN(this%unit,FILE=GetFilename(this),FORM=fformat,STATUS="OLD", &
            ACTION="READ",POSITION="APPEND",IOSTAT=this%error)
#endif
    CASE(REPLACE)
#ifdef PARALLEL
       CALL MPI_File_delete(GetFilename(this),MPI_INFO_NULL,this%error)
       CALL MPI_File_open(MPI_COMM_WORLD,GetFilename(this),IOR(MPI_MODE_WRONLY,&
            MPI_MODE_CREATE),MPI_INFO_NULL,this%handle,this%error)
#else
       OPEN(this%unit,FILE=GetFilename(this),FORM=fformat,STATUS="REPLACE",&
            ACTION="WRITE",POSITION="REWIND",IOSTAT=this%error)
#endif
    CASE(APPEND)
#ifdef PARALLEL
       CALL MPI_File_open(MPI_COMM_WORLD,GetFilename(this),IOR(MPI_MODE_RDWR,&
            MPI_MODE_APPEND),MPI_INFO_NULL,this%handle,this%error)       
       ! opening in append mode doesn't seem to work for pvfs2, hence ...
       offset = 0
       CALL MPI_File_seek(this%handle,offset,MPI_SEEK_END,this%error)
       CALL MPI_File_sync(this%handle,this%error)
#else
       OPEN(this%unit,FILE=GetFilename(this),FORM=fformat,STATUS="OLD",&
            ACTION="READWRITE",POSITION="APPEND",IOSTAT=this%error)
#endif
    CASE DEFAULT
       CALL Error(this,"OpenFile","Unknown access mode.")
    END SELECT
  END SUBROUTINE OpenFile


  SUBROUTINE RewindFile(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    this%offset = 0
    CALL MPI_File_seek(this%handle,this%offset,MPI_SEEK_SET,this%error)
#else
    REWIND (UNIT=this%unit,IOSTAT=this%error)
#endif
    IF (this%error.NE.0) CALL Error(this,"RewindFile",&
         "Cannot rewind data file")
  END SUBROUTINE RewindFile


  SUBROUTINE OpenFile_gnuplot(this,action)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    INTEGER          :: action
    !------------------------------------------------------------------------!
    INTENT(IN)       :: action
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    CALL OpenFile(this,action,ASCII)
  END SUBROUTINE OpenFile_gnuplot


  SUBROUTINE CloseFile_gnuplot(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    CALL MPI_File_close(this%handle,this%error)
#else
    CLOSE(this%unit,IOSTAT=this%error)
#endif
  END SUBROUTINE CloseFile_gnuplot


  SUBROUTINE WriteHeader_gnuplot(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    IF (GetRank(this).EQ.0) THEN
#ifdef PARALLEL
       CALL MPI_File_write(this%handle,header_string,HLEN, &
            MPI_CHARACTER,this%status,this%error)
#else
       WRITE (this%unit,FMT='(A)',IOSTAT=this%error) header_string(1:HLEN-1)
#endif
    END IF
  END SUBROUTINE WriteHeader_gnuplot


  SUBROUTINE ReadHeader_gnuplot(this,success)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    LOGICAL          :: success
    !------------------------------------------------------------------------!
    INTENT(OUT)      :: success
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    IF (GetRank(this).EQ.0) THEN
    END IF
#else
#endif
    success = .FALSE.
  END SUBROUTINE ReadHeader_gnuplot


  SUBROUTINE WriteTimestamp_gnuplot(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    REAL             :: time
    !------------------------------------------------------------------------!
    INTENT(IN)       :: time
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    IF (GetRank(this).EQ.0) THEN
#ifdef PARALLEL
#else
#endif
    END IF
!!$    CALL Warning(this,"WriteTimestamp_gnuplot",&
!!$         "function is not implemented")
  END SUBROUTINE WriteTimestamp_gnuplot


  SUBROUTINE ReadTimestamp_gnuplot(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    REAL             :: time
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: time
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    IF (GetRank(this).EQ.0) THEN
#ifdef PARALLEL
#else
#endif
    END IF
    time = 0.0
!!$    CALL Warning(this,"ReadTimestamp_gnuplot",&
!!$         "function is not implemented")
  END SUBROUTINE ReadTimestamp_gnuplot


  SUBROUTINE WriteDataset_gnuplot(this,Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    INTEGER          :: i,j,k,l
#ifdef PARALLEL
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
    INTEGER          :: request
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)       :: Mesh,Physics
    INTENT(INOUT)    :: this,Timedisc
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    ! be sure to write at the end by getting the offset from the file's size
    CALL MPI_File_get_size(this%handle,offset,this%error)
    ! very importan
    CALL MPI_Barrier(MPI_COMM_WORLD,this%error)
    ! write _one_ line feed at the beginning of each time step
    IF (GetRank(this).EQ.0) THEN
       CALL MPI_File_write_at(this%handle, offset, LF, 1, MPI_CHARACTER, &
            this%status, this%error)
    END IF
    ! add the initial line feed and the general offset (depends on Mesh%IMIN)
    offset = offset + 1
    ! create the file view
    CALL MPI_File_set_view(this%handle,offset,this%basictype,this%filetype, &
         'native',MPI_INFO_NULL,this%error)
#else
    ! write _one_ line feed at the beginning of each time step
    WRITE (this%unit,FMT='(A)',ADVANCE='NO') LF
#endif
    ! trim the data for gnuplot output
    WHERE (ABS(Timedisc%pvar(:,:,:)).LT.(MAX(TINY(1.0),1.0D-99)))
       Timedisc%pvar(:,:,:) = 0.0E+00
    END WHERE
    DO i=Mesh%IMIN,Mesh%IMAX
       DO j=Mesh%JMIN,Mesh%JMAX
          ! write positions to line buffer
          l = 0
          IF (Mesh%INUM.GT.1) THEN
             l = l + 1
             WRITE (this%linebuf((l-1)*FLEN+1:l*FLEN),TRIM(this%fmtstr)) &
                  Mesh%bcenter(i,j,1), RECSEP
          END IF
          IF (Mesh%JNUM.GT.1) THEN
             l = l + 1
             WRITE (this%linebuf((l-1)*FLEN+1:l*FLEN),TRIM(this%fmtstr)) &
                  Mesh%bcenter(i,j,2), RECSEP
          END IF

          ! write variables to line buffer
          DO k=l+1,l+Physics%vnum-1
             WRITE (this%linebuf((k-1)*FLEN+1:k*FLEN),TRIM(this%fmtstr)) &
                  Timedisc%pvar(i,j,k-l), RECSEP
          END DO

          IF ((j.EQ.Mesh%JNUM).AND.((Mesh%JNUM.GT.1).OR.(Mesh%INUM.EQ.i))) THEN
             ! finish the block
             WRITE (this%linebuf((this%cols-1)*FLEN+1:this%cols*FLEN),TRIM(this%fmtstr)) &
                  Timedisc%pvar(i,j,Physics%vnum), BLKSEP
          ELSE
             ! finish the line
             WRITE (this%linebuf((this%cols-1)*FLEN+1:this%cols*FLEN),TRIM(this%fmtstr)) &
                  Timedisc%pvar(i,j,Physics%vnum), LINSEP
          END IF

#ifdef PARALLEL
          ! write line buffer to output buffer
          DO k=1,this%linelen
             this%outbuf(k,j) = this%linebuf(k:k)
          END DO
#else
          ! write line buffer to output file
          WRITE (this%unit,FMT=TRIM(this%linefmt),ADVANCE='YES') this%linebuf(1:this%linelen-1)
#endif
       END DO
#ifdef PARALLEL
       !*****************************************************************!
       ! This collective call doesn't work for pvfs2 -> bug in ROMIO ?
!!$       CALL MPI_File_write_all(this%handle,this%binout,this%bufsize,&  
!!$            this%basictype, this%status, this%error)
       !*****************************************************************!
       ! so we use these two commands instead
       CALL MPI_File_iwrite(this%handle,this%outbuf,this%bufsize,this%basictype,&
            request,this%error)
       CALL MPI_Wait(request,this%status,this%error)
#endif
    END DO
  END SUBROUTINE WriteDataset_gnuplot


  SUBROUTINE ReadDataset_gnuplot(this,Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Timedisc
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
  END SUBROUTINE ReadDataset_gnuplot


  SUBROUTINE Error_gnuplot(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(INOUT) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL CloseFile_gnuplot(this)
    CALL Error_fileio(this,modproc,msg)
  END SUBROUTINE Error_gnuplot


  SUBROUTINE CloseFileIO_gnuplot(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    DEALLOCATE(this%outbuf)
#endif
  END SUBROUTINE CloseFileIO_gnuplot

END MODULE fileio_gnuplot
