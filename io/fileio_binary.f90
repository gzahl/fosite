!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fileio_binary.f90                                                 #
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
! module for BINARY file I/O
!----------------------------------------------------------------------------!
MODULE fileio_binary
  USE fileio_gnuplot, CloseFile_binary => CloseFile_gnuplot
  USE geometry_common, ONLY : Geometry_TYP, GetType
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP, GetType
  USE timedisc_common, ONLY : Timedisc_TYP
  IMPLICIT NONE
#ifdef PARALLEL
  include 'mpif.h'
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
  ! size of header data fields
  INTEGER, PARAMETER :: HISIZE = 10              ! integer data              !
  INTEGER, PARAMETER :: HRSIZE = 10              ! real data                 !
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       FileIO_TYP, &
       ! constants
       ! methods
       InitFileIO_binary, &
       OpenFile_binary, &
       CloseFile_binary, &
       WriteHeader_binary, &
       ReadHeader_binary, &
       WriteTimestamp_binary,&
       ReadTimestamp_binary,&
       WriteDataset_binary, &
       ReadDataset_binary, &
       CloseFileIO_binary
  !--------------------------------------------------------------------------!

CONTAINS
  
  SUBROUTINE InitFileIO_binary(this,Mesh,Physics,fmt,filename,stoptime,dtwall,&
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
#ifdef PARALLEL
    INTEGER, DIMENSION(2) :: gsizes,lsizes,indices
    INTEGER           :: lb,extent
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,fmt,filename,stoptime,count,fcycles
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL InitFileIO(this,Mesh,Physics,fmt,"binary",filename,"bin",stoptime,&
         dtwall,count,fcycles)
    ! allocate memory
    ALLOCATE(this%header%idata(HISIZE),this%header%rdata(HRSIZE),&
         ! for output buffer
         this%binout(1:this%cols,Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX),&
#ifdef PARALLEL
         ! and displacements
         this%disp(Mesh%IMIN:Mesh%IMAX),&
#endif
         STAT=this%error)
    IF (this%error.NE.0) THEN
       CALL Error(this,"InitFileIO_binary","Unable to allocate memory.")
    END IF
#ifdef PARALLEL
    ! create the basic type: (1 or 2) coordinates + (4 or 5) simulation data
    CALL MPI_Type_contiguous(this%cols,DEFAULT_MPI_REAL,this%basictype,this%error)
    CALL MPI_Type_commit(this%basictype,this%error)

    ! create the data type for the distributed array of
    ! coordinates and simulation data
    gsizes(1) = Mesh%INUM
    gsizes(2) = Mesh%JNUM
    lsizes(1) = Mesh%IMAX-Mesh%IMIN+1
    lsizes(2) = Mesh%JMAX-Mesh%JMIN+1
    indices(1)= Mesh%IMIN-1
    indices(2)= Mesh%JMIN-1
    this%bufsize = lsizes(1) * lsizes(2)
    CALL MPI_Type_create_subarray(2, gsizes, lsizes, indices,MPI_ORDER_FORTRAN,&
         this%basictype,this%filetype,this%error)
    CALL MPI_Type_commit(this%filetype,this%error)
#endif
    ! initialize header data fields
    ! 1. integer data:
    this%header%idata(1) = GetType(Physics)
    this%header%idata(2) = GetType(Mesh%Geometry)
    this%header%idata(3) = Mesh%INUM
    this%header%idata(4) = Mesh%JNUM
    this%header%idata(5:HISIZE) = 0      ! insert additional parameters here !
    ! 2. real data:
    this%header%rdata(1) = Mesh%xmin
    this%header%rdata(2) = Mesh%xmax
    this%header%rdata(3) = Mesh%ymin
    this%header%rdata(4) = Mesh%ymax
    this%header%rdata(5:HRSIZE) = 0.0    ! insert additional parameters here !
  END SUBROUTINE InitFileIO_binary


  SUBROUTINE OpenFile_binary(this,action)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    INTEGER          :: action
    !------------------------------------------------------------------------!
    INTENT(IN)       :: action
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    CALL OpenFile(this,action,BIN)
  END SUBROUTINE OpenFile_binary


  SUBROUTINE SeekBackwards(this,count,error)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    INTEGER, OPTIONAL :: count
    INTEGER, OPTIONAL :: error
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER           :: k
    INTEGER(KIND=4)   :: size1,size2  ! force 4 byte integer
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
#endif
    INTEGER           :: i,j,n,err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: count
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: error
    !------------------------------------------------------------------------!
    IF (PRESENT(count)) THEN
       n = count
    ELSE
       n = 1
    END IF
    err=0
#ifdef PARALLEL
    ! get the current position of the individual file pointer
    ! each process may return a different result
    CALL MPI_File_get_position(this%handle,this%offset,err)
    IF (this%error.Eq.0) CALL MPI_File_get_byte_offset(this%handle,this%offset,&
         offset,err)
#endif
    ! loop over unformatted data sets
    DO i=1,n
       ! abort if reading fails
       IF (err.NE.0) EXIT
#ifdef PARALLEL
       ! read both subsequent and preceding size information
       ! of the data record
       size1 = 0
       DO j=1,2
          size2 = size1
          ! compute address of the size information record
          offset = offset - (4+size1)
          IF (offset.LT.0) THEN
             err=-1
             EXIT
          ELSE
             CALL MPI_File_read_at(this%handle,offset,size1,1,MPI_INTEGER4,this%status,err)
          END IF
       END DO
       IF ((err.NE.0).OR.(size1.NE.size2)) err=-1
#else
       ! one backspace
       BACKSPACE (UNIT=this%unit,IOSTAT=err)
#endif
    END DO
#ifdef PARALLEL
    IF (err.EQ.0) THEN
       this%offset = offset
       CALL MPI_File_seek(this%handle,offset,MPI_SEEK_SET,err)
    END IF
#endif
    IF (PRESENT(error)) error=err
  END SUBROUTINE SeekBackwards


  SUBROUTINE WriteHeader_binary(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    !------------------------------------------------------------------------!
    INTEGER(KIND=4)  :: size  ! force 4 byte integer
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    IF (GetRank(this).EQ.0) THEN
#ifdef PARALLEL
       ! size in bytes of the data (for compatibility with standard Fortran i/o)
       size = HISIZE*this%intext + HRSIZE*this%realext
       CALL MPI_File_write(this%handle,size,1,MPI_INTEGER4,this%status,this%error)
       ! integer data
       IF (this%error.EQ.0) CALL MPI_File_write(this%handle,this%header%idata,&
            HISIZE,MPI_INTEGER,this%status,this%error)
       ! real data
       IF (this%error.EQ.0) CALL MPI_File_write(this%handle,this%header%rdata,&
            HRSIZE,DEFAULT_MPI_REAL,this%status,this%error)
       ! write the size again to finish the i/o of the unformatted data
       IF (this%error.EQ.0) CALL MPI_File_write(this%handle,size,1,MPI_INTEGER4,&
            this%status,this%error)
#else
       WRITE (UNIT=this%unit,IOSTAT=this%error) this%header%idata, &
            this%header%rdata
#endif
       IF (this%error.NE.0) CALL Error(this,"WriteHeader_binary",&
            "Cannot write data file header")
    END IF
  END SUBROUTINE WriteHeader_binary


  SUBROUTINE ReadHeader_binary(this,success)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    LOGICAL          :: success
    !------------------------------------------------------------------------!
    INTEGER(KIND=4)  :: size1, size2    ! force 4 byte integer
    INTEGER, DIMENSION(HISIZE) :: idata
    REAL, DIMENSION(HRSIZE)    :: rdata
    !------------------------------------------------------------------------!
    INTENT(OUT)      :: success
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    success = .TRUE.
    this%error = 0
    IF (GetRank(this).EQ.0) THEN
    ! read the data
#ifdef PARALLEL
       size2 = this%intext*HISIZE + this%realext*HRSIZE
       ! size of the header data in bytes
       CALL MPI_File_read(this%handle,size1,1,MPI_INTEGER4,this%status,this%error)
       IF ((this%error.NE.0).OR.(size1.NE.size2)) this%error=-1
       ! integer data
       IF (this%error.EQ.0) CALL MPI_File_read(this%handle,idata,HISIZE,&
            MPI_INTEGER,this%status,this%error)
       ! real data
       IF (this%error.EQ.0) CALL MPI_File_read(this%handle,rdata,HRSIZE,&
            DEFAULT_MPI_REAL,this%status,this%error)
       ! size of the header data in bytes
       IF (this%error.EQ.0) CALL MPI_File_read(this%handle,size2,1,MPI_INTEGER4,&
            this%status,this%error)
       IF ((this%error.NE.0).OR.(size1.NE.size2)) THEN
          CALL Warning(this,"ReadHeader_binary","header size mismatch")
          success = .FALSE.
          RETURN
       END IF
#else
       READ (UNIT=this%unit,IOSTAT=this%error) idata, rdata
#endif
       IF (this%error.NE.0) THEN
          CALL Warning(this,"ReadHeader_binary",&
               "File header corrupt")
          success = .FALSE. 
          RETURN
       END IF
       ! check header data
       IF (idata(1).NE.this%header%idata(1)) THEN
          CALL Warning(this,"ReadHeader_binary","physics mismatch")
          success = .FALSE.
       END IF
       IF (idata(2).NE.this%header%idata(2)) THEN
          CALL Warning(this,"ReadHeader_binary","geometry mismatch")
          success = .FALSE.
       END IF
       IF (ALL(idata(3:4).NE.this%header%idata(3:4))) THEN
          CALL Warning(this,"ReadHeader_binary","resolution mismatch")
          success = .FALSE.
       END IF
       IF (ALL(rdata(1:4).NE.this%header%rdata(1:4))) THEN
          CALL Warning(this,"ReadHeader_binary","computational domain mismatch")
          success = .FALSE.
       END IF
    END IF
  END SUBROUTINE ReadHeader_binary


  SUBROUTINE WriteTimestamp_binary(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    REAL             :: time
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER(KIND=4)  :: size  ! force 4 byte integer
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)       :: time
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    IF (GetRank(this).EQ.0) THEN
#ifdef PARALLEL
       size = this%realext
       ! write size of data (for compatibility with standard Fortran I/O)
       CALL MPI_File_write(this%handle,size,1,MPI_INTEGER4,this%status,this%error)
       ! write the time stamp
       IF (this%error.EQ.0) CALL MPI_File_write(this%handle,time,1,&
            DEFAULT_MPI_REAL,this%status,this%error)
       ! write size of data (for compatibility with standard Fortran I/O)
       IF (this%error.EQ.0) CALL MPI_File_write(this%handle,size,1,&
            MPI_INTEGER4,this%status,this%error)
#else
       ! time stamp
       WRITE (this%unit,IOSTAT=this%error) time
#endif
       IF (this%error.NE.0) CALL Error(this,"WriteTimestamp_binary", &
            "writing time stamp failed")
    END IF
  END SUBROUTINE WriteTimestamp_binary


  SUBROUTINE ReadTimestamp_binary(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    REAL             :: time
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER(KIND=4) :: size
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
#endif
    !------------------------------------------------------------------------!
    INTENT(OUT)      :: time
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    IF (GetRank(this).EQ.0) THEN
       ! return a negative time if reading fails
       time=-1.0
       ! seek backwards one data record
       CALL SeekBackwards(this,1,this%error)
       IF (this%error.EQ.0) THEN
#ifdef PARALLEL
          offset = this%offset + 4
          CALL MPI_File_read_at(this%handle,offset,time,1,DEFAULT_MPI_REAL,&
               this%status,this%error)
#else
          READ (this%unit,IOSTAT=this%error) time
#endif
       END IF
       IF (this%error.NE.0) CALL Warning(this,"ReadTimestamp_binary",&
            "can't read the time stamp")
    END IF
#ifdef PARALLEL
    ! send the time stamp to all other processes
    CALL MPI_Bcast(time,1,DEFAULT_MPI_REAL,0,MPI_COMM_WORLD,this%error)
#endif
  END SUBROUTINE ReadTimestamp_binary


  SUBROUTINE WriteDataset_binary(this,Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    INTEGER          :: i,j,kx,ky
#ifdef PARALLEL
    INTEGER(KIND=4)  :: size  ! force 4 byte integer
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
    INTEGER          :: request
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)       :: Mesh,Physics,Timedisc
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    IF (Mesh%INUM.GT.1) THEN
       kx=1
    ELSE
       kx=0
    END IF
    IF (Mesh%JNUM.GT.1) THEN
       ky=kx+1
    ELSE
       ky=kx
    END IF
    DO j=Mesh%JMIN,Mesh%JMAX
       DO i=Mesh%IMIN,Mesh%IMAX
          ! copy coordinates
          IF (kx.GT.0) this%binout(kx,i,j) = Mesh%bcenter(i,j,1)
          IF (ky.GT.kx) this%binout(ky,i,j) = Mesh%bcenter(i,j,2)
          ! copy data
          this%binout(ky+1:ky+Physics%vnum,i,j) = Timedisc%pvar(i,j,1:Physics%vnum)
       END DO
    END DO

    ! write data
#ifdef PARALLEL
    CALL MPI_File_get_size(this%handle,this%offset,this%error)
    CALL MPI_Barrier(MPI_COMM_WORLD,this%error)
    ! compute the number of bytes in real numbers we are going to
    ! write with _all_ processes in total 
    size=(Mesh%INUM*Mesh%JNUM*this%cols) * this%realext
    IF (GetRank(this).EQ.0) THEN
       ! write size information for compatiblity with standard Fortran I/O
       CALL MPI_File_write_at(this%handle,this%offset,size,1,MPI_INTEGER4,&
            this%status,this%error)
    END IF
    ! skip the size information
    this%offset = this%offset + 4
    CALL MPI_File_set_view(this%handle,this%offset,this%basictype,&
         this%filetype, 'native', MPI_INFO_NULL, this%error)
    !*****************************************************************!
    ! This collective call doesn't work for pvfs2 -> bug in ROMIO ?
!!$    CALL MPI_File_write_all(this%handle,this%binout,this%bufsize,&  
!!$         this%basictype, this%status, this%error)
    !*****************************************************************!
    CALL MPI_File_iwrite(this%handle,this%binout,this%bufsize,this%basictype,&
         request,this%error)
    CALL MPI_Wait(request,this%status,this%error)
    this%offset = this%offset + size
    CALL MPI_File_set_view(this%handle,this%offset,MPI_INTEGER4,&
         MPI_INTEGER4, 'native', MPI_INFO_NULL, this%error)
    IF (GetRank(this).EQ.0) THEN
       ! write size information again at the end
       CALL MPI_File_write(this%handle,size,1,MPI_INTEGER4,&
            this%status,this%error)
    END IF
#else
    WRITE (this%unit) this%binout(:,:,:)
#endif
  END SUBROUTINE WriteDataset_binary


  SUBROUTINE ReadDataset_binary(this,Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k,k0
#ifdef PARALLEL
    INTEGER(KIND=4)  :: size1, size2  ! force 4 byte integer
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
    INTEGER          :: request
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: this,Timedisc
    !------------------------------------------------------------------------!
    ! point on data record
    CALL SeekBackwards(this,2,this%error)
    IF (this%error.EQ.0) THEN
#ifdef PARALLEL
       ! compute the number of bytes in real numbers we are going to
       ! read with _all_ processes in total
       size1=(Mesh%INUM*Mesh%JNUM*this%cols) * this%realext
       IF (GetRank(this).EQ.0) THEN
          ! check extent of data in file
          CALL MPI_File_read_at(this%handle,this%offset,size2,1,&
               MPI_INTEGER4,this%status,this%error)
          IF ((this%error.NE.0).OR.(size1.NE.size2)) CALL Error(this,&
               "ReadDataset_binary","data size mismatch")
       END IF
       this%offset = this%offset + 4
       CALL MPI_File_set_view(this%handle,this%offset,&
            this%basictype,this%filetype, 'native', MPI_INFO_NULL, this%error)
       !*****************************************************************!
       ! This collective call doesn't work for pvfs2 -> bug in ROMIO ?
!!$       CALL MPI_File_write_all(this%handle,this%binout,this%bufsize,&  
!!$            this%basictype, this%status, this%error)
       !*****************************************************************!
       CALL MPI_File_iread(this%handle,this%binout,this%bufsize,this%basictype,&
            request,this%error)
       CALL MPI_Wait(request,this%status,this%error)
#else
       READ (this%unit,IOSTAT=this%error) this%binout(:,:,:)
#endif
    END IF
    IF (this%error.NE.0) CALL Error(this,"ReadDataset_binary",&
         "reading data set failed")
    ! copy data
    k0 = this%cols-Physics%vnum
    FORALL (i=Mesh%IMIN:Mesh%IMAX,j=Mesh%JMIN:Mesh%JMAX,k=1:Physics%vnum) &
       Timedisc%pvar(i,j,k) = this%binout(k0+k,i,j)
  END SUBROUTINE ReadDataset_binary


  SUBROUTINE CloseFileIO_binary(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%header%idata,this%header%rdata,&
#ifdef PARALLEL
         this%disp,&
#endif
         this%binout)
  END SUBROUTINE CloseFileIO_binary

END MODULE fileio_binary
