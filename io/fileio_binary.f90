!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fileio_binary.f90                                                 #
!#                                                                           #
!# Copyright (C) 2008-2014                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling <sperling@astrophysik.uni-kiel.de>                         #
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
!> \author Tobias Illenseer
!! \author Björn Sperling
!!
!! \brief I/O for binary file format
!!
!! This module implements file I/O, which writes array data and a
!! short header with basic information to a file as binary data.
!! It is possible to select which data arrays should be written. 
!!
!! Specification: [header],[data],[bflux],[timestamp],[[data],[bflux],....]
!! - header : (4 + 10 * sizeof(INTEGER) + 10 * sizeof(REAL) + 4) bytes
!! - data : (4 + sizeof(REAL) * INUM * JNUM * (2+VNUM) + 4) bytes
!! - bflux : (4 + sizeof(REAL) * 4 * VNUM + 4) bytes
!! - timestamp: (4 + sizeof(REAL) + 4) bytes
!!
!! the leading and trailing 4 bytes are caused by the Fortran output
!!
!! \extends fileio_gnuplot
!! \ingroup fileio
!----------------------------------------------------------------------------!
MODULE fileio_binary
  USE fileio_gnuplot, CloseFile_binary => CloseFile_gnuplot
  USE geometry_common, ONLY : Geometry_TYP, GetType
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP, GetType
  USE timedisc_common, ONLY : Timedisc_TYP
  USE fluxes_generic, ONLY : Fluxes_TYP, GetBoundaryFlux
  USE common_dict
#ifdef PARALLEL
#ifdef HAVE_MPI_MOD
  USE mpi
#endif
#endif
  IMPLICIT NONE
#ifdef PARALLEL
#ifdef HAVE_MPIF_H
  include 'mpif.h'
#endif
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
  ! size of header data fields
  INTEGER, PARAMETER :: HISIZE  = 10      !< number of integer data in header    
  INTEGER, PARAMETER :: HRSIZE  = 10      !< number of real data in header
  INTEGER, PARAMETER :: MAXCOLS = 40      !< max number of data fields
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
  !> \public Constructor for the binary file I/O 
  !!
  !! Initilizes the file I/O type, filename, stoptime, number of outputs, 
  !! number of files, unit number, config as a dict
  SUBROUTINE InitFileIO_binary(this,Mesh,Physics,IO,fmt,fpath,filename,stoptime,dtwall,&
       count,fcycles,unit)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this          !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh          !< \param [in] Mesh mesh type
    TYPE(Physics_TYP) :: Physics       !< \param [in] Physics Physics type
    TYPE(Dict_TYP),POINTER :: IO       !< \param [in] IO Dictionary for I/O 
    INTEGER           :: fmt           !< \param [in] fmt fileio type number
    CHARACTER(LEN=*)  :: fpath         !< \param [in] fpath
    CHARACTER(LEN=*)  :: filename      !< \param [in] filename
    REAL              :: stoptime      !< \param [in] stoptime
    INTEGER           :: dtwall        !< \param [in] dtwall wall clock time
    INTEGER           :: count         !< \param [in] count number of outputs
    INTEGER           :: fcycles       !< \param [in] fcycles file cycle number
    INTEGER, OPTIONAL :: unit          !< \param [in] unit fileio unit number
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: node
    REAL,DIMENSION(:,:),POINTER :: dummy2
    REAL,DIMENSION(:,:,:),POINTER :: dummy3
    INTEGER           :: err
#ifdef PARALLEL
    INTEGER, DIMENSION(2) :: gsizes,lsizes,indices
    INTEGER           :: lb,extent
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,fmt,fpath,filename,stoptime,count,fcycles,unit
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL InitFileIO(this,Mesh,Physics,fmt,"binary",fpath,filename,"bin",stoptime, &
         dtwall,count,fcycles,.FALSE.,unit)

    ALLOCATE(this%output(MAXCOLS),STAT=err)
    IF (this%error.NE.0) THEN
       CALL Error(this,"InitFileIO_binary","Unable to allocate memory.")
    END IF

    CALL GetAttr(IO,"/mesh/bary_curv/value",dummy3)
    ! pointer to sub-array: set lower bounds
    ! use special function for bound remapping
    IF (Mesh%INUM.EQ.1) THEN
      ! use format string as temp
      this%output(1)%val => remap_bounds2(Mesh%IGMIN,Mesh%JGMIN,dummy3(:,:,2))
      this%cols = 1
    ELSE IF (Mesh%JNUM.EQ.1) THEN 
      this%output(1)%val => remap_bounds2(Mesh%IGMIN,Mesh%JGMIN,dummy3(:,:,1))
      this%cols = 1
    ELSE
      this%output(1)%val => remap_bounds2(Mesh%IGMIN,Mesh%JGMIN,dummy3(:,:,1))
      this%output(2)%val => remap_bounds2(Mesh%IGMIN,Mesh%JGMIN,dummy3(:,:,2))
      this%cols = 2
    END IF

    ! set output-pointer
    node => IO
    CALL GetOutputPointer(this,Mesh,node,this%cols)

    ! allocate memory
    ALLOCATE(this%header%idata(HISIZE),this%header%rdata(HRSIZE), &
         ! for output buffer
         this%binout(1:this%cols,Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX), &
         this%bflux(Physics%VNUM,4), &
#ifdef PARALLEL
         ! and displacements
         this%disp(Mesh%IMIN:Mesh%IMAX), &
#endif
         STAT=err)
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
    this%header%idata(5) = this%cols     
    this%header%idata(6:HISIZE) = 0      ! insert additional parameters here !
    ! 2. real data:
    this%header%rdata(1) = Mesh%xmin
    this%header%rdata(2) = Mesh%xmax
    this%header%rdata(3) = Mesh%ymin
    this%header%rdata(4) = Mesh%ymax
    this%header%rdata(5:HRSIZE) = 0.0    ! insert additional parameters here !
  END SUBROUTINE InitFileIO_binary

  !> Sets pointer to rank 2 array with a given lower bound
  !!
  !! This function sets a pointer to a rank 2 array with new
  !! lower boundaries. This is necessary because the lower bound of a array 
  !! will be set to zero when it is saved in a dictionary. 
  !! \todo reuse remap_bounds2 function from gnuplot 
  !! \return pointer to data array with new bounds
  FUNCTION remap_bounds2(lb1,lb2,array) RESULT(ptr)
    INTEGER, INTENT(IN) :: lb1      !< \param [in] lb1 new lower bound
    INTEGER, INTENT(IN) :: lb2      !< \param [in] lb2 new lower bound
    !> \param [in] array data array
    REAL, DIMENSION(lb1:,lb2:), INTENT(IN), TARGET :: array
  !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:), POINTER                  :: ptr
  !------------------------------------------------------------------------!
    ptr => array 
  END FUNCTION

  !> Creates a list of all data arrays which will be written to file
  !!
  !! Therefore it ignores all arrays of coordinates and checks if the data 
  !! arrays are of the same dimension as the mesh.
  !! \todo reuse GetOutputPointer subroutine from gnuplot 
  RECURSIVE SUBROUTINE GetOutputPointer(this,Mesh,node,k)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)       :: this  !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)         :: Mesh  !< \param [in] mesh mesh type
    TYPE(Dict_TYP),POINTER :: node  !< \param [in,out] node pointer to (sub-)dict
    INTEGER                :: k     !< \param [in,out] k number of data arrays
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: dir
    REAL,DIMENSION(:,:),POINTER :: dummy2
    REAL,DIMENSION(:,:,:),POINTER :: dummy3
    REAL,DIMENSION(:,:,:,:),POINTER :: dummy4
    INTEGER                :: dim3,dim4,i
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    INTENT(INOUT)     :: this,k
    !------------------------------------------------------------------------!
    DO WHILE(ASSOCIATED(node))
      IF(GetDataType(node) .EQ. DICT_DIR .AND. GetKey(node) .NE. "mesh") THEN
      ! recursion
        CALL GetAttr(node,GetKey(node),dir)
        CALL GetOutputPointer(this,Mesh,dir,k)
      ELSE IF (GetKey(node) .EQ. "value") THEN
      ! value found!
        SELECT CASE(GetDataType(node))
        CASE(DICT_REAL_TWOD)
          CALL GetAttr(node,GetKey(node),dummy2)
          IF (LBOUND(dummy2,DIM=1) .EQ. Mesh%IGMIN .AND.&
             UBOUND(dummy2,DIM=1) .EQ. Mesh%IGMAX .AND.&
             LBOUND(dummy2,DIM=2) .EQ. Mesh%JGMIN .AND.&
             UBOUND(dummy2,DIM=2) .EQ. Mesh%JGMAX) THEN
            k = k+1
            IF (k .GT. MAXCOLS) &
              CALL Error(this,"GetOutputPointer_binary","reached MAXCOLS")
            this%output(k)%val=> dummy2
          END IF
        CASE(DICT_REAL_THREED)
          CALL GetAttr(node,GetKey(node),dummy3)
          IF (LBOUND(dummy3,DIM=1) .EQ. Mesh%IGMIN .AND.&
             UBOUND(dummy3,DIM=1) .EQ. Mesh%IGMAX .AND.&
             LBOUND(dummy3,DIM=2) .EQ. Mesh%JGMIN .AND.&
             UBOUND(dummy3,DIM=2) .EQ. Mesh%JGMAX) THEN
            dim3 = SIZE(dummy3, DIM = 3)
            IF (k+dim3 .GT. MAXCOLS)&
              CALL Error(this,"GetOutputPointer_binary","reached MAXCOLS")
            DO i=k+1, k+dim3
              this%output(i)%val => remap_bounds2(Mesh%IGMIN,Mesh%JGMIN,dummy3(:,:,i-k))
            END DO
            k = k+dim3  
          END IF     
        CASE(DICT_REAL_FOURD)
          CALL GetAttr(node,GetKey(node),dummy4)
          IF (LBOUND(dummy4,DIM=1) .EQ. Mesh%IGMIN .AND.&
             UBOUND(dummy4,DIM=1) .EQ. Mesh%IGMAX .AND.&
             LBOUND(dummy4,DIM=2) .EQ. Mesh%JGMIN .AND.&
             UBOUND(dummy4,DIM=2) .EQ. Mesh%JGMAX) THEN
            dim3 = SIZE(dummy4, DIM = 3)
            dim4 = SIZE(dummy4, DIM = 4)
            IF (k+dim3*dim4 .GT. MAXCOLS)&
              CALL Error(this,"GetOutputPointer_binary","reached MAXCOLS")
            DO i=k+1, k+dim3*dim4
              this%output(k)%val => remap_bounds2(Mesh%IGMIN,Mesh%JGMIN,&
                  dummy4(:,:,(i-k-1)/dim4+1,mod((i-k-1),dim4)+1))
            END DO
            k = dim3*dim4
          END IF
        CASE DEFAULT
          !do nothing (wrong type)
        END SELECT
      END IF
      node=>GetNext(node)
    END DO
  END SUBROUTINE GetOutputPointer

  !> \public Specific routine to open a file for binary I/O
  !!
  SUBROUTINE OpenFile_binary(this,action)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this    !< \param [in,out] this fileio type
    INTEGER          :: action  !< \param [in] action mode of file access
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
    INTEGER           :: j,k
    INTEGER(KIND=4)   :: size1,size2  ! force 4 byte integer
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
#endif
    INTEGER           :: i,n,err
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

  !> \public Writes a short file format description as a header
  !!
  SUBROUTINE WriteHeader_binary(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this  !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    INTEGER(KIND=4)  :: size  !< force 4 byte integer
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


  !> \public Reads header from file and checks if it is consistent with simulation
  !!
  SUBROUTINE ReadHeader_binary(this,success)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this     !< \param [in,out] this fileio type
    LOGICAL          :: success  !< \param [out] success consistent
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER(KIND=4)  :: size1, size2    ! force 4 byte integer
#endif
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
       IF (idata(5).NE.this%header%idata(5)) THEN
          CALL Warning(this,"ReadHeader_binary","number of data fields mismatch")
          success = .FALSE.
       END IF
       IF (ALL(rdata(1:4).NE.this%header%rdata(1:4))) THEN
          CALL Warning(this,"ReadHeader_binary","computational domain mismatch")
          success = .FALSE.
       END IF
    END IF
  END SUBROUTINE ReadHeader_binary

  !> \public Writes timestamp to file
  !!
  SUBROUTINE WriteTimestamp_binary(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this  !< \param [in,out] this fileio type
    REAL             :: time  !< \param [in] time timestamp
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


  !> \public Reads timestamp to file
  !!
  SUBROUTINE ReadTimestamp_binary(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this  !< \param [in,out] this fileio type
    REAL             :: time  !< \param [out] time timestamp
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

  !> \public Writes all desired data arrays to a file 
  !!
  SUBROUTINE WriteDataset_binary(this,Mesh,Physics,Fluxes,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this      !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh      !< \param [in] mesh mesh type
    TYPE(Physics_TYP) :: Physics   !< \param [in] physics physics type
    TYPE(Fluxes_TYP)  :: Fluxes    !< \param [in] fluxes fluxes type
    TYPE(Timedisc_TYP):: Timedisc  !< \param [in] timedisc timedisc type
    !------------------------------------------------------------------------!
    INTEGER          :: i,j,k
#ifdef PARALLEL
    INTEGER(KIND=4)  :: size  ! force 4 byte integer
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
    INTEGER          :: request
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)       :: Mesh,Physics,Fluxes,Timedisc
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    ! prepare data for output
    DO k=1,this%cols
      DO j=Mesh%JMIN,Mesh%JMAX
        DO i=Mesh%IMIN,Mesh%IMAX
          ! copy data
          this%binout(k,i,j) = this%output(k)%val(i,j)
        END DO
      END DO
    END DO

    ! write data
#ifdef PARALLEL
    CALL MPI_File_get_size(this%handle,this%offset,this%error)
    CALL MPI_Barrier(MPI_COMM_WORLD,this%error)
    ! compute the number of bytes in real numbers we are going to
    ! write with _all_ processes in total 
    size = (Mesh%INUM*Mesh%JNUM*this%cols) * this%realext
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

    ! get boundary fluxes across all 4 boundaries
    DO i=1,4
       ! in PARALLEL mode the result is sent to process with rank 0
       this%bflux(:,i) = GetBoundaryFlux(Fluxes,Mesh,Physics,i)
    END DO
    ! write boundary fluxes
#ifdef PARALLEL
    ! compute size of boundary flux data array in bytes
    i = 4 * Physics%VNUM          ! number of real numbers
    size = i * this%realext       ! number of bytes
    ! skip 4 bytes of the size information of the main data set
    this%offset = this%offset + 4
    ! write size information for compatiblity with standard Fortran I/O
    IF (GetRank(this).EQ.0) THEN
       CALL MPI_File_write(this%handle,size,1,MPI_INTEGER4, &
            this%status,this%error)
    END IF
    this%offset = this%offset + 4
    CALL MPI_File_set_view(this%handle,this%offset,DEFAULT_MPI_REAL,&
         DEFAULT_MPI_REAL, 'native', MPI_INFO_NULL, this%error)
    ! write boundary fluxes
    IF (GetRank(this).EQ.0) THEN
       CALL MPI_File_write(this%handle,this%bflux,i,DEFAULT_MPI_REAL, &
            this%status,this%error)
    END IF
    this%offset = this%offset + size
    CALL MPI_File_set_view(this%handle,this%offset,MPI_INTEGER4,&
         MPI_INTEGER4, 'native', MPI_INFO_NULL, this%error)
    ! write size information again at the end
    IF (GetRank(this).EQ.0) THEN
       CALL MPI_File_write(this%handle,size,1,MPI_INTEGER4, &
            this%status,this%error)
    END IF
#else
    WRITE (this%unit) this%bflux(:,:)    
#endif
  END SUBROUTINE WriteDataset_binary

  !> \public Reads all data arrays from file 
  !!
  SUBROUTINE ReadDataset_binary(this,Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this      !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh      !< \param [in] mesh mesh type
    TYPE(Physics_TYP) :: Physics   !< \param [in] physics physics type
    TYPE(Timedisc_TYP):: Timedisc  !< \param [in] timedisc timedisc type
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
    ! point on data record (skip time stamp and boundary fluxes)
    CALL SeekBackwards(this,3,this%error)
    ! ***********************************************************************!
    ! FIXME: read boundary fluxes not implemented yet
    ! ***********************************************************************!    
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
       ! This collective call doesn't work for older pvfs2 versions 
       ! (maybe a bug in ROMIO) but the performance is better.
       CALL MPI_File_write_all(this%handle,this%binout,this%bufsize,&  
            this%basictype, this%status, this%error)
       !*****************************************************************!
       ! use non-collective calls only if collective call doesn't work
!       CALL MPI_File_iread(this%handle,this%binout,this%bufsize,this%basictype,&
!            request,this%error)
!       CALL MPI_Wait(request,this%status,this%error)
#else
       READ (this%unit,IOSTAT=this%error) this%binout(:,:,:)
#endif
    END IF
    IF (this%error.NE.0) CALL Error(this,"ReadDataset_binary",&
         "reading data set failed")

    ! copy data
    k0 = 0
    IF (Mesh%INUM .GT. 1) k0= 1
    IF (Mesh%JNUM .GT. 1) k0= k0+1

    DO k=k0+1,this%cols
      DO j=Mesh%JMIN,Mesh%JMAX
        DO i=Mesh%IMIN,Mesh%IMAX
          ! copy data
          this%output(k)%val(i,j) = this%binout(k,i,j)
        END DO
      END DO
    END DO
  END SUBROUTINE ReadDataset_binary

  !> \public Closes the file I/O
  !!
  SUBROUTINE CloseFileIO_binary(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this   !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    DEALLOCATE(this%header%idata,this%header%rdata,&
#ifdef PARALLEL
         this%disp,&
#endif
         this%binout,this%bflux)
  END SUBROUTINE CloseFileIO_binary

END MODULE fileio_binary
