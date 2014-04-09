!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fileio_common.f90                                                 #
!#                                                                           #
!# Copyright (C) 2006-2010                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
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
! basic module for file I/O
!----------------------------------------------------------------------------!
MODULE fileio_common
  USE common_types, &
       GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Initialized_common => Initialized, Info_common => Info, &
       Warning_common => Warning, Error_common => Error
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
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
  INTERFACE GetType
     MODULE PROCEDURE GetFileIOType, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetFileIOName, GetName_common
  END INTERFACE
  INTERFACE GetRank
     MODULE PROCEDURE GetFileIORank, GetRank_common
  END INTERFACE
  INTERFACE GetNumProcs
     MODULE PROCEDURE GetFileIONumProcs, GetNumProcs_common
  END INTERFACE
  INTERFACE Initialized
     MODULE PROCEDURE FileIOInitialized, Initialized_common
  END INTERFACE
  INTERFACE Info
     MODULE PROCEDURE FileIOInfo, Info_common
  END INTERFACE
  INTERFACE Warning
     MODULE PROCEDURE FileIOWarning, Warning_common
  END INTERFACE
  INTERFACE Error
     MODULE PROCEDURE FileIOError_rank0, FileIOError_rankX, Error_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  ! length of file names and extensions
  INTEGER, PARAMETER :: FEXTLEN  = 4               ! extension length        !
  INTEGER, PARAMETER :: FNAMLEN  = 256-(FEXTLEN+1) ! file name length        !
  ! handling multiple data files with time step in their names
  INTEGER, PARAMETER      :: MAXCYCLES = 10000 ! max. number of data files   !
  INTEGER, PARAMETER      :: FCYCLEN   = 4     ! num. of digits in file names!
  CHARACTER(LEN=32), SAVE :: cycfmt            ! format string for cycles    !
  !--------------------------------------------------------------------------!
  ! data type for file header data
  TYPE Header_TYP
     INTEGER, DIMENSION(:), POINTER :: idata
     REAL, DIMENSION(:), POINTER    :: rdata
  END TYPE Header_TYP
  ! basic data type for file i/o
  TYPE FileIO_TYP
     TYPE(Common_TYP)       :: format              ! i/o file format         !
     TYPE(Header_TYP)       :: header              ! file header             !
     CHARACTER(LEN=256)     :: linebuf             ! buffer for character i/o!
     CHARACTER(LEN=FNAMLEN) :: filename            ! file name w/o extension !
     CHARACTER(LEN=FEXTLEN) :: extension           ! file name extension     !
     CHARACTER(LEN=64)      :: fmtstr              ! format string           !
     CHARACTER(LEN=64)      :: linefmt             ! output line format str. !
     CHARACTER(LEN=12)      :: realfmt             ! real format str. for vtk!
     CHARACTER(LEN=14)      :: endianness          ! endianness str. for vtk !
     INTEGER                :: cols                ! no. of output columns   !
     INTEGER                :: linelen             ! length of a line        !
     INTEGER                :: error               ! i/o error code          !
     INTEGER                :: step                ! counter for output steps!
     INTEGER                :: count               ! number of output steps  !
     INTEGER                :: cycles              ! number of output files  !
     INTEGER                :: dtwall              ! wall clock time diff.   !
     INTEGER                :: ioffset             ! VTK: appended data offset!
     REAL                   :: stoptime            ! end of data output      !
     REAL                   :: time                ! output time             !
     REAL, DIMENSION(:,:,:), POINTER :: &
                               vtktemp,vtktemp2    ! vtk temp. data          !
     REAL, DIMENSION(:,:,:), POINTER :: binout     ! binary data output buf. !
     REAL, DIMENSION(:,:) , POINTER  :: bflux      ! bound. flux data output !               
#ifdef HAVE_NETCDF
     INTEGER                :: ncid                ! file id                 !
     INTEGER                :: ncfmt               ! file format             !
     INTEGER                :: rank                ! 1D or 2D data           !
#endif
#ifdef PARALLEL
     CHARACTER, DIMENSION(:,:), POINTER :: outbuf  ! output buffer           !
     INTEGER                :: bufsize             ! size of output buffer   !
     INTEGER                :: handle              ! MPI file handle         !
     INTEGER                :: blocknum            ! no. of blocks in data f.!
     INTEGER                :: basictype           ! data type for points    !
     INTEGER                :: filetype            ! data type for i/o       !
     INTEGER(KIND=MPI_ADDRESS_KIND) :: realext,intext ! data type extent     !
     INTEGER(KIND=MPI_OFFSET_KIND) :: offset       ! skip header bytes       !
     INTEGER, DIMENSION(:), POINTER :: disp        ! array of displacements  !
     INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status ! MPI i/o status record   !
     LOGICAL                :: sepfiles            ! one or multiple files in parallel mode!
#endif
     INTEGER                :: unit                ! i/o unit                !
  END TYPE FileIO_TYP
  !--------------------------------------------------------------------------!
  INTEGER, PARAMETER :: READONLY = 1               ! file access modes       !
  INTEGER, PARAMETER :: READEND  = 2
  INTEGER, PARAMETER :: REPLACE  = 3
  INTEGER, PARAMETER :: APPEND   = 4
  ! file status
  INTEGER, PARAMETER :: FILE_EXISTS = B'00000001'
  ! file formats
  CHARACTER(LEN=9), PARAMETER :: ASCII  = "formatted"
  CHARACTER(LEN=11), PARAMETER :: BIN    = "unformatted"
  !--------------------------------------------------------------------------!
  INTEGER, SAVE :: lastunit = 10
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       FileIO_TYP, Header_TYP, &
       ! constants
       READONLY, READEND, REPLACE, APPEND, &
       ASCII, BIN, &
       FILE_EXISTS, &
#ifdef PARALLEL
       DEFAULT_MPI_REAL,&
#endif
       ! methods
       InitFileIO, &
       GetFilename, &
       GetFilestatus, &
       GetType, &
       GetName, &
       GetRank, &
       GetNumProcs, &
       Initialized, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitFileIO(this,fmt,fmt_name,fname,fext,fcycles,sepfiles,unit)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    INTEGER          :: fmt
    CHARACTER(LEN=*) :: fmt_name
    CHARACTER(LEN=*) :: fname
    CHARACTER(LEN=*) :: fext
    INTEGER          :: fcycles
    LOGICAL          :: sepfiles
    INTEGER, OPTIONAL:: unit
    !------------------------------------------------------------------------!
    CHARACTER(LEN=FCYCLEN+2)  :: rankstr  ! same length as fcycle 
    !------------------------------------------------------------------------!
    INTENT(IN)       :: fmt,fmt_name,fname,fext,fcycles,sepfiles,unit
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%format,fmt,fmt_name)
    ! check cycles
    IF (fcycles.GT.MAXCYCLES) THEN
       CALL Error(this,"InitFileIO","file cycles exceed limits")
    END IF
    ! check file name
    IF (fcycles.GT.0) THEN
       IF (LEN_TRIM(fname).GT.(FNAMLEN-(FCYCLEN+1))) &
          CALL Error(this,"InitFileIO","file name too long")
    ELSE IF (LEN_TRIM(fname).GT.FNAMLEN) THEN
       CALL Error(this,"InitFileIO","file name too long")
    END IF
    ! check file name extension
    IF (LEN_TRIM(fext).GT.FEXTLEN) THEN
       CALL Error(this,"InitFileIO","file name extension too long")
    END IF
    ! format string for writing file names with explicit time step
    WRITE (cycfmt, "('(A,I',I1,'.',I1,',A)')") FCYCLEN,FCYCLEN

    this%filename = fname
#ifdef PARALLEL
    this%sepfiles = sepfiles
    IF (this%sepfiles) THEN
       ! check file name
       IF (fcycles.GT.0) THEN
          IF (LEN_TRIM(fname).GT.(FNAMLEN-(2*FCYCLEN+3))) &
             CALL Error(this,"InitFileIO","file name too long")
       ELSE IF (LEN_TRIM(fname).GT.FNAMLEN-(FCYCLEN+2)) THEN
          CALL Error(this,"InitFileIO","file name too long")
       END IF
       ! sets rankstr
       WRITE (rankstr, FMT=TRIM(cycfmt))"-r",GetRank(this)
       this%filename = fname//rankstr
    END IF
#endif
    this%extension= fext
    this%cycles   = fcycles
    this%error    = 0

    IF (PRESENT(unit)) THEN
       this%unit = unit
    ELSE
       ! this ensures that no unit number is assigned twice
       this%unit = lastunit + 1
       lastunit  = this%unit
    END IF
  END SUBROUTINE InitFileIO


  FUNCTION GetFilename(this) RESULT (fname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this
    CHARACTER(LEN=256) :: fname
    !------------------------------------------------------------------------!
    CHARACTER(LEN=FCYCLEN+2)  :: cycstr
    !------------------------------------------------------------------------!
    IF (this%cycles.GT.0) THEN
       ! generate a file name with time step
       WRITE (cycstr, FMT=TRIM(cycfmt)) "_", MODULO(this%step,this%cycles), "."
       WRITE (fname,"(A,A,A)") TRIM(this%filename),TRIM(cycstr), TRIM(this%extension)
    ELSE
       ! file name + extension
       WRITE (fname,"(A)") TRIM(this%filename)// "." // TRIM(this%extension)
    END IF
  END FUNCTION GetFilename


  FUNCTION GetFilestatus(this) RESULT(fstatus)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this
    INTEGER :: fstatus
    !------------------------------------------------------------------------!
    LOGICAL :: success
    !------------------------------------------------------------------------!
    fstatus = 0
    INQUIRE (FILE=GetFilename(this),EXIST=success)
    IF (success) fstatus = IOR(fstatus,FILE_EXISTS)
  END FUNCTION GetFilestatus

 
  PURE FUNCTION GetFileIOType(this) RESULT(fmt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this
    INTEGER :: fmt
    !------------------------------------------------------------------------!
    fmt = GetType_common(this%format)
  END FUNCTION GetFileIOType


  PURE FUNCTION GetFileIOName(this) RESULT(fn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: fn
    !------------------------------------------------------------------------!
    fn = GetName_common(this%format)
  END FUNCTION GetFileIOName


  PURE FUNCTION GetFileIORank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%format)
  END FUNCTION GetFileIORank


  PURE FUNCTION GetFileIONumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this
    INTEGER :: p
    !------------------------------------------------------------------------!
    p = GetNumProcs_common(this%format)
  END FUNCTION GetFileIONumProcs


  PURE FUNCTION FileIOInitialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this
    LOGICAL :: i
    !------------------------------------------------------------------------!
    i = Initialized_common(this%format)
  END FUNCTION FileIOInitialized


  SUBROUTINE FileIOInfo(this,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: msg
    !------------------------------------------------------------------------!
    CALL Info_common(this%format,msg)
  END SUBROUTINE FileIOInfo


  SUBROUTINE FileIOWarning(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Warning_common(this%format,modproc,msg)
  END SUBROUTINE FileIOWarning


  SUBROUTINE FileIOError_rank0(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN)    :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Error_common(this%format,modproc,msg)
  END SUBROUTINE FileIOError_rank0


  SUBROUTINE FileIOError_rankX(this,modproc,msg,rank)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN)    :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    INTEGER, INTENT(IN)           :: rank
    !------------------------------------------------------------------------!
    CALL Error_common(this%format,modproc,msg,rank)
  END SUBROUTINE FileIOError_rankX


END MODULE fileio_common
