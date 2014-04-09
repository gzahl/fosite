!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fileio_common.f90                                                 #
!#                                                                           #
!# Copyright (C) 2006-2014                                                   #
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
!> \defgroup fileio file I/O
!! \{
!! \brief Family of file I/O modules
!!
!! This is the family of file I/O modules. The generic interface routines are
!! defined in the module \link fileio_generic \endlink. The basic file I/O
!! data type and common basic subroutines and functions are defined in
!! \link fileio_common \endlink. Any other module of the family supports
!! a different file format.
!! \}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!!
!! \brief basic module for file I/O
!!
!! \extends common_types
!! \ingroup fileio
!----------------------------------------------------------------------------!
MODULE fileio_common
  USE common_types, &
       GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Initialized_common => Initialized, Info_common => Info, &
       Warning_common => Warning, Error_common => Error
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  USE common_dict
#ifdef HAVE_HDF5_MOD
  USE hdf5
#endif
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
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
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
     MODULE PROCEDURE FileIOError, Error_common
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  ! Private Attributes section starts here:
  !> \name Private Attributes
  !>#### file name and extension lengths
  INTEGER, PARAMETER :: FEXTLEN  = 4               !< file name extension length
  INTEGER, PARAMETER :: FNAMLEN  = 128-(FEXTLEN+1) !< file name length
  !> \name
  !!#### handling multiple files in parallel mode
  INTEGER, PARAMETER :: FMEXTLEN = 3               !< file extension length
  CHARACTER(LEN=FMEXTLEN+2), SAVE :: fmextstr = "" !< multi extension string
  !> \name
  !!#### handling multiple data files with time step in their names
  INTEGER, PARAMETER      :: MAXCYCLES = 10000     !< max. number of data files
  INTEGER, PARAMETER      :: FCYCLEN   = 4         !< num. of digits in file names
  CHARACTER(LEN=32), SAVE :: cycfmt                !< format string for cycles
  !--------------------------------------------------------------------------!
  !> data type for file header data
  TYPE Header_TYP
     INTEGER, DIMENSION(:), POINTER :: idata
     REAL, DIMENSION(:), POINTER    :: rdata
  END TYPE Header_TYP
  !> output-pointer (binary,gnuplot,vtk)
  TYPE Output_TYP
    REAL, DIMENSION(:,:),POINTER :: val
    CHARACTER(LEN=128)           :: key
  END TYPE Output_TYP
  !> basic data type for file i/o
  TYPE FileIO_TYP
     !> \name Variables
     TYPE(Common_TYP)       :: format      !< i/o file format
     CHARACTER(LEN=512)     :: linebuf     !< buffer for character i/o
     CHARACTER(LEN=FNAMLEN) :: filename    !< file name without extension
     CHARACTER(LEN=FNAMLEN) :: path        !< file path without filename 
     CHARACTER(LEN=FEXTLEN) :: extension   !< file name extension
     INTEGER                :: unit        !< i/o unit
     INTEGER                :: error       !< i/o error code
     INTEGER                :: step        !< counter for output steps
     INTEGER                :: count       !< number of output steps
     INTEGER                :: cycles      !< number of output files
     INTEGER                :: dtwall      !< wall clock time difference
     INTEGER                :: realsize    !< byte size of real numbers
     REAL                   :: stoptime    !< final simulation time for data output
     REAL                   :: time        !< output time
     !> \name
     !!#### VTK specific variables
     CHARACTER(LEN=512)     :: linebuf2    !< buffer for character i/o
     CHARACTER(LEN=32)      :: buf         !< buffer for character i/o
     CHARACTER(LEN=12)      :: realfmt     !< real format string
     CHARACTER(LEN=14)      :: endianness  !< endianness string
     INTEGER                :: ioffset     !< appended data offset
     !> \name
     !!#### GNUPLOT and BINARY specific variables
     TYPE(Header_TYP)       :: header      !< \public binary file header
     CHARACTER(LEN=64)      :: fmtstr      !< format string
     CHARACTER(LEN=64)      :: linefmt     !< output line format string
     INTEGER                :: COLS        !< number of output columns
     INTEGER                :: MAXCOLS     !< upper limit for output cols
                                           !! (MAXCOLS < LEN(linebuf)/FLEN)
     INTEGER                :: DECS        !< decimal places for real number output
     INTEGER                :: FLEN        !< output field length
     INTEGER                :: linelen     !< length of a line
#ifdef HAVE_NETCDF
     !> \name
     !!#### netCDF specific variables
     INTEGER                :: ncid        !< file id
     INTEGER                :: ncfmt       !< file format
     INTEGER                :: rank        !< 1D or 2D data
     !< \public
#endif
#ifdef HAVE_HDF5_MOD
     !> \name
     !!#### HDF specific variables
     INTEGER(HID_T)         :: fid         !< file id
     INTEGER(HID_T)         :: xferid      !< xfer id
#endif
     TYPE(Output_TYP),DIMENSION(:), POINTER :: &
                               output      !< list of output fields
     REAL, DIMENSION(:,:) , POINTER  :: &
                               bflux       !< boundary flux output buffer
     REAL, DIMENSION(:,:,:), POINTER :: &
                               vtktemp, &  !< VTK temporary data 
                               vtktemp2, & !< VTK temporary data
                               binout      !< binary data output buffer
#ifdef PARALLEL
     !> \name Variables in Parallel Mode
     LOGICAL                :: multfiles   !< spread files across nodes
     INTEGER                :: handle      !< MPI file handle
     INTEGER                :: bufsize     !< output data buffer size
     INTEGER                :: basictype   !< data type for points
     INTEGER                :: filetype    !< data type for data i/o
     INTEGER, DIMENSION(MPI_STATUS_SIZE) :: &
                               status      !< MPI i/o status record
     !> \name
     !!#### GNUPLOT and BINARY specific variables
     INTEGER                :: blocknum    !< number of output blocks
     INTEGER(KIND=MPI_OFFSET_KIND) :: &
                               offset      !< skip header bytes
     INTEGER(KIND=MPI_ADDRESS_KIND) :: &
                               realext, &  !< real data type extent
                               intext      !< integer data type extent
     !> \name
     !!#### VTK specific variables
     CHARACTER(LEN=64), DIMENSION(:), POINTER :: &
                               extent      !< extent of all pieces in VTK
     CHARACTER, DIMENSION(:,:), POINTER :: &
                               outbuf      !< output buffer
     !!#### XDMF specific variables
     INTEGER                :: meshbufsize !< size of mesh output buffer
     INTEGER                :: meshtype    !< data type for mesh i/o
     INTEGER                :: memtype     !< data type for memory 
     INTEGER, DIMENSION(:), POINTER :: &
                               disp        !< array of displacements
#endif
  END TYPE FileIO_TYP
  !--------------------------------------------------------------------------!
  !> \name Public Attributes
  !! #### file status and access modes
  INTEGER, PARAMETER :: FILE_EXISTS = B'00000001' !< file status for existing files
  INTEGER, PARAMETER :: READONLY = 1       !< readonly access
  INTEGER, PARAMETER :: READEND  = 2       !< readonly access at end
  INTEGER, PARAMETER :: REPLACE  = 3       !< read/write access replacing file
  INTEGER, PARAMETER :: APPEND   = 4       !< read/write access at end
  !> \name
  !! #### file formats
  CHARACTER(LEN=9), PARAMETER :: ASCII = "formatted"   !< for ASCII data
  CHARACTER(LEN=11), PARAMETER :: BIN  = "unformatted" !< for BINARY data
  !> \}
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       FileIO_TYP, Header_TYP, Output_TYP, &
       ! constants
       READONLY, READEND, REPLACE, APPEND, &
       ASCII, BIN, &
       FILE_EXISTS, &
#ifdef PARALLEL
       DEFAULT_MPI_REAL,&
#endif
       ! methods
       InitFileIO, &
       CloseFileIO, &
       MakeMultstr, &
       GetBasename, &
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

  !> \public Constructor of common fileio class.
  !!
  !! Initializes the file format, format name, file name with extension,
  !! number of file cycles, number of files in parallel mode and the
  !! i/o unit number.
  SUBROUTINE InitFileIO(this,fmt,fmt_name,fpath,fname,fext,fcycles,&
      multfiles,unit)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this      !< \param [in,out] this fileio type
    INTEGER          :: fmt       !< \param [in] fmt file format
    CHARACTER(LEN=*) :: fmt_name  !< \param [in] fmt_name file format name
    CHARACTER(LEN=*) :: fpath     !< \param [in] fpath file path
    CHARACTER(LEN=*) :: fname     !< \param [in] fname file name
    CHARACTER(LEN=*) :: fext      !< \param [in] fext file name extension
    INTEGER          :: fcycles   !< \param [in] fcycles number of file cycles
    LOGICAL          :: multfiles !< \param [in] multfiles spread file in parallel i/o
    INTEGER          :: unit      !< \param [in] unit fortran i/o unit number
    !------------------------------------------------------------------------!
    INTENT(IN)       :: fmt,fmt_name,fpath,fname,fext,fcycles,multfiles,unit
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%format,fmt,fmt_name)
    ! check cycles
    IF (fcycles.GT.MAXCYCLES) THEN
       CALL Error(this,"InitFileIO","file cycles exceed limits")
    END IF
    ! check file name
    IF (fcycles.GT.0) THEN
       IF (LEN_TRIM(fpath)+LEN_TRIM(fname).GT.(FNAMLEN-(FCYCLEN+1))) &
          CALL Error(this,"InitFileIO","file name too long")
    ELSE IF (LEN_TRIM(fpath)+LEN_TRIM(fname).GT.FNAMLEN) THEN
       CALL Error(this,"InitFileIO","file name too long")
    END IF

    ! check file name extension
    IF (LEN_TRIM(fext).GT.FEXTLEN) THEN
       CALL Error(this,"InitFileIO","file name extension too long")
    END IF
    ! format string for writing file names with explicit time step
    WRITE (cycfmt, "('(A,I',I1,'.',I1,',A)')") FCYCLEN,FCYCLEN

#ifdef PARALLEL
    this%multfiles = multfiles
    IF (this%multfiles) THEN
       ! check file name
       IF (fcycles.GT.0) THEN
          IF (LEN_TRIM(fpath)+LEN_TRIM(fname).GT.(FNAMLEN-(FMEXTLEN+2+FCYCLEN+1))) &
             CALL Error(this,"InitFileIO","file name too long")
       ELSE IF (LEN_TRIM(fpath)+LEN_TRIM(fname).GT.FNAMLEN-(FMEXTLEN+2)) THEN
          CALL Error(this,"InitFileIO","file name too long")
       END IF
       fmextstr = MakeMultstr(this)
    END IF
#endif

    this%path     = fpath
    this%filename = fname
    this%extension= fext
    this%cycles   = fcycles
    this%error    = 0
    this%unit     = unit
  END SUBROUTINE InitFileIO


  !> \public Destructor of common fileio class.
  SUBROUTINE CloseFileIO(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(INOUT) :: this !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    CALL CloseCommon(this%format)
  END SUBROUTINE CloseFileIO

  !> \public Get a file label (multiples files in parallel mode)
  !! without filenumber => use GetRank; with filenumber < 0 => empty label
  !! \result file label
  FUNCTION MakeMultstr(this,fn) RESULT (multstr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this !< \param [in,out] this fileio type
    INTEGER, OPTIONAL, INTENT(IN):: fn       !< \param [in] fn number of file
    CHARACTER(LEN=FMEXTLEN+2)    :: multstr
    !------------------------------------------------------------------------!
    INTEGER            :: fn_l
    CHARACTER(LEN=32)  :: mextfmt
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    IF (PRESENT(fn)) THEN
      fn_l = fn
    ELSE
      fn_l = GetRank(this)
    END IF
    ! with fn < 0 you can suppress a label
    IF (this%multfiles .AND. fn_l .GE. 0) THEN
        WRITE (mextfmt, "('(A,I',I1,'.',I1,')')") FMEXTLEN,FMEXTLEN
        
        WRITE (multstr, mextfmt) "-r", fn_l
    ELSE 
      multstr = ""
    END IF
#else
    multstr = ""
#endif
   
  END FUNCTION MakeMultstr


  !> \public Get the current file name without path 
  !! e.g. important for vtk (pvts files)
  !! \result current file name
  FUNCTION GetBasename(this,fn) RESULT (fname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this !< \param [in] this fileio type
    INTEGER, OPTIONAL, INTENT(IN):: fn   !< \param [in] fn number of file 
    CHARACTER(LEN=256)           :: fname
    !------------------------------------------------------------------------!
    CHARACTER(LEN=FCYCLEN+2)  :: cycstr
    !------------------------------------------------------------------------!
    IF (this%cycles.GT.0) THEN
       ! generate a file name with time step
       WRITE (cycstr, FMT=TRIM(cycfmt)) "_", MODULO(this%step,this%cycles), "."
       WRITE (fname,"(A,A,A,A)")  TRIM(this%filename),&
              TRIM(MakeMultstr(this,fn)),TRIM(cycstr),TRIM(this%extension)
    ELSE
       ! file name + extension
       WRITE (fname,"(A,A,A,A)") TRIM(this%filename),&
              TRIM(MakeMultstr(this,fn)),".",TRIM(this%extension)
    END IF
  END FUNCTION GetBasename

  !> \public Get the current file name
  !! \result current file name
  FUNCTION GetFilename(this,fn) RESULT (fname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this !< \param [in] this fileio type
    INTEGER, OPTIONAL, INTENT(IN):: fn   !< \param [in] fn number of file 
    CHARACTER(LEN=256) :: fname
    !------------------------------------------------------------------------!
    fname = TRIM(this%path) // TRIM(GetBasename(this,fn))
  END FUNCTION GetFilename

  !> \public Get the file status
  !!
  !! Checks if the file exists.
  !! \result file status
  FUNCTION GetFilestatus(this) RESULT(fstatus)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this !< \param [in] this fileio type
    INTEGER :: fstatus
    !------------------------------------------------------------------------!
    LOGICAL :: success
    !------------------------------------------------------------------------!
    fstatus = 0
    INQUIRE (FILE=GetFilename(this),EXIST=success)
    IF (success) fstatus = IOR(fstatus,FILE_EXISTS)
  END FUNCTION GetFilestatus

 
  !> \public Get the file type; 
  !! overloads \b GetType from \link common_types::GetType \endlink
  !! \result file type
  PURE FUNCTION GetFileIOType(this) RESULT(fmt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this !< \param [in] this fileio type
    INTEGER :: fmt
    !------------------------------------------------------------------------!
    fmt = GetType_common(this%format)
  END FUNCTION GetFileIOType


  !> \public Get the file type name; 
  !! overloads \b GetName from \link common_types::GetName \endlink
  !! \result file type name 
  PURE FUNCTION GetFileIOName(this) RESULT(fn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this !< \param [in] this fileio type
    CHARACTER(LEN=32) :: fn
    !------------------------------------------------------------------------!
    fn = GetName_common(this%format)
  END FUNCTION GetFileIOName


  !> \public Get the MPI rank;
  !! overloads \b GetRank from \link common_types::GetRank \endlink
  !! \return MPI rank
  PURE FUNCTION GetFileIORank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this !< \param [in] this fileio type
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%format)
  END FUNCTION GetFileIORank


  !> \public Get the total number of MPI processes;
  !! overloads \b GetNumProcs from \link common_types::GetNumProcs \endlink
  !! \return number of MPI processes
  PURE FUNCTION GetFileIONumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this !< \param [in] this fileio type
    INTEGER :: p
    !------------------------------------------------------------------------!
    p = GetNumProcs_common(this%format)
  END FUNCTION GetFileIONumProcs


  !> \public Query initialization status;
  !! overloads \b Initialized from \link common_types::Initialized \endlink
   PURE FUNCTION FileIOInitialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this !< \param [in] this fileio type
    LOGICAL :: i
    !------------------------------------------------------------------------!
    i = Initialized_common(this%format)
  END FUNCTION FileIOInitialized


  !> \public Print information on standard output;
  !! overloads \b Info from \link common_types::Info \endlink
  SUBROUTINE FileIOInfo(this,msg,rank,node_info,tostderr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this !< \param [in] this fileio type
    CHARACTER(LEN=*),  INTENT(IN) :: msg  !< \param [in] msg info message
    INTEGER, OPTIONAL, INTENT(IN) :: rank !< \param [in] rank MPI rank
    LOGICAL, OPTIONAL, INTENT(IN) :: node_info !< \param [in] node_info enable rank output
    LOGICAL, OPTIONAL, INTENT(IN) :: tostderr  !< \param [in] tostderr enable STDERR output
    !------------------------------------------------------------------------!
    CALL Info_common(this%format,msg,rank,node_info,tostderr)
  END SUBROUTINE FileIOInfo


  !> \public Print warning message on standard error;
  !! overloads \b Warning from \link common_types::Warning \endlink
  SUBROUTINE FileIOWarning(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this !< \param [in] this fileio type
    CHARACTER(LEN=*),  INTENT(IN) :: modproc !< \param [in] modproc name of module procedure
    CHARACTER(LEN=*),  INTENT(IN) :: msg !< \param [in] msg warning message
    !------------------------------------------------------------------------!
    CALL Warning_common(this%format,modproc,msg)
  END SUBROUTINE FileIOWarning


  !> \public Print error message on standard error and terminate the program;
  !! overloads \b Error from \link common_types::Error \endlink
  SUBROUTINE FileIOError(this,modproc,msg,rank,node_info)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(IN) :: this !< \param [in] this fileio type
    CHARACTER(LEN=*),  INTENT(IN) :: modproc !< \param [in] modproc name of module procedure
    CHARACTER(LEN=*),  INTENT(IN) :: msg !< \param [in] msg warning message
    INTEGER, OPTIONAL, INTENT(IN) :: rank !< \param [in] rank MPI rank
    LOGICAL, OPTIONAL, INTENT(IN) :: node_info !< \param [in] node_info enable rank output
    !------------------------------------------------------------------------!
    CALL Error_common(this%format,modproc,msg,rank,node_info)
  END SUBROUTINE FileIOError


END MODULE fileio_common
