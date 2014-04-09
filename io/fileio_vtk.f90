!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fileio_vtk.f90                                                    #
!#                                                                           #
!# Copyright (C) 2010-2014                                                   #
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
!> \author Björn Sperling
!!
!! \brief I/O for VTK files as XML format (vtkStructuredGrid)
!!
!! This module implements VTK file I/O to write vts files (vtkStructuredGrid) 
!! in XML syntax. In each file is exactly one timestep. 
!! In case of parallel computation there is one file per job (and timestep) 
!! and only one global container file (pvts) which groups all vts files.
!! 
!! see http://www.vtk.org for more details
!!  
!!
!! \extends fileio_common
!! \ingroup fileio
!----------------------------------------------------------------------------!
#ifdef FORTRAN_STREAMS
#define HAVE_VTK
#elif NECSX9
#define HAVE_VTK
#define NOSTREAM
#elif NECSX8
#define HAVE_VTK
#define NOSTREAM
#endif
MODULE fileio_vtk
  USE fileio_common
  USE geometry_common, ONLY : Geometry_TYP, GetName, GetType
  USE geometry_generic
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP, GetName, GetType
  USE fluxes_generic, ONLY : Fluxes_TYP, GetBoundaryFlux
  USE timedisc_common, ONLY : Timedisc_TYP
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
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       FileIO_TYP, &
       ! constants

       ! methods
       InitFileio_vtk, &
       OpenFile_vtk, &
       CloseFile_vtk, &
       WriteHeader_vtk, &
       ReadHeader_vtk, &
       WriteTimestamp_vtk, &
       ReadTimestamp_vtk,&
       WriteDataset_vtk, &
       GetPrecision_vtk, &
       GetEndianness_vtk, &
       CloseFileio_vtk
  !--------------------------------------------------------------------------!
   INTEGER, PARAMETER      :: MAXCOMP  = 9   !< max. of allowed components 
                                             !! 9 is a tensor (rank 2, dim 3)
   INTEGER(KIND = 4)       :: N_Byte         !< blocksize (4 byte!!!)
   INTEGER                 :: indent         !< indent of current paragraph (xml)
   INTEGER, PARAMETER      :: MAXCOLS  = 40  !< max of different output arrays
   INTEGER, PARAMETER      :: MAXKEY   = 64  !< max length of keyname
   CHARACTER, PARAMETER    :: LF = ACHAR(10) !< line feed    
   !> names of fluxes
   CHARACTER(LEN=15),DIMENSION(4),PARAMETER  :: fluxkey = (/'bflux_WEST ', &
                                                            'bflux_EAST ', &
                                                            'bflux_SOUTH', &
                                                            'bflux_NORTH' /)
!----------------------------------------------------------------------------------------------------------------------------------

CONTAINS
  !> \public Constructor for the VTK file I/O 
  !!
  !! Initilizes the file I/O type, filename, stoptime, number of outputs, 
  !! number of files, unit number, config as a dict
  SUBROUTINE InitFileio_vtk(this,Mesh,Physics,IO,fmt,fpath,filename,stoptime,dtwall,&
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
    INTEGER           :: err,k,n,i
    REAL              :: ftime
    INTEGER, DIMENSION(:), POINTER :: sendbuf,recvbuf
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,fmt,fpath,filename,stoptime,count,fcycles,unit
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL InitFileIO(this,fmt,"VTK",fpath,filename,'vts',fcycles,.TRUE.,unit)
       this%stoptime = stoptime
       this%dtwall   = dtwall
       this%time     = 0.
       this%count    = count
       this%step     = 0
       this%ioffset  = 0
#ifndef HAVE_VTK
    CALL Error(this,"InitFileIO","VTK support disabled. Check config.log")
#endif
    IF (fcycles .NE. count+1) CALL Error(this,'InitFileIO','VTK need filecycles = count+1')

		CALL GetPrecision_vtk(this, this%realsize)

    ! save size of floats to string: realfmt
    write(this%linebuf,'(I4)',IOSTAT=this%error) 8*this%realsize
    write(this%realfmt,'(A,A,A)',IOSTAT=this%error)'"Float',trim(AdjustL(this%linebuf)),'"'

    IF (BIT_SIZE(N_BYTE) .NE. 4*8) &
          CALL Error(this, "InitFileio_vtk", "INTERGER(KIND=4)::N_Byte failed")

		CALL GetEndianness_vtk(this,this%endianness,'"LittleEndian"','"BigEndian"')

    ALLOCATE(this%vtktemp(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,3),&
             this%vtktemp2(MAXCOMP,Mesh%IMIN:Mesh%IMAX+1,Mesh%JMIN:Mesh%JMAX+1),&
             this%output(MAXCOLS),&
             this%bflux(Physics%VNUM,4),&
             STAT = err)
    IF (err.NE.0) &
         CALL Error(this, "InitFileio_vtk", "Unable to allocate memory.")

#ifdef PARALLEL
    ALLOCATE(sendbuf(4),STAT = err)
    IF (err.NE.0) &
      CALL Error(this, "InitFileio_vtk", "Unable to allocate memory.")

    ! allocate memory for receive buffers and this%extent;
    ! both arrays are only needed on the receiver node (rank 0)
    ! but parallel profiling using scalasca with mpich3 crashes
    ! if recvbuf is unallocated hence we allocate a minimal array
    ! of size 1 on all nodes except for the receiver node;
    ! remark: according to the MPI-2 standard recvbuf can be undefined
    IF (GetRank(this) .EQ. 0 ) THEN
      ALLOCATE(this%extent(0:GetNumProcs(this)-1),&
               recvbuf(4*GetNumProcs(this)),&
               STAT = err)
    ELSE
      ! see comment above
      ALLOCATE(recvbuf(1),&
               STAT = err)
    END IF
    IF (err.NE.0) &
         CALL Error(this, "InitFileio_vtk", "Unable to allocate memory.")

    sendbuf(1) = Mesh%IMIN
    sendbuf(2) = Mesh%IMAX+1 
    sendbuf(3) = Mesh%JMIN 
    sendbuf(4) = Mesh%JMAX+1  
 
    CALL MPI_Gather(sendbuf, 4, MPI_INTEGER, &
               recvbuf, 4, MPI_INTEGER, &
               0, MPI_COMM_WORLD, this%error)

    IF (GetRank(this) .EQ. 0 ) THEN
      DO i=0,GetNumProcs(this)-1
        write(this%extent(i),fmt='(6(I7))',IOSTAT=this%error)&
          recvbuf(i*4+1),recvbuf(i*4+2),recvbuf(i*4+3),recvbuf(i*4+4),1,1
      END DO
    END IF
    DEALLOCATE(sendbuf,recvbuf,STAT=err)
#endif

    ! set mesh
    this%output(1)%key ="/mesh/corners/value"
    this%cols = 1
    ! set outputlist
    node => IO
    CALL GetOutputlist(this,Mesh,node,this%cols)
 
!    ! write a pvd-file for paraview, default: enabled
!    CALL RequireKey(config, "ParaviewFile", 1)
!    CALL GetAttr(config, "ParaviewFile", pvd)
    
!    IF ((pvd .EQ. 1) .AND. (GetRank(this).EQ.0)) &
     IF (GetRank(this).EQ.0) &
       CALL WriteParaviewFile(this)
  END SUBROUTINE InitFileIO_vtk

  SUBROUTINE WriteParaviewFile(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    !------------------------------------------------------------------------!
    INTEGER           :: i,k
    REAL              :: ftime
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! write a pvd-file: this is a "master" file of all timesteps of all ranks
#ifdef HAVE_VTK
      OPEN(this%unit, FILE=TRIM(this%filename)//'.pvd', &
           STATUS     = 'REPLACE',      &
#ifndef NOSTREAM
           ACCESS     = 'STREAM' ,   &
#else
           FORM='UNFORMATTED',&
#endif
           ACTION     = 'WRITE',        &
           POSITION   = 'REWIND',       &
           IOSTAT     = this%error)
       IF (this%error.NE. 0) CALL Error(this,"InitFileIO_vtk","Can't open pvd-file")
#endif
       WRITE(this%unit, IOSTAT=this%error)'<?xml version="1.0"?>'//LF &
             //'<VTKFile type="Collection" version="0.1" byte_order=' &
             //this%endianness//'>'//LF &
             //repeat(' ',2)//'<Collection>'//LF
       IF (this%error.NE. 0) CALL Error(this,"InitFileIO_vtk","Can't write pvd-file")
    
       ftime = 0.0
       DO k=0,this%cycles-1 
         ftime = this%stoptime/(this%cycles-1)*k
#ifdef PARALLEL       
         DO i=0,GetNumProcs(this)-1
           WRITE(this%linebuf,fmt='(A,E11.5,A,I4.4,A)',IOSTAT=this%error)&
              repeat(' ',4)//'<DataSet timestep="',&
              ftime,'" part="', i ,'" file="'//TRIM(GetBasename(this,i))//'"/>' // LF
           WRITE(this%unit, IOSTAT=this%error) TRIM(this%linebuf)
         END DO
#else
         WRITE(this%linebuf,fmt='(A,E11.5,A,I4.4,A)',IOSTAT=this%error)repeat(' ',4)//&
            '<DataSet timestep="',ftime, &
            '" part="0" file="'//TRIM(GetBasename(this,i))//'"/>' // LF
         WRITE(this%unit, IOSTAT=this%error) TRIM(this%linebuf)
#endif
       END DO
       WRITE(this%unit, IOSTAT=this%error)repeat(' ',2)//'</Collection>'//LF &
              //'</VTKFile>'//LF
       CALL CloseFile_vtk(this)

  END SUBROUTINE WriteParaviewFile


  !> \public Determines precision of real numbers in bytes
  !!
  !! Determines the precision (aka size) of a real number.
  !! Single precision (4 bytes), double precision (8 bytes) and 
  !! quad precision (16 bytes) are possible results.
  SUBROUTINE GetPrecision_vtk(this, realsize)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this      !< \param [in,out] this fileio type
    INTEGER           :: realsize  !< \param [out] realsize size of real (byte)
    !------------------------------------------------------------------------!
    CHARACTER(LEN=4)  :: cTIPO4
    CHARACTER(LEN=8)  :: cTIPO8
    CHARACTER(LEN=16) :: cTIPO16
    REAL              :: rTIPO1, rTIPO2
    INTEGER           :: k,err,iTIPO
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: realsize
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!

    !endianness
    k = BIT_SIZE(iTIPO)/8
    !HOW big is a REAL?
    rTIPO1 = ACOS(0.0)
    cTIPO4 = transfer(rTIPO1,cTIPO4)
    rTIPO2 = transfer(cTIPO4,rTIPO2)
    IF (rTIPO2 == rTIPO1) THEN
       realsize = 4
    ELSE
       cTIPO8 = transfer(rTIPO1,cTIPO8)
       rTIPO2 = transfer(cTIPO8,rTIPO2)
       IF (rTIPO2 == rTIPO1) THEN
          realsize = 8
       ELSE
          cTIPO16 = transfer(rTIPO1,cTIPO16)
          rTIPO2 = transfer(cTIPO16,rTIPO2)
          IF (rTIPO2 == rTIPO1) THEN
             realsize = 16
          ELSE
             CALL Error(this, "GetRealsize_vtk", "Could not estimate size of float type")
          END IF
       END IF
    END IF

  END SUBROUTINE GetPrecision_vtk

  !> \public Determines the endianness of the system
  !!
  !! Determines the the endianess of the system (big or little endian) 
  SUBROUTINE GetEndianness_vtk(this, res, littlestr, bigstr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this         !< \param [in,out] this fileio type
    CHARACTER(LEN=*)  :: res          !< \param [out] res result string
    CHARACTER(LEN=*)  :: littlestr    !< \param [in] littlestr little endian str
    CHARACTER(LEN=*)  :: bigstr       !< \param [in] bigstr big endian str
    !------------------------------------------------------------------------!
    INTEGER           :: k,err,iTIPO
    CHARACTER, POINTER:: cTIPO(:)
    !------------------------------------------------------------------------!
    INTENT(IN)        :: littlestr, bigstr
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: res
    !------------------------------------------------------------------------!

    !endianness
    k = BIT_SIZE(iTIPO)/8
    ALLOCATE(cTIPO(k),STAT = err)
       IF (err.NE.0) &
         CALL Error(this, "GetEndianness_vtk", "Unable to allocate memory.")
    cTIPO(1)='A'
    !cTIPO(2:k-1) = That's of no importance.
    cTIPO(k)='B'

    iTIPO = transfer(cTIPO, iTIPO)
    DEALLOCATE(cTIPO)
    !Test of 'B'=b'01000010' ('A'=b'01000001')
    IF (BTEST(iTIPO,1)) THEN
       write(res,'(A)',IOSTAT=this%error)bigstr
    ELSE
       write(res,'(A)',IOSTAT=this%error)littlestr
    END IF
  END SUBROUTINE GetEndianness_vtk

  !> Creates a list of all data arrays which will be written to file
  !!
  !! Therefore it ignores all arrays with coordinates and checks if the data 
  !! arrays are of the dimension of the mesh.
  RECURSIVE SUBROUTINE GetOutputlist(this,Mesh,node,k,prefix)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)       :: this   !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)         :: Mesh   !< \param [in] mesh mesh type
    TYPE(Dict_TYP),POINTER :: node   !< \param [in,out] node pointer to (sub-)dict
    INTEGER                :: k      !< \param [in,out] k number of data arrays
    !> \param [in,out] prefix namespace (path) to sub-dict
    CHARACTER(LEN=*),OPTIONAL :: prefix        
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: dir
    CHARACTER(LEN=MAXKEY)  :: key
    REAL,DIMENSION(:,:),POINTER :: dummy2
    REAL,DIMENSION(:,:,:),POINTER :: dummy3
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    INTENT(INOUT)     :: this,k,prefix
    !------------------------------------------------------------------------!
    DO WHILE(ASSOCIATED(node))
      IF(GetDataType(node) .EQ. DICT_DIR) THEN
      ! recursion
        IF(PRESENT(prefix)) THEN
           key = TRIM(prefix)//'/'//TRIM(GetKey(node))
        ELSE
           key = '/'//TRIM(GetKey(node))
        END IF
        CALL GetAttr(node,GetKey(node),dir)
        CALL GetOutputlist(this,Mesh,dir,k,key)
      ELSE IF (GetKey(node) .EQ. "value") THEN
         IF (k+1 .GT. MAXCOLS) CALL Error(this,"GetOutputlist_vtk","reached MAXCOLS")
         IF(PRESENT(prefix)) THEN
           this%output(k+1)%key = TRIM(prefix)//'/'//TRIM(GetKey(node))
         ELSE
           this%output(k+1)%key = '/'//TRIM(GetKey(node))
         END IF
         ! check dim of array
         SELECT CASE(GetDataType(node,TRIM(GetKey(node))))
         CASE(DICT_REAL_TWOD)
           CALL GetAttr(node,TRIM(GetKey(node)),dummy2)
           ! if not => reset k
           IF (.NOT.(LBOUND(dummy2,1) .NE. Mesh%IGMIN .OR. &
               LBOUND(dummy2,2) .NE. Mesh%JGMIN .OR. &
               UBOUND(dummy2,1) .NE. Mesh%IGMAX .OR. &
               UBOUND(dummy2,2) .NE. Mesh%JGMAX)) k = k+1
         CASE(DICT_REAL_THREED)
           CALL GetAttr(node,TRIM(GetKey(node)),dummy3)
           ! if not => reset k
           IF (.NOT.(LBOUND(dummy3,1) .NE. Mesh%IGMIN .OR. &
               LBOUND(dummy3,2) .NE. Mesh%JGMIN .OR. &
               UBOUND(dummy3,1) .NE. Mesh%IGMAX .OR. &
               UBOUND(dummy3,2) .NE. Mesh%JGMAX)) k = k+1
         END SELECT
      END IF
      node=>GetNext(node)
    END DO
  END SUBROUTINE GetOutputlist

  !> \public Specific routine to open a file for vtk I/O
  !!
  SUBROUTINE OpenFile_vtk(this,action)
     IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this    !< \param [in,out] this fileio type
    INTEGER          :: action  !< \param [in] action mode of file access
    !------------------------------------------------------------------------!
    INTENT(IN)       :: action
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
#ifdef HAVE_VTK
    SELECT CASE(action)
    CASE(READONLY)
       OPEN(this%unit, FILE=GetFilename(this), &
         STATUS     = 'OLD',          &
#ifndef NOSTREAM 
         ACCESS     = 'STREAM' ,   &
#else
         FORM='UNFORMATTED',&
#endif
         action     = 'READ',         &
         POSITION   = 'REWIND',       &
         iostat     = this%error)
    CASE(READEND)
       open(this%unit, FILE=GetFilename(this), &
         STATUS     = 'OLD',          &
#ifndef NOSTREAM
         ACCESS     = 'STREAM' ,   &
#else
         FORM='UNFORMATTED',&
#endif
         action     = 'READ',         &
         POSITION   = 'APPEND',       &
         iostat     = this%error)
    CASE(REPLACE)
       open(this%unit, FILE=GetFilename(this), &
         STATUS     = 'REPLACE',      &
#ifndef NOSTREAM
         ACCESS     = 'STREAM' ,   &
#else
         FORM='UNFORMATTED',&
#endif
         action     = 'WRITE',        &
         POSITION   = 'REWIND',       &
         iostat     = this%error)
#ifdef PARALLEL
    ! open pvts-file
      IF (GetRank(this).EQ.0) THEN
        this%extension='pvts'
        OPEN(this%unit+100, FILE=GetFilename(this,-1), &
           STATUS     = 'REPLACE',      &
#ifndef NOSTREAM
           ACCESS     = 'STREAM' ,   &
#else
           FORM='UNFORMATTED',&
#endif
           ACTION     = 'WRITE',        &
           POSITION   = 'REWIND',       &
           IOSTAT     = this%error)
        this%extension='vts'
      END IF
#endif
    CASE(APPEND)
       open(this%unit, FILE=GetFilename(this), &
         STATUS     = 'OLD',          &
#ifndef NOSTREAM
         ACCESS     = 'STREAM' ,   &
#else
         FORM='UNFORMATTED',&
#endif
         action     = 'READWRITE',    &
         POSITION   = 'APPEND',       &
         iostat     = this%error)
#ifdef PARALLEL
    ! open pvts-file
      IF (GetRank(this).EQ.0) THEN
        this%extension='pvts'
        OPEN(this%unit+100, FILE=GetFilename(this,-1), &
           STATUS     = 'OLD',      &
#ifndef NOSTREAM
           ACCESS     = 'STREAM' ,   &
#else
           FORM='UNFORMATTED',&
#endif
           ACTION     = 'READWRITE',        &
           POSITION   = 'APPEND',       &
           IOSTAT     = this%error)
        this%extension='vts'
      END IF
#endif

    CASE DEFAULT
       CALL Error(this,"OpenFile","Unknown access mode.")
    END SELECT
    IF (this%error.NE. 0) CALL Error(this,"OpenFile_vtk","Can't open file")
#endif
  END SUBROUTINE OpenFile_vtk

  !> \public Specific routine to close a file for vtk I/O
  !!
  SUBROUTINE CloseFile_vtk(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!

    CLOSE(this%unit,IOSTAT=this%error)
    IF(this%error .NE. 0) CALL ERROR(this, "CloseFileIO_vtk", "Can't close file")

#ifdef PARALLEL
    IF (GetRank(this) .EQ. 0 ) THEN
      CLOSE(this%unit+100,IOSTAT=this%error)
      IF(this%error .NE. 0) CALL ERROR(this, "CloseFileIO_vtk", "Can't close pvts file")
    END IF
#endif

  END SUBROUTINE CloseFile_vtk

  !> Extract a subkey of a string
  !! 
  !! Extract a subkey from a string (key) of type "abcd/SUBKEY/efgh"
  FUNCTION GetSubKey(key) RESULT(subkey)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CHARACTER(LEN=*) :: key         
    !------------------------------------------------------------------------!
    CHARACTER(LEN=MAXKEY) :: subkey 
    !------------------------------------------------------------------------!
    INTENT(IN)       :: key
    !------------------------------------------------------------------------!
    ! format of key like: "abcd/SUBKEY/efgh"
    subkey = key(1:SCAN(key,"/",.TRUE.)-1)
    subkey = subkey(SCAN(subkey,"/",.TRUE.)+1:)
  END FUNCTION GetSubKey

  !> \public Writes XML header to file
  !!
  SUBROUTINE WriteHeader_vtk(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this       !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh       !< \param [in] mesh mesh type
    TYPE(Physics_TYP) :: Physics    !< \param [in] physics physics type
    TYPE(Dict_TYP),POINTER :: IO    !< \param [in,out] IO I/O dictionary
    TYPE(Dict_TYP),POINTER :: config!< \param [in,out] config config dictionary
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    this%ioffset = 0

    WRITE(this%unit, IOSTAT=this%error)&
        '<?xml version="1.0"?>' &
        //LF//'<VTKFile type="StructuredGrid" version="0.1" byte_order=' & 
        //this%endianness//'>'

#ifdef PARALLEL
    IF (GetRank(this) .EQ. 0 ) THEN
       WRITE(this%unit+100, IOSTAT=this%error)&
        '<?xml version="1.0"?>' &
        //LF//'<VTKFile type="PStructuredGrid" version="0.1" byte_order=' & 
        //this%endianness//'>'
    END IF
#endif


  END SUBROUTINE WriteHeader_vtk

  !> \public Reads the header (not yet implemented)
  !!
  SUBROUTINE ReadHeader_vtk(this,success)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    LOGICAL          :: success
    !------------------------------------------------------------------------!
    INTENT(OUT)      :: success
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    success = .FALSE.
  END SUBROUTINE ReadHeader_vtk

  !> \public Writes the timestep (not yet implemented)
  !!
  SUBROUTINE WriteTimestamp_vtk(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    REAL             :: time
    !------------------------------------------------------------------------!
    INTEGER          :: i
    INTENT(IN)       :: time
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
  END SUBROUTINE WriteTimestamp_vtk

  !> \public Reads the timestep (not yet implemented)
  !!
  SUBROUTINE ReadTimestamp_vtk(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    REAL             :: time
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: time
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    time = 0.0
  END SUBROUTINE ReadTimestamp_vtk

  !> \public Writes all desired data arrays to a file 
  !!
  SUBROUTINE WriteDataset_vtk(this,Mesh,Physics,Fluxes,Timedisc,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this      !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh      !< \param [in] mesh mesh type
    TYPE(Physics_TYP) :: Physics   !< \param [in] physics physics type
    TYPE(Fluxes_TYP)  :: Fluxes    !< \param [in] fluxes fluxes type
    TYPE(Timedisc_TYP):: Timedisc  !< \param [in] timedisc timedisc type
    TYPE(Dict_TYP),POINTER :: IO   !< \param [in,out] IO I/O dictionary
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k,m, spos, epos, n
    REAL,DIMENSION(:,:),POINTER :: dummy2
    REAL,DIMENSION(:,:,:),POINTER :: dummy3
    REAL,DIMENSION(:,:,:,:),POINTER :: dummy4
    TYPE(Dict_TYP),POINTER :: node
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Fluxes,Timedisc
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    WRITE(this%linebuf,fmt='(6(I7))',IOSTAT=this%error)&
         Mesh%IMIN,Mesh%IMAX+1,Mesh%JMIN,Mesh%JMAX+1,1,1
    WRITE(this%linebuf2,fmt='(ES14.7)',IOSTAT=this%error)&
         Timedisc%time

    WRITE(this%unit, IOSTAT=this%error)&
          LF//repeat(' ',2)//'<StructuredGrid WholeExtent="' &
                //TRIM(this%linebuf)//'">' &
        //LF//repeat(' ',4)//'<FieldData>' &
        //LF//repeat(' ',6)//'<DataArray type="Float64" Name="TIME"' &
                           //' NumberOfTuples="1" format="ascii">' &
        //LF//repeat(' ',8)//TRIM(this%linebuf2) &
        //LF//repeat(' ',6)//'</DataArray>' &
        //LF//repeat(' ',4)//'</FieldData>' &
        //LF//repeat(' ',4)//'<Piece Extent="'//TRIM(this%linebuf)//'">' &
        //LF//repeat(' ',6)//'<Points>' &
        //LF//repeat(' ',8)//'<DataArray type='//this%realfmt &
                //' NumberOfComponents="3" Name="Point"' &
                //' format="appended" offset="0">' &
        //LF//repeat(' ',8)//'</DataArray>' &
        //LF//repeat(' ',6)//'</Points>' &
        //LF//repeat(' ',6)//'<CellData>'&
        //LF

#ifdef PARALLEL
    IF (GetRank(this) .EQ. 0 ) THEN
       WRITE(this%linebuf,fmt='(6(I7))',IOSTAT=this%error)&
         1,Mesh%INUM+1,1,Mesh%JNUM+1,1,1
       WRITE(this%unit+100, IOSTAT=this%error)&
          LF//repeat(' ',2)//'<PStructuredGrid WholeExtent="' &
                //TRIM(this%linebuf)//'">' &
        //LF//repeat(' ',4)//'<PFieldData>' &
        //LF//repeat(' ',6)//'<DataArray type="Float64" Name="TIME"' &
                           //' NumberOfTuples="1" format="ascii">' &
        //LF//repeat(' ',8)//TRIM(this%linebuf2) &
        //LF//repeat(' ',6)//'</DataArray>' &
        //LF//repeat(' ',4)//'</PFieldData>' &
        //LF//repeat(' ',4)//'<PPoints>' &
        //LF//repeat(' ',6)//'<DataArray type='//this%realfmt &
                //' NumberOfComponents="3" Name="Point"/>' &
        //LF//repeat(' ',4)//'</PPoints>' &
        //LF//repeat(' ',4)//'<PCellData>'&
        //LF
    END IF
#endif


    N_Byte = 3*(Mesh%IMAX-Mesh%IMIN+2)*(Mesh%JMAX-Mesh%JMIN+2)*this%realsize
    this%ioffset = this%ioffset + N_Byte + 4 !sizeof(N_Byte) must be 4!!!!!  
 
    N_Byte = (Mesh%IMAX-Mesh%IMIN+1)*(Mesh%JMAX-Mesh%JMIN+1)*this%realsize
    DO k = 2, this%cols
       SELECT CASE(GetDataType(IO,TRIM(this%output(k)%key)))
       CASE(DICT_REAL_TWOD) 
         n = 1
       CASE(DICT_REAL_THREED)
         CALL GetAttr(IO,TRIM(this%output(k)%key),dummy3)
         n = SIZE(dummy3,DIM=3)
       CASE(DICT_REAL_FOURD)
         CALL GetAttr(IO,TRIM(this%output(k)%key),dummy4)
         n = SIZE(dummy4,DIM=3)*SIZE(dummy4,DIM=4)
       CASE DEFAULT
          CALL ERROR(this, "WriteDataset_vtk", "wrong type of output data")
       END SELECT

       WRITE(this%linebuf,fmt='(I8)',IOSTAT=this%error)this%ioffset
       WRITE(this%buf,fmt='(I3)',IOSTAT=this%error)n
       WRITE(this%unit,IOSTAT=this%error)&
           repeat(' ',6)//'<DataArray type='//this%realfmt//&
              ' NumberOfComponents="',TRIM(this%buf),&
              '" Name="'//TRIM(GetSubKey(this%output(k)%key))//&
              '" format="appended" offset="',TRIM(this%linebuf),'"/>'//LF

#ifdef PARALLEL
    IF (GetRank(this) .EQ. 0 ) THEN
      WRITE(this%unit+100,IOSTAT=this%error)&
           repeat(' ',6)//'<DataArray type='//this%realfmt//&
              ' NumberOfComponents="',TRIM(this%buf),&
              '" Name="'//TRIM(GetSubKey(this%output(k)%key))//'"/>'//LF
    END IF
#endif

       IF (k < this%cols) THEN 
          this%ioffset = this%ioffset + n*N_Byte + 4 !sizeof(N_Byte) must be 4!!
       END IF
    END DO

    WRITE(this%unit,IOSTAT=this%error)&
             repeat(' ',6)//'</CellData>' &
       //LF//repeat(' ',4)//'</Piece>' &
       //LF//repeat(' ',2)//'</StructuredGrid>' &
       //LF//repeat(' ',2)//'<GlobalData>'

#ifdef PARALLEL
    IF (GetRank(this) .EQ. 0 ) THEN
      WRITE(this%unit+100, IOSTAT=this%error)&
         repeat(' ',4)//'</PCellData>'//LF

      DO i=0,GetNumProcs(this)-1
        WRITE(this%linebuf,fmt='(A)',IOSTAT=this%error)&
           repeat(' ',4)//'<Piece Extent="'//TRIM(this%extent(i))&
           //'" Source="'//TRIM(GetBasename(this,i))//'"/>' // LF
        WRITE(this%unit+100, IOSTAT=this%error) TRIM(this%linebuf)
      END DO

      WRITE(this%unit+100, IOSTAT=this%error)&
        repeat(' ',2)//'</PStructuredGrid>' &
        //LF//'</VTKFile>'//LF
       
      CLOSE(this%unit+100,IOSTAT=this%error)
      IF(this%error .NE. 0) CALL ERROR(this, "WriteDataset_vtk", "Can't close pvts file")
    END IF
#endif


    ! get boundary fluxes across all 4 boundaries
    DO i=1,4
       ! in PARALLEL mode the result is sent to process with rank 0
       this%bflux(:,i) = GetBoundaryFlux(Fluxes,Mesh,Physics,i)
    
       WRITE(this%unit,IOSTAT=this%error)&
              LF//repeat(' ',4)//'<DataArray type='//this%realfmt//&
              'Name="'//TRIM(fluxkey(i))//'" format="ascii" >' &
              //LF//repeat(' ',6)
       DO k=1,Physics%VNUM
           WRITE(this%linebuf((k-1)*20+1:k*20),fmt='(E16.9,A)',IOSTAT=this%error)&
                this%bflux(k,i),repeat(' ',4)
       END DO

       WRITE(this%unit,IOSTAT=this%error)TRIM(this%linebuf(1:Physics%VNUM*20)) &
              //LF//repeat(' ',4)//'</DataArray>'

    END DO
      
    !size of mesh
    N_Byte = 3*(Mesh%IMAX-Mesh%IMIN+2)*(Mesh%JMAX-Mesh%JMIN+2)*this%realsize
    WRITE(this%unit,IOSTAT=this%error)&
       LF//repeat(' ',2)//'</GlobalData>' &
       //LF//'<AppendedData encoding="raw">' &
       //LF//'_',N_Byte

    !get mesh corners
    CALL GetAttr(IO,TRIM(this%output(1)%key),dummy4)
    DO j=Mesh%JMIN,Mesh%JMAX+1
       DO i=Mesh%IMIN,Mesh%IMAX+1
          ! only use SOUTH/WEST corner '1'
          this%vtktemp2(1,i,j) = dummy4(i,j,1,1)
          this%vtktemp2(2,i,j) = dummy4(i,j,1,2)
          ! set z-position to zero (vtk need 3D coords)
          this%vtktemp2(3,i,j) = 0.0
       END DO
    END DO

    ! write mesh
    WRITE(this%unit,IOSTAT=this%error)this%vtktemp2(1:3,:,:)

    DO k = 2, this%cols
       SELECT CASE(GetDataType(IO,TRIM(this%output(k)%key)))
       CASE(DICT_REAL_TWOD)
         CALL GetAttr(IO,TRIM(this%output(k)%key),dummy2)
         N_Byte  = (Mesh%IMAX-Mesh%IMIN+1)*(Mesh%JMAX-Mesh%JMIN+1)*this%realsize
         WRITE(this%unit,IOSTAT=this%error) N_Byte, &
            dummy2(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)
       CASE(DICT_REAL_THREED)
         CALL GetAttr(IO,TRIM(this%output(k)%key),dummy3)
         n = SIZE(dummy3,DIM=3)
         N_Byte  = n*(Mesh%IMAX-Mesh%IMIN+1)*(Mesh%JMAX-Mesh%JMIN+1)*this%realsize

         ! copy for c-like output
         DO j=Mesh%JMIN,Mesh%JMAX
           DO i=Mesh%IMIN,Mesh%IMAX
                this%vtktemp2(1,i,j)=dummy3(i,j,1)
                this%vtktemp2(2,i,j)=dummy3(i,j,2)
                DO m=3, n
                   this%vtktemp2(m,i,j)=dummy3(i,j,m)
                END DO
           END DO
         END DO

         WRITE(this%unit,IOSTAT=this%error) N_Byte, &
            this%vtktemp2(1:n,Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)
       CASE DEFAULT
          CALL ERROR(this, "WriteDataset_vtk", "wrong type of output data")
       END SELECT
    END DO
    WRITE(this%unit,IOSTAT=this%error)LF//repeat(' ',2)//&
         '</AppendedData>'//LF//'</VTKFile>'//LF
      
  END SUBROUTINE WriteDataset_vtk

  !> \public Closes the file I/O
  !!
  SUBROUTINE CloseFileIO_vtk(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this   !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%vtktemp,this%vtktemp2,this%output,this%bflux) 
  END SUBROUTINE CloseFileIO_vtk
END MODULE fileio_vtk
