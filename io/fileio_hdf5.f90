!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fileio_hdf5.f90                                                   #
!#                                                                           #
!# Copyright (C) 2012                                                        #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
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
!! \brief I/O for HDF5 files 
!!
!! \extends fileio_common
!! \ingroup fileio
!----------------------------------------------------------------------------!
MODULE fileio_hdf5
#ifdef HAVE_HDF5_MOD
  USE hdf5
  USE h5lt
#endif
#ifdef PARALLEL
#ifdef HAVE_MPI_MOD
  USE mpi
#endif
#endif
  USE fileio_common
  USE geometry_common, ONLY : Geometry_TYP, GetName, GetType
  USE geometry_generic, ONLY : Convert2Cartesian
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP, GetName, GetType
! Fixme remove hack to load bccsound and fcsound arrays
  USE physics_generic, ONLY : EULER2D_ISOTHERM
  USE timedisc_common, ONLY : Timedisc_TYP
  USE fluxes_common, ONLY : Fluxes_TYP
  USE fluxes_generic, ONLY : GetBoundaryFlux
  USE common_dict
  IMPLICIT NONE
#ifdef PARALLEL
#ifdef HAVE_MPIF_H
  INCLUDE 'mpif.h'
#endif
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER   :: DEFAULT_REAL
  INTEGER   :: DEFAULT_INT
#ifdef HAVE_HDF5_MOD
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE WriteNode_hdf5
    MODULE PROCEDURE WriteNode_hdf5_1, WriteNode_hdf5_2, WriteNode_hdf5_3, &
        WriteNode_hdf5_4, WriteNode_hdf5_5, WriteNode_hdf5_6, &
        WriteNode_hdf5_7, WriteNode_hdf5_8, WriteNode_hdf5_9, &
        WriteNode_hdf5_10
  END INTERFACE WriteNode_hdf5
  !> \endcond
#endif
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       FileIO_TYP, &
       ! constants
       ! methods
#ifdef HAVE_HDF5_MOD
       OpenFile_hdf5, &
       CloseFile_hdf5, &
       WriteHeader_hdf5, &
       ReadHeader_hdf5, &
       WriteTimestamp_hdf5, &
       ReadTimestamp_hdf5, &
       WriteDataset_hdf5, &
       ReadDataset_hdf5, &
       GetConfig_hdf5, &
       CloseFileio_hdf5, &
#endif
       InitFileio_hdf5
  !--------------------------------------------------------------------------!

CONTAINS
  !> \public Constructor for the xdmf file I/O 
  !!
  !! Initilizes the file I/O type, filename, stoptime, number of outputs, 
  !! number of files, unit number, config as a dict
  SUBROUTINE InitFileio_hdf5(this,Mesh,Physics,fmt,fpath,filename,stoptime,dtwall,&
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
    REAL              :: dummy
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,fmt,fpath,filename,stoptime,count,fcycles,unit
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
#ifndef HAVE_HDF5_MOD
    CALL Error(this, "InitFileio_hdf5", "No hdf5 fortran module available. "&
        //"Set Flag --with-hdf5 and make sure hdf5 has been compiled with "&
        //"--enable-fortran.")
#endif

    CALL InitFileIO(this,fmt,"HDF5",fpath,filename,"h5",fcycles,.FALSE.,unit)
    this%count = count
    this%stoptime = stoptime
#ifdef HAVE_HDF5_MOD
    CALL h5open_f(this%error)
    SELECT CASE (SELECTED_REAL_KIND(PRECISION(dummy)))
    CASE(4)
        DEFAULT_REAL = H5T_NATIVE_REAL
    CASE(8)
        DEFAULT_REAL = H5T_NATIVE_DOUBLE
    CASE DEFAULT
        CALL Error(this,"OpenFile_hdf5","Cannot determine real type.")
    END SELECT
    DEFAULT_INT = H5T_NATIVE_INTEGER
#ifdef PARALLEL
    IF(MOD(Mesh%INUM,Mesh%dims(1))>1.0E-10 &
       .OR.MOD(Mesh%JNUM,Mesh%dims(2))>1.0E-10) THEN
       CALL Error(this, "OpenFile_hdf5","Each of INUM and JNUM have to be a"&
            //" multiple of the MPI domain decomposition components")
    END IF
#endif     
#endif
  END SUBROUTINE InitFileIO_hdf5

#ifdef HAVE_HDF5_MOD
  !> \public Specific routine to open a file for hdf5 I/O
  !!
  SUBROUTINE OpenFile_hdf5(this,action)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this    !< \param [in,out] this fileio type
    INTEGER          :: action  !< \param [in] action mode of file access
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER,PARAMETER:: comm = MPI_COMM_WORLD
    INTEGER,PARAMETER:: info = MPI_INFO_NULL
#endif
    INTEGER(HID_T)   :: plistid
    !------------------------------------------------------------------------!
    INTENT(IN)       :: action
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plistid, this%error)
    CALL h5pset_fapl_mpio_f(plistid, comm, info, this%error)
#else
    plistid = H5P_DEFAULT_F
#endif
    SELECT CASE(action)
    CASE(READONLY,READEND)
        CALL h5fopen_f(GetFilename(this), &
                       H5F_ACC_RDONLY_F, &
                       this%fid, &
                       this%error, &
                       access_prp = plistid)
    CASE(REPLACE)
        CALL h5fcreate_f(GetFilename(this), &
                         H5F_ACC_TRUNC_F, &
                         this%fid, &
                         this%error, &
                         access_prp = plistid)
    CASE(APPEND)
        CALL h5fopen_f(GetFilename(this), &
                       H5F_ACC_RDWR_F, &
                       this%fid, &
                       this%error, &
                       access_prp = plistid)
    CASE DEFAULT
       CALL Error(this,"OpenFile_hdf5","Unknown access mode.")
    END SELECT
#ifdef PARALLEL
    CALL h5pclose_f(plistid, this%error)
#endif
#ifdef PARALLEL
    ! Create porperty list for collective dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, this%xferid, this%error)
    CALL h5pset_dxpl_mpio_f(this%xferid, H5FD_MPIO_COLLECTIVE_F, this%error)
    ! For independent write use:
    !CALL h5pset_dxpl_mpio_f(this%xferid, H5FD_MPIO_INDEPENDENT_F, this%error)
#else
    this%xferid = H5P_DEFAULT_F
#endif
  END SUBROUTINE OpenFile_hdf5

  !> \public Closes the file I/O
  !!
  SUBROUTINE CloseFile_hdf5(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this  !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    CALL h5pclose_f(this%xferid, this%error)
#endif
    CALL h5fclose_f(this%fid, this%error)
  END SUBROUTINE CloseFile_hdf5


  !> \public Writes header to file
  !!
  SUBROUTINE WriteHeader_hdf5(this,Mesh,Header)
    IMPLICIT NONE
    !----------------------------------- -------------------------------------!
    TYPE(FileIO_TYP)  :: this        !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh        !< \param [in,out] mesh mesh type
    TYPE(Dict_TYP),POINTER :: Header !< \param [in,out] header dict type
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    INTENT(IN)        :: Mesh
    !------------------------------------------------------------------------!
    CALL WriteDict_hdf5(this, Mesh, Header)
  END SUBROUTINE WriteHeader_hdf5


  !> \public Reads the header
  !!
  SUBROUTINE ReadHeader_hdf5(this,Header,success)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this        !< \param [in,out] this fileio type
    TYPE(Dict_TYP),POINTER :: Header !< \param [in,out] header dict type
    LOGICAL           :: success     !< \param [out] success
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: success
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    NULLIFY(Header)
    CALL ReadDict_hdf5(this, Header)
    success = .True.
  END SUBROUTINE ReadHeader_hdf5

  !> \public Writes timestamp to file
  !!
  !! This is implemented in WriteDataset, because Mesh is needed
  SUBROUTINE WriteTimestamp_hdf5(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    REAL             :: time
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    INTENT(IN)       :: time
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    ! This is implemented in WriteDataset, because Mesh is needed
  END SUBROUTINE WriteTimestamp_hdf5

  !> \public Reads the timestep (not yet implemented)
  !!
  SUBROUTINE ReadTimestamp_hdf5(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    REAL             :: time
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    INTENT(OUT)      :: time
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    time = 0.0
  END SUBROUTINE ReadTimestamp_hdf5

  !> \public Writes all desired data arrays to a file 
  !!
  SUBROUTINE WriteDataset_hdf5(this,Mesh,Physics,Timedisc,Fluxes,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this      !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh      !< \param [in] mesh mesh type
    TYPE(Physics_TYP) :: Physics   !< \param [in] physics physics type
    TYPE(Fluxes_TYP)  :: Fluxes    !< \param [in] fluxes fluxes type
    TYPE(Timedisc_TYP):: Timedisc  !< \param [in] timedisc timedisc type
    TYPE(Dict_TYP),POINTER :: IO   !< \param [in,out] IO I/O dictionary
    !------------------------------------------------------------------------!
    CHARACTER(LEN=MAX_CHAR_LEN) :: path = ""
    TYPE(Dict_TYP),POINTER :: tmp
    REAL, DIMENSION(:,:),POINTER :: bflux
    INTEGER           :: k,status
    CHARACTER(LEN=32) :: str = ""
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Timedisc,Fluxes
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    WRITE(path, "(I4.4)")this%step
    tmp => Dict(TRIM(path) / IO)
    CALL WriteDict_hdf5(this,Mesh,tmp)
    CALL WriteNode_hdf5(this,Mesh,TRIM(path)//"/timedisc/time",Timedisc%time)
    DEALLOCATE(tmp,stat=status)
    ! nesh: Error code 195 means deallocation of a unassociated pointer.
    ! Checking with ASSOCIATED says it is indeed associated. So ignore error
    ! 195
    IF((status.NE.0).AND.(status.NE.195)) THEN
      WRITE(str, *)status
      CALL Error(this,"WriteDataset_hdf5", "status = "//TRIM(str))
    END IF
  END SUBROUTINE WriteDataset_hdf5

  !> \public Reads all data arrays from file 
  !!
  SUBROUTINE ReadDataset_hdf5(this,Mesh,Physics,Timedisc,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this      !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh      !< \param [in] mesh mesh type
    TYPE(Physics_TYP) :: Physics   !< \param [in] physics physics type
    TYPE(Timedisc_TYP):: Timedisc  !< \param [in] timedisc timedisc type
    TYPE(Dict_TYP),POINTER :: IO   !< \param [in,out] IO I/O dictionary
    !------------------------------------------------------------------------!
    REAL,DIMENSION(:,:),POINTER &
                      :: ptr
    REAL,DIMENSION(:,:,:),POINTER &
                      :: ptr2
    CHARACTER(LEN=MAX_CHAR_LEN) :: path, keyIO, key
    INTEGER           :: i
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    INTENT(INOUT)     :: this,Timedisc,Physics
    !------------------------------------------------------------------------!
    !CALL ReadDict_hdf5(this, config, Mesh=Mesh)
    !CALL PrintDict(config)
    path = "/timedisc/"
    key = ""
    keyIO = ""
    DO i = 1,Physics%VNUM
        keyIO = TRIM(path)//TRIM(Physics%pvarname(i))//"/value"
        key = "/0000"//TRIM(keyIO)
        CALL GetAttr(IO, keyIO, ptr)
        CALL ReadNode_hdf5_7(this, Mesh, key, ptr)
        
        ! This is a workaround. Should already be written by assigning the data
        ! to ptr. But it doesn't work at the moment. So let us copy them again.
        Timedisc%pvar(:,:,i) = ptr
    END DO

    NULLIFY(ptr)
    
    IF(GetType(Physics).EQ.EULER2D_ISOTHERM) THEN
        ! Small hack to load bccsound and fcsound arrays.
        ! It is planned to load everything in the IO array.
        keyIO = "/physics/bccsound/value"
        key = "/0000"//TRIM(keyIO)
        CALL GetAttr(IO, keyIO, ptr)
        CALL ReadNode_hdf5_7(this, Mesh, key, ptr)
        ! workaround, see above
        Physics%bccsound(:,:) = ptr

        keyIO = "/physics/fcsound/value"
        key = "/0000"//TRIM(keyIO)
        CALL GetAttr(IO, keyIO, ptr2)
        CALl ReadNode_hdf5_8(this, Mesh, key, ptr2)

        ! workaround, see above
        Physics%fcsound(:,:,:) = ptr2

        ! Copy last cells in calc domain into boundaries
        Physics%bccsound(Mesh%IGMIN,:) = Physics%bccsound(Mesh%IMIN,:)
        Physics%bccsound(Mesh%IMIN-1,:) = Physics%bccsound(Mesh%IMIN,:)
        Physics%bccsound(Mesh%IMAX+1,:) = Physics%bccsound(Mesh%IMAX,:)
        Physics%bccsound(Mesh%IGMAX,:) = Physics%bccsound(Mesh%IMAX,:)
        Physics%bccsound(:,Mesh%JGMIN) = Physics%bccsound(:,Mesh%JMIN)
        Physics%bccsound(:,Mesh%JMIN-1) = Physics%bccsound(:,Mesh%JMIN)
        Physics%bccsound(:,Mesh%JMAX+1) = Physics%bccsound(:,Mesh%JMAX)
        Physics%bccsound(:,Mesh%JGMAX) = Physics%bccsound(:,Mesh%JMAX)

        Physics%fcsound(Mesh%IGMIN,:,:) = Physics%fcsound(Mesh%IMIN,:,:)
        Physics%fcsound(Mesh%IMIN-1,:,:) = Physics%fcsound(Mesh%IMIN,:,:)
        Physics%fcsound(Mesh%IMAX+1,:,:) = Physics%fcsound(Mesh%IMAX,:,:)
        Physics%fcsound(Mesh%IGMAX,:,:) = Physics%fcsound(Mesh%IMAX,:,:)
        Physics%fcsound(:,Mesh%JGMIN,:) = Physics%fcsound(:,Mesh%JMIN,:)
        Physics%fcsound(:,Mesh%JMIN-1,:) = Physics%fcsound(:,Mesh%JMIN,:)
        Physics%fcsound(:,Mesh%JMAX+1,:) = Physics%fcsound(:,Mesh%JMAX,:)
        Physics%fcsound(:,Mesh%JGMAX,:) = Physics%fcsound(:,Mesh%JMAX,:)

        ! known edges get known values
        Physics%fcsound(Mesh%IMIN-1,:,2) = Physics%fcsound(Mesh%IMIN,:,1)
        Physics%fcsound(Mesh%IMAX+1,:,1) = Physics%fcsound(Mesh%IMAX,:,2)
        Physics%fcsound(:,Mesh%JMIN-1,3) = Physics%fcsound(:,Mesh%JMIN,4)
        Physics%fcsound(:,Mesh%JMAX+1,4) = Physics%fcsound(:,Mesh%JMAX,3)

    END IF

    !CALL Error(this, "ReadDataset_hdf5", "Operation not supported")
  END SUBROUTINE ReadDataset_hdf5

  SUBROUTINE ReadNode_hdf5_7(this, Mesh, key, val)
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)        :: this
    TYPE(Mesh_TYP)          :: Mesh
    CHARACTER(LEN=*)        :: key
    REAL,DIMENSION(:,:),POINTER :: val
    !------------------------------------------------------------------------!
    INTEGER,PARAMETER       :: rank = 2
    INTEGER(HSIZE_T),DIMENSION(rank) :: cdims
    INTEGER(HID_T)          :: filespace, memspace, dsetid
    !------------------------------------------------------------------------!
    INTENT(INOUT)           :: this
    INTENT(IN)              :: key, Mesh
    !------------------------------------------------------------------------!
    cdims = SHAPE(val)
    CALL TouchDataset_hdf5(this, Mesh, key, dsetid, filespace, memspace, &
                           rank, cdims)
    CALL h5dread_f(dsetid, DEFAULT_REAL, &
                   val(Mesh%IMIN:Mesh%IMAX,&
                       Mesh%JMIN:Mesh%JMAX),&
                   cdims, &
                   this%error, &
                   file_space_id = filespace, &
                   mem_space_id = memspace, &
                   xfer_prp = this%xferid)    
    
    CALL CloseDataset_hdf5(this, dsetid, filespace, memspace)

  END SUBROUTINE ReadNode_hdf5_7

  SUBROUTINE ReadNode_hdf5_8(this, Mesh, key, val)
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)        :: this
    TYPE(Mesh_TYP)          :: Mesh
    CHARACTER(LEN=*)        :: key
    REAL,DIMENSION(:,:,:),POINTER :: val
    !------------------------------------------------------------------------!
    INTEGER,PARAMETER       :: rank = 3
    INTEGER(HSIZE_T),DIMENSION(rank) :: cdims
    INTEGER(HID_T)          :: filespace, memspace, dsetid
    !------------------------------------------------------------------------!
    INTENT(INOUT)           :: this
    INTENT(IN)              :: key, Mesh
    !------------------------------------------------------------------------!
    cdims = SHAPE(val)
    CALL TouchDataset_hdf5(this, Mesh, key, dsetid, filespace, memspace, &
                           rank, cdims)
    CALL h5dread_f(dsetid, DEFAULT_REAL, &
                   val(Mesh%IMIN:Mesh%IMAX,&
                       Mesh%JMIN:Mesh%JMAX,&
                       :),&
                   cdims, &
                   this%error, &
                   file_space_id = filespace, &
                   mem_space_id = memspace, &
                   xfer_prp = this%xferid)    
   
    CALL CloseDataset_hdf5(this, dsetid, filespace, memspace)

  END SUBROUTINE ReadNode_hdf5_8

  SUBROUTINE ReadData_hdf5(this, Mesh, config, key)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)        :: this
    TYPE(Mesh_TYP)          :: Mesh
    TYPE(dict_TYP),POINTER  :: config
    CHARACTER(LEN=*)        :: key
    !------------------------------------------------------------------------!
    INTEGER                 :: dtype, rank
    INTEGER(SIZE_T)         :: tsize
    REAL,DIMENSION(:),ALLOCATABLE :: val5
    REAL,DIMENSION(:,:),ALLOCATABLE :: val7
    REAL,DIMENSION(:,:,:),ALLOCATABLE :: val8
    REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: val9
    INTEGER(HSIZE_T),DIMENSION(:),ALLOCATABLE &
                            :: cdims, dims
    INTEGER(HID_T)    :: filespace, memspace, dsetid
    !------------------------------------------------------------------------!
    INTENT(INOUT)           :: this
    INTENT(IN)              :: key, Mesh
    !------------------------------------------------------------------------!
    CALL h5ltget_dataset_ndims_f(this%fid, key, rank, this%error)
    IF(rank.LT.1) &
        CALL Error(this, "ReadData_hdf5", "rank<1 cannot be handled.")
    ALLOCATE(cdims(rank),dims(rank))
    CALL h5ltget_dataset_info_f(this%fid, key, cdims, dtype, tsize, &
        this%error)
    IF(rank.GT.1) THEN
        ! Only read the domain specific block
        cdims(1) = Mesh%IMAX-Mesh%IMIN+1+2*Mesh%GNUM
        cdims(2) = Mesh%JMAX-Mesh%JMIN+1+2*Mesh%GNUM
    END IF
    !print *,"key, dtype: ", TRIM(key), ", ", dtype
    SELECT CASE(dtype)
    !CASE(H5T_NATIVE_INTEGER)
    !CASE(H5T_NATIVE_FLOAT)
    CASE(1) !H5T_NATIVE_DOUBLE
        dims = cdims
        CALL TouchDataset_hdf5(this, Mesh, key, dsetid, filespace, memspace, &
                               rank, dims)
        SELECT CASE(rank)
        CASE(1)
           ALLOCATE(val5(cdims(1)))
           CALL h5dread_f(dsetid, DEFAULT_REAL, &
                          val5, &
                          cdims, this%error, &
                          file_space_id = filespace, &
                          mem_space_id = memspace, &
                          xfer_prp = this%xferid)
           CALL SetAttr(config, key, val5)
           DEALLOCATE(val5)
        CASE(2)
           !ALLOCATE(val7(cdims(1),cdims(2)))
           ALLOCATE(val7(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX))
           CALL h5dread_f(dsetid, DEFAULT_REAL, &
                          val7(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX), &
                          cdims, this%error, &
                          file_space_id = filespace, &
                          mem_space_id = memspace, &
                          xfer_prp = this%xferid)
           CALL SetAttr(config, key, val7)
           DEALLOCATE(val7)
        CASE(3)
           !ALLOCATE(val8(cdims(1),cdims(2),cdims(3)))
           ALLOCATE(val8(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,cdims(3)))
           CALL h5dread_f(dsetid, DEFAULT_REAL, &
                          val8(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:), &
                          cdims, this%error, &
                          file_space_id = filespace, &
                          mem_space_id = memspace, &
                          xfer_prp = this%xferid)
            CALL SetAttr(config, key, val8)
            DEALLOCATE(val8)
        CASE(4)
            !ALLOCATE(val9(cdims(1),cdims(2),cdims(3),cdims(4)))
            ALLOCATE(val9(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,cdims(3),cdims(4)))
            CALL h5dread_f(dsetid, DEFAULT_REAL, &
                           val9(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:,:), &
                           cdims, this%error, &
                           file_space_id = filespace, &
                           mem_space_id = memspace, &
                           xfer_prp = this%xferid)
            CALL SetAttr(config, key, val9)
            DEALLOCATE(val9)
        CASE DEFAULT
            CALL Error(this, "ReadData_hdf5", "Reading dataset with"&
                //" rank>4 is not supported")
        END SELECT
        CALL CloseDataset_hdf5(this, dsetid, filespace, memspace)
    CASE DEFAULT
        CALL Warning(this, "ReadData_hdf5", "Unknown dataset type.")
    END SELECT
    DEALLOCATE(cdims,dims)
  END SUBROUTINE ReadData_hdf5


  SUBROUTINE ReadAttr_hdf5(this, config, grpid, path, key)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)        :: this
    TYPE(dict_TYP),POINTER  :: config
    CHARACTER(LEN=*)        :: key, path
    !------------------------------------------------------------------------!
    CHARACTER(LEN=MAX_CHAR_LEN) :: name, str
    INTEGER                 :: type, rank
    INTEGER(SIZE_T)         :: tsize
    INTEGER(HID_T)          :: grpid
    INTEGER(HSIZE_T),DIMENSION(:),ALLOCATABLE :: dims
    INTEGER,DIMENSION(:),ALLOCATABLE :: val1
    REAL,DIMENSION(:),ALLOCATABLE :: val2
    CHARACTER(LEN=MAX_CHAR_LEN)   :: val3
    !------------------------------------------------------------------------!
    INTENT(INOUT)           :: this
    INTENT(IN)              :: key, path, grpid
    !------------------------------------------------------------------------!
    CALL h5ltget_attribute_ndims_f(grpid, ".", key, rank, this%error)
    IF(rank.GT.1) THEN
        ! rank > 1: should be a dataset. We don't use attributes with
        ! rank>1. 
        ! rank = 0: A string!
        CALL Warning(this, "ReadAttr_hdf5", "rank<>1 not supported.")
    ELSE
        name = TRIM(path) // "/" // TRIM(key)
        ALLOCATE(dims(1))
        CALL h5ltget_attribute_info_f(grpid, ".", key, dims, type, tsize, &
            this%error)
        IF((rank.EQ.0).AND.(type.NE.3)) &
            CALL Error(this, "ReadAttr_hdf5", "rank==0 and type<>3? Unknown"&
                //" state.")
        SELECT CASE(type)
        CASE(0) !H5T_INTEGER_F)
            ALLOCATE(val1(dims(1)))
            IF(dims(1).NE.1) &
                CALL Error(this, "ReadAttr_hdf5", "1D int arrays not "&
                    //"supported.")
            CALL h5ltget_attribute_int_f(grpid, ".", key, val1, this%error)
            !print *,"key, val: ", TRIM(key), ", ", val1
            CALL SetAttr(config, name, val1(1))
            DEALLOCATE(val1)
!        CASE(H5T_FLOAT_F)
!           !CALL h5ltget_attribute_float_f(grpid, ".", key, this%error, val2)
        CASE(1) !H5T_DOUBLE_F)
            ALLOCATE(val2(dims(1)))
            CALL h5ltget_attribute_double_f(grpid, ".", key, val2, this%error)
            !print *,"key, val: ", TRIM(key), ", ", val2
            IF(dims(1).EQ.1) THEN
                CALL SetAttr(config, name, val2(1))
            ELSE
                CALL SetAttr(config, name, val2)
            END IF
            DEALLOCATE(val2)
        CASE(3) !H5T_STRING_F
            val3 = ""
            CALL h5ltget_attribute_string_f(grpid, ".", key, val3, this%error)
           !print *,"key, val: ", TRIM(key), ", ", TRIM(val3)
            ! Remove trailing '\0' important!!!
            val3 = val3(1:LEN_TRIM(val3)-1)
            CALL SetAttr(config, name, TRIM(val3))
        CASE DEFAULT
!             print *,"Unknown - key, type: ", TRIM(key), ", ", type
             CALL Warning(this, "ReadAttr_hdf5", "Unknown attribute type.")
        END SELECT
        DEALLOCATE(dims)
     END IF
  END SUBROUTINE ReadAttr_hdf5


  RECURSIVE SUBROUTINE WriteDict_hdf5(this, Mesh, config, path)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)        :: this
    TYPE(Mesh_TYP)          :: Mesh
    TYPE(dict_TYP),POINTER  :: config
    CHARACTER(LEN=*),OPTIONAL :: path
    !------------------------------------------------------------------------!
    CHARACTER(LEN=MAX_CHAR_LEN) :: str, key
    TYPE(dict_TYP),POINTER   :: node
    !------------------------------------------------------------------------!
    INTENT(INOUT)            :: this
    !------------------------------------------------------------------------!
    IF(PRESENT(path)) THEN
        str = path
    ELSE
        str = "/"
    ENDIF
    node => config
    DO WHILE(ASSOCIATED(node))
        key = TRIM(str)//TRIM(GetKey(node))
        SELECT CASE(GetDataType(node))
        CASE(1)
            CALL WriteNode_hdf5(this, Mesh, key, GetVal1(node))
        CASE(2)
            CALL WriteNode_hdf5(this, Mesh, key, GetVal2(node))
        CASE(3)
            CALL WriteNode_hdf5(this, Mesh, key, GetVal3(node))
        CASE(4)
            CALL WriteNode_hdf5(this, Mesh, key, GetVal4(node))
        CASE(5)
            CALL WriteNode_hdf5(this, Mesh, key, GetVal5(node))
        CASE(6) !DICT_DIR
            CALL WriteNode_hdf5(this, Mesh, key)
            key = TRIM(key)//"/"
            CALL WriteDict_hdf5(this, Mesh, GetVal6(node), key)
        CASE(7)
            CALL WriteNode_hdf5(this, Mesh, key, GetVal7(node))
        CASE(8)
            CALL WriteNode_hdf5(this, Mesh, key, GetVal8(node))
        CASE(9)
            CALL WriteNode_hdf5(this, Mesh, key, GetVal9(node))
        CASE(10)
            CALL WriteNode_hdf5(this, Mesh, key, GetVal10(node))
        CASE DEFAULT
            CALL Error(this, "WriteDict_hdf5", "node has a unknown type")
        END SELECT
        node => GetNext(node)
    END DO
  END SUBROUTINE WriteDict_hdf5

  RECURSIVE SUBROUTINE ReadDict_hdf5(this, config, path, Mesh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)        :: this
    TYPE(Mesh_TYP),OPTIONAL :: Mesh
    TYPE(dict_TYP),POINTER  :: config
    CHARACTER(LEN=*),OPTIONAL :: path
    !------------------------------------------------------------------------!
    CHARACTER(LEN=MAX_CHAR_LEN) :: str, key, name
    INTEGER                 :: nmembers, i, type, dtype, rank
    INTEGER(HSIZE_T)        :: j, tsize
    INTEGER(HID_T)          :: grpid
    !------------------------------------------------------------------------!
    INTENT(INOUT)           :: this
    INTENT(IN)              :: path, Mesh
    !------------------------------------------------------------------------!
    IF(PRESENT(path)) THEN
        str = path
    ELSE
        str = "/"
    ENDIF
    CALL h5gn_members_f(this%fid, str, nmembers, this%error)
    DO i=0, nmembers-1
        CALL h5gget_obj_info_idx_f(this%fid, str, i, key, type, this%error)
        !print *,TRIM(key),type
        IF(PRESENT(path)) THEN
            key = TRIM(str)//"/"//TRIM(key)
        ELSE
            key = "/"//TRIM(key)
        END IF
        SELECT CASE(type)
        !CASE(0) !H5G_LINK_F
        CASE(0) !H5G_GROUP_F
            CALL ReadDict_hdf5(this, config, key, Mesh)
        CASE(1) !H5G_DATASET_F
            IF(PRESENT(Mesh)) &
                CALL ReadData_hdf5(this, Mesh, config, key)    
        !CASE(3) !H5G_TYPE_F
        CASE DEFAULT
            CALL Warning(this, "ReadDict_hdf5", "Unknown group member type.")
        END SELECT
    END DO
    CALL h5gopen_f(this%fid, str, grpid, this%error)
    CALL h5aget_num_attrs_f(grpid, nmembers, this%error)
    !print *,"Attr no: ",nmembers
    DO j=0,nmembers-1
        CALL h5aget_name_by_idx_f(grpid, ".", H5_INDEX_NAME_F, &
            H5_ITER_NATIVE_F, j, key, this%error)
        CALL ReadAttr_hdf5(this, config, grpid, str, key)
    END DO
    CALL h5gclose_f(grpid, this%error)

  END SUBROUTINE ReadDict_hdf5

  SUBROUTINE SplitPath(key, path, name)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CHARACTER(LEN=*)  :: key
    CHARACTER(LEN=MAX_CHAR_LEN) :: path, name
    !------------------------------------------------------------------------!
    INTEGER           :: k
    !------------------------------------------------------------------------!
    INTENT(IN)        :: key
    INTENT(INOUT)     :: path, name
    !------------------------------------------------------------------------!
    k = SCAN(key, '/', .TRUE.)
    IF(k.NE.0) THEN
        path = key(1:k-1)
    ELSE
        path = '/'
    ENDIF
    name = key(k+1:)
    !print *,"key, path, name: ", TRIM(key), ", ", TRIM(path), ", ", TRIM(name)
  END SUBROUTINE SplitPath

  SUBROUTINE WriteNode_hdf5_1(this,Mesh,key,val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    CHARACTER(LEN=*)  :: key
    INTEGER           :: val
    !------------------------------------------------------------------------!
    INTEGER,DIMENSION(1) :: data
    CHARACTER(LEN=MAX_CHAR_LEN) :: path, name
    INTEGER(SIZE_T)   :: size = 1
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    INTENT(IN)        :: key, val, Mesh
    !------------------------------------------------------------------------!
    CALL SplitPath(key, path, name)
    data(1) = val
    CALL h5ltset_attribute_int_f(this%fid, path, name, data, size, this%error)
  END SUBROUTINE WriteNode_hdf5_1

  SUBROUTINE WriteNode_hdf5_2(this,Mesh,key,val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    CHARACTER(LEN=*)  :: key
    REAL              :: val
    !------------------------------------------------------------------------!
    REAL,DIMENSION(1) :: data
    CHARACTER(LEN=MAX_CHAR_LEN) :: path, name
    INTEGER(SIZE_T )  :: size = 1
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    INTENT(IN)        :: key, val, Mesh
    !------------------------------------------------------------------------!
    CALL SplitPath(key, path, name)
    data(1) = val
    SELECT CASE (SELECTED_REAL_KIND(PRECISION(val)))
    CASE(4)
        !CALL h5ltset_attribute_float_f(this%fid, path, name, data, size, this%error)
    CASE(8)
        CALL h5ltset_attribute_double_f(this%fid, path, name, data, size, this%error)
    CASE DEFAULT
        CALL Error(this,"WriteNode_hdf5_2","Cannot determine real type.")
    END SELECT
  END SUBROUTINE WriteNode_hdf5_2

  SUBROUTINE WriteNode_hdf5_3(this,Mesh,key,val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    CHARACTER(LEN=*)  :: key
    CHARACTER(LEN=*)  :: val
    !------------------------------------------------------------------------!
    CHARACTER(LEN=MAX_CHAR_LEN) :: data
    CHARACTER(LEN=MAX_CHAR_LEN) :: path, name
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    INTENT(IN)        :: key, val, Mesh
    !------------------------------------------------------------------------!
    CALL SplitPath(key, path, name)
    data = val
    CALL h5ltset_attribute_string_f(this%fid, path, name, data, this%error)
  END SUBROUTINE WriteNode_hdf5_3

  SUBROUTINE WriteNode_hdf5_4(this,Mesh,key,val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    CHARACTER(LEN=*)  :: key
    LOGICAL           :: val
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    INTENT(IN)        :: key, val, Mesh
    !------------------------------------------------------------------------!
    CALL Error(this, "WriteNode_hdf5_4", "Writing of logical types not "&
        //"implemented.")
  END SUBROUTINE WriteNode_hdf5_4

  SUBROUTINE TouchDataset_hdf5(this,Mesh,key,dsetid,filespace,memspace,rank, &
                               cdims)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    CHARACTER(LEN=*)  :: key
    INTEGER           :: rank
    INTEGER(HID_T)    :: dsetid, filespace, memspace
    INTEGER(HSIZE_T),DIMENSION(:) &
                      :: cdims
    !------------------------------------------------------------------------!
    LOGICAL           :: exists
    INTEGER(HID_T)    :: plistid
    INTEGER(HSIZE_T),DIMENSION(rank) &
                      :: offset,count,stride,block,dims
    INTEGER(HSIZE_T)  :: nblocks
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this, dsetid, filespace, memspace, cdims
    INTENT(IN)        :: Mesh, key, rank
    !------------------------------------------------------------------------!
    dims = cdims
    ! standard values (used for 1D)
    offset(:) = 0
    count = cdims
    block(:) = 1
    stride(:) = 1


    ! modify for rank>1 for parallel write
    IF(rank.GT.1) THEN
        IF(.NOT.((cdims(1).EQ.(Mesh%IMAX-Mesh%IMIN+1+2*Mesh%GNUM)).AND.&
                 (cdims(2).EQ.(Mesh%JMAX-Mesh%JMIN+1+2*Mesh%GNUM)))) THEN
            !print *,cdims
            !CALL Error(this, "TouchDataset_hdf5", "First and second dimension "//&
            !                 "have to be in grid dimensions (i,j)")
            dims(:) = cdims(:)
            offset(:) = 0
            count(:) = 1
            block(:) = cdims(:)
        ELSE
            cdims(1) = Mesh%IMAX-Mesh%IMIN+1
            cdims(2) = Mesh%JMAX-Mesh%JMIN+1

            dims(1) = Mesh%INUM
            dims(2) = Mesh%JNUM
            ! e.g. dims = (/ Mesh%INUM, Mesh%JNUM, cdims(3), cdims(4) /)
            offset(1) = Mesh%IMIN-1
            offset(2) = Mesh%JMIN-1
            ! e.g. offset = (/ Mesh%IMIN-1, Mesh%JMIN-1, 0, 0 /)
            count(1) = 1
            count(2) = 1
            ! e.g. count = (/ 1, 1, cdims(3), cdims(4) /)
            block(1) = Mesh%IMAX-Mesh%IMIN+1
            block(2) = Mesh%JMAX-Mesh%JMIN+1
            ! e.g. block = (/ Mesh%IMAX-Mesh%IMIN+1, Mesh%JMAX-Mesh%JMIN+1, 1, 1 /)
        END IF
    END IF

    !print *,GetRank(this)," dims: ",dims
    !print *,GetRank(this)," offset: ", offset
    !print *,GetRank(this)," count: ", count
    !print *,GetRank(this)," block: ", block

    CALL h5screate_simple_f(rank, cdims, memspace, this%error)

    ! Does the Dataset already exist?
    CALL h5lexists_f(this%fid, key, exists, this%error)
    IF(exists) THEN
        CALL h5dopen_f(this%fid, key, dsetid, this%error)
    ELSE
        ! Create Dataspace
        CALL h5screate_simple_f(rank, dims, filespace, this%error)

        ! Create chunked Dataset with default properties
        CALL h5pcreate_f(H5P_DATASET_CREATE_F, plistid, this%error)
        CALL h5pset_chunk_f(plistid, rank, cdims, this%error)
        CALL h5dcreate_f(this%fid, key, DEFAULT_REAL, filespace, &
                        dsetid, this%error, plistid)
        CALL h5pclose_f(plistid, this%error)
        CALL h5sclose_f(filespace, this%error)
    END IF
       
    
    CALL h5dget_space_f(dsetid, filespace, this%error)

    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
                               offset, count, &
                               this%error, &
    !                           stride = stride, &
                               block = block)
    CALL h5sget_select_hyper_nblocks_f(filespace,nblocks,this%error)
    !print *,GetRank(this)," nblocks: ",nblocks
    cdims = dims

  END SUBROUTINE TouchDataset_hdf5

  SUBROUTINE CloseDataset_hdf5(this, dsetid, filespace, memspace)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    INTEGER(HID_T)    :: dsetid, filespace, memspace
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this, dsetid, filespace, memspace
    !------------------------------------------------------------------------!
    
    CALL h5sclose_f(filespace, this%error)
    CALL h5sclose_f(memspace, this%error)
    CALL h5dclose_f(dsetid, this%error)

  END SUBROUTINE CloseDataset_hdf5


  SUBROUTINE WriteNode_hdf5_5(this,Mesh,key,val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    CHARACTER(LEN=*)  :: key
    REAL,DIMENSION(:) :: val
    !------------------------------------------------------------------------!
    INTEGER,PARAMETER :: rank = 1
    INTEGER(HSIZE_T),DIMENSION(rank) &
                      :: cdims
    INTEGER(HID_T)    :: filespace, memspace, dsetid
    LOGICAL           :: exists
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    INTENT(IN)        :: key, val, Mesh
    !------------------------------------------------------------------------!
    cdims = SHAPE(val)
    CALL TouchDataset_hdf5(this, Mesh, key, dsetid, filespace, memspace, &
                           rank, cdims)

    ! Write collectively
    CALL h5dwrite_f(dsetid, DEFAULT_REAL, val, cdims, this%error, &
                    file_space_id = filespace, &
                    mem_space_id = memspace, &
                    xfer_prp = this%xferid)

    CALL CloseDataset_hdf5(this, dsetid, filespace, memspace)

!    Simple way, but no xfer_prp == no parallel write possible    
!    CALL h5ltmake_dataset_f(this%fid, key, 1, dims, DEFAULT_REAL, &
!        val, this%error)
  END SUBROUTINE WriteNode_hdf5_5

  SUBROUTINE WriteNode_hdf5_6(this,Mesh,key)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    CHARACTER(LEN=*)  :: key
    !------------------------------------------------------------------------!
    INTEGER           :: grpid
    LOGICAL           :: link_exists
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    INTENT(IN)        :: key, Mesh
    !------------------------------------------------------------------------!
    CALL h5lexists_f(this%fid, key, link_exists, this%error)
    ! Only create group if not already existing
    IF(.NOT.link_exists) THEN
        CALL h5gcreate_f(this%fid, key, grpid, this%error)
        CALL h5gclose_f(grpid, this%error)
    END IF
  END SUBROUTINE WriteNode_hdf5_6

  SUBROUTINE WriteNode_hdf5_7(this,Mesh,key,val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    CHARACTER(LEN=*)  :: key
    REAL,DIMENSION(:,:),POINTER :: val
    !------------------------------------------------------------------------!
    INTEGER,PARAMETER :: rank = 2
    INTEGER(HSIZE_T),DIMENSION(rank) &
                      :: cdims, dims
    INTEGER(HID_T)    :: filespace, memspace, dsetid
    LOGICAL           :: exists
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    INTENT(IN)        :: Mesh, key
    !------------------------------------------------------------------------!
    cdims = SHAPE(val)
    dims = cdims
    CALL TouchDataset_hdf5(this, Mesh, key, dsetid, filespace, memspace, &
                           rank, cdims)
    ! Write collectively
    IF(.NOT.((dims(1).EQ.(Mesh%IMAX-Mesh%IMIN+1+2*Mesh%GNUM)).AND.&
             (dims(2).EQ.(Mesh%JMAX-Mesh%JMIN+1+2*Mesh%GNUM)))) THEN
      !print *,dims
      !print *,val
      CALL h5dwrite_f(dsetid, DEFAULT_REAL, &
                      val, & 
                      cdims, this%error, &
                      file_space_id = filespace, &
                      mem_space_id = memspace)
    ELSE
      CALL h5dwrite_f(dsetid, DEFAULT_REAL, &
                      val(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX), & 
                      cdims, this%error, &
                      file_space_id = filespace, &
                      mem_space_id = memspace, &
                      xfer_prp = this%xferid)
    END IF
    CALL CloseDataset_hdf5(this, dsetid, filespace, memspace)
  END SUBROUTINE WriteNode_hdf5_7

  SUBROUTINE WriteNode_hdf5_8(this,Mesh,key,val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    CHARACTER(LEN=*)  :: key
    REAL,DIMENSION(:,:,:),POINTER :: val
    !------------------------------------------------------------------------!
    INTEGER,PARAMETER :: rank = 3
    INTEGER(HSIZE_T),DIMENSION(rank) &
                      :: cdims
    INTEGER(HID_T)    :: filespace, memspace, dsetid
    LOGICAL           :: exists
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    INTENT(IN)        :: Mesh, key
    !------------------------------------------------------------------------!
    cdims = SHAPE(val)
    CALL TouchDataset_hdf5(this, Mesh, key, dsetid, filespace, memspace, &
                           rank, cdims)
    !print *,GetRank(this)," Write dims: ", cdims
    ! Write collectively
    CALL h5dwrite_f(dsetid, DEFAULT_REAL, &
                    val(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:), &
                    cdims, this%error, &
                    file_space_id = filespace, &
                    mem_space_id = memspace, &
                    xfer_prp = this%xferid)
    !print *,"Greetings from rank ",GetRank(this)

    CALL CloseDataset_hdf5(this, dsetid, filespace, memspace)
  END SUBROUTINE WriteNode_hdf5_8

  SUBROUTINE WriteNode_hdf5_9(this,Mesh,key,val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    CHARACTER(LEN=*)  :: key
    REAL,DIMENSION(:,:,:,:),POINTER :: val
    !------------------------------------------------------------------------!
    INTEGER,PARAMETER :: rank = 4
    INTEGER(HSIZE_T),DIMENSION(rank) &
                      :: cdims
    INTEGER(HID_T)    :: filespace, memspace, dsetid
    LOGICAL           :: exists
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    INTENT(IN)        :: Mesh, key
    !------------------------------------------------------------------------!
    cdims = SHAPE(val)
    CALL TouchDataset_hdf5(this, Mesh, key, dsetid, filespace, memspace, &
                           rank, cdims)

    ! Write collectively
    CALL h5dwrite_f(dsetid, DEFAULT_REAL, &
                    val(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:,:), &
                    cdims, this%error, &
                    file_space_id = filespace, &
                    mem_space_id = memspace, &
                    xfer_prp = this%xferid)

    CALL CloseDataset_hdf5(this, dsetid, filespace, memspace)
  END SUBROUTINE WriteNode_hdf5_9

  SUBROUTINE WriteNode_hdf5_10(this,Mesh,key,val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    CHARACTER(LEN=*)  :: key
    INTEGER,DIMENSION(:) :: val
    !------------------------------------------------------------------------!
    INTEGER,PARAMETER :: rank = 1
    INTEGER(HSIZE_T),DIMENSION(rank) &
                      :: cdims
    INTEGER(HID_T)    :: filespace, memspace, dsetid
    LOGICAL           :: exists
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    INTENT(IN)        :: Mesh, key, val
    !------------------------------------------------------------------------!
    cdims = SHAPE(val)
    CALL TouchDataset_hdf5(this, Mesh, key, dsetid, filespace, memspace, &
                           rank, cdims)

    ! Write collectively
    CALL h5dwrite_f(dsetid, DEFAULT_INT, val, cdims, this%error, &
                    file_space_id = filespace, &
                    mem_space_id = memspace, &
                    xfer_prp = this%xferid)

    CALL CloseDataset_hdf5(this, dsetid, filespace, memspace)
  END SUBROUTINE WriteNode_hdf5_10

  FUNCTION GetConfig_hdf5(this, filename) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    CHARACTER(LEN=*)  :: filename
    TYPE(Dict_TYP),POINTER :: res
    !------------------------------------------------------------------------!
    INTEGER           :: k
    CHARACTER(LEN=MAX_CHAR_LEN) :: str
    LOGICAL           :: success 
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    INTENT(IN)        :: filename
    !------------------------------------------------------------------------!
    CALL h5open_f(this%error)
    NULLIFY(res)
    str = filename
    k = SCAN(str, ".", .TRUE.)
    IF(k.NE.0) THEN
        this%extension = str(k+1:)
        str = str(1:k-1)
    ELSE
        this%extension = "h5"
    END IF
    this%filename = TRIM(str)
    this%cycles = 0
    CALL OpenFile_hdf5(this, READONLY)
    CALL ReadHeader_hdf5(this, res, success)
    CALL CloseFile_hdf5(this)    
    CALL CloseFileIO_hdf5(this)
  END FUNCTION GetConfig_hdf5

  !> \public Closes the file I/O
  !!
  SUBROUTINE CloseFileIO_hdf5(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    CALL h5close_f(this%error)
  END SUBROUTINE CloseFileIO_hdf5
#endif


END MODULE fileio_hdf5
