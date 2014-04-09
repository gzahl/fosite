!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fileio_generic.f90                                                #
!#                                                                           #
!# Copyright (C) 2008-2014                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
!# Manuel Jung      <mjung@astrophysik.uni-kiel.de>                          #
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
!> \addtogroup fileio
!! - general parameters of fileio group as key-values
!! \key{count,INTEGER,number of output steps,1}
!! \key{filecycles,INTEGER,number of data files (=0 => one file and append data),count+1}
!! \key{stoptime,REAL,stop time for output,Timedisc%stoptime}
!! \key{fileformat,INTEGER,type of fileio}
!! \key{filepath,CHARACTER,file path,""}
!! \key{filename,CHARACTER,file name}
!! \key{dtwall,INTEGER,wall clock time between successive outputs,3600}
!! \key{unit,INTEGER,fortran unit number for I/O,lastunit+1}
#ifdef HAVE_NETCDF
!! \key{ncfmt,INTEGER,netcdf format type}
#endif
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer 
!! \author Björn Sperling
!! \author Manuel Jung
!!
!! \brief Generic file I/O module
!!
!! This module provides the generic interface routines to all file I/O
!! modules.
!!
!! \ingroup fileio
!----------------------------------------------------------------------------!
MODULE fileio_generic
  USE fileio_gnuplot, InitFileIO_common => InitFileIO, &
        CloseFileIO_common => CloseFileIO, &
        OpenFile_basic => OpenFile
  USE fileio_binary
  USE fileio_netcdf
  USE fileio_vtk
  USE fileio_npy
  USE fileio_hdf5
  USE fileio_xdmf
  USE fluxes_common, ONLY : Fluxes_TYP
  USE mesh_common, ONLY : Mesh_TYP
  USE timedisc_common, ONLY : Timedisc_TYP
  USE physics_generic, ONLY : Physics_TYP, Convert2Conservative
  USE sources_common, ONLY : Sources_TYP
  USE common_dict
  USE fluxes_generic, ONLY : GetBoundaryFlux
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! file formats
  INTEGER, PARAMETER :: BINARY  = 1
  INTEGER, PARAMETER :: GNUPLOT = 2
  INTEGER, PARAMETER :: NETCDF  = 3
  INTEGER, PARAMETER :: VTK     = 4
  INTEGER, PARAMETER :: NPY     = 5
  INTEGER, PARAMETER :: HDF     = 6
  INTEGER, PARAMETER :: XDMF    = 7
  !--------------------------------------------------------------------------!
  INTEGER, SAVE :: lastunit = 10
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       FileIO_TYP, &
       ! constants
       BINARY, GNUPLOT, NETCDF, VTK, NPY, HDF, XDMF, &
       DTCAUSE_FILEIO, &
#ifdef HAVE_NETCDF
       NF90_NOCLOBBER, NF90_SHARE, NF90_64BIT_OFFSET, &
#ifdef HAVE_HDF5
       NF90_FORMAT_CLASSIC, NF90_FORMAT_64BIT, &
       NF90_FORMAT_NETCDF4, NF90_FORMAT_NETCDF4_CLASSIC, &
       NF90_CLASSIC_MODEL, NF90_NETCDF4, NF90_HDF5, &
#endif
#endif
       ! methods
       InitFileIO, &
       WriteHeader, &
       ReadHeader, &
       WriteTimestamp, &
       ReadTimestamp, &
       WriteDataset, &
       ReadDataset, &
       AdjustTimestep, &
       CloseFileIO, &
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
 
  !> \public Generic constructor for file I/O 
  !!
  SUBROUTINE InitFileIO(this,Mesh,Physics,Timedisc,Sources,config,IO,global_config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this            !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh            !< \param [in] Mesh mesh type
    TYPE(Physics_TYP) :: Physics         !< \param [in] Physics physics type
    TYPE(Timedisc_TYP):: Timedisc        !< \param [in] Timedisc timedisc type
    TYPE(Sources_TYP),POINTER :: Sources !< \param [in] Sources sources type
    
    TYPE(Dict_TYP),POINTER :: config!< \param [in] config dict with I/O configuration
    TYPE(Dict_TYP),POINTER :: IO    !< \param [in] IO dict with pointers to I/O arrays
    !> \param [in] global_config dict with global configuration
    TYPE(Dict_TYP),POINTER,OPTIONAL :: global_config
    !------------------------------------------------------------------------!
    LOGICAL           :: success
    CHARACTER(LEN=32) :: timestamp
    INTEGER           :: i,fstatus
    INTEGER           :: count_def, fcycles_def, dtwall_def, ncfmt_def
    REAL              :: stoptime_def
    REAL              :: time,new_time
    INTEGER           :: fileformat, unit
    CHARACTER(LEN=MAX_CHAR_LEN) :: filename, fpath
    TYPE(Dict_TYP),POINTER :: oldconfig => null()
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    INTENT(INOUT)     :: this,Timedisc,Physics
    !------------------------------------------------------------------------!
    CALL RequireKey(config, "fileformat")
    fpath = ""
    CALL RequireKey(config, "filepath", fpath)
    CALL RequireKey(config, "filename")

    ! wall clock time between successive outputs
    ! this is mainly intended for log file outputs
    CALL RequireKey(config, "dtwall", 3600) ! default is one hour
    
    ! number of output steps
    CALL RequireKey(config, "count", 1)

    CALL GetAttr(config, "count", count_def)

    ! number of data files
    ! fcycles = 0     : one data file, append data
    ! fcycles = 1     : one data file, overwrite data
    ! fcycles = X > 1 : X data files
    CALL RequireKey(config, "filecycles", count_def+1)

    ! stop time for output defaults to simulation stop time
    CALL RequireKey(config, "stoptime", Timedisc%stoptime)

    CALL GetAttr(config, "fileformat", fileformat)
    CALL GetAttr(config, "filepath", fpath)
    CALL GetAttr(config, "filename", filename)
    CALL GetAttr(config, "dtwall", dtwall_def)
    CALL GetAttr(config, "filecycles", fcycles_def)
    CALL GetAttr(config, "stoptime", stoptime_def)

    CALL RequireKey(config, "unit", lastunit+1)  
    CALL GetAttr(config, "unit", unit)
    lastunit = unit

    ! initialize file object
    SELECT CASE(fileformat)
    CASE(BINARY)
       CALL InitFileIO_binary(this,Mesh,Physics,IO,fileformat,fpath,filename,stoptime_def,&
            dtwall_def,count_def,fcycles_def,unit)
    CASE(GNUPLOT)
       CALL InitFileIO_gnuplot(this,Mesh,Physics,IO,fileformat,fpath,filename,stoptime_def,&
            dtwall_def,count_def,fcycles_def,unit,global_config)
    CASE(NETCDF)
#ifdef HAVE_NETCDF
#ifdef HAVE_HDF5
       CALL RequireKey(config, "ncfmt", NF90_NETCDF4)
#else
       CALL RequireKey(config, "ncfmt", NF90_NOCLOBBER)
#endif
#else
       CALL RequireKey(config, "ncfmt")
#endif
       CALL GetAttr(config, "ncfmt", ncfmt_def)
       CALL InitFileIO_netcdf(this,Mesh,Physics,fileformat,fpath,filename,stoptime_def,&
            dtwall_def,count_def,fcycles_def,ncfmt_def,unit)
       !TODO: support of custom output fields      
       CALL Warning(this,"InitFileIO_netcdf", "No support of custom output fields")
    CASE(VTK)
       CALL InitFileIO_vtk(this,Mesh,Physics,IO,fileformat,fpath,filename,stoptime_def,&
            dtwall_def,count_def,fcycles_def,unit)
    CASE(NPY)
       CALL InitFileIO_npy(this,Mesh,Physics,fileformat,fpath,filename,stoptime_def,&
            dtwall_def,count_def,fcycles_def,unit)
       !TODO: support of custom output fields      
       CALL Warning(this,"InitFileIO_npy", "No support of custom output fields")
    CASE(HDF)
       CALL InitFileIO_hdf5(this,Mesh,Physics,fileformat,fpath,filename,stoptime_def,&
            dtwall_def,count_def,fcycles_def,unit)
    CASE(XDMF)
       CALL InitFileIO_xdmf(this,Mesh,Physics,IO,fileformat,fpath,filename,stoptime_def,&
            dtwall_def,count_def,fcycles_def,unit)
    CASE DEFAULT
       CALL Error(this,"InitFileIO","Unknown file format.")
    END SELECT

    SELECT CASE(fileformat)
    CASE(XDMF)
      success = .TRUE.
    CASE DEFAULT
      ! find the most recent data file
      DO i=0,this%cycles
         this%step = i
         fstatus = GetFilestatus(this)
         IF (IAND(fstatus,FILE_EXISTS).GT.0) THEN
            CALL ReadTimestamp(this,time)
            IF (time.GE.this%time) THEN
               this%time = time
               CYCLE
            END IF
         END IF
         this%step = MAX(0,i-1)
         EXIT
      END DO
  
      ! read and check file header
      success = .FALSE.
      fstatus = GetFilestatus(this)
      oldconfig => Dict("leer" / "leer")
      IF (IAND(fstatus,FILE_EXISTS).GT.0) CALL ReadHeader(this,Mesh,Physics,oldconfig,success)
  
! \todo compare oldconfig and global_config => maybe success = .FALSE.
   
      ! read the data if the file is ok and the data is newer
      IF (success.AND.(this%time.GE.Timedisc%time)) THEN
         CALL ReadDataset(this,Mesh,Physics,Timedisc,IO)
         ! set new simulation time
         Timedisc%time  = this%time
         Timedisc%dtmin = Timedisc%stoptime-Timedisc%time
      ELSE
         success = .FALSE.
      END IF
    END SELECT
      
    ! compute the (actual) output time
    time = ABS(this%stoptime) / this%count
    !print *,"timetimetime ", this%stoptime, this%count, time
    this%time = time*FLOOR(Timedisc%time/time)

    ! print some information
    CALL Info(this," FILEIO---> file name:         " // TRIM(GetFilename(this)))
    IF (success) THEN
       WRITE (timestamp,'(ES10.4)') Timedisc%time
       CALL Info(this,"            time stamp:        " // TRIM(timestamp))
    END IF
    ! time for next output
    IF (Timedisc%time.GT.0.0) CALL IncTime(this)
    IF (ASSOCIATED(oldconfig)) CALL DeleteDict(oldconfig)
  END SUBROUTINE InitFileIO
 

  !> \public Generic routine to open a file
  !!
  SUBROUTINE OpenFile(this,action)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this      !< \param [in,out] this fileio type
    INTEGER          :: action    !< \param [in] action 
    !------------------------------------------------------------------------!
    INTENT(IN)       :: action
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    IF (.NOT.Initialized(this)) &
         CALL Error(this,"OpenFile","file module uninitialized")
    SELECT CASE(GetType(this))
    CASE(BINARY)
       CALL OpenFile_binary(this,action)
    CASE(GNUPLOT)
       CALL OpenFile_gnuplot(this,action)
#ifdef HAVE_NETCDF
    CASE(NETCDF)
       CALL OpenFile_netcdf(this,action)
#endif
    CASE(VTK)
       CALL OpenFile_vtk(this,action)
    CASE(NPY)
       CALL OpenFile_npy(this,action)
#ifdef HAVE_HDF5_MOD
    CASE(HDF)
       CALL OpenFile_hdf5(this,action)
#endif
    CASE(XDMF)
       CALL OpenFile_xdmf(this,action)
    END SELECT
  END SUBROUTINE OpenFile

  !> \public Generic routine to close the a file
  !!
  SUBROUTINE CloseFile(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this    !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    IF (.NOT.Initialized(this)) &
         CALL Error(this,"CloseFile","file module uninitialized")
    SELECT CASE(GetType(this))
    CASE(BINARY)
       CALL CloseFile_binary(this)
    CASE(GNUPLOT)
       CALL CloseFile_gnuplot(this)
#ifdef HAVE_NETCDF
    CASE(NETCDF)
       CALL CloseFile_netcdf(this)
#endif
    CASE(VTK)
       CALL CloseFile_vtk(this)
    CASE(NPY)
       CALL CloseFile_npy(this)
#ifdef HAVE_HDF5_MOD
    CASE(HDF)
       CALL CloseFile_hdf5(this)
#endif
    CASE(XDMF)
       CALL CloseFile_xdmf(this)
    END SELECT
  END SUBROUTINE CloseFile

  !> \public Generic routine to write a header to a file
  !!
  SUBROUTINE WriteHeader(this,Mesh,Physics,Header,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this      !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh      !< \param [in] Mesh mesh type
    TYPE(Physics_TYP) :: Physics   !< \param [in] Physics physics type
    !> \param [in] Header dict with header
    TYPE(Dict_TYP),POINTER :: Header,&
                              IO !< \param [in] IO dict with I/O arrays
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! create new empty data file
    CALL OpenFile(this,REPLACE)
    SELECT CASE(GetType(this))
    CASE(BINARY)
       CALL WriteHeader_binary(this)
    CASE(GNUPLOT)
       CALL WriteHeader_gnuplot(this)
#ifdef HAVE_NETCDF
    CASE(NETCDF)
       CALL WriteHeader_netcdf(this,Mesh,Physics)
#endif
    CASE(VTK)
       CALL WriteHeader_vtk(this,Mesh,Physics,Header,IO)
    CASE(NPY)
       CALL WriteHeader_npy(this)
#ifdef HAVE_HDF5_MOD
    CASE(HDF)
       CALL WriteHeader_hdf5(this,Mesh,Header)
#endif
    CASE(XDMF)
       CALL WriteHeader_xdmf(this)
    END SELECT
    CALL CloseFile(this)
  END SUBROUTINE WriteHeader

  !> \public Generic routine to read a header from a file
  !!
  SUBROUTINE ReadHeader(this,Mesh,Physics,Header,success)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this      !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh      !< \param [in] Mesh mesh type
    TYPE(Physics_TYP) :: Physics   !< \param [in] Physics physics type
    !> \param [in] Header dict with header
    TYPE(Dict_TYP),POINTER :: Header
    LOGICAL           :: success   !< \param [out] success
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(OUT)       :: success
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL OpenFile(this,READONLY)
    SELECT CASE(GetType(this))
    CASE(BINARY)
       CALL ReadHeader_binary(this,success)
    CASE(GNUPLOT)
       CALL ReadHeader_gnuplot(this,success)
#ifdef HAVE_NETCDF
    CASE(NETCDF)
       CALL ReadHeader_netcdf(this,Mesh,Physics,success)
#endif
    CASE(VTK)
       CALL ReadHeader_vtk(this,success)
    CASE(NPY)
       CALL ReadHeader_npy(this,success)
#ifdef HAVE_HDF5_MOD
    CASE(HDF)
       CALL ReadHeader_hdf5(this,Header,success)
#endif
    CASE(XDMF)
       CALL ReadHeader_xdmf(this,success)
    END SELECT
    CALL CloseFile(this)
  END SUBROUTINE ReadHeader

  !> \public Generic routine to write a timestamp to a file
  !!
  SUBROUTINE WriteTimestamp(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this    !< \param [in,out] this fileio type
    REAL             :: time    !< \param [in] time
    !------------------------------------------------------------------------!
    INTENT(IN)       :: time
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    CALL OpenFile(this,APPEND)
    SELECT CASE(GetType(this))
    CASE(BINARY)
       CALL WriteTimestamp_binary(this,time)
    CASE(GNUPLOT)
       CALL WriteTimestamp_gnuplot(this,time)
#ifdef HAVE_NETCDF
    CASE(NETCDF)
       CALL WriteTimestamp_netcdf(this,time)
#endif
    CASE(VTK)
       CALL WriteTimestamp_vtk(this,time)
    CASE(NPY)
       CALL WriteTimestamp_npy(this,time)
#ifdef HAVE_HDF5_MOD
    CASE(HDF)
       CALL WriteTimestamp_hdf5(this,time)
#endif
    CASE(XDMF)
       CALL WriteTimestamp_xdmf(this,time)
    END SELECT
    CALL CloseFile(this)
  END SUBROUTINE WriteTimestamp

  !> \public Generic routine to read a timestamp from a file
  !!
  SUBROUTINE ReadTimestamp(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this    !< \param [in,out] this fileio type
    REAL             :: time    !< \param [out] time
    !------------------------------------------------------------------------!
    INTENT(OUT)      :: time
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    ! open datafile in r/w mode and seek to the end
    CALL OpenFile(this,READEND)
    SELECT CASE(GetType(this))
    CASE(BINARY)
       CALL ReadTimestamp_binary(this,time)
    CASE(GNUPLOT)
       CALL ReadTimestamp_gnuplot(this,time)
#ifdef HAVE_NETCDF
    CASE(NETCDF)
       CALL ReadTimestamp_netcdf(this,time)
#endif
    CASE(VTK)
       CALL ReadTimestamp_vtk(this,time)
    CASE(NPY)
       CALL ReadTimestamp_npy(this,time)
#ifdef HAVE_HDF5_MOD
    CASE(HDF)
       CALL ReadTimestamp_hdf5(this,time)
#endif
    CASE(XDMF)
       CALL ReadTimestamp_xdmf(this,time)
    END SELECT
    CALL CloseFile(this)
  END SUBROUTINE ReadTimestamp

  !> \public Generic routine to write a dataset to a file
  !!
  SUBROUTINE WriteDataset(this,Mesh,Physics,Fluxes,Timedisc,Header,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this      !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh      !< \param [in] Mesh mesh type
    TYPE(Physics_TYP) :: Physics   !< \param [in] Physics physics type
    TYPE(Fluxes_TYP)  :: Fluxes    !< \param [in] Fluxes fluxes type
    TYPE(Timedisc_TYP):: Timedisc  !< \param [in] Timedisc timedisc type
    !> \param [in] Header dict with header
    TYPE(Dict_TYP),POINTER :: Header,&
                              IO !< \param [in] IO dict with I/O arrays
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Fluxes
    INTENT(INOUT)     :: this,Timedisc
    !------------------------------------------------------------------------!
    INTEGER           :: k, ierror
    !------------------------------------------------------------------------!
    ! calculate boundary fluxes, if they were requested for write
    IF(ASSOCIATED(Timedisc%bflux)) THEN
        DO k=1,4
            Timedisc%bflux(:,k) = GetBoundaryFlux(Fluxes,Mesh,Physics,k)
        END DO
#ifdef PARALLEL        
        CALL MPI_BCAST(Timedisc%bflux, Physics%VNUM*4, DEFAULT_MPI_REAL, 0, &
                Mesh%comm_cart, ierror)
#endif
    END IF
    ! write the header if either this is the first data set we write or
    ! each data set is written into a new file
    IF ((this%step.EQ.0).OR.(this%cycles.GT.0)) THEN
       CALL WriteHeader(this,Mesh,Physics,Header,IO)
    END IF
    CALL OpenFile(this,APPEND)
    SELECT CASE(GetType(this))
    CASE(BINARY)
       CALL WriteDataset_binary(this,Mesh,Physics,Fluxes,Timedisc)
    CASE(GNUPLOT)
       CALL WriteDataset_gnuplot(this,Mesh)
#ifdef HAVE_NETCDF
    CASE(NETCDF)
       CALL WriteDataset_netcdf(this,Mesh,Physics,Timedisc)
#endif
    CASE(VTK)
       CALL WriteDataset_vtk(this,Mesh,Physics,Fluxes,Timedisc,IO)
    CASE(NPY)
       CALL WriteDataset_npy(this,Mesh,Physics,Timedisc)
#ifdef HAVE_HDF5_MOD
    CASE(HDF)
       CALL WriteDataset_hdf5(this,Mesh,Physics,Timedisc,Fluxes,IO)
#endif
    CASE(XDMF)
       CALL WriteDataset_xdmf(this,Mesh,Physics,Fluxes,Timedisc,IO)
    END SELECT
    CALL CloseFile(this)
    ! append the time stamp
    CALL WriteTimestamp(this,Timedisc%time)
    CALL IncTime(this)
  END SUBROUTINE WriteDataset

  !> \public Generic routine to read a dataset from a file
  !!
  SUBROUTINE ReadDataset(this,Mesh,Physics,Timedisc,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this      !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh      !< \param [in] Mesh mesh type
    TYPE(Physics_TYP) :: Physics   !< \param [in] Physics physics type
    TYPE(Timedisc_TYP):: Timedisc  !< \param [in] Timedisc timedisc type
    TYPE(Dict_TYP),POINTER :: IO   !< \param [in] IO dict with I/O arrays
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    INTENT(INOUT)     :: this,Timedisc,Physics
    !------------------------------------------------------------------------!
    CALL OpenFile(this,APPEND)
    SELECT CASE(GetType(this))
    CASE(BINARY)
       CALL ReadDataset_binary(this,Mesh,Physics,Timedisc)
    CASE(GNUPLOT)
       CALL ReadDataset_gnuplot(this,Mesh,Physics,Timedisc)
#ifdef HAVE_NETCDF
    CASE(NETCDF)
       CALL ReadDataset_netcdf(this,Mesh,Physics,Timedisc)
#endif
    CASE(VTK)
       CALL Error(this,"ReadDataset","function not supported")
    CASE(NPY)
       CALL Error(this,"ReadDataset","function not supported")
#ifdef HAVE_HDF5_MOD
    CASE(HDF)
       CALL ReadDataset_hdf5(this,Mesh,Physics,Timedisc,IO)
#endif
    CASE(XDMF)
       CALL ReadDataset_xdmf(this,Mesh,Physics,Timedisc)
    END SELECT
    CALL CloseFile(this)
    ! calculate conservative variables
    CALL Convert2Conservative(Physics,Mesh,Mesh%IMIN,Mesh%IMAX,Mesh%JMIN,Mesh%JMAX, &
         Timedisc%pvar,Timedisc%cvar)
  END SUBROUTINE ReadDataset

  !> \public Generic deconstructor of the file I/O
  !!
  SUBROUTINE CloseFileIO(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this      !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    IF (.NOT.Initialized(this)) &
        CALL Error(this,"CloseFileIO","not initialized")
    SELECT CASE(GetType(this))
    CASE(BINARY)
       CALL CloseFileIO_binary(this)
    CASE(GNUPLOT)
       CALL CloseFileIO_gnuplot(this)
#ifdef HAVE_NETCDF
    CASE(NETCDF)
       CALL CloseFileIO_netcdf(this)
#endif
    CASE(VTK)
       CALL CloseFileIO_vtk(this)
    CASE(NPY)
       CALL CloseFileIO_npy(this)
#ifdef HAVE_HDF5_MOD
    CASE(HDF)
       CALL CloseFileIO_hdf5(this)
#endif
    CASE(XDMF)
       CALL CloseFileIO_xdmf(this)
    END SELECT
    CALL CloseFileIO_common(this)    
  END SUBROUTINE CloseFileIO
    
END MODULE fileio_generic
