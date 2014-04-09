!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fileio_generic.f90                                                #
!#                                                                           #
!# Copyright (C) 2008-2010                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
! generic module for file I/O 
!----------------------------------------------------------------------------!
MODULE fileio_generic
  USE fileio_gnuplot, InitFileIO_common => InitFileIO, OpenFile_basic => OpenFile
  USE fileio_binary
  USE fileio_netcdf
  USE fileio_vtk
  USE fluxes_common, ONLY : Fluxes_TYP
  USE mesh_common, ONLY : Mesh_TYP
  USE timedisc_common, ONLY : Timedisc_TYP
  USE physics_generic, ONLY : Physics_TYP, Convert2Conservative
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! file formats
  INTEGER, PARAMETER :: BINARY  = 1
  INTEGER, PARAMETER :: GNUPLOT = 2
  INTEGER, PARAMETER :: NETCDF  = 3
  INTEGER, PARAMETER :: VTK     = 4
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       FileIO_TYP, &
       ! constants
       BINARY, GNUPLOT, NETCDF, VTK, &
#ifdef HAVE_NETCDF
#ifdef HAVE_HDF5
       NF90_CLASSIC_MODEL, NF90_NETCDF4, &
#else
       NF90_FORMAT_CLASSIC, NF90_FORMAT_64BIT, &       
#endif
#endif
       ! methods
       InitFileIO, &
       WriteHeader, &
       ReadHeader, &
       WriteTimestamp, &
       ReadTimestamp, &
       WriteDataset, &
       AdjustTimestep, &
       CloseFileIO, &
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

  SUBROUTINE InitFileIO(this,Mesh,Physics,Timedisc,fileformat,filename,&
       stoptime,dtwall,count,filecycles,ncfmt,unit)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Timedisc_TYP):: Timedisc
    INTEGER           :: fileformat
    CHARACTER(LEN=*)  :: filename
    REAL, OPTIONAL    :: stoptime
    INTEGER, OPTIONAL :: dtwall
    INTEGER, OPTIONAL :: count
    INTEGER, OPTIONAL :: filecycles
    INTEGER, OPTIONAL :: ncfmt
    INTEGER, OPTIONAL :: unit
    !------------------------------------------------------------------------!
    LOGICAL           :: success
    CHARACTER(LEN=32) :: timestamp
    INTEGER           :: i,fstatus
    INTEGER           :: count_def, fcycles_def, dtwall_def, ncfmt_def
    REAL              :: stoptime_def
    REAL              :: time,new_time
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,fileformat,filename,stoptime,dtwall,&
                         count,filecycles,ncfmt,unit
    INTENT(INOUT)     :: this,Timedisc
    !------------------------------------------------------------------------!
    ! wall clock time between successive outputs
    ! this is mainly intended for log file outputs
    IF (PRESENT(dtwall)) THEN
       dtwall_def = dtwall
    ELSE
       dtwall_def = 3600 ! default is one hour
    END IF
    ! number of output steps
    IF (PRESENT(count)) THEN
       count_def = count
    ELSE
       count_def = 1
    END IF

    ! number of data files
    ! fcycles = 0     : one data file, append data
    ! fcycles = 1     : one data file, overwrite data
    ! fcycles = X > 1 : X data files
    IF (PRESENT(filecycles)) THEN
       fcycles_def = filecycles
    ELSE ! default
       fcycles_def = count_def+1
    END IF

    ! stop time for output defaults to simulation stop time
    IF (PRESENT(stoptime)) THEN
       stoptime_def = stoptime
    ELSE
       stoptime_def = Timedisc%stoptime
    END IF

    ! initialize file object
    SELECT CASE(fileformat)
    CASE(BINARY)
       CALL InitFileIO_binary(this,Mesh,Physics,fileformat,filename,stoptime_def,&
            dtwall_def,count_def,fcycles_def,unit)
    CASE(GNUPLOT)
       CALL InitFileIO_gnuplot(this,Mesh,Physics,fileformat,filename,stoptime_def,&
            dtwall_def,count_def,fcycles_def,unit)
    CASE(NETCDF)
#ifdef HAVE_NETCDF
#ifdef PARALLEL
#ifdef HAVE_HDF5
       ncfmt_def = NF90_NETCDF4
#else
       CALL Error(this,"InitFileIO","HDF5 required for parallel NetCDF i/o")
#endif
#else
       IF (PRESENT(ncfmt)) THEN
          ncfmt_def = ncfmt
       ELSE
#ifdef HAVE_HDF5
          ncfmt_def = NF90_NETCDF4
#else
          ncfmt_def = NF90_FORMAT_CLASSIC
#endif
       END IF
#endif
       CALL InitFileIO_netcdf(this,Mesh,Physics,fileformat,filename,stoptime_def,&
            dtwall_def,count_def,fcycles_def,ncfmt_def,unit)
#else
       CALL Error(this,"InitFileIO","NetCDF support disabled")
#endif
    CASE(VTK)
       CALL InitFileIO_vtk(this,Mesh,Physics,fileformat,filename,stoptime_def,&
            dtwall_def,count_def,fcycles_def,unit)
    CASE DEFAULT
       CALL Error(this,"InitFileIO","Unknown file format.")
    END SELECT

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
    IF (IAND(fstatus,FILE_EXISTS).GT.0) CALL ReadHeader(this,Mesh,Physics,success)
 
    ! read the data if the file is ok and the data is newer
    IF (success.AND.(this%time.GT.Timedisc%time)) THEN
       CALL ReadDataset(this,Mesh,Physics,Timedisc)
       ! set new simulation time
       Timedisc%time  = this%time
       Timedisc%dtmin = Timedisc%stoptime-Timedisc%time
    ELSE
       success = .FALSE.
    END IF
    
    ! compute the (actual) output time
    time = ABS(this%stoptime) / this%count
    this%time = time*FLOOR(Timedisc%time/time)

    ! print some information
    CALL Info(this," FILEIO---> file name:         " // TRIM(GetFilename(this)))
    IF (success) THEN
       WRITE (timestamp,'(ES10.4)') Timedisc%time
       CALL Info(this,"            time stamp:        " // TRIM(timestamp))
    END IF
    ! time for next output
    IF (Timedisc%time.GT.0.0) CALL IncTime(this)
  END SUBROUTINE InitFileIO

  
  SUBROUTINE OpenFile(this,action)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    INTEGER          :: action
    !------------------------------------------------------------------------!
    INTENT(IN)       :: action
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
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
    END SELECT
  END SUBROUTINE OpenFile


  SUBROUTINE CloseFile(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
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
    END SELECT
  END SUBROUTINE CloseFile


  SUBROUTINE WriteHeader(this,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
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
       CALL WriteHeader_vtk(this,Mesh,Physics)
    END SELECT
    CALL CloseFile(this)
  END SUBROUTINE WriteHeader

  
  SUBROUTINE ReadHeader(this,Mesh,Physics,success)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    LOGICAL           :: success
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(OUT)       :: success
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL OpenFile(this,READONLY)
    CALL RewindFile(this)
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
    END SELECT
    CALL CloseFile(this)
  END SUBROUTINE ReadHeader

  
  SUBROUTINE WriteTimestamp(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    REAL             :: time
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
    END SELECT
    CALL CloseFile(this)
  END SUBROUTINE WriteTimestamp

  
  SUBROUTINE ReadTimestamp(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    REAL             :: time
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
    END SELECT
    CALL CloseFile(this)
  END SUBROUTINE ReadTimestamp

  
  SUBROUTINE WriteDataset(this,Mesh,Physics,Fluxes,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Fluxes
    INTENT(INOUT)     :: this,Timedisc
    !------------------------------------------------------------------------!
    ! write the header if either this is the first data set we write or
    ! each data set is written into a new file
    IF ((this%step.EQ.0).OR.(this%cycles.GT.0)) THEN
       CALL WriteHeader(this,Mesh,Physics)
    END IF
    CALL OpenFile(this,APPEND)
    SELECT CASE(GetType(this))
    CASE(BINARY)
       CALL WriteDataset_binary(this,Mesh,Physics,Fluxes,Timedisc)
    CASE(GNUPLOT)
       CALL WriteDataset_gnuplot(this,Mesh,Physics,Timedisc)
#ifdef HAVE_NETCDF
    CASE(NETCDF)
       CALL WriteDataset_netcdf(this,Mesh,Physics,Timedisc)
#endif
    CASE(VTK)
       CALL WriteDataset_vtk(this,Mesh,Physics,Timedisc)
    END SELECT
    CALL CloseFile(this)
    ! append the time stamp
    CALL WriteTimestamp(this,Timedisc%time)
    CALL IncTime(this)
  END SUBROUTINE WriteDataset


  SUBROUTINE ReadDataset(this,Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: this,Timedisc
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
    END SELECT
    CALL CloseFile(this)
    ! calculate conservative variables
    CALL Convert2Conservative(Physics,Mesh,Mesh%IMIN,Mesh%IMAX,Mesh%JMIN,Mesh%JMAX, &
         Timedisc%pvar,Timedisc%cvar)
  END SUBROUTINE ReadDataset


  SUBROUTINE CloseFileIO(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
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
    END SELECT
  END SUBROUTINE CloseFileIO
    
END MODULE fileio_generic
