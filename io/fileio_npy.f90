!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fileio_npy.f90                                                    #
!#                                                                           #
!# Copyright (C) 2011                                                        #
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
!! \brief I/O for numpy *.npy binary file output
!!
!! \extends fileio_common
!! \ingroup fileio
!----------------------------------------------------------------------------!
MODULE fileio_npy
  USE fileio_common
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  USE timedisc_common, ONLY : Timedisc_TYP
  USE sources_common, ONLY : Sources_TYP
#ifdef HAVE_NPY
  USE fnpy
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       FileIO_TYP, &
       ! constants
       ! methods
       InitFileIO_npy, &
       CloseFileIO_npy, &
       WriteHeader_npy, &
       ReadHeader_npy, &
       WriteTimestamp_npy,&
       ReadTimestamp_npy,&
       WriteDataset_npy, &
       ReadDataset_npy, &
       OpenFile_npy, &
       CloseFile_npy
  !--------------------------------------------------------------------------!

CONTAINS
  !> \public Constructor for the numpy binary file I/O 
  !!
  !! Initilizes the file I/O type, filename, stoptime, number of outputs, 
  !! number of files, unit number, config as a dict
  SUBROUTINE InitFileIO_npy(this,Mesh,Physics,fmt,fpath,filename,stoptime,dtwall,&
       count,fcycles,unit)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this          !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh          !< \param [in] Mesh mesh type
    TYPE(Physics_TYP) :: Physics       !< \param [in] Physics Physics type
    INTEGER           :: fmt           !< \param [in] fmt fileio type number
    CHARACTER(LEN=*)  :: fpath         !< \param [in] fpath
    CHARACTER(LEN=*)  :: filename      !< \param [in] filename
    REAL              :: stoptime      !< \param [in] stoptime
    INTEGER           :: dtwall        !< \param [in] dtwall wall clock time
    INTEGER           :: count         !< \param [in] count number of outputs
    INTEGER           :: fcycles       !< \param [in] fcycles file cycle number
    INTEGER, OPTIONAL :: unit          !< \param [in] unit fileio unit number
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,fmt,fpath,filename,stoptime,dtwall,&
                     count,fcycles,unit
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!

    CALL InitFileIO(this,fmt,"NPY",fpath,filename,"npy",fcycles,.TRUE.,unit)
    this%stoptime = stoptime
    this%dtwall   = dtwall
    this%time     = 0.
    this%count    = count
    this%step     = 0

#ifndef HAVE_NPY
    CALL Error(this, "InitFileIO_npy", "configure with --with-npy=[path]")
#endif 
   END SUBROUTINE InitFileIO_npy

  !> \public Specific routine to open a file for numpy binary I/O
  !!
  SUBROUTINE OpenFile(this,action,fformat)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    INTEGER          :: action
    CHARACTER(LEN=*) :: fformat
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    INTENT(IN)       :: action,fformat
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    SELECT CASE(action)
    CASE(READONLY)
    CASE(READEND)
    CASE(REPLACE)
    CASE(APPEND)
    CASE DEFAULT
       CALL ERROR(this,"OpenFile","Unknown access mode.")
    END SELECT
  END SUBROUTINE OpenFile

  SUBROUTINE OpenFile_npy(this,action)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    INTEGER          :: action
    !------------------------------------------------------------------------!
    INTENT(IN)       :: action
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
  END SUBROUTINE OpenFile_npy


  SUBROUTINE CloseFile_npy(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
  END SUBROUTINE CloseFile_npy
  
  SUBROUTINE WriteHeader_npy(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
  END SUBROUTINE WriteHeader_npy


  SUBROUTINE ReadHeader_npy(this,success)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    LOGICAL          :: success
    !------------------------------------------------------------------------!
    INTENT(OUT)      :: success
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    success = .FALSE.
  END SUBROUTINE ReadHeader_npy


  SUBROUTINE WriteTimestamp_npy(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    REAL             :: time
    !------------------------------------------------------------------------!
    INTENT(IN)       :: time
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
  END SUBROUTINE WriteTimestamp_npy


  SUBROUTINE ReadTimestamp_npy(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    REAL             :: time
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: time
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    time = 0.0
  END SUBROUTINE ReadTimestamp_npy


  SUBROUTINE SaveDouble2D_npy(filename, data)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CHARACTER(LEN=*)  :: filename
    REAL, DIMENSION(:,:) &
                      :: data
    !------------------------------------------------------------------------!
    INTENT(IN)        :: filename, data
    !------------------------------------------------------------------------!
#ifdef HAVE_NPY    
    CALL save_double(filename, SHAPE(data), data)
#endif
    END SUBROUTINE SaveDouble2D_npy

  SUBROUTINE SaveDouble3D_npy(filename, data)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CHARACTER(LEN=*)  :: filename
    REAL, DIMENSION(:,:,:) &
                      :: data
    !------------------------------------------------------------------------!
    INTENT(IN)        :: filename, data
    !------------------------------------------------------------------------!
    
#ifdef HAVE_NPY
    CALL save_double(filename, SHAPE(data), data)
#endif
    END SUBROUTINE SaveDouble3D_npy

  !> \public Writes all desired data arrays to a file 
  !!
  SUBROUTINE WriteDataset_npy(this,Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this      !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh      !< \param [in] mesh mesh type
    TYPE(Physics_TYP) :: Physics   !< \param [in] physics physics type
    TYPE(Timedisc_TYP):: Timedisc  !< \param [in,out] timedisc timedisc type
    !------------------------------------------------------------------------!
    !CHARACTER(LEN=*)  :: filename
    INTEGER           :: IMIN,IMAX,JMIN,JMAX
    !------------------------------------------------------------------------!
    INTENT(IN)       :: Mesh,Physics
    INTENT(INOUT)    :: this,Timedisc
    !------------------------------------------------------------------------!
    ! trim the data for output
!    WHERE (ABS(Timedisc%pvar(:,:,:)).LT.(MAX(TINY(Timedisc%pvar),1.0D-99)))
!       Timedisc%pvar(:,:,:) = 0.0E+00
!    END WHERE

    !filename = GetFilename(this)

!    IMIN = Mesh%IGMIN
!    IMAX = Mesh%IGMAX
!    JMIN = Mesh%JGMIN
!    JMAX = Mesh%JGMAX

    IMIN = Mesh%IMIN
    IMAX = Mesh%IMAX
    JMIN = Mesh%JMIN
    JMAX = Mesh%JMAX
    
    CALL SaveDouble3D_npy(TRIM(GetFilename(this)) // "_pvar.npy", &
                          Timedisc%pvar(IMIN:IMAX, &
                                        JMIN:JMAX, &
                                        :))
    !CALL SaveDouble3D_npy(TRIM(GetFilename(this)) // "_accel.npy", &
    !                      Physics%sources%accel(IMIN:IMAX, &
    !                                            JMIN:JMAX, &
    !                                            1:2))
!FIXME:
!    IF(ASSOCIATED(Physics%sources)) THEN
!        IF(ASSOCIATED(Physics%sources%poisson%phi)) THEN
!            CALL SaveDouble2D_npy(TRIM(GetFilename(this)) // "_phi.npy", &
!                                  Physics%sources%poisson%phi(IMIN:IMAX, &
!                                                              JMIN:JMAX))
!        ENDIF
!    ENDIF
  END SUBROUTINE WriteDataset_npy


  !> \public Reads the data arrays from file (not yet implemented)
  !!
  SUBROUTINE ReadDataset_npy(this,Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this      !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh      !< \param [in] mesh mesh type
    TYPE(Physics_TYP) :: Physics   !< \param [in] physics physics type
    TYPE(Timedisc_TYP):: Timedisc  !< \param [in] timedisc timedisc type
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Timedisc
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
  END SUBROUTINE ReadDataset_npy


  SUBROUTINE Error_npy(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP), INTENT(INOUT) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    IF (Initialized(this)) &
         CALL CloseFile_npy(this)
    CALL ERROR(this,modproc,msg)
  END SUBROUTINE Error_npy

  !> \public Closes the file I/O
  !!
  SUBROUTINE CloseFileIO_npy(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    !------------------------------------------------------------------------!
  END SUBROUTINE CloseFileIO_npy

END MODULE fileio_npy
