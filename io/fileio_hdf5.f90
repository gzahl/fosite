!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fileio_hdf5.f90                                                   #
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
! module for HDF5 I/O
!----------------------------------------------------------------------------!
MODULE fileio_hdf5
  USE fileio_gnuplot, CloseFile_gnuplot => CloseFile, &
    CloseFileio_hdf5 => CloseFileIO_gnuplot
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: header_bytes = 25
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       FileIO_TYP, &
       ! constants
       ! methods
       InitFileio_hdf5, &
       WriteHeader_hdf5, &
       CloseFileio_hdf5
  !--------------------------------------------------------------------------!

CONTAINS
  
  SUBROUTINE InitFileio_hdf5(this,Mesh,Physics,fmt,filename,stoptime,dtwall,&
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
    INTENT(IN)    :: Mesh,Physics,fmt,filename,stoptime,count,fcycles
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL InitFileIO(this,Mesh,Physics,fmt,"hdf5",filename,"hdf",stoptime,&
         dtwall,count,fcycles)
  END SUBROUTINE InitFileio_hdf5


  SUBROUTINE OpenFile_hdf5(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
  END SUBROUTINE OpenFile_hdf5


  SUBROUTINE WriteHeader_hdf5(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
  END SUBROUTINE WriteHeader_hdf5

END MODULE fileio_hdf5
