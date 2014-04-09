!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fileio_xdmf.f90                                                   #
!#                                                                           #
!# Copyright (C) 2013                                                        #
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
!> \author Manuel Jung
!!
!! \brief module for XDMF file I/O
!!
!! The xdmf file format carries light data in a xml file and can store heavy
!! data in binary or hdf5 files. This implementation stores heavy data in
!! binary files. To write the binary files, we need one of the following:
!! - a fortran compiler with f2003 Stream IO
!! - a MPI build
!! - on NEC sx9: Set the runtime enviroment variable F_NORCW=5555 (or to
!!   another value)
!! specifications: http://www.xdmf.org/index.php/XDMF_Model_and_Format
!!
!! \extends fileio_vtk
!! \ingroup fileio
!----------------------------------------------------------------------------!
MODULE fileio_xdmf
  USE fileio_gnuplot, CloseFile_xdmf => CloseFile_gnuplot
  USE fileio_vtk, ONLY : GetPrecision_xdmf => GetPrecision_vtk, &
      GetEndianness_xdmf => GetEndianness_vtk, OpenFile_vtk
  USE geometry_common, ONLY : Geometry_TYP, GetType
  USE mesh_common, ONLY : Mesh_TYP
  USE mesh_generic, ONLY : BIANGLESPHERICAL
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
  CHARACTER, PARAMETER    :: LF = ACHAR(10)        !< line feed
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       FileIO_TYP, &
       ! constants
       ! methods
       InitFileIO_xdmf, &
       OpenFile_xdmf, &
       CloseFile_xdmf, &
       WriteHeader_xdmf, &
       ReadHeader_xdmf, &
       WriteTimestamp_xdmf,&
       ReadTimestamp_xdmf,&
       WriteDataset_xdmf, &
       ReadDataset_xdmf, &
       CloseFileIO_xdmf
  !--------------------------------------------------------------------------!

CONTAINS
  !> \public Constructor for the xdmf file I/O 
  !!
  !! Initilizes the file I/O type, filename, stoptime, number of outputs, 
  !! number of files, unit number, config as a dict
  SUBROUTINE InitFileIO_xdmf(this,Mesh,Physics,IO,fmt,fpath,filename,stoptime,dtwall,&
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
#ifdef PARALLEL
    INTEGER, DIMENSION(2) :: gsizes,lsizes,indices,memsizes
    INTEGER           :: lb,extent
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,fmt,filename,stoptime,count,fcycles,unit
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL InitFileIO(this,Mesh,Physics,fmt,"xdmf",fpath,filename,"bin",stoptime, &
         dtwall,count,fcycles,.FALSE.,unit)

    CALL GetEndianness_xdmf(this, this%endianness, '"Little"', '"Big"')
    CALL GetPrecision_xdmf(this, this%realsize)
    IF(this%realsize.GT.8) &
      CALL Error(this,"WriteXMF_xmf","Only single and double precision are allowed")

    WRITE(this%realfmt,'(A,I1,A)',IOSTAT=this%error) '"',this%realsize,'"'

    CALL WriteXMF_xdmf(this,Mesh,IO)

#ifdef PARALLEL
    CALL MPI_Type_contiguous(1,DEFAULT_MPI_REAL,this%basictype,this%error)
    CALL MPI_Type_commit(this%basictype,this%error)

    ! create the data type for the distributed array of
    ! coordinates and simulation data
    gsizes(1) = Mesh%INUM
    gsizes(2) = Mesh%JNUM
    lsizes(1) = Mesh%IMAX-Mesh%IMIN+1
    lsizes(2) = Mesh%JMAX-Mesh%JMIN+1
    indices(1)= Mesh%IMIN-1
    indices(2)= Mesh%JMIN-1
    this%bufsize = PRODUCT(lsizes)
    CALL MPI_Type_create_subarray(2, gsizes, lsizes, indices,MPI_ORDER_FORTRAN,&
         DEFAULT_MPI_REAL,this%filetype,this%error)
    CALL MPI_Type_commit(this%filetype,this%error)

    memsizes(:) = lsizes(:) + Mesh%GNUM*2
    indices(:) = Mesh%GNUM
    CALL MPI_Type_create_subarray(2, memsizes, lsizes, indices, &
           MPI_ORDER_FORTRAN, DEFAULT_MPI_REAL, this%memtype, this%error)
    CALL MPI_Type_commit(this%memtype,this%error)

    ! create the data type for the distributed array of
    ! mesh corner positions
    gsizes(1) = Mesh%INUM+1
    gsizes(2) = Mesh%JNUM+1
    IF(Mesh%INUM.EQ.Mesh%IMAX) THEN
      lsizes(1) = Mesh%IMAX-Mesh%IMIN+2
    ELSE
      lsizes(1) = Mesh%IMAX-Mesh%IMIN+1
    END IF      
    IF(Mesh%JNUM.EQ.Mesh%JMAX) THEN
      lsizes(2) = Mesh%JMAX-Mesh%JMIN+2
    ELSE
      lsizes(2) = Mesh%JMAX-Mesh%JMIN+1
    END IF
    indices(1)= Mesh%IMIN-1
    indices(2)= Mesh%JMIN-1
    this%meshbufsize = lsizes(1) * lsizes(2)
    CALL MPI_Type_create_subarray(2, gsizes, lsizes, indices, MPI_ORDER_FORTRAN,&
         this%basictype,this%meshtype,this%error)
    CALL MPI_Type_commit(this%meshtype,this%error)
#endif

  END SUBROUTINE InitFileIO_xdmf


  !> Sets pointer to rank 2 array with a given lower bound
  !!
  !! This function sets a pointer to a rank 2 array with new
  !! lower boundaries. This is necessary because the lower bound of a array 
  !! will be set to zero when it is saved in a dictionary. 
  !! \todo reuse remap_bounds2 function from gnuplot 
  !! \return pointer to data array with new bounds
 FUNCTION remap_bounds2(lb1,lb2,array) RESULT(ptr)
  !------------------------------------------------------------------------!
    INTEGER, INTENT(IN) :: lb1      !< \param [in] lb1 new lower bound
    INTEGER, INTENT(IN) :: lb2      !< \param [in] lb2 new lower bound
    !> \param [in] array data array
    REAL, DIMENSION(lb1:,lb2:), INTENT(IN), TARGET :: array
  !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:), POINTER                  :: ptr
  !------------------------------------------------------------------------!
    ptr => array 
  END FUNCTION

  !> Writes description of data item in xml syntax
  !!
  SUBROUTINE WriteDataItem_xdmf(this,dims,filename,offset)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this           !< \param [in,out] this fileio type
    CHARACTER(LEN=*) :: dims           !< \param [in] dims 
    CHARACTER(LEN=*) :: filename       !< \param [in] filename 
    INTEGER          :: offset         !< \param [in] offset 
    !------------------------------------------------------------------------!
    CHARACTER(LEN=16):: seek
    !------------------------------------------------------------------------!
    INTENT(IN)       :: dims,filename,offset
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    WRITE(this%unit, IOSTAT=this%error)&
           '<DataItem Dimensions=' // TRIM(dims) // ' ' &
        // 'NumberType="Float" Precision=' // TRIM(this%realfmt) // LF &
        // 'Format="Binary" Endian=' // TRIM(this%endianness)
    WRITE(seek,'(I16)',IOSTAT=this%error) offset*this%realsize
    WRITE(this%unit, IOSTAT=this%error)&
       LF // 'Seek="' // TRIM(ADJUSTL(seek)) // '"'
    WRITE(this%unit, IOSTAT=this%error)&
           '>' // LF &
        // TRIM(filename) // LF &
        // '</DataItem>' // LF

  END SUBROUTINE WriteDataItem_xdmf

  !> Writes reference to an data item
  !!
  SUBROUTINE WriteReference_xdmf(this,dims,name,i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this        !< \param [in,out] this fileio type
    CHARACTER(LEN=*) :: dims        !< \param [in] dims
    CHARACTER(LEN=*) :: name        !< \param [in] name
    INTEGER,OPTIONAL :: i           !< \param [in] i
    !------------------------------------------------------------------------!
    CHARACTER(LEN=16):: no
    !------------------------------------------------------------------------!
    INTENT(IN)       :: dims,name
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    IF(PRESENT(i)) THEN
      WRITE(no,'(I4)') i
    ELSE
      no = "1"
    END IF
    WRITE(this%unit, IOSTAT=this%error)&
       '<DataItem Dimensions=' // TRIM(dims) // ' ' &
    // 'Reference="/Xdmf/Domain/Grid/Grid/Attribute[@Name=' // "'" &
    // TRIM(name) // "'" // ']/DataItem[' // TRIM(no) // ']"/>' // LF

  END SUBROUTINE WriteReference_xdmf  

  !> Writes description of data item in xml syntax
  !!
  RECURSIVE SUBROUTINE WriteAttribute_xdmf(this,Mesh,name,dims,filename,offset,ref)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this       !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)   :: Mesh       !< \param [in] mesh mesh type
    CHARACTER(LEN=*) :: name       !< \param [in] name
    CHARACTER(LEN=*) :: filename   !< \param [in] filename
    INTEGER, DIMENSION(:) :: dims  !< \param [in] dims
    INTEGER          :: offset     !< \param [in,out] offset
    LOGICAL          :: ref        !< \param [in] ref
    !------------------------------------------------------------------------!
    CHARACTER(LEN=8) :: inum,jnum
    CHARACTER(LEN=32):: type,center
    INTEGER          :: i
    !------------------------------------------------------------------------!
    INTENT(IN)       :: Mesh,name,filename
    INTENT(INOUT)    :: this,offset
    !------------------------------------------------------------------------!
    SELECT CASE(SIZE(dims))
    CASE(1)
      type = "Scalar"
      center = "Grid"
    CASE(2)
      type = "Scalar"
      center = "Cell"
    CASE(3)
      type = "Vector"
      center = "Cell"
    END SELECT 
    SELECT CASE(SIZE(dims))
    CASE(1,2)
      WRITE(this%unit, IOSTAT=this%error)&
           '<Attribute Name="' // TRIM(name) // '" '&
           // 'AttributeType="' // TRIM(type) // '" ' &
           // 'Center="' // TRIM(center) // '">' // LF
      IF(.NOT.ref) THEN
        CALL WriteDataItem_xdmf(this, GetDimsStr(Mesh,dims), TRIM(filename), offset)
        offset = offset + GetSize(Mesh,dims)
      ELSE
        CALL WriteReference_xdmf(this, GetDimsStr(Mesh,dims), TRIM(name))
      END IF
      WRITE(this%unit, IOSTAT=this%error)&
        '</Attribute>' // LF
    CASE(3)
      CALL WriteAttribute_xdmf(this,Mesh,TRIM(name)//"_x",dims(1:2),filename,offset,ref)
      CALL WriteAttribute_xdmf(this,Mesh,TRIM(name)//"_y",dims(1:2),filename,offset,ref)
    END SELECT
  END SUBROUTINE WriteAttribute_xdmf

  !> Writes the mesh to file
  !!
  SUBROUTINE WriteMeshXML_xdmf(this,Mesh,filename,offset,ref)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this        !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)   :: Mesh        !< \param [in] mesh mesh type
    CHARACTER(LEN=*) :: filename    !< \param [in] filename
    INTEGER          :: offset      !< \param [in,out] offset
    LOGICAL          :: ref         !< \param [in] ref
    !------------------------------------------------------------------------!
    CHARACTER(LEN=64):: dims
    CHARACTER(LEN=8) :: inum,jnum
    INTEGER          :: i
    !------------------------------------------------------------------------!
    INTENT(IN)       :: Mesh,filename,ref
    INTENT(INOUT)    :: this, offset
    !------------------------------------------------------------------------!
    WRITE(inum,'(I8)') Mesh%INUM+1
    WRITE(jnum,'(I8)') Mesh%JNUM+1
    WRITE(dims,'(A,A,A)') TRIM(ADJUSTL(jnum)),' ',TRIM(ADJUSTL(inum))

    WRITE(this%unit, IOSTAT=this%error)&
           '<Topology TopologyType="2DSMesh" ' &
        // 'NumberOfElements="' // TRIM(dims) // '"/>' // LF &
        // '<Geometry GeometryType="X_Y_Z">' // LF
      IF(.NOT.ref) THEN
        DO i=1,3
          CALL WriteDataItem_xdmf(this, '"'//TRIM(dims)//'"', &
                                  TRIM(filename),offset)
          offset = offset + (Mesh%INUM+1)*(Mesh%JNUM+1)
        END DO
      ELSE
        WRITE(this%unit, IOSTAT=this%error)&
              '<DataItem Reference="/Xdmf/Domain/Grid/Grid/Geometry/DataItem[1]"/>' // LF &
           // '<DataItem Reference="/Xdmf/Domain/Grid/Grid/Geometry/DataItem[2]"/>' // LF &
           // '<DataItem Reference="/Xdmf/Domain/Grid/Grid/Geometry/DataItem[3]"/>' // LF
      END IF

      WRITE(this%unit, IOSTAT=this%error)&
           '</Geometry>' // LF

  END SUBROUTINE WriteMeshXML_xdmf

  !> Main routine to write all data to xmf file
  !!
  SUBROUTINE WriteXMF_xdmf(this,Mesh,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this        !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)   :: Mesh        !< \param [in] mesh mesh type
    TYPE(Dict_TYP),POINTER :: IO    !< \param [in] io I/O dictionary
    !------------------------------------------------------------------------!
    INTEGER          :: i, offset
    REAL             :: ftime
    CHARACTER(LEN=4) :: step
    CHARACTER(LEN=32):: time
    CHARACTER(LEN=256):: filename
    TYPE(Dict_TYP),POINTER :: meshIO => Null()
    !------------------------------------------------------------------------!
    INTENT(IN)       :: Mesh
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!

    IF (GetRank(this).EQ.0) THEN
      ! write a xmf-file
      OPEN(this%unit, FILE=TRIM(this%filename)//'.xmf', &
           STATUS     = 'REPLACE',      &
           ACTION     = 'WRITE',        &
#ifdef FORTRAN_STREAMS
           ACCESS     = 'STREAM',       &
#else
           FORM       = 'UNFORMATTED',  &
#endif
           POSITION   = 'REWIND',       &
           IOSTAT     = this%error)
       IF (this%error.NE. 0) CALL Error(this,"WriteXMF_xdmf","Can't open xmf-file")

      WRITE(this%unit, IOSTAT=this%error)&
           '<?xml version="1.0" ?>' // LF &
        // '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>' // LF &
        // '<Xdmf Version="2.0">' // LF &
        // '<Domain>' // LF &
        // '<Grid Name="mesh" GridType="Collection" CollectionType="Temporal">' // LF 

! this works in paraview, but not in visit :(. So we have to use the simpler method for now
!        // '<Time TimeType="List">' // LF &
!        // '<DataItem Format="XML" NumberType="Float" Dimensions="11">' // LF &
!        // '0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0' // LF &
!        // '</DataItem>' // LF &
!        // '</Time>' // LF 
!      CALL WriteDataItem_xdmf(this, '"1"', TRIM(filename)//'_'//step//'.bin')
!      WRITE(this%unit, IOSTAT=this%error)&      

      CALL GetAttr(IO,"mesh",meshIO)

      DO i=0,this%count
        offset = 0
        ftime = this%stoptime/(this%count)*i
        WRITE(step,'(I4.4)') i
        WRITE(time,'(ES16.9)') ftime
        WRITE(filename, '(A)') TRIM(this%filename)//"_"//step//".bin"

        WRITE(this%unit, IOSTAT=this%error)&
             '<Grid Name="step' // step // '" GridType="Uniform">' // LF &
          // '<Time Value="' // TRIM(time) // '" />' // LF

        ! Write time again as Attribute
        CALL WriteAttribute_xdmf(this,Mesh,"/timedisc/time",(/1/),filename,offset,.FALSE.)

        CALL WriteAttributes_xdmf(this,Mesh,IO,filename,offset,.FALSE.)

        IF(i.EQ.0) THEN    ! Write mesh arrays only into file 0, because they are constant
          CALL WriteMeshXML_xdmf(this,Mesh,filename,offset,.FALSE.)
          CALL WriteAttributes_xdmf(this,Mesh,meshIO,filename,offset,.FALSE.,"/mesh")
        ELSE               ! Reference file 0 in other timesteps
          CALL WriteMeshXML_xdmf(this,Mesh,filename,offset,.TRUE.)
          CALL WriteAttributes_xdmf(this,Mesh,meshIO,filename,offset,.TRUE.,"/mesh")
        END IF

        !WRITE(this%unit, IOSTAT=this%error)&
        !     '<Attribute Name="/timedisc/velocity" AttributeType="Vector" Center="Node">' // LF &
        !  // '<DataItem Dimensions="384 128" Function="JOIN($0, $1, 0)" ItemType="Function">' // LF &
        !  // '<DataItem Dimensions="384 128" Reference="/Xdmf/Domain/Grid/Grid/Attribute[@Name=' &
        !    // "'" // '/timedisc/xvelocity' // "'" // ']/DataItem"/>' // LF &
        !  // '<DataItem Dimensions="384 128" Reference="/Xdmf/Domain/Grid/Grid/Attribute[@Name=' &
        !    // "'" // '/timedisc/yvelocity' // "'" // ']/DataItem"/>' // LF &
        !  // '</DataItem>' // LF &
        !  // '</Attribute>' // LF


        WRITE(this%unit, IOSTAT=this%error)&
             '</Grid>' // LF
      END DO

      WRITE(this%unit, IOSTAT=this%error)&
           '</Grid>' // LF


      WRITE(this%unit, IOSTAT=this%error)&
           '</Domain>' // LF &
        // '</Xdmf>' // LF

      IF(this%error.NE.0) CALL Error(this,"WriteXMF_xmf","Can't write xmf-file")
      CLOSE(this%unit,IOSTAT=this%error)
      IF(this%error.NE.0) CALL Error(this, "WriteXMF_xdmf", "Can't close xmf-file")
   END IF

  END SUBROUTINE WriteXMF_xdmf

  !> Writes attributes....
  !!
  RECURSIVE SUBROUTINE WriteAttributes_xdmf(this,Mesh,config,filename,offset,ref,path)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)       :: this        !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)         :: Mesh        !< \param [in] mesh mesh type
    TYPE(Dict_TYP),POINTER :: config      !< \param [in] config config dictionary
    CHARACTER(LEN=*)       :: filename    !< \param [in] filename 
    INTEGER                :: offset      !< \param [in,out] offset
    LOGICAL                :: ref         !< \param [in] ref
    !! \todo intent of path is what?
    CHARACTER(LEN=*),OPTIONAL & 
                           :: path        !< \param [in] path 

    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: dir
    CHARACTER(LEN=MAX_CHAR_LEN) :: str, key
    TYPE(dict_TYP),POINTER   :: node
    REAL,POINTER           :: realptr
    REAL,DIMENSION(:,:),POINTER :: ptr2
    REAL,DIMENSION(:,:,:),POINTER :: ptr3
    REAL,DIMENSION(:,:,:,:),POINTER :: ptr4
    INTEGER                :: dims2(2),dims3(3),dims4(4),dims1(1)
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh, filename, ref
    INTENT(INOUT)     :: this, offset
    !------------------------------------------------------------------------!
    IF(PRESENT(path)) THEN
        str = path
    ELSE
        str = ""
    ENDIF
    node => config
    DO WHILE(ASSOCIATED(node))
      key = TRIM(str)//"/"//TRIM(GetKey(node))
      SELECT CASE(GetDataType(node))
      CASE(DICT_DIR)
        IF(GetKey(node).NE."mesh") THEN
          key = TRIM(key)
          CALL GetAttr(node,GetKey(node),dir)
          CALL WriteAttributes_xdmf(this, Mesh, dir, filename, offset, ref, key)
        END IF
      CASE(DICT_REAL)
        CALL GetAttr(node,"value",realptr)
        dims1(1) = 1
        CALL WriteAttribute_xdmf(this,Mesh,TRIM(str),dims1,filename,offset,ref)
      CASE(DICT_REAL_TWOD)
        CALL GetAttr(node,"value",ptr2)
        dims2 = SHAPE(ptr2)
        IF((dims2(1).EQ.(Mesh%IGMAX-Mesh%IGMIN+1)).AND.(dims2(2).EQ.(Mesh%JGMAX-Mesh%JGMIN+1))) THEN
          CALL WriteAttribute_xdmf(this,Mesh,str,dims2,filename,offset,ref)
        ELSE
          CALL WriteAttribute_xdmf(this,Mesh,str,dims2,filename,offset,ref)
        END IF
      CASE(DICT_REAL_THREED)
        CALL GetAttr(node,"value",ptr3)
        !dims3 = SHAPE(ptr3(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:))
        dims3 = SHAPE(ptr3)
        CALL WriteAttribute_xdmf(this,Mesh,TRIM(str),dims3,filename,offset,ref)
        
        !CALL WriteAttribute_xdmf(this,Mesh,TRIM(str)//"_x","Scalar",GetDimsStr(Mesh,dims3(1:2)),filename,offset,ref)
        !offset = offset + GetSize(Mesh,dims3(1:2))
        !CALL WriteAttribute_xdmf(this,Mesh,TRIM(str)//"_y","Scalar",GetDimsStr(Mesh,dims3(1:2)),filename,offset,ref)
        !offset = offset + GetSize(Mesh,dims3(1:2))
      CASE(DICT_REAL_FOURD)
        CALL GetAttr(node,"value",ptr4)
        dims4 = SHAPE(ptr4(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:,:))
        !CALL WriteAttribute_xdmf(this,Mesh,str,"Tensor",GetDimsStr(Mesh,dims4),filename,offset,ref)
        !offset = offset + GetSize(Mesh,dims4)
      END SELECT
      node=>GetNext(node)
    END DO
  END SUBROUTINE WriteAttributes_xdmf


  FUNCTION GetSize(Mesh,dims) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)         :: Mesh
    INTEGER,DIMENSION(:)   :: dims
    INTEGER                :: res
    !------------------------------------------------------------------------!
    INTEGER                :: i, l
    !------------------------------------------------------------------------!
    INTENT(IN)             :: Mesh,dims
    !------------------------------------------------------------------------!
    l = SIZE(dims)
    IF(l.GE.2) THEN
      IF((dims(1).EQ.(Mesh%IGMAX-Mesh%IGMIN+1)).AND.(dims(2).EQ.(Mesh%JGMAX-Mesh%JGMIN+1))) THEN
        res = Mesh%INUM*Mesh%JNUM
      ELSE
        res = dims(1)*dims(2)
      END IF
      DO i=3,l
        res = res * dims(i)
      END DO
    ELSE
      res = dims(1)
    END IF
  END FUNCTION GetSize

  
  FUNCTION GetDimsStr(Mesh,dims) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)         :: Mesh
    INTEGER,DIMENSION(:)   :: dims
    CHARACTER(LEN=128)     :: res
    !------------------------------------------------------------------------!
    INTEGER                :: i, l
    CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE &
                           :: buf
    !------------------------------------------------------------------------!
    INTENT(IN)             :: Mesh,dims
    !------------------------------------------------------------------------!
    l = SIZE(dims)
    ALLOCATE(buf(l))
    IF(l.GE.2) THEN
      IF((dims(1).EQ.(Mesh%IGMAX-Mesh%IGMIN+1)).AND.(dims(2).EQ.(Mesh%JGMAX-Mesh%JGMIN+1))) THEN
        WRITE(buf(1),'(I8)') Mesh%JNUM
        WRITE(buf(2),'(I8)') Mesh%INUM
      ELSE
        WRITE(buf(1),'(I8)') dims(2)
        WRITE(buf(2),'(I8)') dims(1)
      END IF
      DO i=3,l
        WRITE(buf(i),'(I8)') dims(i)
      END DO
      SELECT CASE(l)
      CASE(2)
        WRITE(res,'(A,A,A,A,A)') '"',TRIM(ADJUSTL(buf(1))),&
          ' ',TRIM(ADJUSTL(buf(2))),'"'
      CASE(3)
        WRITE(res,'(A,A,A,A,A,A,A)') '"',TRIM(ADJUSTL(buf(1))),&
          ' ',TRIM(ADJUSTL(buf(2))),' ',TRIM(ADJUSTL(buf(3))),'"'
      CASE(4)
        WRITE(res,'(A,A,A,A,A,A,A,A,A)') '"',TRIM(ADJUSTL(buf(1))),&
          ' ',TRIM(ADJUSTL(buf(2))),' ',TRIM(ADJUSTL(buf(3))), &
          ' ',TRIM(ADJUSTL(buf(4))),'"'
      END SELECT
    ELSE
      WRITE(buf(1),'(I8)') dims(1)
      WRITE(res,'(A,A,A)') '"',TRIM(ADJUSTL(buf(1))),'"'
    END IF
    DEALLOCATE(buf)
  END FUNCTION GetDimsStr

  !> Writes mesh to a file
  !!
  SUBROUTINE WriteMesh_xdmf(this,Mesh,offset)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this      !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)   :: Mesh      !< \param [in] mesh mesh type
#ifdef PARALLEL
    !> \param [in,out] offset
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset  
#else
    INTEGER          :: offset     !< \param [in,out] offset
#endif
    !------------------------------------------------------------------------!
    REAL,DIMENSION(:,:,:),ALLOCATABLE :: buf
    REAL,DIMENSION(:,:),POINTER :: buf2
    INTEGER          :: imax,jmax,i,err,dim
#ifdef PARALLEL
    INTEGER          :: request
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)       :: Mesh
    INTENT(INOUT)    :: this, offset
    !------------------------------------------------------------------------!
    IF(Mesh%INUM.EQ.Mesh%IMAX) THEN
      imax = Mesh%IMAX+1
    ELSE
      imax = Mesh%IMAX
    END IF
    IF(Mesh%JNUM.EQ.Mesh%JMAX) THEN
      jmax = Mesh%JMAX+1
    ELSE
      jmax = Mesh%JMAX
    END IF

    ALLOCATE(&
      buf2(Mesh%IMIN:imax,Mesh%JMIN:jmax),&
      STAT=this%error)
    IF (this%error.NE.0) &
       CALL Error(this,"WriteMesh_xdmf","Unable to allocate memory.")

    DO i=1,3
      SELECT CASE(i)
      CASE(1)
        buf2(:,:) = Mesh%ccart(Mesh%IMIN:imax,Mesh%JMIN:jmax,1,1)
      CASE(2)
        buf2(:,:) = Mesh%ccart(Mesh%IMIN:imax,Mesh%JMIN:jmax,1,2)
      CASE(3)
        IF(GetType(Mesh%geometry).EQ.BIANGLESPHERICAL) THEN
          buf2(:,:) = Mesh%ccart(Mesh%IMIN:imax,Mesh%JMIN:jmax,1,3)
        ELSE
          buf2(:,:) = 0.
        END IF
      END SELECT
#ifndef PARALLEL
      WRITE (this%unit) buf2
#else
      CALL MPI_File_set_view(this%handle,offset,DEFAULT_MPI_REAL,&
        this%meshtype, 'native', MPI_INFO_NULL, this%error)
      CALL MPI_File_write_all(this%handle,buf2,this%meshbufsize,DEFAULT_MPI_REAL,&
        this%status, this%error)
      offset = offset + (Mesh%INUM+1)*(Mesh%JNUM+1)*this%realsize
#endif
    END DO

    DEALLOCATE(buf2)
  END SUBROUTINE WriteMesh_xdmf 

  !> \public Specific routine to open a file for xdmf I/O
  !!
  SUBROUTINE OpenFile_xdmf(this,action)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this    !< \param [in,out] this fileio type
    INTEGER          :: action  !< \param [in] action mode of file access
    !------------------------------------------------------------------------!
    INTENT(IN)       :: action
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    CALL OpenFile(this,action,BIN)
#else
    CALL OpenFile_vtk(this,action)
#endif
  END SUBROUTINE OpenFile_xdmf


  !> \public do nothing
  !!
  SUBROUTINE WriteHeader_xdmf(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this  !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
  END SUBROUTINE WriteHeader_xdmf

  !> \public Reads the header of a file (not yet implemented)
  !!
  SUBROUTINE ReadHeader_xdmf(this,success)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this     !< \param [in,out] this fileio type
    LOGICAL          :: success  !< \param [out] success
    !------------------------------------------------------------------------!
    INTENT(OUT)      :: success
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    success = .FALSE.
  END SUBROUTINE ReadHeader_xdmf

  !> \public Writes the timestep (not yet implemented)
  !!
  SUBROUTINE WriteTimestamp_xdmf(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this     !< \param [in,out] this fileio type
    REAL             :: time     !< \param [in] time
    !------------------------------------------------------------------------!
    INTENT(IN)       :: time
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
  END SUBROUTINE WriteTimestamp_xdmf

  !> \public Writes the timestep (not yet implemented)
  !!
  SUBROUTINE ReadTimestamp_xdmf(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this     !< \param [in,out] this fileio type
    REAL             :: time     !< \param [out] time
    !------------------------------------------------------------------------!
    INTENT(OUT)      :: time
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
  END SUBROUTINE ReadTimestamp_xdmf

  !> Writes data attributes to a file
  !!
  RECURSIVE SUBROUTINE WriteDataAttributes_xdmf(this,Mesh,config,offset,path)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this         !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh         !< \param [in] mesh mesh type
    TYPE(Dict_TYP),POINTER :: config  !< \param [in] config dict of configuration
    CHARACTER(LEN=*),OPTIONAL &
                           :: path    !< \param [in] path 
#ifdef PARALLEL
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
#else
    INTEGER                :: offset  !< \param [in,out] offset 
#endif
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: dir
    CHARACTER(LEN=MAX_CHAR_LEN) :: str, key
    TYPE(dict_TYP),POINTER   :: node
    REAL,POINTER           :: realptr
    REAL,DIMENSION(:,:),POINTER :: ptr2
    REAL,DIMENSION(:,:,:),POINTER :: ptr3
    REAL,DIMENSION(:,:,:,:),POINTER :: ptr4
    INTEGER                :: dims2(2),dims3(3),dims4(4)
    REAL,DIMENSION(2,Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) :: buf
    INTEGER                :: i
#ifdef PARALLEL
    INTEGER                :: request,j,k
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh, path
    INTENT(INOUT)     :: this, offset
    !------------------------------------------------------------------------!
    IF(PRESENT(path)) THEN
        str = path
    ELSE
        str = ""
    ENDIF
    node => config
    DO WHILE(ASSOCIATED(node))
      key = TRIM(str)//"/"//TRIM(GetKey(node))
      SELECT CASE(GetDataType(node))
      CASE(DICT_DIR)
        IF(GetKey(node).NE."mesh") THEN
          key = TRIM(key)
          CALL GetAttr(node,GetKey(node),dir)
          CALL WriteDataAttributes_xdmf(this, Mesh, dir, offset, key)
        END IF
      CASE(DICT_REAL)
        CALL GetAttr(node,"value",realptr)
#ifndef PARALLEL
        WRITE(this%unit) realptr
#else
        CALL MPI_File_write_all(this%handle,realptr,1,DEFAULT_MPI_REAL, &
             this%status,this%error)
#endif        
        offset = offset + this%realsize
      CASE(DICT_REAL_TWOD)
!        print *,"writing: ", TRIM(key), " offset: ",offset*8
        CALL GetAttr(node,"value",ptr2)
        dims2 = SHAPE(ptr2)
        IF((dims2(1).EQ.(Mesh%IGMAX-Mesh%IGMIN+1)).AND.(dims2(2).EQ.(Mesh%JGMAX-Mesh%JGMIN+1))) THEN
#ifndef PARALLEL
          WRITE(this%unit) ptr2(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)
#else
          CALL MPI_File_set_view(this%handle,offset,DEFAULT_MPI_REAL,&
               this%filetype, 'native', MPI_INFO_NULL, this%error)
          CALL MPI_File_write_all(this%handle,ptr2,1,this%memtype,&
                                  this%status,this%error)
          offset = offset + GetSize(Mesh,dims2)*this%realsize
#endif
        ELSE
#ifndef PARALLEL
          WRITE(this%unit) ptr2(:,:)
#else
          !CALL MPI_File_write_all(this%handle,ptr2(:,:),PRODUCT(dims2),&
          !                        DEFAULT_MPI_REAL,this%status,this%error)
          IF(GetRank(this).EQ.0) THEN
            CALL MPI_File_write_at(this%handle, offset, ptr2, PRODUCT(dims2), &
                                   DEFAULT_MPI_REAL,this%status, this%error)
          END IF
#endif
          offset = offset + PRODUCT(dims2)*this%realsize
        END IF
      CASE(DICT_REAL_THREED)
        CALL GetAttr(node,"value",ptr3)
#ifndef PARALLEL
        WRITE(this%unit) ptr3(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:)
#else
        dims3 = SHAPE(ptr3)
        DO i=1,dims3(3)
          CALL MPI_File_set_view(this%handle,offset,DEFAULT_MPI_REAL,&
               this%filetype, 'native', MPI_INFO_NULL, this%error)
          CALL MPI_File_write_all(this%handle,&
                                  ptr3(:,:,i), &
                                  1,this%memtype,this%status,this%error)
          offset = offset + GetSize(Mesh,dims3(1:2))*this%realsize
        END DO
#endif
      CASE(DICT_REAL_FOURD)
        CALL GetAttr(node,"value",ptr4)
#ifndef PARALLEL
        !WRITE(this%unit) ptr4(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:,:)
#else
        dims4 = SHAPE(ptr4)
        !offset = offset + GetSize(Mesh,dims4)*this%realsize
#endif
      END SELECT
      node=>GetNext(node)
    END DO
  END SUBROUTINE WriteDataAttributes_xdmf

  !> \public Writes all desired data arrays to a file 
  !!
  SUBROUTINE WriteDataset_xdmf(this,Mesh,Physics,Fluxes,Timedisc,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this      !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh      !< \param [in] mesh mesh type
    TYPE(Physics_TYP) :: Physics   !< \param [in] physics physics type
    TYPE(Fluxes_TYP)  :: Fluxes    !< \param [in] fluxes fluxes type
    TYPE(Timedisc_TYP):: Timedisc  !< \param [in] timedisc timedisc type
    TYPE(Dict_TYP),POINTER :: IO   !< \param [in,out] IO I/O dictionary
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: meshIO
#ifdef PARALLEL
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
#else
    INTEGER           :: offset
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Fluxes,Timedisc,IO
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    CALL MPI_File_get_size(this%handle,offset,this%error)
    CALL MPI_Barrier(MPI_COMM_WORLD,this%error)
#else
    offset = 0
#endif

#ifndef PARALLEL
    WRITE(this%unit) Timedisc%time
#else
    IF (GetRank(this).EQ.0) THEN
      CALL MPI_File_write(this%handle,Timedisc%time,1,DEFAULT_MPI_REAL,&
        this%status,this%error)
    END IF
#endif
    offset = offset + this%realsize

    ! write data
    CALL WriteDataAttributes_xdmf(this,Mesh,IO,offset)

    ! only write mesh arrays in first data file
    IF((this%step.EQ.0).OR.(this%cycles.EQ.0)) THEN
      CALL WriteMesh_xdmf(this,Mesh,offset)
      CALL GetAttr(IO,"mesh",meshIO)
      CALL WriteDataAttributes_xdmf(this,Mesh,meshIO,offset,"/mesh")
    END IF
    
  END SUBROUTINE WriteDataset_xdmf

  !> \public Reads the data arrays from file (not yet implemented)
  !!
  SUBROUTINE ReadDataset_xdmf(this,Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP)  :: this      !< \param [in,out] this fileio type
    TYPE(Mesh_TYP)    :: Mesh      !< \param [in] mesh mesh type
    TYPE(Physics_TYP) :: Physics   !< \param [in] physics physics type
    TYPE(Timedisc_TYP):: Timedisc  !< \param [in,out] timedisc timedisc type
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: this,Timedisc
    !------------------------------------------------------------------------!
  END SUBROUTINE ReadDataset_xdmf

  !> \public Closes the file I/O
  !!
  SUBROUTINE CloseFileIO_xdmf(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(FileIO_TYP) :: this
    !------------------------------------------------------------------------!
  END SUBROUTINE CloseFileIO_xdmf

END MODULE fileio_xdmf
