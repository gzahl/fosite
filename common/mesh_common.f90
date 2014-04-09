!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: mesh_common.f90                                                   #
!#                                                                           #
!# Copyright (C) 2006-2008                                                   #
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
! basic mesh module
!----------------------------------------------------------------------------!
MODULE mesh_common
  USE common_types, GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Info_common => Info, Warning_common => Warning, Error_common => Error
  USE geometry_common, ONLY : Geometry_TYP
  IMPLICIT NONE
#ifdef PARALLEL
  include 'mpif.h'
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: NDIMS = 2
  !--------------------------------------------------------------------------!
  INTERFACE GetRank
     MODULE PROCEDURE GetMeshRank, GetRank_common
  END INTERFACE
  INTERFACE GetNumProcs
     MODULE PROCEDURE GetMeshNumProcs, GetNumProcs_common
  END INTERFACE
  INTERFACE Info
     MODULE PROCEDURE MeshInfo, Info_common
  END INTERFACE
  INTERFACE Warning
     MODULE PROCEDURE MeshWarning, Warning_common
  END INTERFACE
  INTERFACE Error
     MODULE PROCEDURE MeshError_rank0, MeshError_rankX, Error_common
  END INTERFACE
  !--------------------------------------------------------------------------!
  ! mesh data structure
  TYPE Mesh_TYP
     TYPE(Common_TYP)           :: mtype           ! mesh type               !
     TYPE(Geometry_TYP)         :: Geometry        ! geometrical properties  !
     INTEGER                    :: GNUM            ! num. ghost cells        !
     INTEGER                    :: INUM,JNUM       ! resolution              !
     INTEGER                    :: IMIN,IMAX       ! min. & max. index in x- !
     INTEGER                    :: JMIN,JMAX       ! and y-direction         !
     INTEGER                    :: IGMIN,IGMAX     ! same with ghost cells   !
     INTEGER                    :: JGMIN,JGMAX     !                         !
#ifdef PARALLEL
     INTEGER                    :: MAXINUM,MAXJNUM ! max. of local INUM,JNUM !
     INTEGER                    :: comm_cart       ! cartesian communicator  !
     INTEGER, DIMENSION(4)      :: comm_boundaries ! comm. for bound. procs. !
     INTEGER, DIMENSION(4)      :: rank0_boundaries! map rank0 -> world rank !
     INTEGER, DIMENSION(4)      :: neighbor        ! ranks of neighbor proc. !
     INTEGER, DIMENSION(NDIMS)  :: dims            ! dimensions of cart comm !
     INTEGER, DIMENSION(NDIMS)  :: mycoords        ! par. proc coordinates   !
!!$     INTEGER                 :: weblocks, snblocks ! boundary data handles   !
#endif
     REAL                       :: xmin, xmax      ! comput. domain in x-    !
     REAL                       :: ymin, ymax      ! and y-direction         !
     REAL                       :: dx,dy           ! width of cells in x-&   !
     REAL                       :: invdx, invdy    ! y-direction; inverse    !
     REAL, POINTER              :: center(:,:,:)   ! cell geometr. centers   !
     REAL, POINTER              :: bcenter(:,:,:)  ! cell bary centers       !
     REAL, POINTER              :: bccart(:,:,:)   ! cartesian bary centers  !
     REAL, POINTER              :: fpos(:,:,:,:)   ! face centered positions !
     REAL, POINTER              :: cpos(:,:,:,:)   ! corner positions        !
     REAL, POINTER              :: fhx(:,:,:)      ! face centered scale     !
     REAL, POINTER              :: fhy(:,:,:)      !  factors                !
     REAL, POINTER              :: fhz(:,:,:)      !                         !
     REAL, POINTER              :: chx(:,:,:)      ! corner scale factors    !
     REAL, POINTER              :: chy(:,:,:)      !                         !
     REAL, POINTER              :: chz(:,:,:)      !                         !
     REAL, POINTER              :: bhx(:,:)        ! bary center scale       !
     REAL, POINTER              :: bhy(:,:)        !  factors                !
     REAL, POINTER              :: bhz(:,:)        !                         !
     REAL, POINTER              :: cyxy(:,:,:)     ! commutator coefficients !
     REAL, POINTER              :: cxyx(:,:,:)     ! for geometrical sources !
     REAL, POINTER              :: czxz(:,:,:)
     REAL, POINTER              :: czyz(:,:,:)
     REAL, POINTER              :: sqrtg(:,:)      ! sqrt(det(g_ij))         !
     REAL, POINTER              :: weights(:,:,:,:)! interpolation weights   !
     REAL, POINTER              :: dlx(:,:)        ! line elements           !
     REAL, POINTER              :: dly(:,:)        !   at cell centers       !
     REAL, POINTER              :: dAx(:,:,:)      ! surface elements        !
     REAL, POINTER              :: dAy(:,:,:)      !   on cell faces         !
     REAL, POINTER              :: dAxdy(:,:,:)    ! surface elements divided!
     REAL, POINTER              :: dAydx(:,:,:)    !   by dx or dy           !
     REAL, POINTER              :: volume(:,:)     ! cell volumes            !
     REAL, POINTER              :: dxdV(:,:)       ! inverse volume elements !
     REAL, POINTER              :: dydV(:,:)       ! multiplied with dx or dy!
  END TYPE Mesh_TYP
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Mesh_TYP, &
#ifdef PARALLEL
       DEFAULT_MPI_REAL, &
#endif
       ! methods
       InitMesh, &
       GetRank, &
       GetNumProcs, &
       Info, &
       Warning, &
       Error, &
       CloseMesh
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitMesh(this,mtype,mname,inum,jnum,xmin,xmax,ymin,ymax)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    INTEGER           :: mtype
    CHARACTER(LEN=*)  :: mname
    INTEGER           :: inum,jnum
    REAL              :: xmin,xmax,ymin,ymax    
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
#ifdef PARALLEL
    INTEGER           :: maxinum, maxjnum
#endif    
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: inum,jnum,xmin,xmax,ymin,ymax
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!

    CALL InitCommon(this%mtype,mtype,mname)

    ! total resolution
    this%INUM = inum
    this%JNUM = jnum

    ! number of ghost rows/columns
    this%GNUM = 2

    ! coordinate domain
    this%xmin = xmin
    this%xmax = xmax
    this%ymin = ymin
    this%ymax = ymax

    ! coordinate differences in each direction
    this%dx = (xmax - xmin) / inum
    this%dy = (ymax - ymin) / jnum
    
    ! inverse coordinate differences
    this%invdx = 1./this%dx
    this%invdy = 1./this%dy

    ! set index ranges
#ifdef PARALLEL
    CALL InitMesh_parallel(this)
    CALL MPI_Barrier(MPI_COMM_WORLD,err)
    maxinum = this%IMAX - this%IMIN + 1
    CALL MPI_Allreduce(maxinum,this%MAXINUM,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,err)
    maxjnum = this%JMAX - this%JMIN + 1
    CALL MPI_Allreduce(maxjnum,this%MAXJNUM,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,err)
#else
    this%IMIN   = 1
    this%IMAX   = this%INUM
    this%JMIN   = 1
    this%JMAX   = this%JNUM
#endif

    ! index ranges with ghost cells
    this%IGMIN = this%IMIN - this%GNUM
    this%IGMAX = this%IMAX + this%GNUM
    this%JGMIN = this%JMIN - this%GNUM
    this%JGMAX = this%JMAX + this%GNUM

    ! allocate memory for all pointers that are independent of fluxtype
    ALLOCATE(this%center(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2), &
         this%bcenter(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2), &
         this%bccart(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2), &
         this%fpos(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4,2), &
         this%bhx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         this%bhy(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         this%bhz(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         this%fhx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4), &
         this%fhy(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4), &
         this%fhz(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4), &
         this%volume(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         this%dxdV(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         this%dydV(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         this%dlx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         this%dly(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         STAT=err)
    IF (err.NE.0) THEN
       CALL Error(this,"InitMesh_common","Unable to allocate memory!")
    END IF

    ! calculate geometrical cell centers
    FORALL (i=this%IGMIN:this%IGMAX,j=this%JGMIN:this%JGMAX)
       this%center(i,j,1) = this%xmin + (i - 0.5)*this%dx
       this%center(i,j,2) = this%ymin + (j - 0.5)*this%dy
    END FORALL

    ! calculate face centered positions
    FORALL (i=this%IGMIN:this%IGMAX,j=this%JGMIN:this%JGMAX)
       this%fpos(i,j,1,1) = this%center(i,j,1) - 0.5*this%dx  ! western x coord
       this%fpos(i,j,1,2) = this%center(i,j,2)                ! western y coord
       this%fpos(i,j,2,1) = this%center(i,j,1) + 0.5*this%dx  ! eastern x coord
       this%fpos(i,j,2,2) = this%center(i,j,2)                ! eastern y coord
       this%fpos(i,j,3,1) = this%center(i,j,1)                ! southern x coord
       this%fpos(i,j,3,2) = this%center(i,j,2) - 0.5*this%dy  ! southern y coord
       this%fpos(i,j,4,1) = this%center(i,j,1)                ! northern x coord
       this%fpos(i,j,4,2) = this%center(i,j,2) + 0.5*this%dy  ! northern y coord
    END FORALL
  END SUBROUTINE InitMesh


#ifdef PARALLEL
  SUBROUTINE InitMesh_parallel(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP) :: this
    !------------------------------------------------------------------------!
    LOGICAL, DIMENSION(NDIMS), PARAMETER :: periods = .FALSE.
    INTEGER        :: num,rem
    INTEGER        :: ierror
    INTEGER        :: i,j
    INTEGER        :: worldgroup,newgroup
    INTEGER, DIMENSION(NDIMS) :: coords
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ranks
    !------------------------------------------------------------------------!
    INTENT(INOUT)  :: this
    !------------------------------------------------------------------------!

    ! create a cartesian topology of processes
    ! 1. balance number of processes per direction
    this%dims(1)=GetNumProcs(this)
    this%dims(2)=1
    CALL CalculateDecomposition(this%INUM,this%JNUM,this%dims(1),this%dims(2))
    IF (this%dims(2).LE.0) THEN
       CALL Error(this,"InitMesh_parallel","Domain decomposition algorithm failed.")
    END IF

    ! 2. create the cartesian communicator
    ! IMPORTANT: disable reordering of nodes
    CALL MPI_Cart_create(MPI_COMM_WORLD,NDIMS,this%dims,periods,.TRUE., &
         this%comm_cart,ierror)

    ! 3. inquire and save the own position
    CALL MPI_Cart_get(this%comm_cart,NDIMS,this%dims,periods,this%mycoords,ierror)

    ! subdivide the mesh and set mesh indices
    ! x-direction
    rem = MOD(this%INUM,this%dims(1))            ! remainder
    num = (this%INUM-rem) / this%dims(1)         ! fraction
    IF (this%mycoords(1).LT.rem) THEN
       ! the first (rem-1) nodes get one more to account for the remainder
       this%IMIN = this%mycoords(1) * (num + 1) + 1
       this%IMAX = this%IMIN + num
    ELSE
       this%IMIN = rem * (num+1) + (this%mycoords(1) - rem) * num + 1
       this%IMAX = this%IMIN + num - 1
    END IF
    ! y-direction
    rem = MOD(this%JNUM,this%dims(2))            ! remainder
    num = (this%JNUM-rem) / this%dims(2)         ! fraction
    IF (this%mycoords(2).LT.rem) THEN
       ! the first (rem-1) nodes get one more to account for the remainder
       this%JMIN = this%mycoords(2) * (num + 1) + 1
       this%JMAX = this%JMIN + num
    ELSE
       this%JMIN = rem * (num+1) + (this%mycoords(2) - rem) * num + 1
       this%JMAX = this%JMIN + num - 1
    END IF

    ! create communicators for all boundaries
    ALLOCATE(ranks(MAX(this%dims(1),this%dims(2))))
    ranks = 0
    CALL MPI_Comm_group(MPI_COMM_WORLD,worldgroup,ierror)
    ! western boundary
    coords(1) = 0
    DO j=0,this%dims(2)-1
       coords(2) = j
       CALL MPI_Cart_rank(this%comm_cart,coords,i,ierror)
       ranks(j+1)=i
    END DO
    CALL MPI_Group_incl(worldgroup,this%dims(2),ranks,newgroup,ierror)
    CALL MPI_Comm_create(this%comm_cart,newgroup,this%comm_boundaries(1),ierror)
    CALL MPI_Group_translate_ranks(newgroup,1,0,worldgroup,&
         this%rank0_boundaries(1),ierror)
    CALL MPI_Group_free(newgroup,ierror)
    ! eastern boundary
    coords(1) = this%dims(1)-1
    DO j=0,this%dims(2)-1
       coords(2) = j
       CALL MPI_Cart_rank(this%comm_cart,coords,i,ierror)
       ranks(j+1)=i
    END DO
    CALL MPI_Group_incl(worldgroup,this%dims(2),ranks,newgroup,ierror)
    CALL MPI_Comm_create(this%comm_cart,newgroup,this%comm_boundaries(2),ierror)
    CALL MPI_Group_translate_ranks(newgroup,1,0,worldgroup,&
         this%rank0_boundaries(2),ierror)
    CALL MPI_Group_free(newgroup,ierror)
    ! southern boundary
    coords(2) = 0
    DO i=0,this%dims(1)-1
       coords(1) = i
       CALL MPI_Cart_rank(this%comm_cart,coords,j,ierror)
       ranks(i+1)=j
    END DO
    CALL MPI_Group_incl(worldgroup,this%dims(1),ranks,newgroup,ierror)
    CALL MPI_Comm_create(this%comm_cart,newgroup,this%comm_boundaries(3),ierror)
    CALL MPI_Group_translate_ranks(newgroup,1,0,worldgroup, &
         this%rank0_boundaries(3),ierror)
    CALL MPI_Group_free(newgroup,ierror)
    ! northern boundary
    coords(2) = this%dims(2)-1
    DO i=0,this%dims(1)-1
       coords(1) = i
       CALL MPI_Cart_rank(this%comm_cart,coords,j,ierror)
       ranks(i+1)=j
    END DO
    CALL MPI_Group_incl(worldgroup,this%dims(1),ranks,newgroup,ierror)
    CALL MPI_Comm_create(this%comm_cart,newgroup,this%comm_boundaries(4),ierror)
    CALL MPI_Group_translate_ranks(newgroup,1,0,worldgroup, &
         this%rank0_boundaries(4),ierror)
    CALL MPI_Group_free(newgroup,ierror)
    CALL MPI_Group_free(worldgroup,ierror)
    DEALLOCATE(ranks)
  END SUBROUTINE InitMesh_parallel


  ! return the best partitioning of processes
  ! pj x pk for a given mesh with resolution nj x nk
  ! pj : total number of processes (<10000)
  ! pk = 1 initially
  ! pk return 0, for erroneous input
  SUBROUTINE CalculateDecomposition(nj,nk,pj,pk)
    USE factors
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER, INTENT(IN)    :: nj,nk
    INTEGER, INTENT(INOUT) :: pj,pk
    !------------------------------------------------------------------------!
    INTEGER :: pf,p1,p2,pp,pn,tmp,bl,bl1
    !------------------------------------------------------------------------!
    ! return immediatly for malformed input
    IF (((pj.LT.2).AND.(pk.NE.1)).OR.(pj.GT.10000) &
        .OR.(nj*nk.LT.pj)) THEN
        pk = 0
        RETURN
    END IF
    p1=pj
    p2=pk
    pp=pj
    pn=pj
    pf=1
    tmp=0
    bl1=MAX(nj*pk+nk*pj,nj*pj+nk*pk)
    DO
       IF (pf.GT.tmp) THEN
          tmp=pf
          p1=pn/pf
          p2=pf
          bl = Balance(nj,nk,p1,p2)
          IF (bl.LT.bl1) THEN
             bl1 = bl
             pj = p1
             pk = p2
          END IF
       END IF
       IF (pp.EQ.1) EXIT
       pf=Reduce(pp)
    END DO

  CONTAINS
    
    RECURSIVE FUNCTION Balance(nj,nk,pj,pk) RESULT(bl)
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      INTEGER, INTENT(IN)    :: nj,nk
      INTEGER, INTENT(INOUT) :: pj,pk
      INTEGER :: bl
      !------------------------------------------------------------------------!
      INTEGER :: n1,n2,p1,p2,boundlen
      INTEGER :: primfac,bl1,bl2
      !------------------------------------------------------------------------!
      boundlen(n1,n2,p1,p2) = n1*(p2-1) + n2*(p1-1)
      !------------------------------------------------------------------------!
      bl = boundlen(nj,nk,pj,pk)
      IF (pj.EQ.1) RETURN
      primfac = GetFactor(pj)
      p1 = pj/primfac
      p2 = pk*primfac
      bl1 = Balance(nj,nk,p1,p2)
      IF (bl.GT.bl1) THEN
         pj = p1
         pk = p2
         bl = bl1
      END IF
    END FUNCTION Balance
    
  END SUBROUTINE CalculateDecomposition
#endif


  PURE FUNCTION GetMeshRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP), INTENT(IN) :: this
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%mtype)
  END FUNCTION GetMeshRank


  PURE FUNCTION GetMeshNumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP), INTENT(IN) :: this
    INTEGER :: p
    !------------------------------------------------------------------------!
    p = GetNumProcs_common(this%mtype)
  END FUNCTION GetMeshNumProcs


  SUBROUTINE MeshInfo(this,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: msg
    !------------------------------------------------------------------------!
    CALL Info_common(this%mtype,msg)
  END SUBROUTINE MeshInfo


  SUBROUTINE MeshWarning(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Warning_common(this%mtype,modproc,msg)
  END SUBROUTINE MeshWarning


  SUBROUTINE MeshError_rank0(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP), INTENT(IN)    :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Error_common(this%mtype,modproc,msg)
  END SUBROUTINE MeshError_rank0


  SUBROUTINE MeshError_rankX(this,modproc,msg,rank)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP), INTENT(IN)    :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    INTEGER, INTENT(IN)           :: rank
    !------------------------------------------------------------------------!
    CALL Error_common(this%mtype,modproc,msg,rank)
  END SUBROUTINE MeshError_rankX


  SUBROUTINE CloseMesh(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%center,this%bcenter,this%bccart,this%fpos, &
         this%bhx,this%bhy,this%bhz, &
         this%fhx,this%fhy,this%fhz, &
         this%volume,this%dxdV,this%dydV,this%dlx,this%dly)
  END SUBROUTINE CloseMesh

END MODULE mesh_common
