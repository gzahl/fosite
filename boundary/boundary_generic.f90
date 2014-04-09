!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_generic.f90                                              #
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
!> \addtogroup boundary
!! \key{western,INTEGER,boundary condition at western boundary}
!! \key{eastern,INTEGER,boundary condition at eastern boundary}
!! \key{southern,INTEGER,boundary condition at southern boundary}
!! \key{northern,INTEGER,boundary condition at northern boundary}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!!
!! \brief Generic boundary module
!!
!! This module provides the generic interface routines to all boundary
!! modules.
!!
!----------------------------------------------------------------------------!
MODULE boundary_generic
  USE timedisc_common, ONLY : Timedisc_TYP
  USE mesh_common, ONLY : Mesh_TYP, GetRank, Initialized, Error
  USE fluxes_common, ONLY : Fluxes_TYP
  USE reconstruction_common, ONLY : Reconstruction_TYP, PrimRecon
  USE boundary_nogradients, InitBoundary_common => InitBoundary, &
       CloseBoundary_common => CloseBoundary
  USE boundary_periodic
  USE boundary_reflecting, InitBoundary_common1 => InitBoundary,&
       CloseBoundary_common1 => CloseBoundary
  USE boundary_axis
  USE boundary_folded
  USE boundary_fixed
  USE boundary_extrapolation
  USE boundary_noh
  USE boundary_noslip
  USE boundary_absorbing
  USE boundary_custom
  USE boundary_farfield
  USE physics_generic, ONLY : Physics_TYP, Initialized, &
       Convert2Primitive, Convert2Conservative, &
       EULER2D_ISOIAMT, EULER2D_IAMT, GetType
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
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE CloseBoundary
     MODULE PROCEDURE CloseBoundary_one, CloseBoundary_all
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
#ifdef PARALLEL
  INTEGER, PARAMETER :: NONE            = 0
#endif
  INTEGER, PARAMETER :: NO_GRADIENTS    = 1
  INTEGER, PARAMETER :: PERIODIC        = 2
  INTEGER, PARAMETER :: REFLECTING      = 3
  INTEGER, PARAMETER :: AXIS            = 4
  INTEGER, PARAMETER :: FOLDED          = 5
  INTEGER, PARAMETER :: FIXED           = 6
  INTEGER, PARAMETER :: EXTRAPOLATION   = 7
  INTEGER, PARAMETER :: NOH2D           = 8
  INTEGER, PARAMETER :: NOH3D           = 9
  INTEGER, PARAMETER :: NOSLIP          = 10
  INTEGER, PARAMETER :: CUSTOM          = 11
  INTEGER, PARAMETER :: FARFIELD        = 12
  INTEGER, PARAMETER :: ABSORBING       = 13
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Boundary_TYP, &
       ! constants
       WEST, EAST, SOUTH, NORTH, &
       NO_GRADIENTS, PERIODIC, REFLECTING, AXIS, FOLDED, FIXED, EXTRAPOLATION, &
       NOH2D, NOH3D, NOSLIP, CUSTOM, FARFIELD, ABSORBING, &
       CUSTOM_NOGRAD, CUSTOM_PERIOD, CUSTOM_REFLECT, CUSTOM_REFLNEG, &
       CUSTOM_EXTRAPOL, CUSTOM_FIXED, CUSTOM_LOGEXPOL, &
       CUSTOM_OUTFLOW, CUSTOM_KEPLER, CUSTOM_ANGKEPLER, CUSTOM_POISSON, &
#ifdef PARALLEL
       NONE,&
#endif
       ! methods
       InitBoundary, &
       CenterBoundary, &
       CloseBoundary, &
       GetType, &
       GetName, &
       GetDirection, &
       GetDirectionName, &
       GetRank, &
       GetNumProcs, &
       Initialized, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitBoundary_one(this,Mesh,Physics,btype,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: btype,dir
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER, PARAMETER :: strlen = 32
    CHARACTER(LEN=strlen) :: sendbuf
    CHARACTER(LEN=strlen) :: recvbuf
    INTEGER            :: status(MPI_STATUS_SIZE)
    INTEGER            :: ierror
    INTEGER            :: i
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,btype,dir
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF(((GetType(Physics).EQ.EULER2D_ISOIAMT).OR.&
        (GetType(Physics).EQ.EULER2D_IAMT)).AND.&
       ((dir.EQ.NORTH).OR.(dir.EQ.SOUTH)).AND.&
       (.NOT.((btype.EQ.PERIODIC) &
#ifdef PARALLEL
              .OR.(btype.EQ.NONE) &
#endif
              ))) &
      CALL Error(this, "InitBoundary_one", "All IAMT Physics need periodic" & 
        // " boundary conditions in NORTH/SOUTH direction")

    ! set boundary properties
    SELECT CASE(btype)
#ifdef PARALLEL
    CASE(NONE)
       CALL InitBoundary_common(this,btype,"none",dir)
#endif
    CASE(NO_GRADIENTS)
       CALL InitBoundary_nogradients(this,btype,dir)
    CASE(PERIODIC)
       CALL InitBoundary_periodic(this,btype,dir)
    CASE(REFLECTING)
       CALL InitBoundary_reflecting(this,Physics,btype,dir)
    CASE(AXIS)
       CALL InitBoundary_axis(this,Physics,btype,dir)
    CASE(FOLDED)
       CALL InitBoundary_folded(this,Mesh,Physics,btype,dir)
    CASE(FIXED)
       CALL InitBoundary_fixed(this,Mesh,Physics,btype,dir)
    CASE(EXTRAPOLATION)
       CALL InitBoundary_extrapolation(this,Physics,btype,dir)
    CASE(NOH2D)
       CALL InitBoundary_noh(this,Mesh,Physics,btype,dir,2)
    CASE(NOH3D)
       CALL InitBoundary_noh(this,Mesh,Physics,btype,dir,3)
    CASE(NOSLIP)
       CALL InitBoundary_noslip(this,Mesh,Physics,btype,dir)
    CASE(CUSTOM)
       CALL InitBoundary_custom(this,Mesh,Physics,btype,dir)
    CASE(FARFIELD)
       CALL InitBoundary_farfield(this,Mesh,Physics,btype,dir)
    CASE(ABSORBING)
       CALL InitBoundary_absorbing(this,Mesh,Physics,btype,dir)
    CASE DEFAULT
       CALL Error(this, "InitBoundary_one", "Unknown boundary condition.")
    END SELECT

    ! print some information
#ifdef PARALLEL
    ! send boundary information to the rank 0 process;
    ! we only need this to synchronize the output
    IF (GetRank(this) .EQ. 0 .AND. GetRank(this).EQ.Mesh%rank0_boundaries(dir)) THEN
       ! print output without communication 
#endif
       CALL Info(this, " BOUNDARY-> condition:         " //  TRIM(GetDirectionName(this)) &
            // " " // TRIM(GetName(this)), GetRank(this))
#ifdef PARALLEL
    ELSE IF (GetRank(this).EQ.Mesh%rank0_boundaries(dir)) THEN
       ! send info to root
       sendbuf = TRIM(GetDirectionName(this))//" "//TRIM(GetName(this))
       CALL MPI_SEND(sendbuf,strlen,MPI_CHARACTER,0, &
                     0,MPI_COMM_WORLD,ierror)
    ELSE IF (GetRank(this).EQ.0) THEN
       ! receive input from rank0_boundaries(dir)
       CALL MPI_RECV(recvbuf,strlen,MPI_CHARACTER,Mesh%rank0_boundaries(dir),&
                     MPI_ANY_TAG,MPI_COMM_WORLD,status,ierror)
       CALL Info(this, " BOUNDARY-> condition:         " // TRIM(recvbuf),GetRank(this))
    END IF
#endif
  END SUBROUTINE InitBoundary_one


  SUBROUTINE InitBoundary(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), DIMENSION(4) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Dict_TYP),POINTER &
                       :: config,IO
    INTEGER            :: western, eastern, southern, northern
    !------------------------------------------------------------------------!
    INTEGER            :: new_western, new_eastern, new_southern, new_northern
#ifdef PARALLEL
    INTEGER            :: comm_old, dir
    INTEGER            :: sizeofreal, ignum, jgnum, twoslices
    INTEGER            :: ierr
    LOGICAL, DIMENSION(SIZE(Mesh%dims)) :: periods = .FALSE.
    LOGICAL, DIMENSION(SIZE(Mesh%dims)) :: remain_dims = .FALSE.
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Physics
    INTENT(INOUT) :: this,Mesh
    !------------------------------------------------------------------------!
    IF (.NOT.Initialized(Physics).OR..NOT.Initialized(Mesh)) &
         CALL Error(this(WEST),"InitBoundary","physics and/or mesh module uninitialized")    

    CALL RequireKey(config, "western")
    CALL RequireKey(config, "eastern")
    CALL RequireKey(config, "southern")
    CALL RequireKey(config, "northern")

    CALL GetAttr(config, "western", western)
    CALL GetAttr(config, "eastern", eastern)
    CALL GetAttr(config, "southern", southern)
    CALL GetAttr(config, "northern", northern)

    new_western = western
    new_eastern = eastern
    new_southern = southern
    new_northern = northern

#ifdef PARALLEL
    ! define connections
    IF (Mesh%mycoords(1).NE.0)  new_western = NONE
    IF (Mesh%mycoords(1).NE.Mesh%dims(1)-1)  new_eastern = NONE
    IF (Mesh%mycoords(2).NE.0)  new_southern = NONE
    IF (Mesh%mycoords(2).NE.Mesh%dims(2)-1)  new_northern = NONE
#endif
    
    ! initialize every boundary
    ! IMPORTANT: do this before anything else
    CALL InitBoundary_one(this(WEST),Mesh,Physics,new_western,WEST)
    CALL InitBoundary_one(this(EAST),Mesh,Physics,new_eastern,EAST)
    CALL InitBoundary_one(this(SOUTH),Mesh,Physics,new_southern,SOUTH)
    CALL InitBoundary_one(this(NORTH),Mesh,Physics,new_northern,NORTH)
    
#ifdef PARALLEL
    ! check periodicity
    IF (western.EQ.PERIODIC.AND.eastern.EQ.PERIODIC) THEN
       periods(1) = .TRUE.
    ELSE IF (western.EQ.PERIODIC.NEQV.eastern.EQ.PERIODIC) THEN
       CALL Error(this(WEST),"InitBoundary", &
            "Opposite boundary should be periodic.")
    END IF
    IF (southern.EQ.PERIODIC.AND.northern.EQ.PERIODIC) THEN
       periods(2) = .TRUE.
    ELSE IF (southern.EQ.PERIODIC.NEQV.northern.EQ.PERIODIC) THEN
       CALL Error(this(SOUTH),"InitBoundary", &
            "Opposite boundary should be periodic.")
    END IF

    ! create new cartesian communicator using Mesh%comm_cart
    ! and account for the periodicity
    ! IMPORTANT: disable reordering of nodes
    comm_old = Mesh%comm_cart
    CALL MPI_Cart_create(comm_old,SIZE(Mesh%dims),Mesh%dims,periods,.FALSE., &
         Mesh%comm_cart,ierr)

    ! save ranks of neighbor processes
    CALL MPI_Cart_shift(Mesh%comm_cart,0,1,Mesh%neighbor(WEST),Mesh%neighbor(EAST),ierr)
    CALL MPI_Cart_shift(Mesh%comm_cart,1,1,Mesh%neighbor(SOUTH),Mesh%neighbor(NORTH),ierr)

    ! create communicators for every column and row of the cartesian
    !	topology (used eg. for fargo shifts)
    remain_dims = (/ .FALSE., .TRUE. /)
    CALL MPI_Cart_Sub(Mesh%comm_cart,remain_dims,Mesh%Icomm,ierr)
    remain_dims = (/ .TRUE., .FALSE. /)
    CALL MPI_Cart_Sub(Mesh%comm_cart,remain_dims,Mesh%Jcomm,ierr)

    ! allocate memory for boundary data buffers
    ALLOCATE(this(WEST)%sendbuf(Mesh%GNUM,Mesh%JMIN:Mesh%JMAX,Physics%VNUM), &
         this(WEST)%recvbuf(Mesh%GNUM,Mesh%JMIN:Mesh%JMAX,Physics%VNUM), &
         this(EAST)%sendbuf(Mesh%GNUM,Mesh%JMIN:Mesh%JMAX,Physics%VNUM), &
         this(EAST)%recvbuf(Mesh%GNUM,Mesh%JMIN:Mesh%JMAX,Physics%VNUM), &
         this(SOUTH)%sendbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%GNUM,Physics%VNUM), &
         this(SOUTH)%recvbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%GNUM,Physics%VNUM), &
         this(NORTH)%sendbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%GNUM,Physics%VNUM), &
         this(NORTH)%recvbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%GNUM,Physics%VNUM), &
         STAT=ierr)
    IF (ierr.NE.0) THEN
       CALL Error(this(WEST),"InitBoundary", &
            "Unable to allocate memory for data buffers.")
    END IF
    DO dir=WEST,NORTH
      this(dir)%recvbuf = 0.
      this(dir)%sendbuf = 0.
    END DO
#endif

  END SUBROUTINE InitBoundary


  SUBROUTINE CenterBoundary(this,Mesh,Fluxes,Physics,time,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this(4)
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Fluxes_TYP)   :: Fluxes
    TYPE(Physics_TYP)  :: Physics
    REAL               :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                       :: pvar, cvar
    !------------------------------------------------------------------------!
    INTEGER       :: i,j
#ifdef PARALLEL
    INTEGER       :: req(4)
    INTEGER       :: ierr
#ifdef MPI_USE_SENDRECV
    INTEGER       :: status(MPI_STATUS_SIZE)
#else
    INTEGER       :: status(MPI_STATUS_SIZE,4)
#endif
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,Fluxes,time
    INTENT(INOUT) :: this,pvar,cvar
    !------------------------------------------------------------------------!
    CALL Convert2Primitive(Physics,Mesh,Mesh%IMIN,Mesh%IMAX,Mesh%JMIN, &
             Mesh%JMAX,cvar,pvar)
#ifdef PARALLEL
    ! NOTE: if you want to use MPI_Sendrecv instead of nonblocking
    ! MPI_Irecv and  MPI_Issend for exchange of ghost cell data, 
    ! you must add -DMPI_USE_SENDRECV to the compile command

    ! initiate western/eastern MPI communication
#ifdef MPI_USE_SENDRECV
    ! send boundary data to western and receive from eastern neighbor
    IF (Mesh%neighbor(WEST).NE.MPI_PROC_NULL) &
         this(WEST)%sendbuf(:,:,:) = pvar(Mesh%IMIN:Mesh%IMIN+Mesh%GNUM-1, &
                                          Mesh%JMIN:Mesh%JMAX,1:Physics%VNUM)
    CALL MPI_Sendrecv(this(WEST)%sendbuf,Mesh%GNUM*(Mesh%JMAX-Mesh%JMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(WEST),10+WEST,this(EAST)%recvbuf, &
         Mesh%GNUM*(Mesh%JMAX-Mesh%JMIN+1)*Physics%VNUM,DEFAULT_MPI_REAL,Mesh%neighbor(EAST), &
         MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
    IF (Mesh%neighbor(EAST).NE.MPI_PROC_NULL) &
         pvar(Mesh%IMAX+1:Mesh%IGMAX,Mesh%JMIN:Mesh%JMAX,1:Physics%VNUM) = this(EAST)%recvbuf(:,:,:)
#else
    ! receive boundary data from eastern neighbor
    CALL MPI_Irecv(this(EAST)%recvbuf,Mesh%GNUM*(Mesh%JMAX-Mesh%JMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(EAST),10+WEST,Mesh%comm_cart,req(1),ierr)
    ! fill send buffer if western neighbor exists
    IF (Mesh%neighbor(WEST).NE.MPI_PROC_NULL) &
         this(WEST)%sendbuf(:,:,:) = pvar(Mesh%IMIN:Mesh%IMIN+Mesh%GNUM-1, &
                                          Mesh%JMIN:Mesh%JMAX,1:Physics%VNUM)
    ! send boundary data to western neighbor
    CALL MPI_Issend(this(WEST)%sendbuf,Mesh%GNUM*(Mesh%JMAX-Mesh%JMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(WEST),10+WEST,Mesh%comm_cart,req(2),ierr)
#endif
#endif
    ! set physical boundary conditions at western boundaries
    IF (Mesh%INUM.GT.1) CALL SetBoundaryData(WEST)
#ifdef PARALLEL
#ifdef MPI_USE_SENDRECV
    ! send boundary data to eastern and receive from western neighbor
    IF (Mesh%neighbor(EAST).NE.MPI_PROC_NULL) &
         this(EAST)%sendbuf(:,:,:) = pvar(Mesh%IMAX-Mesh%GNUM+1:Mesh%IMAX, &
                                          Mesh%JMIN:Mesh%JMAX,1:Physics%VNUM)
    CALL MPI_Sendrecv(this(EAST)%sendbuf,Mesh%GNUM*(Mesh%JMAX-Mesh%JMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(EAST),10+EAST,this(WEST)%recvbuf, &
         Mesh%GNUM*(Mesh%JMAX-Mesh%JMIN+1)*Physics%VNUM,DEFAULT_MPI_REAL,Mesh%neighbor(WEST), &
         MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
    IF (Mesh%neighbor(WEST).NE.MPI_PROC_NULL) &
         pvar(Mesh%IGMIN:Mesh%IMIN-1,Mesh%JMIN:Mesh%JMAX,1:Physics%VNUM) = this(WEST)%recvbuf(:,:,:)
#else
    ! receive boundary data from western neighbor
    CALL MPI_Irecv(this(WEST)%recvbuf,Mesh%GNUM*(Mesh%JMAX-Mesh%JMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(WEST),10+EAST,Mesh%comm_cart,req(3),ierr)
    ! fill send buffer if eastern neighbor exists
    IF (Mesh%neighbor(EAST).NE.MPI_PROC_NULL) &
         this(EAST)%sendbuf(:,:,:) = pvar(Mesh%IMAX-Mesh%GNUM+1:Mesh%IMAX, &
                                          Mesh%JMIN:Mesh%JMAX,1:Physics%VNUM)
    ! send boundary data to eastern neighbor
    CALL MPI_Issend(this(EAST)%sendbuf,Mesh%GNUM*(Mesh%JMAX-Mesh%JMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(EAST),10+EAST,Mesh%comm_cart,req(4),ierr)
#endif
#endif
    ! set physical boundary conditions at eastern boundaries
    IF (Mesh%INUM.GT.1) CALL SetBoundaryData(EAST)
#ifdef PARALLEL
#ifndef MPI_USE_SENDRECV
    ! wait for unfinished MPI communication
    CALL MPI_Waitall(4,req,status,ierr)
    ! copy data from recieve buffers into ghosts cells
    IF (Mesh%neighbor(WEST).NE.MPI_PROC_NULL) &
         pvar(Mesh%IGMIN:Mesh%IMIN-1,Mesh%JMIN:Mesh%JMAX,1:Physics%VNUM) = this(WEST)%recvbuf(:,:,:)
    IF (Mesh%neighbor(EAST).NE.MPI_PROC_NULL) &
         pvar(Mesh%IMAX+1:Mesh%IGMAX,Mesh%JMIN:Mesh%JMAX,1:Physics%VNUM) = this(EAST)%recvbuf(:,:,:)
#endif
#endif

#ifdef PARALLEL
    ! initiate southern/northern MPI communication
#ifdef MPI_USE_SENDRECV
    ! send boundary data to southern and receive from northern neighbor
    IF (Mesh%neighbor(SOUTH).NE.MPI_PROC_NULL) &
         this(SOUTH)%sendbuf(:,:,:) = pvar(Mesh%IGMIN:Mesh%IGMAX, &
                                           Mesh%JMIN:Mesh%JMIN+Mesh%GNUM-1,1:Physics%VNUM)
    CALL MPI_Sendrecv(this(SOUTH)%sendbuf,Mesh%GNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(SOUTH),10+SOUTH,this(NORTH)%recvbuf, &
         Mesh%GNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*Physics%VNUM,DEFAULT_MPI_REAL,Mesh%neighbor(NORTH), &
         MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
    IF (Mesh%neighbor(NORTH).NE.MPI_PROC_NULL) &
         pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMAX+1:Mesh%JGMAX,1:Physics%VNUM) = this(NORTH)%recvbuf(:,:,:)
#else
    ! receive boundary data from northern neighbor
    CALL MPI_Irecv(this(NORTH)%recvbuf,Mesh%GNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*Physics%vnum, &
         DEFAULT_MPI_REAL,Mesh%neighbor(NORTH),10+SOUTH,Mesh%comm_cart,req(1),ierr)
    ! fill send buffer if southern neighbor exists
    IF (Mesh%neighbor(SOUTH).NE.MPI_PROC_NULL) &
         this(SOUTH)%sendbuf(:,:,:) = pvar(Mesh%IGMIN:Mesh%IGMAX, &
                                           Mesh%JMIN:Mesh%JMIN+Mesh%GNUM-1,1:Physics%VNUM)
    ! send boundary data to southern neighbor
    CALL MPI_Issend(this(SOUTH)%sendbuf,Mesh%GNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*Physics%vnum, &
         DEFAULT_MPI_REAL,Mesh%neighbor(SOUTH),10+SOUTH,Mesh%comm_cart,req(2),ierr)
#endif
#endif
    ! set physical boundary conditions at southern boundaries
    IF (Mesh%JNUM.GT.1) CALL SetBoundaryData(SOUTH)
#ifdef PARALLEL
#ifdef MPI_USE_SENDRECV
    ! send boundary data to northern and receive from southern neighbor
    IF (Mesh%neighbor(NORTH).NE.MPI_PROC_NULL) &
         this(NORTH)%sendbuf(:,:,:) = pvar(Mesh%IGMIN:Mesh%IGMAX, &
                                           Mesh%JMAX-Mesh%GNUM+1:Mesh%JMAX,1:Physics%VNUM)

    CALL MPI_Sendrecv(this(NORTH)%sendbuf,Mesh%GNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(NORTH),10+NORTH,this(SOUTH)%recvbuf, &
         Mesh%GNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*Physics%VNUM,DEFAULT_MPI_REAL,Mesh%neighbor(SOUTH), &
         MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
    IF (Mesh%neighbor(SOUTH).NE.MPI_PROC_NULL) &
         pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JMIN-1,1:Physics%VNUM) = this(SOUTH)%recvbuf(:,:,:)
#else
    ! receive boundary data from southern neighbor
    CALL MPI_Irecv(this(SOUTH)%recvbuf,Mesh%GNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*Physics%vnum, &
         DEFAULT_MPI_REAL,Mesh%neighbor(SOUTH),10+NORTH,Mesh%comm_cart,req(3),ierr)
    ! fill send buffer if northern neighbor exists
    IF (Mesh%neighbor(NORTH).NE.MPI_PROC_NULL) &
         this(NORTH)%sendbuf(:,:,:) = pvar(Mesh%IGMIN:Mesh%IGMAX, &
                                           Mesh%JMAX-Mesh%GNUM+1:Mesh%JMAX,1:Physics%VNUM)
    ! send boundary data to northern neighbor
    CALL MPI_Issend(this(NORTH)%sendbuf,Mesh%GNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*Physics%vnum, &
         DEFAULT_MPI_REAL,Mesh%neighbor(NORTH),10+NORTH,Mesh%comm_cart,req(4),ierr)
#endif
#endif
    ! set physical boundary conditions at northern boundaries
    IF (Mesh%JNUM.GT.1) CALL SetBoundaryData(NORTH)
#ifdef PARALLEL
#ifndef MPI_USE_SENDRECV
    ! wait for unfinished MPI communication
    CALL MPI_Waitall(4,req,status,ierr)
    IF (Mesh%neighbor(SOUTH).NE.MPI_PROC_NULL) &
         pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JMIN-1,1:Physics%VNUM) = this(SOUTH)%recvbuf(:,:,:)
    IF (Mesh%neighbor(NORTH).NE.MPI_PROC_NULL) &
         pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMAX+1:Mesh%JGMAX,1:Physics%VNUM) = this(NORTH)%recvbuf(:,:,:)
#endif
#endif
    ! FIXME: implement MPI communication to exchange corner values
    ! this is a quick hack to set defined boundary values in
    ! the corners outside the computational domain;
    ! this is also necessary, because we need some of these values in the
    ! viscosity module
    ! Check if it is a real corner
    IF (Mesh%JNUM.EQ.1) THEN
      ! for 1D simulations along the x-direction copy internal data
      ! into the ghost cells at the y-boundaries
      DO j=1,Mesh%GNUM
        pvar(:,Mesh%JMIN-j,:) = pvar(:,Mesh%JMIN,:)
        pvar(:,Mesh%JMAX+j,:) = pvar(:,Mesh%JMAX,:)
      END DO
    ELSE IF (Mesh%INUM.EQ.1) THEN
      ! for 1D simulations along the y-direction copy internal data
      ! into the ghost cells at the x-boundaries
      DO i=1,Mesh%GNUM
        pvar(Mesh%IMIN-i,:,:) = pvar(Mesh%IMIN,:,:)
        pvar(Mesh%IMAX+i,:,:) = pvar(Mesh%IMAX,:,:)
      END DO
    ELSE IF ((Mesh%INUM.GT.1) .AND. (Mesh%JNUM.GT.1) .AND.&
      (.NOT.((GetType(this(NORTH)).EQ.PERIODIC).OR.(GetType(this(SOUTH)).EQ.PERIODIC) &
#ifdef PARALLEL
         .OR.(GetType(this(NORTH)).EQ.NONE).OR.(GetType(this(SOUTH)).EQ.NONE) &
#endif
      ))) THEN
      DO j=1,Mesh%GNUM
        DO i=j+1,Mesh%GNUM
          ! copy data left of the diagonal within the corner areas at
          ! the western boundary
          ! south west
          pvar(Mesh%IMIN-i,Mesh%JMIN-j,:) = pvar(Mesh%IMIN-i,Mesh%JMIN-j+1,:)
          ! north west
          pvar(Mesh%IMIN-i,Mesh%JMAX+j,:) = pvar(Mesh%IMIN-i,Mesh%JMAX+j-1,:)
          ! copy data right of the diagonal within the corner areas at
          ! the eastern boundary
          ! south east
          pvar(Mesh%IMAX+i,Mesh%JMIN-j,:) = pvar(Mesh%IMAX+i,Mesh%JMIN-j+1,:)
          ! north east
          pvar(Mesh%IMAX+i,Mesh%JMAX+j,:) = pvar(Mesh%IMAX+i,Mesh%JMAX+j-1,:)

          ! copy data below the diagonal within the corner areas at
          ! the southern boundary (transposed problem i<->j)
          ! south west
          pvar(Mesh%IMIN-j,Mesh%JMIN-i,:) = pvar(Mesh%IMIN-j+1,Mesh%JMIN-i,:)
          ! south east
          pvar(Mesh%IMAX+j,Mesh%JMIN-i,:) = pvar(Mesh%IMAX+j-1,Mesh%JMIN-i,:)
          ! copy data above the diagonal within the corner areas at
          ! the northern boundary (transposed problem i<->j)
          ! north west
          pvar(Mesh%IMIN-j,Mesh%JMAX+i,:) = pvar(Mesh%IMIN-j+1,Mesh%JMAX+i,:)
          ! north east
          pvar(Mesh%IMAX+j,Mesh%JMAX+i,:) = pvar(Mesh%IMAX+j-1,Mesh%JMAX+i,:)
        END DO
      END DO
      ! set the diagonal data 
      DO i=1,Mesh%GNUM
        ! south west
        pvar(Mesh%IMIN-i,Mesh%JMIN-i,:) = 0.5 * (pvar(Mesh%IMIN-i,Mesh%JMIN,:) &
             + pvar(Mesh%IMIN,Mesh%JMIN-i,:))
        ! south east
        pvar(Mesh%IMAX+i,Mesh%JMIN-i,:) = 0.5 * (pvar(Mesh%IMAX+i,Mesh%JMIN,:) &
             + pvar(Mesh%IMAX,Mesh%JMIN-i,:))
        ! north west
        pvar(Mesh%IMIN-i,Mesh%JMAX+i,:) = 0.5 * (pvar(Mesh%IMIN-i,Mesh%JMAX,:) &
             + pvar(Mesh%IMIN,Mesh%JMAX+i,:))
        ! north east
        pvar(Mesh%IMAX+i,Mesh%JMAX+i,:) = 0.5 * (pvar(Mesh%IMAX+i,Mesh%JMAX,:) &
             + pvar(Mesh%IMAX,Mesh%JMAX+i,:))
      END DO
    END IF
    ! convert primitive variables in ghost cells
    CALL Convert2Conservative(Physics,Mesh,Mesh%IGMIN,Mesh%IMIN-1,&
         Mesh%JGMIN,Mesh%JGMAX,pvar,cvar)
    CALL Convert2Conservative(Physics,Mesh,Mesh%IMAX+1,Mesh%IGMAX,&
         Mesh%JGMIN,Mesh%JGMAX,pvar,cvar)
    CALL Convert2Conservative(Physics,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
         Mesh%JGMIN,Mesh%JMIN-1,pvar,cvar)
    CALL Convert2Conservative(Physics,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
         Mesh%JMAX+1,Mesh%JGMAX,pvar,cvar)

  CONTAINS

    ! set boundary data for (real) physical boundaries in direction "dir"
    SUBROUTINE SetBoundaryData(dir)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: dir
!CDIR IEXPAND
      SELECT CASE(GetType(this(dir)))
      CASE(NO_GRADIENTS)
         CALL CenterBoundary_nogradients(this(dir),Mesh,Physics,pvar)
      CASE(PERIODIC)
         ! do nothing in parallel version, because periodicity is
         ! handled via MPI communication
#ifndef PARALLEL
         CALL CenterBoundary_periodic(this(dir),Mesh,Physics,pvar)
#endif
      CASE(REFLECTING)
         CALL CenterBoundary_reflecting(this(dir),Mesh,Physics,pvar)
      CASE(AXIS)
         CALL CenterBoundary_axis(this(dir),Mesh,Physics,pvar)
      CASE(FOLDED)
         CALL CenterBoundary_folded(this(dir),Mesh,Physics,pvar)
      CASE(FIXED)
         CALL CenterBoundary_fixed(this(dir),Mesh,Physics,pvar)
      CASE(EXTRAPOLATION)
         CALL CenterBoundary_extrapolation(this(dir),Mesh,Physics,pvar)
      CASE(NOH2D,NOH3D)
         CALL CenterBoundary_noh(this(dir),Mesh,Physics,time,pvar)
      CASE(NOSLIP)
         CALL CenterBoundary_noslip(this(dir),Mesh,Physics,pvar)
      CASE(CUSTOM)
         CALL CenterBoundary_custom(this(dir),Mesh,Physics,pvar)
      CASE(FARFIELD)
         CALL CenterBoundary_farfield(this(dir),Mesh,Physics,pvar)
      CASE(ABSORBING)
         CALL CenterBoundary_absorbing(this(dir),Mesh,Physics,cvar,pvar)
      END SELECT
    END SUBROUTINE SetBoundaryData

  END SUBROUTINE CenterBoundary


  SUBROUTINE CloseBoundary_one(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF (.NOT.Initialized(this)) &
        CALL Error(this,"CloseBoundary_one","not initialized")
#ifdef PARALLEL
    ! deallocate MPI send/recv buffers 
    DEALLOCATE(this%sendbuf,this%recvbuf)
#endif
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(NO_GRADIENTS,PERIODIC,EXTRAPOLATION)
       ! do nothing
    CASE(REFLECTING)
       CALL CloseBoundary_reflecting(this)
    CASE(AXIS)
       CALL CloseBoundary_axis(this)
    CASE(FOLDED)
       CALL CloseBoundary_folded(this)
    CASE(FIXED)
       CALL CloseBoundary_fixed(this)
    CASE(NOSLIP)
       CALL CloseBoundary_noslip(this)
    CASE(CUSTOM)
       CALL CloseBoundary_custom(this)
    CASE(FARFIELD)
       CALL CloseBoundary_farfield(this)
    CASE(ABSORBING)
       CALL CloseBoundary_absorbing(this)
    END SELECT
    CALL CloseBoundary_common(this)
  END SUBROUTINE CloseBoundary_one


  SUBROUTINE CloseBoundary_all(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), DIMENSION(4) :: this
    !------------------------------------------------------------------------!
    INTEGER       :: dir
    !------------------------------------------------------------------------!
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    ! loop over all boundaries
    DO dir=1,4
       CALL CloseBoundary_one(this(dir))
    END DO
  END SUBROUTINE CloseBoundary_all

END MODULE boundary_generic
