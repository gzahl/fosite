!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# program file: fosite.f90                                                  #
!#                                                                           #
!# Copyright (C) 2006-2012                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
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
! Module fosite
!----------------------------------------------------------------------------!
MODULE fosite
  USE fosite_common
  USE sources_generic, CalcTimestep_sources => CalcTimestep
  USE physics_generic, ONLY : InitPhysics, ClosePhysics
  USE fluxes_generic, ONLY : InitFluxes, CloseFluxes, GetBoundaryFlux
  USE reconstruction_generic, ONLY : InitReconstruction
  USE mesh_generic, ONLY : InitMesh, CloseMesh
  USE boundary_generic, ONLY : InitBoundary, WEST, EAST, SOUTH, NORTH
  USE sources_generic, ONLY : InitSources, CloseSources
  USE timedisc_generic, ONLY : InitTimedisc, CloseTimedisc, SolveODE, &
                               CalcTimestep
  USE fileio_generic, ONLY : InitFileIO, CloseFileIO, Initialized, &
                             WriteDataset, AdjustTimestep

  USE integration
#ifdef PARALLEL
  USE common_types
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
  INTEGER, PARAMETER :: simtype = 1
  CHARACTER(LEN=32), PARAMETER :: simname = "fosite"  
  INTEGER, PARAMETER :: MAXLEN = 500
  CHARACTER(MAXLEN)  :: buffer
  !--------------------------------------------------------------------------!
  PUBLIC               :: Fosite_TYP, &
                          InitFosite, &
                          RunFosite, &
                          StepFosite, &
                          CloseFosite, &
                          InitPhysics, &
                          InitFluxes, &
                          InitReconstruction, &
                          InitMesh, &
                          InitBoundary, &
                          InitSources, &
                          InitTimedisc, &
                          InitFileIO
  !--------------------------------------------------------------------------!


CONTAINS

  SUBROUTINE InitFosite(this)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Fosite_TYP)   :: this
    !--------------------------------------------------------------------------!
    INTENT(INOUT)      :: this
    !--------------------------------------------------------------------------!
#ifdef PARALLEL
    ! initialize MPI library for parallel execution
    CALL MPI_Init(this%ierror)
#endif
    CALL InitSim(this, simtype, simname)

    IF (GetRank(this).EQ.0) THEN
        ! print some information
        WRITE(buffer, "(A)")&
            "+---------------------------------------------------------+"
        CALL Info(this, buffer)
        WRITE(buffer, "(A)")&
            "|          Solution of 2D advection problems              |"
        CALL Info(this, buffer)
        WRITE(buffer, "(A)")&
            "+---------------------------------------------------------+"
        CALL Info(this, buffer)
        CALL Info(this, "Initializing simulation:")
    END IF

    this%iter = 1
  
  END SUBROUTINE InitFosite



  SUBROUTINE FirstStepFosite(this)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Fosite_TYP)   :: this
    !--------------------------------------------------------------------------!
    INTENT(INOUT)      :: this
    !--------------------------------------------------------------------------!

#ifdef PARALLEL
    this%wall_time = MPI_Wtime()
    CALL MPI_Allreduce(this%wall_time,this%start_time,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
        this%Mesh%comm_cart,this%ierror)
#else
    CALL CPU_TIME(this%start_time)
#endif

    ! make sure that the initial data is written to the log file
    this%wall_time = this%start_time
    this%log_time  = this%wall_time

    IF (GetRank(this).EQ.0) THEN
        CALL Info(this,&
            "===================================================================")
        CALL Info(this, "Starting calculation...")
    END IF

    ! store old values
    this%Timedisc%cold(:,:,:) = this%Timedisc%cvar(:,:,:)
    this%Timedisc%pold(:,:,:) = this%Timedisc%pvar(:,:,:)

    ! store initial data
    IF (this%Timedisc%time.EQ.0.0) THEN
        CALL WriteDataset(this%Datafile,this%Mesh,this%Physics,this%Fluxes,&
                          this%Timedisc)
        IF (GetRank(this).EQ.0) CALL PrintInfo(this,0,0.0,0.0,this%Timedisc%n_adj)
    END IF

  END SUBROUTINE FirstStepFosite

  SUBROUTINE RunFosite(this)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Fosite_TYP)   :: this
    !--------------------------------------------------------------------------!
    INTENT(INOUT)      :: this
    !--------------------------------------------------------------------------!

    ! main loop
    DO WHILE(this%iter.LE.this%Timedisc%maxiter)
        IF(StepFosite(this)) EXIT
    END DO

  END SUBROUTINE RunFosite
     
  FUNCTION StepFosite(this) RESULT(break)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Fosite_TYP)   :: this
    LOGICAL            :: break
    !--------------------------------------------------------------------------!
    INTENT(INOUT)      :: this
    !--------------------------------------------------------------------------!

    break = .FALSE.

    IF (this%iter == 1) &
        CALL FirstStepFosite(this)

    ! finish simulation if stop time is reached
    IF (ABS(this%Timedisc%stoptime-this%Timedisc%time)&
            .LE.1.0E-05*this%Timedisc%stoptime) THEN
            break = .TRUE.
            RETURN
    END IF

    ! calculate timestep
    CALL CalcTimestep(this%Timedisc,this%Mesh,this%Physics,this%Fluxes)

    ! adjust timestep for output
    ! and calculate the wall clock time
    CALL AdjustTimestep(this%Datafile,this%Timedisc%time,this%Timedisc%dt)
#ifdef PARALLEL
    CALL MPI_Allreduce(this%Timedisc%dt,this%dt_all,1,DEFAULT_MPI_REAL,MPI_MIN,&
         this%Mesh%comm_cart,this%ierror)
    this%Timedisc%dt = this%dt_all
    this%run_time = MPI_Wtime()
    CALL MPI_Allreduce(this%run_time,this%wall_time,1,MPI_DOUBLE_PRECISION,&
                       MPI_MIN,this%Mesh%comm_cart,this%ierror)
#else
    CALL CPU_TIME(this%wall_time)
#endif

    ! advance the solution in time
    CALL SolveODE(this%Timedisc,this%Mesh,this%Physics,this%Fluxes)

    ! write output to data file
    IF (ABS(this%Datafile%time-this%Timedisc%time)&
            .LE.1.0E-5*this%Datafile%time) THEN
       CALL WriteDataset(this%Datafile,this%Mesh,this%Physics,this%Fluxes,&
                         this%Timedisc)
       CALL PrintInfo(this, this%iter,this%Timedisc%time,this%Timedisc%dtmin,&
                      this%Timedisc%n_adj)
       ! reset dt_min and n_adj
       this%Timedisc%dtmin = this%Timedisc%stoptime
       this%Timedisc%n_adj = 0
    END IF

    ! write output to log file
    IF (Initialized(this%Logfile).AND.(this%wall_time.GE.this%log_time)) THEN
       CALL WriteDataset(this%Logfile,this%Mesh,this%Physics,this%Fluxes,&
                         this%Timedisc)
       this%log_time = this%wall_time + this%Logfile%dtwall
    END IF

    this%iter = this%iter + 1

  END FUNCTION StepFosite

  SUBROUTINE CloseFosite(this)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Fosite_TYP)   :: this
    !--------------------------------------------------------------------------!
    INTENT(INOUT)      :: this
    !--------------------------------------------------------------------------!

#ifdef PARALLEL
    this%run_time = MPI_Wtime()
    CALL MPI_Allreduce(this%run_time,this%end_time,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
        this%Mesh%comm_cart,this%ierror)
#else
    CALL CPU_TIME(this%end_time)
#endif

    CALL PrintBoundaryFluxes(this)
    CALL PrintSummary(this)

    CALL CloseFileIO(this%Datafile)
    IF (Initialized(this%Logfile)) CALL CloseFileIO(this%Logfile)
    CALL CloseTimedisc(this%Timedisc)

    IF (ASSOCIATED(this%Physics%sources)) &
        CALL CloseSources(this%Physics%sources,this%Fluxes)
    CALL ClosePhysics(this%Physics)
    CALL CloseFluxes(this%Fluxes)
    CALL CloseMesh(this%Mesh)

#ifdef PARALLEL
    CALL MPI_Finalize(this%ierror)
#endif

    CALL CloseSim(this)
  END SUBROUTINE CloseFosite


  SUBROUTINE PrintInfo(this,i,t,d,na)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)   :: this
    INTEGER     :: i, na, datetime(8)
    REAL        :: t, d
    !------------------------------------------------------------------------!
    INTENT(IN)  :: i,t,d
    !------------------------------------------------------------------------!
    INTENT(INOUT)      :: this
    !--------------------------------------------------------------------------!
    IF (GetRank(this).EQ.0) THEN
       CALL date_and_time(values = datetime)
       !CALL IDATE(today)
       !CALL ITIME(now)
       WRITE(buffer,"(I2.2,A,I2.2,A,I4.4,A,I2.2,A,I2.2,A,I2.2,A,I8,A,ES11.3,A,ES11.3,A,I5)")&
             datetime(3), ".", datetime(2), ".", datetime(1), " ",&
             datetime(5), ":", datetime(6), ":", datetime(7),&
             " n ", i, "  time ", t, "  min dt ", d, "  adj ", na
       CALL Info(this, buffer)
    END IF
  END SUBROUTINE PrintInfo

  SUBROUTINE PrintBoundaryFluxes(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)  :: this
    INTEGER :: k
    REAL, DIMENSION(this%Physics%VNUM,4) :: bflux
    !------------------------------------------------------------------------!
    INTENT(INOUT)      :: this
    !--------------------------------------------------------------------------!
    DO k=1,4
       bflux(:,k) = GetBoundaryFlux(this%Fluxes,this%Mesh,this%Physics,k)
    END DO
    IF (GetRank(this).EQ.0) THEN
       CALL Info(this, "-------------------------------------------------------------------")
       CALL Info(this, "total boundary fluxes:")
       CALL Info(this, "                      west        east        south       north")
       DO k=1,this%Physics%VNUM
          WRITE(buffer,"(T2,A,T21,4(ES12.3))")TRIM(this%Physics%cvarname(k)), &
               bflux(k,WEST), bflux(k,EAST), bflux(k,SOUTH), bflux(k,NORTH)   
          CALL Info(this, buffer)
       END DO
    END IF
  END SUBROUTINE PrintBoundaryFluxes

  SUBROUTINE PrintSummary(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)   :: this
    !--------------------------------------------------------------------------!
    INTENT(INOUT)      :: this
    !--------------------------------------------------------------------------!
    IF (GetRank(this).EQ.0) THEN
       CALL Info(this, "===================================================================")
       IF (this%iter.LT.this%Timedisc%maxiter) THEN
          CALL Info(this, "calculation finished correctly.")
       ELSE
          CALL Info(this, "too many iterations, aborting!")
       END IF
       this%run_time = this%end_time - this%start_time
       WRITE(buffer,"(A,F10.2,A)")" main loop runtime: ", this%run_time, " sec."
       CALL Info(this, buffer)
    END IF
  END SUBROUTINE PrintSummary

END MODULE fosite
