!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# program file: main.f90                                                    #
!#                                                                           #
!# Copyright (C) 2006 - 2010                                                 #
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
! Main program file
!----------------------------------------------------------------------------!
PROGRAM fosite
  USE physics_generic
  USE sources_generic, CalcTimestep_sources => CalcTimestep
  USE mesh_generic
  USE fluxes_generic
  USE boundary_generic
  USE fileio_generic
  USE reconstruction_generic
  USE timedisc_generic
  USE init
  USE integration
  IMPLICIT NONE
#ifdef PARALLEL
    include 'mpif.h'
#endif
  !--------------------------------------------------------------------------!
  TYPE(Mesh_TYP)       :: Mesh
  TYPE(Fluxes_TYP)     :: Fluxes
  TYPE(Physics_TYP)    :: Physics
  TYPE(FileIO_TYP)     :: Datafile
  TYPE(Timedisc_TYP)   :: Timedisc
  TYPE(FileIO_TYP)     :: Logfile
  
  INTEGER              :: n
  INTEGER              :: myrank
  DOUBLE PRECISION     :: wall_time           ! wall clock elapsed time      !
  DOUBLE PRECISION     :: log_time            ! time for next log output     !
  DOUBLE PRECISION     :: start_time          ! system clock start time      !
  DOUBLE PRECISION     :: end_time            ! system clock end time        !
  DOUBLE PRECISION     :: run_time            ! = end_time - start_time      !
#ifdef PARALLEL
  INTEGER              :: ierror
  REAL                 :: dt_all              ! min timestep of all processes!
#endif
  !--------------------------------------------------------------------------!

#ifdef PARALLEL
  ! initialize MPI library for parallel execution
  CALL MPI_Init(ierror)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierror)     
#else
    myrank = 0
#endif

  CALL InitIntegration

  IF (myrank.EQ.0) THEN
     ! print some information
     PRINT "(A)", "+---------------------------------------------------------+"
     PRINT "(A)", "|          Solution of 2D advection problems              |"
     PRINT "(A)", "+---------------------------------------------------------+"
     PRINT *, "Initializing simulation:"
  END IF
  
  ! setup simulation
  CALL InitProgram(Mesh,Physics,Fluxes,Timedisc,Datafile,Logfile)

  ! allocate memory for physics and fluxes modules
  CALL MallocPhysics(Physics,Mesh)
  CALL MallocFluxes(Fluxes,Mesh,Physics)

#ifdef PARALLEL
  wall_time = MPI_Wtime()
  CALL MPI_AllReduce(wall_time,start_time,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
       Mesh%comm_cart,ierror)
#else
  CALL CPU_TIME(start_time)
#endif

  ! make sure that the initial data is written to the log file
  wall_time = start_time
  log_time  = wall_time

  IF (myrank.EQ.0) THEN
     PRINT *, "==================================================================="
     PRINT *, "Starting calculation..."
  END IF

  ! store initial data
  IF (Timedisc%time.EQ.0.0) THEN
     CALL WriteDataset(Datafile,Mesh,Physics,Fluxes,Timedisc)
     IF (myrank.EQ.0) CALL PrintInfo(0,0.0,0.0,Timedisc%n_adj)
  END IF

  ! main loop
  DO n=1,Timedisc%maxiter
     ! finish simulation if stop time is reached
     IF (ABS(Timedisc%stoptime-Timedisc%time).LE.1.0E-05*Timedisc%stoptime) EXIT

     ! calculate timestep
     CALL CalcTimestep(Timedisc,Mesh,Physics,Fluxes)

     ! adjust timestep for output
     ! and calculate the wall clock time
     CALL AdjustTimestep(Datafile,Timedisc%time,Timedisc%dt)
#ifdef PARALLEL
     CALL MPI_Allreduce(Timedisc%dt,dt_all,1,DEFAULT_MPI_REAL,MPI_MIN,&
          Mesh%comm_cart,ierror)
     Timedisc%dt = dt_all
     run_time = MPI_Wtime()
     CALL MPI_Allreduce(run_time,wall_time,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
          Mesh%comm_cart,ierror)
#else
     CALL CPU_TIME(wall_time)
#endif

     ! advance the solution in time
     CALL SolveODE(Timedisc,Mesh,Physics,Fluxes)

     ! write output to data file
     IF (ABS(Datafile%time-Timedisc%time).LE.1.0E-5*Datafile%time) THEN
        CALL WriteDataset(Datafile,Mesh,Physics,Fluxes,Timedisc)
        CALL PrintInfo(n,Timedisc%time,Timedisc%dtmin,Timedisc%n_adj)
        ! reset dt_min and n_adj
        Timedisc%dtmin = Timedisc%stoptime
        Timedisc%n_adj = 0
     END IF

     ! write output to log file
     IF (Initialized(Logfile).AND.(wall_time.GE.log_time)) THEN
        CALL WriteDataset(Logfile,Mesh,Physics,Fluxes,Timedisc)
        log_time = wall_time + Logfile%dtwall
     END IF

  END DO

#ifdef PARALLEL
  run_time = MPI_Wtime()
  CALL MPI_AllReduce(run_time,end_time,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       Mesh%comm_cart,ierror)
#else
  CALL CPU_TIME(end_time)
#endif

  CALL PrintBoundaryFluxes(Physics)
  CALL PrintSummary

  CALL CloseFileIO(Datafile)
  CALL CloseFileIO(Logfile)
  CALL CloseTimedisc(Timedisc)

  IF (ASSOCIATED(Physics%sources)) CALL CloseSources(Physics%sources,Fluxes)
  CALL ClosePhysics(Physics)
  CALL CloseMesh(Mesh,Fluxes)
  CALL CloseFluxes(Fluxes)
  CALL CloseIntegration

#ifdef PARALLEL
  CALL MPI_Finalize(ierror)
#endif

CONTAINS

  SUBROUTINE PrintInfo(i,t,d,na)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER     :: i, na
    REAL        :: t, d
    !------------------------------------------------------------------------!
    INTENT(IN)  :: i,t,d
    !------------------------------------------------------------------------!
    IF (myrank.EQ.0) THEN
       PRINT "(A,I8,A,ES11.3,A,ES11.3,A,I5)", " Iteration ", i, &
            "  time ", t, "  min dt ", d, "  adj ", na
    END IF
  END SUBROUTINE PrintInfo

  SUBROUTINE PrintBoundaryFluxes(Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    INTEGER :: k
    REAL, DIMENSION(Physics%VNUM,4) :: bflux
    !------------------------------------------------------------------------!
    DO k=1,4
       bflux(:,k) = GetBoundaryFlux(Fluxes,Mesh,Physics,k)
    END DO
    IF (myrank.EQ.0) THEN
       PRINT *, "-------------------------------------------------------------------"
       PRINT *, "total boundary fluxes:"
       PRINT *, "                      west        east        south       north"
       DO k=1,Physics%VNUM
          PRINT "(T2,A,T21,4(ES12.3))", TRIM(Physics%cvarname(k)), &
               bflux(k,WEST), bflux(k,EAST), bflux(k,SOUTH), bflux(k,NORTH)   
       END DO
    END IF
  END SUBROUTINE PrintBoundaryFluxes

  SUBROUTINE PrintSummary
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    IF (myrank.EQ.0) THEN
       PRINT *, "==================================================================="
       IF (n.LT.Timedisc%maxiter) THEN
          PRINT *, "calculation finished correctly."
       ELSE
          PRINT *, "too many iterations, aborting!"
       END IF
       run_time = end_time - start_time
       PRINT "(A,F10.2,A)", " main loop runtime: ", run_time, " sec."
    END IF
  END SUBROUTINE PrintSummary

END PROGRAM fosite
