!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# program file: fosite.f90                                                  #
!#                                                                           #
!# Copyright (C) 2006-2014                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Manuel Jung      <mjung@astrophysik.uni-kiel.de>                          #
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
!> Module fosite
!----------------------------------------------------------------------------!
MODULE fosite
  USE fosite_common
  USE sources_generic, CalcTimestep_sources => CalcTimestep
  USE physics_generic, ONLY : InitPhysics, ClosePhysics, MaxWaveSpeeds
  USE fluxes_generic, ONLY : InitFluxes, CloseFluxes, GetBoundaryFlux
  USE reconstruction_generic, ONLY : InitReconstruction
  USE mesh_generic, ONLY : InitMesh, CloseMesh
  USE boundary_generic, ONLY : InitBoundary, WEST, EAST, SOUTH, NORTH, &
                               CenterBoundary
  USE sources_generic, ONLY : InitSources, CloseSources
  USE timedisc_generic, ONLY : InitTimedisc, CloseTimedisc, SolveODE, &
                               CalcTimestep, DTCAUSE_CFL, DTCAUSE_ERRADJ, &
                               DTCAUSE_FILEIO, DTCAUSE_SMALLERR
  USE fileio_generic, ONLY : InitFileIO, CloseFileIO, Initialized, &
                             WriteDataset, AdjustTimestep
  USE timedisc_rkfehlberg, ONLY: ComputeRHS_rkfehlberg
  USE integration
  USE common_dict
#ifdef PARALLEL
  USE common_types, ONLY: DEFAULT_MPI_REAL,DEFAULT_MPI_2REAL 
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
                          SetupFosite, &
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
                          InitFileIO, &
                          GetType, &
                          GetName, &
                          GetRank, &
                          GetNumProcs, &
                          Initialized, &
                          Info, Warning, Error
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
! initialize MPI library for parallel execution, if Fosite is not initialized
    IF(.NOT.Initialized(this)) &
      CALL MPI_Init(this%ierror)
#endif
    ! if Fosite is already initialized, close it, but do not finalize MPI
    IF(Initialized(this)) &
      CALL CloseFosite(this,.FALSE.)

    CALL InitSim(this, simtype, simname)

    IF (GetRank(this).EQ.0) THEN
        ! print some information
        WRITE(buffer, "(A)")&
            "+---------------------------------------------------------+"
        CALL Info(this, buffer)
        WRITE(buffer, "(A1,A29,A28,A1)")&
            "|",TRIM(simname),"","|"
        CALL Info(this, buffer)
        WRITE(buffer, "(A1,A35,A22,A1)")&
            "|",TRIM(VERSION),"","|"
        CALL Info(this, buffer)
        WRITE(buffer, "(A)")&
            "|          Solution of 2D advection problems              |"
        CALL Info(this, buffer)
        WRITE(buffer, "(A)")&
            "+---------------------------------------------------------+"
        CALL Info(this, buffer)
        CALL Info(this, "Initializing simulation:")
    END IF
    CALL InitDict()
    this%iter = 0
  END SUBROUTINE InitFosite

  SUBROUTINE SetupFosite(this)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Fosite_TYP)        :: this
    TYPE(Dict_TYP),POINTER  :: dir, IOdir
    !--------------------------------------------------------------------------!
    INTENT(INOUT)           :: this
    !--------------------------------------------------------------------------!
    IF (.NOT.Initialized(this)) &
      CALL Error(this,"SetupFosite","Sim is uninitialized")

    CALL CopyHierarchy(this%config,this%IO)

    CALL RequireKey(this%config, "mesh")
    CALL GetAttr(this%config, "mesh", dir)
    CALL RequireKey(this%IO, "mesh")
    CALL GetAttr(this%IO, "mesh", IOdir)
    IF(ASSOCIATED(dir)) THEN
       CALL InitMesh(this%Mesh, dir, IOdir)
    END IF

    CALL RequireKey(this%config, "physics")
    CALL GetAttr(this%config, "physics", dir)
    CALL RequireKey(this%IO, "physics")
    CALL GetAttr(this%IO, "physics", IOdir)
    IF(ASSOCIATED(dir)) THEN
        CALL InitPhysics(this%Physics, this%Mesh, dir, IOdir)
    END IF
    
    CALL RequireKey(this%config, "fluxes")
    CALL GetAttr(this%config, "fluxes", dir)
    CALL RequireKey(this%IO, "fluxes")
    CALL GetAttr(this%IO, "fluxes", IOdir)
    IF(ASSOCIATED(dir)) THEN
        CALL InitFluxes(this%Fluxes, this%Mesh, this%Physics, dir, IOdir)
    END IF

    CALL RequireKey(this%config, "boundary")
    CALL GetAttr(this%config, "boundary", dir)
    CALL RequireKey(this%IO, "boundary")
    CALL GetAttr(this%IO, "boundary", IOdir)
    IF(ASSOCIATED(dir)) THEN
        CALL InitBoundary(this%Timedisc%boundary, this%Mesh, this%Physics, dir,IOdir)
    END IF

    CALL GetAttr(this%config, "sources", dir)
    CALL GetAttr(this%IO, "sources", IOdir)
    IF(ASSOCIATED(dir)) THEN
        IF(GetDataType(dir).EQ.DICT_DIR) &
            CALL InitSources(this%Physics%sources, this%Mesh, this%Fluxes, this%Physics,&
                    this%Timedisc%boundary,dir,IOdir)
    END IF

    CALL RequireKey(this%config, "timedisc")
    CALL GetAttr(this%config, "timedisc", dir)
    CALL RequireKey(this%IO, "timedisc")
    CALL GetAttr(this%IO, "timedisc", IOdir)
    IF(ASSOCIATED(dir)) THEN
        CALL InitTimedisc(this%Timedisc,this%Mesh,this%Physics,dir,IOdir)
    END IF

    CALL RequireKey(this%config, "datafile")
    CALL GetAttr(this%config, "datafile", dir)
    IF(ASSOCIATED(dir)) THEN
        CALL InitFileIO(this%Datafile, this%Mesh, this%Physics, this%Timedisc,&
                        this%Physics%sources, dir,this%IO,this%config)
    END IF
    
    CALL GetAttr(this%config, "logfile", dir)
    IF(ASSOCIATED(dir)) THEN
        CALL InitFileIO(this%Logfile, this%Mesh, this%Physics, this%Timedisc,&
                        this%Physics%sources,dir,this%IO)
    END IF

  END SUBROUTINE SetupFosite


  SUBROUTINE FirstStepFosite(this)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Fosite_TYP)   :: this
    INTEGER            :: datetime(8)
    REAL, DIMENSION(this%Mesh%IGMIN:this%Mesh%IGMAX, &
                    this%Mesh%JGMIN:this%Mesh%JGMAX,this%Physics%vnum) &
                       :: trash
    REAL, DIMENSION(this%Mesh%IGMIN:this%Mesh%IGMAX,this%Mesh%JGMIN:this%Mesh%JGMAX,2) &
                       :: amax
    !--------------------------------------------------------------------------!
    INTENT(INOUT)      :: this
    !--------------------------------------------------------------------------!

#ifdef PARALLEL
    this%start_time = MPI_Wtime()
#else
    CALL CPU_TIME(this%start_time)
#endif

    CALL SYSTEM_CLOCK(this%start_count)

    ! make sure that the initial data is written to the log file
    this%wall_time = this%start_time
    this%log_time  = this%wall_time

    IF (GetRank(this).EQ.0) THEN
        CALL Info(this,&
            "===================================================================")
        CALL Info(this, "Starting calculation...")
        CALL date_and_time(values = datetime)
        WRITE(buffer,"(A,I2.2,A,I2.2,A,I4.4,A,I2.2,A,I2.2,A,I2.2)")&
              "Time: ", &
              datetime(3), ".", datetime(2), ".", datetime(1), " ",&
              datetime(5), ":", datetime(6), ":", datetime(7)
       CALL Info(this, buffer)
       CALL Info(this,&
            "step     time        n           t     min(dt) due to    adj.")
       CALL Info(this,&
            "-------------------------------------------------------------------")
    END IF

    ! initialize ghost cell data
    CALL CenterBoundary(this%Timedisc%boundary,this%Mesh,this%Fluxes,&
        this%Physics,this%Timedisc%time,this%Timedisc%pvar,this%Timedisc%cvar)

    CALL MaxWaveSpeeds(this%Physics,this%Mesh,0.0,this%Timedisc%pvar(:,:,:),amax)
    IF(MINVAL(this%Physics%bccsound).LE.0.) &
      CALL Error(this,"FirstStepFosite","Illegal speed of sound value less than 0.")
    IF(this%Physics%csiso.GT.0.) THEN
      IF(ANY(this%Physics%bccsound.NE.this%Physics%csiso)) THEN
        CALL Error(this,"FirstStepFosite","isothermal sound speed set, but "&
        // "arrays bccsound and/or fcsound have been overwritten.")
      END IF
    END IF

    ! store old values
    this%Timedisc%cold(:,:,:) = this%Timedisc%cvar(:,:,:)
    this%Timedisc%pold(:,:,:) = this%Timedisc%pvar(:,:,:)

    ! update all sources for first output
    CALL ExternalSources(this%Physics%Sources,this%Mesh,this%Fluxes, &
               this%Physics,0.0,1.0,this%Timedisc%pvar,this%Timedisc%cvar,trash)

    CALL ComputeRHS_rkfehlberg(this%Timedisc,this%Mesh,this%Physics,this%Fluxes,&
      0.,0.,this%Timedisc%pvar,this%Timedisc%cvar,this%Timedisc%rhs)

    ! store initial data
    IF (this%Timedisc%time.EQ.0.0) THEN
        CALL WriteDataset(this%Datafile,this%Mesh,this%Physics,this%Fluxes,&
                          this%Timedisc,this%config,this%IO)
        IF (GetRank(this).EQ.0) CALL PrintInfo(this,0,0,0.0,0.0,0,this%Timedisc%n_adj)
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
#ifdef PARALLEL
    REAL,DIMENSION(2)  :: dt_buf
#endif
    !--------------------------------------------------------------------------!
    INTENT(INOUT)      :: this
    !--------------------------------------------------------------------------!

    break = .FALSE.

    IF (this%iter.EQ.0) &
        CALL FirstStepFosite(this)

    ! finish simulation if stop time is reached
    IF (ABS(this%Timedisc%stoptime-this%Timedisc%time)&
            .LE.1.0E-05*this%Timedisc%stoptime) THEN
            break = .TRUE.
            RETURN
    END IF

    ! calculate timestep
    this%Timedisc%dt = CalcTimestep(this%Timedisc,this%Mesh,this%Physics,&
                                    this%Fluxes,this%Timedisc%time,this%Timedisc%dtcause)

#ifdef PARALLEL
    ! In Fortran MPI_MINLOC is only able to have two values of the same kind!
    dt_buf(1) = this%Timedisc%dt
    dt_buf(2) = this%Timedisc%dtcause
    CALL MPI_Allreduce(MPI_IN_PLACE,dt_buf,1,DEFAULT_MPI_2REAL,MPI_MINLOC,&
         this%Mesh%comm_cart,this%ierror)
    this%Timedisc%dt = dt_buf(1)
    this%Timedisc%dtcause = dt_buf(2)

    this%run_time = MPI_Wtime() - this%log_time
    CALL MPI_Allreduce(this%run_time,this%wall_time,1,MPI_DOUBLE_PRECISION,&
                       MPI_MIN,this%Mesh%comm_cart,this%ierror)
#else
    CALL CPU_TIME(this%wall_time)
    this%wall_time = this%wall_time - this%log_time
#endif

    ! adjust timestep for output
    ! and calculate the wall clock time
    CALL AdjustTimestep(this%Datafile,this%Timedisc%time,&
                        this%Timedisc%dt,this%Timedisc%dtcause)

#ifdef PARALLEL
    ! In Fortran MPI_MINLOC is only able to have two values of the same kind!
    dt_buf(1) = this%Timedisc%dt
    dt_buf(2) = this%Timedisc%dtcause
    CALL MPI_Allreduce(MPI_IN_PLACE,dt_buf,1,DEFAULT_MPI_2REAL,MPI_MINLOC,&
         this%Mesh%comm_cart,this%ierror)
    this%Timedisc%dt = dt_buf(1)
    this%Timedisc%dtcause = dt_buf(2)
#endif

    ! advance the solution in time
    CALL SolveODE(this%Timedisc,this%Mesh,this%Physics,this%Fluxes,&
                  this%Datafile,this%iter,this%config,this%IO)

    ! write output to data file
    IF ((ABS(this%Datafile%time-this%Timedisc%time)&
            .LE.1.0E-5*this%Datafile%time).OR.&
        this%Timedisc%break) THEN
       CALL WriteDataset(this%Datafile,this%Mesh,this%Physics,this%Fluxes,&
                         this%Timedisc,this%config,this%IO)
       CALL PrintInfo(this, this%Datafile%step-1, this%iter,this%Timedisc%time,&
               this%Timedisc%dtmin,this%Timedisc%dtmincause,this%Timedisc%n_adj)

       ! Stop program if a break is requested from SolveODE
       IF(this%Timedisc%break.AND.GetRank(this).EQ.0) &
         CALL Error(this,"SolveODE", "Time step too small, aborting.",0,.FALSE.)

       ! reset dt_min,dtmincause and n_adj
       this%Timedisc%dtmin = this%Timedisc%stoptime
       this%Timedisc%dtmincause = -99
       this%Timedisc%n_adj = 0
       !IF(GetRank(this).EQ.0) &
       !  WRITE(*,"(A,ES10.4,A,ES10.4)") "dtmean: ", this%Timedisc%dtmean, " +- ",&
       !    SQRT(this%Timedisc%dtstddev/(this%Timedisc%dtaccept-1))/this%Timedisc%dtmean
       this%Timedisc%dtmean = 0.
       this%Timedisc%dtstddev = 0.
       this%Timedisc%dtaccept = 0
       IF(this%Timedisc%write_error) &
         this%Timedisc%error = 0.
    END IF

    ! write output to log file
    IF (Initialized(this%Logfile)) THEN
      IF(this%wall_time.GE.this%Logfile%dtwall) THEN
        CALL WriteDataset(this%Logfile,this%Mesh,this%Physics,this%Fluxes,&
                          this%Timedisc,this%config,this%IO)
        this%log_time = this%wall_time + this%log_time
      END IF
    END IF

  END FUNCTION StepFosite

  SUBROUTINE CloseFosite(this,finalize_)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Fosite_TYP)   :: this
    LOGICAL,OPTIONAL   :: finalize_
    !--------------------------------------------------------------------------!
    LOGICAL            :: finalize = .TRUE.
    !--------------------------------------------------------------------------!
    INTENT(INOUT)      :: this
    !--------------------------------------------------------------------------!

#ifdef PARALLEL
    this%end_time = MPI_Wtime() - this%start_time
    CALL MPI_Allreduce(this%end_time,this%run_time,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
        this%Mesh%comm_cart,this%ierror)
#else
    CALL CPU_TIME(this%end_time)
    this%run_time = this%end_time - this%start_time
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
    CALL DeleteDict(this%config)
    CALL DeleteDict(this%IO)

#ifdef PARALLEL
    IF(PRESENT(finalize_)) &
      finalize = finalize_
    IF(finalize) &
      CALL MPI_Finalize(this%ierror)
#endif

    CALL CloseSim(this)
  END SUBROUTINE CloseFosite


  SUBROUTINE PrintInfo(this,step,i,t,d,dc,na)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)   :: this
    INTEGER            :: step, i, dc, na
    CHARACTER(LEN=9)   :: dtcause
    REAL               :: t, d
    !------------------------------------------------------------------------!
    INTEGER            :: drt, c, c_rate, c_max
    DOUBLE PRECISION   :: time 
    !------------------------------------------------------------------------!
    INTENT(IN)         :: i,t,d
    INTENT(INOUT)      :: this
    !--------------------------------------------------------------------------!
    IF (GetRank(this).EQ.0) THEN

       SELECT CASE (dc)
        ! positive values represent source terms
        CASE(1:)
          WRITE(dtcause, "(A,I2.2,A)") " S",dc," "
        CASE(DTCAUSE_CFL)
          WRITE(dtcause, "(A)") " cfl "
        CASE(DTCAUSE_ERRADJ)
           WRITE(dtcause, "(A)") " err_adj "
        CASE(DTCAUSE_SMALLERR)
           WRITE(dtcause, "(A)") " err "
        CASE(DTCAUSE_FILEIO)
           ! output by this reason is suppressed by default 
           WRITE(dtcause, "(A)") " fileio "
        CASE DEFAULT
           WRITE(dtcause, "(A,I3.2,A)") " ?", dc, "? "
       END SELECT

       CALL SYSTEM_CLOCK(c,c_rate, c_max)
       ! overflow every 24days => assumption: max one per simulation 
       drt = (c - this%start_count)/c_rate
       IF (drt .LT. 0) THEN
         drt = c_max/c_rate + drt
       END IF
       WRITE(buffer,"(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I8,A,ES11.3,A,ES11.3,A,I5)")&
             step, " ", drt/3600, ":", mod(drt,3600)/60, ":", mod(drt,60),&
             " ", i, " ", t, " ", d, dtcause, na
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
       WRITE(buffer,"(A,F10.2,A)")" main loop runtime: ", this%run_time, " sec."
       CALL Info(this, buffer)
    END IF
  END SUBROUTINE PrintSummary

END MODULE fosite
