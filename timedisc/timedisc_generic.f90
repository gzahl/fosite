!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: timedisc_generic.f90                                              #
!#                                                                           #
!# Copyright (C) 2007-2012                                                   #
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
!> \addtogroup timedisc
!! - general parameters of timedisc group as key-values
!! \key{method,INTEGER,time integration method}
!! \key{stoptime,REAL,physical stop time of simulation}
!! \key{cfl,REAL,CFL number,0.4}
!! \key{dtlimit,REAL,time step minimum,EPSILON(dtlimit)*stoptime}
!! \key{dtmax,REAL,time step maximum in units of [CFL timestep] 
!!   (Used in Dumka), 5}
!! \key{maxiter,INTEGER,maximum iterations,HUGE(maxiter)}
!! \key{tol_rel,REAL,relative tolerance for adaptive step size control,0.01}
!! \key{tol_abs,REAL\, DIMENSION(Physics%VNUM), absolute tolerance for adaptive 
!!   step size control, (/0.001\,0.001\,../)}
!! \key{fargo,INTEGER,enable(=1) fargo timestepping for polar geometries,0}
!! \key{beta,REAL,time step friction parameter for PI-Controller,0.08}

!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Björn Sperling
!! \author Manuel Jung
!!
!! \brief generic subroutines for time discretization
!!
!! \ingroup timedisc
!----------------------------------------------------------------------------!
MODULE timedisc_generic
  USE timedisc_modeuler, CloseTimedisc_common => CloseTimedisc
  USE timedisc_rkfehlberg
  USE timedisc_cashkarp
  USE timedisc_ssprk
  USE timedisc_dumka
  USE timedisc_dormand_prince
  USE timedisc_multistep
  USE boundary_generic
  USE mesh_generic
  USE physics_generic
  USE fluxes_generic
  USE reconstruction_generic, ONLY: CalculateSlopes
  USE sources_generic, CalcTimestep_sources => CalcTimestep
  USE fileio_generic
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
  INTEGER, PARAMETER :: MODIFIED_EULER = 1
  INTEGER, PARAMETER :: RK_FEHLBERG    = 2
  INTEGER, PARAMETER :: CASH_KARP      = 3
  INTEGER, PARAMETER :: DUMKA          = 4
  INTEGER, PARAMETER :: DORMAND_PRINCE = 5
  INTEGER, PARAMETER :: MULTISTEP      = 6
  INTEGER, PARAMETER :: SSPRK          = 7
  !--------------------------------------------------------------------------!
  INTEGER, PARAMETER :: DTCAUSE_CFL    =  0  ! smallest ts due to cfl cond.  !
  INTEGER, PARAMETER :: DTCAUSE_ERRADJ = -1  ! smallest ts due to err adj.   !
  PUBLIC :: &
       ! types
       Timedisc_TYP, &
       ! constants
       MODIFIED_EULER, RK_FEHLBERG, CASH_KARP, DUMKA, DORMAND_PRINCE, MULTISTEP, &
       SSPRK, &
       DTCAUSE_CFL,DTCAUSE_ERRADJ,DTCAUSE_FILEIO,DTCAUSE_SMALLERR, &
       ! methods 
       InitTimedisc, &
       CalcTimestep, &
       SolveODE, &
       CloseTimedisc, &
       GetOrder, &
       GetCFL, &
       GetCentrifugalVelocity, &
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

  SUBROUTINE InitTimedisc(this,Mesh,Physics,config,IO,dumpfile)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fileio_TYP), OPTIONAL :: Dumpfile
    TYPE(Dict_TYP), POINTER &
                       :: config,IO
    !------------------------------------------------------------------------!
    INTEGER            :: err
    CHARACTER(LEN=32)   :: order_str,cfl_str,stoptime_str,dtmax_str,beta_str
    CHARACTER(LEN=32)  :: info_str 
    REAL               :: tol_abs_def(Physics%VNUM)
    INTEGER            :: method
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,Physics,Dumpfile
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    IF (.NOT.Initialized(Physics).OR..NOT.Initialized(Mesh)) &
         CALL Error(this,"InitTimedisc","physics and/or mesh module uninitialized")

    CALL RequireKey(config, "method")
    CALL RequireKey(config, "stoptime")
 
    CALL GetAttr(config, "method", method)
    CALL GetAttr(config, "stoptime", this%stoptime)
 
    ! set default values
    ! CFL number
    CALL RequireKey(config, "cfl", 0.4)
    
    ! time step minimum
    CALL RequireKey(config, "dtlimit", EPSILON(this%dtlimit)*this%stoptime)
    
    ! time step maximum in units of [CFL timestep] (Used in Dumka)
    CALL RequireKey(config, "dtmax", 5.0)    

    ! maximum iterations
    CALL RequireKey(config, "maxiter", HUGE(this%maxiter))
    
    ! relative tolerance for adaptive step size control
    CALL RequireKey(config, "tol_rel", 0.01)    ! 1%
    
    ! absolute tolerance for adaptive step size control
    tol_abs_def(:) = 0.001
    CALL RequireKey(config, "tol_abs", tol_abs_def(:))

    ! enable fargo timestepping for polar geometries
    CALL RequireKey(config, "fargo", 0) ! 0 = disable, 1 = enable

    ! time step friction parameter for PI-Controller
    CALL RequireKey(config, "beta", 0.08)

    CALL RequireKey(config, "output/"//TRIM(Physics%pvarname(1)), 1)

    ! allocate memory for data structures needed in all timedisc modules
    ALLOCATE(this%pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
         this%cvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
         this%pold(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
         this%cold(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
         this%ptmp(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
         this%ctmp(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
         this%geo_src(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
         this%src(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
         this%rhs(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
         this%rhsold(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
         this%xfluxdy(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
         this%yfluxdx(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
         this%amax(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
         this%tol_abs(Physics%VNUM), &
         this%shift(Mesh%IGMIN:Mesh%IGMAX), &
         this%dtmean, this%dtstddev, &
         STAT = err)
    IF (err.NE.0) THEN
       CALL Error(this,"InitTimedisc", "Unable to allocate memory.")
    END IF

    CALL GetAttr(config, "cfl", this%cfl)
    CALL GetAttr(config, "dtlimit", this%dtlimit)
    CALL GetAttr(config, "dtmax", this%dtmax)
    CALL GetAttr(config, "maxiter", this%maxiter)
    CALL GetAttr(config, "fargo", this%fargo)
    CALL GetAttr(config, "tol_rel", this%tol_rel)
    CALL GetAttr(config, "tol_abs", this%tol_abs)
    CALL GetAttr(config, "fargo", this%fargo)
    CALL GetAttr(config, "beta", this%beta)

    ! call individual constructors
    SELECT CASE(method)
    CASE(MODIFIED_EULER)
       CALL InitTimedisc_modeuler(this,config)
    CASE(RK_FEHLBERG)
       CALL InitTimedisc_rkfehlberg(this,Mesh,Physics,config)
    CASE(CASH_KARP)
       CALL InitTimedisc_cashkarp(this,Mesh,Physics,config)
    CASE(DUMKA)
       CALL InitTimedisc_dumka(this,Mesh,Physics,config)
    CASE(DORMAND_PRINCE)
       CALL InitTimedisc_dormand_prince(this,Mesh,Physics,config)
    CASE(MULTISTEP)
       CALL InitTimedisc_multistep(this,Mesh,Physics,config)
    CASE(SSPRK)
       CALL InitTimedisc_ssprk(this,Mesh,Physics,config)
    CASE DEFAULT
       CALL Error(this,"InitTimedisc", "Unknown ODE solver.")
    END SELECT

    ! check dumpfile
    IF (PRESENT(dumpfile)) THEN
       IF(.NOT.Initialized(Dumpfile)) &
            CALL Warning(this,"InitTimedisc", &
            "dump file uninitialized, dumping disabled")
    END IF

    ! initialize all variables
    this%pvar = 0.
    this%cvar = 0.
    this%pold = 0.
    this%cold = 0.
    this%ctmp = 0.
    this%src = 0.
    this%geo_src = 0.
    this%rhs = 0.
    this%xfluxdy = 0.
    this%yfluxdx = 0.
    this%amax = 0.
    this%shift = 0.
    this%break = .FALSE.

    CALL SetOutput(this,Mesh,Physics,config,IO)    

    ! check if fargo can be used, if requested
    IF(this%fargo.EQ.1) THEN
      SELECT CASE(GetType(Physics))
      CASE(EULER2D_ISOIAMT,EULER2D_IAMT)
        ! do nothing
      CASE DEFAULT
        this%fargo = 0
        CALL Warning(this,"InitTimedisc","fargo has been disabled, because the physics are not supported.")
      END SELECT
      SELECT CASE(GetType(Mesh%geometry))
      CASE(POLAR,TANPOLAR,LOGPOLAR,SINHPOLAR)
        ! do nothing
      CASE DEFAULT
        this%fargo = 0
        CALL Warning(this,"InitTimedisc","fargo has been disabled, because the physics are not supported.")
      END SELECT
#ifdef PARALLEL
      IF(this%fargo.EQ.1) THEN
        ALLOCATE(this%buf(Physics%VNUM,1:Mesh%MINJNUM), &
                 STAT = err)
        IF (err.NE.0) THEN
          CALL Error(this,"InitTimedisc", "Unable to allocate memory.")
        END IF
      END IF
#endif
    END IF

    ! print some information
    WRITE (order_str, '(I0)') GetOrder(this)
    WRITE (cfl_str, '(F4.2)') GetCFL(this)
    WRITE (stoptime_str, '(ES10.4)') this%stoptime
    WRITE (dtmax_str, '(ES10.4)') this%dtmax
    WRITE (beta_str, '(ES10.4)') this%beta
    CALL Info(this," TIMEDISC-> ODE solver:        " //TRIM(GetName(this)))
    CALL Info(this,"            order:             " //TRIM(order_str))
    CALL Info(this,"            CFL number:        " //TRIM(cfl_str))
    CALL Info(this,"            dtmax:             " //TRIM(dtmax_str))
    CALL Info(this,"            stoptime:          " //TRIM(stoptime_str))
    CALL Info(this,"            beta:              " //TRIM(beta_str))
    IF(this%fargo.EQ.1) &
      CALL Info(this,"            fargo:             " // "enabled")
    ! adaptive step size control
    IF (this%tol_rel.LT.1.0) THEN
       WRITE (info_str,'(ES7.1)') this%tol_rel*100
       CALL Info(this,"            step size control: enabled")
       CALL Info(this,"            rel. precision:    "//TRIM(info_str)//" %")
    ELSE
       WRITE (info_str,'(A)') "disabled"       
    END IF

  END SUBROUTINE InitTimedisc

 SUBROUTINE SetOutput(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !-----------------------------------------------------------------------!
    TYPE(Timedisc_TYP)   :: this
    TYPE(Mesh_TYP)       :: Mesh
    TYPE(Physics_TYP)    :: Physics
    TYPE(Dict_TYP),POINTER  :: config,IO
    !------------------------------------------------------------------------!
    INTEGER              :: valwrite,i,err
    CHARACTER(LEN=60)    :: key
    !------------------------------------------------------------------------!
    INTENT(IN)           :: Mesh,Physics
    INTENT(INOUT)        :: this
    !------------------------------------------------------------------------! 
    valwrite = 0
    IF(HasKey(config, "output/error")) &
      CALL GetAttr(config, "output/error", valwrite)
    IF(valwrite.EQ.1) THEN
      this%write_error = .TRUE.
      ALLOCATE(this%error(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
        STAT = err)
      IF (err.NE.0) &
        CALL Error(this,"SetOutput_timedisc", "Unable to allocate memory.")
      this%error = 0.
    ELSE
      this%write_error = .FALSE.
    END IF

!    CALL AddField(IO, "dtmean", this%dtmean, Dict("name" / "dtmean"))
!    CALL AddField(IO, "dtstddev", this%dtstddev, Dict("name" / "dtstddev"))

    DO i=1, Physics%VNUM
      !prim
      key = TRIM(Physics%pvarname(i))
      valwrite = 0
      CALL RequireKey(config, "output/" // TRIM(key), 1)
      CALL GetAttr(config, "output/" // TRIM(key), valwrite)
      ! second argument is important if pvarname is used twice
      IF (valwrite.EQ.1 .AND. (.NOT.HasKey(IO, TRIM(key)))) THEN
        CALL AddField(IO, &
                      TRIM(key), &
                      remap_bounds(Mesh,this%pvar(:,:,i)), &
                      Dict("name" / TRIM(key)))
      END IF

      !cons
      key = TRIM(Physics%cvarname(i))
      valwrite = 0
      IF (HasKey(config, "output/" // TRIM(key))) &
        CALL GetAttr(config, "output/" // TRIM(key), valwrite)
      ! second argument is important if pvarname is used twice
      IF (valwrite.EQ.1 .AND. (.NOT.HasKey(IO, TRIM(key)))) THEN
           CALL AddField(IO, &
                      TRIM(key), &
                      remap_bounds(Mesh,this%cvar(:,:,i)), &
                      Dict("name" / TRIM(key)))
      END IF
      
      key = TRIM(Physics%cvarname(i))
      IF(this%write_error) THEN
        CALL AddField(IO, &
                    "error_" // TRIM(key), &
                    remap_bounds(Mesh,this%error(:,:,i)), &
                    Dict("name" / ("error_" // TRIM(key))))
      END IF

      ! write geometrical sources 
      valwrite = 0
      IF (HasKey(config, "output/" // "geometrical_sources")) &
        CALL GetAttr(config, "output/" // "geometrical_sources", valwrite)
      IF (valwrite.EQ.1) THEN
           CALL AddField(IO, &
                      TRIM(key)//"_geo_src", &
                      remap_bounds(Mesh,this%geo_src(:,:,i)), &
                      Dict("name" / (TRIM(key)//"_geo_src")))
      END IF

      ! write external sources 
      valwrite = 0
      IF (HasKey(config, "output/" // "external_sources")) &
        CALL GetAttr(config, "output/" // "external_sources", valwrite)
      IF (valwrite.EQ.1) THEN
           CALL AddField(IO, &
                      TRIM(key)//"_src", &
                      remap_bounds(Mesh,this%src(:,:,i)), &
                      Dict("name" / (TRIM(key)//"_src")))
      END IF

      ! write right hand side
      valwrite = 0
      IF (HasKey(config, "output/" // "rhs")) &
        CALL GetAttr(config, "output/" // "rhs", valwrite)
      IF (valwrite.EQ.1) THEN
           CALL AddField(IO, &
                      TRIM(key)//"_rhs", &
                      remap_bounds(Mesh,this%rhs(:,:,i)), &
                      Dict("name" / (TRIM(key)//"_rhs")))
      END IF

      ! write fluxes 
      ! ATTENTION: this are the numerical fluxes devided by dy or dx respectively
      valwrite = 0
      IF (HasKey(config, "output/" // "fluxes")) &
        CALL GetAttr(config, "output/" // "fluxes", valwrite)
      IF (valwrite.EQ.1) THEN
           CALL AddField(IO, &
                      TRIM(key)//"_xfluxdy", &
                      remap_bounds(Mesh,this%xfluxdy(:,:,i)), &
                      Dict("name" / (TRIM(key)//"_xfluxdy")))
           CALL AddField(IO, &
                      TRIM(key)//"_yfluxdx", &
                      remap_bounds(Mesh,this%yfluxdx(:,:,i)), &
                      Dict("name" / (TRIM(key)//"_yfluxdx")))
      END IF

    END DO
    ! write bfluxes
    valwrite = 0
    IF(HasKey(config, "output/" // "bflux")) &
        CALL GetAttr(config, "output/" // "bflux", valwrite)
    IF(valwrite.EQ.1) THEN
        ALLOCATE(this%bflux(Physics%VNUM,4))
        CALL AddField(IO, &
                      "bflux", &
                      remap_bounds(this%bflux), &
                      Dict("name" / "bflux"))
    ELSE
        NULLIFY(this%bflux)
    END IF
  END SUBROUTINE SetOutput


  REAL FUNCTION CalcTimestep(this,Mesh,Physics,Fluxes,time,dtcause) RESULT(dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    REAL               :: time
    INTEGER            :: dtcause
    !------------------------------------------------------------------------!    
    INTEGER            :: i, err
    REAL               :: wi
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh, time
    INTENT(INOUT)      :: this,Physics,dtcause
    !------------------------------------------------------------------------!   
    IF(this%fargo.EQ.1) THEN
      DO i=Mesh%IGMIN,Mesh%IGMAX
        wi = SUM(this%pvar(i,Mesh%JMIN:Mesh%JMAX,Physics%YVELOCITY))
#ifdef PARALLEL
        IF(Mesh%dims(2).GT.1) THEN
          CALL MPI_AllReduce(MPI_IN_PLACE, wi, 1, DEFAULT_MPI_REAL, MPI_SUM, &
                             Mesh%Icomm, err)
        END IF
#endif
        Physics%w(i,Mesh%JGMIN:Mesh%JGMAX,:) &
          = wi / (Mesh%bhy(i,Mesh%JMIN) * Mesh%JNUM)
      END DO
    END IF

!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(MODIFIED_EULER)
      dt = CalcTimestep_modeuler(this,Mesh,Physics,Fluxes,time,dtcause)
    CASE(RK_FEHLBERG)
      dt = CalcTimestep_rkfehlberg(this,Mesh,Physics,Fluxes,time,dtcause)
    CASE(CASH_KARP)
      dt = CalcTimestep_cashkarp(this,Mesh,Physics,Fluxes,time,dtcause)
    CASE(DUMKA)
      dt = CalcTimestep_dumka(this,Mesh,Physics,Fluxes,time,dtcause)
    CASE(DORMAND_PRINCE)
      dt = CalcTimestep_dormand_prince(this,Mesh,Physics,Fluxes,time,dtcause)
    CASE(MULTISTEP)
      dt = CalcTimestep_multistep(this,Mesh,Physics,Fluxes,time,dtcause)
    CASE(SSPRK)
      dt = CalcTimestep_ssprk(this,Mesh,Physics,Fluxes,time,dtcause)
    END SELECT
  END FUNCTION CalcTimestep


  SUBROUTINE SolveODE(this,Mesh,Physics,Fluxes,Datafile,iter,config,IO)
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    TYPE(FileIO_TYP)   :: Datafile
    INTEGER            :: iter
    TYPE(Dict_TYP),POINTER :: config,IO
    !------------------------------------------------------------------------!
    INTEGER            :: i,j
#ifdef PARALLEL
    REAL               :: err_all,dt_all
    INTEGER            :: ierror
#endif
    REAL               :: err,dtold,dt,time,dtmeanold
    INTEGER            :: M,k
    REAL               :: f
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh
    INTENT(INOUT)      :: this,Physics,Fluxes,iter
    !------------------------------------------------------------------------!
    time = this%time
    dt   = this%dt
    IF (dt.LT.this%dtmin .AND. this%dtcause .NE. DTCAUSE_FILEIO) THEN
      ! only save dtmin if the reasion is not the fileio
      this%dtmin = dt
      this%dtmincause = this%dtcause
    END IF
    timestep: DO WHILE (time+dt.LE.this%time+this%dt)
      dtold = dt
!CDIR IEXPAND
      SELECT CASE(GetType(this))
      CASE(MODIFIED_EULER)
        CALL SolveODE_modeuler(this,Mesh,Physics,Fluxes,time,dt,err)
      CASE(RK_FEHLBERG)
        CALL SolveODE_rkfehlberg(this,Mesh,Physics,Fluxes,time,dt,err)
      CASE(CASH_KARP)
        CALL SolveODE_cashkarp(this,Mesh,Physics,Fluxes,time,dt,err)
      CASE(DUMKA)
        CALL SolveODE_dumka(this,Mesh,Physics,Fluxes,time,dt,err)
      CASE(DORMAND_PRINCE)
        CALL SolveODE_dormand_prince(this,Mesh,Physics,Fluxes,time,dt,err)
      CASE(MULTISTEP)
        CALL SolveODE_multistep(this,Mesh,Physics,Fluxes,time,dt,err)
      CASE(SSPRK)
        CALL SolveODE_ssprk(this,Mesh,Physics,Fluxes,time,dt,err)
      END SELECT
      ! check truncation error and restart if necessary
      IF (err.LT.1.0) THEN
        time=time+dtold
        this%dtaccept = this%dtaccept + 1
        dtmeanold = this%dtmean
        this%dtmean = this%dtmean + (dtold - this%dtmean)/this%dtaccept
        this%dtstddev = this%dtstddev + (dtold - dtmeanold)*(dtold-this%dtmean)
        this%cold(:,:,:) = this%cvar(:,:,:)
        this%pold(:,:,:) = this%pvar(:,:,:)
        Fluxes%bxfold(:,:,:) = Fluxes%bxflux(:,:,:)
        Fluxes%byfold(:,:,:) = Fluxes%byflux(:,:,:)
        iter = iter + 1
!!$          PRINT '(A,4(A,ES12.6))'," Horray!"," t=",time," err=",err,&
!!$               " dtold=",dtold," dt=",dt
      ELSE
        this%cvar(:,:,:) = this%cold(:,:,:)
        this%pvar(:,:,:) = this%pold(:,:,:)
        Fluxes%bxflux(:,:,:) = Fluxes%bxfold(:,:,:)
        Fluxes%byflux(:,:,:) = Fluxes%byfold(:,:,:)
        ! count adjustments for information
        this%n_adj = this%n_adj + 1
        this%dtcause = DTCAUSE_ERRADJ
        ! only save dtmin if the reasion is not the fileio (automatically satisfied)
        IF (dt.LT.this%dtmin) THEN
          this%dtmin = dt
          this%dtmincause = this%dtcause
        END IF
!!$         PRINT '(A,4(A,ES12.6))'," Argggh!"," t=",time," err=",err,&
!!$               " dtold=",dtold," dt=",dt
      END IF
      IF (dt.LT.this%dtlimit) THEN
        this%break = .TRUE.
        ! Do not attempt to fargo shift anymore
        this%fargo = 0
        EXIT timestep
      END IF
!!$       ! check data
!!$       bad_data = CheckData(Physics,Mesh,this%pvar,this%pold)
!!$       IF (bad_data.NE.0) THEN
!!$          IF ((this%dt * 0.5).LT.this%dtlimit) THEN
!!$             CALL Print_Checkdata(this,Mesh,Physics)
!!$             CALL Error(this,"SolveODE", "Time step to small, aborting.",GetRank(this))
!!$          END IF
!!$       END IF
!!$#ifdef PARALLEL
!!$       CALL MPI_Allreduce(bad_data,bad_data_all,1,MPI_LOGICAL,MPI_LOR,Mesh%comm_cart,ierror)
!!$       bad_data = bad_data_all
!!$#endif
!!$
    END DO timestep

    this%dt = time - this%time
    this%time  = time
    this%dtold = dt

    IF(this%fargo.EQ.1) &
      CALL FargoAdvection(this,Fluxes,Mesh,Physics)

  END SUBROUTINE SolveODE


  SUBROUTINE FargoAdvection(this,Fluxes,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP)   :: this
    TYPE(Fluxes_TYP)     :: Fluxes
    TYPE(Mesh_TYP)       :: Mesh
    TYPE(Physics_TYP)    :: Physics
    !------------------------------------------------------------------------!
    REAL                 :: dely,f
    INTEGER              :: i,j,k,ierror
#ifdef PARALLEL
    INTEGER              :: status(MPI_STATUS_SIZE)
#endif
    !------------------------------------------------------------------------!
    f = this%dt * Mesh%invdy
    DO i=Mesh%IGMIN,Mesh%IGMAX
#ifdef PARALLEL
      IF(Mesh%dims(2).GT.1) THEN
        this%shift(i) = MIN(NINT(Physics%w(i,Mesh%JMIN,1) * f), Mesh%MINJNUM)
      ELSE
        this%shift(i) = NINT(Physics%w(i,Mesh%JMIN,1) * f)
      END IF
#else
      this%shift(i) = NINT(Physics%w(i,Mesh%JMIN,1) * f)
#endif
    END DO

    ! advect with residual velocity
    CALL CalculateSlopes(Fluxes%Reconstruction,Mesh,Physics,this%cvar)
    DO i=Mesh%IMIN,Mesh%IMAX
      dely = (Physics%w(i,Mesh%JMIN,1)*this%dt-this%shift(i)*Mesh%dy) * Mesh%invdy
      DO k=1,Physics%VNUM
        DO j=Mesh%JMIN,Mesh%JMAX
          IF(dely.GE.0.) THEN
            this%cvar(i,j,k) = this%cold(i,j,k) &
              - dely & 
                * (this%cold(i,j,k)-this%cold(i,j-1,k) &
                   + 0.5 * Fluxes%reconstruction%yslopes(i,j,k) * Mesh%dy &
                   *( 1. - dely) &
                   - 0.5 * Fluxes%reconstruction%yslopes(i,j-1,k) * Mesh%dy&
                   *( 1. - dely))
          ELSE
            this%cvar(i,j,k) = this%cold(i,j,k) &
              - dely & 
                * (this%cold(i,j+1,k)-this%cold(i,j,k) &
                   - 0.5 * Fluxes%reconstruction%yslopes(i,j+1,k) * Mesh%dy&
                   *( 1. + dely) &
                   + 0.5 * Fluxes%reconstruction%yslopes(i,j,k) * Mesh%dy&
                   *( 1. + dely))
          END IF
        END DO
      END DO
    END DO
    Physics%w(:,:,:) = 0.


#ifdef PARALLEL
    ! We only need to do something, if we (also) are dealing with domain decomposition in
    ! the second (phi) direction
    IF(Mesh%dims(2).GT.1) THEN
      DO i=Mesh%IGMIN,Mesh%IGMAX
        IF(this%shift(i).GT.0) THEN
          DO k=1,Physics%VNUM
            this%buf(k,1:this%shift(i)) = this%cvar(i,Mesh%JMAX-this%shift(i)+1:Mesh%JMAX,k)
          END DO
          CALL MPI_Sendrecv_replace(&
            this%buf,&
            this%shift(i)*Physics%VNUM, &
            DEFAULT_MPI_REAL, &
            Mesh%neighbor(NORTH), 130, &
            Mesh%neighbor(SOUTH), 130, &
            Mesh%comm_cart, status, ierror)
          DO k=1,Physics%VNUM
            this%cvar(i,Mesh%JMAX-this%shift(i)+1:Mesh%JMAX,k) = this%buf(k,1:this%shift(i))
          END DO
        ELSE IF(this%shift(i).LT.0) THEN
          DO k=1,Physics%VNUM
            this%buf(k,1:-this%shift(i)) = this%cvar(i,Mesh%JMIN:Mesh%JMIN-this%shift(i)-1,k)
          END DO
          CALL MPI_Sendrecv_replace(&
            this%buf,&
            -this%shift(i)*Physics%VNUM, &
            DEFAULT_MPI_REAL, &
            Mesh%neighbor(SOUTH), 150, &
            Mesh%neighbor(NORTH), 150, &
            Mesh%comm_cart, status, ierror)
          DO k=1,Physics%VNUM
            this%cvar(i,Mesh%JMIN:Mesh%JMIN-this%shift(i)-1,k) = this%buf(k,1:-this%shift(i))
          END DO
        END IF
      END DO
    END IF
#endif
    DO k=1,Physics%VNUM
      this%cvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMIN:Mesh%JMAX,k) &
        = CSHIFT(this%cvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMIN:Mesh%JMAX,k),-this%shift,2)
    END DO
    CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,this%time,this%pvar,this%cvar)

    this%cold(:,:,:) = this%cvar(:,:,:)
    this%pold(:,:,:) = this%pvar(:,:,:)
  END SUBROUTINE FargoAdvection


  FUNCTION GetCentrifugalVelocity(this,Mesh,Physics,Fluxes,&
                       dir_omega_,accel_,centrot) RESULT(velo)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP):: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM) &
                      :: velo
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2),OPTIONAL &
                      :: accel_
    REAL, DIMENSION(2), OPTIONAL :: centrot
    REAL, DIMENSION(3):: dir_omega_
    !------------------------------------------------------------------------!
    REAL, DIMENSION(3):: dir_omega
    REAL              :: omega2
    INTEGER           :: k,i,j
    REAL              :: rotoemga
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
                      :: bccart, bcposvec,  accel
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: tmp
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,centrot,dir_omega_,accel_
    INTENT(INOUT)     :: this,Physics,Fluxes
    !------------------------------------------------------------------------!

    IF(PRESENT(accel_)) THEN
      accel = accel_
    ELSE
      ! Works only for centrot = 0
      IF(PRESENT(centrot)) &
        CALL Error(this,"GetCentrifugalVelocity","You are not allowed to "&
          //"define centrot without accel.")
      ! This may not work for physics with unusual conservative variables.
      ! We assume
      !   conservative momentum = density * velocity
      ! This is not true for the second component of physics with
      ! angular momentum transport, but that component should be zero.
      SELECT CASE(GetType(Physics))
      CASE(EULER2D,EULER2D_ISOTHERM,EULER3D_ROTSYM,EULER3D_ROTAMT,&
           EULER3D_ROTSYMSGS,EULER2D_SGS,EULER3D_ROTAMTSGS,EULER2D_ISOIAMT, &
           EULER2D_IAMT, EULER2D_IAMROT,EULER2D_ISOIAMROT)
        ! do nothing
      CASE DEFAULT
        CALL Error(this,"GetCentrifugalVelocity","It is unknown, if the "&
          //"selected physics module works with this routine.")
      END SELECT
      CALL ComputeRHS_rkfehlberg(this,Mesh,Physics,Fluxes,&
        this%time,0.,this%pvar,this%cvar,this%rhs)
      DO k=Physics%XMOMENTUM,Physics%YMOMENTUM
        accel(:,:,k-Physics%XMOMENTUM+1) = -1. * this%rhs(:,:,k) &
                                           / this%pvar(:,:,Physics%DENSITY)
      END DO
    END IF

    dir_omega = dir_omega_

    ! be sure: |dir_omega| == 1
    omega2 = SUM(dir_omega(:)*dir_omega(:))
    ! omega must not be the zero vector
    IF (omega2 .EQ. 0.0) &
        CALL Error(this,"GetCentrifugalVelocity", &
           "omega must not be the zero vector")
    ! norm must be one
    IF (omega2 .NE. 1.0) dir_omega(:) = dir_omega(:) / SQRT(omega2)

    IF ((Physics%DIM .EQ. 2) .AND. &
        ((dir_omega(1) .NE. 0.0) .OR. (dir_omega(2) .NE. 0.0))) &
        CALL Error(this,"GetCentrifugalVelocity", &
           "the direction of omega should be (0,0,+-1) in case of two dimensions")

!attention: centrot is in cartesian coor
    IF (present(centrot)) THEN
      ! translate the position vector to the center of rotation
      bccart(:,:,1) = Mesh%bccart(:,:,1) - centrot(1)
      bccart(:,:,2) = Mesh%bccart(:,:,2) - centrot(2)

      ! compute curvilinear components of translated position vectors
      CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,bccart,bcposvec)
    ELSE
      bcposvec = Mesh%bposvec
    END IF


    ! compute distance to axis of rotation (It is automatically fulfilled in 2D.)
    IF (Physics%DIM .GT. 2) THEN
      tmp(:,:) =  bcposvec(:,:,1)*dir_omega(1)&
                    + bcposvec(:,:,2)*dir_omega(2)
      bcposvec(:,:,1) = bcposvec(:,:,1) - tmp(:,:)*dir_omega(1)
      bcposvec(:,:,2) = bcposvec(:,:,2) - tmp(:,:)*dir_omega(2)
    END IF

    ! compute omega = SQRT(-dot(g,r)/|r|**2)
    tmp(:,:) = SQRT(MAX(0.0,-SUM(accel(:,:,1:2)*bcposvec(:,:,:),DIM=3))&
                    / SUM(bcposvec(:,:,:)*bcposvec(:,:,:),DIM=3))

    ! v / |omega| = dir_omega x r
    velo(:,:,:) = CROSS_PRODUCT(Mesh,Physics,dir_omega,bcposvec)
    ! v = |omega| * dir_omega x r
    DO k=1,Physics%DIM
      velo(:,:,k) = tmp(:,:)*velo(:,:,k)
    END DO

  CONTAINS
    FUNCTION CROSS_PRODUCT(Mesh,Physics,a,b) RESULT(cp)
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      TYPE(Mesh_TYP)    :: Mesh
      TYPE(Physics_TYP) :: Physics
      REAL, DIMENSION(3):: a
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
                        :: b
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM) &
                        :: cp
      !------------------------------------------------------------------------!
      INTENT(IN)        :: Mesh,Physics,a,b
      !------------------------------------------------------------------------!
      IF (Physics%DIM .EQ. 3) THEN
        ! cp(:,:,1) = a(2)*b(:,:,3) - a(3)*b(:,:,2)
        ! cp(:,:,2) = a(3)*b(:,:,1) - a(1)*b(:,:,3)
        ! due to symmetry => axis of rotation only along 1,2 => a(3) == 0 and b(,3) == 0
        cp(:,:,1:2) = 0.
        cp(:,:,3) = a(1)*b(:,:,2) - a(2)*b(:,:,1)
      ELSE
        ! => a(1) = a(2) = 0 and a(3) = 1 or -1
        cp(:,:,1) = - a(3)*b(:,:,2)
        cp(:,:,2) = a(3)*b(:,:,1)
      END IF
     END FUNCTION CROSS_PRODUCT
  END FUNCTION GetCentrifugalVelocity



  SUBROUTINE CloseTimedisc(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP)   :: this
    !------------------------------------------------------------------------!
    IF (.NOT.Initialized(this)) &
        CALL Error(this,"CloseTimedisc","not initialized")
    ! call boundary destructor
    CALL CloseBoundary(this%Boundary)

    ! call individual destructors
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(MODIFIED_EULER)
       CALL CloseTimedisc_modeuler(this)
    CASE(RK_FEHLBERG)
       CALL CloseTimedisc_rkfehlberg(this)
    CASE(CASH_KARP)
       CALL CloseTimedisc_cashkarp(this)
    CASE(DUMKA)
       CALL CloseTimedisc_dumka(this)
    CASE(DORMAND_PRINCE)
       CALL CloseTimedisc_dormand_prince(this)
    CASE(MULTISTEP)
       CALL CloseTimedisc_multistep(this)
    CASE(SSPRK)
       CALL CloseTimedisc_ssprk(this)
    END SELECT

    DEALLOCATE(this%pvar,this%cvar,this%pold,this%cold,this%ptmp,this%ctmp, &
         this%geo_src,this%src,this%rhs,this%rhsold,&
         this%xfluxdy,this%yfluxdx,this%amax,this%tol_abs,this%shift,&
         this%dtmean,this%dtstddev)
#ifdef PARALLEL
    IF(this%fargo.EQ.1) &
      DEALLOCATE(this%buf)
#endif
    IF(ASSOCIATED(this%bflux)) &
        DEALLOCATE(this%bflux)
    IF(ASSOCIATED(this%error)) &
      DEALLOCATE(this%error)
    CALL CloseTimedisc_common(this)
  END SUBROUTINE CloseTimedisc

END MODULE timedisc_generic
