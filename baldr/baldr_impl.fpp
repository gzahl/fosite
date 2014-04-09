!#############################################################################
!#                                                                           #
!# baldr - 2D hydrodynamical simulation program                              #
!# module: baldr_impl.fpp                                                    #
!#                                                                           #
!# Copyright (C) 2011,2012                                                   #
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

#define REAL DOUBLE PRECISION

MODULE baldr
  USE fosite_common, ONLY : Fosite_TYP
  USE fosite, ONLY : InitFosite, StepFosite, RunFosite, CloseFosite, &
                     InitPhysics, InitFluxes, InitReconstruction, InitMesh, &
                     InitBoundary, InitSources, InitTimedisc, InitFileIO
  USE physics_generic, ONLY : Convert2Conservative
  USE sources_generic, ONLY : GetSourcesPointer
  USE geometry_generic, ONLY : Convert2Curvilinear
  USE boundary_generic, ONLY : NORTH, SOUTH, WEST, EAST
  USE sources_generic, ONLY : POISSON
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  TYPE(Fosite_TYP),ALLOCATABLE &
                         :: Sim
  !f2py intent(hide)     :: Sim 
  !--------------------------------------------------------------------------!
  ! flags from common/physics_common.f90
  INTEGER                :: VNUM
  INTEGER                :: DENSITY
  INTEGER                :: PRESSURE, ENERGY
  INTEGER                :: XVELOCITY, XMOMENTUM
  INTEGER                :: YVELOCITY, YMOMENTUM
  INTEGER                :: ZVELOCITY, ZMOMENTUM
  
  !--------------------------------------------------------------------------!
  INTEGER                  :: IMIN,IMAX,JMIN,JMAX,INUM,JNUM,IGMIN,IGMAX,&
                              JGMIN,JGMAX,GNUM
  REAL, DIMENSION(:,:,:), ALLOCATABLE &
                           :: curv_coords
  REAL, DIMENSION(:,:), ALLOCATABLE &
                           :: radius


  !--------------------------------------------------------------------------!


#define GET0D(varname, var, type)                                             @@\
  SUBROUTINE get##varname(varname)                                            @@\
    IMPLICIT NONE                                                             @@\
  !--------------------------------------------------------------------------!@@\
    type              :: varname                                              @@\
    !f2py intent(out) :: varname                                              @@\
  !--------------------------------------------------------------------------!@@\
  varname = var                                                               @@\
    END SUBROUTINE get##varname

#define SET0D(varname, var, type)                                             @@\
  SUBROUTINE set##varname(value)                                              @@\
    IMPLICIT NONE                                                             @@\
  !--------------------------------------------------------------------------!@@\
    type                :: value                                              @@\
    !f2py intent(in)    :: value                                              @@\
  !--------------------------------------------------------------------------!@@\
  var = value                                                                 @@\
    END SUBROUTINE set##varname

#define GET1D(varname, var)                                                   @@\
  SUBROUTINE get##varname(m,varname)                                          @@\
    IMPLICIT NONE                                                             @@\
  !--------------------------------------------------------------------------!@@\
    INTEGER             :: m                                                  @@\
    REAL, DIMENSION(m)&                                                       @@\
                        :: varname                                            @@\
    !f2py intent(out), depend(m) :: varname                                   @@\
  !--------------------------------------------------------------------------!@@\
  varname = var                                                               @@\
    END SUBROUTINE get##varname

#define SET1D(varname, var)                                                   @@\
  SUBROUTINE set##varname(m,value)                                            @@\
    IMPLICIT NONE                                                             @@\
  !--------------------------------------------------------------------------!@@\
    INTEGER             :: m                                                  @@\
    REAL, DIMENSION(m)&                                                       @@\
                        :: value                                              @@\
    !f2py intent(hide)  :: m                                                  @@\
    !f2py intent(in), depend(m) :: value                                      @@\
  !--------------------------------------------------------------------------!@@\
  var = value                                                                 @@\
    END SUBROUTINE set##varname

#define GET2D(varname, var)                                                   @@\
  SUBROUTINE get##varname(m,n,varname)                                        @@\
    IMPLICIT NONE                                                             @@\
  !--------------------------------------------------------------------------!@@\
    INTEGER             :: m, n                                               @@\
    REAL, DIMENSION(m,n)&                                                     @@\
                        :: varname                                            @@\
    !f2py intent(out), depend(m, n) :: varname                                @@\
  !--------------------------------------------------------------------------!@@\
  varname = var                                                               @@\
    END SUBROUTINE get##varname

#define SET2D(varname, var)                                                   @@\
  SUBROUTINE set##varname(m,n,value)                                          @@\
    IMPLICIT NONE                                                             @@\
  !--------------------------------------------------------------------------!@@\
    INTEGER             :: m, n                                               @@\
    REAL, DIMENSION(m,n)&                                                     @@\
                        :: value                                              @@\
    !f2py intent(hide)  :: m, n                                               @@\
    !f2py intent(in), depend(m, n) :: value                                   @@\
  !--------------------------------------------------------------------------!@@\
  var = value                                                                 @@\
    END SUBROUTINE set##varname

#define GET3D(varname, var)                                                   @@\
  SUBROUTINE get##varname(m,n,o,res)                                          @@\
    IMPLICIT NONE                                                             @@\
  !--------------------------------------------------------------------------!@@\
    INTEGER             :: m, n, o                                            @@\
    REAL, DIMENSION(m,n,o)&                                                   @@\
                        :: res                                                @@\
    !f2py intent(out), depend(m, n, o) :: res                                 @@\
  !--------------------------------------------------------------------------!@@\
  res = var                                                                   @@\
    END SUBROUTINE get##varname

#define SET3D(varname, var)                                                   @@\
  SUBROUTINE set##varname(m,n,o,value)                                        @@\
    IMPLICIT NONE                                                             @@\
  !--------------------------------------------------------------------------!@@\
    INTEGER             :: m, n, o                                            @@\
    REAL, DIMENSION(m,n,o)&                                                   @@\
                        :: value                                              @@\
    !f2py intent(hide)  :: m, n, o                                            @@\
    !f2py intent(in), depend(m, n, o) :: value                                @@\
  !--------------------------------------------------------------------------!@@\
  var = value                                                                 @@\
    END SUBROUTINE set##varname

#define GETSHAPE(varname, var, dim)                                           @@\
  SUBROUTINE shape##varname(res)                                              @@\
    IMPLICIT NONE                                                             @@\
  !--------------------------------------------------------------------------!@@\
    INTEGER, DIMENSION(dim)&                                                  @@\
                        :: res                                                @@\
    !f2py intent(out)   :: res                                                @@\
  !--------------------------------------------------------------------------!@@\
  res = SHAPE(var)                                                            @@\
    END SUBROUTINE shape##varname

#define GETSHAPE0D(varname, var)                                              @@\
  SUBROUTINE shape##varname(res)                                              @@\
    IMPLICIT NONE                                                             @@\
  !--------------------------------------------------------------------------!@@\
    INTEGER, DIMENSION(1)&                                                    @@\
                        :: res                                                @@\
    !f2py intent(out)   :: res                                                @@\
  !--------------------------------------------------------------------------!@@\
  res = 1                                                                     @@\
    END SUBROUTINE shape##varname

#define ADD0D(varname, var, type)                                             @@\
  GET0D(varname, var, type)                                                   @@\
  SET0D(varname, var, type)                                                   @@\
  GETSHAPE0D(varname, var)

#define ADD1D(varname, var)                                                   @@\
  GET1D(varname, var)                                                         @@\
  SET1D(varname, var)                                                         @@\
  GETSHAPE(varname, var, 1)

#define ADD2D(varname, var)                                                   @@\
  GET2D(varname, var)                                                         @@\
  SET2D(varname, var)                                                         @@\
  GETSHAPE(varname, var, 2)

#define ADD3D(varname, var)                                                   @@\
  GET3D(varname, var)                                                         @@\
  SET3D(varname, var)                                                         @@\
  GETSHAPE(varname, var, 3)


  !--------------------------------------------------------------------------!
  CONTAINS
  !--------------------------------------------------------------------------!

  SUBROUTINE InitBaldr()
    IMPLICIT NONE
  !--------------------------------------------------------------------------!

    IF(.NOT.ALLOCATED(Sim)) THEN
        ALLOCATE(Sim)
        CALL InitFosite(Sim)
    ELSE
        print *,"Fosite already allocated. Can't initialize baldr"
    END IF

    END SUBROUTINE InitBaldr

  SUBROUTINE StepBaldr(halt)
    IMPLICIT NONE
  !--------------------------------------------------------------------------!
    LOGICAL      :: halt
    !f2py intent(out) halt
  !--------------------------------------------------------------------------!

    halt = StepFosite(Sim)

    END SUBROUTINE StepBaldr

  SUBROUTINE Convert2ConservativeBaldr()
  !--------------------------------------------------------------------------!
    ! transform to conservative variables
    CALL Convert2Conservative(Sim%Physics,Sim%Mesh,&
                              Sim%Timedisc%pvar,Sim%Timedisc%cvar)

  END SUBROUTINE Convert2ConservativeBaldr

  SUBROUTINE RunBaldr()
    IMPLICIT NONE
  !--------------------------------------------------------------------------!

    CALL RunFosite(Sim)

    END SUBROUTINE RunBaldr

  SUBROUTINE CloseBaldr()
    IMPLICIT NONE
  !--------------------------------------------------------------------------!

    IF(ALLOCATED(Sim)) THEN
        CALL CloseFosite(Sim)

        IF(ALLOCATED(curv_coords)) DEALLOCATE(curv_coords)
        IF(ALLOCATED(radius)) DEALLOCATE(radius)

        IF(ALLOCATED(Sim)) DEALLOCATE(Sim)
    END IF

    END SUBROUTINE CloseBaldr


  SUBROUTINE InitPhysicsBaldr(problem,units,gamma,mu,cs, &
                             rhomin,pmin,dpmax)
    IMPLICIT NONE
  !--------------------------------------------------------------------------!
    INTEGER                         :: problem
    INTEGER, OPTIONAL               :: units
    REAL, OPTIONAL      :: gamma,mu,cs,rhomin,pmin,dpmax
    !f2py integer,intent(in)        :: problem
    !f2py integer,intent(in),optional  :: units = 1    ! SI
    !f2py double precision,intent(in),optional  :: gamma = 1.4
    !f2py double precision,intent(in),optional  :: mu = 0.029
    !f2py double precision,intent(in),optional  :: cs = 343.0
    !f2py double precision,intent(in),optional  :: rhomin = 1.0E-30
    !f2py double precision,intent(in),optional  :: pmin = 1.0E-30
    !f2py double precision,intent(in),optional  :: dpmax = 1.0
  !--------------------------------------------------------------------------!

    CALL InitPhysics(Sim%Physics,&
                     Sim%Mesh, &
                     problem = problem, &
                     units = units, &
                     gamma = gamma, &
                     mu = mu, &
                     cs = cs, &
                     rhomin = rhomin, &
                     pmin = pmin, &
                     dpmax = dpmax)


    VNUM = Sim%Physics%VNUM
    DENSITY = Sim%Physics%DENSITY-1
    PRESSURE = Sim%Physics%PRESSURE-1
    ENERGY = Sim%Physics%ENERGY-1
    XVELOCITY = Sim%Physics%XVELOCITY-1
    XMOMENTUM = Sim%Physics%XMOMENTUM-1
    YVELOCITY = Sim%Physics%YVELOCITY-1
    YMOMENTUM = Sim%Physics%YMOMENTUM-1
    ZVELOCITY = Sim%Physics%ZVELOCITY-1
    ZMOMENTUM = Sim%Physics%ZMOMENTUM-1
    
    END SUBROUTINE InitPhysicsBaldr

  SUBROUTINE InitFluxesBaldr(order,variables,limiter,theta)
    IMPLICIT NONE
  !--------------------------------------------------------------------------!
    INTEGER                             :: order
    LOGICAL                             :: variables
    INTEGER, OPTIONAL                   :: limiter
    REAL, OPTIONAL                      :: theta
    !f2py integer,intent(in)            :: order
    !f2py logical,intent(in)            :: variables
    !f2py integer,intent(in),optional   :: limiter = 1 ! MINMOD
    !f2py double precision,intent(in),optional      :: theta = 1.0
  !--------------------------------------------------------------------------!

    CALL InitFluxes(Sim%Fluxes, &
                    Sim%Mesh, &
                    Sim%Physics, &
                    order, &
                    variables, &
                    limiter = limiter, &
                    theta = theta)
    END SUBROUTINE InitFluxesBaldr

  SUBROUTINE InitMeshBaldr(geometry,meshtype,inum,jnum,xmin,xmax,ymin,ymax,gparam)
    IMPLICIT NONE
  !--------------------------------------------------------------------------!
    INTEGER                         :: geometry,meshtype
    INTEGER                         :: inum,jnum
    REAL                            :: xmin,xmax,ymin,ymax
    REAL, OPTIONAL                  :: gparam
    !f2py integer,intent(in)        :: geometry,meshtype
    !f2py integer,intent(in)        :: inum,jnum
    !f2py double precision,intent(in)           :: xmin,xmax,ymin,ymax
    !f2py double precision,intent(in),optional  :: gparam = 1.0
  !--------------------------------------------------------------------------!
      

    CALL InitMesh(Sim%Mesh, &
                  meshtype, &
                  geometry, &
                  inum, jnum, &
                  xmin, xmax, ymin, ymax, &
                  gparam = gparam)

    IGMIN = 0
    JGMIN = 0
    IMIN = IGMIN + Sim%Mesh%GNUM
    JMIN = JGMIN + Sim%Mesh%GNUM
    IMAX = IMIN + Sim%Mesh%INUM - 1
    JMAX = JMIN + Sim%Mesh%JNUM - 1
    IGMAX = IMAX + Sim%Mesh%GNUM
    JGMAX = JMAX + Sim%Mesh%GNUM
    INUM = Sim%Mesh%INUM
    JNUM = Sim%Mesh%JNUM
    GNUM = Sim%Mesh%GNUM

    ALLOCATE(radius(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,&
                    Sim%Mesh%JGMIN:Sim%Mesh%JGMAX))

    radius = SQRT(Sim%Mesh%bccart(:,:,1)**2 + Sim%Mesh%bccart(:,:,2)**2)

    ALLOCATE(curv_coords(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX, &
                         Sim%Mesh%JGMIN:Sim%Mesh%JGMAX,1:2))

    CALL Convert2Curvilinear(Sim%Mesh%Geometry,&
                             Sim%Mesh%bccart,&
                             curv_coords)

    END SUBROUTINE InitMeshBaldr

  SUBROUTINE InitBoundaryBaldr(western,eastern,southern,northern)
    IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTEGER                   :: western,eastern,southern,northern
  !f2py integer,intent(in)  :: western,eastern,southern,northern
  !--------------------------------------------------------------------------!
  
  CALL InitBoundary(Sim%Timedisc%boundary, &
                    Sim%Mesh, &
                    Sim%Physics, &
                    western, &
                    eastern, &
                    southern, &
                    northern)

    END SUBROUTINE InitBoundaryBaldr

  SUBROUTINE InitSourcesBaldr(stype,potential,vismodel,mass,mdot,rin,rout,&
                             dynconst,bulkconst,cvis,xaccel,yaccel,solver,&
                             maxresidnorm,maxmult,bndrytype,relaxtype,npre,&
                             npost,minres,nmaxcycle,omega,mass1,mass2,&
                             excentricity,semimayoraxis,green,sigma,x,y)
    IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTEGER           :: stype
  INTEGER, OPTIONAL :: potential,vismodel,solver,maxmult,bndrytype,relaxtype,&
                       npre,npost,minres,nmaxcycle,&
                       green
  REAL, OPTIONAL    :: mass,mdot,rin,rout,dynconst,bulkconst,cvis,&
                       xaccel,yaccel,maxresidnorm,&
                       omega,mass1,mass2,excentricity,semimayoraxis,&
                       sigma,x,y
  !f2py integer,intent(in)          :: stype
  !f2py integer,intent(in),optional :: potential,vismodel,solver,maxmult
  !f2py integer,intent(in),optional :: bndrytype,relaxtype,npre,npost
  !f2py integer,intent(in),optional :: minres,nmaxcycle
  !f2py double precision,intent(in),optional    :: mass,mdot,rin,rout,dynconst,bulkconst
  !f2py double precision,intent(in),optional    :: cvis,xaccel,yaccel,maxresidnorm
  !--------------------------------------------------------------------------!
  ! TODO: Define default values for optional arguments.

  CALL InitSources(Sim%Physics%sources,&
                   Sim%Mesh,&
                   Sim%Fluxes,&
                   Sim%Physics,&
                   Sim%Timedisc%Boundary,&
                   stype,&
                   potential = potential,&
                   vismodel = vismodel,&
                   mass = mass,&
                   mdot = mdot,&
                   rin = rin,&
                   rout = rout,&
                   dynconst = dynconst,&
                   bulkconst = bulkconst,&
                   cvis = cvis,&
                   xaccel = xaccel,&
                   yaccel = yaccel,&
                   solver = solver,&
                   maxresidnorm = maxresidnorm,&
                   maxmult = maxmult,&
                   bndrytype = bndrytype,&
                   relaxtype = relaxtype,&
                   npre = npre,&
                   npost = npost,&
                   minres = minres,&
                   nmaxcycle = nmaxcycle)


  END SUBROUTINE InitSourcesBaldr

  SUBROUTINE InitTimediscBaldr(method,order,stoptime,cfl,dtlimit,maxiter,&
                              tol_rel)!,&
                              !tol_abs)
    IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTEGER                       :: method
  INTEGER                       :: order
  REAL              :: stoptime
  REAL, OPTIONAL    :: cfl,dtlimit,tol_rel!,tol_abs(VNUM)
  INTEGER, OPTIONAL             :: maxiter
  !--------------------------------------------------------------------------!
  INTENT(IN)                    :: method,order,stoptime,cfl,dtlimit,tol_rel,&
                                   maxiter!,tol_abs  
  !--------------------------------------------------------------------------!
  !f2py integer,intent(in)                      :: method
  !f2py integer,intent(in)                      :: order
  !f2py double precision,intent(in)             :: stoptime
  !f2py double precision,intent(in),optional    :: cfl = 0.4
  !f2py double precision,intent(in),optional    :: dtlimit = -1.0
  !f2py double precision,intent(in),optional    :: maxiter = -1.0
  !f2py double precision,intent(in),optional    :: tol_rel = 0.01
  !!f2py double precision,intent(in),optional    :: tol_abs
  !--------------------------------------------------------------------------!

  !IF(dtlimit.EQ.-1.0) dtlimit = stoptime
  !IF(maxiter.EQ.-1.0) maxiter = HUGE(maxiter)
  ! tol_abs(:) = 0.001 if not defined...

  CALL InitTimedisc(Sim%Timedisc,&
                    Sim%Mesh,& 
                    Sim%Physics,&
                    method,&
                    order,&
                    stoptime,&
                    cfl,&
                    dtlimit,&
                    maxiter,&
                    tol_rel)!,&
                    !tol_abs)

  END SUBROUTINE InitTimediscBaldr

  SUBROUTINE InitDatafileBaldr(fileformat,filename,stoptime,dtwall,count,&
                              filecycles,ncfmt)
                              !,unit)
    IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTEGER               :: fileformat
  CHARACTER(LEN=*)      :: filename
  REAL, OPTIONAL        :: stoptime
  INTEGER, OPTIONAL     :: dtwall, count, filecycles, ncfmt!, unit
  !f2py integer,intent(in)          :: fileformat
  !f2py character,intent(in)        :: filename
  !f2py double precision,intent(in),optional    :: stoptime = -1
  !f2py integer,intent(in),optional :: dtwall = 3600
  !f2py integer,intent(in),optional :: count = 1
  !f2py integer,intent(in),optional :: filecycles = -1
  !f2py integer,intent(in),optional :: ncfmt = 0 
  !!f2py integer,intent(in),optional :: unit = 
  !--------------------------------------------------------------------------!

    IF(filecycles==-1) filecycles = count+1
    IF(stoptime==-1) stoptime = Sim%Timedisc%stoptime
   
    CALL InitFileIO(Sim%Datafile, &
                    Sim%Mesh, &
                    Sim%Physics, &
                    Sim%Timedisc, &
                    fileformat, &
                    filename, &
                    stoptime, &
                    dtwall, &
                    count, &
                    filecycles, &
                    ncfmt)!, &
                    !unit)

    END SUBROUTINE InitDatafileBaldr

  SUBROUTINE InitLogfileBaldr(fileformat,filename,stoptime,dtwall,count,&
                             filecycles,ncfmt)
                             !,unit)
    IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTEGER               :: fileformat
  CHARACTER(LEN=*)      :: filename
  REAL, OPTIONAL        :: stoptime
  INTEGER, OPTIONAL     :: dtwall, count, filecycles, ncfmt!, unit
  !f2py integer,intent(in)          :: fileformat
  !f2py character,intent(in)        :: filename
  !f2py double precision,intent(in),optional    :: stoptime = -1
  !f2py integer,intent(in),optional :: dtwall = 3600
  !f2py integer,intent(in),optional :: count = 1
  !f2py integer,intent(in),optional :: filecycles = -1
  !f2py integer,intent(in),optional :: ncfmt = 0 
  !!f2py integer,intent(in),optional :: unit = 
  !--------------------------------------------------------------------------!

    IF(filecycles==-1) filecycles = count+1
    IF(stoptime==-1) stoptime = Sim%Timedisc%stoptime
   
    CALL InitFileIO(Sim%Logfile, &
                    Sim%Mesh, &
                    Sim%Physics, &
                    Sim%Timedisc, &
                    fileformat, &
                    filename, &
                    stoptime, &
                    dtwall, &
                    count, &
                    filecycles, &
                    ncfmt)!, &
                    !unit)

    END SUBROUTINE InitLogfileBaldr

  
  SUBROUTINE NoOutboundBaldr
    IMPLICIT NONE
  !--------------------------------------------------------------------------!
  
    Sim%Physics%sources%outbound = 0

    END SUBROUTINE NoOutboundBaldr

#ifdef FFTW
    ADD2D(phi,Sim%Physics%sources%poisson%phi)
    ADD3D(accel,Sim%Physics%sources%poisson%accel)
    ADD3D(fi,Sim%Physics%sources%poisson%FI)
#endif
    ADD3D(pvar,Sim%Timedisc%pvar)
    ADD2D(dlx,Sim%Mesh%dlx)
    ADD2D(dly,Sim%Mesh%dly)
    ADD2D(bhx,Sim%Mesh%bhx)
    ADD2D(bhy,Sim%Mesh%bhy)
    ADD2D(bhz,Sim%Mesh%bhz)
    ADD3D(fhx,Sim%Mesh%fhx)
    ADD3D(fhy,Sim%Mesh%fhy)
    ADD3D(fhz,Sim%Mesh%fhz)
    ADD2D(volume,Sim%Mesh%volume)
    ADD0D(inum,Sim%Mesh%INUM,INTEGER)
    ADD0D(jnum,Sim%Mesh%JNUM,INTEGER)
    ADD0D(gnum,Sim%Mesh%GNUM,INTEGER)
    ADD3D(east,Sim%Timedisc%Boundary(EAST)%data)
    ADD3D(west,Sim%Timedisc%Boundary(WEST)%data)
    ADD3D(south,Sim%Timedisc%Boundary(SOUTH)%data)
    ADD3D(north,Sim%Timedisc%Boundary(NORTH)%data)
    ADD2D(cbeast,Sim%Timedisc%Boundary(EAST)%cbtype)
    ADD2D(cbwest,Sim%Timedisc%Boundary(WEST)%cbtype)
    ADD2D(cbsouth,Sim%Timedisc%Boundary(SOUTH)%cbtype)
    ADD2D(cbnorth,Sim%Timedisc%Boundary(NORTH)%cbtype)
    ADD3D(cart_coords,Sim%Mesh%bccart)
#ifdef HAVE_FFTW
    ADD3D(fi,Sim%Physics%sources%poisson%FI)
#endif
    ADD0D(gn,Sim%Physics%Constants%GN,REAL)

    ADD0D(time,Sim%Timedisc%time,REAL)
    ADD0D(dtmin,Sim%Timedisc%dtmin,REAL)
    ADD0D(dt,Sim%Timedisc%dt,REAL)
    ADD0D(iter,Sim%iter,INTEGER)
    

END MODULE baldr
