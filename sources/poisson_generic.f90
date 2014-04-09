!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: poisson_generic.f90                                               #
!#                                                                           #
!# Copyright (C) 2011                                                        #
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
! generic poisson solver module providing functionaly common
! to all source terms
!----------------------------------------------------------------------------!
MODULE poisson_generic
  USE poisson_multigrid
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  USE boundary_common, ONLY : Boundary_TYP
  USE sources_pointmass
  USE physics_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! flags for source terms
  INTEGER, PARAMETER :: MULTIGRID    = 1
!   INTEGER, PARAMETER :: MULTIPOL     = 2
!   INTEGER, PARAMETER :: FFT          = 3
  CHARACTER(LEN=32), PARAMETER :: poisson_name = "Poisson"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Poisson_TYP, &
       Grid_TYP, &
       ! constants
       MULTIGRID, &
       RED_BLACK_GAUSS_SEIDEL,BLOCK_GAUSS_SEIDEL,GAUSS_SEIDEL,&
       SPHERMULTEXPAN, CYLINMULTEXPAN, &
       ! methods
       InitPoisson, &
       PoissonSource, &
       ClosePoisson, &
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

  SUBROUTINE InitPoisson(this,Mesh,Physics,Boundary,stype,solver,maxresidnorm,maxmult,&
       bndrytype,relaxtype,npre,npost,minres,nmaxcycle)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER  :: this
    TYPE(Mesh_TYP)             :: Mesh
    TYPE(Physics_TYP)          :: Physics
    TYPE(Boundary_TYP), DIMENSION(4) :: Boundary
    INTEGER           :: stype
    INTEGER, OPTIONAL :: solver,maxmult,bndrytype,relaxtype,npre,npost,minres,nmaxcycle
    REAL, OPTIONAL    :: maxresidnorm
    !------------------------------------------------------------------------!
    INTEGER           :: solver_def,maxmult_def,bndrytype_def,relaxtype_def,npre_def,&
                         npost_def,minres_def,nmaxcycle_def
    REAL              :: maxresidnorm_def
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Boundary,stype,solver,maxresidnorm,maxmult,&
                         bndrytype,relaxtype,npre,npost,minres,nmaxcycle
    !------------------------------------------------------------------------!
    IF (.NOT.Initialized(Physics).OR..NOT.Initialized(Mesh)) &
         CALL Error(this,"InitPoisson","physics and/or mesh module uninitialized")
    CALL InitSources(this,stype,poisson_name)

    IF (PRESENT(solver)) THEN
          solver_def = solver
    ELSE
          solver_def = MULTIGRID
    END IF

    SELECT CASE(solver_def)
    CASE(MULTIGRID)
       ! number of multipol moments (in case of spherical grids)
       IF (PRESENT(maxmult)) THEN
          maxmult_def = maxmult
       ELSE
          maxmult_def = 5
       END IF
       ! accuracy of multigrid solver
       IF (PRESENT(maxresidnorm)) THEN
          maxresidnorm_def = maxresidnorm
       ELSE
          maxresidnorm_def = 1.0E-5
       END IF
       ! type of multipol expansion (spherical, cylindrical)
       IF (PRESENT(bndrytype)) THEN
          bndrytype_def = bndrytype
       ELSE
          bndrytype_def = CYLINMULTEXPAN
       END IF
       ! type of relaxation method
       IF (PRESENT(relaxtype)) THEN
          relaxtype_def = relaxtype
       ELSE
          relaxtype_def = BLOCK_GAUSS_SEIDEL
       END IF  
       ! number of pre smoothings
       IF (PRESENT(npre)) THEN
          npre_def = npre
       ELSE
          npre_def = 1
       END IF  
       ! number of post smoothings
       IF (PRESENT(npost)) THEN
          npost_def = npost
       ELSE
          npost_def = 1
       END IF  
       ! resolution of coarest grid
       IF (PRESENT(minres)) THEN
          minres_def = minres
       ELSE
          minres_def = 3
       END IF  
       ! max iteration steps of mg solver 
       IF (PRESENT(nmaxcycle)) THEN
          nmaxcycle_def = nmaxcycle
       ELSE
          nmaxcycle_def = 100
       END IF  
       CALL InitPoisson_multigrid(this%poisson,Mesh,Physics,Boundary,solver_def, &
            maxmult_def,maxresidnorm_def,bndrytype_def,relaxtype_def, &
            npre_def,npost_def,minres_def,nmaxcycle_def)
     CASE DEFAULT
       CALL Error(this,"InitPoisson","unknown Poisson type")
     END SELECT

     ! allocate common memory for all Poisson
     CALL MallocPoisson(this%poisson,Mesh,Physics)
  END SUBROUTINE InitPoisson

  SUBROUTINE MallocPoisson(this,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Poisson_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ALLOCATE(this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
             this%phi(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1:Mesh%JMAX+1), &  
             STAT = err)
    IF (err.NE.0) CALL Error(this, "MallocPoisson", &
         "Unable to allocate memory.")
    ! set gravitational acceleration
    this%accel(:,:,:) = 0.
    this%phi(:,:) = 0.
  END SUBROUTINE MallocPoisson

  SUBROUTINE PoissonSource(this,Mesh,Physics,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Poisson_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,pvar,cvar
    INTENT(INOUT)     :: sterm
    !------------------------------------------------------------------------!
    ! call specific subroutine
    SELECT CASE(GetType(this))
    CASE(MULTIGRID)
       CALL CalcPotential_multigrid(this,Mesh,Physics,pvar)
    CASE DEFAULT
       CALL Error(this,"PoissonSource", "unknown poisson term")
    END SELECT

    DO i = Mesh%IMIN,Mesh%IMAX
      DO j = Mesh%JMIN,Mesh%JMAX
        ! g(x) = - grad(phi(x)) 
        this%accel(i,j,2) = 0.5*(this%phi(i,j-1)-this%phi(i,j+1))/Mesh%dly(i,j)
        this%accel(i,j,1) = 0.5*(this%phi(i-1,j)-this%phi(i+1,j))/Mesh%dlx(i,j)
      END DO
    END DO

    !no acceleration in boundary cells 
    this%accel(Mesh%IGMIN:Mesh%IMIN-1,:,:) = 0.0
    this%accel(Mesh%IMAX+1:Mesh%IGMAX,:,:) = 0.0
    this%accel(:,Mesh%JGMIN:Mesh%JMIN-1,:) = 0.0
    this%accel(:,Mesh%JMAX+1:Mesh%JGMAX,:) = 0.0
    ! gravitational source terms
    CALL ExternalSources(Physics,Mesh,this%accel,pvar,cvar,sterm)
  END SUBROUTINE PoissonSource

  SUBROUTINE ClosePoisson(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Poisson_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%accel,this%phi)
    SELECT CASE(GetType(this))
    CASE(MULTIGRID)
       CALL ClosePoisson_multigrid(this)
    END SELECT
  END SUBROUTINE ClosePoisson
END MODULE poisson_generic
