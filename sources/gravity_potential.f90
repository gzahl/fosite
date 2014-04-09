!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: poisson_spectral.f90                                              #
!#                                                                           #
!# Copyright (C) 2011                                                        #
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
!> \addtogroup gravity
!! - parameters of \link gravity_potential \endlink as key-values
!! \todo description of n,s0 and sdelta 
!! \key{n,INTEGER,number of components?}
!! \key{s0,REAL,???,1.0}
!! \key{sdelta,REAL,???,0.0}
!! \key{switchon,REAL,soft switch on,0.}
!----------------------------------------------------------------------------!
!> \author Manuel Jung
!!
!! \brief evaluate gravitational forces by finite differences of constant
!! potentials
!!
!! This source module allows to define \f$n\f$ constant gravitational
!! potentials \f$\Phi_i\f$, \f$i\in n\f$. The gravitational acceleration is
!! then calculated by evaluation of
!! \f[
!!   \mathbf{a} = -\nabla \left(\sum_{i\in n} \Phi_i \right)
!! \f]
!!
!! Also the user is able to provide a linear scaling function, so non
!! axisymmetric potentials can be switched on smoothly.
!! It is defined as:
!! \f[
!!   \mathbf{f}(t)=\mathbf{s}_0 + \mathbf{s}_\Delta \cdot\text{min}\left(t/\tau,1\right)
!! \f]
!!
!! \extends gravity_common
!! \ingroup gravity
!----------------------------------------------------------------------------!
MODULE gravity_potential
  USE gravity_common
  USE mesh_generic
  USE physics_generic
  USE boundary_generic
  USE functions
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER &
                    :: solver_name  = "potential"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Gravity_TYP, &
       Grid_TYP, &
       ! constants
       ! methods
       InitGravity_potential, &
       GetAccelGravity_potential, &
       CloseGravity_potential
  !--------------------------------------------------------------------------!
  CONTAINS

  SUBROUTINE InitGravity_potential(this,Mesh,Physics,Boundary,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Boundary_TYP), DIMENSION(4) :: Boundary
    TYPE(Dict_TYP), POINTER :: config, IO
    !------------------------------------------------------------------------!
    INTEGER           :: solver,stype
    INTEGER           :: err, valwrite, i
    CHARACTER(LEN=32) :: info_str, no
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Boundary
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "gtype", solver)
    CALL InitGravity(this,solver,solver_name)

    CALL RequireKey(config, "n")
    CALL GetAttr(config, "n", this%n)

    ALLOCATE(this%mphi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%n), &
             this%s0(this%n), &
             this%sdelta(this%n), &
             this%lastfac(this%n), &
             this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM), &
             STAT = err)
    IF (err.NE.0) CALL Error(this, "InitGravity_potential", &
         "Unable to allocate memory.")

    this%accel(:,:,:) = 0.
    
    this%s0(:) = 1.
    CALL RequireKey(config, "s0", this%s0(:)) 
    CALL GetAttr(config, "s0", this%s0)

    this%sdelta(:) = 0.
    CALL RequireKey(config, "sdelta", this%sdelta(:))
    CALL GetAttr(config, "sdelta", this%sdelta)

    CALL RequireKey(config, "switchon", 0.)
    CALL GetAttr(config, "switchon", this%switchon)
    
    CALL Info(this, " POISSON--> constant potentials")

    WRITE (info_str, '(ES10.4)') this%switchon
    CALL Info(this, " POISSON--> switch on:         " // TRIM(info_str))

!    valwrite = 0
!    IF (HasKey(config, "output/mpotential")) &
!      CALL GetAttr(config, "output/mpotential", valwrite)
!    IF (valwrite .EQ. 1) THEN
!       DO i = 1,this%n
!         WRITE(no, "(I1)")i 
!         CALL AddField(IO, &
!                 "mpotential_" // TRIM(no), &
!                 remap_bounds(Mesh,this%mphi(:,:,:,i)), &
!                 Dict("name" / ("potential_" // TRIM(no))))
!       END DO
!    END IF

    this%lastfac(:) = 0.
  END SUBROUTINE InitGravity_potential


  FUNCTION GetAccelGravity_potential(this,Mesh,Physics,time) RESULT(ac)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Gravity_TYP), POINTER :: this
   TYPE(Mesh_TYP)       :: Mesh
   TYPE(Physics_TYP)    :: Physics
   REAL                 :: time
   REAL, DIMENSION(:,:,:), POINTER &
                        :: ac
   !------------------------------------------------------------------------!
   INTEGER              :: i, j, k
   REAL, DIMENSION(this%n) :: fac
   !------------------------------------------------------------------------!
   INTENT(IN)           :: Mesh,Physics,time
   !------------------------------------------------------------------------!

   fac = this%s0 &
      + MAX(MIN(time/(this%switchon+TINY(this%switchon)),1.),0.)*this%sdelta

   IF((SUM(ABS(this%lastfac-fac)).GT.TINY(time))) THEN
     this%accel(:,:,:) = 0.
     DO k = 1,this%n
       DO i = Mesh%IGMIN,Mesh%IGMAX
         DO j = Mesh%JGMIN,Mesh%JGMAX
           ! g(x) = - grad(phi(x)) 
           this%accel(i,j,1) = this%accel(i,j,1) &
            + fac(k)*(this%mphi(i,j,1,k)-this%mphi(i,j,2,k))/Mesh%dlx(i,j)
           this%accel(i,j,2) = this%accel(i,j,2) &
            + fac(k)*(this%mphi(i,j,3,k)-this%mphi(i,j,4,k))/Mesh%dly(i,j)
         END DO
       END DO
     END DO
   END IF
   this%lastfac = fac
   ac => this%accel
  END FUNCTION GetAccelGravity_potential


  SUBROUTINE CloseGravity_potential(this)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Gravity_TYP), POINTER :: this
   !------------------------------------------------------------------------!
   DEALLOCATE(this%mphi,this%s0,this%sdelta,this%accel)
  END SUBROUTINE CloseGravity_potential

END MODULE gravity_potential
