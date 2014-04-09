!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: gravity_monopol.f90                                               #
!#                                                                           #
!# Copyright (C) 2007-2010,2014                                              #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
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
!> \author Tobias Illenseer
!! \author Björn Sperling
!!
!! \brief source terms module for gravitational acceleration due to ...
!!
!! \extends gravity_common
!! \ingroup gravity
!----------------------------------------------------------------------------!
MODULE gravity_monopol
  USE gravity_common
  USE physics_generic
  USE mesh_generic
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: gravity_name = "monopol approx. for gravity"
  REAL, PARAMETER              :: SQRTTWOPI = SQRT(2.0*PI)
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Gravity_TYP, &
       ! constants
       ! methods
       InitGravity_monopol, &
       GetAccelGravity_monopol, &
       CloseGravity_monopol
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGravity_monopol(this,Mesh,Physics,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Dict_TYP),POINTER :: config
    INTEGER           :: gtype
    !------------------------------------------------------------------------!
    INTEGER           :: err
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "gtype", gtype)
    CALL InitGravity(this,gtype,gravity_name)

#ifdef PARALLEL
    CALL Error(this,"InitGravity_gravmonopol", &
         "parallel monopol approximation not supported yet")
#endif

    ALLOCATE(this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM), &
         this%radius(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%radius3(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%gposvecr3(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
         this%cellmass(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%enclmass(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%tmp(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         STAT = err)
    IF (err.NE.0) CALL Error(this,"InitGravity_gravmonopol", "Unable allocate memory!")

    ! compute distance to the origin for each cell
    this%radius(:,:) = SQRT(Mesh%bccart(:,:,1)**2 + Mesh%bccart(:,:,2)**2) 
    this%radius3(:,:) = this%radius(:,:)**3

    ! f_grav / mass = -GN*posvec/r**3
    this%gposvecr3(:,:,1) = -Physics%constants%GN *  Mesh%bposvec(:,:,1) &
         / (this%radius3(:,:)+TINY(this%radius(:,:)))
    this%gposvecr3(:,:,2) = this%gposvecr3(:,:,1) * Mesh%bposvec(:,:,2) / Mesh%bposvec(:,:,1)

    ! reset gravitational acceleration
    this%accel(:,:,:)  = 0.0

    this%cs => GetSoundSpeeds(Physics)

  END SUBROUTINE InitGravity_monopol


  PURE SUBROUTINE UpdateGravMonopol(this,Mesh,Physics,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,i0,j0
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,time,pvar
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    this%accel(:,:,:) = 0.0
    ! compute cell masses
!CDIR OUTERUNROLL=8
       DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
          DO i=Mesh%IMIN,Mesh%IMAX
             this%cellmass(i,j) = pvar(i,j,Physics%DENSITY) * Mesh%volume(i,j)
          END DO
       END DO
       ! compute cartesian acceleration
!CDIR OUTERUNROLL=8
       DO j0=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
          DO i0=Mesh%IMIN,Mesh%IMAX
          ! compute enclosed mass
             this%enclmass(i0,j0) = 0.0
             DO j=Mesh%JMIN,Mesh%JMAX
                DO i=Mesh%IMIN,Mesh%IMAX
                   IF (this%radius(i,j).LT.this%radius(i0,j0)) &
                        this%enclmass(i0,j0) = this%enclmass(i0,j0) + this%cellmass(i,j)
                END DO
             END DO
             ! compute components of gravitational acceleration
             this%accel(i0,j0,1:2) = this%enclmass(i0,j0) * this%gposvecr3(i0,j0,1:2)
          END DO
       END DO

  END SUBROUTINE UpdateGravMonopol

  FUNCTION GetAccelGravity_monopol(this,Mesh,Physics,time,pvar) RESULT(ac)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    REAL, DIMENSION(:,:,:),&
      POINTER         :: ac
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: gravptr
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,time,pvar
    INTENT(INOUT)     :: Physics
    !------------------------------------------------------------------------!
    CALL UpdateGravMonopol(this,Mesh,Physics,time,pvar)
    ac => this%accel
  END FUNCTION GetAccelGravity_monopol


  SUBROUTINE CloseGravity_monopol(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%accel,this%radius,this%radius3,this%gposvecr3, &
               this%cellmass,this%enclmass,this%tmp)
    CALL CloseGravity(this)
  END SUBROUTINE CloseGravity_monopol

END MODULE gravity_monopol
