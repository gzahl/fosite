!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: reconstruction_linear.f90                                         #
!#                                                                           #
!# Copyright (C) 2007-2008                                                   #
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
! module for linear (first order) reconstruction
!----------------------------------------------------------------------------!
MODULE reconstruction_linear
  USE reconstruction_constant
  USE common_types, ONLY : Common_TYP, InitCommon
  USE mesh_common, ONLY : Mesh_TYP
  USE physics_common, ONLY : Physics_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: MINMOD   = 1
  INTEGER, PARAMETER :: MONOCENT = 2
  INTEGER, PARAMETER :: SWEBY    = 3
  INTEGER, PARAMETER :: SUPERBEE = 4
  INTEGER, PARAMETER :: OSPRE    = 5
  CHARACTER(LEN=32), PARAMETER  :: recontype_name = "linear"  
  CHARACTER(LEN=32), DIMENSION(5), PARAMETER :: limitertype_name = (/ &
         "minmod  ", "mc      ", "sweby   ", "superbee", "ospre   " /)
  REAL, PARAMETER :: TINY = 1.0E-30 
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! constants
       MINMOD, MONOCENT, SWEBY, SUPERBEE, OSPRE, &
       ! methods
       InitReconstruction_linear, &
       MallocReconstruction_linear, &
       GetLimiterName, &
       PrimRecon, &
       CalculateSlopes_linear, &
       CalculateStates_linear, &
       CloseReconstruction_linear
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitReconstruction_linear(this,rtype,pc,ltype,lparam)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP) :: this
    INTEGER                  :: rtype,ltype
    REAL                     :: lparam
    LOGICAL                  :: pc
    !------------------------------------------------------------------------!
    INTENT(IN)               :: rtype,ltype,pc,lparam
    INTENT(INOUT)            :: this
    !------------------------------------------------------------------------!
    CALL InitReconstruction(this,rtype,recontype_name,pc)
    ! limiter settings
    CALL InitCommon(this%limiter,ltype,limitertype_name(ltype))
    this%limiter_param= lparam
  END SUBROUTINE InitReconstruction_linear


  SUBROUTINE MallocReconstruction_linear(this,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP) :: this
    TYPE(Mesh_TYP)           :: Mesh
    TYPE(Physics_TYP)        :: Physics
    !------------------------------------------------------------------------!
    INTEGER       :: err
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    ! allocate memory for all arrays used in reconstruction_linear
    ALLOCATE(this%xslopes(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
         this%yslopes(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
         this%temp1(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
         this%temp2(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
         this%temp3(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
         this%temp4(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
         STAT = err)
    IF (err.NE.0) THEN
       CALL Error(this, "MallocReconstruction_linear",  "Unable to allocate memory.")
    END IF
    ! zero the slopes
    this%xslopes(:,:,:) = 0.0
    this%yslopes(:,:,:) = 0.0
  END SUBROUTINE MallocReconstruction_linear


  PURE FUNCTION GetLimiter(this) RESULT(lt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP), INTENT(IN) :: this
    INTEGER :: lt
    !------------------------------------------------------------------------!
    lt = GetType(this%limiter)
  END FUNCTION GetLimiter


  PURE FUNCTION GetLimiterName(this) RESULT(ln)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: ln
    !------------------------------------------------------------------------!
    ln = GetName(this%limiter)
  END FUNCTION GetLimiterName


  PURE SUBROUTINE CalculateSlopes_linear(this,Mesh,Physics,rvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP) :: this
    TYPE(Mesh_TYP)           :: Mesh
    TYPE(Physics_TYP)        :: Physics
    REAL    :: rvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,rvar
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!

    ! choose limiter & check for 1D case ... not beautiful but currently best
    SELECT CASE(GetLimiter(this))
    CASE(MINMOD)
       ! calculate slopes in x-direction
       IF (Mesh%INUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR NOUNROLL
             DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
                DO i=Mesh%IMIN-1,Mesh%IMAX+1
                   this%xslopes(i,j,k) = Mesh%invdx * minmod2_limiter(&
                      rvar(i,j,k) - rvar(i-1,j,k), rvar(i+1,j,k) - rvar(i,j,k))
                END DO
             END DO
          END DO
       END IF
       ! calculate slopes in y-direction
       IF (Mesh%JNUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR NOUNROLL
             DO j=Mesh%JMIN-1,Mesh%JMAX+1
!CDIR NODEP
                DO i=Mesh%IGMIN,Mesh%IGMAX
                   this%yslopes(i,j,k) = Mesh%invdy * minmod2_limiter(&
                      rvar(i,j,k) - rvar(i,j-1,k), rvar(i,j+1,k) - rvar(i,j,k))
                END DO
             END DO
          END DO
       END IF
    CASE(MONOCENT)
       ! calculate slopes in x-direction
       IF (Mesh%INUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR NOUNROLL
             DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
                DO i=Mesh%IMIN-1,Mesh%IMAX+1
                   this%xslopes(i,j,k) = Mesh%invdx * minmod3_limiter(&
                      this%limiter_param*(rvar(i,j,k) - rvar(i-1,j,k)),&
                      this%limiter_param*(rvar(i+1,j,k) - rvar(i,j,k)),&
                      0.5*(rvar(i+1,j,k) - rvar(i-1,j,k)))
                END DO
             END DO
          END DO
       END IF
       ! calculate slopes in y-direction
       IF (Mesh%JNUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR NOUNROLL
             DO j=Mesh%JMIN-1,Mesh%JMAX+1
!CDIR NODEP
                DO i=Mesh%IGMIN,Mesh%IGMAX
                   this%yslopes(i,j,k) = Mesh%invdy * minmod3_limiter(&
                      this%limiter_param*(rvar(i,j,k)- rvar(i,j-1,k)),&
                      this%limiter_param*(rvar(i,j+1,k) - rvar(i,j,k)),&
                      0.5*(rvar(i,j+1,k) - rvar(i,j-1,k)))
                END DO
             END DO
          END DO
       END IF
    CASE(SWEBY)
       ! calculate slopes in x-direction
       IF (Mesh%INUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR NOUNROLL
             DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
                DO i=Mesh%IMIN-1,Mesh%IMAX+1
                   this%xslopes(i,j,k) = Mesh%invdx * sweby_limiter(&
                      rvar(i,j,k) - rvar(i-1,j,k), rvar(i+1,j,k) - rvar(i,j,k),&
                      this%limiter_param)
                END DO
             END DO
          END DO
       END IF
       ! calculate slopes in y-direction
       IF (Mesh%JNUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR NOUNROLL
             DO j=Mesh%JMIN-1,Mesh%JMAX+1
!CDIR NODEP
                DO i=Mesh%IGMIN,Mesh%IGMAX
                   this%yslopes(i,j,k) = Mesh%invdy * sweby_limiter(&
                      rvar(i,j,k) - rvar(i,j-1,k), rvar(i,j+1,k) - rvar(i,j,k),&
                      this%limiter_param)
                END DO
             END DO
          END DO
       END IF
    CASE(SUPERBEE)
       ! calculate slopes in x-direction
       IF (Mesh%INUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR NOUNROLL
             DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
                DO i=Mesh%IMIN-1,Mesh%IMAX+1
                   this%xslopes(i,j,k) = Mesh%invdx * sweby_limiter(&
                      rvar(i,j,k) - rvar(i-1,j,k), rvar(i+1,j,k) - rvar(i,j,k), 2.0)
                END DO
             END DO
          END DO
       END IF
       ! calculate slopes in y-direction
       IF (Mesh%JNUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR NOUNROLL
             DO j=Mesh%JMIN-1,Mesh%JMAX+1
!CDIR NODEP
                DO i=Mesh%IGMIN,Mesh%IGMAX
                   this%yslopes(i,j,k) = Mesh%invdy * sweby_limiter(&
                      rvar(i,j,k) - rvar(i,j-1,k), rvar(i,j+1,k) - rvar(i,j,k), 2.0)
                END DO
             END DO
          END DO
       END IF
    CASE(OSPRE)
       ! calculate slopes in x-direction
       IF (Mesh%INUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR NOUNROLL
             DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
                DO i=Mesh%IMIN-1,Mesh%IMAX+1
                   this%xslopes(i,j,k) = Mesh%invdx * ospre_limiter(&
                      rvar(i,j,k) - rvar(i-1,j,k), rvar(i+1,j,k) - rvar(i,j,k))
                END DO
             END DO
          END DO
       END IF
       ! calculate slopes in y-direction
       IF (Mesh%JNUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR NOUNROLL
             DO j=Mesh%JMIN-1,Mesh%JMAX+1
!CDIR NODEP
                DO i=Mesh%IGMIN,Mesh%IGMAX
                   this%yslopes(i,j,k) = Mesh%invdy * ospre_limiter(&
                      rvar(i,j,k) - rvar(i,j-1,k), rvar(i,j+1,k) - rvar(i,j,k))
                END DO
             END DO
          END DO
       END IF
    END SELECT
    
    CONTAINS
      
      ELEMENTAL FUNCTION minmod2_limiter(arg1,arg2) RESULT(limarg)
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        REAL             :: limarg
        REAL, INTENT(IN) :: arg1, arg2
        !--------------------------------------------------------------------!
        IF (arg1*arg2.GT.0) THEN
           limarg = SIGN(1.0,arg1) * MIN(ABS(arg1),ABS(arg2))
        ELSE
           limarg = 0.
        END IF
      END FUNCTION minmod2_limiter
      
      ELEMENTAL FUNCTION minmod3_limiter(arg1,arg2,arg3) RESULT(limarg)
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        REAL             :: limarg
        REAL, INTENT(IN) :: arg1, arg2, arg3
        !--------------------------------------------------------------------!
        IF (((arg1*arg2).GT.0).AND.((arg2*arg3).GT.0)) THEN
           limarg = SIGN(1.0,arg1) * MIN(ABS(arg1),ABS(arg2),ABS(arg3))
        ELSE
           limarg = 0.
        END IF
      END FUNCTION minmod3_limiter
      
      ELEMENTAL FUNCTION sweby_limiter(arg1,arg2,param) RESULT(limarg)
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        REAL             :: limarg
        REAL, INTENT(IN) :: arg1, arg2, param
        !--------------------------------------------------------------------!
        IF (arg1*arg2.GT.0) THEN
           limarg = SIGN(1.0,arg1) * MAX(MIN(param*ABS(arg1),ABS(arg2)), &
                MIN(ABS(arg1),param*ABS(arg2)))
        ELSE
           limarg = 0.
        END IF
      END FUNCTION sweby_limiter
      
      ELEMENTAL FUNCTION ospre_limiter(arg1,arg2) RESULT(limarg)
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        REAL             :: limarg
        REAL, INTENT(IN) :: arg1, arg2
        !--------------------------------------------------------------------!
        limarg = 1.5*arg1*arg2*(arg1 + arg2) / (arg1*(arg1 + 0.5*arg2) &
             + arg2*(arg2 + 0.5*arg1) + TINY)
      END FUNCTION ospre_limiter
      
  END SUBROUTINE CalculateSlopes_linear
  

  PURE SUBROUTINE CalculateStates_linear(this,Mesh,Physics,npos,pos0,pos,rvar,rstates)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP) :: this
    TYPE(Mesh_TYP)           :: Mesh
    TYPE(Physics_TYP)        :: Physics
    INTEGER                  :: npos
    REAL :: pos0(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2)
    REAL :: pos(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,npos,2)
    REAL :: rvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)
    REAL :: rstates(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,npos,Physics%vnum)
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k,n
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,Physics,npos,pos0,pos,rvar
    INTENT(OUT) :: rstates
    !------------------------------------------------------------------------!
    ! reconstruct states at positions pos 
    FORALL (i=Mesh%IGMIN:Mesh%IGMAX, j=Mesh%JGMIN:Mesh%JGMAX, n=1:npos, k=1:Physics%vnum)
       rstates(i,j,n,k) = reconstruct(rvar(i,j,k),this%xslopes(i,j,k), &
            this%yslopes(i,j,k),pos0(i,j,1),pos0(i,j,2),pos(i,j,n,1),pos(i,j,n,2))
    END FORALL

    CONTAINS

      ELEMENTAL FUNCTION reconstruct(cvar0,xslope0,yslope0,x0,y0,x,y) RESULT(rstate)
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        REAL :: rstate
        REAL, INTENT(IN) :: cvar0,xslope0,yslope0,x0,y0,x,y 
        !--------------------------------------------------------------------!
        rstate = cvar0 + xslope0*(x-x0) + yslope0*(y-y0)
      END FUNCTION reconstruct
      
  END SUBROUTINE CalculateStates_linear


  SUBROUTINE CloseReconstruction_linear(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP)  :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)             :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%xslopes,this%yslopes,&
         this%temp1,this%temp2,this%temp3,this%temp4)
  END SUBROUTINE CloseReconstruction_linear
  
END MODULE reconstruction_linear
