!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: reconstruction_linear.f90                                         #
!#                                                                           #
!# Copyright (C) 2007-2012                                                   #
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
!! \brief module for linear (first order) reconstruction
!!
!! \extends reconstruction_common
!! \ingroup reconstruction
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
  INTEGER, PARAMETER :: PP       = 6
  INTEGER, PARAMETER :: VANLEER  = 7
  INTEGER, PARAMETER :: NOLIMIT  = 8
  CHARACTER(LEN=32), PARAMETER  :: recontype_name = "linear"  
  CHARACTER(LEN=32), DIMENSION(8), PARAMETER :: limitertype_name = (/ &
         "minmod    ", "mc        ", "sweby     ", "superbee  ", "ospre     ", &
         "pp        ", "van leer  ", "no limit  "/)
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! constants
       MINMOD, MONOCENT, SWEBY, SUPERBEE, OSPRE, PP, VANLEER, NOLIMIT, &
       ! methods
       InitReconstruction_linear, &
       GetLimiterName, &
       PrimRecon, &
       CalculateSlopes_linear, &
       CalculateStates_linear, &
       CloseReconstruction_linear
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitReconstruction_linear(this,Mesh,Physics,rtype,pc,ltype,lparam)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP) :: this
    TYPE(Mesh_TYP)           :: Mesh
    TYPE(Physics_TYP)        :: Physics
    INTEGER                  :: rtype,ltype
    REAL                     :: lparam
    INTEGER                  :: pc
    !------------------------------------------------------------------------!
    INTEGER                  :: err
    !------------------------------------------------------------------------!
    INTENT(IN)               :: Mesh,Physics,rtype,ltype,pc,lparam
    INTENT(INOUT)            :: this
    !------------------------------------------------------------------------!
    CALL InitReconstruction(this,rtype,recontype_name,pc)
    ! limiter settings
    SELECT CASE(ltype)
    CASE(MINMOD,MONOCENT,SWEBY,SUPERBEE,OSPRE,PP,VANLEER,NOLIMIT)
       CALL InitCommon(this%limiter,ltype,limitertype_name(ltype))
       this%limiter_param= lparam
    CASE DEFAULT
       CALL Error(this, "InitReconstruction_linear", "Unknown limiter")
    END SELECT

    ! allocate memory for all arrays used in reconstruction_linear
    ALLOCATE(this%xslopes(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
         this%yslopes(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
         STAT = err)
    IF (err.NE.0) THEN
       CALL Error(this, "InitReconstruction_linear",  "Unable to allocate memory.")
    END IF
    ! zero the slopes
    this%xslopes(:,:,:) = 0.0
    this%yslopes(:,:,:) = 0.0
  END SUBROUTINE InitReconstruction_linear


  PURE FUNCTION GetLimiter(this) RESULT(lt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP), INTENT(IN) :: this
    INTEGER :: lt
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    lt = GetType(this%limiter)
  END FUNCTION GetLimiter


  PURE FUNCTION GetLimiterName(this) RESULT(ln)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: ln
    !------------------------------------------------------------------------!
!CDIR IEXPAND
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

    ! choose limiter & check for 1D case
!CDIR IEXPAND
    SELECT CASE(GetLimiter(this))
    CASE(MINMOD)
       ! calculate slopes in x-direction
       IF (Mesh%INUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR UNROLL=8
             DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
                DO i=Mesh%IGMIN+1,Mesh%IGMAX-1
                   this%xslopes(i,j,k) = Mesh%invdx * minmod2_limiter(&
                      rvar(i,j,k) - rvar(i-1,j,k), rvar(i+1,j,k) - rvar(i,j,k))
                END DO
             END DO
          END DO
       END IF
       ! calculate slopes in y-direction
       IF (Mesh%JNUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR COLLAPSE
             DO j=Mesh%JGMIN+1,Mesh%JGMAX-1
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
!CDIR UNROLL=8
             DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
                DO i=Mesh%IGMIN+1,Mesh%IGMAX-1
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
!CDIR COLLAPSE
             DO j=Mesh%JGMIN+1,Mesh%JGMAX-1
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
!CDIR UNROLL=8
             DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
                DO i=Mesh%IGMIN+1,Mesh%IGMAX-1
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
!CDIR COLLAPSE
             DO j=Mesh%JGMIN+1,Mesh%JGMAX-1
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
!CDIR UNROLL=8
             DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
                DO i=Mesh%IGMIN+1,Mesh%IGMAX-1
                   this%xslopes(i,j,k) = Mesh%invdx * sweby_limiter(&
                      rvar(i,j,k) - rvar(i-1,j,k), rvar(i+1,j,k) - rvar(i,j,k), 2.0)
                END DO
             END DO
          END DO
       END IF
       ! calculate slopes in y-direction
       IF (Mesh%JNUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR COLLAPSE
             DO j=Mesh%JGMIN+1,Mesh%JGMAX-1
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
!CDIR UNROLL=8
             DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR COLLAPSE
                DO i=Mesh%IGMIN+1,Mesh%IGMAX-1
                   this%xslopes(i,j,k) = Mesh%invdx * ospre_limiter(&
                      rvar(i,j,k) - rvar(i-1,j,k), rvar(i+1,j,k) - rvar(i,j,k))
                END DO
             END DO
          END DO
       END IF
       ! calculate slopes in y-direction
       IF (Mesh%JNUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR COLLAPSE
             DO j=Mesh%JGMIN+1,Mesh%JGMAX-1
                DO i=Mesh%IGMIN,Mesh%IGMAX
                   this%yslopes(i,j,k) = Mesh%invdy * ospre_limiter(&
                      rvar(i,j,k) - rvar(i,j-1,k), rvar(i,j+1,k) - rvar(i,j,k))
                END DO
             END DO
          END DO
       END IF
       CASE(PP)
          ! calculate slopes in both-directions
          DO k=1,Physics%VNUM
!CDIR UNROLL=8
             DO j=Mesh%JGMIN+1,Mesh%JGMAX-1
!CDIR NODEP
                DO i=Mesh%IGMIN+1,Mesh%IGMAX-1
                   CALL pp_limiter(this%xslopes(i,j,k),&
                                   this%yslopes(i,j,k),&
                      rvar(i-1,j+1,k)-rvar(i,j,k),&
                      rvar(i  ,j+1,k)-rvar(i,j,k),&
                      rvar(i+1,j+1,k)-rvar(i,j,k),&
                      rvar(i-1,j  ,k)-rvar(i,j,k),&
                      1.0E-10,&
                      rvar(i+1,j  ,k)-rvar(i,j,k),&
                      rvar(i-1,j-1,k)-rvar(i,j,k),&
                      rvar(i  ,j-1,k)-rvar(i,j,k),&
                      rvar(i+1,j-1,k)-rvar(i,j,k))
                END DO
             END DO
          END DO
    CASE(VANLEER)
       ! calculate slopes in x-direction
       IF (Mesh%INUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR UNROLL=8
             DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR COLLAPSE
                DO i=Mesh%IGMIN+1,Mesh%IGMAX-1
                   this%xslopes(i,j,k) = Mesh%invdx * vanleer_limiter(&
                      rvar(i,j,k) - rvar(i-1,j,k), rvar(i+1,j,k) - rvar(i,j,k))
                END DO
             END DO
          END DO
       END IF
       ! calculate slopes in y-direction
       IF (Mesh%JNUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR COLLAPSE
             DO j=Mesh%JGMIN+1,Mesh%JGMAX-1
                DO i=Mesh%IGMIN,Mesh%IGMAX
                   this%yslopes(i,j,k) = Mesh%invdy * vanleer_limiter(&
                      rvar(i,j,k) - rvar(i,j-1,k), rvar(i,j+1,k) - rvar(i,j,k))
                END DO
             END DO
          END DO
       END IF
    CASE(NOLIMIT)
       ! calculate slopes in x-direction
       IF (Mesh%INUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR UNROLL=8
             DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
                DO i=Mesh%IGMIN+1,Mesh%IGMAX-1
                   this%xslopes(i,j,k) = Mesh%invdx * nolimit_limiter(&
                      rvar(i,j,k) - rvar(i-1,j,k), rvar(i+1,j,k) - rvar(i,j,k))
                END DO
             END DO
          END DO
       END IF
       ! calculate slopes in y-direction
       IF (Mesh%JNUM.GT.1) THEN
          DO k=1,Physics%VNUM
!CDIR COLLAPSE
             DO j=Mesh%JGMIN+1,Mesh%JGMAX-1
                DO i=Mesh%IGMIN,Mesh%IGMAX
                   this%yslopes(i,j,k) = Mesh%invdy * nolimit_limiter(&
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
        IF (SIGN(1.0,arg1)*SIGN(1.0,arg2).GT.0) THEN
           limarg = SIGN(MIN(ABS(arg1),ABS(arg2)),arg1)
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
        IF (((SIGN(1.0,arg1)*SIGN(1.0,arg2)).GT.0).AND.&
            ((SIGN(1.0,arg2)*SIGN(1.0,arg3)).GT.0)) THEN
           limarg = SIGN(MIN(ABS(arg1),ABS(arg2),ABS(arg3)),arg1)
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
        IF (SIGN(1.0,arg1)*SIGN(1.0,arg2).GT.0) THEN
           limarg = SIGN(MAX(MIN(param*ABS(arg1),ABS(arg2)), &
                MIN(ABS(arg1),param*ABS(arg2))),arg1)
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
             + arg2*(arg2 + 0.5*arg1) + TINY(1.0))
      END FUNCTION ospre_limiter

      ELEMENTAL SUBROUTINE pp_limiter(xslope,yslope,arg1,arg2,arg3,arg4,param,arg6,arg7,arg8,arg9) 
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        REAL, INTENT(OUT):: xslope,yslope
        REAL, INTENT(IN) :: arg1,arg2,arg3,arg4,param,arg6,arg7,arg8,arg9
        !--------------------------------------------------------------------!
        REAL             :: Vmin,Vmax,V
        !--------------------------------------------------------------------!
        xslope = (arg6-arg4)*0.5
        yslope = (arg2-arg8)*0.5
        Vmin = MIN(arg1,arg2,arg3,arg4,-param,arg6,arg7,arg8,arg9)
        Vmax = MAX(arg1,arg2,arg3,arg4,+param,arg6,arg7,arg8,arg9)
        V = 2.0*MIN(ABS(Vmin),ABS(Vmax))/(ABS(xslope)+ABS(yslope))
        V = MIN(1.0,V)
        xslope = V*xslope*Mesh%invdx
        yslope = V*yslope*Mesh%invdy
      END SUBROUTINE pp_limiter
      
      ELEMENTAL FUNCTION vanleer_limiter(arg1,arg2) RESULT(limarg)
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        REAL             :: limarg
        REAL, INTENT(IN) :: arg1, arg2
        !--------------------------------------------------------------------!
        REAL             :: a1,a2
        !--------------------------------------------------------------------!
        a1 = ABS(arg1)
        a2 = ABS(arg2)
        limarg = (arg1 * a2 + arg2 * a1)/(a1 + a2 + TINY(a1))
      END FUNCTION vanleer_limiter

      ELEMENTAL FUNCTION nolimit_limiter(arg1,arg2) RESULT(limarg)
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        REAL             :: limarg
        REAL, INTENT(IN) :: arg1, arg2
        !--------------------------------------------------------------------!
        limarg = 0.5*(arg1+arg2)
      END FUNCTION nolimit_limiter

  END SUBROUTINE CalculateSlopes_linear
  

  PURE SUBROUTINE CalculateStates_linear(this,Mesh,Physics,npos,dx,dy,rvar,rstates)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP) :: this
    TYPE(Mesh_TYP)           :: Mesh
    TYPE(Physics_TYP)        :: Physics
    INTEGER                  :: npos
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,npos) &
                             :: dx,dy
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                             :: rvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,npos,Physics%VNUM) &
                             :: rstates
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k,n
    !------------------------------------------------------------------------!
    INTENT(IN)  :: this,Mesh,Physics,npos,dx,dy,rvar
    INTENT(OUT) :: rstates
    !------------------------------------------------------------------------!
    ! reconstruct states at positions pos 
    DO k=1,Physics%VNUM
       DO n=1,npos
!CDIR COLLAPSE
          DO j=Mesh%JGMIN,Mesh%JGMAX
             DO i=Mesh%IGMIN,Mesh%IGMAX
!CDIR IEXPAND
                rstates(i,j,n,k) = reconstruct(rvar(i,j,k),this%xslopes(i,j,k), &
                     this%yslopes(i,j,k),dx(i,j,n),dy(i,j,n))
             END DO
          END DO
       END DO
    END DO
 
    CONTAINS

      ELEMENTAL FUNCTION reconstruct(cvar0,xslope0,yslope0,dx,dy) RESULT(rstate)
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        REAL :: rstate
        REAL, INTENT(IN) :: cvar0,xslope0,yslope0,dx,dy 
        !--------------------------------------------------------------------!
        rstate = cvar0 + xslope0*dx + yslope0*dy
      END FUNCTION reconstruct
      
  END SUBROUTINE CalculateStates_linear


  SUBROUTINE CloseReconstruction_linear(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Reconstruction_TYP)  :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)             :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%xslopes,this%yslopes)
  END SUBROUTINE CloseReconstruction_linear
  
END MODULE reconstruction_linear
