!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: gravity_spectral.f90                                              #
!#                                                                           #
!# Copyright (C) 2011-2014                                                   #
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
!! - parameters of \link gravity_spectral \endlink as key-values
!! \key{green,INTEGER, type of Green-function,1}
!! \key{sigma,REAL,standard deviation,0.05}
!! \key{output/potential,INTEGER,enable(=1) output of grav. potential}
!----------------------------------------------------------------------------!
!> \author Manuel Jung
!!
!! \brief poisson solver via spectral methods and direct integration
!!
!! \extends gravity_common
!! \ingroup gravity
!----------------------------------------------------------------------------!
MODULE gravity_spectral
  USE gravity_common
  USE mesh_generic
  USE physics_generic
  USE boundary_generic
  USE functions
  USE common_dict
#ifdef PARALLEL
  USE common_types, ONLY : DEFAULT_MPI_COMPLEX
#endif
#ifdef HAVE_FFTW
  USE fftw
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
  CHARACTER(LEN=32), PARAMETER &
                    :: solver_name  = "spectral"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Gravity_TYP, &
       Grid_TYP, &
       ! constants
       ! methods
       InitGravity_spectral, &
       GetAccelGravity_spectral, &
       CalcPotential_spectral, &
       CloseGravity_spectral
  !--------------------------------------------------------------------------!
  CONTAINS

  SUBROUTINE InitGravity_spectral(this,Mesh,Physics,Boundary,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP),POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Boundary_TYP), DIMENSION(4) :: Boundary
    TYPE(Dict_TYP), POINTER :: config,IO
    INTEGER           :: solver
    !------------------------------------------------------------------------!
    INTEGER           :: err, valwrite, i
    CHARACTER(LEN=32) :: info_str
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Boundary
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "gtype", solver)
    CALL InitGravity(this,solver,solver_name)
#ifndef HAVE_FFTW
    CALL Error(this,"InitGravity_spectral", &
         "Mandatory requirement fftw has not been enabled. "//&
         "Please add --with-fftw=$FFTWDIR or similar to configure call.")
#endif
#ifndef HAVE_ISO_C_BINDING
    CALL Error(this,"InitGravity_spectral", &
         "No ISO_C_BINDINGs are available for this compiler, but they are "//&
         "a Mandatory requirement for fftw.")
#endif
#ifdef HAVE_FFTW
#ifdef PARALLEL
    ! Check domains a annular rings
    IF(.NOT.(Mesh%dims(2).EQ.1)) &
      CALL Error(this,"InitGravity_spectral", &
                 "Only domains shaped as annular rings are allowed (N x 1 decompositions).")
#endif
    
    IF(.NOT.(GetType(Boundary(NORTH)).EQ.PERIODIC .AND. &
             GetType(Boundary(SOUTH)).EQ.PERIODIC)) THEN
      CALL Error(this,"InitGravity_spectral", &
                 "The boundary conditions in north and south direction have " // &
                 "to be periodic! Don't forget: This kind of " // &
                 "self-gravitation only works for polar-like coordinate " // &
                 "systems.")
    END IF
       

    CALL RequireKey(config, "green", 1) 
    CALL GetAttr(config, "green", this%green)
    
    CALL RequireKey(config, "sigma", 0.05) 
    CALL GetAttr(config, "sigma", this%sigma)
                

    IF(.NOT.(MOD(Mesh%JNUM,2)==0)) THEN
      CALL Error(this,"InitGravity_spectral", &
                 "The spectral poisson solver needs an even number of cells " // &
                 "in the phi direction due to the discrete cosinus transform.")
    END IF


! If this is the innermost domain include 1 ghost cell,
! if this is the outermost domain include also 1 ghost cell.
#ifdef PARALLEL
    IF(Mesh%IMAX.EQ.Mesh%INUM) THEN
      this%IMAX = Mesh%IMAX+1
    ELSE
      this%IMAX = Mesh%IMAX
    END IF
#else 
    !always innermost and outermost domain
    this%IMAX = Mesh%IMAX+1    
#endif
    
    this%INUM = this%IMAX - Mesh%IMIN + 1

    CALL Info(this, " POISSON--> spectral method")
    CALL Info(this, " POISSON--> Initializing")

    WRITE (info_str, '(I8)') this%green
    CALL Info(this, " POISSON--> green-fn type:     " // TRIM(info_str))
    WRITE (info_str, '(ES8.2)') this%sigma
    CALL Info(this, " POISSON--> sigma:             " // TRIM(info_str))
            
    ALLOCATE(this%row(Mesh%JNUM), &
             this%Frow(Mesh%JNUM/2+1), &
             this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
             this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM), &
#ifdef PARALLEL
             this%sbuf1(Mesh%JNUM), &
             this%rbuf1(Mesh%JNUM), &
             this%sbuf2(Mesh%JNUM), &
             this%rbuf2(Mesh%JNUM), &
             this%displ(0:GetNumProcs(this)-1),&
             this%num(0:GetNumProcs(this)-1),&
#endif
             STAT=err)

    this%p_FI = fftw_alloc_real(INT(Mesh%JNUM/2 * (this%INUM)*(Mesh%INUM), C_SIZE_T))
    CALL C_F_POINTER(this%p_FI, this%FI, &
                     [Mesh%JNUM/2,(this%INUM),Mesh%INUM])

    IF (err.NE.0) &
        CALL Error(this,"InitGravity_spectral","Memory allocation failed.")

    this%phi(:,:) = 0.
    this%accel(:,:,:) = 0.

    this%cs => GetSoundSpeeds(Physics)
    
    ! Create plans for fftw
    
    ! Use FFTW_MEASURE for calculating the fastest plan, but this
    ! costs some extra seconds 
    this%plan_r2c = fftw_plan_dft_r2c_1d(Mesh%JNUM, this%row, this%Frow, FFTW_MEASURE)
    this%plan_c2r = fftw_plan_dft_c2r_1d(Mesh%JNUM, this%Frow, this%row, FFTW_MEASURE)

    valwrite = 0
    IF (HasKey(config, "output/potential")) CALL GetAttr(config, "output/potential", valwrite)
    IF (valwrite .EQ. 1) THEN
       CALL AddField(IO, &
               "potential", &
               this%phi, &
               Dict("name" / "potential"))
    END IF

#ifdef PARALLEL
    this%displ(GetRank(this)) = Mesh%IMIN-1
    this%num(GetRank(this)) = Mesh%IMAX-Mesh%IMIN+1
    CALL MPI_AllGather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                       this%displ, 1, MPI_INTEGER, MPI_COMM_WORLD, this%error)
    CALL MPI_AllGather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                       this%num, 1, MPI_INTEGER, MPI_COMM_WORLD, this%error)
#endif

    CALL PrecomputeI_spectral(this, Mesh)

    CALL Info(this, " POISSON--> .. done initializing")
#endif

  END SUBROUTINE InitGravity_spectral

  SUBROUTINE CalcPotential_spectral(this,Mesh,Physics,pvar)
    IMPLICIT NONE

    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,m
    COMPLEX, DIMENSION(1:Mesh%INUM, &
                       Mesh%JMIN:Mesh%JMAX/2+1) &
                      :: Fdensity
    COMPLEX, DIMENSION(Mesh%IMIN:this%IMAX, &
                       Mesh%JMIN:Mesh%JMAX/2+1) &
                      :: Fphi
#ifdef PARALLEL
    COMPLEX, DIMENSION(1:Mesh%INUM) :: buf
    INTEGER           :: status(MPI_STATUS_SIZE)
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,pvar
    !------------------------------------------------------------------------!
#ifdef HAVE_FFTW
    ! Fourier transform the density with respect to Phi
    DO i=Mesh%IMIN, Mesh%IMAX
      this%row = pvar(i,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY)
      CALL fftw_execute_dft_r2c(this%plan_r2c, this%row, this%Frow)
      Fdensity(i,:) = this%Frow * Mesh%dlx(i,Mesh%JMAX)
    END DO
  
#ifdef PARALLEL
    ! Distribute Fdensity to all processes
    DO m=Mesh%JMIN,Mesh%JMAX/2+1
      buf(Mesh%IMIN:Mesh%IMAX) = Fdensity(Mesh%IMIN:Mesh%IMAX,m)
      CALL MPI_AllGatherV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,&
                          buf, this%num, this%displ, DEFAULT_MPI_COMPLEX, &
                          MPI_COMM_WORLD, this%error)
      Fdensity(:,m) = buf(:)
    END DO
#endif

    ! Integrate in radial direction by numerical quadrature
    ! Fphi_m(t,r) = int_rmin^rmax Fdensity_m(t,r') * I_m(r,r') dr'
    DO i=Mesh%IMIN, this%IMAX
      ! integrate
      DO m=Mesh%JMIN, Mesh%JMAX/2
        FPhi(i,m) = SUM(Fdensity(:,m) * this%FI(m,i-Mesh%IMIN+1,:))
      END DO

      ! Nyquist has to be zero, because FI does not include this frequency and
      ! therefore sets it to zero.
      FPhi(i,Mesh%JMAX/2+1)=0.0
        
    END DO

    ! Inverse fourier transform Fphi with respect to Phi
    DO i=Mesh%IMIN, this%IMAX
      this%Frow = Fphi(i,Mesh%JMIN:Mesh%JMAX/2+1)
      CALL fftw_execute_dft_c2r(this%plan_c2r, this%Frow, this%row)
      this%phi(i,Mesh%JMIN:Mesh%JMAX) = this%row * Physics%Constants%GN / (Mesh%JNUM)**2
    END DO

#ifdef PARALLEL
    ! send boundary data to western and receive from eastern neighbor
    IF (Mesh%neighbor(WEST).NE.MPI_PROC_NULL) &
         this%sbuf1 = this%phi(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX)
    !CALL MPI_Sendrecv_replace(&
    !  this%row,Mesh%JMAX-Mesh%JMIN+1,DEFAULT_MPI_REAL,&
    !  Mesh%neighbor(WEST),10+WEST,&
    !  Mesh%neighbor(WEST),MPI_ANY_TAG,&
    !  Mesh%comm_cart,status,this%error)
    CALL MPI_Sendrecv(this%sbuf1,Mesh%JNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(WEST),10+WEST,this%rbuf1, &
         Mesh%JNUM,DEFAULT_MPI_REAL,Mesh%neighbor(EAST), &
         MPI_ANY_TAG,Mesh%comm_cart,status,this%error)
    IF (Mesh%neighbor(EAST).NE.MPI_PROC_NULL) &
         this%phi(Mesh%IMAX+1,Mesh%JMIN:Mesh%JMAX) = this%rbuf1
    ! send boundary data to western and receive from eastern neighbor
    IF (Mesh%neighbor(EAST).NE.MPI_PROC_NULL) &
         this%sbuf2 = this%phi(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)
    !CALL MPI_Sendrecv_replace(&
    !  this%row,Mesh%JMAX-Mesh%JMIN+1,DEFAULT_MPI_REAL,&
    !  Mesh%neighbor(EAST),10+EAST,&
    !  Mesh%neighbor(WEST),MPI_ANY_TAG,&
    !  Mesh%comm_cart,status,this%error)
    CALL MPI_Sendrecv(this%sbuf2,Mesh%JNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(EAST),10+EAST,this%rbuf2, &
         Mesh%JNUM,DEFAULT_MPI_REAL,Mesh%neighbor(WEST), &
         MPI_ANY_TAG,Mesh%comm_cart,status,this%error)
    IF (Mesh%neighbor(WEST).NE.MPI_PROC_NULL) &
         this%phi(Mesh%IMIN-1,Mesh%JMIN:Mesh%JMAX) = this%rbuf2
#endif
    !print *,"rank: ", GetRank(this), "  ", this%phi(:,Mesh%JMIN)
    !CALL Error(this,"ASD","das")
    this%phi(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN-1) = this%phi(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMAX)
    this%phi(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMAX+1) = this%phi(Mesh%IMIN-1:Mesh%IMAX+1,Mesh%JMIN)

#endif
  END SUBROUTINE CalcPotential_spectral

! This is the Green's function
ELEMENTAL FUNCTION GreenFunction(r0, r1, phi, green, sigma, eps) RESULT(G)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   REAL              :: r0, r1, phi
   REAL              :: G
   REAL              :: RSQUARE
   REAL              :: eps
   INTEGER           :: green
   REAL              :: sigma
   !------------------------------------------------------------------------!
   INTENT(IN)        :: r0, r1, phi, green, sigma, eps
   !------------------------------------------------------------------------!

   
   ! This is integral, which has to be calculated for a general Z(r,z)
   !
   ! G = int_-oo^oo -Z(r1,z1) / sqrt(r**2 + r'**2 - 2*r*r'*cos(phi) +
   ! epsilon**2 + z'**2) dz1

   
   ! For a vertical Gaussian density distribution Z(r,z)
   !
   ! Z(r,z) = (2*pi*(H(r))**2)**-0.5 * exp(-z**2 / (2*(H(r))**2))
   !
   ! this integral can be evaluated analyticly, which results in:
   !
   ! G(r,r',phi) = - exp(R**2/4) * K_0(R**2/4) / (sqrt(2*pi) * H(r'))
   !
   ! with R**2 = r**2 + r'**2 - 2*r*r'*cos(phi) + epsilon**2)/(H(r'))**2
   ! and K_0 the modified Bessel function of the second kind.

   SELECT CASE(green)
     CASE(1)  ! Razor sharp disc
       RSQUARE = (r0*r0 + r1*r1 - 2.0*r0*r1*COS(phi) + eps**2)
       G = -1.0/SQRT(RSQUARE)
     CASE(2)  ! gaussian spheres
       RSQUARE = (r0*r0 + r1*r1 - 2.0*r0*r1*COS(phi) + eps**2)&
                 / (sigma*sigma)
       G = -1.0 *  Bessel_K0e(RSQUARE/4.0) &
           / SQRT(2.0*PI*sigma*sigma) 
     CASE(3)  ! Used for the orbiting cylinder example
       RSQUARE = (r0*r0 + r1*r1 - 2.0*r0*r1*COS(phi) + eps**2)
       G = LOG(RSQUARE)
     CASE DEFAULT
       ! should never happen
       G = 0.0
   END SELECT
   END FUNCTION GreenFunction


SUBROUTINE PrecomputeI_spectral(this, Mesh)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Gravity_TYP), POINTER :: this
   TYPE(Mesh_TYP)    :: Mesh
   !------------------------------------------------------------------------!
#ifdef HAVE_FFTW
   TYPE(C_PTR)        :: plan_r2r
#endif
   INTEGER            :: i0, i1, j
   REAL               :: r0, eps
   REAL, PARAMETER    :: epsmax = 4.0
   REAL, DIMENSION(Mesh%JMIN:Mesh%JMAX/2) &
                      :: phi
   REAL, DIMENSION(1:Mesh%INUM+1) :: r
#ifdef HAVE_FFTW
   REAL(C_DOUBLE), DIMENSION(:,:,:), POINTER &
                      :: I_tmp
   TYPE(C_PTR)        :: p
#endif
   INTEGER     :: rank, howmany, istride, idist, ostride, odist
   INTEGER, DIMENSION(1) &
                      :: fftkind, n
#ifdef PARALLEL
   INTEGER, DIMENSION(0:GetNumProcs(this)-1) :: num
#endif
   !------------------------------------------------------------------------!
   INTENT(IN)         :: Mesh
   !------------------------------------------------------------------------!
#ifdef HAVE_FFTW
    n(1)    = Mesh%JNUM/2
    !n(1)    = Mesh%JMAX/2-Mesh%JMIN+1       ! == Mesh%JNUM/2
    rank    = 1
    howmany = this%INUM * Mesh%INUM
    !inembed =  
    !istride = howmany
    !idist   = 1 
    istride = 1
    idist = n(1)
    ostride = istride
    odist   = idist
    fftkind(1) = FFTW_REDFT10

    p = fftw_alloc_real(INT(n(1)*howmany, C_SIZE_T))
    CALL C_F_POINTER(p, I_tmp, [Mesh%JNUM/2, this%INUM, Mesh%INUM])
    
    plan_r2r = fftw_plan_many_r2r(rank, n, howmany, &
                                  I_tmp, n, &
                                  istride, idist, &
                                  this%FI, n, &
                                  ostride, odist, &
                                  fftkind, &
                                  IOR(FFTW_MEASURE, FFTW_DESTROY_INPUT) &
                                  )
    ! Precompute the fourier transform with respect to Phi-Phi' of
    ! I(r,r',Phi-Phi') = 2 * pi * r' * G(r,r',Phi-Phi')
    ! when G is the softened Green's function.

    ! In the case of a gaussian density distribution in the z direction,
    ! we have:
    !
    ! Z(r,z) = (2*pi*H(r))**0.5 * exp( -z**2 / (2*(H(r))**2) )
    !
    ! and therefore as Green's function:
    !
    ! G(r,r',Phi-Phi') = -(exp(R**2/4)*K_0(R**2/4))/(sqrt(2*pi)*H(r'))
    !
    ! with R**2 = (r**2 + r'**2 - 2*r*r'*cos(Phi-Phi') + epsilon**2) /
    ! (H(r'))**2
    !
    ! (epsilon is a small softening parameter)

    phi = ATAN2(Mesh%bccart(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX/2,2),&
                Mesh%bccart(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX/2,1))
    WHERE(phi.LT.0.0) &
      phi = phi + 2.0*PI

    r(Mesh%IMIN:this%IMAX) &
      = SQRT(Mesh%bccart(Mesh%IMIN:this%IMAX,Mesh%JMIN,1)**2 + Mesh%bccart(Mesh%IMIN:this%IMAX,Mesh%JMIN,2)**2)
#ifdef PARALLEL
    num(GetRank(this)) = this%INUM
    CALL MPI_AllGather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                       num, 1, MPI_INTEGER, MPI_COMM_WORLD, this%error)
    CALL MPI_AllGatherV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                        r, num, this%displ, &
                        DEFAULT_MPI_REAL, MPI_COMM_WORLD, this%error)
#endif

    DO i0=Mesh%IMIN,this%IMAX
      r0 = r(i0) - 0.5 * Mesh%dlx(i0,Mesh%JMIN)

      DO i1=1,Mesh%INUM
        I_tmp(:,i0-Mesh%IMIN+1,i1) = 2.0 * PI * r(i1) &
                         * GreenFunction(r0,r(i1),phi,this%green,this%sigma,0.)
      END DO
    END DO

    CALL fftw_execute_r2r(plan_r2r, I_tmp, this%FI)

    ! Free plan_r2r
    CALL fftw_destroy_plan(plan_r2r)
    CALL fftw_free(p)
#endif
  END SUBROUTINE PrecomputeI_spectral

  FUNCTION GetAccelGravity_spectral(this,Mesh,Physics,pvar) RESULT(ac)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Gravity_TYP), POINTER :: this
   TYPE(Mesh_TYP)       :: Mesh
   TYPE(Physics_TYP)    :: Physics
   REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                        :: pvar
   !------------------------------------------------------------------------!
   INTEGER              :: i, j
   REAL, DIMENSION(:,:,:), POINTER &
                        :: ac
   !------------------------------------------------------------------------!
   INTENT(IN)           :: Mesh,Physics,pvar
   INTENT(INOUT)        :: this
   !------------------------------------------------------------------------!
   ! calc potential first
   CALL CalcPotential_spectral(this,Mesh,Physics,pvar)

   ! Maybe Physics%VNUM is greater than 2 => set all to zero
   this%accel(:,:,:) = 0.
   DO j = Mesh%JMIN,Mesh%JMAX
     DO i = Mesh%IMIN,Mesh%IMAX
       this%accel(i,j,1) = -1.0*(this%phi(i+1,j)-this%phi(i,j))/Mesh%dlx(i,j)
       this%accel(i,j,2) = -1.0*(this%phi(i+1,j+1)+this%phi(i,j+1) &
                        -this%phi(i+1,j-1)-this%phi(i,j-1))&
                         /(4.0*Mesh%dly(i,j))
     END DO
   END DO
  ac => this%accel
  END FUNCTION GetAccelGravity_spectral

  SUBROUTINE CloseGravity_spectral(this)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Gravity_TYP), POINTER :: this
   !------------------------------------------------------------------------!
#ifdef HAVE_FFTW

    ! Destroy plans
    CALL fftw_destroy_plan(this%plan_r2c)
    CALL fftw_destroy_plan(this%plan_c2r)

    ! Free memomry
    DEALLOCATE(&
#ifdef PARALLEL
               this%sbuf1, &
               this%rbuf1, &
               this%sbuf2, &
               this%rbuf2, &
               this%displ,this%num,&
#endif
               this%phi, &
               this%accel, &
               this%row, &
               this%Frow &
               )
    CALL fftw_free(this%p_FI)

#endif
    END SUBROUTINE CloseGravity_spectral

END MODULE gravity_spectral
