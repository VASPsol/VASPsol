!**********************************************************************
! Copyright 2013-2015 Kiran Mathew and Richard Hennig
!
!Licensed under the Apache License, Version 2.0 (the "License");
!you may not use this file except in compliance with the License.
!You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
!Unless required by applicable law or agreed to in writing, software
!distributed under the License is distributed on an "AS IS" BASIS,
!WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!See the License for the specific language governing permissions and
!limitations under the License.
!**********************************************************************

#include "symbol.inc"



#define sol_io IF(IU0>=0) WRITE(IU0,*)

#define debugsol
#ifdef debugsol
#define array_sum_rc(x) debug_tmp =  SUM_GRID_RC(CMPLX(ABS(x),KIND=q)); dsol_io #x,' :  ',debug_tmp
#define array_sum_rl(x) debug_tmp =  SUM_GRID_RL(ABS(x)); dsol_io #x,' :  ',debug_tmp
#define array_min_rl(x) debug_tmp = -MAX_GRID_RL(-x); dsol_io #x,' :  ',debug_tmp
#define array_max_rl(x) debug_tmp =  MAX_GRID_RL(x); dsol_io #x,' :  ',debug_tmp
#define ifdsol
#else
#define array_sum_rc(x)
#define array_sum_rl(x)
#define array_min_rl(x)
#define array_max_rl(x)
#define ifdsol !
#endif

#define ABORT      CALL TRACEBACKQQ; CALL M_exit(); CALL ABORT

#define dsol_io ifdsol sol_io









MODULE solvation_grid

 USE prec
 USE base
 USE constant
 USE mpimy
 USE poscar,  ONLY : type_info
 USE lattice, ONLY : latt
 USE mgrid,   ONLY : grid_3d
 USE IFCORE

 IMPLICIT NONE
 
 PRIVATE
 
 PUBLIC :: INIT_SOL_GRID, SETGRID, GRIDC, LATT_CUR, T_INFO
 PUBLIC :: CONV, XCORR, GRAD, DIV, LAPL
 PUBLIC :: MINIMIZE_L, L_op, L0_op, L0_inv_op
 PUBLIC :: F_soft_sphere, F_gaussian, F_exp
 PUBLIC :: D_PROD_RC, SUM_GRID_RC, SUM_GRID_RL, MAX_GRID_RL, FFT3D_RC2RL, FFT3D_RL2RC, AVG_GRID_RC
 PUBLIC :: GET_G0, SET_G0
 PUBLIC :: WRITE_TO_FILE_RC, WRITE_TO_FILE_RL
 PUBLIC :: GSQU, R, debug_tmp
 
 SAVE

 TYPE (grid_3d),   POINTER :: GRIDC
 TYPE (latt),      POINTER :: LATT_CUR
 TYPE (type_info), POINTER :: T_INFO

 REAL(q), ALLOCATABLE :: G(:,:), GSQU(:), R(:)
 
 INTEGER :: IU0, IND_G0
 
 REAL(q) :: debug_tmp
 
 
 INTERFACE GRAD
   PROCEDURE GRAD1, GRAD3
 END INTERFACE
 
 INTERFACE DIV
   PROCEDURE DIV1, DIV3
 END INTERFACE
 
 INTERFACE LAPL
   PROCEDURE LAPL_inpl, LAPL1, LAPL3
 END INTERFACE
   
 
 
CONTAINS

!>=======================SUBROUTINE INIT_SOL_GRID===========
!!
!! set pointers to VASP structures
!!
!!====================================================================
SUBROUTINE INIT_SOL_GRID(GRIDC1, LATT_CUR1, T_INFO1, IO)
 
 TYPE (grid_3d),     TARGET, INTENT(IN) :: GRIDC1
 TYPE (latt),        TARGET, INTENT(IN) :: LATT_CUR1
 TYPE (type_info),   TARGET, INTENT(IN) :: T_INFO1
 TYPE (in_struct),           INTENT(IN) :: IO

 GRIDC    => GRIDC1
 LATT_CUR => LATT_CUR1
 T_INFO   => T_INFO1
 
 IU0 = IO%IU0
 
 ALLOCATE(G(GRIDC%MPLWV,3),GSQU(GRIDC%MPLWV), R(DIMREAL(GRIDC%MPLWV)))

END SUBROUTINE INIT_SOL_GRID


!>=====================SUBROUTINE SETGRID=============================
!!
!!set up grid vectors in reciprocal space stored in G(:,:)
!!
!!====================================================================
SUBROUTINE SETGRID
 
 REAL(q) :: X1, X2, X3, X, Y, Z
 INTEGER :: NC, I, N1, N2, N3
 
 G = 0.
 I = 0
 IND_G0 = -1
 DO NC=1,GRIDC%RC%NCOL
   N2= GRIDC%RC%I2(NC)
   N3= GRIDC%RC%I3(NC)
   IF (N2==1 .AND. N3==1) THEN
      IND_G0 = I + 1
   ENDIF
   DO N1=1,GRIDC%RC%NROW
      I = I +1
      !GX,GY,GZ
      G(I,1) = (GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*&
           & LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
      G(I,2) = (GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*&
           & LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
      G(I,3) = (GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*&
           & LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))
   ENDDO
 ENDDO
 GSQU = G(:,1)**2+G(:,2)**2+G(:,3)**2
 
 
 R = 0.
 I = 0
 DO NC=1,GRIDC%RL%NCOL
   N1= GRIDC%RL%I2(NC)
   N2= GRIDC%RL%I3(NC)
   X1 = (N1-1._q)/GRIDC%NGX
   X2 = (N2-1._q)/GRIDC%NGY
   X1 = MODULO(X1+0.5_q,1._q) - 0.5_q
   X2 = MODULO(X2+0.5_q,1._q) - 0.5_q
   DO N3=1,GRIDC%RL%NROW
      X3 = (N3-1._q)/GRIDC%NGZ
      X3 = MODULO(X3+0.5_q,1._q) - 0.5_q
      I = I +1
      !X,Y,Z
      X = (X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3))
      Y = (X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3))
      Z = (X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3))
      R(I) = SQRT(X**2 + Y**2 + Z**2)
   ENDDO
 ENDDO
 
END SUBROUTINE SETGRID


!>=====================SUBROUTINE MINIMIZE_L============================
!!
!!routine for solving the LPB equation using conjugate gradient
!!
!!======================================================================
SUBROUTINE MINIMIZE_L(resid, chi, ekappa2, w, dphi, tol, A_LPB, NCGSOL)
 
 COMPLEX(q), INTENT(IN)    :: resid(:)
 REAL(q),    INTENT(IN)    :: chi(:,:,:) 
 REAL(q),    INTENT(IN)    :: ekappa2(:)
 COMPLEX(q), INTENT(IN)    :: w(:)
 COMPLEX(q), INTENT(OUT)   :: dphi(:)
 REAL(q),    INTENT(IN)    :: tol
 REAL(q),    INTENT(INOUT) :: A_LPB
 INTEGER,    INTENT(OUT)   :: NCGSOL
 
 !local variables
 COMPLEX(q), ALLOCATABLE :: r(:), z(:), p(:), Lp(:), Lp0(:)
 REAL(q),    ALLOCATABLE :: PRECOND(:)
 
 REAL(q) :: A_0, alpha,beta,gamma,lambda0,lambda,rMr,rMr_old,pLp,rms,r0,p0
 INTEGER :: iter
 LOGICAL :: LION
 
 ALLOCATE(r(GRIDC%MPLWV),z(GRIDC%MPLWV),p(GRIDC%MPLWV),Lp(GRIDC%MPLWV),Lp0(GRIDC%MPLWV), PRECOND(GRIDC%MPLWV))
 
 !*********************************************************
 ! start of Pre-conditioned Conjugate Gradient
 !**********************************************************
 
 !dsol_io '  Starting LPB CG loop'
 
 ! construct preconditioner
 WHERE (GSQU > 0._q)
   PRECOND = 1._q/(TPI**2*GSQU) * EDEPS/LATT_CUR%OMEGA
 ELSEWHERE
   PRECOND = 0._q
 ENDWHERE
 
 ! initialize dphi, residual, and energy
 dphi = 0._q
 r = resid
 A_0 = A_LPB
 r0 = 0._q
 
 z = PRECOND * r
 rms = sqrt(D_PROD_RC(z,z))
 
 !ifdsol IF(IU0>=0) WRITE(IU0,'(3X,5X,A4,2A24,3A15)') 'iter', 'A_solv', 'rms', 'alpha', 'beta', 'gamma'
 ifdsol IF(IU0>=0) WRITE(IU0,'(3X,A5,4X,2E24.15)') 'LPB: ', A_LPB, rms
 
 ! set the G=0 component of dphi to zero out the G=0 component of residual
 LION = SIZE(ekappa2) > 0
 IF (LION) THEN
   CALL GET_G0(r, r0)
   CALL FFT3D_RL2RC(ekappa2, Lp0)
   Lp0 = Lp0 * LATT_CUR%OMEGA/EDEPS
   CALL GET_G0(Lp0, lambda0)
   alpha = r0/lambda0
   CALL SET_G0(dphi, alpha)
   r = r - alpha*Lp0
   z = PRECOND * r
   rms = sqrt(D_PROD_RC(z,z))
   A_LPB = A_0 + 0.5_q*alpha*r0
   ifdsol IF(IU0>=0) WRITE(IU0,'(3X,A5,I4,2E24.15,4E15.6)') 'LPB: ', 0, A_LPB, rms, alpha, 0._q, lambda0, r0
 ENDIF
 
 
 !-------------------------------------------------------------
 !CG loop
 !-------------------------------------------------------------
 iter = 1
 beta = 0._q
 DO
   
   rMr = D_PROD_RC(r,z)
   
   ! CG search direction
   IF (iter > 1) THEN
      beta = rMr/rMr_old
      p = z + beta*p
   ELSE
      p = z;
   ENDIF
   
   ! set the G=0 component of p to zero out the G=0 component of Lp
   IF (LION) THEN
      lambda = D_PROD_RC(z,Lp0)
      CALL GET_G0(p, p0)
      p0 = p0 - lambda/lambda0
      CALL SET_G0(p, p0)
   ENDIF
   CALL GET_G0(r, r0)

   ! compute Lp = L(p)
   CALL L_op(p, chi, ekappa2, w, Lp)
   
   ! compute step size
   ! alpha = r.M.r/p.L.p
   pLp = D_PROD_RC(p,Lp)
   alpha = rMr/pLp
   gamma = pLp/D_PROD_RC(p,p)
   
   ! take step along search direction and compute new residual
   dphi = dphi + alpha * p
   r = r - alpha*Lp
   ! new preconditioned residual
   z = PRECOND * r
   
   ! rms of the residual
   rms = sqrt(D_PROD_RC(z,z))
   ! energy
   A_LPB = A_0 + 0.5_q*D_PROD_RC(dphi,resid+r)
   
   ifdsol IF(IU0>=0) WRITE(IU0,'(3X,A5,I4,2E24.15,4E15.6)') 'LPB: ', iter, A_LPB, rms, alpha, beta, gamma, r0
   
   IF (rms <= tol .AND. iter >= 4) EXIT
   
   iter = iter + 1
   IF (iter > 200) THEN
      IF (IU0>=0) WRITE(IU0,*) 'ERROR: MINIMIZE_L: failed to converge in 200 iterations'
      STOP
   ENDIF
   
   rMr_old = rMr
  
 ENDDO
 !-------------------------------------------------------------
 !end of CG loop
 !-------------------------------------------------------------
 NCGSOL = iter
 
 !dsol_io '  LPB CG loop finished'
 
END SUBROUTINE MINIMIZE_L


!>======================SUBROUTINE L_op===============================
!!
!!subroutine that evaluates the fourier transform of 
!!Lp = L(dphi) = -(div(epsilon(grad(dphi))) - ekappa2*dphi)/EDEPS
!!epsilon is the nonlocal dielectric tensor operator defined by
!!  epsilon(v) = v + conv(w, chi*xcorr(w,v))
!!Lp is scaled by OMEGA to have units of charge
!!
!!====================================================================
SUBROUTINE L_op(dphi, chi, ekappa2, w, Lp)
 
 COMPLEX(q), INTENT(IN)  :: dphi(:)
 REAL(q),    INTENT(IN)  :: chi(:,:,:), ekappa2(:)
 COMPLEX(q), INTENT(IN)  :: w(:)
 COMPLEX(q), INTENT(OUT) :: Lp(:)
 
 COMPLEX(q) :: CWORK(GRIDC%MPLWV)
 REAL(q)    :: RWORK(DIMREAL(GRIDC%MPLWV))
 
 ! compute conv(w, div(chi.xcorr(w,grad(dphi)))) + lapl(dphi)
 IF (SIZE(w)>0) THEN
   CWORK = conjg(w)*dphi
   CALL LAPL(CWORK,Lp,chi)
   Lp = w*Lp
   CALL LAPL(dphi,CWORK)
   Lp = Lp + CWORK
 ELSE
   CALL LAPL(dphi,Lp,1._q+chi)
 ENDIF
 
 ! add -ekappa2*dphi
 IF (SIZE(ekappa2)>0) THEN
   CALL FFT3D_RC2RL(dphi,RWORK)
   RWORK = -ekappa2*RWORK
   CALL FFT3D_RL2RC(RWORK,CWORK)
   Lp = Lp + CWORK
 ENDIF
 
 ! scale Lp
 Lp = -LATT_CUR%OMEGA/EDEPS * Lp

END SUBROUTINE L_op


!>==========================SUBROUTINE CONV===========================
!!
!!
!!====================================================================
SUBROUTINE CONV(phi, w)

 REAL(q),    INTENT(INOUT) :: phi(:)
 COMPLEX(q), INTENT(IN)    :: w(:)
 
 COMPLEX(q) :: CWORK(GRIDC%MPLWV)
 
 CALL FFT3D_RL2RC(phi,CWORK)
 CWORK = w*CWORK
 CALL FFT3D_RC2RL(CWORK,phi)

END SUBROUTINE CONV


!>==========================SUBROUTINE XCORR==========================
!!
!!
!!====================================================================
SUBROUTINE XCORR(phi, w)

 REAL(q),    INTENT(INOUT) :: phi(:)
 COMPLEX(q), INTENT(IN)    :: w(:)
 
 COMPLEX(q) :: CWORK(GRIDC%MPLWV)
 
 CALL FFT3D_RL2RC(phi,CWORK)
 CWORK = CONJG(w)*CWORK
 CALL FFT3D_RC2RL(CWORK,phi)

END SUBROUTINE XCORR


!>==========================SUBROUTINE GRAD===========================
!!
!! phi, input scalar field in reciprocal space
!! gradphi(:,1:3) contains grad(phi) in real space
!! gradphi(:,4) contains |grad(phi)| in real space
!!
!!====================================================================
SUBROUTINE GRAD3(phi, gradphi)

 COMPLEX(q), INTENT(IN)  :: phi(:)
 REAL(q),    INTENT(OUT) :: gradphi(:,:)
 
!local variables
 COMPLEX(q) :: CWORK(GRIDC%MPLWV)
 INTEGER    :: I
 
 gradphi = 0._q
 DO I = 1,3
   CWORK = CITPI*G(:,I)*phi
   CALL FFT3D_RC2RL(CWORK,gradphi(:,I))
 ENDDO
 IF (SIZE(gradphi,2) >= 4) gradphi(:,4) = SQRT(SUM(gradphi(:,1:3)**2,2))

END SUBROUTINE GRAD3


!>==========================SUBROUTINE GRAD===========================
!!
!! phi, input scalar field in reciprocal space
!! gradphi(:,1:3) contains grad(phi) in real space
!! gradphi(:,4) contains |grad(phi)| in real space
!!
!!====================================================================
SUBROUTINE GRAD1(phi, gradphi)

 COMPLEX(q), INTENT(IN)  :: phi(:)
 REAL(q),    INTENT(OUT) :: gradphi(:)
 
!local variables
 REAL(q)    :: RWORK(DIMREAL(GRIDC%MPLWV))
 COMPLEX(q) :: CWORK(GRIDC%MPLWV)
 INTEGER    :: I
 
 gradphi = 0._q
 DO I = 1,3
   CWORK = CITPI*G(:,I)*phi
   CALL FFT3D_RC2RL(CWORK,RWORK)
   gradphi = gradphi + RWORK**2
 ENDDO
 gradphi = SQRT(gradphi)

END SUBROUTINE GRAD1


!>==========================SUBROUTINE DIV3===========================
!!
!! calculate the divergence of a vector field in real space
!! A, the input in real space
!! chi, an optional scalar or tensor field in real space
!!  that multiplies the vector field
!! divA, the output in rec space
!!
!!====================================================================
SUBROUTINE DIV3(A,divA,chi)

 REAL(q),    INTENT(IN)  :: A(:,:)
 COMPLEX(q), INTENT(OUT) :: divA(:)
 REAL(q),    INTENT(IN)  :: chi(:,:,:)
 
 !local variables
 COMPLEX(q) :: CWORK(GRIDC%MPLWV)
 REAL(q)    :: RWORK(DIMREAL(GRIDC%MPLWV))
 INTEGER    :: I, J
 
 IF (SIZE(chi,2) == 1 .AND. SIZE(chi,3) == 1) THEN
   CALL DIV1(A,divA,chi(:,1,1))
   RETURN
 ENDIF

 divA = 0.
 DO I = 1,3
   RWORK = 0.
   DO J = 1,3
      RWORK = RWORK + chi(:,I,J)*A(:,J)
   ENDDO
   CALL FFT3D_RL2RC(RWORK,CWORK)
   divA = divA + CITPI*G(:,I)*CWORK
 ENDDO

END SUBROUTINE DIV3


! allow DIV to be called with a scalar chi
SUBROUTINE DIV1(A,divA,chi)

 REAL(q),    INTENT(IN)  :: A(:,:)
 COMPLEX(q), INTENT(OUT) :: divA(:)
 REAL(q),    INTENT(IN)  :: chi(:)
 
 OPTIONAL :: chi
 
 COMPLEX(q) :: CWORK(GRIDC%MPLWV)
 REAL(q)    :: RWORK(DIMREAL(GRIDC%MPLWV))
 INTEGER    :: I
 
 divA = 0.
 DO I = 1,3
   IF (PRESENT(chi)) THEN
      RWORK = chi*A(:,I)
      CALL FFT3D_RL2RC(RWORK,CWORK)
   ELSE
      CALL FFT3D_RL2RC(A(:,I),CWORK)
   ENDIF
   divA = divA + CITPI*G(:,I)*CWORK
 ENDDO
 
END SUBROUTINE DIV1


!>======================SUBROUTINE LAPL==========================
!!
!! on input, phi contains scalar field in rec space
!! takes optional scalar or tensor metric chi in real space
!! on output, phi contains lapl(phi) in rec space
!!
!!===============================================================
SUBROUTINE LAPL3(phi,lapl_phi,chi)

 COMPLEX(q), INTENT(IN)  :: phi(:)
 COMPLEX(q), INTENT(OUT) :: lapl_phi(:)
 REAL(q),    INTENT(IN)  :: chi(:,:,:)
 
 COMPLEX(q) :: CWORK(GRIDC%MPLWV)
 REAL(q)    :: RWORK(DIMREAL(GRIDC%MPLWV))
 REAL(q)    :: RWORKV(DIMREAL(GRIDC%MPLWV),3)
 INTEGER    :: I, J
 
 IF (SIZE(chi,2) == 1 .AND. SIZE(chi,3) == 1) THEN
   CALL LAPL1(phi,lapl_phi,chi(:,1,1))
   RETURN
 ENDIF
 
 DO I = 1,3
   CWORK = CITPI*G(:,I)*phi
   CALL FFT3D_RC2RL(CWORK,RWORKV(:,I))
 ENDDO
 
 lapl_phi = 0._q
 DO I = 1,3
   RWORK = 0._q
   DO J = 1,3
      RWORK = RWORK + chi(:,I,J)*RWORKV(:,J)
   ENDDO
   CALL FFT3D_RL2RC(RWORK,CWORK)
   lapl_phi = lapl_phi + CITPI*G(:,I)*CWORK
 ENDDO
 
END SUBROUTINE LAPL3


SUBROUTINE LAPL1(phi,lapl_phi,chi)

 COMPLEX(q), INTENT(IN)  :: phi(:)
 COMPLEX(q), INTENT(OUT) :: lapl_phi(:)
 REAL(q),    INTENT(IN)  :: chi(:)
 
 OPTIONAL :: chi
 
 COMPLEX(q) :: CWORK(GRIDC%MPLWV)
 REAL(q)    :: RWORK(DIMREAL(GRIDC%MPLWV))
 INTEGER    :: I
 
 IF(PRESENT(chi)) THEN
   lapl_phi = 0._q
   DO I = 1,3
      CWORK = CITPI*G(:,I)*phi
      CALL FFT3D_RC2RL(CWORK,RWORK)
      RWORK = RWORK*chi
      CALL FFT3D_RL2RC(RWORK,CWORK)
      lapl_phi = lapl_phi + CITPI*G(:,I)*CWORK
   ENDDO
 ELSE
   lapl_phi = -TPI**2*GSQU * phi
 ENDIF
 
END SUBROUTINE LAPL1


SUBROUTINE LAPL_inpl(phi)
 COMPLEX(q), INTENT(INOUT)  :: phi(:)
 phi = -TPI**2*GSQU * phi
END SUBROUTINE LAPL_inpl


!>======================SUBROUTINE L0_op=========================
!!
!! on input, phi contains phi in rec space
!! on output, L0_phi contains L0(phi) in rec space
!!
!!===============================================================
SUBROUTINE L0_op(phi,L0_phi)

 COMPLEX(q), INTENT(IN)  :: phi(:)
 COMPLEX(q), INTENT(OUT) :: L0_phi(:)
 
 L0_phi = phi * (TPI**2*GSQU) * LATT_CUR%OMEGA/EDEPS
 
END SUBROUTINE L0_op


!>======================SUBROUTINE L0_inv_op=====================
!!
!! on input, L0_phi contains L0(phi) in rec space
!! on output, phi contains phi in rec space
!!
!!===============================================================
SUBROUTINE L0_inv_op(L0_phi,phi)

 COMPLEX(q), INTENT(IN)  :: L0_phi(:)
 COMPLEX(q), INTENT(OUT) :: phi(:)
 
 WHERE (GSQU > 0._q)
   phi = L0_phi/(TPI**2*GSQU) * EDEPS/LATT_CUR%OMEGA
 ELSEWHERE
   phi = 0._q
 ENDWHERE
 
END SUBROUTINE L0_inv_op


!>======================SUBROUTINE F_soft_sphere=================
!!
!! construct real space soft sphere function on rec space grid
!!
!!===============================================================
SUBROUTINE F_soft_sphere(F, R0, sigma)

 COMPLEX(q), INTENT(OUT) :: F(:)
 REAL(q),    INTENT(IN)  :: R0, sigma
 
 REAL(q) :: RWORK(DIMREAL(GRIDC%MPLWV))
 
 RWORK = (R - R0)/(sqrt(2._q)*sigma)
 RWORK = 0.5_q - 0.25_q*RWORK
 WHERE (RWORK <= 0._q)
   RWORK = 0._q
 ELSEWHERE (1._q <= RWORK)
   RWORK = 1._q
 ELSEWHERE
   RWORK = RWORK - sin(TPI*RWORK)/TPI
 ENDWHERE
 
 RWORK = RWORK*GRIDC%NPLWV/SUM_GRID_RL(RWORK)
 CALL FFT3D_RL2RC(RWORK,F)
 
END SUBROUTINE F_soft_sphere


!>======================SUBROUTINE F_gaussian====================
!!
!! construct real space gaussian function on rec space grid
!!
!!===============================================================
SUBROUTINE F_gaussian(F, sigma)

 COMPLEX(q), INTENT(OUT) :: F(:)
 REAL(q),    INTENT(IN)  :: sigma
 
 REAL(q) :: RWORK(DIMREAL(GRIDC%MPLWV))
 
 RWORK = EXP(-0.5_q*(R/sigma)**2) / (sigma*SQRT(TPI))**3
 RWORK = RWORK*GRIDC%NPLWV/SUM_GRID_RL(RWORK)
 CALL FFT3D_RL2RC(RWORK,F)
 
END SUBROUTINE F_gaussian


!>======================SUBROUTINE F_exp=========================
!!
!! construct real space exp function on rec space grid
!!
!!===============================================================
SUBROUTINE F_exp(F, R_c, sigma)

 COMPLEX(q), INTENT(OUT) :: F(:)
 REAL(q),    INTENT(IN)  :: R_c, sigma
 
 REAL(q), PARAMETER :: CUTOFF=100._q
 REAL(q) :: RWORK(DIMREAL(GRIDC%MPLWV))
 
 RWORK = 1._q/(1._q/CUTOFF + EXP((R-R_c)/sigma)) / (4._q*TPI*sigma**3) * 4._q/(2._q+R_c/sigma)
 RWORK = RWORK * GRIDC%NPLWV/LATT_CUR%OMEGA
 CALL FFT3D_RL2RC(RWORK,F)
 
END SUBROUTINE F_exp


!>==========================SUBROUTINE D_PROD_RC====================
!!
!! helper routine for dot product over scalar fields in rec space
!! with parallel storage
!! assumes the fields are real-valued in real space
!!
!!==================================================================
REAL(q) FUNCTION D_PROD_RC(ARRAY1,ARRAY2)

 COMPLEX(q),    INTENT(IN)  :: ARRAY1(:)
 COMPLEX(q),    INTENT(IN)  :: ARRAY2(:)
 
 D_PROD_RC = SUM_GRID_RC(ARRAY1*CONJG(ARRAY2))

END FUNCTION D_PROD_RC


!>==========================SUBROUTINE SUM_GRID_RC==================
!!
!! helper routine for sum over scalar fields in real space
!! with parallel storage
!! assumes the fields are real-valued in real space
!!
!!==================================================================
REAL(q) FUNCTION SUM_GRID_RC(ARRAY)

 COMPLEX(q), INTENT(IN)  :: ARRAY(:)
 
 REAL(q) :: sum_val
 
 REAL(q) :: FACTM
 INTEGER :: I, NC, N1, N2, N3, N1UNB, N2UNB, N3UNB
 
 ! don't include unbalanced G-vectors
 N1UNB = GRIDC%NGX/2+1
 N2UNB = GRIDC%NGY/2+1
 N3UNB = GRIDC%NGZ/2+1
 
 sum_val = 0.
 I = 0
 col: DO NC=1,GRIDC%RC%NCOL
   N2= GRIDC%RC%I2(NC)
   N3= GRIDC%RC%I3(NC)
   row: DO N1=1,GRIDC%RC%NROW
      I=I+1
      !IF (N1==N1UNB .OR. N2==N2UNB .OR. N3==N3UNB) CYCLE
#ifdef NGZhalf
      IF (N3==1 .OR. N3==N3UNB) THEN
         FACTM = 1._q
      ELSE
         FACTM = 2._q
      ENDIF
#elif defined(NGXhalf)
      IF (N1==1 .OR. N1==N1UNB) THEN
         FACTM = 1._q
      ELSE
         FACTM = 2._q
      ENDIF
#else
      FACTM = 1._q
#endif
      sum_val = sum_val + FACTM * REAL(ARRAY(I),KIND=q)
   ENDDO row
 ENDDO col
 CALLMPI( M_sum_d(GRIDC%COMM,sum_val,1))
 SUM_GRID_RC = sum_val

END FUNCTION SUM_GRID_RC


!>==========================SUBROUTINE SUM_GRID_RL==================
!!
!! helper routine for sum over scalar fields in real space
!! with parallel storage
!!
!!==================================================================
REAL(q) FUNCTION SUM_GRID_RL(ARRAY)

 REAL(q), INTENT(IN) :: ARRAY(:)
 
 REAL(q) :: sum_val
 
 sum_val = sum(ARRAY(1:GRIDC%RL%NP))
 CALLMPI( M_sum_d(GRIDC%COMM,sum_val,1))
 SUM_GRID_RL = sum_val

END FUNCTION SUM_GRID_RL


!>==========================SUBROUTINE MAX_GRID_RL==================
!!
!! helper routine for sum over scalar fields in real space
!! with parallel storage
!!
!!==================================================================
REAL(q) FUNCTION MAX_GRID_RL(ARRAY)

 REAL(q), INTENT(IN) :: ARRAY(:)
 
 REAL(q) :: max_val
 
 !max_val = maxval(ARRAY(1:GRIDC%RL%NP))
 max_val = maxval(ARRAY)
 CALLMPI( M_max_d(GRIDC%COMM,max_val,1))
 MAX_GRID_RL = max_val

END FUNCTION MAX_GRID_RL


!>=====================SUBROUTINE SETGRID=============================
!!
!!set up grid vectors in reciprocal space stored in G(:,:)
!!
!!====================================================================
SUBROUTINE AVG_GRID_RC(ARRAY, ARRAY_AVG, IDIR)
 
 COMPLEX(q), INTENT(IN)  :: ARRAY(:)
 COMPLEX(q), INTENT(OUT) :: ARRAY_AVG(:)
 INTEGER,    INTENT(IN)  :: IDIR
 
 INTEGER :: NC, N1, N2, N3, I
 
 ARRAY_AVG = 0._q
 
 I = 0
 SELECT CASE(IDIR)
   CASE(1)
      DO NC=1,GRIDC%RC%NCOL
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)
         IF (N2==1 .AND. N3==1) ARRAY_AVG(I+1:I+GRIDC%RC%NROW) = ARRAY(I+1:I+GRIDC%RC%NROW)
         I = I + GRIDC%RC%NROW
      ENDDO
   CASE(2)
      DO NC=1,GRIDC%RC%NCOL
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)
         IF (N3==1) ARRAY_AVG(I+1) = ARRAY(I+1)
         I = I + GRIDC%RC%NROW
      ENDDO
   CASE(3)
      DO NC=1,GRIDC%RC%NCOL
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)
         IF (N2==1) ARRAY_AVG(I+1) = ARRAY(I+1)
         I = I + GRIDC%RC%NROW
      ENDDO
   CASE DEFAULT
      IF (IU0>=0) WRITE(IU0,*) 'ERROR: AVG_GRID_RC: invalid direction'
      STOP
   END SELECT
   
END SUBROUTINE AVG_GRID_RC
 


!>==========================SUBROUTINE FFT3D_RC2RL==================
!!
!! helper routine to FFT a scalar field from rec to real space
!! assumes the field is real-valued in real space
!!
!!==================================================================
SUBROUTINE FFT3D_RC2RL(GRIDQ_RC,GRIDQ_RL)

 COMPLEX(q), INTENT(IN)  :: GRIDQ_RC(:)
 REAL(q),    INTENT(OUT) :: GRIDQ_RL(:)
 
 COMPLEX(q) :: CWORK(GRIDC%MPLWV)
 
 ! fft to real space
 CWORK = GRIDQ_RC
 !CALL SETUNB(CWORK,GRIDC)
 CALL FFT3D(CWORK,GRIDC,1)
 ! unpack complex array into real array
 CALL UNPACK_C2R(CWORK,GRIDQ_RL)

END SUBROUTINE FFT3D_RC2RL


!>==========================SUBROUTINE FFT3D_RL2RC==================
!!
!! helper routine to FFT a scalar field from real to rec space
!! assumes the field is real-valued in real space
!!
!!==================================================================
SUBROUTINE FFT3D_RL2RC(GRIDQ_RL,GRIDQ_RC)

 REAL(q),       INTENT(IN)  :: GRIDQ_RL(:)
 COMPLEX(q),    INTENT(OUT) :: GRIDQ_RC(:)

 CALL PACK_R2C(GRIDQ_RL,GRIDQ_RC)
 CALL FFT3D(GRIDQ_RC,GRIDC,-1)
 !CALL SETUNB(GRIDQ_RC,GRIDC)
 GRIDQ_RC = GRIDQ_RC/GRIDC%NPLWV

END SUBROUTINE FFT3D_RL2RC


!>==========================SUBROUTINE PACK_R2C======================
!!
!!  helper routine to pack real array(RA) to a complex array(CB)
!!
!!===================================================================
SUBROUTINE PACK_R2C(RA,CB)

 REAL(q),       INTENT(IN)  :: RA(:)
 COMPLEX(q),    INTENT(OUT) :: CB(:)
 
 INTEGER(q) :: I,J

 CB = 0.
 
 J = 1
 DO I=1,GRIDC%MPLWV
#ifdef realmode
    CB(I) = CMPLX(RA(J),RA(J+1),KIND=q)
    J = J + 2
#else
    CB(I) = CMPLX(RA(I),KIND=q)
#endif
 ENDDO

END SUBROUTINE PACK_R2C


!>==========================SUBROUTINE UNPACK_C2R====================
!!
!!  helper routine to unpack complex array(CA) to real array(RB)
!!
!!===================================================================
SUBROUTINE UNPACK_C2R(CA,RB)

 COMPLEX(q),     INTENT(IN)  :: CA(:)
 REAL(q),        INTENT(OUT) :: RB(:)
 
 INTEGER(q) :: I ,J
 
 RB = 0.
 
 J = 1
 DO I=1,GRIDC%MPLWV
#ifdef realmode
    RB(J) = REAL(CA(I),KIND=q)
    RB(J+1) = AIMAG(CA(I))
    J = J + 2
#else
    RB(I) = REAL(CA(I),KIND=q)
#endif
 ENDDO

END SUBROUTINE UNPACK_C2R


!>==========================SUBROUTINE GET_G0========================
!!
!!  retrieve G=0 component of array in rec space
!!
!!===================================================================
SUBROUTINE GET_G0(C, C0)
 
 COMPLEX(q), INTENT(IN) ::  C(:)
 REAL(q),    INTENT(OUT) :: C0
 
 C0 = 0._q
 IF (IND_G0 > 0) C0 = C(IND_G0)
 
 CALLMPI( M_sum_d( GRIDC%COMM, C0, 1))
 
END SUBROUTINE GET_G0


!>==========================SUBROUTINE GET_G0========================
!!
!!  set G=0 component of array in rec space
!!
!!===================================================================

SUBROUTINE SET_G0(C, C0)
 
 COMPLEX(q), INTENT(INOUT) ::  C(:)
 REAL(q),    INTENT(IN)    :: C0
 
 IF (IND_G0 > 0) C(IND_G0) = C0
 
END SUBROUTINE SET_G0


!>==========================SUBROUTINE WRITE_TO_FILE_RC==================
!!
!! write rec space quantities to file
!! writes poscar first and then A
!! could be any quantity defined on RC grid
!!
!!====================================================================
SUBROUTINE WRITE_TO_FILE_RC(FNAME, A)

 USE poscar, ONLY : OUTPOS
 USE fileio, ONLY : OUTCHG
 
 COMPLEX(q),      INTENT(IN) :: A(:)
 CHARACTER(*),    INTENT(IN) :: FNAME

 CHARACTER(40) :: SZNAM
 INTEGER :: IU
 
 IU = 468
 SZNAM = 'VASPsol'

 IF (IU0>0) THEN
     OPEN(IU, FILE = FNAME, STATUS = 'UNKNOWN')
     REWIND IU
     CALL OUTPOS(IU, .FALSE., SZNAM, T_INFO, LATT_CUR%SCALE, LATT_CUR%A, .FALSE., T_INFO%POSION)
 ENDIF
 
 CALL OUTCHG(GRIDC, IU, .TRUE., A)

 IF (IU0>0) CLOSE(IU)

END SUBROUTINE WRITE_TO_FILE_RC


!>==========================SUBROUTINE WRITE_TO_FILE_RL==================
!!
!! write real space quantities to file
!! writes poscar first and then A
!! could be any quantity defined on RL grid
!!
!!====================================================================
SUBROUTINE WRITE_TO_FILE_RL(FNAME, A)

 USE poscar, ONLY : OUTPOS
 USE fileio, ONLY : OUTPOT
 
 REAL(q),         INTENT(IN) :: A(:)
 CHARACTER(*),    INTENT(IN) :: FNAME
 
 REAL(q)    :: RWORK(DIMREAL(GRIDC%MPLWV))
 COMPLEX(q) :: CWORK(GRIDC%MPLWV)
 CHARACTER(40) :: SZNAM
 INTEGER :: IU
 
 IU = 468
 SZNAM = 'VASPsol'

 IF (IU0>0) THEN
     OPEN(IU, FILE = FNAME, STATUS = 'UNKNOWN')
     REWIND IU
     CALL OUTPOS(IU, .FALSE., SZNAM, T_INFO, LATT_CUR%SCALE, LATT_CUR%A, .FALSE., T_INFO%POSION)
 ENDIF
 
 RWORK = 0._q
 WHERE (ABS(A) > 1e-50_q) RWORK = A
 CALL PACK_R2C(RWORK,CWORK)
 CALL OUTPOT(GRIDC, IU, .TRUE., CWORK)

 IF (IU0>0) CLOSE(IU)

END SUBROUTINE WRITE_TO_FILE_RL




END MODULE solvation_grid









MODULE solvation_lpcm
 
 USE prec
 USE base
 USE constant 
 USE mpimy 
 USE IFCORE
 USE solvation_grid
 
 IMPLICIT NONE
 
 PRIVATE

 PUBLIC :: LPCM_solver, SET_SOL_WRITE_LPCM
 PUBLIC :: LPCM_READER, LPCM_INIT, LPCM_WRITER, XML_WRITE_LPCM
 
 SAVE
 
 INTEGER :: IU0, IU6
 
 COMPLEX(q), ALLOCATABLE :: phi(:), phi_solv(:)
 REAL(q),    ALLOCATABLE :: S(:), Sprime(:)
 
 REAL(q) :: eps = 1E-10
 
 REAL(q) :: sigma_k = 0.6_q
 REAL(q) :: nc_k = 0.0025_q
 REAL(q) :: tau = 5.25E-4_q
 REAL(q) :: eb_k = 78.4_q
 REAL(q) :: lambda_d_k = 0._q
 
 REAL(q) :: ekappa2_b
 LOGICAL :: LION
 LOGICAL :: LWRITE=.FALSE.

CONTAINS


!>=====================SUBROUTINE LPCM_solver=========================
!!
!! solve the LPB or NLPB equations for phi
!! also compute the correction to the KS potential and free energy
!!
!!====================================================================
SUBROUTINE LPCM_solver(phi_sol, n_e, q_sol, tol, V_corr, n_solv, A_corr)

 USE ini
 
 COMPLEX(q), INTENT(IN)  :: phi_sol(:)
 COMPLEX(q), INTENT(IN)  :: n_e(:)
 REAL(q),    INTENT(IN)  :: q_sol, tol
 COMPLEX(q), INTENT(OUT) :: V_corr(:)
 COMPLEX(q), INTENT(OUT) :: n_solv(:)
 REAL(q),    INTENT(OUT) :: A_corr
 
 COMPLEX(q) :: dphi(:), resid(:), n_sol(:), n_b(:), n_ion(:)
 REAL(q)    :: ekappa2(:), chi(:,:,:)
 REAL(q)    :: A_cav, A_solv, A_sol
 INTEGER    :: NCGSOL
 
 ALLOCATABLE :: dphi, resid, n_sol, n_b, n_ion, chi, ekappa2
 
 ALLOCATE(n_sol(GRIDC%MPLWV), n_b(GRIDC%MPLWV), n_ion(GRIDC%MPLWV))
 ALLOCATE(dphi(GRIDC%MPLWV), resid(GRIDC%MPLWV))
 ALLOCATE(chi(DIMREAL(GRIDC%MPLWV),1,1))
 IF (LION) THEN
   ALLOCATE(ekappa2(DIMREAL(GRIDC%MPLWV)))
 ELSE
   ALLOCATE(ekappa2(0))
 ENDIF
 
 ! create solvent cavity S and compute cavity formation potential and energy
 dsol_io 'Creating dielectric and ionic cavities'
 CALL CREATECAVITY(n_e, V_corr, A_cav)
 CALL STOP_TIMING("SOL",IU6,"CAVITY")
 
 ! construct chi and ekappa2
 chi(:,1,1) = (eb_k-1._q)*S
 IF (LION) ekappa2 = ekappa2_b*S
 
 ! compute initial residual and solvation free energy
 phi = phi_sol + phi_solv
 CALL L0_op(phi_sol, n_sol)
 CALL SET_G0(n_sol, q_sol)
 CALL L_op(phi, chi, ekappa2, [COMPLEX(q)::], resid)
 resid = n_sol - resid
 A_solv = 0.5_q*D_PROD_RC(phi_solv, n_sol) + 0.5_q*D_PROD_RC(phi, resid)
 
 ! solve LPB equations for phi
 dsol_io ' Solving linear dielectric/ionic screening model'
 ifdsol IF(IU0>=0) WRITE(IU0,'(2X,6X,A4,2A24)') 'iter', 'A_solv', 'rms'
 CALL MINIMIZE_L(resid, chi, ekappa2, [COMPLEX(q)::], dphi, tol, A_solv, NCGSOL)
 phi = phi + dphi
 phi_solv = phi_solv + dphi
 dsol_io ' Finished solving LPB equations'
 CALL STOP_TIMING("SOL",IU6,"MIN_L")
 
 !Computes V_solv and adds to V_corr
 !also computes n_b and n_ion
 CALL SOL_pot(V_corr, n_b, n_ion)
 CALL STOP_TIMING("SOL",IU6,"SOL_pot")
 
 IF (LWRITE) CALL WRITE_TO_FILE_RC('VSOLV', V_corr)
 
 !add electrostatic potential correction (phi_solv) to V_corr
 V_corr = V_corr + phi_solv
 !compute solvation free energy (A_solv) and add A_cav
 A_corr = A_solv + A_cav
 !compute solvent charge density n_solv
 n_solv = n_b + n_ion
 
 dsol_io ' A_corr:  ', A_corr
 
 
 IF (LWRITE) THEN
   CALL WRITE_TO_FILE_RC('RHOB', n_b)
   CALL WRITE_TO_FILE_RC('RHOION', n_ion)
   CALL WRITE_TO_FILE_RC('PHI', phi_solv)
   CALL WRITE_TO_FILE_RC('PHISOLV', phi_solv)
   CALL WRITE_TO_FILE_RL('S', S)
   CALL STOP_TIMING("SOL",IU6,"WRITE")
 ENDIF
 

END SUBROUTINE LPCM_solver


!>=======================SUBROUTINE CREATECAVITY======================
!!
!!initializes epsilon in the grid
!!calculates the cavitation energy and cavitation potential
!!Acav and Vcav
!!
!!====================================================================
SUBROUTINE CREATECAVITY(n_e, Vcav, Acav)

 COMPLEX(q), INTENT(IN)  :: n_e(:)
 COMPLEX(q), INTENT(OUT) :: Vcav(:)
 REAL(q),    INTENT(OUT) :: Acav
 
 COMPLEX(q) :: CWORK(:)
 REAL(q)    :: RWORK(:), gradS(:,:), qs, qv, n_shift
 INTEGER    :: I
 
 ALLOCATABLE :: CWORK, RWORK, gradS
 
 ALLOCATE(CWORK(GRIDC%MPLWV), RWORK(DIMREAL(GRIDC%MPLWV)), gradS(DIMREAL(GRIDC%MPLWV),4))

 
 !----------------------------------------------------------------------------
 !SHAPE_FUNC returns shape function and its derivative wrt the charge density
 !charge density must be in real space (unscaled), result in real space
 !----------------------------------------------------------------------------
 
 ! cavity shape function
 CWORK = n_e/LATT_CUR%OMEGA
 CALL FFT3D_RC2RL(CWORK,RWORK)
 CALL SHAPE_FUNC(RWORK, S, Sprime)
 
 ! compute grad(S) in gradS
 ! gradS(:,4) contains |grad(S)|
 CALL FFT3D_RL2RC(S,CWORK)
 CALL GRAD(CWORK, gradS)
 
 !-------------------------------------------------------------
 !cavitation energy
 !Acav = tau * \int_dr^3 |grad(S)|
 !-------------------------------------------------------------
 qs = SUM_GRID_RL(gradS(:,4)) * LATT_CUR%OMEGA/GRIDC%NPLWV
 qv = SUM_GRID_RL(S)/GRIDC%NPLWV
 qv = LATT_CUR%OMEGA*(1._q-qv)
 Acav = tau*qs
 
 dsol_io ' area:   ', qs
 dsol_io ' volume: ', qv
 dsol_io ' Acav:   ', Acav

 ! compute cavity contribution to potential
 ! Vcav = -tau * div((grad(S))/|grad(S)|)
 RWORK = 1._q/(gradS(:,4)+eps)
 ! returns div(gradS/RWORK) in CWORK
 CALL DIV(gradS, CWORK, RWORK)
 CALL FFT3D_RC2RL(CWORK,RWORK)
 RWORK = -tau*Sprime*RWORK
 CALL FFT3D_RL2RC(RWORK,Vcav)
 

END SUBROUTINE CREATECAVITY


!>=====================SUBROUTINE SOL_pot=================================
 ! calculate LPB correction to KS potential
 ! also computes bound and counterion charges
!!========================================================================
SUBROUTINE SOL_pot(V_solv, n_b, n_ion)
 
 COMPLEX(q), INTENT(INOUT) :: V_solv(:)
 COMPLEX(q), INTENT(OUT)   :: n_b(:), n_ion(:)
 
 COMPLEX(q), ALLOCATABLE :: CWORK(:)
 REAL(q),    ALLOCATABLE :: RWORK(:), VWORK(:), RWORKV(:,:)
 
 ALLOCATE(CWORK(GRIDC%MPLWV), RWORK(DIMREAL(GRIDC%MPLWV)))
 ALLOCATE(VWORK(DIMREAL(GRIDC%MPLWV)), RWORKV(DIMREAL(GRIDC%MPLWV),4))
 
 ! compute V_diel in VWORK
 CALL GRAD(phi, RWORKV)
 VWORK = -0.5_q/EDEPS * (eb_k-1._q) * RWORKV(:,4)**2
 ! compute bound charge in n_b
 CALL DIV(RWORKV, n_b, S)
 n_b = n_b * (eb_k-1._q) * LATT_CUR%OMEGA/EDEPS
 
 IF (LION) THEN
   ! add V_ion to VWORK
   CALL FFT3D_RC2RL(phi, RWORK)
   VWORK = VWORK - 0.5_q/EDEPS * ekappa2_b * RWORK**2
   ! compute ionic charge in n_ion
   RWORK = -S*RWORK
   CALL FFT3D_RL2RC(RWORK, n_ion)
   n_ion = n_ion * ekappa2_b * LATT_CUR%OMEGA/EDEPS
 ENDIF
 
 VWORK = Sprime * VWORK
 CALL FFT3D_RL2RC(VWORK,CWORK)
 V_solv = V_solv + CWORK
 
END SUBROUTINE SOL_pot


!>======================SUBROUTINE SHAPE_FUNC=========================
!!
!!calculates the cavity shape function and its first derivative wrt n_e 
!!in realspace
!!
!!====================================================================
SUBROUTINE SHAPE_FUNC(n_e, S, Sprime)

 REAL(q), INTENT(IN)  :: n_e(:)
 REAL(q), INTENT(OUT) :: S(:), Sprime(:)

 !local variables
 REAL(q) :: X(:), Z(:), C
 
 ALLOCATABLE :: X, Z
 
 ALLOCATE(X(DIMREAL(GRIDC%MPLWV)), Z(DIMREAL(GRIDC%MPLWV)))
 
 C = 1._q/(sqrt(2._q)*sigma_k)
 X = MAX(n_e,eps)/nc_k
 Z = C*log(X)
 
 S = 0.5_q*ERFC(Z)
 Sprime = -C/(nc_k*sqrt(PI)) * EXP(-Z**2)/X
 
END SUBROUTINE SHAPE_FUNC


!******************** SUBROUTINE LPCM_READER *************************
!
! Reads in the solvation model parameters
!
!**********************************************************************
SUBROUTINE LPCM_READER(IO)

 USE vaspxml

 TYPE (in_struct), INTENT(IN) :: IO
 
 INTEGER      :: IDUM, N, IERR, IU5
 REAL(q)      :: RDUM
 COMPLEX(q)   :: CDUM
 LOGICAL      :: LOPEN,LDUM
 CHARACTER(1) :: CHARAC
 
 IU0 = IO%IU0
 IU6 = IO%IU6
 IU5 = IO%IU5

 LOPEN = .FALSE.

 OPEN(UNIT=IU5,FILE='INCAR',STATUS='OLD')      

 
 !---------------------------------------------------------------
 !sigma_k
 !--------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'SIGMA_K','=','#',';','F',IDUM,sigma_k,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''SIGMA_K''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('SIGMA_K','F',IDUM,sigma_k,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !nc_k
 !---------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'NC_K','=','#',';','F',IDUM,nc_k,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''NC_K''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('NC_K','F',IDUM,nc_k,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !tau
 !---------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'TAU','=','#',';','F',IDUM,tau,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''TAU''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('TAU','F',IDUM,tau,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !eb_k
 !---------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'EB_K','=','#',';','F',IDUM,eb_k,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''EB_K''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('EB_K','F',IDUM,eb_k,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !lambda_d_k = debye screening length in angstrom
 !--------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'LAMBDA_D_K','=','#',';','F',IDUM,lambda_d_k,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''LAMBDA_D_K''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 
 
 CLOSE(IU5)
 
END SUBROUTINE LPCM_READER


!******************** SUBROUTINE LPCM_INIT ****************************
!      
! initialize LPCM module
!
!***********************************************************************
SUBROUTINE LPCM_INIT
 
 ! intialize real and rec space grid data
 CALL SETGRID
 
 IF (lambda_d_k > 0._q) THEN
   ekappa2_b = eb_k/lambda_d_k**2
   LION = .TRUE.
 ELSE
   ekappa2_b = 0._q
   LION = .FALSE.
 ENDIF
 
 ALLOCATE(phi(GRIDC%MPLWV), phi_solv(GRIDC%MPLWV), S(DIMREAL(GRIDC%MPLWV)), Sprime(DIMREAL(GRIDC%MPLWV)))
 
 phi_solv = 0._q

END SUBROUTINE LPCM_INIT


!******************** SUBROUTINE LPCM_WRITER **************************
!      
! writes the solvation model parameters to the OUTCAR file
!
!***********************************************************************
SUBROUTINE LPCM_WRITER(IO)
 
 TYPE (in_struct), INTENT(IN) :: IO
 
 IF (IO%IU6>=0) THEN
   WRITE(IO%IU6,100) eb_k, sigma_k, nc_k, tau, lambda_d_k
 ENDIF

100 FORMAT( &
            ' Linear solvation model:'/ &
            '   EB_K        =',F10.6,  '     relative permittivity of the bulk solvent' /&
            '   SIGMA_K     =',F10.6,  '     width of the dielectric cavity' /&
            '   NC_K        =',F10.6,  '     cutoff charge density' /&
            '   TAU         =',F10.6,  '     cavity surface tension in eV/Angstrom**2' /&        
            '   LAMBDA_D_K  =',F10.6,  '     Debye length in Angstrom' / )


END SUBROUTINE LPCM_WRITER



!******************** SUBROUTINE XML_WRITE_LPCM ***********************
!      
!
! writes the solvation model parameters to vasprun.xml
!
!
!***********************************************************************
SUBROUTINE XML_WRITE_LPCM

 USE vaspxml
 
 INTEGER      :: IDUM
 REAL(q)      :: RDUM
 COMPLEX(q)   :: CDUM
 LOGICAL      :: LDUM
 CHARACTER(1) :: CHARAC
 
 CALL XML_TAG("separator","solvation model")
 
 CALL XML_INCAR('SIGMA_K','F',IDUM,sigma_k,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('NC_K','F',IDUM,nc_k,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('EB_K','F',IDUM,eb_k,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('TAU','F',IDUM,tau,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('LAMBDA_D_K','F',IDUM,lambda_d_k,CDUM,LDUM,CHARAC,1)
 
 CALL XML_CLOSE_TAG

END SUBROUTINE XML_WRITE_LPCM


!>=======================SUBROUTINE SET_SOL_WRITE=========================
!!
!! write ionic and bound charges to files
!!
!!====================================================================
SUBROUTINE SET_SOL_WRITE_LPCM
 LWRITE = .TRUE.
END SUBROUTINE SET_SOL_WRITE_LPCM


END MODULE solvation_lpcm









MODULE solvation_nlpcm
 
 USE prec
 USE base
 USE constant 
 USE mpimy 
 USE IFCORE
 USE solvation_grid
 
 IMPLICIT NONE
 
 PRIVATE

 PUBLIC :: NLPCM_solver, SET_SOL_WRITE_NLPCM, qion_max
 PUBLIC :: NLPCM_READER, NLPCM_INIT, NLPCM_WRITER, XML_WRITE_NLPCM
 
 SAVE
 
 INTEGER :: IU0
 
 REAL(q) :: A_sol, A_solv
 
 COMPLEX(q), ALLOCATABLE :: phi(:), phi_solv(:), n_b(:), n_ion(:)
 REAL(q),    ALLOCATABLE :: E_loc(:,:), P(:,:), lambda_diel(:), lambda_ion(:), lambda_cav(:)
 COMPLEX(q), ALLOCATABLE :: w_solv(:), w_cav(:), w_ion(:), w_diel(:), w(:), w_a(:), w_b(:)
 REAL(q),    ALLOCATABLE :: S_cav(:), S_diel(:), S_ion(:)
 REAL(q),    ALLOCATABLE :: Xcav(:,:)
 
 REAL(q) :: qion_max = 0._q
 
 LOGICAL :: LNLDIEL=.TRUE.
 LOGICAL :: LNLION=.TRUE.
 LOGICAL :: LNLTEST=.FALSE.
 
 INTEGER :: I_nloc_sol=1
 
 REAL(q) :: eps = 1.E-10_q
 REAL(q) :: n_min = 1.E-4_q
 REAL(q) :: a_max = 100._q
 
 REAL(q) :: sigma_k = 0.6_q
 REAL(q) :: a_k = 0.125_q
 REAL(q) :: nc_k = 0.015_q
 REAL(q) :: tau = 8.79E-4_q
 REAL(q) :: r_solv = 1.4_q
 REAL(q) :: r_cav = 0._q
 REAL(q) :: r_diel = 1.0_q
 REAL(q) :: r_ion = 0._q
 REAL(q) :: r_b = 0._q
 REAL(q) :: eb_k = 78.4_q
 REAL(q) :: lambda_d_k = 0._q
 REAL(q) :: SolTemp = 298._q
 REAL(q) :: ZION = 1._q
 REAL(q) :: d_ion = -1._q
 REAL(q) :: c_molar = 0._q
 REAL(q) :: n_mol = 0.0335_q
 REAL(q) :: p_mol = 0.50_q
 REAL(q) :: epsilon_inf = 1.78_q
 
 REAL(q) :: invBETA
 REAL(q) :: PBETA, alpha0_rot, alpha_pol, invalpha_sic
 REAL(q) :: ZBETA, n_max, alpha0_ion, theta_b
 LOGICAL :: LION
 LOGICAL :: LWRITE=.FALSE.

CONTAINS


!>=====================SUBROUTINE NLPCM_solver========================
!!
!! solve the NLPB equations for phi
!! also compute the correction to the KS potential and free energy
!!
!!====================================================================
SUBROUTINE NLPCM_solver(phi_sol, n_e, q_sol, tol, V_corr, n_solv, A_corr)

 COMPLEX(q), INTENT(IN)  :: phi_sol(:)
 COMPLEX(q), INTENT(IN)  :: n_e(:)
 REAL(q),    INTENT(IN)  :: q_sol, tol
 COMPLEX(q), INTENT(OUT) :: V_corr(:)
 COMPLEX(q), INTENT(OUT) :: n_solv(:)
 REAL(q),    INTENT(OUT) :: A_corr
 
 COMPLEX(q) :: dphi(:), resid(:)
 REAL(q)    :: ekappa2(:), chi(:,:,:)
 REAL(q)    :: tol_L, rms, A_LPB, A_cav, phi_G0, f_ion, q_ion, dphi0
 INTEGER    :: iter, NCGSOL, NCGSOL_L
 
 COMPLEX(q) :: CWORK(GRIDC%MPLWV)
 REAL(q)    :: RWORKV(DIMREAL(GRIDC%MPLWV),3)
 
 ALLOCATABLE :: dphi, resid, ekappa2, chi
 
 ALLOCATE(dphi(GRIDC%MPLWV), resid(GRIDC%MPLWV))
 IF (LNLDIEL) THEN
   ALLOCATE(chi(DIMREAL(GRIDC%MPLWV),3,3))
 ELSE
   ALLOCATE(chi(DIMREAL(GRIDC%MPLWV),1,1))
 ENDIF
 IF (LION) THEN
   ALLOCATE(ekappa2(DIMREAL(GRIDC%MPLWV)))
 ELSE
   ALLOCATE(ekappa2(0))
 ENDIF
 
 ! create solvent cavity S and compute cavity formation potential and energy
 dsol_io 'Creating dielectric and ionic cavities'
 CALL CREATECAVITY(n_e, A_cav)
 
 ! check whether the electrolyte region is large enough to accommodate the necessary counter charge
 CALL FFT3D_RL2RC(S_ion, CWORK)
 CALL GET_G0(CWORK, f_ion)
 qion_max = f_ion * n_max*LATT_CUR%OMEGA
 dsol_io ' qion_max:  ', qion_max
 IF (abs(q_sol) > qion_max*0.99_q .AND. LNLION) THEN
   IF (IU0>=0) THEN
      WRITE(IU0,*) 'ERROR: NLPCM_solver: your ''vacuum'' region is not large enough to accommodate the necessary counter charge'
      WRITE(IU0,*) 'The counter charge is exceeding 99% of the saturation value'
      WRITE(IU0,*) 'Maximum charge: ', qion_max
      WRITE(IU0,*) 'Surface/solute charge: ', q_sol
   ENDIF
   STOP
 ENDIF
 
 !phi = phi_sol + phi_solv
 phi_solv = phi - phi_sol
 
 ! solve NLPB equations for phi
 dsol_io ' Solving nonlinear dielectric/ionic screening model'
 dsol_io ' Starting outer loop'
 ifdsol IF(IU0>=0) WRITE(IU0,'(2X,6X,A4,2A24)') 'iter', 'A_solv', 'rms'
 NCGSOL = 0
 iter = 0
 ! calculate E_loc, n_b, n_ion, lambda, for initial phi
 CALL Update_NLPB(q_sol)
 ! outer loop that solves the NLPB equations using Newton's method
 DO
   ! compute residual
   ! resid = q_sol + n_b(phi) + n_ion(phi) - L0(phi_solv)
   CALL NLPB_resid(resid, q_sol, rms)
   ifdsol IF(IU0>=0) WRITE(IU0,'(2X,A6,I4,2E24.15)') 'NLPB: ', iter, A_solv, rms
   IF (rms < tol .AND. iter>=1) EXIT
   iter = iter + 1
   IF (iter > 100) THEN
      IF (IU0>=0) WRITE(IU0,*) 'ERROR: NLPCM_solver: outer loop failed to converge in 100 iterations'
      STOP
   ENDIF
   ! compute linearized response functions chi and ekappa2
   CALL SET_LPB(chi, ekappa2)
   dphi = 0._q
   CALL GET_G0(n_ion, q_ion)
   ! adjust the G=0 component of phi if the counter charge is too close to saturation
   IF (abs(q_ion) > qion_max*0.99_q .AND. LNLION) THEN
      dphi0 = -D_PROD_RC(CWORK, phi)/f_ion
      dsol_io ' Adjusting G=0 component of phi, q_ion = ', q_ion, ' qion_max = ', qion_max, ' dphi0 = ', dphi0
      CALL SET_G0(dphi, dphi0)
   ELSE
      A_LPB = A_solv
      !tol_L = MIN(rms/10._q,tol)
      tol_L = rms/10._q
      ! solve L(dphi) = resid
      CALL MINIMIZE_L(resid, chi, ekappa2, w_b, dphi, tol_L, A_LPB, NCGSOL_L)
      NCGSOL = NCGSOL + NCGSOL_L
   ENDIF
   ! update phi and E_loc and return free energy
   CALL UPDATE_PHI(dphi, resid, q_sol)
 ENDDO
 dsol_io ' Finished solving NLPB equations'
 
 !Computes V_solv and adds to V_corr
 CALL SOL_pot(V_corr)
 
 IF (LWRITE) CALL WRITE_TO_FILE_RC('VSOLV', V_corr)
 
 !add electrostatic potential correction (phi_solv) to V_corr
 V_corr = V_corr + phi_solv
 !compute solvation free energy (A_solv) and add A_cav
 A_corr = A_solv + A_cav
 !compute solvent charge density n_solv
 n_solv = n_b + n_ion
 
 IF (LWRITE) THEN
   !CALL GET_G0(phi, phi_G0)
   !CALL SET_G0(phi_sol, phi_G0)
   CALL WRITE_TO_FILE_RC('PHI', phi)
   !CALL WRITE_TO_FILE_RC('PHISOL', phi_sol)
   CALL WRITE_TO_FILE_RC('PHISOLV', phi_solv)
   CALL WRITE_TO_FILE_RL('ELOC', E_loc(:,3))
   CALL WRITE_TO_FILE_RL('P', P(:,3))
   CALL WRITE_TO_FILE_RC('RHOB', n_b)
   CALL WRITE_TO_FILE_RC('RHOION', n_ion)
   !write Ez
   !CALL GRAD(phi, RWORKV)
   !RWORKV = -RWORKV
   !CALL WRITE_TO_FILE_RL('Ez', RWORKV(:,3))
   !write Dz
   !CALL L0_inv_op(n_ion, CWORK)
   !CWORK = CWORK + phi_sol
   !CALL GRAD(CWORK, RWORKV)
   !RWORKV = -RWORKV
   !CALL WRITE_TO_FILE_RL('Dz', RWORKV(:,3))
   !write Dz-Ez
   !CALL L0_inv_op(n_b, CWORK)
   !CALL GRAD(CWORK, RWORKV)
   !RWORKV = -RWORKV
   !CALL WRITE_TO_FILE_RL('Dz-Ez', RWORKV(:,3))
 ENDIF
 
END SUBROUTINE NLPCM_solver
 

!>=====================SUBROUTINE UPDATE_PHI============================
!!
!!
!!
!!======================================================================
SUBROUTINE UPDATE_PHI(dphi, resid, q_sol)
 
 COMPLEX(q), INTENT(IN) :: dphi(:)
 COMPLEX(q), INTENT(IN) :: resid(:)
 REAL(q),    INTENT(IN) :: q_sol
 
 REAL(q) :: m, alpha, A_old
 INTEGER :: iter
 
 !! derivative test
 !REAL(q) :: m_0, m_diel, m_ion, Aold_0, Aold_diel, Aold_ion, A_0, A_diel, A_ion
 !COMPLEX(q) :: CWORK(GRIDC%MPLWV)
 !! derivative test
 
 REAL(q), PARAMETER :: gamma = 0.5_q
 REAL(q), PARAMETER :: c = 0.1_q
 
 !dsol_io '  Starting line search'
 
 !! derivative test
 ! compute L0(phi) in CWORK
 !CALL L0_op(phi,CWORK)
 !Aold_0 = D_PROD_RC(phi, n_sol) - 0.5_q*D_PROD_RC(phi, CWORK)
 !Aold_diel = SUM_GRID_RL(S_diel*lambda_diel) * LATT_CUR%OMEGA/GRIDC%NPLWV
 !Aold_ion = SUM_GRID_RL(S_ion*lambda_ion) * LATT_CUR%OMEGA/GRIDC%NPLWV
 !CWORK = n_sol - CWORK
 !m_0 = D_PROD_RC(CWORK,dphi)
 !m_diel = D_PROD_RC(n_b,dphi)
 !m_ion = D_PROD_RC(n_ion,dphi)
 !! derivative test
 
 A_old = A_solv
 m = D_PROD_RC(resid,dphi)
 alpha = 1._q
 
 phi = phi + dphi
 phi_solv = phi_solv + dphi
 
 iter = 0
 !ifdsol IF(IU0>=0) WRITE(IU0,'(3X,5X,A4,24X,3A24)') 'iter', 'A_solv', 'dA', 'min dA'
 
 DO
   CALL Update_NLPB(q_sol)
   ifdsol IF(IU0>=0) WRITE(IU0,'(3X,A5,I4,3E24.15,2F15.6)') ' LS: ', iter, A_solv, A_solv-A_old, alpha*m
   IF (A_solv >= A_old + alpha*c*m) EXIT
   iter = iter + 1
   IF (iter > 10) EXIT
   phi = phi - (1._q-gamma)*alpha*dphi
   phi_solv = phi_solv - (1._q-gamma)*alpha*dphi
   alpha = gamma*alpha
 ENDDO
 
 !! derivative test
 !CALL L0_op(phi,CWORK)
 !A_0 = D_PROD_RC(phi, n_sol) - 0.5_q*D_PROD_RC(phi, CWORK)
 !A_diel = n_mol*SUM_GRID_RL(S_diel*lambda_diel) * LATT_CUR%OMEGA/GRIDC%NPLWV
 !A_ion = n_max*SUM_GRID_RL(S_ion*lambda_ion) * LATT_CUR%OMEGA/GRIDC%NPLWV
 !ifdsol IF(IU0>=0) WRITE(IU0,'(12X,6E15.6)') A_0-Aold_0, alpha*m_0, A_diel-Aold_diel, alpha*m_diel, A_ion-Aold_ion, alpha*m_ion
 !! derivative test
 
 !dsol_io '  Line search finished'

END SUBROUTINE UPDATE_PHI


!>=====================SUBROUTINE NLPB_resid==============================
 ! calculate NLPB residual
!!========================================================================
SUBROUTINE NLPB_resid(resid, q_sol, rms)

 COMPLEX(q), INTENT(OUT) :: resid(:)
 REAL(q),    INTENT(IN)  :: q_sol
 REAL(q),    INTENT(OUT) :: rms
 
 OPTIONAL :: rms
 
 COMPLEX(q) :: CWORK(GRIDC%MPLWV)

 ! compute L0(phi_solv) - q_sol in resid
 CALL L0_op(phi_solv,resid)
 CALL SET_G0(resid, -q_sol)
 ! resid = q_sol + n_b + n_ion - L0(phi_solv)
 resid = n_b + n_ion - resid
 ! compute rms of resid
 IF (PRESENT(rms)) THEN
   CALL L0_inv_op(resid, CWORK)
   rms = 0._q
   IF (LION) THEN
      CALL GET_G0(resid, rms)
      dsol_io 'resid(G=0) = ', rms
   ENDIF
   rms = sqrt( rms**2 + D_PROD_RC(CWORK,CWORK) )
 ENDIF

END SUBROUTINE NLPB_resid


!>=======================SUBROUTINE CREATECAVITY======================
!!
!!initializes S_cav, S_diel, and S_ion on the real space grid
!!
!!====================================================================
SUBROUTINE CREATECAVITY(n_e, A_cav)

 COMPLEX(q), INTENT(IN)  :: n_e(:)
 REAL(q),    INTENT(OUT) :: A_cav
 
 COMPLEX(q) :: CWORK(:)
 REAL(q)    :: RWORK(:), SWORK(:), gradS(:,:)
 
 REAL(q) :: qs, qv
 
 ALLOCATABLE :: CWORK, RWORK, SWORK, gradS
 
 ALLOCATE(CWORK(GRIDC%MPLWV), SWORK(DIMREAL(GRIDC%MPLWV)), RWORK(DIMREAL(GRIDC%MPLWV)))
 ALLOCATE(gradS(DIMREAL(GRIDC%MPLWV),4))
 
 
 ! determine vdW cavity
 ! compute xcav_vdW
 CWORK = n_e/LATT_CUR%OMEGA / nc_k 
 CALL FFT3D_RC2RL(CWORK,RWORK)
 !IF (LWRITE) CALL WRITE_TO_FILE_RL('NE', RWORK)
 RWORK = MAX(RWORK,n_min)
 Xcav(:,1) = LOG(RWORK)
 !IF (LWRITE) CALL WRITE_TO_FILE_RL('XVDW', Xcav(:,1))
 ! compute S_vdW in SWORK
 CALL SHAPE_FUNC(Xcav(:,1), SWORK)
 
 IF (LWRITE) CALL WRITE_TO_FILE_RL('SVDW', SWORK)
 
 ! determine ionic cavity
 ! S_vdW is stored in SWORK
 IF (SIZE(w_ion) > 0) THEN
   ! compute xcav_ion
   RWORK = 1._q - SWORK
   CALL CONV(RWORK, w_ion)
   RWORK = MAX(RWORK,n_min)
   Xcav(:,3) = LOG(RWORK)
   !IF (LWRITE) CALL WRITE_TO_FILE_RL('XION', Xcav(:,3))
   ! compute S_ion
   CALL SHAPE_FUNC(Xcav(:,3), S_ion)
 ELSE
   S_ion = SWORK
 ENDIF
 
 IF (LWRITE) CALL WRITE_TO_FILE_RL('SION', S_ion)
 
 ! determine solvent accessible cavity
 ! S_vdW is stored in SWORK
 IF (SIZE(w_solv) > 0) THEN
   ! compute xcav_vdW
   RWORK = 1._q - SWORK
   CALL CONV(RWORK, w_solv)
   RWORK = MAX(RWORK,n_min)
   Xcav(:,2) = LOG(RWORK)
   !IF (LWRITE) CALL WRITE_TO_FILE_RL('XSOLV', Xcav(:,2))
   ! compute S_solv in SWORK
   CALL SHAPE_FUNC(Xcav(:,2), SWORK)
 ELSE
   SWORK = SWORK
 ENDIF
 
 IF (LWRITE) CALL WRITE_TO_FILE_RL('SSOLV', SWORK)
 
 ! determine solvent excluded cavity
 ! S_solv is stored in SWORK
 IF (SIZE(w_cav) > 0) THEN
   ! compute xcav_cav
   RWORK = SWORK
   CALL CONV(RWORK, w_cav)
   RWORK = MAX(RWORK,n_min)
   Xcav(:,4) = LOG(RWORK)
   !IF (LWRITE) CALL WRITE_TO_FILE_RL('XCAV', Xcav(:,4))
   ! compute S_cav
   CALL SHAPE_FUNC(Xcav(:,4), S_cav)
   S_cav = 1._q - S_cav
 ELSE
   S_cav = SWORK
 ENDIF
 
 IF (LWRITE) CALL WRITE_TO_FILE_RL('SCAV', S_cav)
 
 ! determine dielectric cavity
 ! S_solv is stored in SWORK
 IF (SIZE(w_diel) > 0) THEN
   ! compute xcav_diel
   RWORK = SWORK
   CALL CONV(RWORK, w_diel)
   RWORK = MAX(RWORK,n_min)
   Xcav(:,5) = LOG(RWORK)
   !IF (LWRITE) CALL WRITE_TO_FILE_RL('XDIEL', Xcav(:,5))
   ! compute S_diel
   CALL SHAPE_FUNC(Xcav(:,5), S_diel)
   S_diel = 1._q - S_diel
 ELSE
   S_diel = SWORK
 ENDIF
 
 IF (LWRITE) CALL WRITE_TO_FILE_RL('SDIEL', S_diel)
 
 
 !-------------------------------------------------------------
 !cavitation energy
 !Acav = tau * \int_dr^3 |grad(S_cav)|
 !-------------------------------------------------------------
 
 ! compute grad(S_cav) in gradS
 CALL FFT3D_RL2RC(S_cav, CWORK)
 CALL GRAD(CWORK, gradS)
 ! compute lambda_cav in RWORK
 RWORK = -tau/n_mol / (gradS(:,4)+eps)
 CALL DIV(gradS, CWORK, RWORK)
 CALL FFT3D_RC2RL(CWORK,lambda_cav)
 
 qs = SUM_GRID_RL(gradS(:,4)) * LATT_CUR%OMEGA/GRIDC%NPLWV
 qv = SUM_GRID_RL(S_cav)/GRIDC%NPLWV
 qv = LATT_CUR%OMEGA*(1._q-qv)
 A_cav = tau*qs
 
 dsol_io ' area:   ', qs
 dsol_io ' volume: ', qv
 dsol_io ' A_cav:  ', A_cav
 
 
END SUBROUTINE CREATECAVITY


!>=====================SUBROUTINE SOL_pot=================================
 ! calculate NLPB correction to KS potential
!!========================================================================
SUBROUTINE SOL_pot(V_solv)
 
 COMPLEX(q), INTENT(OUT) :: V_solv(:)
 
 REAL(q),    ALLOCATABLE :: RWORK(:)
 COMPLEX(q), ALLOCATABLE :: CWORK(:)
 
 ! ! derivative test
 ! REAL(q), ALLOCATABLE, SAVE :: Xcav_old(:,:)
 ! REAL(q), ALLOCATABLE, SAVE :: S_cav_old(:), S_diel_old(:), S_ion_old(:)
 ! REAL(q), ALLOCATABLE, SAVE :: n(:), n_old(:), V_cav(:), V_diel(:), V_ion(:)
 ! REAL(q) :: dE1, dE2
 ! LOGICAL, SAVE :: initialized = .FALSE.
 ! ! derivative test
 
 ALLOCATE(RWORK(DIMREAL(GRIDC%MPLWV)), CWORK(GRIDC%MPLWV))
 
 ! ! derivative test
 ! If (.NOT. initialized) THEN
   ! ALLOCATE(S_diel_old(DIMREAL(GRIDC%MPLWV)), S_ion_old(DIMREAL(GRIDC%MPLWV)), S_cav_old(DIMREAL(GRIDC%MPLWV)))
   ! ALLOCATE(Xcav_old(DIMREAL(GRIDC%MPLWV),5))
   ! ALLOCATE(V_cav(DIMREAL(GRIDC%MPLWV)), V_diel(DIMREAL(GRIDC%MPLWV)), V_ion(DIMREAL(GRIDC%MPLWV)))
   ! ALLOCATE(n(DIMREAL(GRIDC%MPLWV)), n_old(DIMREAL(GRIDC%MPLWV)))
   ! Xcav_old = 0._q
   ! S_cav_old = 0._q
   ! S_diel_old = 0._q
   ! S_ion_old = 0._q
   ! n_old = 0._q
   ! initialized = .TRUE.
 ! ENDIF
 ! ! derivative test
 
 ! compute cavity contribution to potential
 RWORK = n_mol*lambda_cav
 IF (SIZE(w_cav) > 0) THEN
   ! compute V_cav_3
   CALL SHAPE_POT(Xcav(:,4), RWORK)
   ! compute lambda_cav_2
   CALL XCORR(RWORK, w_cav)
   RWORK = -RWORK
 ENDIF
 IF (SIZE(w_solv) > 0) THEN
   ! compute V_cav_2
   CALL SHAPE_POT(Xcav(:,2), RWORK)
   ! compute lambda_cav_1
   CALL XCORR(RWORK, w_solv)
   RWORK = -RWORK
 ENDIF
 ! compute V_cav
 CALL SHAPE_POT(Xcav(:,1), RWORK)
 !V_cav = RWORK
 ! add to V_solv
 CALL FFT3D_RL2RC(RWORK,V_solv)
 !IF (LWRITE) CALL WRITE_TO_FILE_RL('VCAV', RWORK)
 
 ! contribution from S_diel
 RWORK = n_mol*lambda_diel
 IF (SIZE(w_diel) > 0) THEN
   ! compute V_diel_3
   CALL SHAPE_POT(Xcav(:,5), RWORK)
   ! compute lambda_diel_2
   CALL XCORR(RWORK, w_diel)
   RWORK = -RWORK
 ENDIF
 IF (SIZE(w_solv) > 0) THEN
   ! compute V_diel_2
   CALL SHAPE_POT(Xcav(:,2), RWORK)
   ! compute lambda_diel_1
   CALL XCORR(RWORK, w_solv)
   RWORK = -RWORK
 ENDIF
 ! compute V_diel
 CALL SHAPE_POT(Xcav(:,1), RWORK)
 !V_diel = RWORK
 ! add to V_solv
 CALL FFT3D_RL2RC(RWORK,CWORK)
 V_solv = V_solv + CWORK
 !IF (LWRITE) CALL WRITE_TO_FILE_RL('VDIEL', RWORK)
 
 ! contribution from S_ion
 RWORK = n_max*lambda_ion
 IF (SIZE(w_ion) > 0) THEN
   ! compute V_ion_2
   CALL SHAPE_POT(Xcav(:,3), RWORK)
   ! compute lambda_ion_1
   CALL XCORR(RWORK, w_ion)
   RWORK = -RWORK
 ENDIF
 ! compute V_ion
 CALL SHAPE_POT(Xcav(:,1), RWORK)
 !V_ion = RWORK
 ! add to V_solv
 CALL FFT3D_RL2RC(RWORK,CWORK)
 V_solv = V_solv + CWORK
 !IF (LWRITE) CALL WRITE_TO_FILE_RL('VION', RWORK)
 
 V_solv = V_solv/nc_k
 
 !Xcav_old = Xcav
 
 
 ! ! derivative test
 ! n = exp(Xcav(:,1))
 ! dE1  = SUM_GRID_RL((S_cav-S_cav_old)*n_mol*lambda_cav) * LATT_CUR%OMEGA/GRIDC%NPLWV
 ! dE2  = SUM_GRID_RL((n-n_old)*V_cav) * LATT_CUR%OMEGA/GRIDC%NPLWV
 ! ifdsol IF(IU0>=0) WRITE(IU0,'(12X,2E15.6)') dE1, dE2
 ! dE1  = SUM_GRID_RL((S_diel-S_diel_old)*n_mol*lambda_diel) * LATT_CUR%OMEGA/GRIDC%NPLWV
 ! dE2  = SUM_GRID_RL((n-n_old)*V_diel) * LATT_CUR%OMEGA/GRIDC%NPLWV
 ! ifdsol IF(IU0>=0) WRITE(IU0,'(12X,2E15.6)') dE1, dE2
 ! dE1  = SUM_GRID_RL((S_ion-S_ion_old)*n_max*lambda_ion) * LATT_CUR%OMEGA/GRIDC%NPLWV
 ! dE2  = SUM_GRID_RL((n-n_old)*V_ion) * LATT_CUR%OMEGA/GRIDC%NPLWV
 ! ifdsol IF(IU0>=0) WRITE(IU0,'(12X,2E15.6)') dE1, dE2
 ! S_cav_old = S_cav
 ! S_diel_old = S_diel
 ! S_ion_old = S_ion
 ! n_old = n
 ! ! derivative test
 
 
END SUBROUTINE SOL_pot



!>=====================SUBROUTINE Update_NLPB=============================
 ! computes:
 !    phi and E_loc
 !    n_b and n_ion
 !    lambda_diel and lambda_ion
 !    A_solv
!!========================================================================
SUBROUTINE Update_NLPB(q_sol)

 REAL(q), INTENT(IN) :: q_sol
 
 REAL(q),    ALLOCATABLE :: RWORK(:), RWORKV(:,:), X(:), n_work(:)
 COMPLEX(q), ALLOCATABLE :: CWORK(:)
 
 INTEGER :: I
 
 ALLOCATE(RWORK(DIMREAL(GRIDC%MPLWV)), RWORKV(DIMREAL(GRIDC%MPLWV),3), CWORK(GRIDC%MPLWV))
 ALLOCATE(X(DIMREAL(GRIDC%MPLWV)))
 ALLOCATE(n_work(DIMREAL(GRIDC%MPLWV)))
 
 ! calculate local field
 CALL SOLVE_E_LOC
 
 !--------------------------
 ! ionic contribution
 !--------------------------
 
 ! compute X = ZBETA*phi
 CALL FFT3D_RC2RL(phi, RWORK)
 X = ZBETA*RWORK
 ! compute x*g_ion(x) in n_work
 IF (LNLION .AND. theta_b > 0._q) THEN
   WHERE (ABS(X) > 100._q)
      n_work = SIGN(1._q,X)
      lambda_ion = ( LOG(0.5_q*theta_b) + ABS(X) )
   ELSEWHERE
      WHERE (ABS(X) < sqrt(theta_b)*2e-4_q)
         n_work = theta_b*X
         lambda_ion = 0.5_q*theta_b*X**2
      ELSEWHERE
         RWORK = 1._q + theta_b*(COSH(X)-1._q)
         n_work = theta_b*SINH(X)/RWORK
         lambda_ion = LOG(RWORK)
      ENDWHERE
   ENDWHERE
 ELSEIF (LNLION) THEN
   X =  MAX(-100._q,MIN(X,100._q))     ! let's avoid forming black holes
   n_work = SINH(X)
   WHERE (ABS(X) < 2e-4_q)
      lambda_ion = 0.5_q*X**2
   ELSEWHERE
      lambda_ion = COSH(X)-1._q
   ENDWHERE
 ELSE
   n_work = X
   lambda_ion = 0.5_q*X**2
 ENDIF
 
 ! compute n_ion and fft to rec space
 !n_work = -n_max*alpha0_ion/EDEPS/ZBETA * S_ion*n_work * LATT_CUR%OMEGA
 n_work = -n_max*invBETA*ZBETA * S_ion*n_work * LATT_CUR%OMEGA
 CALL FFT3D_RL2RC(n_work, n_ion)
 lambda_ion = -invBETA*lambda_ion
 
 
 !--------------------------
 ! dielectric contribution
 !--------------------------
 
 ! rotational part
 ! compute X = PBETA*|E_loc|
 X = PBETA*E_loc(:,4)
 ! compute g_rot(x) in n_work
 ! add lambda_rot to lambda_diel
 IF (LNLDIEL) THEN
   WHERE (X < 2e-4_q)
      n_work = 1._q
      lambda_diel = -invBETA/6._q * X**2
   ELSEWHERE (X > 100._q)
      RWORK  = 1._q/X
      n_work = 3._q*(1._q - RWORK)*RWORK
      lambda_diel = -invBETA*(X - LOG(2*X))
   ELSEWHERE
      RWORK  = 1._q/X
      n_work = 3._q*(1._q/TANH(X) - RWORK)*RWORK
      lambda_diel = -invBETA*LOG(SINH(X)*RWORK)
   ENDWHERE
 ELSE
   n_work = 1._q
   lambda_diel = -invBETA/6._q * X**2
 ENDIF
 
 !IF (LWRITE) CALL WRITE_TO_FILE_RL('g_rot', n_work)
 
 ! compute |p|/|E_loc| in n_work
 n_work = alpha0_rot/EDEPS * n_work + alpha_pol/EDEPS
 
 RWORK = n_mol*S_diel/(1._q/(n_work*EDEPS) - invalpha_sic)
 !IF (LWRITE) CALL WRITE_TO_FILE_RL('chi', RWORK)
 
 ! add linear and sic parts to lambda_diel in lambda_diel
 ! lambda_pol = 0.5*EDEPS/alpha_pol*p_pol**2
 ! p_pol = alpha_pol/EDEPS * |E_loc|
 lambda_diel = lambda_diel - 0.5_q*alpha_pol/EDEPS * E_loc(:,4)**2
 ! lambda_sic = 0.5*EDEPS*invalpha_sic*p**2
 lambda_diel = lambda_diel + 0.5_q*EDEPS*invalpha_sic * (n_work*E_loc(:,4))**2
 
 ! compute |P|/|E_loc| in n_work
 n_work = n_mol*S_diel*n_work
 
 ! compute bound charge
 ! n_b = -div(P)
 DO I = 1,3
   P(:,I) = n_work*E_loc(:,I)
   CALL CONV(P(:,I), w_b)
 ENDDO
 CALL DIV(E_loc,n_b,n_work)
 n_b = -w_b*n_b * LATT_CUR%OMEGA
 
 
 !--------------------------
 ! solvation free energy
 !--------------------------
 ! compute q_sol * <phi_solv>
 CALL GET_G0(phi_solv, A_solv)
 A_solv = q_sol * A_solv
 ! compute L0(phi) in CWORK
 CALL L0_op(phi_solv,CWORK)
 A_solv = A_solv - 0.5_q*D_PROD_RC(phi_solv, CWORK) &
      & + n_mol*SUM_GRID_RL(S_diel*lambda_diel) * LATT_CUR%OMEGA/GRIDC%NPLWV &
      & + n_max*SUM_GRID_RL(S_ion*lambda_ion) * LATT_CUR%OMEGA/GRIDC%NPLWV
 

END SUBROUTINE Update_NLPB


!>======================SUBROUTINE SET_LPB==========================
!!
!!Chi = n_mol*S_solv*alpha_3*[(1-alpha_3)^-1]
!!Nalpha is a 3*3 matrix
!!
!!==================================================================
SUBROUTINE SET_LPB(chi, ekappa2)
 
 REAL(q), INTENT(OUT) :: chi(:,:,:)
 REAL(q), INTENT(OUT) :: ekappa2(:)

 !local variables
 REAL(q), ALLOCATABLE :: chi_par(:), chi_perp(:), RWORK(:), RWORK2(:)
 
 INTEGER :: I,J
 
 ALLOCATE(chi_par(DIMREAL(GRIDC%MPLWV)), chi_perp(DIMREAL(GRIDC%MPLWV)))
 ALLOCATE(RWORK(DIMREAL(GRIDC%MPLWV)), RWORK2(DIMREAL(GRIDC%MPLWV)))
 
 
 !--------------------------
 ! compute ekappa2
 !--------------------------
 
 IF (LNLION) THEN
 
   ! compute X = ZBETA*phi in RWORK
   CALL FFT3D_RC2RL(phi, RWORK)
   RWORK = ZBETA*RWORK
   WHERE (ABS(RWORK) > 100._q)
      ekappa2 = 0._q
   ELSEWHERE
      WHERE (ABS(RWORK) < 2e-4_q)
         RWORK2 = 0.5_q*RWORK**2
      ELSEWHERE
         RWORK2 = COSH(RWORK)-1._q
      ENDWHERE
      ekappa2 = (1._q + (1._q-theta_b)*RWORK2)/(1._q + theta_b*RWORK2)**2
   ENDWHERE
   ! scale ekappa2 and multiply by cavity functions
   ekappa2 = n_max*alpha0_ion * S_ion*ekappa2
   
 ELSEIF (LION) THEN
 
   ekappa2 = n_max*alpha0_ion * S_ion
   
 ENDIF
 
 !--------------------------
 ! compute chi
 !--------------------------
 
 IF (LNLDIEL) THEN
 
   ! store x = PBETA*|E_loc| in RWORK
   RWORK = PBETA*E_loc(:,4)
   
   ! store parallel component of chi in chi_par
   ! store perpendicular component in chi_perp
   WHERE (RWORK < 2e-4_q)
      chi_par  = 1._q
      chi_perp = 1._q
   ELSEWHERE (RWORK > 100._q)
      chi_par  = 3._q * 1._q/(RWORK)**2
      chi_perp = 3._q * (1._q - 1._q/RWORK)/RWORK
   ELSEWHERE
      chi_par  = 3._q * (1._q/RWORK**2 - 1._q/SINH(RWORK)**2)
      chi_perp = 3._q * (1._q/TANH(RWORK) - 1._q/RWORK)/RWORK
   ENDWHERE
   chi_par  = alpha_pol + alpha0_rot*chi_par
   chi_perp = alpha_pol + alpha0_rot*chi_perp
   chi_par  = n_mol*S_diel/(1._q/chi_par  - invalpha_sic)
   chi_perp = n_mol*S_diel/(1._q/chi_perp - invalpha_sic)
   
   IF (ANY(chi_par < 0._q)) THEN
      IF (IU0>=0) WRITE(IU0,*) 'Internal ERROR: SET_LPB: chi_par < 0'
      STOP
   ENDIF
   IF (ANY(chi_perp < 0._q)) THEN
      IF (IU0>=0) WRITE(IU0,*) 'Internal ERROR: SET_LPB: chi_perp < 0'
      STOP
   ENDIF
   
   ! store 1/|E_loc|^2 in RWORK
   WHERE (RWORK < 2e-4_q)
      RWORK = 0._q
   ELSEWHERE
      RWORK = 1._q/E_loc(:,4)**2
   ENDWHERE
   
   ! compute chi
   DO I = 1,3
      ! anisotropic component
      DO J = 1,3
         RWORK2  = E_loc(:,I)*E_loc(:,J)*RWORK
         chi(:,I,J) = (chi_par - chi_perp)*RWORK2
      ENDDO
      ! isotropic component
      chi(:,I,I) = chi(:,I,I) + chi_perp
   ENDDO
   
 ELSE
   
   RWORK = alpha0_rot + alpha_pol
   chi(:,1,1) = n_mol*S_diel/(1._q/RWORK - invalpha_sic)
   
 ENDIF
 
END SUBROUTINE SET_LPB


!>======================SUBROUTINE SOLVE_E_LOC========================
!! Compute E_loc from phi
!!====================================================================
 SUBROUTINE SOLVE_E_LOC

 REAL(q), ALLOCATABLE :: x0(:), f_loc(:)
 
 INTEGER :: I
 REAL(q) :: a, b, tol
 REAL(q), PARAMETER :: tol_grot = 2e-8_q     ! machine precision of g_rot
 
 ALLOCATE(f_loc(DIMREAL(GRIDC%MPLWV)))
 
 a = 1._q/(1._q - alpha_pol*invalpha_sic)                     !lower limit of f_loc
 b = 1._q/(1._q - (alpha_pol+alpha0_rot)*invalpha_sic)        !upper limit of f_loc
 
 ! machine precision of F(f_loc)
 tol = alpha0_rot*invalpha_sic*b*tol_grot
 ! ensure that F(b) < 0 within machine precision
 b = b*(1._q+2._q*tol)
 
 ! compute -grad(x_corr(w,phi)) and store in E_loc
 CALL GRAD(-conjg(w_b)*phi, E_loc)
 
 ! solve for f_loc
 IF (LNLDIEL) THEN
   ALLOCATE(x0(DIMREAL(GRIDC%MPLWV)))
   x0 = PBETA*E_loc(:,4)
   CALL fzero_vec(F_floc, a, b, tol, f_loc, x0)
 ELSE
   f_loc = 1._q/(1._q - (alpha_pol+alpha0_rot)*invalpha_sic)
 ENDIF
 ! scale E_loc by f_loc
 DO I = 1,4
   E_loc(:,I) = f_loc*E_loc(:,I)
 END DO

 
CONTAINS


ELEMENTAL REAL(q) FUNCTION F_floc(f_loc, x0)
 REAL(q), INTENT(IN) :: f_loc, x0
 REAL(q) :: x, g_rot
 x = f_loc*x0
 ! g_rot has relative fp error of 1e-8
 IF (x < 2e-4_q) THEN
   g_rot = 1._q
 ELSE
   g_rot = 3._q*(x-TANH(x))/(x**2*TANH(x))
 ENDIF
 F_floc = 1._q /(1._q - (g_rot*alpha0_rot + alpha_pol)*invalpha_sic) - f_loc  
END FUNCTION F_floc
 
 
END SUBROUTINE SOLVE_E_LOC


!>======================SUBROUTINE fzero_vec==================================
!! Vectorized root finding suroutine to find zeros of an elemental function
!! Uses Brent's method - adapted from:
!!       https://people.math.sc.edu/Burkardt/f_src/brent/brent.f90
!!============================================================================
SUBROUTINE fzero_vec(f, a0, b0, tol, x, x0)
 
 USE prec
 
 REAL(q), INTENT(IN)  :: a0, b0, tol, x0(:)
 REAL(q), INTENT(OUT) :: x(:)
 
 INTERFACE
   ELEMENTAL REAL(q) FUNCTION f(x, x0)
   USE prec
      REAL(q), INTENT(IN) :: x, x0
   END FUNCTION f
 END INTERFACE
 
 REAL(q), ALLOCATABLE :: a(:), b(:), c(:), d(:), e(:), m(:), p(:), r(:), s1(:), s2(:), s3(:), fa(:), fb(:), fc(:)
 REAL(q), ALLOCATABLE :: tol_abs(:)
 
 REAL(q) :: tmp
 INTEGER :: ITER
 
 ALLOCATE(a(DIMREAL(GRIDC%MPLWV)), b(DIMREAL(GRIDC%MPLWV)), c(DIMREAL(GRIDC%MPLWV)))
 ALLOCATE(d(DIMREAL(GRIDC%MPLWV)), e(DIMREAL(GRIDC%MPLWV)), m(DIMREAL(GRIDC%MPLWV)))
 ALLOCATE(p(DIMREAL(GRIDC%MPLWV)), r(DIMREAL(GRIDC%MPLWV)))
 ALLOCATE(s1(DIMREAL(GRIDC%MPLWV)), s2(DIMREAL(GRIDC%MPLWV)), s3(DIMREAL(GRIDC%MPLWV)))
 ALLOCATE(fa(DIMREAL(GRIDC%MPLWV)), fb(DIMREAL(GRIDC%MPLWV)), fc(DIMREAL(GRIDC%MPLWV)))
 ALLOCATE(tol_abs(DIMREAL(GRIDC%MPLWV)))
 
 ! early exit if possible
 IF ( 0.5_q*(b0-a0) <= tol*ABS(b0) ) THEN
   x = b0
   RETURN
 ENDIF
 
 a  = a0
 b  = b0
 fa = f(a, x0)
 fb = f(b, x0)
 
 !IF ( ANY( ((fa>0._q).AND.(fb>0._q)) .OR. ((fa<0._q).AND.(fb<0._q)) ) ) THEN
 !  ABORT
 !ENDIF
 
 c = b; fc = fb
 
 ! a is the previous value of b
 ! b is the current guess for the root
 ! c is the contrapoint
 ! d is the next step size
 ! e is the last or second-to-last step size
 ! m is the step size to the midpoint
 
 ITER = 0
 DO
   ! Ensure that c is on the opposite side of the zero from b
   WHERE ( (fb > 0._q) .EQV. (fc > 0._q) )
      c = a; fc = fa
      d = b - a
      e = d
   ENDWHERE
   ! Ensure that b is the best result so far 
   WHERE ( ABS(fc) < ABS(fb) )
      a = b; fa = fb
      b = c; fb = fc
      c = a; fc = fa
   ENDWHERE
   
   ! Convergence test and possible exit
   tol_abs = tol*ABS(b)
   m = 0.5_q * (c - b)
   IF ( ALL( ABS(m) <= tol_abs .OR. fb == 0._q ) ) THEN
      EXIT
   ENDIF
   
   ITER = ITER + 1
   IF (ITER > 100) THEN
      IF (IU0>=0) WRITE(IU0,*) 'ERROR: fzero_vec: failed to converge in 100 iterations'
      STOP
   ENDIF
   
   ! Choose bisection or interpolation
   WHERE ( ABS(e) >= tol_abs .AND. ABS(fa) > ABS(fb) )
      s1 = fb/fa
      ! linear interpolation
      WHERE ( a == c )
         p = (c - b)*s1
         r = 1._q - s1
      ! inverse quadratic interpolation
      ELSEWHERE
         s2 = fa/fc
         s3 = fb/fc
         p = s1*((c - b)*s2*(s2 - s3) - (b - a)*(s3 - 1._q))
         r = (s2 - 1._q)*(s3 - 1._q)*(s1 - 1._q)
      ENDWHERE
      WHERE ( 0._q < p )
         r = -r
      ELSEWHERE
         p = -p
      ENDWHERE
      ! accept interpolation
      ! sgn(d) = sgn(m)
      ! |d| < 1.5*|m| - 0.5*tol_abs
      ! |d| < 0.5*|e|
      WHERE ( 2._q*p < MIN(3._q*m*r - tol_abs*ABS(r), ABS(e*r)) )
         e = d
         d = p/r
      ! interpolation failed, use bisection
      ELSEWHERE
         d = m
         e = d
      ENDWHERE
   ! bounds decreasing too slowly or f not decreasing, use bisection
   ELSEWHERE
      d = m
      e = d
   ENDWHERE
   
   ! Next point
   a = b; fa = fb
   b = b + SIGN( MAX(ABS(d), tol_abs), m )
   fb = f(b, x0)
   
 ENDDO
 
 x = b
 
END SUBROUTINE fzero_vec


!>======================SUBROUTINE SHAPE_FUNC=========================
!!
!! calculates the cavity shape function in real and recip space
!!
!!====================================================================
SUBROUTINE SHAPE_FUNC(X, S)

 REAL(q), INTENT(IN)  :: X(:)
 REAL(q), INTENT(OUT) :: S(:)
 
 S = ERFC(X/(SQRT(2._q)*sigma_k))/2._q
 
END SUBROUTINE SHAPE_FUNC


!>======================SUBROUTINE SHAPE_POT=========================
!!
!! calculates the potential due to changes in the cavity shape
!! V contains nlambda at start of function, and the potential at exit
!!
!!====================================================================
SUBROUTINE SHAPE_POT(X, V, X_old, TAG)

 REAL(q), INTENT(IN)    :: X(:), X_old(:)
 REAL(q), INTENT(INOUT) :: V(:)
 
 CHARACTER(*), OPTIONAL, INTENT(IN) :: TAG
 
 OPTIONAL :: X_old
 
 REAL(q) :: Sprime(DIMREAL(GRIDC%MPLWV))
 
 LOGICAL :: ldebug
 REAL(q) :: dE1, dE2
 REAL(q), ALLOCATABLE :: S(:), S_old(:)
 
 ldebug = PRESENT(X_old)
 
 IF (ldebug) ALLOCATE(S(DIMREAL(GRIDC%MPLWV)), S_old(DIMREAL(GRIDC%MPLWV)))
 
 ! compute S' in Sprime
 Sprime = -EXP(-X**2/(2._q*sigma_k**2))/(SQRT(TPI)*sigma_k)
 
 IF (ldebug) THEN
   S_old = ERFC(X_old/(SQRT(2._q)*sigma_k))/2._q
   S = ERFC(X/(SQRT(2._q)*sigma_k))/2._q
   dE1 = SUM_GRID_RL((S-S_old)*V) * LATT_CUR%OMEGA/GRIDC%NPLWV
   dE2 = SUM_GRID_RL((X-X_old)*Sprime*V) * LATT_CUR%OMEGA/GRIDC%NPLWV
 ENDIF
 
 IF (ldebug .AND. LWRITE .AND. PRESENT(TAG)) THEN
   CALL WRITE_TO_FILE_RL('S'//TAG,S)
   CALL WRITE_TO_FILE_RL('S_old'//TAG,S_old)
   CALL WRITE_TO_FILE_RL('X'//TAG,X)
   CALL WRITE_TO_FILE_RL('X_old'//TAG,X_old)
   CALL WRITE_TO_FILE_RL('Sprime'//TAG,Sprime)
 ENDIF
 
 ! compute V = Sprime*lambda
 V = Sprime*V
 ! scale by exp(-x)
 V = V*EXP(-X)
 
 ifdsol IF(IU0>=0 .AND. ldebug) WRITE(IU0,'(2E15.6)') dE1, dE2
 
 
END SUBROUTINE SHAPE_POT


!>=====================SUBROUTINE SET_SOLV_FUNC===============================
!!
!! calculates the solvent, dielectric, and ionic shape functions in rec space
!! these are stored in the module arrays w_solv, w_diel, and w_ion
!!
!!============================================================================

SUBROUTINE SET_SOLV_FUNC
 
 IF (r_solv > 0._q) THEN
   ALLOCATE(w_solv(GRIDC%MPLWV))
   CALL F_exp(w_solv, r_solv, a_k/sigma_k)
 ELSE
   ALLOCATE(w_solv(0))
 ENDIF
 
 IF (r_cav > 0._q) THEN
   ALLOCATE(w_cav(GRIDC%MPLWV))
   CALL F_exp(w_cav, r_cav, a_k/sigma_k)
 ELSE
   ALLOCATE(w_cav(0))
 ENDIF
 
 IF (r_ion > 0._q) THEN
   ALLOCATE(w_ion(GRIDC%MPLWV))
   CALL F_exp(w_ion, r_ion, a_k/sigma_k)
 ELSE
   ALLOCATE(w_ion(0))
 ENDIF
 
 IF (r_diel > 0._q) THEN
   ALLOCATE(w_diel(GRIDC%MPLWV))
   CALL F_exp(w_diel, r_diel, a_k/sigma_k)
 ELSE
   ALLOCATE(w_diel(0))
 ENDIF
 
 !ALLOCATE(w(GRIDC%MPLWV))
 !CALL F_gaussian(w, a_k)
 
 ALLOCATE(w_b(GRIDC%MPLWV))
 IF (r_b > 0._q) THEN
   IF (I_nloc_sol == 1) THEN
      CALL F_gaussian(w_b, r_b)
   ELSEIF (I_nloc_sol == 2) THEN
      CALL F_soft_sphere(w_b, r_b, 4._q*a_k)
   ENDIF
 ELSE
   IF (I_nloc_sol == 1) THEN
      CALL F_gaussian(w_b, a_k)
   ELSEIF (I_nloc_sol == 2) THEN
      w_b = 1._q
   ENDIF
 ENDIF

END SUBROUTINE SET_SOLV_FUNC


!******************** SUBROUTINE NLPCM_READER *************************
!
! Reads in the solvation model parameters
!
!**********************************************************************
SUBROUTINE NLPCM_READER(IO)

 USE vaspxml

 TYPE (in_struct), INTENT(IN) :: IO
 
 INTEGER      :: IDUM, N, IERR, IU5
 REAL(q)      :: RDUM
 COMPLEX(q)   :: CDUM
 LOGICAL      :: LOPEN,LDUM
 CHARACTER(1) :: CHARAC
 
 IU0 = IO%IU0
 IU5 = IO%IU5

 LOPEN = .FALSE.

 OPEN(UNIT=IU5,FILE='INCAR',STATUS='OLD')      

 !-----------------------------------------------------------------
 ! LNLDIEL, nonlinear dielectric effect on/off
 !-----------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'LNLDIEL','=','#',';','L',IDUM,RDUM,CDUM,LNLDIEL,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading item ''LNLDIEL'' from file INCAR.'
   STOP
 ENDIF
 CALL XML_INCAR('LNLDIEL','L',IDUM,RDUM,CDUM,LNLDIEL,CHARAC,N)

 !-----------------------------------------------------------------
 ! LNLION, nonlinear ionic effect on/off
 !-----------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'LNLION','=','#',';','L',IDUM,RDUM,CDUM,LNLION,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading item ''LNLION'' from file INCAR.'
   STOP
 ENDIF
 CALL XML_INCAR('LNLION','L',IDUM,RDUM,CDUM,LNLION,CHARAC,N)
 
 !---------------------------------------------------------------
 !sigma_k
 !--------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'SIGMA_K','=','#',';','F',IDUM,sigma_k,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''SIGMA_K''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('SIGMA_K','F',IDUM,sigma_k,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !a_k
 !--------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'A_K','=','#',';','F',IDUM,a_k,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''A_K''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('A_K','F',IDUM,a_k,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !nc_k
 !---------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'NC_K','=','#',';','F',IDUM,nc_k,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''NC_K''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('NC_K','F',IDUM,nc_k,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !r_solv
 !---------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'R_SOLV','=','#',';','F',IDUM,r_solv,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''R_SOLV''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('R_SOLV','F',IDUM,r_solv,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !r_cav
 !---------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'R_CAV','=','#',';','F',IDUM,r_cav,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''R_CAV''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('R_CAV','F',IDUM,r_cav,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !r_diel
 !---------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'R_DIEL','=','#',';','F',IDUM,r_diel,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''R_DIEL''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('R_DIEL','F',IDUM,r_diel,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !r_ion
 !---------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'R_ION','=','#',';','F',IDUM,r_ion,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''R_ION''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('R_ION','F',IDUM,r_ion,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !r_b
 !---------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'R_B','=','#',';','F',IDUM,r_b,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''R_B''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('R_B','F',IDUM,r_b,CDUM,LDUM,CHARAC,N)
 
 !-----------------------------------------------------------------
 ! I_nloc_sol
 !-----------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'I_nloc_sol','=','#',';','I',I_nloc_sol,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading item ''I_nloc_sol'' from file INCAR.'
   STOP
 ENDIF
 CALL XML_INCAR('I_nloc_sol','I',I_nloc_sol,RDUM,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !tau
 !---------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'TAU','=','#',';','F',IDUM,tau,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''TAU''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('TAU','F',IDUM,tau,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !eb_k
 !---------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'EB_K','=','#',';','F',IDUM,eb_k,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''EB_K''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('EB_K','F',IDUM,eb_k,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !lambda_d_k = debye screening length in angstrom
 !--------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'LAMBDA_D_K','=','#',';','F',IDUM,lambda_d_k,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''LAMBDA_D_K''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 
 !---------------------------------------------------------------
 !SolTemp
 !--------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'SolTemp','=','#',';','F',IDUM,SolTemp,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''SolTemp''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('SolTemp','F',IDUM,SolTemp,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !d_ion
 !--------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'d_ion','=','#',';','F',IDUM,d_ion,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''d_ion''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('d_ion','F',IDUM,d_ion,CDUM,LDUM,CHARAC,N)     
 
 !---------------------------------------------------------------
 !c_molar
 !--------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'c_molar','=','#',';','F',IDUM,c_molar,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''c_molar''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('c_molar','F',IDUM,c_molar,CDUM,LDUM,CHARAC,N)        
 
 !---------------------------------------------------------------
 !ZION
 !--------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'ZION','=','#',';','F',IDUM,ZION,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''ZION''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('ZION','F',IDUM,ZION,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !n_mol
 !--------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'n_mol','=','#',';','F',IDUM,n_mol,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''n_mol''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('n_mol','F',IDUM,n_mol,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !p_mol
 !--------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'p_mol','=','#',';','F',IDUM,p_mol,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''p_mol''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('p_mol','F',IDUM,p_mol,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !epsilon_inf
 !--------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'epsilon_inf','=','#',';','F',IDUM,epsilon_inf,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''epsilon_inf''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('epsilon_inf','F',IDUM,epsilon_inf,CDUM,LDUM,CHARAC,N)
 
 !---------------------------------------------------------------
 !LNLTEST
 !--------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'LNLTEST','=','#',';','L',IDUM,RDUM,CDUM,LNLTEST,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading ''LNLTEST''from INCAR.'
   sol_io ' Error code = ', IERR
   STOP
 ENDIF
 CALL XML_INCAR('LNLTEST','L',IDUM,RDUM,CDUM,LNLTEST,CHARAC,N)
 
 
 CLOSE(IU5)
 
END SUBROUTINE NLPCM_READER


!******************** SUBROUTINE NLPCM_INIT ****************************
!      
! initialize NLPCM module
!
!***********************************************************************
SUBROUTINE NLPCM_INIT
 
 REAL(q) :: c_ion_b, alpha_lin, f_loc, alpha0_diel
 
 REAL(q), PARAMETER :: MOLAR=6.022e-4_q, lin_fact=1e10_q
 
 ! intialize real and rec space grid data
 CALL SETGRID
 
 IF (LNLTEST) THEN
   SolTemp = SolTemp*lin_fact**2
   p_mol   = p_mol*lin_fact
   ZION    = ZION*lin_fact
 ENDIF
 
 invBETA = BOLKEV*SolTemp
 !p_mol = sqrt(3._q*invBETA*(eb_k-epsilon_inf)/(n_mol*EDEPS)) ! for testing
 PBETA = p_mol/invBETA
 ZBETA = ZION/invBETA
 
 dsol_io
 dsol_io 'invBETA = ', invBETA
 dsol_io 'PBETA   = ', PBETA
 dsol_io 'ZBETA   = ', ZBETA
 dsol_io

 alpha0_rot = EDEPS*invBETA*PBETA**2/3._q
 alpha_pol = alpha0_rot/(eb_k-epsilon_inf)*(epsilon_inf-1._q)
 invalpha_sic = ((eb_k-epsilon_inf)/alpha0_rot - n_mol)/(eb_k-1._q)
 alpha_lin = alpha0_rot + alpha_pol
 f_loc = 1._q/(1._q - alpha_lin*invalpha_sic)
 alpha0_diel = alpha_lin*f_loc
 
 dsol_io 'alpha0_rot   = ', alpha0_rot
 dsol_io 'alpha_pol    = ', alpha_pol
 dsol_io 'invalpha_sic = ', invalpha_sic
 dsol_io 'alpha_lin    = ', alpha_lin
 dsol_io 'f_loc        = ', f_loc
 dsol_io 'alpha0_diel  = ', alpha0_diel
 dsol_io

 ! if alpha0_rot is too small, it will result in significant floating point error
 IF (n_mol*alpha0_rot <= (eb_k-epsilon_inf)*1e-10_q .OR. invalpha_sic < 0._q) THEN
   sol_io 'ERROR: Invalid parameters for dielectric screening model'
   STOP
 ENDIF
 
 alpha0_ion = EDEPS*invBETA*ZBETA**2
 IF (lambda_d_k > 0._q) THEN
   c_ion_b = eb_k/(lambda_d_k**2 * alpha0_ion)
 ELSE
   c_ion_b = 2._q*c_molar*MOLAR
 ENDIF
 
 LION = (c_ion_b > 0._q)
 LNLION = (LNLION .AND. LION)
 IF (d_ion < 0._q) d_ion = 2._q**(5._q/6._q)*R_ION
 IF (LNLION .AND. d_ion > 0._q) THEN
   n_max = 1._q/d_ion**3
   theta_b = c_ion_b/n_max
   alpha0_ion = theta_b*alpha0_ion
   IF (theta_b < 1e-8_q .OR. theta_b > 1._q) THEN
      sol_io 'ERROR: Invalid C_ION_B and/or D_ION for nonlinear ionic screening model'
      STOP
   ENDIF
 ELSE
   n_max = c_ion_b
   theta_b = 0._q
 ENDIF
 
 dsol_io 'alpha0_ion = ', alpha0_ion
 dsol_io 'n_max      = ', n_max
 dsol_io 'theta_b    = ', theta_b
 dsol_io 'LION       = ', LION
 dsol_io
 
 
 ! initialize w_solv, w_ion, w_diel, w
 CALL SET_SOLV_FUNC
 
 
 ALLOCATE(phi(GRIDC%MPLWV), phi_solv(GRIDC%MPLWV), E_loc(DIMREAL(GRIDC%MPLWV),4), P(DIMREAL(GRIDC%MPLWV),3))
 ALLOCATE(S_cav(DIMREAL(GRIDC%MPLWV)), S_ion(DIMREAL(GRIDC%MPLWV)), S_diel(DIMREAL(GRIDC%MPLWV)))
 ALLOCATE(Xcav(DIMREAL(GRIDC%MPLWV),5))
 ALLOCATE(lambda_diel(DIMREAL(GRIDC%MPLWV)), lambda_ion(DIMREAL(GRIDC%MPLWV)), lambda_cav(DIMREAL(GRIDC%MPLWV)))
 ALLOCATE(n_b(GRIDC%MPLWV), n_ion(GRIDC%MPLWV))
 
 phi = 0._q
 phi_solv = 0._q

END SUBROUTINE NLPCM_INIT


!******************** SUBROUTINE NLPCM_WRITER **************************
!      
! writes the solvation model parameters to the OUTCAR file
!
!***********************************************************************
SUBROUTINE NLPCM_WRITER(IO)
 
 TYPE (in_struct), INTENT(IN) :: IO
 
 IF (IO%IU6>=0) THEN
   WRITE(IO%IU6,100) LNLDIEL, LNLION, LNLTEST, I_nloc_sol, &
               &  nc_k, tau, sigma_k, a_k, r_solv, r_cav, r_diel, r_ion, r_b, &
               &  SolTemp, eb_k, epsilon_inf, n_mol, p_mol, c_molar, lambda_d_k, ZION, d_ion
 ENDIF

100 FORMAT( &
            ' Nonlinear solvation model:'/ &
            '   LNLDIEL     =',L6, '         switch on non-linear dielectric mode' / &
            '   LNLION      =',L6, '         switch on non-linear ionic mode' / &
            '   LNLTEST     =',L6, '         test non-linear model in linear regime' / &
            '   I_NLOC_SOL  =',I6, '         type of nonlocal dielectric response' / &
            '   NC_K        =',F10.6,  '     cutoff charge density for cavities (1/Angstrom**3)' /&
            '   TAU         =',F10.6,  '     cavity formation energy (eV/Angstrom**2)' /&
            '   SIGMA_K     =',F10.6,  '     cavity width' /&
            '   A_K         =',F10.6,  '     smoothing length (Angstrom)' /&
            '   R_SOLV      =',F10.6,  '     solvent radius (Angstrom)' / &
            '   R_CAV       =',F10.6,  '     cavity radius (Angstrom)' / &
            '   R_DIEL      =',F10.6,  '     dielectric radius (Angstrom)' / &
            '   R_ION       =',F10.6,  '     ionic radius (Angstrom)' / &
            '   R_B         =',F10.6,  '     bound charge radius (Angstrom)' / &
            '   SOLTEMP     =',F10.6,  '     solvent temperature (K)' / &
            '   EB_K        =',F10.6,  '     relative permittivity of the bulk solvent' /&
            '   EPSILON_INF =',F10.6,  '     optical relative permittivity of bulk solvent' /&         
            '   N_MOL       =',F10.6,  '     molar density of solvent (1/Angstrom**3)' / &
            '   P_MOL       =',F10.6,  '     solvent dipole moment (e*Angstrom)' / &
            '   C_MOLAR     =',F10.6,  '     electrolyte concentration (mol/L)' / &
            '   LAMBDA_D_K  =',F10.6,  '     Debye length (Angstrom)' /&
            '   ZION        =',F10.6,  '     charge of ions in electrolyte (e)' / &
            '   D_ION       =',F10.6,  '     minimum ion packing diameter (Angstrom)' / )


END SUBROUTINE NLPCM_WRITER



!******************** SUBROUTINE XML_WRITE_NLPCM ***********************
!      
!
! writes the solvation model parameters to vasprun.xml
!
!
!***********************************************************************
SUBROUTINE XML_WRITE_NLPCM

 USE vaspxml
 
 INTEGER      :: IDUM
 REAL(q)      :: RDUM
 COMPLEX(q)   :: CDUM
 LOGICAL      :: LDUM
 CHARACTER(1) :: CHARAC
 
 CALL XML_TAG("separator","solvation model")
 
 CALL XML_INCAR('LNLDIEL','L',IDUM,RDUM,CDUM,LNLDIEL,CHARAC,1)
 CALL XML_INCAR('LNLION','L',IDUM,RDUM,CDUM,LNLION,CHARAC,1)
 CALL XML_INCAR('LNLTEST','L',IDUM,RDUM,CDUM,LNLTEST,CHARAC,1)
 CALL XML_INCAR('I_NLOC_SOL','I',I_nloc_sol,RDUM,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('NC_K','F',IDUM,nc_k,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('TAU','F',IDUM,tau,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('SIGMA_K','F',IDUM,sigma_k,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('A_K','F',IDUM,a_k,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('R_SOLV','F',IDUM,r_solv,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('R_CAV','F',IDUM,r_cav,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('R_DIEL','F',IDUM,r_diel,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('R_ION','F',IDUM,r_ion,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('R_B','F',IDUM,r_b,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('SOLTEMP','F',IDUM,SolTemp,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('EB_K','F',IDUM,eb_k,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('EPSILON_INF','F',IDUM,epsilon_inf,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('N_MOL','F',IDUM,n_mol,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('P_MOL','F',IDUM,p_mol,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('C_MOLAR','F',IDUM,c_molar,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('LAMBDA_D_K','F',IDUM,lambda_d_k,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('ZION','F',IDUM,ZION,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('D_ION','F',IDUM,d_ion,CDUM,LDUM,CHARAC,1)
 
 CALL XML_CLOSE_TAG

END SUBROUTINE XML_WRITE_NLPCM


!>=======================SUBROUTINE SET_SOL_WRITE=========================
!!
!! write ionic and bound charges to files
!!
!!====================================================================
SUBROUTINE SET_SOL_WRITE_NLPCM
 LWRITE = .TRUE.
END SUBROUTINE SET_SOL_WRITE_NLPCM


END MODULE solvation_nlpcm









!>==========================MODULE solvation===============================
!!
!! This module glues together the solvation engine module POT_LPCM_K
!!with the rest of the VASP code
!!
!! Contains subroutines that perform solvation I/O, pre-process 
!!the charge densities, call the main subroutine, Vcorrection_lpcm,  
!!from the module POT_LPCM_K and compute the force corrections
!!
!!=====================================================================
MODULE solvation
 
 USE prec
 USE base
 USE constant
 USE mpimy
 USE mgrid,   ONLY : grid_3d
 USE lattice, ONLY : latt
 USE poscar,  ONLY : type_info
 USE solvation_grid
 USE solvation_lpcm
 USE solvation_nlpcm
 
 IMPLICIT NONE

#ifdef jdftx
 EXTERNAL initjdftx, getionsigma, minimizefluid
#endif

 PRIVATE
 
 PUBLIC :: SOL_READER, SOL_INIT, SOL_WRITER, XML_WRITE_SOL, LDO_SOL
 PUBLIC :: SOL_VCORRECTION, SET_RMS_SCF, SET_SOL_WRITE, update_NELECT, check_EFERMI, L_const_pot
 PUBLIC :: Ecorr_sol, V_corr, n_solv
 
 SAVE
 
 INTEGER :: PK_counter = 0
 
 INTEGER :: IU0, IU6
 
 LOGICAL :: LSOL = .FALSE.
 LOGICAL :: LSOL_SCF = .TRUE.
 LOGICAL :: LDO_SOL = .FALSE.
 INTEGER :: ISOL = 1
 INTEGER :: INISOL
 
 REAL(q) :: rms_SCF = 0._q
 REAL(q) :: q_sol = 0._q
 REAL(q) :: ZVSUM = 0._q
 REAL(q) :: EFERMI_ref = 0._q
 REAL(q) :: EFERMI_tol = 0._q
 REAL(q) :: capacitance_init = 1._q
 
 LOGICAL :: LWRITE = .FALSE.
 
 !solvent charge density
 COMPLEX(q), ALLOCATABLE :: n_solv(:)
 !correction to KS potential
 COMPLEX(q), ALLOCATABLE :: V_corr(:)
 !correction to free energy
 REAL(q) :: Ecorr_sol = 0._q
 
 
CONTAINS


!******************** SUBROUTINE SOL_VCORRECTION ***********************
!
! Computes the potential, energy and force corrections due to solvation
!
!***********************************************************************
SUBROUTINE SOL_VCORRECTION(CHTOT, DENCOR, CVHAR)
 
 USE ini
 
 COMPLEX(q),      INTENT(IN) :: CHTOT(:,:)
 RGRID,           INTENT(IN) :: DENCOR(:)
 COMPLEX(q),      INTENT(IN) :: CVHAR(:)
 
 COMPLEX(q) :: n_val(:), n_e(:), phi_sol(:), CWORK(:)
 REAL(q)    :: RWORK(:)
 REAL(q)    :: Ecorr_band, A_corr, tol, tmp
 
 REAL(q), PARAMETER :: tol_init = 1.e-3_q
 
 ALLOCATABLE :: n_val, n_e, phi_sol, CWORK, RWORK
 
 
 ! early exit if possible
 IF (.NOT.LDO_SOL) RETURN
 
 
!-----------------------------------------------------------------------------------------
! if INISOL=1, accurately solve the PCM equation in the first SCF step using tol_init
!   this is used by default when the initial charge density
!   is read from a file or calculated from the wavefunctions (ICHARG=0,1)
! if INISOL=2, turn of solvation during the non-self-consistent phase (first NELMDL steps)
!   this is used by default when the initial charge density
!   is taken from the superposition of atomic charge densities (ICHARG=2,3)
!-----------------------------------------------------------------------------------------
 IF (PK_COUNTER == 0)  THEN
   PK_COUNTER = PK_COUNTER + 1
   IF (INISOL == 2) THEN
      dsol_io ' Solvation turned off during non-self-consistent phase (first NELMDL steps)'
      RETURN
   ENDIF
 ENDIF
 
 CALL SEPERATOR_TIMING(IU6)
 CALL START_TIMING("SOL")
 
 
 dsol_io 'Begin Solvation calculations'
 
 ALLOCATE(n_val(GRIDC%MPLWV), n_e(GRIDC%MPLWV), phi_sol(GRIDC%MPLWV))
 ALLOCATE(CWORK(GRIDC%MPLWV), RWORK(DIMREAL(GRIDC%MPLWV)))
 
 
 ! copy charge to n_val
 ! CHTOT(:,1) contains the total valence charge density
 ! copy using RC_ADD to ensure that unused elements are properly zeroed
 n_val = 0
 CALL RC_ADD(CHTOT(1,1), 1.0_q, CHTOT(1,1), 0.0_q, n_val(1), GRIDC)
 
 ! compute n_e = n_val + n_c
 RWORK = 0
 CALL RL_ADD(DENCOR(1), 1.0_q, DENCOR(1), 0.0_q, RWORK(1), GRIDC)
 CALL FFT3D_RL2RC(RWORK,n_e)
 n_e = n_e + n_val
 !CALL AVG_GRID_RC(n_e, CWORK, 3)
 !n_e = CWORK
 
 ! compute total charge, q_sol
 CALL GET_G0(n_val, q_sol)
 q_sol = q_sol - ZVSUM
 
 ! copy solute Hartree potential to phi_sol
 phi_sol = 0
 CALL RC_ADD(CVHAR(1), 1.0_q, CVHAR(1), 0.0_q, phi_sol(1), GRIDC)
 !CALL AVG_GRID_RC(phi_sol, CWORK, 3)
 !phi_sol = CWORK
 
 
 !------------------------------------------------------------------------------
 ! set tolerance for outer loop, tol is in units of charge
 !------------------------------------------------------------------------------
 IF (rms_SCF > 0._q) THEN
   ! set tolerance for LPB to 1/10 the current rms of the SCF residual
   tol = rms_SCF/10._q
   dsol_io ' rms_SCF = ', rms_SCF, 'tol = ', tol
 ELSE
   ! for first SCF step or constant charge calculations, set LPB tolerance to tol_init
   tol = tol_init
   dsol_io' Setting tol to tol_init = ', tol_init
 ENDIF
 
 CALL STOP_TIMING("SOL",IU6,"INIT")
 
 
 !-------------------------------------------------------------
 !              BEGIN Solvation minimization
 !           Compute Vcorr, the correction potential
 !           and Atot, the total solvation energy
 !-------------------------------------------------------------   
 dsol_io 'Begin minimization'
 SELECT CASE(ISOL)
   CASE (1)
      ! run linear polarizable continuum model
      CALL LPCM_solver(phi_sol, n_e, q_sol, tol, V_corr, n_solv, A_corr)
   CASE (2)
      ! run nonlinear/nonlocal polarizable continuum model
      CALL NLPCM_solver(phi_sol, n_e, q_sol, tol, V_corr, n_solv, A_corr)
   CASE (3)
      ! run JDFT based polarizable continuum model
      IF (IU0>=0) WRITE(IU0,*) 'ERROR: SOL_VCORRECTION: invalid value for ISOL'
      STOP
   CASE (4)
      ! run externally using JDFTx
      IF (IU0>=0) WRITE(IU0,*) 'ERROR: SOL_VCORRECTION: invalid value for ISOL'
      STOP
   CASE DEFAULT
      IF (IU0>=0) WRITE(IU0,*) 'ERROR: SOL_VCORRECTION: invalid value for ISOL'
      STOP
 END SELECT
 dsol_io 'End minimization'
 !-------------------------------------------------------------
 !              END Solvation minimization
 !-------------------------------------------------------------    
 
 
 !-------------------------------------------------------------
 !net correction to the total energy
 !we need to subtract off Ecorr_band = n_val*Vcorr
 ! since this will end up in the band energy
 !-------------------------------------------------------------
 Ecorr_band = D_PROD_RC(n_val,V_corr)
 Ecorr_sol = A_corr - Ecorr_band - q_sol*EFERMI_ref
 
 dsol_io 'Ecorr=', Ecorr_sol
 dsol_io 'A_corr=', A_corr
 dsol_io 'Ecorr_band=', Ecorr_band
 
 dsol_io 'End Solvation calculations'
 
 CALL STOP_TIMING("SOL",IU6,"FINALIZE")
 CALL SEPERATOR_TIMING(IU6)
 
 
END SUBROUTINE SOL_VCORRECTION


!>=======================SUBROUTINE SET_RMS_SCF=======================
!!
!! set module variable rms_SCF based on rms of charge density
!! used to determine convergence criteria for solvent charge densities
!!
!!====================================================================
SUBROUTINE SET_RMS_SCF(CHTOT,CHTOTL,LRESET)
 
 COMPLEX(q), INTENT(IN) :: CHTOT(:,:), CHTOTL(:,:)
 LOGICAL,    INTENT(IN) :: LRESET
 
 COMPLEX(q), ALLOCATABLE :: resid(:), CWORK(:)
 REAL(q) :: rms_SCF_new
 
 IF (.NOT.LDO_SOL) RETURN
 
 ! reintialize real and rec space grid data after ionic step
 !  in case unit cell changed
 IF (LRESET) CALL SETGRID
 
 ALLOCATE(resid(GRIDC%MPLWV), CWORK(GRIDC%MPLWV))
 
 resid = 0.
 CALL RC_ADD(CHTOT(1,1), 1.0_q, CHTOTL(1,1), -1.0_q, resid, GRIDC)
 
 CALL L0_inv_op(resid, CWORK)
 rms_SCF_new = sqrt( D_PROD_RC(CWORK,CWORK) )
 
 IF (rms_SCF > 0._q) THEN
   rms_SCF = min(rms_SCF,rms_SCF_new)
 ELSE
   rms_SCF = rms_SCF_new
 ENDIF

END SUBROUTINE SET_RMS_SCF


!>=======================SUBROUTINE update_NELECT=====================
!!
!! update number of electrons for constant potential calculation
!!
!!====================================================================
SUBROUTINE update_NELECT(NELECT,EFERMI)
 
 REAL(q), INTENT(INOUT) :: NELECT
 REAL(q), INTENT(IN)    :: EFERMI
 
 REAL(q), SAVE :: EFERMI_LAST, NELECT_LAST, capacitance
 LOGICAL, SAVE :: LFIRST = .TRUE.
 
 REAL(q) :: qsol_new
 
 IF (.NOT.L_const_pot()) RETURN
 
 IF (LFIRST) THEN
   capacitance = capacitance_init
   LFIRST = .FALSE.
 ELSEIF (ABS(EFERMI-EFERMI_LAST) > 0.1_q) THEN
   capacitance = (NELECT - NELECT_LAST)/(EFERMI - EFERMI_LAST)
 ENDIF
 
 EFERMI_LAST = EFERMI
 NELECT_LAST = NELECT
 qsol_new = q_sol  - capacitance*(EFERMI - EFERMI_ref)
 IF (qion_max > 0._q .AND. ABS(qsol_new) > qion_max*0.9_q) THEN
   IF (ABS(q_sol) >= qion_max*0.9_q) THEN
      IF (IU0>=0) THEN
         WRITE(IU0,*) 'ERROR: update_NELECT: your ''vacuum'' region is not large enough to accommodate the necessary counter charge'
         WRITE(IU0,*) 'The counter charge is exceeding 90% of the saturation value'
         WRITE(IU0,*) 'Maximum charge: ', qion_max
         WRITE(IU0,*) 'Surface/solute charge: ', q_sol
      ENDIF
      STOP
   ENDIF
   qsol_new = SIGN(MIN(ABS(qsol_new), qion_max*0.9_q), qsol_new)
 ENDIF
 NELECT = ZVSUM + qsol_new
 
 IF (IU6>=0) THEN
   WRITE(IU6,*)
   WRITE(IU6,*) 'Adjusting number of electrons for constant potential calculation'
   WRITE(IU6,*)
   WRITE(IU6,'(3X,A,F18.8,6X,A,F18.8)') 'EFERMI = ', EFERMI,      'EFERMI_ref  = ', EFERMI_ref
   WRITE(IU6,'(3X,A,F18.8,6X,A,F18.8)') 'q_sol  = ', q_sol,       '-q*mu      =  ', -q_sol*EFERMI_ref
   WRITE(IU6,'(3X,A,F18.8,6X,A,F18.8)') 'NELECT = ', NELECT_LAST, 'new NELECT =  ', NELECT
   WRITE(IU6,'(36X,A,F18.8)')                                     'capacitance = ', capacitance
   WRITE(IU6,*)
 ENDIF

END SUBROUTINE update_NELECT


!>=======================FUNCTION check_EFERMI========================
!!
!! check if Fermi level is converged for constant potential calculation
!!
!!====================================================================
PURE FUNCTION check_EFERMI(EFERMI)
 
 LOGICAL :: check_EFERMI
 REAL(q), INTENT(IN) :: EFERMI
 
 check_EFERMI = (ABS(EFERMI-EFERMI_ref) < EFERMI_tol) .OR. (EFERMI_ref >= 0.0_q)

END FUNCTION check_EFERMI


!>=======================FUNCTION L_const_pot=========================
!!
!! Constant potential calculation?
!!
!!====================================================================
PURE FUNCTION L_const_pot
 LOGICAL :: L_const_pot
 L_const_pot = (EFERMI_ref < 0._q)
END FUNCTION L_const_pot


!******************** SUBROUTINE SOL_READER ***************************
!
! Reads in the solvation model parameters
!      
!**********************************************************************
SUBROUTINE SOL_READER(IO)
 
 USE vaspxml

 TYPE (in_struct), INTENT(IN) :: IO
 
 INTEGER      :: IDUM, N, IERR, IU5
 REAL(q)      :: RDUM
 COMPLEX(q)   :: CDUM
 LOGICAL      :: LOPEN,LDUM
 CHARACTER(1) :: CHARAC
 
 IU0 = IO%IU0
 IU6 = IO%IU6
 IU5 = IO%IU5

 LOPEN = .FALSE.

 OPEN(UNIT=IU5,FILE='INCAR',STATUS='OLD')

 !-----------------------------------------------------------------
 ! LSOL, solvation on/off
 !-----------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'LSOL','=','#',';','L',IDUM,RDUM,CDUM,LSOL,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading item ''LSOL'' from file INCAR.'
   STOP
 ENDIF
 
 CALL XML_INCAR('LSOL','L',IDUM,RDUM,CDUM,LSOL,CHARAC,N)
 
 IF (.NOT.LSOL) THEN
   CLOSE(IU5)
   RETURN
 ENDIF
 
 !-----------------------------------------------------------------
 ! LSOL_SCF, self consistent solvation on/off
 !-----------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'LSOL_SCF','=','#',';','L',IDUM,RDUM,CDUM,LSOL_SCF,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading item ''LSOL_SCF'' from file INCAR.'
   STOP
 ENDIF

 !-----------------------------------------------------------------
 ! ISOL, which solvation model
 !-----------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'ISOL','=','#',';','I',ISOL,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading item ''ISOL'' from file INCAR.'
   STOP
 ENDIF
 CALL XML_INCAR('ISOL','I',ISOL,RDUM,CDUM,LDUM,CHARAC,N)
 
 !-----------------------------------------------------------------
 ! EFERMI_ref
 !-----------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'EFERMI_ref','=','#',';','F',IDUM,EFERMI_ref,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading item ''EFERMI_ref'' from file INCAR.'
   STOP
 ENDIF
 CALL XML_INCAR('EFERMI_ref','F',IDUM,EFERMI_ref,CDUM,LDUM,CHARAC,N)
 
 !-----------------------------------------------------------------
 ! EFERMI_tol
 !-----------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'EFERMI_tol','=','#',';','F',IDUM,EFERMI_tol,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading item ''EFERMI_tol'' from file INCAR.'
   STOP
 ENDIF
 CALL XML_INCAR('EFERMI_tol','F',IDUM,EFERMI_tol,CDUM,LDUM,CHARAC,N)
 
 !-----------------------------------------------------------------
 ! capacitance_init
 !-----------------------------------------------------------------
 CALL RDATAB(LOPEN,'INCAR',IU5,'capacitance_init','=','#',';','I',IDUM,capacitance_init,CDUM,LDUM,CHARAC,N,1,IERR)
 IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
   sol_io 'Error reading item ''capacitance_init'' from file INCAR.'
   STOP
 ENDIF
 CALL XML_INCAR('capacitance_init','F',IDUM,capacitance_init,CDUM,LDUM,CHARAC,N)
 
 
 CLOSE(IU5)
 
 
 SELECT CASE (ISOL)
   CASE (1)
      CALL LPCM_READER(IO)
   CASE (2)
      CALL NLPCM_READER(IO)
   CASE DEFAULT
      IF (IU0>=0) WRITE(IU0,*) 'ERROR: SOL_READER: invalid value for ISOL'
      STOP
 END SELECT


END SUBROUTINE SOL_READER


!******************** SUBROUTINE SOL_INIT ******************************
!      
! initialize solvation model
!
!***********************************************************************
SUBROUTINE SOL_INIT(GRIDC, LATT_CUR, T_INFO, P, IO, INFO)

 USE pseudo, ONLY : potcar

 TYPE (grid_3d),   TARGET, INTENT(IN) :: GRIDC
 TYPE (latt),      TARGET, INTENT(IN) :: LATT_CUR
 TYPE (type_info), TARGET, INTENT(IN) :: T_INFO
 TYPE (potcar),            INTENT(IN) :: P(T_INFO%NTYP)
 TYPE (in_struct), TARGET, INTENT(IN) :: IO
 TYPE (info_struct),       INTENT(IN) :: INFO
 
 INTEGER :: NT
 
 ALLOCATE(n_solv(GRIDC%MPLWV),V_corr(GRIDC%MPLWV))
 
 V_corr = 0.
 n_solv = 0.
 
 ! early exit if possible
 IF (.NOT.LSOL) RETURN
 
 ! turn off solvation until the end if not self consistent
 LDO_SOL = LSOL_SCF
 
 CALL INIT_SOL_GRID(GRIDC, LATT_CUR, T_INFO, IO)
 
 ! initial charge density from superposition of atoms
 ! turn solvation off during non-self-consistent phase
 IF ((INFO%INICHG == 2 .OR. INFO%INICHG == 3) .AND. LSOL_SCF) THEN
   INISOL = 2
 ! initial charge density read or computed from wavefunctions
 ! turn solvation on from the start
 ELSE
   INISOL = 1
 ENDIF
 
 ! compute total valence of system
 ZVSUM = 0._q
 DO NT = 1,T_INFO%NTYP
   ZVSUM = ZVSUM + P(NT)%ZVALF * T_INFO%NITYP(NT) * T_INFO%VCA(NT)
 ENDDO
 
 ! set convergence tolerence for fermi level in constant potential calculations
 IF (EFERMI_tol <= 0._q) EFERMI_tol = 10._q*INFO%EDIFF
 
 SELECT CASE (ISOL)
   CASE (1)
      CALL LPCM_INIT
   CASE (2)
      CALL NLPCM_INIT
   CASE DEFAULT
      IF (IU0>=0) WRITE(IU0,*) 'ERROR: SOL_INIT: invalid value for ISOL'
      STOP
 END SELECT
 
END SUBROUTINE SOL_INIT


!******************** SUBROUTINE SOL_WRITER ****************************
!      
! writes the solvation model parameters to the OUTCAR file
!
!***********************************************************************
SUBROUTINE SOL_WRITER(IO)

 TYPE (in_struct), INTENT(IN) :: IO
 
 ! early exit if possible
 IF (.NOT.LSOL) RETURN
 
 IF (IO%IU6>=0) THEN
   WRITE(IO%IU6,100) LSOL, LSOL_SCF, ISOL, EFERMI_ref, capacitance_init
 ENDIF

100 FORMAT( &
   ' Solvation model:'/ &
   '   LSOL             =',L6,    '         switch on solvation model' / &
   '   LSOL_SCF         =',L6,    '         include solvation in SCF' / &
   '   ISOL             =',I6,    '         solvation model type' / &
   '   EFERMI_ref       =',F10.6, '         reference Fermi level' / &
   '   capacitance_init =',F10.6, '         initial capacitance')

 SELECT CASE (ISOL)
   CASE (1)
      CALL LPCM_WRITER(IO)
   CASE (2)
      CALL NLPCM_WRITER(IO)
   CASE DEFAULT
      IF (IU0>=0) WRITE(IU0,*) 'ERROR: SOL_WRITER: invalid value for ISOL'
      STOP
 END SELECT

END SUBROUTINE SOL_WRITER


!******************** SUBROUTINE XML_WRITE_SOL *************************
!      
! writes the solvation model parameters to vasprun.xml
!
!***********************************************************************
SUBROUTINE XML_WRITE_SOL

 USE vaspxml
 
 INTEGER      :: IDUM
 REAL(q)      :: RDUM
 COMPLEX(q)   :: CDUM
 LOGICAL      :: LDUM
 CHARACTER(1) :: CHARAC
 
 CALL XML_TAG("separator","solvation model")
 
 CALL XML_INCAR('LSOL','L',IDUM,RDUM,CDUM,LSOL,CHARAC,1)
 
 IF (.NOT.LSOL) THEN
   CALL XML_CLOSE_TAG
   RETURN
 ENDIF
 
 CALL XML_INCAR('LSOL_SCF','L',IDUM,RDUM,CDUM,LSOL_SCF,CHARAC,1)
 CALL XML_INCAR('ISOL','I',ISOL,RDUM,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('EFERMI_ref','F',IDUM,EFERMI_ref,CDUM,LDUM,CHARAC,1)
 CALL XML_INCAR('capacitance_init','F',IDUM,capacitance_init,CDUM,LDUM,CHARAC,1)
 
 CALL XML_CLOSE_TAG
 
 SELECT CASE (ISOL)
   CASE (1)
      CALL XML_WRITE_LPCM
   CASE (2)
      CALL XML_WRITE_NLPCM
   CASE DEFAULT
      IF (IU0>=0) WRITE(IU0,*) 'ERROR: XML_WRITE_SOL: invalid value for ISOL'
      STOP
 END SELECT

END SUBROUTINE XML_WRITE_SOL


!******************** SUBROUTINE SET_SOL_WRITE *************************
!      
! write ionic and bound charges to files
!
!***********************************************************************
SUBROUTINE SET_SOL_WRITE
 
 LDO_SOL = LSOL
 LWRITE = .TRUE.
 
 SELECT CASE (ISOL)
   CASE (1)
      CALL SET_SOL_WRITE_LPCM
   CASE (2)
      CALL SET_SOL_WRITE_NLPCM
   CASE DEFAULT
      IF (IU0>=0) WRITE(IU0,*) 'ERROR: SET_SOL_WRITE: invalid value for ISOL'
      STOP
 END SELECT

END SUBROUTINE SET_SOL_WRITE





END MODULE solvation















