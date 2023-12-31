!> \brief \b SGTT02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGTT02( TRANS, N, NRHS, DL, D, DU, X, LDX, B, LDB,
!                          RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            LDB, LDX, N, NRHS
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               B( LDB, * ), D( * ), DL( * ), DU( * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGTT02 computes the residual for the solution to a tridiagonal
!> system of equations:
!>    RESID = norm(B - op(A)*X) / (norm(op(A)) * norm(X) * EPS),
!> where EPS is the machine epsilon.
!> The norm used is the 1-norm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER
!>          Specifies the form of the residual.
!>          = 'N':  B - A    * X  (No transpose)
!>          = 'T':  B - A**T * X  (Transpose)
!>          = 'C':  B - A**H * X  (Conjugate transpose = Transpose)
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrices B and X.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] DL
!> \verbatim
!>          DL is REAL array, dimension (N-1)
!>          The (n-1) sub-diagonal elements of A.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The diagonal elements of A.
!> \endverbatim
!>
!> \param[in] DU
!> \verbatim
!>          DU is REAL array, dimension (N-1)
!>          The (n-1) super-diagonal elements of A.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is REAL array, dimension (LDX,NRHS)
!>          The computed solution vectors X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,NRHS)
!>          On entry, the right hand side vectors for the system of
!>          linear equations.
!>          On exit, B is overwritten with the difference B - op(A)*X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
!>          norm(B - op(A)*X) / (norm(op(A)) * norm(X) * EPS)
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
!
!> \ingroup single_lin
!
!  =====================================================================
   SUBROUTINE SGTT02( TRANS, N, NRHS, DL, D, DU, X, LDX, B, LDB, &
                      RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          TRANS
   INTEGER            LDB, LDX, N, NRHS
   REAL               RESID
!     ..
!     .. Array Arguments ..
   REAL               B( LDB, * ), D( * ), DL( * ), DU( * ), &
                      X( LDX, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ONE, ZERO
   PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            J
   REAL               ANORM, BNORM, EPS, XNORM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               SASUM, SLAMCH, SLANGT
   EXTERNAL           LSAME, SASUM, SLAMCH, SLANGT
!     ..
!     .. External Subroutines ..
   EXTERNAL           SLAGTM
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0 or NRHS = 0
!
   RESID = ZERO
   IF( N <= 0 .OR. NRHS == 0 ) &
      RETURN
!
!     Compute the maximum over the number of right hand sides of
!        norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ).
!
   IF( LSAME( TRANS, 'N' ) ) THEN
      ANORM = SLANGT( '1', N, DL, D, DU )
   ELSE
      ANORM = SLANGT( 'I', N, DL, D, DU )
   END IF
!
!     Exit with RESID = 1/EPS if ANORM = 0.
!
   EPS = SLAMCH( 'Epsilon' )
   IF( ANORM <= ZERO ) THEN
      RESID = ONE / EPS
      RETURN
   END IF
!
!     Compute B - op(A)*X and store in B.
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SLAGTM( TRANS, N, NRHS, -ONE, DL, D, DU, X, LDX, ONE, B, &
                LDB )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SLAGTM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   DO J = 1, NRHS
      BNORM = SASUM( N, B( 1, J ), 1 )
      XNORM = SASUM( N, X( 1, J ), 1 )
      IF( XNORM <= ZERO ) THEN
         RESID = ONE / EPS
      ELSE
         RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
      END IF
   ENDDO
!
   RETURN
!
!     End of SGTT02
!
END
                                                                                                                                                                                                                                                                                                            




