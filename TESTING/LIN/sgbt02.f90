!> \brief \b SGBT02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGBT02( TRANS, M, N, KL, KU, NRHS, A, LDA, X, LDX, B,
!                          LDB, RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            KL, KU, LDA, LDB, LDX, M, N, NRHS
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), B( LDB, * ), X( LDX, * ),
!                          RWORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGBT02 computes the residual for a solution of a banded system of
!> equations op(A)*X = B:
!>    RESID = norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ),
!> where op(A) = A or A**T, depending on TRANS, and EPS is the
!> machine epsilon.
!> The norm used is the 1-norm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the form of the system of equations:
!>          = 'N':  A    * X = B  (No transpose)
!>          = 'T':  A**T * X = B  (Transpose)
!>          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns of B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The original matrix A in band storage, stored in rows 1 to
!>          KL+KU+1.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,KL+KU+1).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is REAL array, dimension (LDX,NRHS)
!>          The computed solution vectors for the system of linear
!>          equations.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  If TRANS = 'N',
!>          LDX >= max(1,N); if TRANS = 'T' or 'C', LDX >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,NRHS)
!>          On entry, the right hand side vectors for the system of
!>          linear equations.
!>          On exit, B is overwritten with the difference B - A*X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  IF TRANS = 'N',
!>          LDB >= max(1,M); if TRANS = 'T' or 'C', LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (MAX(1,LRWORK)),
!>          where LRWORK >= M when TRANS = 'T' or 'C'; otherwise, RWORK
!>          is not referenced.
!> \endverbatim
!
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
!>          The maximum over the number of right hand sides of
!>          norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ).
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
   SUBROUTINE SGBT02( TRANS, M, N, KL, KU, NRHS, A, LDA, X, LDX, B, &
                      LDB, RWORK, RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          TRANS
   INTEGER            KL, KU, LDA, LDB, LDX, M, N, NRHS
   REAL               RESID
!     ..
!     .. Array Arguments ..
   REAL               A( LDA, * ), B( LDB, * ), X( LDX, * ), &
                      RWORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO, ONE
   PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I1, I2, J, KD, N1
   REAL               ANORM, BNORM, EPS, TEMP, XNORM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME, SISNAN
   REAL               SASUM, SLAMCH
   EXTERNAL           LSAME, SASUM, SISNAN, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           SGBMV
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if N = 0 pr NRHS = 0
!
   IF( M <= 0 .OR. N <= 0 .OR. NRHS <= 0 ) THEN
      RESID = ZERO
      RETURN
   END IF
!
!     Exit with RESID = 1/EPS if ANORM = 0.
!
   EPS = SLAMCH( 'Epsilon' )
   ANORM = ZERO
   IF( LSAME( TRANS, 'N' ) ) THEN
!
!        Find norm1(A).
!
      KD = KU + 1
      DO J = 1, N
         I1 = MAX( KD+1-J, 1 )
         I2 = MIN( KD+M-J, KL+KD )
         IF( I2 >= I1 ) THEN
            TEMP = SASUM( I2-I1+1, A( I1, J ), 1 )
            IF( ANORM < TEMP .OR. SISNAN( TEMP ) ) ANORM = TEMP
         END IF
      ENDDO
   ELSE
!
!        Find normI(A).
!
      DO I1 = 1, M
         RWORK( I1 ) = ZERO
      ENDDO
      DO J = 1, N
         KD = KU + 1 - J
         DO I1 = MAX( 1, J-KU ), MIN( M, J+KL )
            RWORK( I1 ) = RWORK( I1 ) + ABS( A( KD+I1, J ) )
         ENDDO
      ENDDO
      DO I1 = 1, M
         TEMP = RWORK( I1 )
         IF( ANORM < TEMP .OR. SISNAN( TEMP ) ) ANORM = TEMP
      ENDDO
   END IF
   IF( ANORM <= ZERO ) THEN
      RESID = ONE / EPS
      RETURN
   END IF
!
   IF( LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' ) ) THEN
      N1 = N
   ELSE
      N1 = M
   END IF
!
!     Compute B - op(A)*X
!
   DO J = 1, NRHS
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL SGBMV( TRANS, M, N, KL, KU, -ONE, A, LDA, X( 1, J ), 1, &
                  ONE, B( 1, J ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : SGBMV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   ENDDO
!
!     Compute the maximum over the number of right hand sides of
!        norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ).
!
   RESID = ZERO
   DO J = 1, NRHS
      BNORM = SASUM( N1, B( 1, J ), 1 )
      XNORM = SASUM( N1, X( 1, J ), 1 )
      IF( XNORM <= ZERO ) THEN
         RESID = ONE / EPS
      ELSE
         RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
      END IF
   ENDDO
!
   RETURN
!
!     End of SGBT02
!
END
                                                                                                                                                                                                                                                                                                            




