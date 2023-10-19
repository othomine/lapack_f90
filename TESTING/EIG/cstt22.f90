!> \brief \b CSTT22
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSTT22( N, M, KBAND, AD, AE, SD, SE, U, LDU, WORK,
!                          LDWORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            KBAND, LDU, LDWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               AD( * ), AE( * ), RESULT( 2 ), RWORK( * ),
!      $                   SD( * ), SE( * )
!       COMPLEX            U( LDU, * ), WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSTT22  checks a set of M eigenvalues and eigenvectors,
!>
!>     A U = U S
!>
!> where A is Hermitian tridiagonal, the columns of U are unitary,
!> and S is diagonal (if KBAND=0) or Hermitian tridiagonal (if KBAND=1).
!> Two tests are performed:
!>
!>    RESULT(1) = | U* A U - S | / ( |A| m ulp )
!>
!>    RESULT(2) = | I - U*U | / ( m ulp )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The size of the matrix.  If it is zero, CSTT22 does nothing.
!>          It must be at least zero.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of eigenpairs to check.  If it is zero, CSTT22
!>          does nothing.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] KBAND
!> \verbatim
!>          KBAND is INTEGER
!>          The bandwidth of the matrix S.  It may only be zero or one.
!>          If zero, then S is diagonal, and SE is not referenced.  If
!>          one, then S is Hermitian tri-diagonal.
!> \endverbatim
!>
!> \param[in] AD
!> \verbatim
!>          AD is REAL array, dimension (N)
!>          The diagonal of the original (unfactored) matrix A.  A is
!>          assumed to be Hermitian tridiagonal.
!> \endverbatim
!>
!> \param[in] AE
!> \verbatim
!>          AE is REAL array, dimension (N)
!>          The off-diagonal of the original (unfactored) matrix A.  A
!>          is assumed to be Hermitian tridiagonal.  AE(1) is ignored,
!>          AE(2) is the (1,2) and (2,1) element, etc.
!> \endverbatim
!>
!> \param[in] SD
!> \verbatim
!>          SD is REAL array, dimension (N)
!>          The diagonal of the (Hermitian tri-) diagonal matrix S.
!> \endverbatim
!>
!> \param[in] SE
!> \verbatim
!>          SE is REAL array, dimension (N)
!>          The off-diagonal of the (Hermitian tri-) diagonal matrix S.
!>          Not referenced if KBSND=0.  If KBAND=1, then AE(1) is
!>          ignored, SE(2) is the (1,2) and (2,1) element, etc.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is REAL array, dimension (LDU, N)
!>          The unitary matrix in the decomposition.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U.  LDU must be at least N.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (LDWORK, M+1)
!> \endverbatim
!>
!> \param[in] LDWORK
!> \verbatim
!>          LDWORK is INTEGER
!>          The leading dimension of WORK.  LDWORK must be at least
!>          max(1,M).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (2)
!>          The values computed by the two tests described above.  The
!>          values are currently limited to 1/ulp, to avoid overflow.
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
!> \ingroup complex_eig
!
!  =====================================================================
   SUBROUTINE CSTT22( N, M, KBAND, AD, AE, SD, SE, U, LDU, WORK, &
                      LDWORK, RWORK, RESULT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            KBAND, LDU, LDWORK, M, N
!     ..
!     .. Array Arguments ..
   REAL               AD( * ), AE( * ), RESULT( 2 ), RWORK( * ), &
                      SD( * ), SE( * )
   COMPLEX            U( LDU, * ), WORK( LDWORK, * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, J, K
   REAL               ANORM, ULP, UNFL, WNORM
   COMPLEX            AUKJ
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   REAL               CLANGE, CLANSY, SLAMCH
   EXTERNAL           CLANGE, CLANSY, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMM
!     ..
!     .. Executable Statements ..
!
   RESULT( 1:2 ) = 0.0E+0
   IF( N <= 0 .OR. M <= 0 ) RETURN
!
   UNFL = SLAMCH( 'Safe minimum' )
   ULP = SLAMCH( 'Epsilon' )
!
!     Do Test 1
!
!     Compute the 1-norm of A.
!
   IF( N > 1 ) THEN
      ANORM = ABS( AD( 1 ) ) + ABS( AE( 1 ) )
      DO J = 2, N - 1
         ANORM = MAX( ANORM, ABS( AD( J ) )+ABS( AE( J ) )+ &
                 ABS( AE( J-1 ) ) )
      ENDDO
      ANORM = MAX( ANORM, ABS( AD( N ) )+ABS( AE( N-1 ) ) )
   ELSE
      ANORM = ABS( AD( 1 ) )
   END IF
   ANORM = MAX( ANORM, UNFL )
!
!     Norm of U*AU - S
!
   DO I = 1, M
      DO J = 1, M
         WORK( I, J ) = (0.0E+0,0.0E+0)
         DO K = 1, N
            AUKJ = AD( K )*U( K, J )
            IF( K /= N ) AUKJ = AUKJ + AE( K )*U( K+1, J )
            IF( K /= 1 ) AUKJ = AUKJ + AE( K-1 )*U( K-1, J )
            WORK( I, J ) = WORK( I, J ) + U( K, I )*AUKJ
         ENDDO
      ENDDO
      WORK( I, I ) = WORK( I, I ) - SD( I )
      IF( KBAND == 1 ) THEN
         IF( I /= 1 ) WORK( I, I-1 ) = WORK( I, I-1 ) - SE( I-1 )
         IF( I /= N ) WORK( I, I+1 ) = WORK( I, I+1 ) - SE( I )
      END IF
   ENDDO
!
   WNORM = CLANSY( '1', 'L', M, WORK, M, RWORK )
!
   IF( ANORM > WNORM ) THEN
      RESULT( 1 ) = ( WNORM / ANORM ) / ( M*ULP )
   ELSE
      IF( ANORM < 1.0E+0 ) THEN
         RESULT( 1 ) = ( MIN( WNORM, M*ANORM ) / ANORM ) / ( M*ULP )
      ELSE
         RESULT( 1 ) = MIN( WNORM / ANORM, REAL( M ) ) / ( M*ULP )
      END IF
   END IF
!
!     Do Test 2
!
!     Compute  U*U - I
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CGEMM( 'T', 'N', M, M, N, (1.0E+0,0.0E+0), U, LDU, U, LDU, (0.0E+0,0.0E+0), WORK, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   FORALL (J = 1:M) WORK( J, J ) = WORK( J, J ) - 1.0E+0
!
   RESULT( 2 ) = MIN( REAL( M ), CLANGE( '1', M, M, WORK, M, &
                 RWORK ) ) / ( M*ULP )
!
   RETURN
!
!     End of CSTT22
!
END




