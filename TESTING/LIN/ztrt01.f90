!> \brief \b ZTRT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTRT01( UPLO, DIAG, N, A, LDA, AINV, LDAINV, RCOND,
!                          RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, UPLO
!       INTEGER            LDA, LDAINV, N
!       DOUBLE PRECISION   RCOND, RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( LDA, * ), AINV( LDAINV, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTRT01 computes the residual for a triangular matrix A times its
!> inverse:
!>    RESID = norm( A*AINV - I ) / ( N * norm(A) * norm(AINV) * EPS ),
!> where EPS is the machine epsilon.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A is upper or lower triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A is unit triangular.
!>          = 'N':  Non-unit triangular
!>          = 'U':  Unit triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The triangular matrix A.  If UPLO = 'U', the leading n by n
!>          upper triangular part of the array A contains the upper
!>          triangular matrix, and the strictly lower triangular part of
!>          A is not referenced.  If UPLO = 'L', the leading n by n lower
!>          triangular part of the array A contains the lower triangular
!>          matrix, and the strictly upper triangular part of A is not
!>          referenced.  If DIAG = 'U', the diagonal elements of A are
!>          also not referenced and are assumed to be 1.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] AINV
!> \verbatim
!>          AINV is COMPLEX*16 array, dimension (LDAINV,N)
!>          On entry, the (triangular) inverse of the matrix A, in the
!>          same storage format as A.
!>          On exit, the contents of AINV are destroyed.
!> \endverbatim
!>
!> \param[in] LDAINV
!> \verbatim
!>          LDAINV is INTEGER
!>          The leading dimension of the array AINV.  LDAINV >= max(1,N).
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is DOUBLE PRECISION
!>          The reciprocal condition number of A, computed as
!>          1/(norm(A) * norm(AINV)).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          norm(A*AINV - I) / ( N * norm(A) * norm(AINV) * EPS )
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
!> \ingroup complex16_lin
!
!  =====================================================================
   SUBROUTINE ZTRT01( UPLO, DIAG, N, A, LDA, AINV, LDAINV, RCOND, &
                      RWORK, RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          DIAG, UPLO
   INTEGER            LDA, LDAINV, N
   DOUBLE PRECISION   RCOND, RESID
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   RWORK( * )
   COMPLEX*16         A( LDA, * ), AINV( LDAINV, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            J
   DOUBLE PRECISION   AINVNM, ANORM, EPS
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   DOUBLE PRECISION   DLAMCH, ZLANTR
   EXTERNAL           LSAME, DLAMCH, ZLANTR
!     ..
!     .. External Subroutines ..
   EXTERNAL           ZTRMV
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          DBLE
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0
!
   IF( N <= 0 ) THEN
      RCOND = ONE
      RESID = ZERO
      RETURN
   END IF
!
!     Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.
!
   EPS = DLAMCH( 'Epsilon' )
   ANORM = ZLANTR( '1', UPLO, DIAG, N, N, A, LDA, RWORK )
   AINVNM = ZLANTR( '1', UPLO, DIAG, N, N, AINV, LDAINV, RWORK )
   IF( ANORM <= ZERO .OR. AINVNM <= ZERO ) THEN
      RCOND = ZERO
      RESID = ONE / EPS
      RETURN
   END IF
   RCOND = ( ONE / ANORM ) / AINVNM
!
!     Set the diagonal of AINV to 1 if AINV has unit diagonal.
!
   IF( LSAME( DIAG, 'U' ) ) THEN
      DO J = 1, N
         AINV( J, J ) = ONE
      ENDDO
   END IF
!
!     Compute A * AINV, overwriting AINV.
!
   IF( LSAME( UPLO, 'U' ) ) THEN
      DO J = 1, N
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZTRMV( 'Upper', 'No transpose', DIAG, J, A, LDA, &
                     AINV( 1, J ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZTRMV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
   ELSE
      DO J = 1, N
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZTRMV( 'Lower', 'No transpose', DIAG, N-J+1, A( J, J ), &
                     LDA, AINV( J, J ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZTRMV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
   END IF
!
!     Subtract 1 from each diagonal element to form A*AINV - I.
!
   DO J = 1, N
      AINV( J, J ) = AINV( J, J ) - ONE
   ENDDO
!
!     Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)
!
   RESID = ZLANTR( '1', UPLO, 'Non-unit', N, N, AINV, LDAINV, RWORK )
!
   RESID = ( ( RESID*RCOND ) / DBLE( N ) ) / EPS
!
   RETURN
!
!     End of ZTRT01
!
END
                                                                                                                                                                                                                                                                                                            




