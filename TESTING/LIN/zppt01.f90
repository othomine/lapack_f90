!> \brief \b ZPPT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZPPT01( UPLO, N, A, AFAC, RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            N
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( * ), AFAC( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZPPT01 reconstructs a Hermitian positive definite packed matrix A
!> from its L*L' or U'*U factorization and computes the residual
!>    norm( L*L' - A ) / ( N * norm(A) * EPS ) or
!>    norm( U'*U - A ) / ( N * norm(A) * EPS ),
!> where EPS is the machine epsilon, L' is the conjugate transpose of
!> L, and U' is the conjugate transpose of U.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          Hermitian matrix A is stored:
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows and columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          The original Hermitian matrix A, stored as a packed
!>          triangular matrix.
!> \endverbatim
!>
!> \param[in,out] AFAC
!> \verbatim
!>          AFAC is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          On entry, the factor L or U from the L*L' or U'*U
!>          factorization of A, stored as a packed triangular matrix.
!>          Overwritten with the reconstructed matrix, and then with the
!>          difference L*L' - A (or U'*U - A).
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
!>          If UPLO = 'L', norm(L*L' - A) / ( N * norm(A) * EPS )
!>          If UPLO = 'U', norm(U'*U - A) / ( N * norm(A) * EPS )
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
   SUBROUTINE ZPPT01( UPLO, N, A, AFAC, RWORK, RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            N
   DOUBLE PRECISION   RESID
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   RWORK( * )
   COMPLEX*16         A( * ), AFAC( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, K, KC
   DOUBLE PRECISION   ANORM, EPS, TR
   COMPLEX*16         TC
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   DOUBLE PRECISION   DLAMCH, ZLANHP
   COMPLEX*16         ZDOTC
   EXTERNAL           LSAME, DLAMCH, ZLANHP, ZDOTC
!     ..
!     .. External Subroutines ..
   EXTERNAL           ZHPR, ZSCAL, ZTPMV
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          DBLE, DIMAG
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0
!
   IF( N <= 0 ) THEN
      RESID = ZERO
      RETURN
   END IF
!
!     Exit with RESID = 1/EPS if ANORM = 0.
!
   EPS = DLAMCH( 'Epsilon' )
   ANORM = ZLANHP( '1', UPLO, N, A, RWORK )
   IF( ANORM <= ZERO ) THEN
      RESID = ONE / EPS
      RETURN
   END IF
!
!     Check the imaginary parts of the diagonal elements and return with
!     an error code if any are nonzero.
!
   KC = 1
   IF( LSAME( UPLO, 'U' ) ) THEN
      DO K = 1, N
         IF( DIMAG( AFAC( KC ) ) /= ZERO ) THEN
            RESID = ONE / EPS
            RETURN
         END IF
         KC = KC + K + 1
      ENDDO
   ELSE
      DO K = 1, N
         IF( DIMAG( AFAC( KC ) ) /= ZERO ) THEN
            RESID = ONE / EPS
            RETURN
         END IF
         KC = KC + N - K + 1
      ENDDO
   END IF
!
!     Compute the product U'*U, overwriting U.
!
   IF( LSAME( UPLO, 'U' ) ) THEN
      KC = ( N*( N-1 ) ) / 2 + 1
      DO K = N, 1, -1
!
!           Compute the (K,K) element of the result.
!
         TR = DBLE( ZDOTC( K, AFAC( KC ), 1, AFAC( KC ), 1 ) )
         AFAC( KC+K-1 ) = TR
!
!           Compute the rest of column K.
!
         IF( K > 1 ) THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZTPMV( 'Upper', 'Conjugate', 'Non-unit', K-1, AFAC, &
                        AFAC( KC ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZTPMV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            KC = KC - ( K-1 )
         END IF
      ENDDO
!
!        Compute the difference  L*L' - A
!
      KC = 1
      DO K = 1, N
         DO I = 1, K - 1
            AFAC( KC+I-1 ) = AFAC( KC+I-1 ) - A( KC+I-1 )
         ENDDO
         AFAC( KC+K-1 ) = AFAC( KC+K-1 ) - DBLE( A( KC+K-1 ) )
         KC = KC + K
      ENDDO
!
!     Compute the product L*L', overwriting L.
!
   ELSE
      KC = ( N*( N+1 ) ) / 2
      DO K = N, 1, -1
!
!           Add a multiple of column K of the factor L to each of
!           columns K+1 through N.
!
         IF( K < N )  THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZHPR( 'Lower', N-K, ONE, AFAC( KC+1 ), 1, &
                       AFAC( KC+N-K+1 ) )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZHPR : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ENDIF
!
!           Scale column K by the diagonal element.
!
         TC = AFAC( KC )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZSCAL( N-K+1, TC, AFAC( KC ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZSCAL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
         KC = KC - ( N-K+2 )
      ENDDO
!
!        Compute the difference  U'*U - A
!
      KC = 1
      DO K = 1, N
         AFAC( KC ) = AFAC( KC ) - DBLE( A( KC ) )
         DO I = K + 1, N
            AFAC( KC+I-K ) = AFAC( KC+I-K ) - A( KC+I-K )
         ENDDO
         KC = KC + N - K + 1
      ENDDO
   END IF
!
!     Compute norm( L*U - A ) / ( N * norm(A) * EPS )
!
   RESID = ZLANHP( '1', UPLO, N, AFAC, RWORK )
!
   RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS
!
   RETURN
!
!     End of ZPPT01
!
END
                                                                                                                                                                                                                                                                                                            




