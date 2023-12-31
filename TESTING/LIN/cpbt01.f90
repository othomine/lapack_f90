!> \brief \b CPBT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CPBT01( UPLO, N, KD, A, LDA, AFAC, LDAFAC, RWORK,
!                          RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            KD, LDA, LDAFAC, N
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            A( LDA, * ), AFAC( LDAFAC, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CPBT01 reconstructs a Hermitian positive definite band matrix A from
!> its L*L' or U'*U factorization and computes the residual
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
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of super-diagonals of the matrix A if UPLO = 'U',
!>          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The original Hermitian band matrix A.  If UPLO = 'U', the
!>          upper triangular part of A is stored as a band matrix; if
!>          UPLO = 'L', the lower triangular part of A is stored.  The
!>          columns of the appropriate triangle are stored in the columns
!>          of A and the diagonals of the triangle are stored in the rows
!>          of A.  See CPBTRF for further details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER.
!>          The leading dimension of the array A.  LDA >= max(1,KD+1).
!> \endverbatim
!>
!> \param[in] AFAC
!> \verbatim
!>          AFAC is COMPLEX array, dimension (LDAFAC,N)
!>          The factored form of the matrix A.  AFAC contains the factor
!>          L or U from the L*L' or U'*U factorization in band storage
!>          format, as computed by CPBTRF.
!> \endverbatim
!>
!> \param[in] LDAFAC
!> \verbatim
!>          LDAFAC is INTEGER
!>          The leading dimension of the array AFAC.
!>          LDAFAC >= max(1,KD+1).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
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
!> \ingroup complex_lin
!
!  =====================================================================
   SUBROUTINE CPBT01( UPLO, N, KD, A, LDA, AFAC, LDAFAC, RWORK, &
                      RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            KD, LDA, LDAFAC, N
   REAL               RESID
!     ..
!     .. Array Arguments ..
   REAL               RWORK( * )
   COMPLEX            A( LDA, * ), AFAC( LDAFAC, * )
!     ..
!
!  =====================================================================
!
!
!     .. Parameters ..
   REAL               ZERO, ONE
   PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, J, K, KC, KLEN, ML, MU
   REAL               AKK, ANORM, EPS
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               CLANHB, SLAMCH
   COMPLEX            CDOTC
   EXTERNAL           LSAME, CLANHB, SLAMCH, CDOTC
!     ..
!     .. External Subroutines ..
   EXTERNAL           CHER, CSSCAL, CTRMV
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          AIMAG, MAX, MIN, REAL
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0.
!
   IF( N <= 0 ) THEN
      RESID = ZERO
      RETURN
   END IF
!
!     Exit with RESID = 1/EPS if ANORM = 0.
!
   EPS = SLAMCH( 'Epsilon' )
   ANORM = CLANHB( '1', UPLO, N, KD, A, LDA, RWORK )
   IF( ANORM <= ZERO ) THEN
      RESID = ONE / EPS
      RETURN
   END IF
!
!     Check the imaginary parts of the diagonal elements and return with
!     an error code if any are nonzero.
!
   IF( LSAME( UPLO, 'U' ) ) THEN
      DO J = 1, N
         IF( AIMAG( AFAC( KD+1, J ) ) /= ZERO ) THEN
            RESID = ONE / EPS
            RETURN
         END IF
      ENDDO
   ELSE
      DO J = 1, N
         IF( AIMAG( AFAC( 1, J ) ) /= ZERO ) THEN
            RESID = ONE / EPS
            RETURN
         END IF
      ENDDO
   END IF
!
!     Compute the product U'*U, overwriting U.
!
   IF( LSAME( UPLO, 'U' ) ) THEN
      DO K = N, 1, -1
         KC = MAX( 1, KD+2-K )
         KLEN = KD + 1 - KC
!
!           Compute the (K,K) element of the result.
!
         AKK = REAL( &
            CDOTC( KLEN+1, AFAC( KC, K ), 1, AFAC( KC, K ), 1 ) )
         AFAC( KD+1, K ) = AKK
!
!           Compute the rest of column K.
!
         IF( KLEN > 0 )  THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL CTRMV( 'Upper', 'Conjugate', 'Non-unit', KLEN, &
                        AFAC( KD+1, K-KLEN ), LDAFAC-1, &
                        AFAC( KC, K ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : CTRMV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ENDIF
!
      ENDDO
!
!     UPLO = 'L':  Compute the product L*L', overwriting L.
!
   ELSE
      DO K = N, 1, -1
         KLEN = MIN( KD, N-K )
!
!           Add a multiple of column K of the factor L to each of
!           columns K+1 through N.
!
         IF( KLEN > 0 )  THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL CHER( 'Lower', KLEN, ONE, AFAC( 2, K ), 1, &
                       AFAC( 1, K+1 ), LDAFAC-1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : CHER : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ENDIF
!
!           Scale column K by the diagonal element.
!
         AKK = REAL( AFAC( 1, K ) )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL CSSCAL( KLEN+1, AKK, AFAC( 1, K ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CSSCAL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
      ENDDO
   END IF
!
!     Compute the difference  L*L' - A  or  U'*U - A.
!
   IF( LSAME( UPLO, 'U' ) ) THEN
      DO J = 1, N
         MU = MAX( 1, KD+2-J )
         DO I = MU, KD + 1
            AFAC( I, J ) = AFAC( I, J ) - A( I, J )
         ENDDO
      ENDDO
   ELSE
      DO J = 1, N
         ML = MIN( KD+1, N-J+1 )
         DO I = 1, ML
            AFAC( I, J ) = AFAC( I, J ) - A( I, J )
         ENDDO
      ENDDO
   END IF
!
!     Compute norm( L*L' - A ) / ( N * norm(A) * EPS )
!
   RESID = CLANHB( '1', UPLO, N, KD, AFAC, LDAFAC, RWORK )
!
   RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
!
   RETURN
!
!     End of CPBT01
!
END
                                                                                                                                                                                                                                                                                                            




