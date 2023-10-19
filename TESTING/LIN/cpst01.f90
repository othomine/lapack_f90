!> \brief \b CPST01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CPST01( UPLO, N, A, LDA, AFAC, LDAFAC, PERM, LDPERM,
!                          PIV, RWORK, RESID, RANK )
!
!       .. Scalar Arguments ..
!       REAL               RESID
!       INTEGER            LDA, LDAFAC, LDPERM, N, RANK
!       CHARACTER          UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), AFAC( LDAFAC, * ),
!      $                   PERM( LDPERM, * )
!       REAL               RWORK( * )
!       INTEGER            PIV( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CPST01 reconstructs an Hermitian positive semidefinite matrix A
!> from its L or U factors and the permutation matrix P and computes
!> the residual
!>    norm( P*L*L'*P' - A ) / ( N * norm(A) * EPS ) or
!>    norm( P*U'*U*P' - A ) / ( N * norm(A) * EPS ),
!> where EPS is the machine epsilon, L' is the conjugate transpose of L,
!> and U' is the conjugate transpose of U.
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          The original Hermitian matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N)
!> \endverbatim
!>
!> \param[in] AFAC
!> \verbatim
!>          AFAC is COMPLEX array, dimension (LDAFAC,N)
!>          The factor L or U from the L*L' or U'*U
!>          factorization of A.
!> \endverbatim
!>
!> \param[in] LDAFAC
!> \verbatim
!>          LDAFAC is INTEGER
!>          The leading dimension of the array AFAC.  LDAFAC >= max(1,N).
!> \endverbatim
!>
!> \param[out] PERM
!> \verbatim
!>          PERM is COMPLEX array, dimension (LDPERM,N)
!>          Overwritten with the reconstructed matrix, and then with the
!>          difference P*L*L'*P' - A (or P*U'*U*P' - A)
!> \endverbatim
!>
!> \param[in] LDPERM
!> \verbatim
!>          LDPERM is INTEGER
!>          The leading dimension of the array PERM.
!>          LDAPERM >= max(1,N).
!> \endverbatim
!>
!> \param[in] PIV
!> \verbatim
!>          PIV is INTEGER array, dimension (N)
!>          PIV is such that the nonzero entries are
!>          P( PIV( K ), K ) = 1.
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
!>
!> \param[in] RANK
!> \verbatim
!>          RANK is INTEGER
!>          number of nonzero singular values of A.
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
   SUBROUTINE CPST01( UPLO, N, A, LDA, AFAC, LDAFAC, PERM, LDPERM, &
                      PIV, RWORK, RESID, RANK )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   REAL               RESID
   INTEGER            LDA, LDAFAC, LDPERM, N, RANK
   CHARACTER          UPLO
!     ..
!     .. Array Arguments ..
   COMPLEX            A( LDA, * ), AFAC( LDAFAC, * ), &
                      PERM( LDPERM, * )
   REAL               RWORK( * )
   INTEGER            PIV( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO, ONE
   PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
   COMPLEX            CZERO
   PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
   COMPLEX            TC
   REAL               ANORM, EPS, TR
   INTEGER            I, J, K
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   COMPLEX            CDOTC
   REAL               CLANHE, SLAMCH
   LOGICAL            LSAME
   EXTERNAL           CDOTC, CLANHE, SLAMCH, LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           CHER, CSCAL, CTRMV
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          AIMAG, CONJG, REAL
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
   ANORM = CLANHE( '1', UPLO, N, A, LDA, RWORK )
   IF( ANORM <= ZERO ) THEN
      RESID = ONE / EPS
      RETURN
   END IF
!
!     Check the imaginary parts of the diagonal elements and return with
!     an error code if any are nonzero.
!
   DO J = 1, N
      IF( AIMAG( AFAC( J, J ) ) /= ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
      ENDDO
!
!     Compute the product U'*U, overwriting U.
!
   IF( LSAME( UPLO, 'U' ) ) THEN
!
      IF( RANK < N ) THEN
         DO J = RANK + 1, N
            DO I = RANK + 1, J
               AFAC( I, J ) = CZERO
               ENDDO
            ENDDO
      END IF
!
      DO K = N, 1, -1
!
!           Compute the (K,K) element of the result.
!
         TR = REAL( CDOTC( K, AFAC( 1, K ), 1, AFAC( 1, K ), 1 ) )
         AFAC( K, K ) = TR
!
!           Compute the rest of column K.
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL CTRMV( 'Upper', 'Conjugate', 'Non-unit', K-1, AFAC, &
                     LDAFAC, AFAC( 1, K ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CTRMV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
         ENDDO
!
!     Compute the product L*L', overwriting L.
!
   ELSE
!
      IF( RANK < N ) THEN
         DO J = RANK + 1, N
            DO I = J, N
               AFAC( I, J ) = CZERO
               ENDDO
            ENDDO
      END IF
!
      DO K = N, 1, -1
!           Add a multiple of column K of the factor L to each of
!           columns K+1 through N.
!
         IF( K+1 <= N )  THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL CHER( 'Lower', N-K, ONE, AFAC( K+1, K ), 1, &
                       AFAC( K+1, K+1 ), LDAFAC )
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
         TC = AFAC( K, K )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL CSCAL( N-K+1, TC, AFAC( K, K ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CSCAL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         ENDDO
!
   END IF
!
!        Form P*L*L'*P' or P*U'*U*P'
!
   IF( LSAME( UPLO, 'U' ) ) THEN
!
      DO J = 1, N
         DO I = 1, N
            IF( PIV( I ) <= PIV( J ) ) THEN
               IF( I <= J ) THEN
                  PERM( PIV( I ), PIV( J ) ) = AFAC( I, J )
               ELSE
                  PERM( PIV( I ), PIV( J ) ) = CONJG( AFAC( J, I ) )
               END IF
            END IF
            ENDDO
         ENDDO
!
!
   ELSE
!
      DO J = 1, N
         DO I = 1, N
            IF( PIV( I ) >= PIV( J ) ) THEN
               IF( I >= J ) THEN
                  PERM( PIV( I ), PIV( J ) ) = AFAC( I, J )
               ELSE
                  PERM( PIV( I ), PIV( J ) ) = CONJG( AFAC( J, I ) )
               END IF
            END IF
            ENDDO
         ENDDO
!
   END IF
!
!     Compute the difference  P*L*L'*P' - A (or P*U'*U*P' - A).
!
   IF( LSAME( UPLO, 'U' ) ) THEN
      DO J = 1, N
         DO I = 1, J - 1
            PERM( I, J ) = PERM( I, J ) - A( I, J )
            ENDDO
         PERM( J, J ) = PERM( J, J ) - REAL( A( J, J ) )
         ENDDO
   ELSE
      DO J = 1, N
         PERM( J, J ) = PERM( J, J ) - REAL( A( J, J ) )
         DO I = J + 1, N
            PERM( I, J ) = PERM( I, J ) - A( I, J )
            ENDDO
         ENDDO
   END IF
!
!     Compute norm( P*L*L'P - A ) / ( N * norm(A) * EPS ), or
!     ( P*U'*U*P' - A )/ ( N * norm(A) * EPS ).
!
   RESID = CLANHE( '1', UPLO, N, PERM, LDAFAC, RWORK )
!
   RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
!
   RETURN
!
!     End of CPST01
!
END
                                                                                                                                                                                                                                                                                                            




