!> \brief \b SPPT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SPPT01( UPLO, N, A, AFAC, RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            N
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               A( * ), AFAC( * ), RWORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SPPT01 reconstructs a symmetric positive definite packed matrix A
!> from its L*L' or U'*U factorization and computes the residual
!>    norm( L*L' - A ) / ( N * norm(A) * EPS ) or
!>    norm( U'*U - A ) / ( N * norm(A) * EPS ),
!> where EPS is the machine epsilon.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix A is stored:
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
!>          A is REAL array, dimension (N*(N+1)/2)
!>          The original symmetric matrix A, stored as a packed
!>          triangular matrix.
!> \endverbatim
!>
!> \param[in,out] AFAC
!> \verbatim
!>          AFAC is REAL array, dimension (N*(N+1)/2)
!>          On entry, the factor L or U from the L*L' or U'*U
!>          factorization of A, stored as a packed triangular matrix.
!>          Overwritten with the reconstructed matrix, and then with the
!>          difference L*L' - A (or U'*U - A).
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
!> \ingroup single_lin
!
!  =====================================================================
   SUBROUTINE SPPT01( UPLO, N, A, AFAC, RWORK, RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            N
   REAL               RESID
!     ..
!     .. Array Arguments ..
   REAL               A( * ), AFAC( * ), RWORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO, ONE
   PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, K, KC, NPP
   REAL               ANORM, EPS, T
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               SDOT, SLAMCH, SLANSP
   EXTERNAL           LSAME, SDOT, SLAMCH, SLANSP
!     ..
!     .. External Subroutines ..
   EXTERNAL           SSCAL, SSPR, STPMV
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          REAL
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
   EPS = SLAMCH( 'Epsilon' )
   ANORM = SLANSP( '1', UPLO, N, A, RWORK )
   IF( ANORM <= ZERO ) THEN
      RESID = ONE / EPS
      RETURN
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
         T = SDOT( K, AFAC( KC ), 1, AFAC( KC ), 1 )
         AFAC( KC+K-1 ) = T
!
!           Compute the rest of column K.
!
         IF( K > 1 ) THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL STPMV( 'Upper', 'Transpose', 'Non-unit', K-1, AFAC, &
                        AFAC( KC ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : STPMV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            KC = KC - ( K-1 )
         END IF
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
            CALL SSPR( 'Lower', N-K, ONE, AFAC( KC+1 ), 1, &
                       AFAC( KC+N-K+1 ) )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SSPR : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ENDIF
!
!           Scale column K by the diagonal element.
!
         T = AFAC( KC )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SSCAL( N-K+1, T, AFAC( KC ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SSCAL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
         KC = KC - ( N-K+2 )
      ENDDO
   END IF
!
!     Compute the difference  L*L' - A (or U'*U - A).
!
   NPP = N*( N+1 ) / 2
   DO I = 1, NPP
      AFAC( I ) = AFAC( I ) - A( I )
   ENDDO
!
!     Compute norm( L*U - A ) / ( N * norm(A) * EPS )
!
   RESID = SLANSP( '1', UPLO, N, AFAC, RWORK )
!
   RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
!
   RETURN
!
!     End of SPPT01
!
END
                                                                                                                                                                                                                                                                                                            




