!> \brief \b DGET01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGET01( M, N, A, LDA, AFAC, LDAFAC, IPIV, RWORK,
!                          RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDAFAC, M, N
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * ), AFAC( LDAFAC, * ), RWORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGET01 reconstructs a matrix A from its L*U factorization and
!> computes the residual
!>    norm(L*U - A) / ( N * norm(A) * EPS ),
!> where EPS is the machine epsilon.
!> \endverbatim
!
!  Arguments:
!  ==========
!
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
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The original M x N matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] AFAC
!> \verbatim
!>          AFAC is DOUBLE PRECISION array, dimension (LDAFAC,N)
!>          The factored form of the matrix A.  AFAC contains the factors
!>          L and U from the L*U factorization as computed by DGETRF.
!>          Overwritten with the reconstructed matrix, and then with the
!>          difference L*U - A.
!> \endverbatim
!>
!> \param[in] LDAFAC
!> \verbatim
!>          LDAFAC is INTEGER
!>          The leading dimension of the array AFAC.  LDAFAC >= max(1,M).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices from DGETRF.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          norm(L*U - A) / ( N * norm(A) * EPS )
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
!> \ingroup double_lin
!
!  =====================================================================
   SUBROUTINE DGET01( M, N, A, LDA, AFAC, LDAFAC, IPIV, RWORK, &
                      RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, LDAFAC, M, N
   DOUBLE PRECISION   RESID
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * )
   DOUBLE PRECISION   A( LDA, * ), AFAC( LDAFAC, * ), RWORK( * )
!     ..
!
!  =====================================================================
!
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, J, K
   DOUBLE PRECISION   ANORM, EPS, T
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DDOT, DLAMCH, DLANGE
   EXTERNAL           DDOT, DLAMCH, DLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           DGEMV, DLASWP, DSCAL, DTRMV
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          DBLE, MIN
!     ..
!     .. Executable Statements ..
!
!     Quick exit if M = 0 or N = 0.
!
   IF( M <= 0 .OR. N <= 0 ) THEN
      RESID = ZERO
      RETURN
   END IF
!
!     Determine EPS and the norm of A.
!
   EPS = DLAMCH( 'Epsilon' )
   ANORM = DLANGE( '1', M, N, A, LDA, RWORK )
!
!     Compute the product L*U and overwrite AFAC with the result.
!     A column at a time of the product is obtained, starting with
!     column N.
!
   DO K = N, 1, -1
      IF( K > M ) THEN
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DTRMV( 'Lower', 'No transpose', 'Unit', M, AFAC, &
                     LDAFAC, AFAC( 1, K ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DTRMV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ELSE
!
!           Compute elements (K+1:M,K)
!
         T = AFAC( K, K )
         IF( K+1 <= M ) THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DSCAL( M-K, T, AFAC( K+1, K ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DGEMV( 'No transpose', M-K, K-1, ONE, &
                        AFAC( K+1, 1 ), LDAFAC, AFAC( 1, K ), 1, ONE, &
                        AFAC( K+1, K ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DGEMV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         END IF
!
!           Compute the (K,K) element
!
         AFAC( K, K ) = T + DDOT( K-1, AFAC( K, 1 ), LDAFAC, &
                        AFAC( 1, K ), 1 )
!
!           Compute elements (1:K-1,K)
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DTRMV( 'Lower', 'No transpose', 'Unit', K-1, AFAC, &
                     LDAFAC, AFAC( 1, K ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DTRMV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      END IF
   ENDDO
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLASWP( N, AFAC, LDAFAC, 1, MIN( M, N ), IPIV, -1 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLASWP : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute the difference  L*U - A  and store in AFAC.
!
   DO J = 1, N
      DO I = 1, M
         AFAC( I, J ) = AFAC( I, J ) - A( I, J )
      ENDDO
   ENDDO
!
!     Compute norm( L*U - A ) / ( N * norm(A) * EPS )
!
   RESID = DLANGE( '1', M, N, AFAC, LDAFAC, RWORK )
!
   IF( ANORM <= ZERO ) THEN
      IF( RESID /= ZERO ) &
         RESID = ONE / EPS
   ELSE
      RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS
   END IF
!
   RETURN
!
!     End of DGET01
!
END
                                                                                                                                                                                                                                                                                                            




