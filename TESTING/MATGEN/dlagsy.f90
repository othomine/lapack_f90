!> \brief \b DLAGSY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAGSY( N, K, D, A, LDA, ISEED, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       DOUBLE PRECISION   A( LDA, * ), D( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAGSY generates a real symmetric matrix A, by pre- and post-
!> multiplying a real diagonal matrix D with a random orthogonal matrix:
!> A = U*D*U'. The semi-bandwidth may then be reduced to k by additional
!> orthogonal transformations.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of nonzero subdiagonals within the band of A.
!>          0 <= K <= N-1.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The diagonal elements of the diagonal matrix D.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The generated n by n symmetric matrix A (the full matrix is
!>          stored).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= N.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry, the seed of the random number generator; the array
!>          elements must be between 0 and 4095, and ISEED(4) must be
!>          odd.
!>          On exit, the seed is updated.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup double_matgen
!
!  =====================================================================
   SUBROUTINE DLAGSY( N, K, D, A, LDA, ISEED, WORK, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, K, LDA, N
!     ..
!     .. Array Arguments ..
   INTEGER            ISEED( 4 )
   DOUBLE PRECISION   A( LDA, * ), D( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE, HALF
   PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, J
   DOUBLE PRECISION   ALPHA, TAU, WA, WB, WN
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Subroutines ..
   EXTERNAL           DAXPY, DGEMV, DGER, DLARNV, DSCAL, DSYMV, &
                      DSYR2, XERBLA
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DDOT, DNRM2
   EXTERNAL           DDOT, DNRM2
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX, SIGN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
   INFO = 0
   IF( N < 0 ) THEN
      INFO = -1
   ELSE IF( K < 0 .OR. K > N-1 ) THEN
      INFO = -2
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -5
   END IF
   IF( INFO < 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'DLAGSY', -INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : XERBLA : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      RETURN
   END IF
!
!     initialize lower triangle of A to diagonal matrix
!
   DO J = 1, N
      DO I = J + 1, N
         A( I, J ) = ZERO
      ENDDO
   ENDDO
   DO I = 1, N
      A( I, I ) = D( I )
   ENDDO
!
!     Generate lower triangle of symmetric matrix
!
   DO I = N - 1, 1, -1
!
!        generate random reflection
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLARNV( 3, ISEED, N-I+1, WORK )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      WN = DNRM2( N-I+1, WORK, 1 )
      WA = SIGN( WN, WORK( 1 ) )
      IF( WN == ZERO ) THEN
         TAU = ZERO
      ELSE
         WB = WORK( 1 ) + WA
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DSCAL( N-I, ONE / WB, WORK( 2 ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         WORK( 1 ) = ONE
         TAU = WB / WA
      END IF
!
!        apply random reflection to A(i:n,i:n) from the left
!        and the right
!
!        compute  y := tau * A * u
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DSYMV( 'Lower', N-I+1, TAU, A( I, I ), LDA, WORK, 1, ZERO, &
                  WORK( N+1 ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DSYMV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        compute  v := y - 1/2 * tau * ( y, u ) * u
!
      ALPHA = -HALF*TAU*DDOT( N-I+1, WORK( N+1 ), 1, WORK, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DAXPY( N-I+1, ALPHA, WORK, 1, WORK( N+1 ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DAXPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        apply the transformation as a rank-2 update to A(i:n,i:n)
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DSYR2( 'Lower', N-I+1, -ONE, WORK, 1, WORK( N+1 ), 1, &
                  A( I, I ), LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DSYR2 : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   ENDDO
!
!     Reduce number of subdiagonals to K
!
   DO I = 1, N - 1 - K
!
!        generate reflection to annihilate A(k+i+1:n,i)
!
      WN = DNRM2( N-K-I+1, A( K+I, I ), 1 )
      WA = SIGN( WN, A( K+I, I ) )
      IF( WN == ZERO ) THEN
         TAU = ZERO
      ELSE
         WB = A( K+I, I ) + WA
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DSCAL( N-K-I, ONE / WB, A( K+I+1, I ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         A( K+I, I ) = ONE
         TAU = WB / WA
      END IF
!
!        apply reflection to A(k+i:n,i+1:k+i-1) from the left
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGEMV( 'Transpose', N-K-I+1, K-1, ONE, A( K+I, I+1 ), LDA, &
                  A( K+I, I ), 1, ZERO, WORK, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGEMV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGER( N-K-I+1, K-1, -TAU, A( K+I, I ), 1, WORK, 1, &
                 A( K+I, I+1 ), LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGER : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        apply reflection to A(k+i:n,k+i:n) from the left and the right
!
!        compute  y := tau * A * u
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DSYMV( 'Lower', N-K-I+1, TAU, A( K+I, K+I ), LDA, &
                  A( K+I, I ), 1, ZERO, WORK, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DSYMV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        compute  v := y - 1/2 * tau * ( y, u ) * u
!
      ALPHA = -HALF*TAU*DDOT( N-K-I+1, WORK, 1, A( K+I, I ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DAXPY( N-K-I+1, ALPHA, A( K+I, I ), 1, WORK, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DAXPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        apply symmetric rank-2 update to A(k+i:n,k+i:n)
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DSYR2( 'Lower', N-K-I+1, -ONE, A( K+I, I ), 1, WORK, 1, &
                  A( K+I, K+I ), LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DSYR2 : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
      A( K+I, I ) = -WA
      DO J = K + I + 1, N
         A( J, I ) = ZERO
      ENDDO
   ENDDO
!
!     Store full symmetric matrix
!
   DO J = 1, N
      DO I = J + 1, N
         A( J, I ) = A( I, J )
      ENDDO
   ENDDO
   RETURN
!
!     End of DLAGSY
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

