!> \brief \b CLAGSY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAGSY( N, K, D, A, LDA, ISEED, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       REAL               D( * )
!       COMPLEX            A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAGSY generates a complex symmetric matrix A, by pre- and post-
!> multiplying a real diagonal matrix D with a random unitary matrix:
!> A = U*D*U**T. The semi-bandwidth may then be reduced to k by
!> additional unitary transformations.
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
!>          D is REAL array, dimension (N)
!>          The diagonal elements of the diagonal matrix D.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
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
!>          WORK is COMPLEX array, dimension (2*N)
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
!> \ingroup complex_matgen
!
!  =====================================================================
   SUBROUTINE CLAGSY( N, K, D, A, LDA, ISEED, WORK, INFO )
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
   REAL               D( * )
   COMPLEX            A( LDA, * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, II, J, JJ
   REAL               WN
   COMPLEX            ALPHA, TAU, WA, WB
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Subroutines ..
   EXTERNAL           CAXPY, CGEMV, CGERC, CLACGV, CLARNV, CSCAL, &
                      CSYMV, XERBLA
!     ..
!     .. External Functions ..
   REAL               SCNRM2
   COMPLEX            CDOTC
   EXTERNAL           SCNRM2, CDOTC
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
      CALL XERBLA( 'CLAGSY', -INFO )
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
      A(J+1:N, J ) = (0.0E+0,0.0E+0)
   ENDDO
   FORALL (I = 1:N) A( I, I ) = D( I )
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
      CALL CLARNV( 3, ISEED, N-I+1, WORK )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      WN = SCNRM2( N-I+1, WORK, 1 )
      WA = ( WN / ABS( WORK( 1 ) ) )*WORK( 1 )
      IF( WN == (0.0E+0,0.0E+0) ) THEN
         TAU = (0.0E+0,0.0E+0)
      ELSE
         WB = WORK( 1 ) + WA
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL CSCAL( N-I, (1.0E+0,0.0E+0) / WB, WORK( 2 ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CSCAL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         WORK( 1 ) = (1.0E+0,0.0E+0)
         TAU = REAL( WB / WA )
      END IF
!
!        apply random reflection to A(i:n,i:n) from the left
!        and the right
!
!        compute  y := tau * A * conjg(u)
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CLACGV( N-I+1, WORK, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CLACGV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSYMV( 'Lower', N-I+1, TAU, A( I, I ), LDA, WORK, 1, (0.0E+0,0.0E+0), &
                  WORK( N+1 ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSYMV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CLACGV( N-I+1, WORK, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CLACGV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        compute  v := y - 1/2 * tau * ( u, y ) * u
!
      ALPHA = -(0.5E+0,0.0E+0)*TAU*CDOTC( N-I+1, WORK, 1, WORK( N+1 ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CAXPY( N-I+1, ALPHA, WORK, 1, WORK( N+1 ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CAXPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        apply the transformation as a rank-2 update to A(i:n,i:n)
!
!        CALL CSYR2( 'Lower', N-I+1, -(1.0E+0,0.0E+0), WORK, 1, WORK( N+1 ), 1,
!        $               A( I, I ), LDA )
!
      DO JJ = I, N
         A(JJ:N, JJ ) = A(JJ:N, JJ ) - &
                       WORK( 1-I+JJ:1-I+N )*WORK( N+JJ-I+1 ) - &
                       WORK( N-I+1+JJ:2*N-I+1 )*WORK( JJ-I+1 )
      ENDDO
   ENDDO
!
!     Reduce number of subdiagonals to K
!
   DO I = 1, N - 1 - K
!
!        generate reflection to annihilate A(k+i+1:n,i)
!
      WN = SCNRM2( N-K-I+1, A( K+I, I ), 1 )
      WA = ( WN / ABS( A( K+I, I ) ) )*A( K+I, I )
      IF( WN == (0.0E+0,0.0E+0) ) THEN
         TAU = (0.0E+0,0.0E+0)
      ELSE
         WB = A( K+I, I ) + WA
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL CSCAL( N-K-I, (1.0E+0,0.0E+0) / WB, A( K+I+1, I ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CSCAL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         A( K+I, I ) = (1.0E+0,0.0E+0)
         TAU = REAL( WB / WA )
      END IF
!
!        apply reflection to A(k+i:n,i+1:k+i-1) from the left
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CGEMV( 'Conjugate transpose', N-K-I+1, K-1, (1.0E+0,0.0E+0), &
                  A( K+I, I+1 ), LDA, A( K+I, I ), 1, (0.0E+0,0.0E+0), WORK, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CGEMV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CGERC( N-K-I+1, K-1, -TAU, A( K+I, I ), 1, WORK, 1, &
                  A( K+I, I+1 ), LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CGERC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        apply reflection to A(k+i:n,k+i:n) from the left and the right
!
!        compute  y := tau * A * conjg(u)
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CLACGV( N-K-I+1, A( K+I, I ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CLACGV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSYMV( 'Lower', N-K-I+1, TAU, A( K+I, K+I ), LDA, &
                  A( K+I, I ), 1, (0.0E+0,0.0E+0), WORK, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSYMV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CLACGV( N-K-I+1, A( K+I, I ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CLACGV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        compute  v := y - 1/2 * tau * ( u, y ) * u
!
      ALPHA = -(0.5E+0,0.0E+0)*TAU*CDOTC( N-K-I+1, A( K+I, I ), 1, WORK, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CAXPY( N-K-I+1, ALPHA, A( K+I, I ), 1, WORK, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CAXPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        apply symmetric rank-2 update to A(k+i:n,k+i:n)
!
!        CALL CSYR2( 'Lower', N-K-I+1, -(1.0E+0,0.0E+0), A( K+I, I ), 1, WORK, 1,
!        $               A( K+I, K+I ), LDA )
!
      DO JJ = K + I, N
         A(JJ:N, JJ ) = A(JJ:N, JJ ) - A(JJ:N, I )*WORK( JJ-K-I+1 ) - &
                       WORK( 1-K-I+JJ:1-K-I+N )*A( JJ, I )
      ENDDO
!
      A( K+I, I ) = -WA
      A(K+I+1:N, I ) = (0.0E+0,0.0E+0)
   ENDDO
!
!     Store full symmetric matrix
!
   DO J = 1, N
      A( J,J+1:N) = A(J+1:N, J )
   ENDDO
   RETURN
!
!     End of CLAGSY
!
END

