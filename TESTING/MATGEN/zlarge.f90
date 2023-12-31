!> \brief \b ZLARGE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARGE( N, A, LDA, ISEED, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       COMPLEX*16         A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARGE pre- and post-multiplies a complex general n by n matrix A
!> with a random unitary matrix: A = U*D*U'.
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
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the original n by n matrix A.
!>          On exit, A is overwritten by U*A*U' for some random
!>          unitary matrix U.
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
!>          WORK is COMPLEX*16 array, dimension (2*N)
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
!> \ingroup complex16_matgen
!
!  =====================================================================
   SUBROUTINE ZLARGE( N, A, LDA, ISEED, WORK, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
   INTEGER            ISEED( 4 )
   COMPLEX*16         A( LDA, * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I
   DOUBLE PRECISION   WN
   COMPLEX*16         TAU, WA, WB
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA, ZGEMV, ZGERC, ZLARNV, ZSCAL
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DZNRM2
   EXTERNAL           DZNRM2
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
   INFO = 0
   IF( N < 0 ) THEN
      INFO = -1
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -3
   END IF
   IF( INFO < 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'ZLARGE', -INFO )
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
!     pre- and post-multiply A by random unitary matrix
!
   DO I = N, 1, -1
!
!        generate random reflection
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLARNV( 3, ISEED, N-I+1, WORK )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      WN = DZNRM2( N-I+1, WORK, 1 )
      WA = ( WN / ABS( WORK( 1 ) ) )*WORK( 1 )
      IF( WN == (0.0D+0,0.0D+0) ) THEN
         TAU = (0.0D+0,0.0D+0)
      ELSE
         WB = WORK( 1 ) + WA
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZSCAL( N-I, (1.0D+0,0.0D+0) / WB, WORK( 2 ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZSCAL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         WORK( 1 ) = (1.0D+0,0.0D+0)
         TAU = DBLE( WB / WA )
      END IF
!
!        multiply A(i:n,1:n) by random reflection from the left
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZGEMV( 'Conjugate transpose', N-I+1, N, (1.0D+0,0.0D+0), A( I, 1 ), &
                  LDA, WORK, 1, (0.0D+0,0.0D+0), WORK( N+1 ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZGEMV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZGERC( N-I+1, N, -TAU, WORK, 1, WORK( N+1 ), 1, A( I, 1 ), &
                  LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZGERC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        multiply A(1:n,i:n) by random reflection from the right
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZGEMV( 'No transpose', N, N-I+1, (1.0D+0,0.0D+0), A( 1, I ), LDA, &
                  WORK, 1, (0.0D+0,0.0D+0), WORK( N+1 ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZGEMV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZGERC( N, N-I+1, -TAU, WORK( N+1 ), 1, WORK, 1, A( 1, I ), &
                  LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZGERC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   ENDDO
   RETURN
!
!     End of ZLARGE
!
END

