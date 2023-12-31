*> \brief \b CLARGE
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CLARGE( N, A, LDA, ISEED, WORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            ISEED( 4 )
*       COMPLEX            A( LDA, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CLARGE pre- and post-multiplies a complex general n by n matrix A
*> with a random unitary matrix: A = U*D*U'.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          On entry, the original n by n matrix A.
*>          On exit, A is overwritten by U*A*U' for some random
*>          unitary matrix U.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= N.
*> \endverbatim
*>
*> \param[in,out] ISEED
*> \verbatim
*>          ISEED is INTEGER array, dimension (4)
*>          On entry, the seed of the random number generator; the array
*>          elements must be between 0 and 4095, and ISEED(4) must be
*>          odd.
*>          On exit, the seed is updated.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (2*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complex_matgen
*
*  =====================================================================
      SUBROUTINE CLARGE( N, A, LDA, ISEED, WORK, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      COMPLEX            A( LDA, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ),
     $                   ONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      REAL               WN
      COMPLEX            TAU, WA, WB
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMV, CGERC, CLARNV, CSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, REAL
*     ..
*     .. External Functions ..
      REAL               SCNRM2
      EXTERNAL           SCNRM2
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -3
      END IF
      IF( INFO.LT.0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL XERBLA( 'CLARGE', -INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : XERBLA : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         RETURN
      END IF
*
*     pre- and post-multiply A by random unitary matrix
*
      DO 10 I = N, 1, -1
*
*        generate random reflection
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL CLARNV( 3, ISEED, N-I+1, WORK )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CLARNV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         WN = SCNRM2( N-I+1, WORK, 1 )
         WA = ( WN / ABS( WORK( 1 ) ) )*WORK( 1 )
         IF( WN.EQ.ZERO ) THEN
            TAU = ZERO
         ELSE
            WB = WORK( 1 ) + WA
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
            CALL CSCAL( N-I, ONE / WB, WORK( 2 ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CSCAL : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
            WORK( 1 ) = ONE
            TAU = REAL( WB / WA )
         END IF
*
*        multiply A(i:n,1:n) by random reflection from the left
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL CGEMV( 'Conjugate transpose', N-I+1, N, ONE, A( I, 1 ),
     $               LDA, WORK, 1, ZERO, WORK( N+1 ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CGEMV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL CGERC( N-I+1, N, -TAU, WORK, 1, WORK( N+1 ), 1, A( I, 1 ),
     $               LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CGERC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
*
*        multiply A(1:n,i:n) by random reflection from the right
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL CGEMV( 'No transpose', N, N-I+1, ONE, A( 1, I ), LDA,
     $               WORK, 1, ZERO, WORK( N+1 ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CGEMV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL CGERC( N, N-I+1, -TAU, WORK( N+1 ), 1, WORK, 1, A( 1, I ),
     $               LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CGERC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   10 CONTINUE
      RETURN
*
*     End of CLARGE
*
      END

