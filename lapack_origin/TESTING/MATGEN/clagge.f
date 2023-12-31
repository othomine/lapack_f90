*> \brief \b CLAGGE
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CLAGGE( M, N, KL, KU, D, A, LDA, ISEED, WORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, KL, KU, LDA, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            ISEED( 4 )
*       REAL               D( * )
*       COMPLEX            A( LDA, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CLAGGE generates a complex general m by n matrix A, by pre- and post-
*> multiplying a real diagonal matrix D with random unitary matrices:
*> A = U*D*V. The lower and upper bandwidths may then be reduced to
*> kl and ku by additional unitary transformations.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] KL
*> \verbatim
*>          KL is INTEGER
*>          The number of nonzero subdiagonals within the band of A.
*>          0 <= KL <= M-1.
*> \endverbatim
*>
*> \param[in] KU
*> \verbatim
*>          KU is INTEGER
*>          The number of nonzero superdiagonals within the band of A.
*>          0 <= KU <= N-1.
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is REAL array, dimension (min(M,N))
*>          The diagonal elements of the diagonal matrix D.
*> \endverbatim
*>
*> \param[out] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          The generated m by n matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= M.
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
*>          WORK is COMPLEX array, dimension (M+N)
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
      SUBROUTINE CLAGGE( M, N, KL, KU, D, A, LDA, ISEED, WORK, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      REAL               D( * )
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
      INTEGER            I, J
      REAL               WN
      COMPLEX            TAU, WA, WB
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMV, CGERC, CLACGV, CLARNV, CSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, REAL
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
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 .OR. KL.GT.M-1 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 .OR. KU.GT.N-1 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -7
      END IF
      IF( INFO.LT.0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL XERBLA( 'CLAGGE', -INFO )
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
*     initialize A to diagonal matrix
*
      DO 20 J = 1, N
         DO 10 I = 1, M
            A( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
      DO 30 I = 1, MIN( M, N )
         A( I, I ) = D( I )
   30 CONTINUE
*
*     Quick exit if the user wants a diagonal matrix
*
      IF(( KL .EQ. 0 ).AND.( KU .EQ. 0)) RETURN
*
*     pre- and post-multiply A by random unitary matrices
*
      DO 40 I = MIN( M, N ), 1, -1
         IF( I.LT.M ) THEN
*
*           generate random reflection
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
            CALL CLARNV( 3, ISEED, M-I+1, WORK )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CLARNV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
            WN = SCNRM2( M-I+1, WORK, 1 )
            WA = ( WN / ABS( WORK( 1 ) ) )*WORK( 1 )
            IF( WN.EQ.ZERO ) THEN
               TAU = ZERO
            ELSE
               WB = WORK( 1 ) + WA
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
               CALL CSCAL( M-I, ONE / WB, WORK( 2 ), 1 )
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
*           multiply A(i:m,i:n) by random reflection from the left
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
            CALL CGEMV( 'Conjugate transpose', M-I+1, N-I+1, ONE,
     $                  A( I, I ), LDA, WORK, 1, ZERO, WORK( M+1 ), 1 )
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
            CALL CGERC( M-I+1, N-I+1, -TAU, WORK, 1, WORK( M+1 ), 1,
     $                  A( I, I ), LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CGERC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         END IF
         IF( I.LT.N ) THEN
*
*           generate random reflection
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
*           multiply A(i:m,i:n) by random reflection from the right
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
            CALL CGEMV( 'No transpose', M-I+1, N-I+1, ONE, A( I, I ),
     $                  LDA, WORK, 1, ZERO, WORK( N+1 ), 1 )
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
            CALL CGERC( M-I+1, N-I+1, -TAU, WORK( N+1 ), 1, WORK, 1,
     $                  A( I, I ), LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CGERC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         END IF
   40 CONTINUE
*
*     Reduce number of subdiagonals to KL and number of superdiagonals
*     to KU
*
      DO 70 I = 1, MAX( M-1-KL, N-1-KU )
         IF( KL.LE.KU ) THEN
*
*           annihilate subdiagonal elements first (necessary if KL = 0)
*
            IF( I.LE.MIN( M-1-KL, N ) ) THEN
*
*              generate reflection to annihilate A(kl+i+1:m,i)
*
               WN = SCNRM2( M-KL-I+1, A( KL+I, I ), 1 )
               WA = ( WN / ABS( A( KL+I, I ) ) )*A( KL+I, I )
               IF( WN.EQ.ZERO ) THEN
                  TAU = ZERO
               ELSE
                  WB = A( KL+I, I ) + WA
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                  CALL CSCAL( M-KL-I, ONE / WB, A( KL+I+1, I ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CSCAL : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
                  A( KL+I, I ) = ONE
                  TAU = REAL( WB / WA )
               END IF
*
*              apply reflection to A(kl+i:m,i+1:n) from the left
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
               CALL CGEMV( 'Conjugate transpose', M-KL-I+1, N-I, ONE,
     $                     A( KL+I, I+1 ), LDA, A( KL+I, I ), 1, ZERO,
     $                     WORK, 1 )
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
               CALL CGERC( M-KL-I+1, N-I, -TAU, A( KL+I, I ), 1, WORK,
     $                     1, A( KL+I, I+1 ), LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CGERC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
               A( KL+I, I ) = -WA
            END IF
*
            IF( I.LE.MIN( N-1-KU, M ) ) THEN
*
*              generate reflection to annihilate A(i,ku+i+1:n)
*
               WN = SCNRM2( N-KU-I+1, A( I, KU+I ), LDA )
               WA = ( WN / ABS( A( I, KU+I ) ) )*A( I, KU+I )
               IF( WN.EQ.ZERO ) THEN
                  TAU = ZERO
               ELSE
                  WB = A( I, KU+I ) + WA
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                  CALL CSCAL( N-KU-I, ONE / WB, A( I, KU+I+1 ), LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CSCAL : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
                  A( I, KU+I ) = ONE
                  TAU = REAL( WB / WA )
               END IF
*
*              apply reflection to A(i+1:m,ku+i:n) from the right
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
               CALL CLACGV( N-KU-I+1, A( I, KU+I ), LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CLACGV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
               CALL CGEMV( 'No transpose', M-I, N-KU-I+1, ONE,
     $                     A( I+1, KU+I ), LDA, A( I, KU+I ), LDA, ZERO,
     $                     WORK, 1 )
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
               CALL CGERC( M-I, N-KU-I+1, -TAU, WORK, 1, A( I, KU+I ),
     $                     LDA, A( I+1, KU+I ), LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CGERC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
               A( I, KU+I ) = -WA
            END IF
         ELSE
*
*           annihilate superdiagonal elements first (necessary if
*           KU = 0)
*
            IF( I.LE.MIN( N-1-KU, M ) ) THEN
*
*              generate reflection to annihilate A(i,ku+i+1:n)
*
               WN = SCNRM2( N-KU-I+1, A( I, KU+I ), LDA )
               WA = ( WN / ABS( A( I, KU+I ) ) )*A( I, KU+I )
               IF( WN.EQ.ZERO ) THEN
                  TAU = ZERO
               ELSE
                  WB = A( I, KU+I ) + WA
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                  CALL CSCAL( N-KU-I, ONE / WB, A( I, KU+I+1 ), LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CSCAL : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
                  A( I, KU+I ) = ONE
                  TAU = REAL( WB / WA )
               END IF
*
*              apply reflection to A(i+1:m,ku+i:n) from the right
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
               CALL CLACGV( N-KU-I+1, A( I, KU+I ), LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CLACGV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
               CALL CGEMV( 'No transpose', M-I, N-KU-I+1, ONE,
     $                     A( I+1, KU+I ), LDA, A( I, KU+I ), LDA, ZERO,
     $                     WORK, 1 )
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
               CALL CGERC( M-I, N-KU-I+1, -TAU, WORK, 1, A( I, KU+I ),
     $                     LDA, A( I+1, KU+I ), LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CGERC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
               A( I, KU+I ) = -WA
            END IF
*
            IF( I.LE.MIN( M-1-KL, N ) ) THEN
*
*              generate reflection to annihilate A(kl+i+1:m,i)
*
               WN = SCNRM2( M-KL-I+1, A( KL+I, I ), 1 )
               WA = ( WN / ABS( A( KL+I, I ) ) )*A( KL+I, I )
               IF( WN.EQ.ZERO ) THEN
                  TAU = ZERO
               ELSE
                  WB = A( KL+I, I ) + WA
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                  CALL CSCAL( M-KL-I, ONE / WB, A( KL+I+1, I ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CSCAL : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
                  A( KL+I, I ) = ONE
                  TAU = REAL( WB / WA )
               END IF
*
*              apply reflection to A(kl+i:m,i+1:n) from the left
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
               CALL CGEMV( 'Conjugate transpose', M-KL-I+1, N-I, ONE,
     $                     A( KL+I, I+1 ), LDA, A( KL+I, I ), 1, ZERO,
     $                     WORK, 1 )
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
               CALL CGERC( M-KL-I+1, N-I, -TAU, A( KL+I, I ), 1, WORK,
     $                     1, A( KL+I, I+1 ), LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CGERC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
               A( KL+I, I ) = -WA
            END IF
         END IF
*
         IF (I .LE. N) THEN
            DO 50 J = KL + I + 1, M
               A( J, I ) = ZERO
   50       CONTINUE
         END IF
*
         IF (I .LE. M) THEN
            DO 60 J = KU + I + 1, N
               A( I, J ) = ZERO
   60       CONTINUE
         END IF
   70 CONTINUE
      RETURN
*
*     End of CLAGGE
*
      END

