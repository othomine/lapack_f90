!> \brief \b ZLAGGE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAGGE( M, N, KL, KU, D, A, LDA, ISEED, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, KL, KU, LDA, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       DOUBLE PRECISION   D( * )
!       COMPLEX*16         A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAGGE generates a complex general m by n matrix A, by pre- and post-
!> multiplying a real diagonal matrix D with random unitary matrices:
!> A = U*D*V. The lower and upper bandwidths may then be reduced to
!> kl and ku by additional unitary transformations.
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
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of nonzero subdiagonals within the band of A.
!>          0 <= KL <= M-1.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of nonzero superdiagonals within the band of A.
!>          0 <= KU <= N-1.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (min(M,N))
!>          The diagonal elements of the diagonal matrix D.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The generated m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= M.
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
!>          WORK is COMPLEX*16 array, dimension (M+N)
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
   SUBROUTINE ZLAGGE( M, N, KL, KU, D, A, LDA, ISEED, WORK, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, KL, KU, LDA, M, N
!     ..
!     .. Array Arguments ..
   INTEGER            ISEED( 4 )
   DOUBLE PRECISION   D( * )
   COMPLEX*16         A( LDA, * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, J
   DOUBLE PRECISION   WN
   COMPLEX*16         TAU, WA, WB
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA, ZGEMV, ZGERC, ZLACGV, ZLARNV, ZSCAL
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
   IF( M < 0 ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( KL < 0 .OR. KL > M-1 ) THEN
      INFO = -3
   ELSE IF( KU < 0 .OR. KU > N-1 ) THEN
      INFO = -4
   ELSE IF( LDA < MAX( 1, M ) ) THEN
      INFO = -7
   END IF
   IF( INFO < 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'ZLAGGE', -INFO )
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
!     initialize A to diagonal matrix
!
   A(1:M,1:N) = (0.0D+0,0.0D+0)
   FORALL (I = 1:MIN( M, N )) A( I, I ) = D( I )
!
!     Quick exit if the user wants a diagonal matrix
!
   IF(( KL  ==  0 ).AND.( KU  ==  0)) RETURN
!
!     pre- and post-multiply A by random unitary matrices
!
   DO I = MIN( M, N ), 1, -1
      IF( I < M ) THEN
!
!           generate random reflection
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZLARNV( 3, ISEED, M-I+1, WORK )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZLARNV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         WN = DZNRM2( M-I+1, WORK, 1 )
         WA = ( WN / ABS( WORK( 1 ) ) )*WORK( 1 )
         IF( WN == (0.0D+0,0.0D+0) ) THEN
            TAU = (0.0D+0,0.0D+0)
         ELSE
            WB = WORK( 1 ) + WA
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZSCAL( M-I, (1.0D+0,0.0D+0) / WB, WORK( 2 ), 1 )
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
!           multiply A(i:m,i:n) by random reflection from the left
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZGEMV( 'Conjugate transpose', M-I+1, N-I+1, (1.0D+0,0.0D+0), &
                     A( I, I ), LDA, WORK, 1, (0.0D+0,0.0D+0), WORK( M+1 ), 1 )
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
         CALL ZGERC( M-I+1, N-I+1, -TAU, WORK, 1, WORK( M+1 ), 1, &
                     A( I, I ), LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZGERC : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      END IF
      IF( I < N ) THEN
!
!           generate random reflection
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
!           multiply A(i:m,i:n) by random reflection from the right
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZGEMV( 'No transpose', M-I+1, N-I+1, (1.0D+0,0.0D+0), A( I, I ), &
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
         CALL ZGERC( M-I+1, N-I+1, -TAU, WORK( N+1 ), 1, WORK, 1, &
                     A( I, I ), LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZGERC : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      END IF
   ENDDO
!
!     Reduce number of subdiagonals to KL and number of superdiagonals
!     to KU
!
   DO I = 1, MAX( M-1-KL, N-1-KU )
      IF( KL <= KU ) THEN
!
!           annihilate subdiagonal elements first (necessary if KL = 0)
!
         IF( I <= MIN( M-1-KL, N ) ) THEN
!
!              generate reflection to annihilate A(kl+i+1:m,i)
!
            WN = DZNRM2( M-KL-I+1, A( KL+I, I ), 1 )
            WA = ( WN / ABS( A( KL+I, I ) ) )*A( KL+I, I )
            IF( WN == (0.0D+0,0.0D+0) ) THEN
               TAU = (0.0D+0,0.0D+0)
            ELSE
               WB = A( KL+I, I ) + WA
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZSCAL( M-KL-I, (1.0D+0,0.0D+0) / WB, A( KL+I+1, I ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZSCAL : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               A( KL+I, I ) = (1.0D+0,0.0D+0)
               TAU = DBLE( WB / WA )
            END IF
!
!              apply reflection to A(kl+i:m,i+1:n) from the left
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZGEMV( 'Conjugate transpose', M-KL-I+1, N-I, (1.0D+0,0.0D+0), &
                        A( KL+I, I+1 ), LDA, A( KL+I, I ), 1, (0.0D+0,0.0D+0), &
                        WORK, 1 )
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
            CALL ZGERC( M-KL-I+1, N-I, -TAU, A( KL+I, I ), 1, WORK, &
                        1, A( KL+I, I+1 ), LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZGERC : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            A( KL+I, I ) = -WA
         END IF
!
         IF( I <= MIN( N-1-KU, M ) ) THEN
!
!              generate reflection to annihilate A(i,ku+i+1:n)
!
            WN = DZNRM2( N-KU-I+1, A( I, KU+I ), LDA )
            WA = ( WN / ABS( A( I, KU+I ) ) )*A( I, KU+I )
            IF( WN == (0.0D+0,0.0D+0) ) THEN
               TAU = (0.0D+0,0.0D+0)
            ELSE
               WB = A( I, KU+I ) + WA
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZSCAL( N-KU-I, (1.0D+0,0.0D+0) / WB, A( I, KU+I+1 ), LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZSCAL : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               A( I, KU+I ) = (1.0D+0,0.0D+0)
               TAU = DBLE( WB / WA )
            END IF
!
!              apply reflection to A(i+1:m,ku+i:n) from the right
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZLACGV( N-KU-I+1, A( I, KU+I ), LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZLACGV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZGEMV( 'No transpose', M-I, N-KU-I+1, (1.0D+0,0.0D+0), &
                        A( I+1, KU+I ), LDA, A( I, KU+I ), LDA, (0.0D+0,0.0D+0), &
                        WORK, 1 )
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
            CALL ZGERC( M-I, N-KU-I+1, -TAU, WORK, 1, A( I, KU+I ), &
                        LDA, A( I+1, KU+I ), LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZGERC : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            A( I, KU+I ) = -WA
         END IF
      ELSE
!
!           annihilate superdiagonal elements first (necessary if
!           KU = 0)
!
         IF( I <= MIN( N-1-KU, M ) ) THEN
!
!              generate reflection to annihilate A(i,ku+i+1:n)
!
            WN = DZNRM2( N-KU-I+1, A( I, KU+I ), LDA )
            WA = ( WN / ABS( A( I, KU+I ) ) )*A( I, KU+I )
            IF( WN == (0.0D+0,0.0D+0) ) THEN
               TAU = (0.0D+0,0.0D+0)
            ELSE
               WB = A( I, KU+I ) + WA
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZSCAL( N-KU-I, (1.0D+0,0.0D+0) / WB, A( I, KU+I+1 ), LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZSCAL : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               A( I, KU+I ) = (1.0D+0,0.0D+0)
               TAU = DBLE( WB / WA )
            END IF
!
!              apply reflection to A(i+1:m,ku+i:n) from the right
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZLACGV( N-KU-I+1, A( I, KU+I ), LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZLACGV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZGEMV( 'No transpose', M-I, N-KU-I+1, (1.0D+0,0.0D+0), &
                        A( I+1, KU+I ), LDA, A( I, KU+I ), LDA, (0.0D+0,0.0D+0), &
                        WORK, 1 )
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
            CALL ZGERC( M-I, N-KU-I+1, -TAU, WORK, 1, A( I, KU+I ), &
                        LDA, A( I+1, KU+I ), LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZGERC : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            A( I, KU+I ) = -WA
         END IF
!
         IF( I <= MIN( M-1-KL, N ) ) THEN
!
!              generate reflection to annihilate A(kl+i+1:m,i)
!
            WN = DZNRM2( M-KL-I+1, A( KL+I, I ), 1 )
            WA = ( WN / ABS( A( KL+I, I ) ) )*A( KL+I, I )
            IF( WN == (0.0D+0,0.0D+0) ) THEN
               TAU = (0.0D+0,0.0D+0)
            ELSE
               WB = A( KL+I, I ) + WA
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZSCAL( M-KL-I, (1.0D+0,0.0D+0) / WB, A( KL+I+1, I ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZSCAL : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               A( KL+I, I ) = (1.0D+0,0.0D+0)
               TAU = DBLE( WB / WA )
            END IF
!
!              apply reflection to A(kl+i:m,i+1:n) from the left
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZGEMV( 'Conjugate transpose', M-KL-I+1, N-I, (1.0D+0,0.0D+0), &
                        A( KL+I, I+1 ), LDA, A( KL+I, I ), 1, (0.0D+0,0.0D+0), &
                        WORK, 1 )
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
            CALL ZGERC( M-KL-I+1, N-I, -TAU, A( KL+I, I ), 1, WORK, &
                        1, A( KL+I, I+1 ), LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZGERC : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            A( KL+I, I ) = -WA
         END IF
      END IF
!
      IF (I  <=  N) A(KL+I+1:M, I ) = (0.0D+0,0.0D+0)
!
      IF (I  <=  M) A( I,KU+I+1:N) = (0.0D+0,0.0D+0)
   ENDDO
   RETURN
!
!     End of ZLAGGE
!
END

