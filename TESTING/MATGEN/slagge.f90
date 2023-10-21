!> \brief \b SLAGGE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAGGE( M, N, KL, KU, D, A, LDA, ISEED, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, KL, KU, LDA, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       REAL               A( LDA, * ), D( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAGGE generates a real general m by n matrix A, by pre- and post-
!> multiplying a real diagonal matrix D with random orthogonal matrices:
!> A = U*D*V. The lower and upper bandwidths may then be reduced to
!> kl and ku by additional orthogonal transformations.
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
!>          D is REAL array, dimension (min(M,N))
!>          The diagonal elements of the diagonal matrix D.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
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
!>          WORK is REAL array, dimension (M+N)
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
!> \ingroup real_matgen
!
!  =====================================================================
   SUBROUTINE SLAGGE( M, N, KL, KU, D, A, LDA, ISEED, WORK, INFO )
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
   REAL               A( LDA, * ), D( * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, J
   REAL               TAU, WA, WB, WN
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Subroutines ..
   EXTERNAL           SGEMV, SGER, SLARNV, SSCAL, XERBLA
!     ..
!     .. External Functions ..
   REAL               SNRM2
   EXTERNAL           SNRM2
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
      CALL XERBLA( 'SLAGGE', -INFO )
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
   A(1:M,1:N) = 0.0E+0
   FORALL (I = 1:MIN( M, N )) A( I, I ) = D( I )
!
!     Quick exit if the user wants a diagonal matrix
!
   IF(( KL  ==  0 ).AND.( KU  ==  0)) RETURN
!
!     pre- and post-multiply A by random orthogonal matrices
!
   DO I = MIN( M, N ), 1, -1
      IF( I < M ) THEN
!
!           generate random reflection
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SLARNV( 3, ISEED, M-I+1, WORK )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SLARNV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         WN = SNRM2( M-I+1, WORK, 1 )
         WA = SIGN( WN, WORK( 1 ) )
         IF( WN == 0.0E+0 ) THEN
            TAU = 0.0E+0
         ELSE
            WB = WORK( 1 ) + WA
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SSCAL( M-I, 1.0E+0 / WB, WORK( 2 ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SSCAL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            WORK( 1 ) = 1.0E+0
            TAU = WB / WA
         END IF
!
!           multiply A(i:m,i:n) by random reflection from the left
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SGEMV( 'Transpose', M-I+1, N-I+1, 1.0E+0, A( I, I ), LDA, &
                     WORK, 1, 0.0E+0, WORK( M+1 ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SGEMV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SGER( M-I+1, N-I+1, -TAU, WORK, 1, WORK( M+1 ), 1, &
                    A( I, I ), LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SGER : ',&
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
         CALL SLARNV( 3, ISEED, N-I+1, WORK )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SLARNV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         WN = SNRM2( N-I+1, WORK, 1 )
         WA = SIGN( WN, WORK( 1 ) )
         IF( WN == 0.0E+0 ) THEN
            TAU = 0.0E+0
         ELSE
            WB = WORK( 1 ) + WA
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SSCAL( N-I, 1.0E+0 / WB, WORK( 2 ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SSCAL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            WORK( 1 ) = 1.0E+0
            TAU = WB / WA
         END IF
!
!           multiply A(i:m,i:n) by random reflection from the right
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SGEMV( 'No transpose', M-I+1, N-I+1, 1.0E+0, A( I, I ), &
                     LDA, WORK, 1, 0.0E+0, WORK( N+1 ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SGEMV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SGER( M-I+1, N-I+1, -TAU, WORK( N+1 ), 1, WORK, 1, &
                    A( I, I ), LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SGER : ',&
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
            WN = SNRM2( M-KL-I+1, A( KL+I, I ), 1 )
            WA = SIGN( WN, A( KL+I, I ) )
            IF( WN == 0.0E+0 ) THEN
               TAU = 0.0E+0
            ELSE
               WB = A( KL+I, I ) + WA
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SSCAL( M-KL-I, 1.0E+0 / WB, A( KL+I+1, I ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SSCAL : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               A( KL+I, I ) = 1.0E+0
               TAU = WB / WA
            END IF
!
!              apply reflection to A(kl+i:m,i+1:n) from the left
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGEMV( 'Transpose', M-KL-I+1, N-I, 1.0E+0, &
                        A( KL+I, I+1 ), LDA, A( KL+I, I ), 1, 0.0E+0, &
                        WORK, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGEMV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGER( M-KL-I+1, N-I, -TAU, A( KL+I, I ), 1, WORK, 1, &
                       A( KL+I, I+1 ), LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGER : ',&
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
            WN = SNRM2( N-KU-I+1, A( I, KU+I ), LDA )
            WA = SIGN( WN, A( I, KU+I ) )
            IF( WN == 0.0E+0 ) THEN
               TAU = 0.0E+0
            ELSE
               WB = A( I, KU+I ) + WA
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SSCAL( N-KU-I, 1.0E+0 / WB, A( I, KU+I+1 ), LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SSCAL : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               A( I, KU+I ) = 1.0E+0
               TAU = WB / WA
            END IF
!
!              apply reflection to A(i+1:m,ku+i:n) from the right
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGEMV( 'No transpose', M-I, N-KU-I+1, 1.0E+0, &
                        A( I+1, KU+I ), LDA, A( I, KU+I ), LDA, 0.0E+0, &
                        WORK, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGEMV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGER( M-I, N-KU-I+1, -TAU, WORK, 1, A( I, KU+I ), &
                       LDA, A( I+1, KU+I ), LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGER : ',&
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
            WN = SNRM2( N-KU-I+1, A( I, KU+I ), LDA )
            WA = SIGN( WN, A( I, KU+I ) )
            IF( WN == 0.0E+0 ) THEN
               TAU = 0.0E+0
            ELSE
               WB = A( I, KU+I ) + WA
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SSCAL( N-KU-I, 1.0E+0 / WB, A( I, KU+I+1 ), LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SSCAL : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               A( I, KU+I ) = 1.0E+0
               TAU = WB / WA
            END IF
!
!              apply reflection to A(i+1:m,ku+i:n) from the right
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGEMV( 'No transpose', M-I, N-KU-I+1, 1.0E+0, &
                        A( I+1, KU+I ), LDA, A( I, KU+I ), LDA, 0.0E+0, &
                        WORK, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGEMV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGER( M-I, N-KU-I+1, -TAU, WORK, 1, A( I, KU+I ), &
                       LDA, A( I+1, KU+I ), LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGER : ',&
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
            WN = SNRM2( M-KL-I+1, A( KL+I, I ), 1 )
            WA = SIGN( WN, A( KL+I, I ) )
            IF( WN == 0.0E+0 ) THEN
               TAU = 0.0E+0
            ELSE
               WB = A( KL+I, I ) + WA
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SSCAL( M-KL-I, 1.0E+0 / WB, A( KL+I+1, I ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SSCAL : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               A( KL+I, I ) = 1.0E+0
               TAU = WB / WA
            END IF
!
!              apply reflection to A(kl+i:m,i+1:n) from the left
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGEMV( 'Transpose', M-KL-I+1, N-I, 1.0E+0, &
                        A( KL+I, I+1 ), LDA, A( KL+I, I ), 1, 0.0E+0, &
                        WORK, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGEMV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGER( M-KL-I+1, N-I, -TAU, A( KL+I, I ), 1, WORK, 1, &
                       A( KL+I, I+1 ), LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGER : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            A( KL+I, I ) = -WA
         END IF
      END IF
!
      IF (I  <=  N) A(KL+I+1:M, I ) = 0.0E+0
!
      IF (I  <=  M) A( I,KU+I+1:N) = 0.0E+0
   ENDDO
   RETURN
!
!     End of SLAGGE
!
END

