!> \brief \b SLATM7
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLATM7( MODE, COND, IRSIGN, IDIST, ISEED, D, N,
!                          RANK, INFO )
!
!       .. Scalar Arguments ..
!       REAL               COND
!       INTEGER            IDIST, INFO, IRSIGN, MODE, N, RANK
!       ..
!       .. Array Arguments ..
!       REAL               D( * )
!       INTEGER            ISEED( 4 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SLATM7 computes the entries of D as specified by MODE
!>    COND and IRSIGN. IDIST and ISEED determine the generation
!>    of random numbers. SLATM7 is called by SLATMT to generate
!>    random test matrices.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \verbatim
!>  MODE   - INTEGER
!>           On entry describes how D is to be computed:
!>           MODE = 0 means do not change D.
!>
!>           MODE = 1 sets D(1)=1 and D(2:RANK)=1.0/COND
!>           MODE = 2 sets D(1:RANK-1)=1 and D(RANK)=1.0/COND
!>           MODE = 3 sets D(I)=COND**(-(I-1)/(RANK-1)) I=1:RANK
!>
!>           MODE = 4 sets D(i)=1 - (i-1)/(N-1)*(1 - 1/COND)
!>           MODE = 5 sets D to random numbers in the range
!>                    ( 1/COND , 1 ) such that their logarithms
!>                    are uniformly distributed.
!>           MODE = 6 set D to random numbers from same distribution
!>                    as the rest of the matrix.
!>           MODE < 0 has the same meaning as ABS(MODE), except that
!>              the order of the elements of D is reversed.
!>           Thus if MODE is positive, D has entries ranging from
!>              1 to 1/COND, if negative, from 1/COND to 1,
!>           Not modified.
!>
!>  COND   - REAL
!>           On entry, used as described under MODE above.
!>           If used, it must be >= 1. Not modified.
!>
!>  IRSIGN - INTEGER
!>           On entry, if MODE neither -6, 0 nor 6, determines sign of
!>           entries of D
!>           0 => leave entries of D unchanged
!>           1 => multiply each entry of D by 1 or -1 with probability .5
!>
!>  IDIST  - CHARACTER*1
!>           On entry, IDIST specifies the type of distribution to be
!>           used to generate a random matrix .
!>           1 => UNIFORM( 0, 1 )
!>           2 => UNIFORM( -1, 1 )
!>           3 => NORMAL( 0, 1 )
!>           Not modified.
!>
!>  ISEED  - INTEGER array, dimension ( 4 )
!>           On entry ISEED specifies the seed of the random number
!>           generator. The random number generator uses a
!>           linear congruential sequence limited to small
!>           integers, and so should produce machine independent
!>           random numbers. The values of ISEED are changed on
!>           exit, and can be used in the next call to SLATM7
!>           to continue the same random number sequence.
!>           Changed on exit.
!>
!>  D      - REAL array, dimension ( MIN( M , N ) )
!>           Array to be computed according to MODE, COND and IRSIGN.
!>           May be changed on exit if MODE is nonzero.
!>
!>  N      - INTEGER
!>           Number of entries of D. Not modified.
!>
!>  RANK   - INTEGER
!>           The rank of matrix to be generated for modes 1,2,3 only.
!>           D( RANK+1:N ) = 0.
!>           Not modified.
!>
!>  INFO   - INTEGER
!>            0  => normal termination
!>           -1  => if MODE not in range -6 to 6
!>           -2  => if MODE neither -6, 0 nor 6, and
!>                  IRSIGN neither 0 nor 1
!>           -3  => if MODE neither -6, 0 nor 6 and COND less than 1
!>           -4  => if MODE equals 6 or -6 and IDIST not in range 1 to 3
!>           -7  => if N negative
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
   SUBROUTINE SLATM7( MODE, COND, IRSIGN, IDIST, ISEED, D, N, &
                      RANK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   REAL               COND
   INTEGER            IDIST, INFO, IRSIGN, MODE, N, RANK
!     ..
!     .. Array Arguments ..
   REAL               D( * )
   INTEGER            ISEED( 4 )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ONE
   PARAMETER          ( ONE = 1.0E0 )
   REAL               ZERO
   PARAMETER          ( ZERO = 0.0E0 )
   REAL               HALF
   PARAMETER          ( HALF = 0.5E0 )
!     ..
!     .. Local Scalars ..
   REAL               ALPHA, TEMP
   INTEGER            I
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   REAL               SLARAN
   EXTERNAL           SLARAN
!     ..
!     .. External Subroutines ..
   EXTERNAL           SLARNV, XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, EXP, LOG, REAL
!     ..
!     .. Executable Statements ..
!
!     Decode and Test the input parameters. Initialize flags & seed.
!
   INFO = 0
!
!     Quick return if possible
!
   IF( N == 0 ) &
      RETURN
!
!     Set INFO if an error
!
   IF( MODE < -6 .OR. MODE > 6 ) THEN
      INFO = -1
   ELSE IF( ( MODE /= -6 .AND. MODE /= 0 .AND. MODE /= 6 ) .AND. &
            ( IRSIGN /= 0 .AND. IRSIGN /= 1 ) ) THEN
      INFO = -2
   ELSE IF( ( MODE /= -6 .AND. MODE /= 0 .AND. MODE /= 6 ) .AND. &
            COND < ONE ) THEN
      INFO = -3
   ELSE IF( ( MODE == 6 .OR. MODE == -6 ) .AND. &
            ( IDIST < 1 .OR. IDIST > 3 ) ) THEN
      INFO = -4
   ELSE IF( N < 0 ) THEN
      INFO = -7
   END IF
!
   IF( INFO /= 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'SLATM7', -INFO )
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
!     Compute D according to COND and MODE
!
   IF( MODE /= 0 ) THEN
      GO TO ( 100, 130, 160, 190, 210, 230 )ABS( MODE )
!
!        One large D value:
!
  100    CONTINUE
      DO I = 2, RANK
         D( I ) = ONE / COND
         ENDDO
      DO I = RANK + 1, N
         D( I ) = ZERO
         ENDDO
      D( 1 ) = ONE
      GO TO 240
!
!        One small D value:
!
  130    CONTINUE
      DO I = 1, RANK - 1
         D( I ) = ONE
         ENDDO
      DO I = RANK + 1, N
         D( I ) = ZERO
         ENDDO
      D( RANK ) = ONE / COND
      GO TO 240
!
!        Exponentially distributed D values:
!
  160    CONTINUE
      D( 1 ) = ONE
      IF( N > 1  .AND. RANK > 1 ) THEN
         ALPHA = COND**( -ONE / REAL( RANK-1 ) )
         DO I = 2, RANK
            D( I ) = ALPHA**( I-1 )
            ENDDO
         DO I = RANK + 1, N
            D( I ) = ZERO
            ENDDO
      END IF
      GO TO 240
!
!        Arithmetically distributed D values:
!
  190    CONTINUE
      D( 1 ) = ONE
      IF( N > 1 ) THEN
         TEMP = ONE / COND
         ALPHA = ( ONE-TEMP ) / REAL( N-1 )
         DO I = 2, N
            D( I ) = REAL( N-I )*ALPHA + TEMP
            ENDDO
      END IF
      GO TO 240
!
!        Randomly distributed D values on ( 1/COND , 1):
!
  210    CONTINUE
      ALPHA = LOG( ONE / COND )
      DO I = 1, N
         D( I ) = EXP( ALPHA*SLARAN( ISEED ) )
         ENDDO
      GO TO 240
!
!        Randomly distributed D values from IDIST
!
  230    CONTINUE
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL SLARNV( IDIST, ISEED, N, D )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : SLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
  240    CONTINUE
!
!        If MODE neither -6 nor 0 nor 6, and IRSIGN = 1, assign
!        random signs to D
!
      IF( ( MODE /= -6 .AND. MODE /= 0 .AND. MODE /= 6 ) .AND. &
          IRSIGN == 1 ) THEN
         DO I = 1, N
            TEMP = SLARAN( ISEED )
            IF( TEMP > HALF ) &
               D( I ) = -D( I )
            ENDDO
      END IF
!
!        Reverse if MODE < 0
!
      IF( MODE < 0 ) THEN
         DO I = 1, N / 2
            TEMP = D( I )
            D( I ) = D( N+1-I )
            D( N+1-I ) = TEMP
            ENDDO
      END IF
!
   END IF
!
   RETURN
!
!     End of SLATM7
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        


