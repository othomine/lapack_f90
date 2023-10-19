!> \brief \b DLATM1
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D, N, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IDIST, INFO, IRSIGN, MODE, N
!       DOUBLE PRECISION   COND
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       DOUBLE PRECISION   D( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DLATM1 computes the entries of D(1..N) as specified by
!>    MODE, COND and IRSIGN. IDIST and ISEED determine the generation
!>    of random numbers. DLATM1 is called by DLATMR to generate
!>    random test matrices for LAPACK programs.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] MODE
!> \verbatim
!>          MODE is INTEGER
!>           On entry describes how D is to be computed:
!>           MODE = 0 means do not change D.
!>           MODE = 1 sets D(1)=1 and D(2:N)=1.0/COND
!>           MODE = 2 sets D(1:N-1)=1 and D(N)=1.0/COND
!>           MODE = 3 sets D(I)=COND**(-(I-1)/(N-1))
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
!> \endverbatim
!>
!> \param[in] COND
!> \verbatim
!>          COND is DOUBLE PRECISION
!>           On entry, used as described under MODE above.
!>           If used, it must be >= 1. Not modified.
!> \endverbatim
!>
!> \param[in] IRSIGN
!> \verbatim
!>          IRSIGN is INTEGER
!>           On entry, if MODE neither -6, 0 nor 6, determines sign of
!>           entries of D
!>           0 => leave entries of D unchanged
!>           1 => multiply each entry of D by 1 or -1 with probability .5
!> \endverbatim
!>
!> \param[in] IDIST
!> \verbatim
!>          IDIST is INTEGER
!>           On entry, IDIST specifies the type of distribution to be
!>           used to generate a random matrix .
!>           1 => UNIFORM( 0, 1 )
!>           2 => UNIFORM( -1, 1 )
!>           3 => NORMAL( 0, 1 )
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension ( 4 )
!>           On entry ISEED specifies the seed of the random number
!>           generator. The random number generator uses a
!>           linear congruential sequence limited to small
!>           integers, and so should produce machine independent
!>           random numbers. The values of ISEED are changed on
!>           exit, and can be used in the next call to DLATM1
!>           to continue the same random number sequence.
!>           Changed on exit.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension ( N )
!>           Array to be computed according to MODE, COND and IRSIGN.
!>           May be changed on exit if MODE is nonzero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           Number of entries of D. Not modified.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
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
!> \ingroup double_matgen
!
!  =====================================================================
   SUBROUTINE DLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D, N, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            IDIST, INFO, IRSIGN, MODE, N
   DOUBLE PRECISION   COND
!     ..
!     .. Array Arguments ..
   INTEGER            ISEED( 4 )
   DOUBLE PRECISION   D( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ONE
   PARAMETER          ( ONE = 1.0D0 )
   DOUBLE PRECISION   HALF
   PARAMETER          ( HALF = 0.5D0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I
   DOUBLE PRECISION   ALPHA, TEMP
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLARAN
   EXTERNAL           DLARAN
!     ..
!     .. External Subroutines ..
   EXTERNAL           DLARNV, XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, DBLE, EXP, LOG
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
      CALL XERBLA( 'DLATM1', -INFO )
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
      GO TO ( 10, 30, 50, 70, 90, 110 )ABS( MODE )
!
!        One large D value:
!
10    CONTINUE
      DO I = 1, N
         D( I ) = ONE / COND
      ENDDO
      D( 1 ) = ONE
      GO TO 120
!
!        One small D value:
!
30    CONTINUE
      DO I = 1, N
         D( I ) = ONE
      ENDDO
      D( N ) = ONE / COND
      GO TO 120
!
!        Exponentially distributed D values:
!
50    CONTINUE
      D( 1 ) = ONE
      IF( N > 1 ) THEN
         ALPHA = COND**( -ONE / DBLE( N-1 ) )
         DO I = 2, N
            D( I ) = ALPHA**( I-1 )
         ENDDO
      END IF
      GO TO 120
!
!        Arithmetically distributed D values:
!
70    CONTINUE
      D( 1 ) = ONE
      IF( N > 1 ) THEN
         TEMP = ONE / COND
         ALPHA = ( ONE-TEMP ) / DBLE( N-1 )
         DO I = 2, N
            D( I ) = DBLE( N-I )*ALPHA + TEMP
         ENDDO
      END IF
      GO TO 120
!
!        Randomly distributed D values on ( 1/COND , 1):
!
90    CONTINUE
      ALPHA = LOG( ONE / COND )
      DO I = 1, N
         D( I ) = EXP( ALPHA*DLARAN( ISEED ) )
         ENDDO
      GO TO 120
!
!        Randomly distributed D values from IDIST
!
  110    CONTINUE
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLARNV( IDIST, ISEED, N, D )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
  120    CONTINUE
!
!        If MODE neither -6 nor 0 nor 6, and IRSIGN = 1, assign
!        random signs to D
!
      IF( ( MODE /= -6 .AND. MODE /= 0 .AND. MODE /= 6 ) .AND. &
          IRSIGN == 1 ) THEN
         DO I = 1, N
            TEMP = DLARAN( ISEED )
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
!     End of DLATM1
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        


