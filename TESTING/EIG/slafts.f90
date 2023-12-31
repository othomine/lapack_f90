!> \brief \b SLAFTS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAFTS( TYPE, M, N, IMAT, NTESTS, RESULT, ISEED,
!                          THRESH, IOUNIT, IE )
!
!       .. Scalar Arguments ..
!       CHARACTER*3        TYPE
!       INTEGER            IE, IMAT, IOUNIT, M, N, NTESTS
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       REAL               RESULT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SLAFTS tests the result vector against the threshold value to
!>    see which tests for this matrix type failed to pass the threshold.
!>    Output is to the file given by unit IOUNIT.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \verbatim
!>  TYPE   - CHARACTER*3
!>           On entry, TYPE specifies the matrix type to be used in the
!>           printed messages.
!>           Not modified.
!>
!>  N      - INTEGER
!>           On entry, N specifies the order of the test matrix.
!>           Not modified.
!>
!>  IMAT   - INTEGER
!>           On entry, IMAT specifies the type of the test matrix.
!>           A listing of the different types is printed by SLAHD2
!>           to the output file if a test fails to pass the threshold.
!>           Not modified.
!>
!>  NTESTS - INTEGER
!>           On entry, NTESTS is the number of tests performed on the
!>           subroutines in the path given by TYPE.
!>           Not modified.
!>
!>  RESULT - REAL               array of dimension( NTESTS )
!>           On entry, RESULT contains the test ratios from the tests
!>           performed in the calling program.
!>           Not modified.
!>
!>  ISEED  - INTEGER            array of dimension( 4 )
!>           Contains the random seed that generated the matrix used
!>           for the tests whose ratios are in RESULT.
!>           Not modified.
!>
!>  THRESH - REAL
!>           On entry, THRESH specifies the acceptable threshold of the
!>           test ratios.  If RESULT( K ) > THRESH, then the K-th test
!>           did not pass the threshold and a message will be printed.
!>           Not modified.
!>
!>  IOUNIT - INTEGER
!>           On entry, IOUNIT specifies the unit number of the file
!>           to which the messages are printed.
!>           Not modified.
!>
!>  IE     - INTEGER
!>           On entry, IE contains the number of tests which have
!>           failed to pass the threshold so far.
!>           Updated on exit if any of the ratios in RESULT also fail.
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
!> \ingroup single_eig
!
!  =====================================================================
   SUBROUTINE SLAFTS( TYPE, M, N, IMAT, NTESTS, RESULT, ISEED, &
                      THRESH, IOUNIT, IE )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER*3        TYPE
   INTEGER            IE, IMAT, IOUNIT, M, N, NTESTS
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   INTEGER            ISEED( 4 )
   REAL               RESULT( * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
   INTEGER            K
!     ..
!     .. External Subroutines ..
   EXTERNAL           SLAHD2
!     ..
!     .. Executable Statements ..
!
   IF( M == N ) THEN
!
!     Output for square matrices:
!
      DO K = 1, NTESTS
         IF( RESULT( K ) >= THRESH ) THEN
!
!           If this is the first test to fail, call SLAHD2
!           to print a header to the data file.
!
            IF( IE == 0 ) &
               CALL SLAHD2( IOUNIT, TYPE )
            IE = IE + 1
            IF( RESULT( K ) < 10000.0 ) THEN
               WRITE( IOUNIT, FMT = 9999 )N, IMAT, ISEED, K, &
                  RESULT( K )
 9999             FORMAT( ' Matrix order=', I5, ', type=', I2, &
                     ', seed=', 4( I4, ',' ), ' result ', I3, ' is', &
                     0P, F8.2 )
            ELSE
               WRITE( IOUNIT, FMT = 9998 )N, IMAT, ISEED, K, &
                  RESULT( K )
 9998             FORMAT( ' Matrix order=', I5, ', type=', I2, &
                     ', seed=', 4( I4, ',' ), ' result ', I3, ' is', &
                     1P, E10.3 )
            END IF
         END IF
      ENDDO
   ELSE
!
!     Output for rectangular matrices
!
      DO K = 1, NTESTS
         IF( RESULT( K ) >= THRESH ) THEN
!
!              If this is the first test to fail, call SLAHD2
!              to print a header to the data file.
!
            IF( IE == 0 ) &
               CALL SLAHD2( IOUNIT, TYPE )
            IE = IE + 1
            IF( RESULT( K ) < 10000.0 ) THEN
               WRITE( IOUNIT, FMT = 9997 )M, N, IMAT, ISEED, K, &
                  RESULT( K )
 9997             FORMAT( 1X, I5, ' x', I5, ' matrix, type=', I2, ', s', &
                     'eed=', 3( I4, ',' ), I4, ': result ', I3, &
                     ' is', 0P, F8.2 )
            ELSE
               WRITE( IOUNIT, FMT = 9996 )M, N, IMAT, ISEED, K, &
                  RESULT( K )
 9996             FORMAT( 1X, I5, ' x', I5, ' matrix, type=', I2, ', s', &
                     'eed=', 3( I4, ',' ), I4, ': result ', I3, &
                     ' is', 1P, E10.3 )
            END IF
         END IF
      ENDDO
!
   END IF
   RETURN
!
!     End of SLAFTS
!
END



