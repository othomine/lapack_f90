!> \brief \b ZCHKUNHR_COL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZCHKUNHR_COL( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB,
!                                NBVAL, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NM, NN, NNB, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            MVAL( * ), NBVAL( * ), NVAL( * )
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZCHKUNHR_COL tests:
!>   1) ZUNGTSQR and ZUNHR_COL using ZLATSQR, ZGEMQRT,
!>   2) ZUNGTSQR_ROW and ZUNHR_COL inside ZGETSQRHRT
!>      (which calls ZLATSQR, ZUNGTSQR_ROW and ZUNHR_COL) using ZGEMQRT.
!> Therefore, ZLATSQR (part of ZGEQR), ZGEMQRT (part of ZGEMQR)
!> have to be tested before this test.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] THRESH
!> \verbatim
!>          THRESH is DOUBLE PRECISION
!>          The threshold value for the test ratios.  A result is
!>          included in the output file if RESULT >= THRESH.  To have
!>          every test ratio printed, use THRESH = 0.
!> \endverbatim
!>
!> \param[in] TSTERR
!> \verbatim
!>          TSTERR is LOGICAL
!>          Flag that indicates whether error exits are to be tested.
!> \endverbatim
!>
!> \param[in] NM
!> \verbatim
!>          NM is INTEGER
!>          The number of values of M contained in the vector MVAL.
!> \endverbatim
!>
!> \param[in] MVAL
!> \verbatim
!>          MVAL is INTEGER array, dimension (NM)
!>          The values of the matrix row dimension M.
!> \endverbatim
!>
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER
!>          The number of values of N contained in the vector NVAL.
!> \endverbatim
!>
!> \param[in] NVAL
!> \verbatim
!>          NVAL is INTEGER array, dimension (NN)
!>          The values of the matrix column dimension N.
!> \endverbatim
!>
!> \param[in] NNB
!> \verbatim
!>          NNB is INTEGER
!>          The number of values of NB contained in the vector NBVAL.
!> \endverbatim
!>
!> \param[in] NBVAL
!> \verbatim
!>          NBVAL is INTEGER array, dimension (NNB)
!>          The values of the blocksize NB.
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The unit number for output.
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
!> \ingroup complex16_lin
!
!  =====================================================================
   SUBROUTINE ZCHKUNHR_COL( THRESH, TSTERR, NM, MVAL, NN, NVAL, &
                            NNB, NBVAL, NOUT )
   IMPLICIT NONE
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            TSTERR
   INTEGER            NM, NN, NNB, NOUT
   DOUBLE PRECISION   THRESH
!     ..
!     .. Array Arguments ..
   INTEGER            MVAL( * ), NBVAL( * ), NVAL( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 6 )
!     ..
!     .. Local Scalars ..
   CHARACTER(LEN=3)   PATH
   INTEGER            I, IMB1, INB1, INB2, J, T, M, N, MB1, NB1, &
                      NB2, NFAIL, NERRS, NRUN
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!
!     .. Local Arrays ..
   DOUBLE PRECISION   RESULT( NTESTS )
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAHD, ALASUM, ZERRUNHR_COL, ZUNHR_COL01, &
                      ZUNHR_COL02
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX, MIN
!     ..
!     .. Scalars in Common ..
   LOGICAL            LERR, OK
   CHARACTER(LEN=32)  SRNAMT
   INTEGER            INFOT, NUNIT
!     ..
!     .. Common blocks ..
   COMMON             / INFOC / INFOT, NUNIT, OK, LERR
   COMMON             / SRNAMC / SRNAMT
!     ..
!     .. Executable Statements ..
!
!     Initialize constants
!
   PATH( 1: 1 ) = 'Z'
   PATH( 2: 3 ) = 'HH'
   NRUN = 0
   NFAIL = 0
   NERRS = 0
!
!     Test the error exits
!
   IF( TSTERR ) CALL ZERRUNHR_COL( PATH, NOUT )
   INFOT = 0
!
!     Do for each value of M in MVAL.
!
   DO I = 1, NM
      M = MVAL( I )
!
!        Do for each value of N in NVAL.
!
      DO J = 1, NN
         N = NVAL( J )
!
!           Only for M >= N
!
         IF ( MIN( M, N ) > 0 .AND. M >= N ) THEN
!
!              Do for each possible value of MB1
!
            DO IMB1 = 1, NNB
               MB1 = NBVAL( IMB1 )
!
!                 Only for MB1 > N
!
               IF ( MB1 > N ) THEN
!
!                    Do for each possible value of NB1
!
                  DO INB1 = 1, NNB
                     NB1 = NBVAL( INB1 )
!
!                       Do for each possible value of NB2
!
                     DO INB2 = 1, NNB
                        NB2 = NBVAL( INB2 )
!
                        IF( NB1 > 0 .AND. NB2 > 0 ) THEN
!
!                             Test ZUNHR_COL
!
                           CALL ZUNHR_COL01( M, N, MB1, NB1, &
                                             NB2, RESULT )
!
!                             Print information about the tests that did
!                             not pass the threshold.
!
                           DO T = 1, NTESTS
                              IF( RESULT( T ) >= THRESH ) THEN
                                 IF( NFAIL == 0 .AND. NERRS == 0 ) &
                                 CALL ALAHD( NOUT, PATH )
                                 WRITE( NOUT, FMT = 9999 ) M, N, MB1, &
                                        NB1, NB2, T, RESULT( T )
                                 NFAIL = NFAIL + 1
                              END IF
                           END DO
                           NRUN = NRUN + NTESTS
                        END IF
                     END DO
                  END DO
               END IF
             END DO
         END IF
      END DO
   END DO
!
!     Do for each value of M in MVAL.
!
   DO I = 1, NM
      M = MVAL( I )
!
!        Do for each value of N in NVAL.
!
      DO J = 1, NN
         N = NVAL( J )
!
!           Only for M >= N
!
         IF ( MIN( M, N ) > 0 .AND. M >= N ) THEN
!
!              Do for each possible value of MB1
!
            DO IMB1 = 1, NNB
               MB1 = NBVAL( IMB1 )
!
!                 Only for MB1 > N
!
               IF ( MB1 > N ) THEN
!
!                    Do for each possible value of NB1
!
                  DO INB1 = 1, NNB
                     NB1 = NBVAL( INB1 )
!
!                       Do for each possible value of NB2
!
                     DO INB2 = 1, NNB
                        NB2 = NBVAL( INB2 )
!
                        IF( NB1 > 0 .AND. NB2 > 0 ) THEN
!
!                             Test ZUNHR_COL
!
                           CALL ZUNHR_COL02( M, N, MB1, NB1, &
                                             NB2, RESULT )
!
!                             Print information about the tests that did
!                             not pass the threshold.
!
                           DO T = 1, NTESTS
                              IF( RESULT( T ) >= THRESH ) THEN
                                 IF( NFAIL == 0 .AND. NERRS == 0 ) &
                                 CALL ALAHD( NOUT, PATH )
                                 WRITE( NOUT, FMT = 9998 ) M, N, MB1, &
                                        NB1, NB2, T, RESULT( T )
                                 NFAIL = NFAIL + 1
                              END IF
                           END DO
                           NRUN = NRUN + NTESTS
                        END IF
                     END DO
                  END DO
               END IF
             END DO
         END IF
      END DO
   END DO
!
!     Print a summary of the results.
!
   CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( 'ZUNGTSQR and ZUNHR_COL: M=', I5, ', N=', I5, &
           ', MB1=', I5, ', NB1=', I5, ', NB2=', I5, &
           ' test(', I2, ')=', G12.5 )
 9998 FORMAT( 'ZUNGTSQR_ROW and ZUNHR_COL: M=', I5, ', N=', I5, &
           ', MB1=', I5, ', NB1=', I5, ', NB2=', I5, &
           ' test(', I2, ')=', G12.5 )
   RETURN
!
!     End of ZCHKUNHR_COL
!
END
                                                                                                                                                                                                                                                                                                            



