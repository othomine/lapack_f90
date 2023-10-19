!> \brief \b ZCKLSE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZCKLSE( NN, MVAL, PVAL, NVAL, NMATS, ISEED, THRESH,
!                          NMAX, A, AF, B, BF, X, WORK, RWORK, NIN, NOUT,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, NIN, NMATS, NMAX, NN, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 ), MVAL( * ), NVAL( * ), PVAL( * )
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( * ), AF( * ), B( * ), BF( * ), WORK( * ),
!      $                   X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZCKLSE tests ZGGLSE - a subroutine for solving linear equality
!> constrained least square problem (LSE).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER
!>          The number of values of (M,P,N) contained in the vectors
!>          (MVAL, PVAL, NVAL).
!> \endverbatim
!>
!> \param[in] MVAL
!> \verbatim
!>          MVAL is INTEGER array, dimension (NN)
!>          The values of the matrix row(column) dimension M.
!> \endverbatim
!>
!> \param[in] PVAL
!> \verbatim
!>          PVAL is INTEGER array, dimension (NN)
!>          The values of the matrix row(column) dimension P.
!> \endverbatim
!>
!> \param[in] NVAL
!> \verbatim
!>          NVAL is INTEGER array, dimension (NN)
!>          The values of the matrix column(row) dimension N.
!> \endverbatim
!>
!> \param[in] NMATS
!> \verbatim
!>          NMATS is INTEGER
!>          The number of matrix types to be tested for each combination
!>          of matrix dimensions.  If NMATS >= NTYPES (the maximum
!>          number of matrix types), then all the different types are
!>          generated for testing.  If NMATS < NTYPES, another input line
!>          is read to get the numbers of the matrix types to be used.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry, the seed of the random number generator.  The array
!>          elements should be between 0 and 4095, otherwise they will be
!>          reduced mod 4096, and ISEED(4) must be odd.
!>          On exit, the next seed in the random number sequence after
!>          all the test matrices have been generated.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is DOUBLE PRECISION
!>          The threshold value for the test ratios.  A result is
!>          included in the output file if RESULT >= THRESH.  To have
!>          every test ratio printed, use THRESH = 0.
!> \endverbatim
!>
!> \param[in] NMAX
!> \verbatim
!>          NMAX is INTEGER
!>          The maximum value permitted for M or N, used in dimensioning
!>          the work arrays.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] BF
!> \verbatim
!>          BF is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (5*NMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (NMAX)
!> \endverbatim
!>
!> \param[in] NIN
!> \verbatim
!>          NIN is INTEGER
!>          The unit number for input.
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The unit number for output.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0 :  successful exit
!>          > 0 :  If ZLATMS returns an error code, the absolute value
!>                 of it is returned.
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
!> \ingroup complex16_eig
!
!  =====================================================================
   SUBROUTINE ZCKLSE( NN, MVAL, PVAL, NVAL, NMATS, ISEED, THRESH, &
                      NMAX, A, AF, B, BF, X, WORK, RWORK, NIN, NOUT, &
                      INFO )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, NIN, NMATS, NMAX, NN, NOUT
   DOUBLE PRECISION   THRESH
!     ..
!     .. Array Arguments ..
   INTEGER            ISEED( 4 ), MVAL( * ), NVAL( * ), PVAL( * )
   DOUBLE PRECISION   RWORK( * )
   COMPLEX*16         A( * ), AF( * ), B( * ), BF( * ), WORK( * ), &
                      X( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 7 )
   INTEGER            NTYPES
   PARAMETER          ( NTYPES = 8 )
!     ..
!     .. Local Scalars ..
   LOGICAL            FIRSTT
   CHARACTER          DISTA, DISTB, TYPE
   CHARACTER*3        PATH
   INTEGER            I, IINFO, IK, IMAT, KLA, KLB, KUA, KUB, LDA, &
                      LDB, LWORK, M, MODEA, MODEB, N, NFAIL, NRUN, &
                      NT, P
   DOUBLE PRECISION   ANORM, BNORM, CNDNMA, CNDNMB
!     ..
!     .. Local Arrays ..
   LOGICAL            DOTYPE( NTYPES )
   DOUBLE PRECISION   RESULT( NTESTS )
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAHDG, ALAREQ, ALASUM, DLATB9, ZLARHS, ZLATMS, &
                      ZLSETS
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   PATH( 1: 3 ) = 'LSE'
   INFO = 0
   NRUN = 0
   NFAIL = 0
   FIRSTT = .TRUE.
   CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
   LDA = NMAX
   LDB = NMAX
   LWORK = NMAX*NMAX
!
!     Check for valid input values.
!
   DO IK = 1, NN
      M = MVAL( IK )
      P = PVAL( IK )
      N = NVAL( IK )
      IF( P > N .OR. N > M+P ) THEN
         IF( FIRSTT ) THEN
            WRITE( NOUT, FMT = * )
            FIRSTT = .FALSE.
         END IF
         WRITE( NOUT, FMT = 9997 )M, P, N
      END IF
   ENDDO
   FIRSTT = .TRUE.
!
!     Do for each value of M in MVAL.
!
   DO IK = 1, NN
      M = MVAL( IK )
      P = PVAL( IK )
      N = NVAL( IK )
      IF( P <= N .AND. N <= M+P ) THEN
!
      DO IMAT = 1, NTYPES
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
         IF( .NOT.DOTYPE( IMAT ) ) GO TO 30
!
!           Set up parameters with DLATB9 and generate test
!           matrices A and B with ZLATMS.
!
         CALL DLATB9( PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, &
                      ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB, &
                      DISTA, DISTB )
!
         CALL ZLATMS( M, N, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, &
                      ANORM, KLA, KUA, 'No packing', A, LDA, WORK, &
                      IINFO )
         IF( IINFO /= 0 ) THEN
            WRITE( NOUT, FMT = 9999 )IINFO
            INFO = ABS( IINFO )
            GO TO 30
         END IF
!
         CALL ZLATMS( P, N, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, &
                      BNORM, KLB, KUB, 'No packing', B, LDB, WORK, &
                      IINFO )
         IF( IINFO /= 0 ) THEN
            WRITE( NOUT, FMT = 9999 )IINFO
            INFO = ABS( IINFO )
            GO TO 30
         END IF
!
!           Generate the right-hand sides C and D for the LSE.
!
         CALL ZLARHS( 'ZGE', 'New solution', 'Upper', 'N', M, N, &
                      MAX( M-1, 0 ), MAX( N-1, 0 ), 1, A, LDA, &
                      X( 4*NMAX+1 ), MAX( N, 1 ), X, MAX( M, 1 ), &
                      ISEED, IINFO )
!
         CALL ZLARHS( 'ZGE', 'Computed', 'Upper', 'N', P, N, &
                      MAX( P-1, 0 ), MAX( N-1, 0 ), 1, B, LDB, &
                      X( 4*NMAX+1 ), MAX( N, 1 ), X( 2*NMAX+1 ), &
                      MAX( P, 1 ), ISEED, IINFO )
!
         NT = 2
!
         CALL ZLSETS( M, P, N, A, AF, LDA, B, BF, LDB, X, &
                      X( NMAX+1 ), X( 2*NMAX+1 ), X( 3*NMAX+1 ), &
                      X( 4*NMAX+1 ), WORK, LWORK, RWORK, &
                      RESULT( 1 ) )
!
!           Print information about the tests that did not
!           pass the threshold.
!
         DO I = 1, NT
            IF( RESULT( I ) >= THRESH ) THEN
               IF( NFAIL == 0 .AND. FIRSTT ) THEN
                  FIRSTT = .FALSE.
                  CALL ALAHDG( NOUT, PATH )
               END IF
               WRITE( NOUT, FMT = 9998 )M, P, N, IMAT, I, &
                  RESULT( I )
               NFAIL = NFAIL + 1
            END IF
         ENDDO
         NRUN = NRUN + NT
!
30    CONTINUE
      ENDDO
      ENDIF
   ENDDO
!
!     Print a summary of the results.
!
   CALL ALASUM( PATH, NOUT, NFAIL, NRUN, 0 )
!
 9999 FORMAT( ' ZLATMS in ZCKLSE   INFO = ', I5 )
 9998 FORMAT( ' M=', I4, ' P=', I4, ', N=', I4, ', type ', I2, &
         ', test ', I2, ', ratio=', G13.6 )
 9997 FORMAT( ' *** Invalid input  for LSE:  M = ', I6, ', P = ', I6, &
         ', N = ', I6, ';', / '     must satisfy P <= N <= P+M  ', &
         '(this set of values will be skipped)' )
   RETURN
!
!     End of ZCKLSE
!
END


