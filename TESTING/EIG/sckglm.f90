!> \brief \b SCKGLM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SCKGLM( NN, MVAL, PVAL, NVAL, NMATS, ISEED, THRESH,
!                          NMAX, A, AF, B, BF, X, WORK, RWORK, NIN, NOUT,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, NIN, NMATS, NMAX, NN, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 ), MVAL( * ), NVAL( * ), PVAL( * )
!       REAL               A( * ), AF( * ), B( * ), BF( * ), RWORK( * ),
!      $                   WORK( * ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SCKGLM tests SGGGLM - subroutine for solving generalized linear
!>                       model problem.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER
!>          The number of values of N, M and P contained in the vectors
!>          NVAL, MVAL and PVAL.
!> \endverbatim
!>
!> \param[in] MVAL
!> \verbatim
!>          MVAL is INTEGER array, dimension (NN)
!>          The values of the matrix column dimension M.
!> \endverbatim
!>
!> \param[in] PVAL
!> \verbatim
!>          PVAL is INTEGER array, dimension (NN)
!>          The values of the matrix column dimension P.
!> \endverbatim
!>
!> \param[in] NVAL
!> \verbatim
!>          NVAL is INTEGER array, dimension (NN)
!>          The values of the matrix row dimension N.
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
!>          THRESH is REAL
!>          The threshold value for the test ratios.  A result is
!>          included in the output file if RESID >= THRESH.  To have
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
!>          A is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] BF
!> \verbatim
!>          BF is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is REAL array, dimension (4*NMAX)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (NMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (NMAX*NMAX)
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
!>          > 0 :  If SLATMS returns an error code, the absolute value
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
!> \ingroup single_eig
!
!  =====================================================================
   SUBROUTINE SCKGLM( NN, MVAL, PVAL, NVAL, NMATS, ISEED, THRESH, &
                      NMAX, A, AF, B, BF, X, WORK, RWORK, NIN, NOUT, &
                      INFO )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, NIN, NMATS, NMAX, NN, NOUT
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   INTEGER            ISEED( 4 ), MVAL( * ), NVAL( * ), PVAL( * )
   REAL               A( * ), AF( * ), B( * ), BF( * ), RWORK( * ), &
                      WORK( * ), X( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            NTYPES
   PARAMETER          ( NTYPES = 8 )
!     ..
!     .. Local Scalars ..
   LOGICAL            FIRSTT
   CHARACTER          DISTA, DISTB, TYPE
   CHARACTER*3        PATH
   INTEGER            I, IINFO, IK, IMAT, KLA, KLB, KUA, KUB, LDA, &
                      LDB, LWORK, M, MODEA, MODEB, N, NFAIL, NRUN, P
   REAL               ANORM, BNORM, CNDNMA, CNDNMB, RESID
!     ..
!     .. Local Arrays ..
   LOGICAL            DOTYPE( NTYPES )
!     ..
!     .. External Functions ..
   REAL               SLARND
   EXTERNAL           SLARND
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAHDG, ALAREQ, ALASUM, SGLMTS, SLATB9, SLATMS
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
!     Initialize constants.
!
   PATH( 1: 3 ) = 'GLM'
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
      IF( M > N .OR. N > M+P ) THEN
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
      IF( M <= N .and. N <= M+P ) THEN
!
      DO IMAT = 1, NTYPES
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
         IF (DOTYPE( IMAT ) ) THEN
!
!           Set up parameters with SLATB9 and generate test
!           matrices A and B with SLATMS.
!
         CALL SLATB9( PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, &
                      ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB, &
                      DISTA, DISTB )
!
         CALL SLATMS( N, M, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, &
                      ANORM, KLA, KUA, 'No packing', A, LDA, WORK, &
                      IINFO )
         IF( IINFO /= 0 ) THEN
            WRITE( NOUT, FMT = 9999 )IINFO
            INFO = ABS( IINFO )
            CYCLE
         END IF
!
         CALL SLATMS( N, P, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, &
                      BNORM, KLB, KUB, 'No packing', B, LDB, WORK, &
                      IINFO )
         IF( IINFO /= 0 ) THEN
            WRITE( NOUT, FMT = 9999 )IINFO
            INFO = ABS( IINFO )
            CYCLE
         END IF
!
!           Generate random left hand side vector of GLM
!
         DO I = 1, N
            X( I ) = SLARND( 2, ISEED )
         ENDDO
!
         CALL SGLMTS( N, M, P, A, AF, LDA, B, BF, LDB, X, &
                      X( NMAX+1 ), X( 2*NMAX+1 ), X( 3*NMAX+1 ), &
                      WORK, LWORK, RWORK, RESID )
!
!           Print information about the tests that did not
!           pass the threshold.
!
         IF( RESID >= THRESH ) THEN
            IF( NFAIL == 0 .AND. FIRSTT ) THEN
               FIRSTT = .FALSE.
               CALL ALAHDG( NOUT, PATH )
            END IF
            WRITE( NOUT, FMT = 9998 )N, M, P, IMAT, 1, RESID
            NFAIL = NFAIL + 1
         END IF
         NRUN = NRUN + 1
!
      ENDIF
      ENDDO
   ENDIF
   ENDDO
!
!     Print a summary of the results.
!
   CALL ALASUM( PATH, NOUT, NFAIL, NRUN, 0 )
!
 9999 FORMAT( ' SLATMS in SCKGLM INFO = ', I5 )
 9998 FORMAT( ' N=', I4, ' M=', I4, ', P=', I4, ', type ', I2, &
         ', test ', I2, ', ratio=', G13.6 )
 9997 FORMAT( ' *** Invalid input  for GLM:  M = ', I6, ', P = ', I6, &
         ', N = ', I6, ';', / '     must satisfy M <= N <= M+P  ', &
         '(this set of values will be skipped)' )
   RETURN
!
!     End of SCKGLM
!
END


