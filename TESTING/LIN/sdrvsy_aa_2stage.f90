!> \brief \b SDRVSY_AA_2STAGE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SDRVSY_AA_2STAGE(
!                             DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
!                             A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK,
!                             NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NOUT, NRHS
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NVAL( * )
!       REAL               RWORK( * )
!       REAL               A( * ), AFAC( * ), AINV( * ), B( * ),
!      $                   WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SDRVSY_AA_2STAGE tests the driver routine SSYSV_AA_2STAGE.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DOTYPE
!> \verbatim
!>          DOTYPE is LOGICAL array, dimension (NTYPES)
!>          The matrix types to be used for testing.  Matrices of type j
!>          (for 1 <= j <= NTYPES) are used for testing if DOTYPE(j) =
!>          .TRUE.; if DOTYPE(j) = .FALSE., then type j is not used.
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
!>          The values of the matrix dimension N.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand side vectors to be generated for
!>          each linear system.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is REAL
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
!> \param[in] NMAX
!> \verbatim
!>          NMAX is INTEGER
!>          The maximum value permitted for N, used in dimensioning the
!>          work arrays.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is REAL array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is REAL array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is REAL array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (NMAX*max(2,NRHS))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (NMAX+2*NRHS)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (NMAX)
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
!
!> \ingroup real_lin
!
!  =====================================================================
   SUBROUTINE SDRVSY_AA_2STAGE( &
                            DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                            NMAX, A, AFAC, AINV, B, X, XACT, WORK, &
                            RWORK, IWORK, NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            TSTERR
   INTEGER            NMAX, NN, NOUT, NRHS
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            IWORK( * ), NVAL( * )
   REAL               RWORK( * )
   REAL               A( * ), AFAC( * ), AINV( * ), B( * ), &
                      WORK( * ), X( * ), XACT( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ONE, ZERO
   PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
   INTEGER            NTYPES, NTESTS
   PARAMETER          ( NTYPES = 10, NTESTS = 3 )
   INTEGER            NFACT
   PARAMETER          ( NFACT = 2 )
!     ..
!     .. Local Scalars ..
   LOGICAL            ZEROT
   CHARACTER          DIST, FACT, TYPE, UPLO, XTYPE
   CHARACTER*3        MATPATH, PATH
   INTEGER            I, I1, I2, IFACT, IMAT, IN, INFO, IOFF, IUPLO, &
                      IZERO, J, K, KL, KU, LDA, LWORK, MODE, N, &
                      NB, NBMIN, NERRS, NFAIL, NIMAT, NRUN, NT
   REAL               ANORM, CNDNUM
!     ..
!     .. Local Arrays ..
   CHARACTER          FACTS( NFACT ), UPLOS( 2 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   REAL               RESULT( NTESTS )
!     ..
!     .. External Functions ..
   REAL               SLANSY, SGET06
   EXTERNAL           SLANSY, SGET06
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALADHD, ALAERH, ALASVM, XLAENV, SERRVX, &
                      SLACPY, SLARHS, SLATB4, SLATMS, &
                      SSYSV_AA_2STAGE, SSYT01_AA, SPOT02, &
                      SSYTRF_AA_2STAGE
!     ..
!     .. Scalars in Common ..
   LOGICAL            LERR, OK
   CHARACTER*32       SRNAMT
   INTEGER            INFOT, NUNIT
!     ..
!     .. Common blocks ..
   COMMON             / INFOC / INFOT, NUNIT, OK, LERR
   COMMON             / SRNAMC / SRNAMT
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          CMPLX, MAX, MIN
!     ..
!     .. Data statements ..
   DATA               ISEEDY / 1988, 1989, 1990, 1991 /
   DATA               UPLOS / 'U', 'L' / , FACTS / 'F', 'N' /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
!     Test path
!
   PATH( 1: 1 ) = 'Single precision'
   PATH( 2: 3 ) = 'S2'
!
!     Path to generate matrices
!
   MATPATH( 1: 1 ) = 'Single precision'
   MATPATH( 2: 3 ) = 'SY'
!
   NRUN = 0
   NFAIL = 0
   NERRS = 0
   DO I = 1, 4
      ISEED( I ) = ISEEDY( I )
   ENDDO
!
!     Test the error exits
!
   IF( TSTERR ) &
      CALL SERRVX( PATH, NOUT )
   INFOT = 0
!
!     Set the block size and minimum block size for testing.
!
   NB = 1
   NBMIN = 2
   CALL XLAENV( 1, NB )
   CALL XLAENV( 2, NBMIN )
!
!     Do for each value of N in NVAL
!
   DO IN = 1, NN
      N = NVAL( IN )
      LDA = MAX( N, 1 )
      XTYPE = 'N'
      NIMAT = NTYPES
      IF( N <= 0 ) &
         NIMAT = 1
!
      DO IMAT = 1, NIMAT
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
         IF( .NOT.DOTYPE( IMAT ) ) &
            GO TO 170
!
!           Skip types 3, 4, 5, or 6 if the matrix size is too small.
!
         ZEROT = IMAT >= 3 .AND. IMAT <= 6
         IF( ZEROT .AND. N < IMAT-2 ) &
            GO TO 170
!
!           Do first for UPLO = 'U', then for UPLO = 'L'
!
         DO IUPLO = 1, 2
            UPLO = UPLOS( IUPLO )
!
!              Begin generate the test matrix A.
!
!              Set up parameters with SLATB4 for the matrix generator
!              based on the type of matrix to be generated.
!
           CALL SLATB4( MATPATH, IMAT, N, N, TYPE, KL, KU, ANORM, &
                            MODE, CNDNUM, DIST )
!
!              Generate a matrix with SLATMS.
!
               SRNAMT = 'SLATMS'
               CALL SLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, &
                            CNDNUM, ANORM, KL, KU, UPLO, A, LDA, &
                            WORK, INFO )
!
!                 Check error code from SLATMS and handle error.
!
               IF( INFO /= 0 ) THEN
                  CALL ALAERH( PATH, 'SLATMS', INFO, 0, UPLO, N, N, &
                               -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 160
               END IF
!
!                 For types 3-6, zero one or more rows and columns of
!                 the matrix to test that INFO is returned correctly.
!
               IF( ZEROT ) THEN
                  IF( IMAT == 3 ) THEN
                     IZERO = 1
                  ELSE IF( IMAT == 4 ) THEN
                     IZERO = N
                  ELSE
                     IZERO = N / 2 + 1
                  END IF
!
                  IF( IMAT < 6 ) THEN
!
!                       Set row and column IZERO to zero.
!
                     IF( IUPLO == 1 ) THEN
                        IOFF = ( IZERO-1 )*LDA
                        DO I = 1, IZERO - 1
                           A( IOFF+I ) = ZERO
                        ENDDO
                        IOFF = IOFF + IZERO
                        DO I = IZERO, N
                           A( IOFF ) = ZERO
                           IOFF = IOFF + LDA
                        ENDDO
                     ELSE
                        IOFF = IZERO
                        DO I = 1, IZERO - 1
                           A( IOFF ) = ZERO
                           IOFF = IOFF + LDA
                        ENDDO
                        IOFF = IOFF - IZERO
                        DO I = IZERO, N
                           A( IOFF+I ) = ZERO
                        ENDDO
                     END IF
                  ELSE
                     IOFF = 0
                     IF( IUPLO == 1 ) THEN
!
!                       Set the first IZERO rows and columns to zero.
!
                        DO J = 1, N
                           I2 = MIN( J, IZERO )
                           DO I = 1, I2
                              A( IOFF+I ) = ZERO
                           ENDDO
                           IOFF = IOFF + LDA
                        ENDDO
                        IZERO = 1
                     ELSE
!
!                       Set the first IZERO rows and columns to zero.
!
                        IOFF = 0
                        DO J = 1, N
                           I1 = MAX( J, IZERO )
                           DO I = I1, N
                              A( IOFF+I ) = ZERO
                           ENDDO
                           IOFF = IOFF + LDA
                        ENDDO
                     END IF
                  END IF
               ELSE
                  IZERO = 0
               END IF
!
!                 End generate the test matrix A.
!
!
            DO IFACT = 1, NFACT
!
!                 Do first for FACT = 'F', then for other values.
!
               FACT = FACTS( IFACT )
!
!                 Form an exact solution and set the right hand side.
!
               SRNAMT = 'SLARHS'
               CALL SLARHS( MATPATH, XTYPE, UPLO, ' ', N, N, KL, KU, &
                            NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, &
                            INFO )
               XTYPE = 'C'
!
!                 --- Test SSYSV_AA_2STAGE  ---
!
               IF( IFACT == 2 ) THEN
                  CALL SLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
                  CALL SLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
!
!                    Factor the matrix and solve the system using SSYSV_AA.
!
                  SRNAMT = 'SSYSV_AA_2STAGE '
                  LWORK = MIN(N*NB, 3*NMAX*NMAX)
                  CALL SSYSV_AA_2STAGE( UPLO, N, NRHS, AFAC, LDA, &
                                    AINV, (3*NB+1)*N, &
                                    IWORK, IWORK( 1+N ), &
                                    X, LDA, WORK, LWORK, INFO )
!
!                    Adjust the expected value of INFO to account for
!                    pivoting.
!
                  IF( IZERO > 0 ) THEN
                     J = 1
                     K = IZERO
  100                   CONTINUE
                     IF( J == K ) THEN
                        K = IWORK( J )
                     ELSE IF( IWORK( J ) == K ) THEN
                        K = J
                     END IF
                     IF( J < K ) THEN
                        J = J + 1
                        GO TO 100
                     END IF
                  ELSE
                     K = 0
                  END IF
!
!                    Check error code from SSYSV_AA .
!
                  IF( INFO /= K ) THEN
                     CALL ALAERH( PATH, 'SSYSV_AA', INFO, K, &
                                  UPLO, N, N, -1, -1, NRHS, &
                                  IMAT, NFAIL, NERRS, NOUT )
                     GO TO 120
                  ELSE IF( INFO /= 0 ) THEN
                     GO TO 120
                  END IF
!
!                    Compute residual of the computed solution.
!
                  CALL SLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
                  CALL SPOT02( UPLO, N, NRHS, A, LDA, X, LDA, WORK, &
                               LDA, RWORK, RESULT( 1 ) )
!
!                    Reconstruct matrix from factors and compute
!                    residual.
!
!                     CALL SSY01_AA( UPLO, N, A, LDA, AFAC, LDA,
!     $                                  IWORK, AINV, LDA, RWORK,
!     $                                  RESULT( 2 ) )
!                     NT = 2
                  NT = 1
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                  DO K = 1, NT
                     IF( RESULT( K ) >= THRESH ) THEN
                        IF( NFAIL == 0 .AND. NERRS == 0 ) &
                           CALL ALADHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9999 )'SSYSV_AA ', &
                            UPLO, N, IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
                     ENDDO
                  NRUN = NRUN + NT
  120                CONTINUE
               END IF
!
               ENDDO
!
  160       CONTINUE
            ENDDO
  170    CONTINUE
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I2, &
         ', test ', I2, ', ratio =', G12.5 )
   RETURN
!
!     End of SDRVSY_AA_2STAGE
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
