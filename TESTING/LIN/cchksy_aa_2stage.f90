!> \brief \b CCHKSY_AA_2STAGE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CCHKSY_AA_2STAGE( DOTYPE, NN, NVAL, NNB, NBVAL,
!                             NNS, NSVAL, THRESH, TSTERR, NMAX, A,
!                             AFAC, AINV, B, X, XACT, WORK, RWORK,
!                             IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NNB, NNS, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * )
!       REAL               RWORK( * )
!       COMPLEX            A( * ), AFAC( * ), AINV( * ), B( * ),
!      $                   WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CCHKSY_AA_2STAGE tests CSYTRF_AA_2STAGE, -TRS_AA_2STAGE.
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
!> \param[in] NNS
!> \verbatim
!>          NNS is INTEGER
!>          The number of values of NRHS contained in the vector NSVAL.
!> \endverbatim
!>
!> \param[in] NSVAL
!> \verbatim
!>          NSVAL is INTEGER array, dimension (NNS)
!>          The values of the number of right hand sides NRHS.
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
!>          A is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX array, dimension (NMAX*NSMAX)
!>          where NSMAX is the largest entry in NSVAL.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is COMPLEX array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (NMAX*max(3,NSMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is COMPLEX array, dimension (max(NMAX,2*NSMAX))
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (2*NMAX)
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
!> \ingroup complex_lin
!
!  =====================================================================
   SUBROUTINE CCHKSY_AA_2STAGE( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, &
                         NSVAL, THRESH, TSTERR, NMAX, A, AFAC, AINV, &
                         B, X, XACT, WORK, RWORK, IWORK, NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
   IMPLICIT NONE
!
!     .. Scalar Arguments ..
   LOGICAL            TSTERR
   INTEGER            NN, NNB, NNS, NMAX, NOUT
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * )
   COMPLEX            A( * ), AFAC( * ), AINV( * ), B( * ), &
                      WORK( * ), X( * ), XACT( * )
   REAL               RWORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   COMPLEX            CZERO
   PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
   INTEGER            NTYPES
   PARAMETER          ( NTYPES = 10 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 9 )
!     ..
!     .. Local Scalars ..
   LOGICAL            ZEROT
   CHARACTER          DIST, TYPE, UPLO, XTYPE
   CHARACTER*3        PATH, MATPATH
   INTEGER            I, I1, I2, IMAT, IN, INB, INFO, IOFF, IRHS, &
                      IUPLO, IZERO, J, K, KL, KU, LDA, LWORK, MODE, &
                      N, NB, NERRS, NFAIL, NIMAT, NRHS, NRUN, NT
   REAL               ANORM, CNDNUM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          UPLOS( 2 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   REAL               RESULT( NTESTS )
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAERH, ALAHD, ALASUM, CERRSY, CLACPY, CLARHS, &
                      CLATB4, CLATMS, CSYT02, CSYT01, &
                      CSYTRF_AA_2STAGE, CSYTRS_AA_2STAGE, &
                      XLAENV
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX, MIN
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
!     .. Data statements ..
   DATA               ISEEDY / 1988, 1989, 1990, 1991 /
   DATA               UPLOS / 'U', 'L' /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
!     Test path
!
   PATH( 1: 1 ) = 'Complex precision'
   PATH( 2: 3 ) = 'S2'
!
!     Path to generate matrices
!
   MATPATH( 1: 1 ) = 'Complex precision'
   MATPATH( 2: 3 ) = 'SY'
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
      CALL CERRSY( PATH, NOUT )
   INFOT = 0
!
!     Set the minimum block size for which the block routine should
!     be used, which will be later returned by ILAENV
!
   CALL XLAENV( 2, 2 )
!
!     Do for each value of N in NVAL
!
   DO IN = 1, NN
      N = NVAL( IN )
      IF( N  >  NMAX ) THEN
         NFAIL = NFAIL + 1
         WRITE(NOUT, 9995) 'M ', N, NMAX
         GO TO 180
      END IF
      LDA = MAX( N, 1 )
      XTYPE = 'N'
      NIMAT = NTYPES
      IF( N <= 0 ) &
         NIMAT = 1
!
      IZERO = 0
!
!        Do for each value of matrix type IMAT
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
!
!              Set up parameters with CLATB4 for the matrix generator
!              based on the type of matrix to be generated.
!
            CALL CLATB4( MATPATH, IMAT, N, N, TYPE, KL, KU, &
                         ANORM, MODE, CNDNUM, DIST )
!
!              Generate a matrix with CLATMS.
!
            SRNAMT = 'CLATMS'
            CALL CLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, &
                         CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, &
                         INFO )
!
!              Check error code from CLATMS and handle error.
!
            IF( INFO /= 0 ) THEN
               CALL ALAERH( PATH, 'CLATMS', INFO, 0, UPLO, N, N, -1, &
                            -1, -1, IMAT, NFAIL, NERRS, NOUT )
!
!                    Skip all tests for this generated matrix
!
               GO TO 160
            END IF
!
!              For matrix types 3-6, zero one or more rows and
!              columns of the matrix to test that INFO is returned
!              correctly.
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
!                    Set row and column IZERO to zero.
!
                  IF( IUPLO == 1 ) THEN
                     IOFF = ( IZERO-1 )*LDA
                     DO I = 1, IZERO - 1
                        A( IOFF+I ) = CZERO
                     ENDDO
                     IOFF = IOFF + IZERO
                     DO I = IZERO, N
                        A( IOFF ) = CZERO
                        IOFF = IOFF + LDA
                     ENDDO
                  ELSE
                     IOFF = IZERO
                     DO I = 1, IZERO - 1
                        A( IOFF ) = CZERO
                        IOFF = IOFF + LDA
                     ENDDO
                     IOFF = IOFF - IZERO
                     DO I = IZERO, N
                        A( IOFF+I ) = CZERO
                     ENDDO
                  END IF
               ELSE
                  IF( IUPLO == 1 ) THEN
!
!                       Set the first IZERO rows and columns to zero.
!
                     IOFF = 0
                     DO J = 1, N
                        I2 = MIN( J, IZERO )
                        DO I = 1, I2
                           A( IOFF+I ) = CZERO
                        ENDDO
                        IOFF = IOFF + LDA
                     ENDDO
                     IZERO = 1
                  ELSE
!
!                       Set the last IZERO rows and columns to zero.
!
                     IOFF = 0
                     DO J = 1, N
                        I1 = MAX( J, IZERO )
                        DO I = I1, N
                           A( IOFF+I ) = CZERO
                        ENDDO
                        IOFF = IOFF + LDA
                     ENDDO
                  END IF
               END IF
            ELSE
               IZERO = 0
            END IF
!
!              End generate the test matrix A.
!
!              Do for each value of NB in NBVAL
!
            DO INB = 1, NNB
!
!                 Set the optimal blocksize, which will be later
!                 returned by ILAENV.
!
               NB = NBVAL( INB )
               CALL XLAENV( 1, NB )
!
!                 Copy the test matrix A into matrix AFAC which
!                 will be factorized in place. This is needed to
!                 preserve the test matrix A for subsequent tests.
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL CLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Compute the L*D*L**T or U*D*U**T factorization of the
!                 matrix. IWORK stores details of the interchanges and
!                 the block structure of D. AINV is a work array for
!                 block factorization, LWORK is the length of AINV.
!
               SRNAMT = 'CSYTRF_AA_2STAGE'
               LWORK = MIN(N*NB, 3*NMAX*NMAX)
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL CSYTRF_AA_2STAGE( UPLO, N, AFAC, LDA, &
                                      AINV, (3*NB+1)*N, &
                                      IWORK, IWORK( 1+N ), &
                                      WORK, LWORK, &
                                      INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : CSYTRF_AA_2STAGE : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Adjust the expected value of INFO to account for
!                 pivoting.
!
               IF( IZERO > 0 ) THEN
                  J = 1
                  K = IZERO
  100                CONTINUE
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
!                 Check error code from CSYTRF and handle error.
!
               IF( INFO /= K ) THEN
                  CALL ALAERH( PATH, 'CSYTRF_AA_2STAGE', INFO, K, &
                               UPLO, N, N, -1, -1, NB, IMAT, NFAIL, &
                               NERRS, NOUT )
               END IF
!
!+    TEST 1
!                 Reconstruct matrix from factors and compute residual.
!
!                  CALL CSYT01_AA( UPLO, N, A, LDA, AFAC, LDA, IWORK,
!     $                            AINV, LDA, RWORK, RESULT( 1 ) )
!                  NT = 1
               NT = 0
!
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
               DO K = 1, NT
                  IF( RESULT( K ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9999 )UPLO, N, NB, IMAT, K, &
                        RESULT( K )
                     NFAIL = NFAIL + 1
                  END IF
                  ENDDO
               NRUN = NRUN + NT
!
!                 Skip solver test if INFO is not 0.
!
               IF( INFO /= 0 ) THEN
                  GO TO 140
               END IF
!
!                 Do for each value of NRHS in NSVAL.
!
               DO IRHS = 1, NNS
                  NRHS = NSVAL( IRHS )
!
!+    TEST 2 (Using TRS)
!                 Solve and compute residual for  A * X = B.
!
!                    Choose a set of NRHS random solution vectors
!                    stored in XACT and set up the right hand side B
!
                  SRNAMT = 'CLARHS'
                  CALL CLARHS( MATPATH, XTYPE, UPLO, ' ', N, N, &
                               KL, KU, NRHS, A, LDA, XACT, LDA, &
                               B, LDA, ISEED, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL CLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
                  SRNAMT = 'CSYTRS_AA_2STAGE'
                  LWORK = MAX( 1, 3*N-2 )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL CSYTRS_AA_2STAGE( UPLO, N, NRHS, AFAC, LDA, &
                               AINV, (3*NB+1)*N, IWORK, IWORK( 1+N ), &
                               X, LDA, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CSYTRS_AA_2STAGE : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from CSYTRS and handle error.
!
                  IF( INFO /= 0 ) THEN
                     IF( IZERO == 0 ) THEN
                        CALL ALAERH( PATH, 'CSYTRS_AA_2STAGE', &
                                     INFO, 0, UPLO, N, N, -1, -1, &
                                     NRHS, IMAT, NFAIL, NERRS, NOUT )
                     END IF
                  ELSE
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA &
                                  )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                       Compute the residual for the solution
!
                     CALL CSYT02( UPLO, N, NRHS, A, LDA, X, LDA, &
                                  WORK, LDA, RWORK, RESULT( 2 ) )
!
!
!                       Print information about the tests that did not pass
!                       the threshold.
!
                     DO K = 2, 2
                        IF( RESULT( K ) >= THRESH ) THEN
                           IF( NFAIL == 0 .AND. NERRS == 0 ) &
                              CALL ALAHD( NOUT, PATH )
                           WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, &
                              IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        END IF
                        ENDDO
                  END IF
                  NRUN = NRUN + 1
!
!                 End do for each value of NRHS in NSVAL.
!
                  ENDDO
  140             CONTINUE
               ENDDO
  160       CONTINUE
            ENDDO
  170    CONTINUE
         ENDDO
  180 CONTINUE
      ENDDO
!
!     Print a summary of the results.
!
   CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NB =', I4, ', type ', &
         I2, ', test ', I2, ', ratio =', G12.5 )
 9998 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NRHS=', I3, ', type ', &
         I2, ', test(', I2, ') =', G12.5 )
 9995 FORMAT( ' Invalid input value: ', A4, '=', I6, '; must be <=', &
         I6 )
   RETURN
!
!     End of CCHKSY_AA_2STAGE
!
END
                                                                                                                                                                                                                                                                                                            




