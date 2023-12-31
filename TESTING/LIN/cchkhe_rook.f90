!> \brief \b CCHKHE_ROOK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CCHKHE_ROOK( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL,
!                               THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X,
!                               XACT, WORK, RWORK, IWORK, NOUT )
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
!> CCHKHE_ROOK tests CHETRF_ROOK, -TRI_ROOK, -TRS_ROOK,
!> and -CON_ROOK.
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
!>          RWORK is REAL array, dimension (max(NMAX,2*NSMAX)
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
   SUBROUTINE CCHKHE_ROOK( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL, &
                           THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, &
                           XACT, WORK, RWORK, IWORK, NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            TSTERR
   INTEGER            NMAX, NN, NNB, NNS, NOUT
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * )
   REAL               RWORK( * )
   COMPLEX            A( * ), AFAC( * ), AINV( * ), B( * ), &
                      WORK( * ), X( * ), XACT( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO, ONE
   PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
   REAL               ONEHALF
   PARAMETER          ( ONEHALF = 0.5E+0 )
   REAL               EIGHT, SEVTEN
   PARAMETER          ( EIGHT = 8.0E+0, SEVTEN = 17.0E+0 )
   COMPLEX            CZERO
   PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
   INTEGER            NTYPES
   PARAMETER          ( NTYPES = 10 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 7 )
!     ..
!     .. Local Scalars ..
   LOGICAL            TRFCON, ZEROT
   CHARACTER          DIST, TYPE, UPLO, XTYPE
   CHARACTER*3        PATH, MATPATH
   INTEGER            I, I1, I2, IMAT, IN, INB, INFO, IOFF, IRHS, &
                      IUPLO, IZERO, J, K, KL, KU, LDA, LWORK, MODE, &
                      N, NB, NERRS, NFAIL, NIMAT, NRHS, NRUN, NT
   REAL               ALPHA, ANORM, CNDNUM, CONST, SING_MAX, &
                      SING_MIN, RCOND, RCONDC, STEMP
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          UPLOS( 2 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   REAL               RESULT( NTESTS )
   COMPLEX            BLOCK( 2, 2 ), CDUMMY( 1 )
!     ..
!     .. External Functions ..
   REAL               CLANGE, CLANHE, SGET06
   EXTERNAL           CLANGE, CLANHE, SGET06
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAERH, ALAHD, ALASUM, CERRHE, CGESVD, CGET04, &
                      CLACPY, CLARHS, CLATB4, CLATMS, CPOT02, &
                      CPOT03, CHECON_ROOK, CHET01_ROOK, CHETRF_ROOK, &
                      CHETRI_ROOK, CHETRS_ROOK, XLAENV
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          CONJG, MAX, MIN, SQRT
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
   ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
!
!     Test path
!
   PATH( 1: 1 ) = 'Complex precision'
   PATH( 2: 3 ) = 'HR'
!
!     Path to generate matrices
!
   MATPATH( 1: 1 ) = 'Complex precision'
   MATPATH( 2: 3 ) = 'HE'
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
      CALL CERRHE( PATH, NOUT )
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
            GO TO 260
!
!           Skip types 3, 4, 5, or 6 if the matrix size is too small.
!
         ZEROT = IMAT >= 3 .AND. IMAT <= 6
         IF( ZEROT .AND. N < IMAT-2 ) &
            GO TO 260
!
!           Do first for UPLO = 'U', then for UPLO = 'L'
!
         DO IUPLO = 1, 2
            UPLO = UPLOS( IUPLO )
!
!                 Begin generate the test matrix A.
!
!                 Set up parameters with CLATB4 for the matrix generator
!                 based on the type of matrix to be generated.
!
               CALL CLATB4( MATPATH, IMAT, N, N, TYPE, KL, KU, ANORM, &
                            MODE, CNDNUM, DIST )
!
!                 Generate a matrix with CLATMS.
!
               SRNAMT = 'CLATMS'
               CALL CLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, &
                            CNDNUM, ANORM, KL, KU, UPLO, A, LDA, &
                            WORK, INFO )
!
!                 Check error code from CLATMS and handle error.
!
               IF( INFO /= 0 ) THEN
                  CALL ALAERH( PATH, 'CLATMS', INFO, 0, UPLO, N, N, &
                               -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
!
!                    Skip all tests for this generated matrix
!
                  GO TO 250
               END IF
!
!                 For matrix types 3-6, zero one or more rows and
!                 columns of the matrix to test that INFO is returned
!                 correctly.
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
!                          Set the first IZERO rows and columns to zero.
!
                        IOFF = 0
                        DO J = 1, N
                           I2 = MIN( J, IZERO )
                           DO I = 1, I2
                              A( IOFF+I ) = CZERO
                           ENDDO
                           IOFF = IOFF + LDA
                        ENDDO
                     ELSE
!
!                          Set the last IZERO rows and columns to zero.
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
!                 End generate the test matrix A.
!
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
               LWORK = MAX( 2, NB )*LDA
               SRNAMT = 'CHETRF_ROOK'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL CHETRF_ROOK( UPLO, N, AFAC, LDA, IWORK, AINV, &
                                 LWORK, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : CHETRF_ROOK : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Adjust the expected value of INFO to account for
!                 pivoting.
!
               K = IZERO
               IF( K > 0 ) THEN
  100                CONTINUE
                  IF( IWORK( K ) < 0 ) THEN
                     IF( IWORK( K ) /= -K ) THEN
                        K = -IWORK( K )
                        GO TO 100
                     END IF
                  ELSE IF( IWORK( K ) /= K ) THEN
                     K = IWORK( K )
                     GO TO 100
                  END IF
               END IF
!
!                 Check error code from CHETRF_ROOK and handle error.
!
               IF( INFO /= K) &
                  CALL ALAERH( PATH, 'CHETRF_ROOK', INFO, K, &
                               UPLO, N, N, -1, -1, NB, IMAT, &
                               NFAIL, NERRS, NOUT )
!
!                 Set the condition estimate flag if the INFO is not 0.
!
               IF( INFO /= 0 ) THEN
                  TRFCON = .TRUE.
               ELSE
                  TRFCON = .FALSE.
               END IF
!
!+    TEST 1
!                 Reconstruct matrix from factors and compute residual.
!
               CALL CHET01_ROOK( UPLO, N, A, LDA, AFAC, LDA, IWORK, &
                                 AINV, LDA, RWORK, RESULT( 1 ) )
               NT = 1
!
!+    TEST 2
!                 Form the inverse and compute the residual,
!                 if the factorization was competed without INFO > 0
!                 (i.e. there is no zero rows and columns).
!                 Do it only for the first block size.
!
               IF( INB == 1 .AND. .NOT.TRFCON ) THEN
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL CLACPY( UPLO, N, N, AFAC, LDA, AINV, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  SRNAMT = 'CHETRI_ROOK'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL CHETRI_ROOK( UPLO, N, AINV, LDA, IWORK, WORK, &
                                    INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CHETRI_ROOK : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from CHETRI_ROOK and handle error.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'CHETRI_ROOK', INFO, -1, &
                                  UPLO, N, N, -1, -1, -1, IMAT, &
                                  NFAIL, NERRS, NOUT )
!
!                    Compute the residual for a Hermitian matrix times
!                    its inverse.
!
                  CALL CPOT03( UPLO, N, A, LDA, AINV, LDA, WORK, LDA, &
                               RWORK, RCONDC, RESULT( 2 ) )
                  NT = 2
               END IF
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
!+    TEST 3
!                 Compute largest element in U or L
!
               RESULT( 3 ) = ZERO
               STEMP = ZERO
!
               CONST = ( ( ALPHA**2-ONE ) / ( ALPHA**2-ONEHALF ) ) / &
                       ( ONE-ALPHA )
!
               IF( IUPLO == 1 ) THEN
!
!                 Compute largest element in U
!
                  K = N
  120                CONTINUE
                  IF( K <= 1 ) &
                     GO TO 130
!
                  IF( IWORK( K ) > ZERO ) THEN
!
!                       Get max absolute value from elements
!                       in column k in U
!
                     STEMP = CLANGE( 'M', K-1, 1, &
                             AFAC( ( K-1 )*LDA+1 ), LDA, RWORK )
                  ELSE
!
!                       Get max absolute value from elements
!                       in columns k and k-1 in U
!
                     STEMP = CLANGE( 'M', K-2, 2, &
                             AFAC( ( K-2 )*LDA+1 ), LDA, RWORK )
                     K = K - 1
!
                  END IF
!
!                    STEMP should be bounded by CONST
!
                  STEMP = STEMP - CONST + THRESH
                  IF( STEMP > RESULT( 3 ) ) &
                     RESULT( 3 ) = STEMP
!
                  K = K - 1
!
                  GO TO 120
  130                CONTINUE
!
               ELSE
!
!                 Compute largest element in L
!
                  K = 1
  140                CONTINUE
                  IF( K >= N ) &
                     GO TO 150
!
                  IF( IWORK( K ) > ZERO ) THEN
!
!                       Get max absolute value from elements
!                       in column k in L
!
                     STEMP = CLANGE( 'M', N-K, 1, &
                             AFAC( ( K-1 )*LDA+K+1 ), LDA, RWORK )
                  ELSE
!
!                       Get max absolute value from elements
!                       in columns k and k+1 in L
!
                     STEMP = CLANGE( 'M', N-K-1, 2, &
                             AFAC( ( K-1 )*LDA+K+2 ), LDA, RWORK )
                     K = K + 1
!
                  END IF
!
!                    STEMP should be bounded by CONST
!
                  STEMP = STEMP - CONST + THRESH
                  IF( STEMP > RESULT( 3 ) ) &
                     RESULT( 3 ) = STEMP
!
                  K = K + 1
!
                  GO TO 140
  150                CONTINUE
               END IF
!
!
!+    TEST 4
!                 Compute largest 2-Norm (condition number)
!                 of 2-by-2 diag blocks
!
               RESULT( 4 ) = ZERO
               STEMP = ZERO
!
               CONST = ( ( ALPHA**2-ONE ) / ( ALPHA**2-ONEHALF ) )* &
                       ( ( ONE + ALPHA ) / ( ONE - ALPHA ) )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL CLACPY( UPLO, N, N, AFAC, LDA, AINV, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
               IF( IUPLO == 1 ) THEN
!
!                    Loop backward for UPLO = 'U'
!
                  K = N
  160                CONTINUE
                  IF( K <= 1 ) &
                     GO TO 170
!
                  IF( IWORK( K ) < ZERO ) THEN
!
!                       Get the two singular values
!                       (real and non-negative) of a 2-by-2 block,
!                       store them in RWORK array
!
                     BLOCK( 1, 1 ) = AFAC( ( K-2 )*LDA+K-1 )
                     BLOCK( 1, 2 ) = AFAC( (K-1)*LDA+K-1 )
                     BLOCK( 2, 1 ) = CONJG( BLOCK( 1, 2 ) )
                     BLOCK( 2, 2 ) = AFAC( (K-1)*LDA+K )
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CGESVD( 'N', 'N', 2, 2, BLOCK, 2, RWORK, &
                                  CDUMMY, 1, CDUMMY, 1, &
                                  WORK, 6, RWORK( 3 ), INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CGESVD : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!
                     SING_MAX = RWORK( 1 )
                     SING_MIN = RWORK( 2 )
!
                     STEMP = SING_MAX / SING_MIN
!
!                       STEMP should be bounded by CONST
!
                     STEMP = STEMP - CONST + THRESH
                     IF( STEMP > RESULT( 4 ) ) &
                        RESULT( 4 ) = STEMP
                     K = K - 1
!
                  END IF
!
                  K = K - 1
!
                  GO TO 160
  170                CONTINUE
!
               ELSE
!
!                    Loop forward for UPLO = 'L'
!
                  K = 1
  180                CONTINUE
                  IF( K >= N ) &
                     GO TO 190
!
                  IF( IWORK( K ) < ZERO ) THEN
!
!                       Get the two singular values
!                       (real and non-negative) of a 2-by-2 block,
!                       store them in RWORK array
!
                     BLOCK( 1, 1 ) = AFAC( ( K-1 )*LDA+K )
                     BLOCK( 2, 1 ) = AFAC( ( K-1 )*LDA+K+1 )
                     BLOCK( 1, 2 ) = CONJG( BLOCK( 2, 1 ) )
                     BLOCK( 2, 2 ) = AFAC( K*LDA+K+1 )
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CGESVD( 'N', 'N', 2, 2, BLOCK, 2, RWORK, &
                                  CDUMMY, 1, CDUMMY, 1, &
                                  WORK, 6, RWORK(3), INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CGESVD : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
                     SING_MAX = RWORK( 1 )
                     SING_MIN = RWORK( 2 )
!
                     STEMP = SING_MAX / SING_MIN
!
!                       STEMP should be bounded by CONST
!
                     STEMP = STEMP - CONST + THRESH
                     IF( STEMP > RESULT( 4 ) ) &
                        RESULT( 4 ) = STEMP
                     K = K + 1
!
                  END IF
!
                  K = K + 1
!
                  GO TO 180
  190                CONTINUE
               END IF
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
               DO K = 3, 4
                  IF( RESULT( K ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9999 )UPLO, N, NB, IMAT, K, &
                        RESULT( K )
                     NFAIL = NFAIL + 1
                  END IF
                  ENDDO
               NRUN = NRUN + 2
!
!                 Skip the other tests if this is not the first block
!                 size.
!
               IF( INB > 1 ) &
                  GO TO 240
!
!                 Do only the condition estimate if INFO is not 0.
!
               IF( TRFCON ) THEN
                  RCONDC = ZERO
                  GO TO 230
               END IF
!
!                 Do for each value of NRHS in NSVAL.
!
               DO IRHS = 1, NNS
                  NRHS = NSVAL( IRHS )
!
!                    Begin loop over NRHS values
!
!
!+    TEST 5 ( Using TRS_ROOK)
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
                  SRNAMT = 'CHETRS_ROOK'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL CHETRS_ROOK( UPLO, N, NRHS, AFAC, LDA, IWORK, &
                                    X, LDA, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CHETRS_ROOK : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from CHETRS_ROOK and handle error.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'CHETRS_ROOK', INFO, 0, &
                                  UPLO, N, N, -1, -1, NRHS, IMAT, &
                                  NFAIL, NERRS, NOUT )
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL CLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Compute the residual for the solution
!
                  CALL CPOT02( UPLO, N, NRHS, A, LDA, X, LDA, WORK, &
                               LDA, RWORK, RESULT( 5 ) )
!
!+    TEST 6
!                 Check solution from generated exact solution.
!
                  CALL CGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                               RESULT( 6 ) )
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                  DO K = 5, 6
                     IF( RESULT( K ) >= THRESH ) THEN
                        IF( NFAIL == 0 .AND. NERRS == 0 ) &
                           CALL ALAHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, &
                           IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
                     ENDDO
                  NRUN = NRUN + 2
!
!                 End do for each value of NRHS in NSVAL.
!
                  ENDDO
!
!+    TEST 7
!                 Get an estimate of RCOND = 1/CNDNUM.
!
  230             CONTINUE
               ANORM = CLANHE( '1', UPLO, N, A, LDA, RWORK )
               SRNAMT = 'CHECON_ROOK'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL CHECON_ROOK( UPLO, N, AFAC, LDA, IWORK, ANORM, &
                                 RCOND, WORK, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : CHECON_ROOK : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Check error code from CHECON_ROOK and handle error.
!
               IF( INFO /= 0 ) &
                  CALL ALAERH( PATH, 'CHECON_ROOK', INFO, 0, &
                                UPLO, N, N, -1, -1, -1, IMAT, &
                                NFAIL, NERRS, NOUT )
!
!                 Compute the test ratio to compare values of RCOND
!
               RESULT( 7 ) = SGET06( RCOND, RCONDC )
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
               IF( RESULT( 7 ) >= THRESH ) THEN
                  IF( NFAIL == 0 .AND. NERRS == 0 ) &
                     CALL ALAHD( NOUT, PATH )
                  WRITE( NOUT, FMT = 9997 )UPLO, N, IMAT, 7, &
                     RESULT( 7 )
                  NFAIL = NFAIL + 1
               END IF
               NRUN = NRUN + 1
  240          CONTINUE
               ENDDO
!
  250       CONTINUE
            ENDDO
  260    CONTINUE
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NB =', I4, ', type ', &
         I2, ', test ', I2, ', ratio =', G12.5 )
 9998 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NRHS=', I3, ', type ', &
         I2, ', test ', I2, ', ratio =', G12.5 )
 9997 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ',', 10X, ' type ', I2, &
         ', test ', I2, ', ratio =', G12.5 )
   RETURN
!
!     End of CCHKHE_ROOK
!
END
                                                                                                                                                                                                                                                                                                            




