!> \brief \b CCHKGB
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CCHKGB( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NNS,
!                          NSVAL, THRESH, TSTERR, A, LA, AFAC, LAFAC, B,
!                          X, XACT, WORK, RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            LA, LAFAC, NM, NN, NNB, NNS, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), MVAL( * ), NBVAL( * ), NSVAL( * ),
!      $                   NVAL( * )
!       REAL               RWORK( * )
!       COMPLEX            A( * ), AFAC( * ), B( * ), WORK( * ), X( * ),
!      $                   XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CCHKGB tests CGBTRF, -TRS, -RFS, and -CON
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
!> \param[out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LA)
!> \endverbatim
!>
!> \param[in] LA
!> \verbatim
!>          LA is INTEGER
!>          The length of the array A.  LA >= (KLMAX+KUMAX+1)*NMAX
!>          where KLMAX is the largest entry in the local array KLVAL,
!>                KUMAX is the largest entry in the local array KUVAL and
!>                NMAX is the largest entry in the input array NVAL.
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is COMPLEX array, dimension (LAFAC)
!> \endverbatim
!>
!> \param[in] LAFAC
!> \verbatim
!>          LAFAC is INTEGER
!>          The length of the array AFAC. LAFAC >= (2*KLMAX+KUMAX+1)*NMAX
!>          where KLMAX is the largest entry in the local array KLVAL,
!>                KUMAX is the largest entry in the local array KUVAL and
!>                NMAX is the largest entry in the input array NVAL.
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX array, dimension (NMAX*NSMAX)
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
!>          WORK is COMPLEX array, dimension
!>                      (NMAX*max(3,NSMAX,NMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension
!>                      (NMAX+2*NSMAX)
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
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
!
!> \ingroup complex_lin
!
!  =====================================================================
   SUBROUTINE CCHKGB( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NNS, &
                      NSVAL, THRESH, TSTERR, A, LA, AFAC, LAFAC, B, &
                      X, XACT, WORK, RWORK, IWORK, NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            TSTERR
   INTEGER            LA, LAFAC, NM, NN, NNB, NNS, NOUT
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            IWORK( * ), MVAL( * ), NBVAL( * ), NSVAL( * ), &
                      NVAL( * )
   REAL               RWORK( * )
   COMPLEX            A( * ), AFAC( * ), B( * ), WORK( * ), X( * ), &
                      XACT( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ONE, ZERO
   PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
   INTEGER            NTYPES, NTESTS
   PARAMETER          ( NTYPES = 8, NTESTS = 7 )
   INTEGER            NBW, NTRAN
   PARAMETER          ( NBW = 4, NTRAN = 3 )
!     ..
!     .. Local Scalars ..
   LOGICAL            TRFCON, ZEROT
   CHARACTER          DIST, NORM, TRANS, TYPE, XTYPE
   CHARACTER*3        PATH
   INTEGER            I, I1, I2, IKL, IKU, IM, IMAT, IN, INB, INFO, &
                      IOFF, IRHS, ITRAN, IZERO, J, K, KL, KOFF, KU, &
                      LDA, LDAFAC, LDB, M, MODE, N, NB, NERRS, NFAIL, &
                      NIMAT, NKL, NKU, NRHS, NRUN
   REAL               AINVNM, ANORM, ANORMI, ANORMO, CNDNUM, RCOND, &
                      RCONDC, RCONDI, RCONDO
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          TRANSS( NTRAN )
   INTEGER            ISEED( 4 ), ISEEDY( 4 ), KLVAL( NBW ), &
                      KUVAL( NBW )
   REAL               RESULT( NTESTS )
!     ..
!     .. External Functions ..
   REAL               CLANGB, CLANGE, SGET06
   EXTERNAL           CLANGB, CLANGE, SGET06
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAERH, ALAHD, ALASUM, CCOPY, CERRGE, CGBCON, &
                      CGBRFS, CGBT01, CGBT02, CGBT05, CGBTRF, CGBTRS, &
                      CGET04, CLACPY, CLARHS, CLASET, CLATB4, CLATMS, &
                      XLAENV
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          CMPLX, MAX, MIN
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
   DATA               ISEEDY / 1988, 1989, 1990, 1991 / , &
                      TRANSS / 'N', 'T', 'C' /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   PATH( 1: 1 ) = 'Complex precision'
   PATH( 2: 3 ) = 'GB'
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
      CALL CERRGE( PATH, NOUT )
   INFOT = 0
!
!     Initialize the first value for the lower and upper bandwidths.
!
   KLVAL( 1 ) = 0
   KUVAL( 1 ) = 0
!
!     Do for each value of M in MVAL
!
   DO IM = 1, NM
      M = MVAL( IM )
!
!        Set values to use for the lower bandwidth.
!
      KLVAL( 2 ) = M + ( M+1 ) / 4
!
!        KLVAL( 2 ) = MAX( M-1, 0 )
!
      KLVAL( 3 ) = ( 3*M-1 ) / 4
      KLVAL( 4 ) = ( M+1 ) / 4
!
!        Do for each value of N in NVAL
!
      DO IN = 1, NN
         N = NVAL( IN )
         XTYPE = 'N'
!
!           Set values to use for the upper bandwidth.
!
         KUVAL( 2 ) = N + ( N+1 ) / 4
!
!           KUVAL( 2 ) = MAX( N-1, 0 )
!
         KUVAL( 3 ) = ( 3*N-1 ) / 4
         KUVAL( 4 ) = ( N+1 ) / 4
!
!           Set limits on the number of loop iterations.
!
         NKL = MIN( M+1, 4 )
         IF( N == 0 ) &
            NKL = 2
         NKU = MIN( N+1, 4 )
         IF( M == 0 ) &
            NKU = 2
         NIMAT = NTYPES
         IF( M <= 0 .OR. N <= 0 ) &
            NIMAT = 1
!
         DO IKL = 1, NKL
!
!              Do for KL = 0, (5*M+1)/4, (3M-1)/4, and (M+1)/4. This
!              order makes it easier to skip redundant values for small
!              values of M.
!
            KL = KLVAL( IKL )
            DO IKU = 1, NKU
!
!                 Do for KU = 0, (5*N+1)/4, (3N-1)/4, and (N+1)/4. This
!                 order makes it easier to skip redundant values for
!                 small values of N.
!
               KU = KUVAL( IKU )
!
!                 Check that A and AFAC are big enough to generate this
!                 matrix.
!
               LDA = KL + KU + 1
               LDAFAC = 2*KL + KU + 1
               IF( ( LDA*N ) > LA .OR. ( LDAFAC*N ) > LAFAC ) THEN
                  IF( NFAIL == 0 .AND. NERRS == 0 ) &
                     CALL ALAHD( NOUT, PATH )
                  IF( N*( KL+KU+1 ) > LA ) THEN
                     WRITE( NOUT, FMT = 9999 )LA, M, N, KL, KU, &
                        N*( KL+KU+1 )
                     NERRS = NERRS + 1
                  END IF
                  IF( N*( 2*KL+KU+1 ) > LAFAC ) THEN
                     WRITE( NOUT, FMT = 9998 )LAFAC, M, N, KL, KU, &
                        N*( 2*KL+KU+1 )
                     NERRS = NERRS + 1
                  END IF
                  GO TO 130
               END IF
!
               DO IMAT = 1, NIMAT
!
!                    Do the tests only if DOTYPE( IMAT ) is true.
!
                  IF( .NOT.DOTYPE( IMAT ) ) &
                     GO TO 120
!
!                    Skip types 2, 3, or 4 if the matrix size is too
!                    small.
!
                  ZEROT = IMAT >= 2 .AND. IMAT <= 4
                  IF( ZEROT .AND. N < IMAT-1 ) &
                     GO TO 120
!
                  IF( .NOT.ZEROT .OR. .NOT.DOTYPE( 1 ) ) THEN
!
!                       Set up parameters with CLATB4 and generate a
!                       test matrix with CLATMS.
!
                     CALL CLATB4( PATH, IMAT, M, N, TYPE, KL, KU, &
                                  ANORM, MODE, CNDNUM, DIST )
!
                     KOFF = MAX( 1, KU+2-N )
                     DO I = 1, KOFF - 1
                        A( I ) = ZERO
                     ENDDO
                     SRNAMT = 'CLATMS'
                     CALL CLATMS( M, N, DIST, ISEED, TYPE, RWORK, &
                                  MODE, CNDNUM, ANORM, KL, KU, 'Z', &
                                  A( KOFF ), LDA, WORK, INFO )
!
!                       Check the error code from CLATMS.
!
                     IF( INFO /= 0 ) THEN
                        CALL ALAERH( PATH, 'CLATMS', INFO, 0, ' ', M, &
                                     N, KL, KU, -1, IMAT, NFAIL, &
                                     NERRS, NOUT )
                        GO TO 120
                     END IF
                  ELSE IF( IZERO > 0 ) THEN
!
!                       Use the same matrix for types 3 and 4 as for
!                       type 2 by copying back the zeroed out column.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CCOPY( I2-I1+1, B, 1, A( IOFF+I1 ), 1 )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CCOPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                  END IF
!
!                    For types 2, 3, and 4, zero one or more columns of
!                    the matrix to test that INFO is returned correctly.
!
                  IZERO = 0
                  IF( ZEROT ) THEN
                     IF( IMAT == 2 ) THEN
                        IZERO = 1
                     ELSE IF( IMAT == 3 ) THEN
                        IZERO = MIN( M, N )
                     ELSE
                        IZERO = MIN( M, N ) / 2 + 1
                     END IF
                     IOFF = ( IZERO-1 )*LDA
                     IF( IMAT < 4 ) THEN
!
!                          Store the column to be zeroed out in B.
!
                        I1 = MAX( 1, KU+2-IZERO )
                        I2 = MIN( KL+KU+1, KU+1+( M-IZERO ) )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CCOPY( I2-I1+1, A( IOFF+I1 ), 1, B, 1 )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CCOPY : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
!
                        DO I = I1, I2
                           A( IOFF+I ) = ZERO
                        ENDDO
                     ELSE
                        DO J = IZERO, N
                           DO I = MAX( 1, KU+2-J ), &
                                   MIN( KL+KU+1, KU+1+( M-J ) )
                              A( IOFF+I ) = ZERO
                           ENDDO
                           IOFF = IOFF + LDA
                        ENDDO
                     END IF
                  END IF
!
!                    These lines, if used in place of the calls in the
!                    loop over INB, cause the code to bomb on a Sun
!                    SPARCstation.
!
!                     ANORMO = CLANGB( 'O', N, KL, KU, A, LDA, RWORK )
!                     ANORMI = CLANGB( 'I', N, KL, KU, A, LDA, RWORK )
!
!                    Do for each blocksize in NBVAL
!
                  DO INB = 1, NNB
                     NB = NBVAL( INB )
                     CALL XLAENV( 1, NB )
!
!                       Compute the LU factorization of the band matrix.
!
                     IF( M > 0 .AND. N > 0 )  THEN
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CLACPY( 'Full', KL+KU+1, N, A, LDA, &
                                     AFAC( KL+1 ), LDAFAC )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                     ENDIF
                     SRNAMT = 'CGBTRF'
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CGBTRF( M, N, KL, KU, AFAC, LDAFAC, IWORK, &
                                  INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CGBTRF : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                       Check error code from CGBTRF.
!
                     IF( INFO /= IZERO ) &
                        CALL ALAERH( PATH, 'CGBTRF', INFO, IZERO, &
                                     ' ', M, N, KL, KU, NB, IMAT, &
                                     NFAIL, NERRS, NOUT )
                     TRFCON = .FALSE.
!
!+    TEST 1
!                       Reconstruct matrix from factors and compute
!                       residual.
!
                     CALL CGBT01( M, N, KL, KU, A, LDA, AFAC, LDAFAC, &
                                  IWORK, WORK, RESULT( 1 ) )
!
!                       Print information about the tests so far that
!                       did not pass the threshold.
!
                     IF( RESULT( 1 ) >= THRESH ) THEN
                        IF( NFAIL == 0 .AND. NERRS == 0 ) &
                           CALL ALAHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9997 )M, N, KL, KU, NB, &
                           IMAT, 1, RESULT( 1 )
                        NFAIL = NFAIL + 1
                     END IF
                     NRUN = NRUN + 1
!
!                       Skip the remaining tests if this is not the
!                       first block size or if M .ne. N.
!
                     IF( INB > 1 .OR. M /= N ) &
                        GO TO 110
!
                     ANORMO = CLANGB( 'O', N, KL, KU, A, LDA, RWORK )
                     ANORMI = CLANGB( 'I', N, KL, KU, A, LDA, RWORK )
!
                     IF( INFO == 0 ) THEN
!
!                          Form the inverse of A so we can get a good
!                          estimate of CNDNUM = norm(A) * norm(inv(A)).
!
                        LDB = MAX( 1, N )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CLASET( 'Full', N, N, CMPLX( ZERO ), &
                                     CMPLX( ONE ), WORK, LDB )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CLASET : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                        SRNAMT = 'CGBTRS'
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CGBTRS( 'No transpose', N, KL, KU, N, &
                                     AFAC, LDAFAC, IWORK, WORK, LDB, &
                                     INFO )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CGBTRS : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
!
!                          Compute the 1-norm condition number of A.
!
                        AINVNM = CLANGE( 'O', N, N, WORK, LDB, &
                                 RWORK )
                        IF( ANORMO <= ZERO .OR. AINVNM <= ZERO ) THEN
                           RCONDO = ONE
                        ELSE
                           RCONDO = ( ONE / ANORMO ) / AINVNM
                        END IF
!
!                          Compute the infinity-norm condition number of
!                          A.
!
                        AINVNM = CLANGE( 'I', N, N, WORK, LDB, &
                                 RWORK )
                        IF( ANORMI <= ZERO .OR. AINVNM <= ZERO ) THEN
                           RCONDI = ONE
                        ELSE
                           RCONDI = ( ONE / ANORMI ) / AINVNM
                        END IF
                     ELSE
!
!                          Do only the condition estimate if INFO /= 0.
!
                        TRFCON = .TRUE.
                        RCONDO = ZERO
                        RCONDI = ZERO
                     END IF
!
!                       Skip the solve tests if the matrix is singular.
!
                     IF( TRFCON ) &
                        GO TO 90
!
                     DO IRHS = 1, NNS
                        NRHS = NSVAL( IRHS )
                        XTYPE = 'N'
!
                        DO ITRAN = 1, NTRAN
                           TRANS = TRANSS( ITRAN )
                           IF( ITRAN == 1 ) THEN
                              RCONDC = RCONDO
                              NORM = 'O'
                           ELSE
                              RCONDC = RCONDI
                              NORM = 'I'
                           END IF
!
!+    TEST 2:
!                             Solve and compute residual for op(A) * X = B.
!
                           SRNAMT = 'CLARHS'
                           CALL CLARHS( PATH, XTYPE, ' ', TRANS, N, &
                                        N, KL, KU, NRHS, A, LDA, &
                                        XACT, LDB, B, LDB, ISEED, &
                                        INFO )
                           XTYPE = 'C'
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                           CALL CLACPY( 'Full', N, NRHS, B, LDB, X, &
                                        LDB )
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S2_time)
                           open(file='results.out', unit=10, position = 'append')
                           write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                                 real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                           close(10)
#endif
!
                           SRNAMT = 'CGBTRS'
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                           CALL CGBTRS( TRANS, N, KL, KU, NRHS, AFAC, &
                                        LDAFAC, IWORK, X, LDB, INFO )
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S2_time)
                           open(file='results.out', unit=10, position = 'append')
                           write(10,'(A,F16.10,A)') 'Total time : CGBTRS : ',&
                                 real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                           close(10)
#endif
!
!                             Check error code from CGBTRS.
!
                           IF( INFO /= 0 ) &
                              CALL ALAERH( PATH, 'CGBTRS', INFO, 0, &
                                           TRANS, N, N, KL, KU, -1, &
                                           IMAT, NFAIL, NERRS, NOUT )
!
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                           CALL CLACPY( 'Full', N, NRHS, B, LDB, &
                                        WORK, LDB )
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S2_time)
                           open(file='results.out', unit=10, position = 'append')
                           write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                                 real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                           close(10)
#endif
                           CALL CGBT02( TRANS, M, N, KL, KU, NRHS, A, &
                                        LDA, X, LDB, WORK, LDB, &
                                        RWORK, RESULT( 2 ) )
!
!+    TEST 3:
!                             Check solution from generated exact
!                             solution.
!
                           CALL CGET04( N, NRHS, X, LDB, XACT, LDB, &
                                        RCONDC, RESULT( 3 ) )
!
!+    TESTS 4, 5, 6:
!                             Use iterative refinement to improve the
!                             solution.
!
                           SRNAMT = 'CGBRFS'
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                           CALL CGBRFS( TRANS, N, KL, KU, NRHS, A, &
                                        LDA, AFAC, LDAFAC, IWORK, B, &
                                        LDB, X, LDB, RWORK, &
                                        RWORK( NRHS+1 ), WORK, &
                                        RWORK( 2*NRHS+1 ), INFO )
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S2_time)
                           open(file='results.out', unit=10, position = 'append')
                           write(10,'(A,F16.10,A)') 'Total time : CGBRFS : ',&
                                 real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                           close(10)
#endif
!
!                             Check error code from CGBRFS.
!
                           IF( INFO /= 0 ) &
                              CALL ALAERH( PATH, 'CGBRFS', INFO, 0, &
                                           TRANS, N, N, KL, KU, NRHS, &
                                           IMAT, NFAIL, NERRS, NOUT )
!
                           CALL CGET04( N, NRHS, X, LDB, XACT, LDB, &
                                        RCONDC, RESULT( 4 ) )
                           CALL CGBT05( TRANS, N, KL, KU, NRHS, A, &
                                        LDA, B, LDB, X, LDB, XACT, &
                                        LDB, RWORK, RWORK( NRHS+1 ), &
                                        RESULT( 5 ) )
!
!                             Print information about the tests that did
!                             not pass the threshold.
!
                           DO K = 2, 6
                              IF( RESULT( K ) >= THRESH ) THEN
                                 IF( NFAIL == 0 .AND. NERRS == 0 ) &
                                    CALL ALAHD( NOUT, PATH )
                                 WRITE( NOUT, FMT = 9996 )TRANS, N, &
                                    KL, KU, NRHS, IMAT, K, &
                                    RESULT( K )
                                 NFAIL = NFAIL + 1
                              END IF
                           ENDDO
                           NRUN = NRUN + 5
                        ENDDO
                     ENDDO
!
!+    TEST 7:
!                          Get an estimate of RCOND = 1/CNDNUM.
!
90                   CONTINUE
                     DO ITRAN = 1, 2
                        IF( ITRAN == 1 ) THEN
                           ANORM = ANORMO
                           RCONDC = RCONDO
                           NORM = 'O'
                        ELSE
                           ANORM = ANORMI
                           RCONDC = RCONDI
                           NORM = 'I'
                        END IF
                        SRNAMT = 'CGBCON'
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CGBCON( NORM, N, KL, KU, AFAC, LDAFAC, &
                                     IWORK, ANORM, RCOND, WORK, &
                                     RWORK, INFO )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CGBCON : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
!
!                             Check error code from CGBCON.
!
                        IF( INFO /= 0 ) &
                           CALL ALAERH( PATH, 'CGBCON', INFO, 0, &
                                        NORM, N, N, KL, KU, -1, IMAT, &
                                        NFAIL, NERRS, NOUT )
!
                        RESULT( 7 ) = SGET06( RCOND, RCONDC )
!
!                          Print information about the tests that did
!                          not pass the threshold.
!
                        IF( RESULT( 7 ) >= THRESH ) THEN
                           IF( NFAIL == 0 .AND. NERRS == 0 ) &
                              CALL ALAHD( NOUT, PATH )
                           WRITE( NOUT, FMT = 9995 )NORM, N, KL, KU, &
                              IMAT, 7, RESULT( 7 )
                           NFAIL = NFAIL + 1
                        END IF
                        NRUN = NRUN + 1
                        ENDDO
  110                CONTINUE
                     ENDDO
  120             CONTINUE
                  ENDDO
  130          CONTINUE
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( ' *** In CCHKGB, LA=', I5, ' is too small for M=', I5, &
         ', N=', I5, ', KL=', I4, ', KU=', I4, &
         / ' ==> Increase LA to at least ', I5 )
 9998 FORMAT( ' *** In CCHKGB, LAFAC=', I5, ' is too small for M=', I5, &
         ', N=', I5, ', KL=', I4, ', KU=', I4, &
         / ' ==> Increase LAFAC to at least ', I5 )
 9997 FORMAT( ' M =', I5, ', N =', I5, ', KL=', I5, ', KU=', I5, &
         ', NB =', I4, ', type ', I1, ', test(', I1, ')=', G12.5 )
 9996 FORMAT( ' TRANS=''', A1, ''', N=', I5, ', KL=', I5, ', KU=', I5, &
         ', NRHS=', I3, ', type ', I1, ', test(', I1, ')=', G12.5 )
 9995 FORMAT( ' NORM =''', A1, ''', N=', I5, ', KL=', I5, ', KU=', I5, &
         ',', 10X, ' type ', I1, ', test(', I1, ')=', G12.5 )
!
   RETURN
!
!     End of CCHKGB
!
END
                                                                                                                                                                                                                                                                                                            




