!> \brief \b CDRVGBX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CDRVGB( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, LA,
!                          AFB, LAFB, ASAV, B, BSAV, X, XACT, S, WORK,
!                          RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            LA, LAFB, NN, NOUT, NRHS
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NVAL( * )
!       REAL               RWORK( * ), S( * )
!       COMPLEX            A( * ), AFB( * ), ASAV( * ), B( * ), BSAV( * ),
!      $                   WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CDRVGB tests the driver routines CGBSV, -SVX, and -SVXX.
!>
!> Note that this file is used only when the XBLAS are available,
!> otherwise cdrvgb.f defines this subroutine.
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
!>          The values of the matrix column dimension N.
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
!> \param[out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LA)
!> \endverbatim
!>
!> \param[in] LA
!> \verbatim
!>          LA is INTEGER
!>          The length of the array A.  LA >= (2*NMAX-1)*NMAX
!>          where NMAX is the largest entry in NVAL.
!> \endverbatim
!>
!> \param[out] AFB
!> \verbatim
!>          AFB is COMPLEX array, dimension (LAFB)
!> \endverbatim
!>
!> \param[in] LAFB
!> \verbatim
!>          LAFB is INTEGER
!>          The length of the array AFB.  LAFB >= (3*NMAX-2)*NMAX
!>          where NMAX is the largest entry in NVAL.
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is COMPLEX array, dimension (LA)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] BSAV
!> \verbatim
!>          BSAV is COMPLEX array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is COMPLEX array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL array, dimension (2*NMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension
!>                      (NMAX*max(3,NRHS,NMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension
!>                      (max(2*NMAX,NMAX+2*NRHS))
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
   SUBROUTINE CDRVGB( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, LA, &
                      AFB, LAFB, ASAV, B, BSAV, X, XACT, S, WORK, &
                      RWORK, IWORK, NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            TSTERR
   INTEGER            LA, LAFB, NN, NOUT, NRHS
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            IWORK( * ), NVAL( * )
   REAL               RWORK( * ), S( * )
   COMPLEX            A( * ), AFB( * ), ASAV( * ), B( * ), BSAV( * ), &
                      WORK( * ), X( * ), XACT( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ONE, ZERO
   PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
   INTEGER            NTYPES
   PARAMETER          ( NTYPES = 8 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 7 )
   INTEGER            NTRAN
   PARAMETER          ( NTRAN = 3 )
!     ..
!     .. Local Scalars ..
   LOGICAL            EQUIL, NOFACT, PREFAC, TRFCON, ZEROT
   CHARACTER          DIST, EQUED, FACT, TRANS, TYPE, XTYPE
   CHARACTER*3        PATH
   INTEGER            I, I1, I2, IEQUED, IFACT, IKL, IKU, IMAT, IN, &
                      INFO, IOFF, ITRAN, IZERO, J, K, K1, KL, KU, &
                      LDA, LDAFB, LDB, MODE, N, NB, NBMIN, NERRS, &
                      NFACT, NFAIL, NIMAT, NKL, NKU, NRUN, NT, &
                      N_ERR_BNDS
   REAL               AINVNM, AMAX, ANORM, ANORMI, ANORMO, ANRMPV, &
                      CNDNUM, COLCND, RCOND, RCONDC, RCONDI, RCONDO, &
                      ROLDC, ROLDI, ROLDO, ROWCND, RPVGRW, &
                      RPVGRW_SVXX
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          EQUEDS( 4 ), FACTS( 3 ), TRANSS( NTRAN )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   REAL               RDUM( 1 ), RESULT( NTESTS ), BERR( NRHS ), &
                      ERRBNDS_N( NRHS,3 ), ERRBNDS_C( NRHS, 3 )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               CLANGB, CLANGE, CLANTB, SGET06, SLAMCH, &
                      CLA_GBRPVGRW
   EXTERNAL           LSAME, CLANGB, CLANGE, CLANTB, SGET06, SLAMCH, &
                      CLA_GBRPVGRW
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALADHD, ALAERH, ALASVM, CERRVX, CGBEQU, CGBSV, &
                      CGBSVX, CGBT01, CGBT02, CGBT05, CGBTRF, CGBTRS, &
                      CGET04, CLACPY, CLAQGB, CLARHS, CLASET, CLATB4, &
                      CLATMS, XLAENV, CGBSVXX
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, CMPLX, MAX, MIN
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
   DATA               TRANSS / 'N', 'T', 'C' /
   DATA               FACTS / 'F', 'N', 'E' /
   DATA               EQUEDS / 'N', 'R', 'C', 'B' /
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
      CALL CERRVX( PATH, NOUT )
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
      LDB = MAX( N, 1 )
      XTYPE = 'N'
!
!        Set limits on the number of loop iterations.
!
      NKL = MAX( 1, MIN( N, 4 ) )
      IF( N == 0 ) &
         NKL = 1
      NKU = NKL
      NIMAT = NTYPES
      IF( N <= 0 ) &
         NIMAT = 1
!
      DO IKL = 1, NKL
!
!           Do for KL = 0, N-1, (3N-1)/4, and (N+1)/4. This order makes
!           it easier to skip redundant values for small values of N.
!
         IF( IKL == 1 ) THEN
            KL = 0
         ELSE IF( IKL == 2 ) THEN
            KL = MAX( N-1, 0 )
         ELSE IF( IKL == 3 ) THEN
            KL = ( 3*N-1 ) / 4
         ELSE IF( IKL == 4 ) THEN
            KL = ( N+1 ) / 4
         END IF
         DO IKU = 1, NKU
!
!              Do for KU = 0, N-1, (3N-1)/4, and (N+1)/4. This order
!              makes it easier to skip redundant values for small
!              values of N.
!
            IF( IKU == 1 ) THEN
               KU = 0
            ELSE IF( IKU == 2 ) THEN
               KU = MAX( N-1, 0 )
            ELSE IF( IKU == 3 ) THEN
               KU = ( 3*N-1 ) / 4
            ELSE IF( IKU == 4 ) THEN
               KU = ( N+1 ) / 4
            END IF
!
!              Check that A and AFB are big enough to generate this
!              matrix.
!
            LDA = KL + KU + 1
            LDAFB = 2*KL + KU + 1
            IF( LDA*N > LA .OR. LDAFB*N > LAFB ) THEN
               IF( NFAIL == 0 .AND. NERRS == 0 ) &
                  CALL ALADHD( NOUT, PATH )
               IF( LDA*N > LA ) THEN
                  WRITE( NOUT, FMT = 9999 )LA, N, KL, KU, &
                     N*( KL+KU+1 )
                  NERRS = NERRS + 1
               END IF
               IF( LDAFB*N > LAFB ) THEN
                  WRITE( NOUT, FMT = 9998 )LAFB, N, KL, KU, &
                     N*( 2*KL+KU+1 )
                  NERRS = NERRS + 1
               END IF
               GO TO 130
            END IF
!
            DO IMAT = 1, NIMAT
!
!                 Do the tests only if DOTYPE( IMAT ) is true.
!
               IF( .NOT.DOTYPE( IMAT ) ) &
                  GO TO 120
!
!                 Skip types 2, 3, or 4 if the matrix is too small.
!
               ZEROT = IMAT >= 2 .AND. IMAT <= 4
               IF( ZEROT .AND. N < IMAT-1 ) &
                  GO TO 120
!
!                 Set up parameters with CLATB4 and generate a
!                 test matrix with CLATMS.
!
               CALL CLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, &
                            MODE, CNDNUM, DIST )
               RCONDC = ONE / CNDNUM
!
               SRNAMT = 'CLATMS'
               CALL CLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, &
                            CNDNUM, ANORM, KL, KU, 'Z', A, LDA, WORK, &
                            INFO )
!
!                 Check the error code from CLATMS.
!
               IF( INFO /= 0 ) THEN
                  CALL ALAERH( PATH, 'CLATMS', INFO, 0, ' ', N, N, &
                               KL, KU, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 120
               END IF
!
!                 For types 2, 3, and 4, zero one or more columns of
!                 the matrix to test that INFO is returned correctly.
!
               IZERO = 0
               IF( ZEROT ) THEN
                  IF( IMAT == 2 ) THEN
                     IZERO = 1
                  ELSE IF( IMAT == 3 ) THEN
                     IZERO = N
                  ELSE
                     IZERO = N / 2 + 1
                  END IF
                  IOFF = ( IZERO-1 )*LDA
                  IF( IMAT < 4 ) THEN
                     I1 = MAX( 1, KU+2-IZERO )
                     I2 = MIN( KL+KU+1, KU+1+( N-IZERO ) )
                     DO I = I1, I2
                        A( IOFF+I ) = ZERO
                     ENDDO
                  ELSE
                     DO J = IZERO, N
                        DO I = MAX( 1, KU+2-J ), &
                                MIN( KL+KU+1, KU+1+( N-J ) )
                           A( IOFF+I ) = ZERO
                        ENDDO
                        IOFF = IOFF + LDA
                     ENDDO
                  END IF
               END IF
!
!                 Save a copy of the matrix A in ASAV.
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL CLACPY( 'Full', KL+KU+1, N, A, LDA, ASAV, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
               DO IEQUED = 1, 4
                  EQUED = EQUEDS( IEQUED )
                  IF( IEQUED == 1 ) THEN
                     NFACT = 3
                  ELSE
                     NFACT = 1
                  END IF
!
                  DO IFACT = 1, NFACT
                     FACT = FACTS( IFACT )
                     PREFAC = LSAME( FACT, 'F' )
                     NOFACT = LSAME( FACT, 'N' )
                     EQUIL = LSAME( FACT, 'E' )
!
                     IF( ZEROT ) THEN
                        IF( PREFAC ) &
                           GO TO 100
                        RCONDO = ZERO
                        RCONDI = ZERO
!
                     ELSE IF( .NOT.NOFACT ) THEN
!
!                          Compute the condition number for comparison
!                          with the value returned by SGESVX (FACT =
!                          'N' reuses the condition number from the
!                          previous iteration with FACT = 'F').
!
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CLACPY( 'Full', KL+KU+1, N, ASAV, LDA, &
                                     AFB( KL+1 ), LDAFB )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                        IF( EQUIL .OR. IEQUED > 1 ) THEN
!
!                             Compute row and column scale factors to
!                             equilibrate the matrix A.
!
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                           CALL CGBEQU( N, N, KL, KU, AFB( KL+1 ), &
                                        LDAFB, S, S( N+1 ), ROWCND, &
                                        COLCND, AMAX, INFO )
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S2_time)
                           open(file='results.out', unit=10, position = 'append')
                           write(10,'(A,F16.10,A)') 'Total time : CGBEQU : ',&
                                 real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                           close(10)
#endif
                           IF( INFO == 0 .AND. N > 0 ) THEN
                              IF( LSAME( EQUED, 'R' ) ) THEN
                                 ROWCND = ZERO
                                 COLCND = ONE
                              ELSE IF( LSAME( EQUED, 'C' ) ) THEN
                                 ROWCND = ONE
                                 COLCND = ZERO
                              ELSE IF( LSAME( EQUED, 'B' ) ) THEN
                                 ROWCND = ZERO
                                 COLCND = ZERO
                              END IF
!
!                                Equilibrate the matrix.
!
#ifdef _TIMER
                              call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                              CALL CLAQGB( N, N, KL, KU, AFB( KL+1 ), &
                                           LDAFB, S, S( N+1 ), &
                                           ROWCND, COLCND, AMAX, &
                                           EQUED )
#ifdef _TIMER
                              call system_clock(count_rate=nb_periods_sec,count=S2_time)
                              open(file='results.out', unit=10, position = 'append')
                              write(10,'(A,F16.10,A)') 'Total time : CLAQGB : ',&
                                    real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                              close(10)
#endif
                           END IF
                        END IF
!
!                          Save the condition number of the
!                          non-equilibrated system for use in CGET04.
!
                        IF( EQUIL ) THEN
                           ROLDO = RCONDO
                           ROLDI = RCONDI
                        END IF
!
!                          Compute the 1-norm and infinity-norm of A.
!
                        ANORMO = CLANGB( '1', N, KL, KU, AFB( KL+1 ), &
                                 LDAFB, RWORK )
                        ANORMI = CLANGB( 'I', N, KL, KU, AFB( KL+1 ), &
                                 LDAFB, RWORK )
!
!                          Factor the matrix A.
!
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CGBTRF( N, N, KL, KU, AFB, LDAFB, IWORK, &
                                     INFO )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CGBTRF : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
!
!                          Form the inverse of A.
!
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
                                     AFB, LDAFB, IWORK, WORK, LDB, &
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
                        AINVNM = CLANGE( '1', N, N, WORK, LDB, &
                                 RWORK )
                        IF( ANORMO <= ZERO .OR. AINVNM <= ZERO ) THEN
                           RCONDO = ONE
                        ELSE
                           RCONDO = ( ONE / ANORMO ) / AINVNM
                        END IF
!
!                          Compute the infinity-norm condition number
!                          of A.
!
                        AINVNM = CLANGE( 'I', N, N, WORK, LDB, &
                                 RWORK )
                        IF( ANORMI <= ZERO .OR. AINVNM <= ZERO ) THEN
                           RCONDI = ONE
                        ELSE
                           RCONDI = ( ONE / ANORMI ) / AINVNM
                        END IF
                     END IF
!
                     DO ITRAN = 1, NTRAN
!
!                          Do for each value of TRANS.
!
                        TRANS = TRANSS( ITRAN )
                        IF( ITRAN == 1 ) THEN
                           RCONDC = RCONDO
                        ELSE
                           RCONDC = RCONDI
                        END IF
!
!                          Restore the matrix A.
!
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CLACPY( 'Full', KL+KU+1, N, ASAV, LDA, &
                                     A, LDA )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
!
!                          Form an exact solution and set the right hand
!                          side.
!
                        SRNAMT = 'CLARHS'
                        CALL CLARHS( PATH, XTYPE, 'Full', TRANS, N, &
                                     N, KL, KU, NRHS, A, LDA, XACT, &
                                     LDB, B, LDB, ISEED, INFO )
                        XTYPE = 'C'
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CLACPY( 'Full', N, NRHS, B, LDB, BSAV, &
                                     LDB )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
!
                        IF( NOFACT .AND. ITRAN == 1 ) THEN
!
!                             --- Test CGBSV  ---
!
!                             Compute the LU factorization of the matrix
!                             and solve the system.
!
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                           CALL CLACPY( 'Full', KL+KU+1, N, A, LDA, &
                                        AFB( KL+1 ), LDAFB )
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S2_time)
                           open(file='results.out', unit=10, position = 'append')
                           write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                                 real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                           close(10)
#endif
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
                           SRNAMT = 'CGBSV '
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                           CALL CGBSV( N, KL, KU, NRHS, AFB, LDAFB, &
                                       IWORK, X, LDB, INFO )
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S2_time)
                           open(file='results.out', unit=10, position = 'append')
                           write(10,'(A,F16.10,A)') 'Total time : CGBSV : ',&
                                 real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                           close(10)
#endif
!
!                             Check error code from CGBSV .
!
                           IF( INFO /= IZERO ) &
                              CALL ALAERH( PATH, 'CGBSV ', INFO, &
                                           IZERO, ' ', N, N, KL, KU, &
                                           NRHS, IMAT, NFAIL, NERRS, &
                                           NOUT )
!
!                             Reconstruct matrix from factors and
!                             compute residual.
!
                           CALL CGBT01( N, N, KL, KU, A, LDA, AFB, &
                                        LDAFB, IWORK, WORK, &
                                        RESULT( 1 ) )
                           NT = 1
                           IF( IZERO == 0 ) THEN
!
!                                Compute residual of the computed
!                                solution.
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
                              CALL CGBT02( 'No transpose', N, N, KL, &
                                           KU, NRHS, A, LDA, X, LDB, &
                                           WORK, LDB, RWORK, &
                                           RESULT( 2 ) )
!
!                                Check solution from generated exact
!                                solution.
!
                              CALL CGET04( N, NRHS, X, LDB, XACT, &
                                           LDB, RCONDC, RESULT( 3 ) )
                              NT = 3
                           END IF
!
!                             Print information about the tests that did
!                             not pass the threshold.
!
                           DO K = 1, NT
                              IF( RESULT( K ) >= THRESH ) THEN
                                 IF( NFAIL == 0 .AND. NERRS == 0 ) &
                                    CALL ALADHD( NOUT, PATH )
                                 WRITE( NOUT, FMT = 9997 )'CGBSV ', &
                                    N, KL, KU, IMAT, K, RESULT( K )
                                 NFAIL = NFAIL + 1
                              END IF
                           ENDDO
                           NRUN = NRUN + NT
                        END IF
!
!                          --- Test CGBSVX ---
!
                        IF( .NOT.PREFAC )  THEN
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                           CALL CLASET( 'Full', 2*KL+KU+1, N, &
                                        CMPLX( ZERO ), CMPLX( ZERO ), &
                                        AFB, LDAFB )
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S2_time)
                           open(file='results.out', unit=10, position = 'append')
                           write(10,'(A,F16.10,A)') 'Total time : CLASET : ',&
                                 real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                           close(10)
#endif
                        ENDIF
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CLASET( 'Full', N, NRHS, CMPLX( ZERO ), &
                                     CMPLX( ZERO ), X, LDB )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CLASET : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                        IF( IEQUED > 1 .AND. N > 0 ) THEN
!
!                             Equilibrate the matrix if FACT = 'F' and
!                             EQUED = 'R', 'C', or 'B'.
!
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                           CALL CLAQGB( N, N, KL, KU, A, LDA, S, &
                                        S( N+1 ), ROWCND, COLCND, &
                                        AMAX, EQUED )
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S2_time)
                           open(file='results.out', unit=10, position = 'append')
                           write(10,'(A,F16.10,A)') 'Total time : CLAQGB : ',&
                                 real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                           close(10)
#endif
                        END IF
!
!                          Solve the system and compute the condition
!                          number and error bounds using CGBSVX.
!
                        SRNAMT = 'CGBSVX'
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CGBSVX( FACT, TRANS, N, KL, KU, NRHS, A, &
                                     LDA, AFB, LDAFB, IWORK, EQUED, &
                                     S, S( LDB+1 ), B, LDB, X, LDB, &
                                     RCOND, RWORK, RWORK( NRHS+1 ), &
                                     WORK, RWORK( 2*NRHS+1 ), INFO )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CGBSVX : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
!
!                          Check the error code from CGBSVX.
!
                        IF( INFO /= IZERO ) &
                           CALL ALAERH( PATH, 'CGBSVX', INFO, IZERO, &
                                        FACT // TRANS, N, N, KL, KU, &
                                        NRHS, IMAT, NFAIL, NERRS, &
                                        NOUT )
!
!                          Compare RWORK(2*NRHS+1) from CGBSVX with the
!                          computed reciprocal pivot growth RPVGRW
!
                        IF( INFO /= 0 ) THEN
                           ANRMPV = ZERO
                           DO J = 1, INFO
                              DO I = MAX( KU+2-J, 1 ), &
                                      MIN( N+KU+1-J, KL+KU+1 )
                                 ANRMPV = MAX( ANRMPV, &
                                          ABS( A( I+( J-1 )*LDA ) ) )
                              ENDDO
                           ENDDO
                           RPVGRW = CLANTB( 'M', 'U', 'N', INFO, &
                                    MIN( INFO-1, KL+KU ), &
                                    AFB( MAX( 1, KL+KU+2-INFO ) ), &
                                    LDAFB, RDUM )
                           IF( RPVGRW == ZERO ) THEN
                              RPVGRW = ONE
                           ELSE
                              RPVGRW = ANRMPV / RPVGRW
                           END IF
                        ELSE
                           RPVGRW = CLANTB( 'M', 'U', 'N', N, KL+KU, &
                                    AFB, LDAFB, RDUM )
                           IF( RPVGRW == ZERO ) THEN
                              RPVGRW = ONE
                           ELSE
                              RPVGRW = CLANGB( 'M', N, KL, KU, A, &
                                       LDA, RDUM ) / RPVGRW
                           END IF
                        END IF
                        RESULT( 7 ) = ABS( RPVGRW-RWORK( 2*NRHS+1 ) ) &
                                       / MAX( RWORK( 2*NRHS+1 ), &
                                      RPVGRW ) / SLAMCH( 'E' )
!
                        IF( .NOT.PREFAC ) THEN
!
!                             Reconstruct matrix from factors and
!                             compute residual.
!
                           CALL CGBT01( N, N, KL, KU, A, LDA, AFB, &
                                        LDAFB, IWORK, WORK, &
                                        RESULT( 1 ) )
                           K1 = 1
                        ELSE
                           K1 = 2
                        END IF
!
                        IF( INFO == 0 ) THEN
                           TRFCON = .FALSE.
!
!                             Compute residual of the computed solution.
!
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                           CALL CLACPY( 'Full', N, NRHS, BSAV, LDB, &
                                        WORK, LDB )
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S2_time)
                           open(file='results.out', unit=10, position = 'append')
                           write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                                 real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                           close(10)
#endif
                           CALL CGBT02( TRANS, N, N, KL, KU, NRHS, &
                                        ASAV, LDA, X, LDB, WORK, LDB, &
                                        RWORK( 2*NRHS+1 ), &
                                        RESULT( 2 ) )
!
!                             Check solution from generated exact
!                             solution.
!
                           IF( NOFACT .OR. ( PREFAC .AND. &
                               LSAME( EQUED, 'N' ) ) ) THEN
                              CALL CGET04( N, NRHS, X, LDB, XACT, &
                                           LDB, RCONDC, RESULT( 3 ) )
                           ELSE
                              IF( ITRAN == 1 ) THEN
                                 ROLDC = ROLDO
                              ELSE
                                 ROLDC = ROLDI
                              END IF
                              CALL CGET04( N, NRHS, X, LDB, XACT, &
                                           LDB, ROLDC, RESULT( 3 ) )
                           END IF
!
!                             Check the error bounds from iterative
!                             refinement.
!
                           CALL CGBT05( TRANS, N, KL, KU, NRHS, ASAV, &
                                        LDA, BSAV, LDB, X, LDB, XACT, &
                                        LDB, RWORK, RWORK( NRHS+1 ), &
                                        RESULT( 4 ) )
                        ELSE
                           TRFCON = .TRUE.
                        END IF
!
!                          Compare RCOND from CGBSVX with the computed
!                          value in RCONDC.
!
                        RESULT( 6 ) = SGET06( RCOND, RCONDC )
!
!                          Print information about the tests that did
!                          not pass the threshold.
!
                        IF( .NOT.TRFCON ) THEN
                           DO K = K1, NTESTS
                              IF( RESULT( K ) >= THRESH ) THEN
                                 IF( NFAIL == 0 .AND. NERRS == 0 ) &
                                    CALL ALADHD( NOUT, PATH )
                                 IF( PREFAC ) THEN
                                    WRITE( NOUT, FMT = 9995 ) &
                                       'CGBSVX', FACT, TRANS, N, KL, &
                                       KU, EQUED, IMAT, K, &
                                       RESULT( K )
                                 ELSE
                                    WRITE( NOUT, FMT = 9996 ) &
                                       'CGBSVX', FACT, TRANS, N, KL, &
                                       KU, IMAT, K, RESULT( K )
                                 END IF
                                 NFAIL = NFAIL + 1
                              END IF
                           ENDDO
                           NRUN = NRUN + 7 - K1
                        ELSE
                           IF( RESULT( 1 ) >= THRESH .AND. .NOT. &
                               PREFAC ) THEN
                              IF( NFAIL == 0 .AND. NERRS == 0 ) &
                                 CALL ALADHD( NOUT, PATH )
                              IF( PREFAC ) THEN
                                 WRITE( NOUT, FMT = 9995 )'CGBSVX', &
                                    FACT, TRANS, N, KL, KU, EQUED, &
                                    IMAT, 1, RESULT( 1 )
                              ELSE
                                 WRITE( NOUT, FMT = 9996 )'CGBSVX', &
                                    FACT, TRANS, N, KL, KU, IMAT, 1, &
                                    RESULT( 1 )
                              END IF
                              NFAIL = NFAIL + 1
                              NRUN = NRUN + 1
                           END IF
                           IF( RESULT( 6 ) >= THRESH ) THEN
                              IF( NFAIL == 0 .AND. NERRS == 0 ) &
                                 CALL ALADHD( NOUT, PATH )
                              IF( PREFAC ) THEN
                                 WRITE( NOUT, FMT = 9995 )'CGBSVX', &
                                    FACT, TRANS, N, KL, KU, EQUED, &
                                    IMAT, 6, RESULT( 6 )
                              ELSE
                                 WRITE( NOUT, FMT = 9996 )'CGBSVX', &
                                    FACT, TRANS, N, KL, KU, IMAT, 6, &
                                    RESULT( 6 )
                              END IF
                              NFAIL = NFAIL + 1
                              NRUN = NRUN + 1
                           END IF
                           IF( RESULT( 7 ) >= THRESH ) THEN
                              IF( NFAIL == 0 .AND. NERRS == 0 ) &
                                 CALL ALADHD( NOUT, PATH )
                              IF( PREFAC ) THEN
                                 WRITE( NOUT, FMT = 9995 )'CGBSVX', &
                                    FACT, TRANS, N, KL, KU, EQUED, &
                                    IMAT, 7, RESULT( 7 )
                              ELSE
                                 WRITE( NOUT, FMT = 9996 )'CGBSVX', &
                                    FACT, TRANS, N, KL, KU, IMAT, 7, &
                                    RESULT( 7 )
                              END IF
                              NFAIL = NFAIL + 1
                              NRUN = NRUN + 1
                           END IF
                        END IF

!                    --- Test CGBSVXX ---

!                    Restore the matrices A and B.

!                     write(*,*) 'begin cgbsvxx testing'

#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL CLACPY( 'Full', KL+KU+1, N, ASAV, LDA, A, &
                             LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL CLACPY( 'Full', N, NRHS, BSAV, LDB, B, LDB )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif

                  IF( .NOT.PREFAC )  THEN
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CLASET( 'Full', 2*KL+KU+1, N, &
                                  CMPLX( ZERO ), CMPLX( ZERO ), &
                                  AFB, LDAFB )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CLASET : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                  ENDIF
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL CLASET( 'Full', N, NRHS, &
                               CMPLX( ZERO ), CMPLX( ZERO ), &
                                  X, LDB )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CLASET : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  IF( IEQUED > 1 .AND. N > 0 ) THEN
!
!                       Equilibrate the matrix if FACT = 'F' and
!                       EQUED = 'R', 'C', or 'B'.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CLAQGB( N, N, KL, KU, A, LDA, S, &
                          S( N+1 ), ROWCND, COLCND, AMAX, EQUED )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CLAQGB : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                  END IF
!
!                    Solve the system and compute the condition number
!                    and error bounds using CGBSVXX.
!
                  SRNAMT = 'CGBSVXX'
                  N_ERR_BNDS = 3
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL CGBSVXX( FACT, TRANS, N, KL, KU, NRHS, A, LDA, &
                       AFB, LDAFB, IWORK, EQUED, S, S( N+1 ), B, LDB, &
                       X, LDB, RCOND, RPVGRW_SVXX, BERR, N_ERR_BNDS, &
                       ERRBNDS_N, ERRBNDS_C, 0, ZERO, WORK, &
                       RWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CGBSVXX : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check the error code from CGBSVXX.
!
                  IF( INFO == N+1 ) GOTO 90
                  IF( INFO /= IZERO ) THEN
                     CALL ALAERH( PATH, 'CGBSVXX', INFO, IZERO, &
                                  FACT // TRANS, N, N, -1, -1, NRHS, &
                                  IMAT, NFAIL, NERRS, NOUT )
                     GOTO 90
                  END IF
!
!                    Compare rpvgrw_svxx from CGESVXX with the computed
!                    reciprocal pivot growth factor RPVGRW
!

                  IF ( INFO  >  0 .AND. INFO  <  N+1 ) THEN
                     RPVGRW = CLA_GBRPVGRW(N, KL, KU, INFO, A, LDA, &
                          AFB, LDAFB)
                  ELSE
                     RPVGRW = CLA_GBRPVGRW(N, KL, KU, N, A, LDA, &
                          AFB, LDAFB)
                  ENDIF

                  RESULT( 7 ) = ABS( RPVGRW-rpvgrw_svxx ) / &
                                MAX( rpvgrw_svxx, RPVGRW ) / &
                                SLAMCH( 'E' )
!
                  IF( .NOT.PREFAC ) THEN
!
!                       Reconstruct matrix from factors and compute
!                       residual.
!
                     CALL CGBT01( N, N, KL, KU, A, LDA, AFB, LDAFB, &
                          IWORK, WORK( 2*NRHS+1 ), RESULT( 1 ) )
                     K1 = 1
                  ELSE
                     K1 = 2
                  END IF
!
                  IF( INFO == 0 ) THEN
                     TRFCON = .FALSE.
!
!                       Compute residual of the computed solution.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CLACPY( 'Full', N, NRHS, BSAV, LDB, WORK, &
                                  LDB )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     CALL CGBT02( TRANS, N, N, KL, KU, NRHS, ASAV, &
                                  LDA, X, LDB, WORK, LDB, RWORK, &
                                  RESULT( 2 ) )
!
!                       Check solution from generated exact solution.
!
                     IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED, &
                         'N' ) ) ) THEN
                        CALL CGET04( N, NRHS, X, LDB, XACT, LDB, &
                                     RCONDC, RESULT( 3 ) )
                     ELSE
                        IF( ITRAN == 1 ) THEN
                           ROLDC = ROLDO
                        ELSE
                           ROLDC = ROLDI
                        END IF
                        CALL CGET04( N, NRHS, X, LDB, XACT, LDB, &
                                     ROLDC, RESULT( 3 ) )
                     END IF
                  ELSE
                     TRFCON = .TRUE.
                  END IF
!
!                    Compare RCOND from CGBSVXX with the computed value
!                    in RCONDC.
!
                  RESULT( 6 ) = SGET06( RCOND, RCONDC )
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                  IF( .NOT.TRFCON ) THEN
                     DO K = K1, NTESTS
                        IF( RESULT( K ) >= THRESH ) THEN
                           IF( NFAIL == 0 .AND. NERRS == 0 ) &
                              CALL ALADHD( NOUT, PATH )
                           IF( PREFAC ) THEN
                              WRITE( NOUT, FMT = 9995 )'CGBSVXX', &
                                   FACT, TRANS, N, KL, KU, EQUED, &
                                   IMAT, K, RESULT( K )
                           ELSE
                              WRITE( NOUT, FMT = 9996 )'CGBSVXX', &
                                   FACT, TRANS, N, KL, KU, IMAT, K, &
                                   RESULT( K )
                           END IF
                           NFAIL = NFAIL + 1
                        END IF
                        ENDDO
                     NRUN = NRUN + 7 - K1
                  ELSE
                     IF( RESULT( 1 ) >= THRESH .AND. .NOT.PREFAC ) &
                          THEN
                        IF( NFAIL == 0 .AND. NERRS == 0 ) &
                           CALL ALADHD( NOUT, PATH )
                        IF( PREFAC ) THEN
                           WRITE( NOUT, FMT = 9995 )'CGBSVXX', FACT, &
                                TRANS, N, KL, KU, EQUED, IMAT, 1, &
                                RESULT( 1 )
                        ELSE
                           WRITE( NOUT, FMT = 9996 )'CGBSVXX', FACT, &
                                TRANS, N, KL, KU, IMAT, 1, &
                                RESULT( 1 )
                        END IF
                        NFAIL = NFAIL + 1
                        NRUN = NRUN + 1
                     END IF
                     IF( RESULT( 6 ) >= THRESH ) THEN
                        IF( NFAIL == 0 .AND. NERRS == 0 ) &
                           CALL ALADHD( NOUT, PATH )
                        IF( PREFAC ) THEN
                           WRITE( NOUT, FMT = 9995 )'CGBSVXX', FACT, &
                                TRANS, N, KL, KU, EQUED, IMAT, 6, &
                                RESULT( 6 )
                        ELSE
                           WRITE( NOUT, FMT = 9996 )'CGBSVXX', FACT, &
                                TRANS, N, KL, KU, IMAT, 6, &
                                RESULT( 6 )
                        END IF
                        NFAIL = NFAIL + 1
                        NRUN = NRUN + 1
                     END IF
                     IF( RESULT( 7 ) >= THRESH ) THEN
                        IF( NFAIL == 0 .AND. NERRS == 0 ) &
                           CALL ALADHD( NOUT, PATH )
                        IF( PREFAC ) THEN
                           WRITE( NOUT, FMT = 9995 )'CGBSVXX', FACT, &
                                TRANS, N, KL, KU, EQUED, IMAT, 7, &
                                RESULT( 7 )
                        ELSE
                           WRITE( NOUT, FMT = 9996 )'CGBSVXX', FACT, &
                                TRANS, N, KL, KU, IMAT, 7, &
                                RESULT( 7 )
                        END IF
                        NFAIL = NFAIL + 1
                        NRUN = NRUN + 1
                     END IF
!
                  END IF
!
90                   CONTINUE
                     ENDDO
  100                CONTINUE
                     ENDDO
                  ENDDO
  120          CONTINUE
               ENDDO
  130       CONTINUE
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )
!

!     Test Error Bounds from CGBSVXX

   CALL CEBCHVXX(THRESH, PATH)

 9999 FORMAT( ' *** In CDRVGB, LA=', I5, ' is too small for N=', I5, &
         ', KU=', I5, ', KL=', I5, / ' ==> Increase LA to at least ', &
         I5 )
 9998 FORMAT( ' *** In CDRVGB, LAFB=', I5, ' is too small for N=', I5, &
         ', KU=', I5, ', KL=', I5, / &
         ' ==> Increase LAFB to at least ', I5 )
 9997 FORMAT( 1X, A, ', N=', I5, ', KL=', I5, ', KU=', I5, ', type ', &
         I1, ', test(', I1, ')=', G12.5 )
 9996 FORMAT( 1X, A, '( ''', A1, ''',''', A1, ''',', I5, ',', I5, ',', &
         I5, ',...), type ', I1, ', test(', I1, ')=', G12.5 )
 9995 FORMAT( 1X, A, '( ''', A1, ''',''', A1, ''',', I5, ',', I5, ',', &
         I5, ',...), EQUED=''', A1, ''', type ', I1, ', test(', I1, &
         ')=', G12.5 )
!
   RETURN
!
!     End of CDRVGBX
!
END
                                                                                                                                                                                                                                                                                                            



