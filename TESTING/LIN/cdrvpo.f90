!> \brief \b CDRVPO
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CDRVPO( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
!                          A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK,
!                          RWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NOUT, NRHS
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            NVAL( * )
!       REAL               RWORK( * ), S( * )
!       COMPLEX            A( * ), AFAC( * ), ASAV( * ), B( * ),
!      $                   BSAV( * ), WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CDRVPO tests the driver routines CPOSV and -SVX.
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
!>          A is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is COMPLEX array, dimension (NMAX*NMAX)
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
!>          S is REAL array, dimension (NMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension
!>                      (NMAX*max(3,NRHS))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (NMAX+2*NRHS)
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
   SUBROUTINE CDRVPO( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, &
                      A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK, &
                      RWORK, NOUT )
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
   INTEGER            NVAL( * )
   REAL               RWORK( * ), S( * )
   COMPLEX            A( * ), AFAC( * ), ASAV( * ), B( * ), &
                      BSAV( * ), WORK( * ), X( * ), XACT( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ONE, ZERO
   PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
   INTEGER            NTYPES
   PARAMETER          ( NTYPES = 9 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 6 )
!     ..
!     .. Local Scalars ..
   LOGICAL            EQUIL, NOFACT, PREFAC, ZEROT
   CHARACTER          DIST, EQUED, FACT, TYPE, UPLO, XTYPE
   CHARACTER*3        PATH
   INTEGER            I, IEQUED, IFACT, IMAT, IN, INFO, IOFF, IUPLO, &
                      IZERO, K, K1, KL, KU, LDA, MODE, N, NB, NBMIN, &
                      NERRS, NFACT, NFAIL, NIMAT, NRUN, NT
   REAL               AINVNM, AMAX, ANORM, CNDNUM, RCOND, RCONDC, &
                      ROLDC, SCOND
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          EQUEDS( 2 ), FACTS( 3 ), UPLOS( 2 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   REAL               RESULT( NTESTS )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               CLANHE, SGET06
   EXTERNAL           LSAME, CLANHE, SGET06
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALADHD, ALAERH, ALASVM, CERRVX, CGET04, CLACPY, &
                      CLAIPD, CLAQHE, CLARHS, CLASET, CLATB4, CLATMS, &
                      CPOEQU, CPOSV, CPOSVX, CPOT01, CPOT02, CPOT05, &
                      CPOTRF, CPOTRI, XLAENV
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
   INTRINSIC          CMPLX, MAX
!     ..
!     .. Data statements ..
   DATA               ISEEDY / 1988, 1989, 1990, 1991 /
   DATA               UPLOS / 'U', 'L' /
   DATA               FACTS / 'F', 'N', 'E' /
   DATA               EQUEDS / 'N', 'Y' /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   PATH( 1: 1 ) = 'Complex precision'
   PATH( 2: 3 ) = 'PO'
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
            GO TO 120
!
!           Skip types 3, 4, or 5 if the matrix size is too small.
!
         ZEROT = IMAT >= 3 .AND. IMAT <= 5
         IF( ZEROT .AND. N < IMAT-2 ) &
            GO TO 120
!
!           Do first for UPLO = 'U', then for UPLO = 'L'
!
         DO IUPLO = 1, 2
            UPLO = UPLOS( IUPLO )
!
!              Set up parameters with CLATB4 and generate a test matrix
!              with CLATMS.
!
            CALL CLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, &
                         CNDNUM, DIST )
!
            SRNAMT = 'CLATMS'
            CALL CLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, &
                         CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, &
                         INFO )
!
!              Check error code from CLATMS.
!
            IF( INFO /= 0 ) THEN
               CALL ALAERH( PATH, 'CLATMS', INFO, 0, UPLO, N, N, -1, &
                            -1, -1, IMAT, NFAIL, NERRS, NOUT )
               GO TO 110
            END IF
!
!              For types 3-5, zero one row and column of the matrix to
!              test that INFO is returned correctly.
!
            IF( ZEROT ) THEN
               IF( IMAT == 3 ) THEN
                  IZERO = 1
               ELSE IF( IMAT == 4 ) THEN
                  IZERO = N
               ELSE
                  IZERO = N / 2 + 1
               END IF
               IOFF = ( IZERO-1 )*LDA
!
!                 Set row and column IZERO of A to 0.
!
               IF( IUPLO == 1 ) THEN
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
               IZERO = 0
            END IF
!
!              Set the imaginary part of the diagonals.
!
            CALL CLAIPD( N, A, LDA+1, 0 )
!
!              Save a copy of the matrix A in ASAV.
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL CLACPY( UPLO, N, N, A, LDA, ASAV, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
            DO IEQUED = 1, 2
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
                        GO TO 90
                     RCONDC = ZERO
!
                  ELSE IF( .NOT.LSAME( FACT, 'N' ) ) THEN
!
!                       Compute the condition number for comparison with
!                       the value returned by CPOSVX (FACT = 'N' reuses
!                       the condition number from the previous iteration
!                       with FACT = 'F').
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CLACPY( UPLO, N, N, ASAV, LDA, AFAC, LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     IF( EQUIL .OR. IEQUED > 1 ) THEN
!
!                          Compute row and column scale factors to
!                          equilibrate the matrix A.
!
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CPOEQU( N, AFAC, LDA, S, SCOND, AMAX, &
                                     INFO )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CPOEQU : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                        IF( INFO == 0 .AND. N > 0 ) THEN
                           IF( IEQUED > 1 ) &
                              SCOND = ZERO
!
!                             Equilibrate the matrix.
!
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                           CALL CLAQHE( UPLO, N, AFAC, LDA, S, SCOND, &
                                        AMAX, EQUED )
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S2_time)
                           open(file='results.out', unit=10, position = 'append')
                           write(10,'(A,F16.10,A)') 'Total time : CLAQHE : ',&
                                 real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                           close(10)
#endif
                        END IF
                     END IF
!
!                       Save the condition number of the
!                       non-equilibrated system for use in CGET04.
!
                     IF( EQUIL ) &
                        ROLDC = RCONDC
!
!                       Compute the 1-norm of A.
!
                     ANORM = CLANHE( '1', UPLO, N, AFAC, LDA, RWORK )
!
!                       Factor the matrix A.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CPOTRF( UPLO, N, AFAC, LDA, INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CPOTRF : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                       Form the inverse of A.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CLACPY( UPLO, N, N, AFAC, LDA, A, LDA )
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
                     CALL CPOTRI( UPLO, N, A, LDA, INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CPOTRI : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                       Compute the 1-norm condition number of A.
!
                     AINVNM = CLANHE( '1', UPLO, N, A, LDA, RWORK )
                     IF( ANORM <= ZERO .OR. AINVNM <= ZERO ) THEN
                        RCONDC = ONE
                     ELSE
                        RCONDC = ( ONE / ANORM ) / AINVNM
                     END IF
                  END IF
!
!                    Restore the matrix A.
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL CLACPY( UPLO, N, N, ASAV, LDA, A, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Form an exact solution and set the right hand side.
!
                  SRNAMT = 'CLARHS'
                  CALL CLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU, &
                               NRHS, A, LDA, XACT, LDA, B, LDA, &
                               ISEED, INFO )
                  XTYPE = 'C'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL CLACPY( 'Full', N, NRHS, B, LDA, BSAV, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
                  IF( NOFACT ) THEN
!
!                       --- Test CPOSV  ---
!
!                       Compute the L*L' or U'*U factorization of the
!                       matrix and solve the system.
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
                     SRNAMT = 'CPOSV '
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CPOSV( UPLO, N, NRHS, AFAC, LDA, X, LDA, &
                                 INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CPOSV : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                       Check error code from CPOSV .
!
                     IF( INFO /= IZERO ) THEN
                        CALL ALAERH( PATH, 'CPOSV ', INFO, IZERO, &
                                     UPLO, N, N, -1, -1, NRHS, IMAT, &
                                     NFAIL, NERRS, NOUT )
                        GO TO 70
                     ELSE IF( INFO /= 0 ) THEN
                        GO TO 70
                     END IF
!
!                       Reconstruct matrix from factors and compute
!                       residual.
!
                     CALL CPOT01( UPLO, N, A, LDA, AFAC, LDA, RWORK, &
                                  RESULT( 1 ) )
!
!                       Compute residual of the computed solution.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CLACPY( 'Full', N, NRHS, B, LDA, WORK, &
                                  LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     CALL CPOT02( UPLO, N, NRHS, A, LDA, X, LDA, &
                                  WORK, LDA, RWORK, RESULT( 2 ) )
!
!                       Check solution from generated exact solution.
!
                     CALL CGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                                  RESULT( 3 ) )
                     NT = 3
!
!                       Print information about the tests that did not
!                       pass the threshold.
!
                     DO K = 1, NT
                        IF( RESULT( K ) >= THRESH ) THEN
                           IF( NFAIL == 0 .AND. NERRS == 0 ) &
                              CALL ALADHD( NOUT, PATH )
                           WRITE( NOUT, FMT = 9999 )'CPOSV ', UPLO, &
                              N, IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        END IF
                     ENDDO
                     NRUN = NRUN + NT
70                   CONTINUE
                  END IF
!
!                    --- Test CPOSVX ---
!
                  IF( .NOT.PREFAC )  THEN
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CLASET( UPLO, N, N, CMPLX( ZERO ), &
                                  CMPLX( ZERO ), AFAC, LDA )
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
                               CMPLX( ZERO ), X, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CLASET : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  IF( IEQUED > 1 .AND. N > 0 ) THEN
!
!                       Equilibrate the matrix if FACT='F' and
!                       EQUED='Y'.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CLAQHE( UPLO, N, A, LDA, S, SCOND, AMAX, &
                                  EQUED )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CLAQHE : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                  END IF
!
!                    Solve the system and compute the condition number
!                    and error bounds using CPOSVX.
!
                  SRNAMT = 'CPOSVX'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL CPOSVX( FACT, UPLO, N, NRHS, A, LDA, AFAC, &
                               LDA, EQUED, S, B, LDA, X, LDA, RCOND, &
                               RWORK, RWORK( NRHS+1 ), WORK, &
                               RWORK( 2*NRHS+1 ), INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CPOSVX : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check the error code from CPOSVX.
!
                  IF( INFO /= IZERO ) THEN
                     CALL ALAERH( PATH, 'CPOSVX', INFO, IZERO, &
                                  FACT // UPLO, N, N, -1, -1, NRHS, &
                                  IMAT, NFAIL, NERRS, NOUT )
                     GO TO 90
                  END IF
!
                  IF( INFO == 0 ) THEN
                     IF( .NOT.PREFAC ) THEN
!
!                          Reconstruct matrix from factors and compute
!                          residual.
!
                        CALL CPOT01( UPLO, N, A, LDA, AFAC, LDA, &
                                     RWORK( 2*NRHS+1 ), RESULT( 1 ) )
                        K1 = 1
                     ELSE
                        K1 = 2
                     END IF
!
!                       Compute residual of the computed solution.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CLACPY( 'Full', N, NRHS, BSAV, LDA, WORK, &
                                  LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     CALL CPOT02( UPLO, N, NRHS, ASAV, LDA, X, LDA, &
                                  WORK, LDA, RWORK( 2*NRHS+1 ), &
                                  RESULT( 2 ) )
!
!                       Check solution from generated exact solution.
!
                     IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED, &
                         'N' ) ) ) THEN
                        CALL CGET04( N, NRHS, X, LDA, XACT, LDA, &
                                     RCONDC, RESULT( 3 ) )
                     ELSE
                        CALL CGET04( N, NRHS, X, LDA, XACT, LDA, &
                                     ROLDC, RESULT( 3 ) )
                     END IF
!
!                       Check the error bounds from iterative
!                       refinement.
!
                     CALL CPOT05( UPLO, N, NRHS, ASAV, LDA, B, LDA, &
                                  X, LDA, XACT, LDA, RWORK, &
                                  RWORK( NRHS+1 ), RESULT( 4 ) )
                  ELSE
                     K1 = 6
                  END IF
!
!                    Compare RCOND from CPOSVX with the computed value
!                    in RCONDC.
!
                  RESULT( 6 ) = SGET06( RCOND, RCONDC )
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                  DO K = K1, 6
                     IF( RESULT( K ) >= THRESH ) THEN
                        IF( NFAIL == 0 .AND. NERRS == 0 ) &
                           CALL ALADHD( NOUT, PATH )
                        IF( PREFAC ) THEN
                           WRITE( NOUT, FMT = 9997 )'CPOSVX', FACT, &
                              UPLO, N, EQUED, IMAT, K, RESULT( K )
                        ELSE
                           WRITE( NOUT, FMT = 9998 )'CPOSVX', FACT, &
                              UPLO, N, IMAT, K, RESULT( K )
                        END IF
                        NFAIL = NFAIL + 1
                     END IF
                  ENDDO
                  NRUN = NRUN + 7 - K1
90             CONTINUE
               ENDDO
               ENDDO
  110       CONTINUE
            ENDDO
  120    CONTINUE
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I1, &
         ', test(', I1, ')=', G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5, &
         ', type ', I1, ', test(', I1, ')=', G12.5 )
 9997 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5, &
         ', EQUED=''', A1, ''', type ', I1, ', test(', I1, ') =', &
         G12.5 )
   RETURN
!
!     End of CDRVPO
!
END
                                                                                                                                                                                                                                                                                                            




