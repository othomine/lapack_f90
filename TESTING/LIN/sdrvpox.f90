!> \brief \b SDRVPOX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SDRVPO( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
!                          A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK,
!                          RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NOUT, NRHS
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NVAL( * )
!       REAL               A( * ), AFAC( * ), ASAV( * ), B( * ),
!      $                   BSAV( * ), RWORK( * ), S( * ), WORK( * ),
!      $                   X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SDRVPO tests the driver routines SPOSV, -SVX, and -SVXX.
!>
!> Note that this file is used only when the XBLAS are available,
!> otherwise sdrvpo.f defines this subroutine.
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
!> \param[out] ASAV
!> \verbatim
!>          ASAV is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is REAL array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] BSAV
!> \verbatim
!>          BSAV is REAL array, dimension (NMAX*NRHS)
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
!> \param[out] S
!> \verbatim
!>          S is REAL array, dimension (NMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension
!>                      (NMAX*max(3,NRHS))
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
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
!
!> \ingroup single_lin
!
!  =====================================================================
   SUBROUTINE SDRVPO( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, &
                      A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK, &
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
   REAL               A( * ), AFAC( * ), ASAV( * ), B( * ), &
                      BSAV( * ), RWORK( * ), S( * ), WORK( * ), &
                      X( * ), XACT( * )
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
                      NERRS, NFACT, NFAIL, NIMAT, NRUN, NT, &
                      N_ERR_BNDS
   REAL               AINVNM, AMAX, ANORM, CNDNUM, RCOND, RCONDC, &
                      ROLDC, SCOND, RPVGRW_SVXX
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          EQUEDS( 2 ), FACTS( 3 ), UPLOS( 2 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   REAL               RESULT( NTESTS ), BERR( NRHS ), &
                      ERRBNDS_N( NRHS, 3 ), ERRBNDS_C( NRHS, 3 )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               SGET06, SLANSY
   EXTERNAL           LSAME, SGET06, SLANSY
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALADHD, ALAERH, ALASVM, SERRVX, SGET04, SLACPY, &
                      SLAQSY, SLARHS, SLASET, SLATB4, SLATMS, SPOEQU, &
                      SPOSV, SPOSVX, SPOT01, SPOT02, SPOT05, SPOTRF, &
                      SPOTRI, XLAENV
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX
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
   DATA               FACTS / 'F', 'N', 'E' /
   DATA               EQUEDS / 'N', 'Y' /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   PATH( 1: 1 ) = 'Single precision'
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
!              Set up parameters with SLATB4 and generate a test matrix
!              with SLATMS.
!
            CALL SLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, &
                         CNDNUM, DIST )
!
            SRNAMT = 'SLATMS'
            CALL SLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, &
                         CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, &
                         INFO )
!
!              Check error code from SLATMS.
!
            IF( INFO /= 0 ) THEN
               CALL ALAERH( PATH, 'SLATMS', INFO, 0, UPLO, N, N, -1, &
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
!              Save a copy of the matrix A in ASAV.
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SLACPY( UPLO, N, N, A, LDA, ASAV, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
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
!                       the value returned by SPOSVX (FACT = 'N' reuses
!                       the condition number from the previous iteration
!                       with FACT = 'F').
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SLACPY( UPLO, N, N, ASAV, LDA, AFAC, LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
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
                        CALL SPOEQU( N, AFAC, LDA, S, SCOND, AMAX, &
                                     INFO )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : SPOEQU : ',&
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
                           CALL SLAQSY( UPLO, N, AFAC, LDA, S, SCOND, &
                                        AMAX, EQUED )
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S2_time)
                           open(file='results.out', unit=10, position = 'append')
                           write(10,'(A,F16.10,A)') 'Total time : SLAQSY : ',&
                                 real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                           close(10)
#endif
                        END IF
                     END IF
!
!                       Save the condition number of the
!                       non-equilibrated system for use in SGET04.
!
                     IF( EQUIL ) &
                        ROLDC = RCONDC
!
!                       Compute the 1-norm of A.
!
                     ANORM = SLANSY( '1', UPLO, N, AFAC, LDA, RWORK )
!
!                       Factor the matrix A.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SPOTRF( UPLO, N, AFAC, LDA, INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SPOTRF : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                       Form the inverse of A.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SLACPY( UPLO, N, N, AFAC, LDA, A, LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SPOTRI( UPLO, N, A, LDA, INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SPOTRI : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                       Compute the 1-norm condition number of A.
!
                     AINVNM = SLANSY( '1', UPLO, N, A, LDA, RWORK )
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
                  CALL SLACPY( UPLO, N, N, ASAV, LDA, A, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Form an exact solution and set the right hand side.
!
                  SRNAMT = 'SLARHS'
                  CALL SLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU, &
                               NRHS, A, LDA, XACT, LDA, B, LDA, &
                               ISEED, INFO )
                  XTYPE = 'C'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLACPY( 'Full', N, NRHS, B, LDA, BSAV, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
                  IF( NOFACT ) THEN
!
!                       --- Test SPOSV  ---
!
!                       Compute the L*L' or U'*U factorization of the
!                       matrix and solve the system.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
                     SRNAMT = 'SPOSV '
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SPOSV( UPLO, N, NRHS, AFAC, LDA, X, LDA, &
                                 INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SPOSV : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                       Check error code from SPOSV .
!
                     IF( INFO /= IZERO ) THEN
                        CALL ALAERH( PATH, 'SPOSV ', INFO, IZERO, &
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
                     CALL SPOT01( UPLO, N, A, LDA, AFAC, LDA, RWORK, &
                                  RESULT( 1 ) )
!
!                       Compute residual of the computed solution.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SLACPY( 'Full', N, NRHS, B, LDA, WORK, &
                                  LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     CALL SPOT02( UPLO, N, NRHS, A, LDA, X, LDA, &
                                  WORK, LDA, RWORK, RESULT( 2 ) )
!
!                       Check solution from generated exact solution.
!
                     CALL SGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
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
                           WRITE( NOUT, FMT = 9999 )'SPOSV ', UPLO, &
                              N, IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        END IF
                     ENDDO
                     NRUN = NRUN + NT
70                   CONTINUE
                  END IF
!
!                    --- Test SPOSVX ---
!
                  IF( .NOT.PREFAC )  THEN
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SLASET( UPLO, N, N, ZERO, ZERO, AFAC, LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                  ENDIF
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLASET( 'Full', N, NRHS, ZERO, ZERO, X, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
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
                     CALL SLAQSY( UPLO, N, A, LDA, S, SCOND, AMAX, &
                                  EQUED )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLAQSY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                  END IF
!
!                    Solve the system and compute the condition number
!                    and error bounds using SPOSVX.
!
                  SRNAMT = 'SPOSVX'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SPOSVX( FACT, UPLO, N, NRHS, A, LDA, AFAC, &
                               LDA, EQUED, S, B, LDA, X, LDA, RCOND, &
                               RWORK, RWORK( NRHS+1 ), WORK, IWORK, &
                               INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SPOSVX : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check the error code from SPOSVX.
!
                  IF( INFO /= IZERO ) THEN
                     CALL ALAERH( PATH, 'SPOSVX', INFO, IZERO, &
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
                        CALL SPOT01( UPLO, N, A, LDA, AFAC, LDA, &
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
                     CALL SLACPY( 'Full', N, NRHS, BSAV, LDA, WORK, &
                                  LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     CALL SPOT02( UPLO, N, NRHS, ASAV, LDA, X, LDA, &
                                  WORK, LDA, RWORK( 2*NRHS+1 ), &
                                  RESULT( 2 ) )
!
!                       Check solution from generated exact solution.
!
                     IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED, &
                         'N' ) ) ) THEN
                        CALL SGET04( N, NRHS, X, LDA, XACT, LDA, &
                                     RCONDC, RESULT( 3 ) )
                     ELSE
                        CALL SGET04( N, NRHS, X, LDA, XACT, LDA, &
                                     ROLDC, RESULT( 3 ) )
                     END IF
!
!                       Check the error bounds from iterative
!                       refinement.
!
                     CALL SPOT05( UPLO, N, NRHS, ASAV, LDA, B, LDA, &
                                  X, LDA, XACT, LDA, RWORK, &
                                  RWORK( NRHS+1 ), RESULT( 4 ) )
                  ELSE
                     K1 = 6
                  END IF
!
!                    Compare RCOND from SPOSVX with the computed value
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
                           WRITE( NOUT, FMT = 9997 )'SPOSVX', FACT, &
                              UPLO, N, EQUED, IMAT, K, RESULT( K )
                        ELSE
                           WRITE( NOUT, FMT = 9998 )'SPOSVX', FACT, &
                              UPLO, N, IMAT, K, RESULT( K )
                        END IF
                        NFAIL = NFAIL + 1
                     END IF
                  ENDDO
                  NRUN = NRUN + 7 - K1
!
!                    --- Test SPOSVXX ---
!
!                    Restore the matrices A and B.
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLACPY( 'Full', N, N, ASAV, LDA, A, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLACPY( 'Full', N, NRHS, BSAV, LDA, B, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif

                  IF( .NOT.PREFAC )  THEN
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SLASET( UPLO, N, N, ZERO, ZERO, AFAC, LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                  ENDIF
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLASET( 'Full', N, NRHS, ZERO, ZERO, X, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
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
                     CALL SLAQSY( UPLO, N, A, LDA, S, SCOND, AMAX, &
                                  EQUED )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLAQSY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                  END IF
!
!                    Solve the system and compute the condition number
!                    and error bounds using SPOSVXX.
!
                  SRNAMT = 'SPOSVXX'
                  N_ERR_BNDS = 3
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SPOSVXX( FACT, UPLO, N, NRHS, A, LDA, AFAC, &
                       LDA, EQUED, S, B, LDA, X, &
                       LDA, RCOND, RPVGRW_SVXX, BERR, N_ERR_BNDS, &
                       ERRBNDS_N, ERRBNDS_C, 0, ZERO, WORK, &
                       IWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SPOSVXX : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check the error code from SPOSVXX.
!
                  IF( INFO == N+1 ) GOTO 90
                  IF( INFO /= IZERO ) THEN
                     CALL ALAERH( PATH, 'SPOSVXX', INFO, IZERO, &
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
                        CALL SPOT01( UPLO, N, A, LDA, AFAC, LDA, &
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
                     CALL SLACPY( 'Full', N, NRHS, BSAV, LDA, WORK, &
                                  LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     CALL SPOT02( UPLO, N, NRHS, ASAV, LDA, X, LDA, &
                                  WORK, LDA, RWORK( 2*NRHS+1 ), &
                                  RESULT( 2 ) )
!
!                       Check solution from generated exact solution.
!
                     IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED, &
                         'N' ) ) ) THEN
                        CALL SGET04( N, NRHS, X, LDA, XACT, LDA, &
                                     RCONDC, RESULT( 3 ) )
                     ELSE
                        CALL SGET04( N, NRHS, X, LDA, XACT, LDA, &
                                     ROLDC, RESULT( 3 ) )
                     END IF
!
!                       Check the error bounds from iterative
!                       refinement.
!
                     CALL SPOT05( UPLO, N, NRHS, ASAV, LDA, B, LDA, &
                                  X, LDA, XACT, LDA, RWORK, &
                                  RWORK( NRHS+1 ), RESULT( 4 ) )
                  ELSE
                     K1 = 6
                  END IF
!
!                    Compare RCOND from SPOSVXX with the computed value
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
                           WRITE( NOUT, FMT = 9997 )'SPOSVXX', FACT, &
                              UPLO, N, EQUED, IMAT, K, RESULT( K )
                        ELSE
                           WRITE( NOUT, FMT = 9998 )'SPOSVXX', FACT, &
                              UPLO, N, IMAT, K, RESULT( K )
                        END IF
                        NFAIL = NFAIL + 1
                     END IF
                     ENDDO
                  NRUN = NRUN + 7 - K1
  90               CONTINUE
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

!     Test Error Bounds from SPOSVXX

   CALL SEBCHVXX(THRESH, PATH)

 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I1, &
         ', test(', I1, ')=', G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5, &
         ', type ', I1, ', test(', I1, ')=', G12.5 )
 9997 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5, &
         ', EQUED=''', A1, ''', type ', I1, ', test(', I1, ') =', &
         G12.5 )
   RETURN
!
!     End of SDRVPOX
!
END
                                                                                                                                                                                                                                                                                                            




