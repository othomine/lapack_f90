*> \brief \b CDRVPP
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CDRVPP( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
*                          A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK,
*                          RWORK, NOUT )
*
*       .. Scalar Arguments ..
*       LOGICAL            TSTERR
*       INTEGER            NMAX, NN, NOUT, NRHS
*       REAL               THRESH
*       ..
*       .. Array Arguments ..
*       LOGICAL            DOTYPE( * )
*       INTEGER            NVAL( * )
*       REAL               RWORK( * ), S( * )
*       COMPLEX            A( * ), AFAC( * ), ASAV( * ), B( * ),
*      $                   BSAV( * ), WORK( * ), X( * ), XACT( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CDRVPP tests the driver routines CPPSV and -SVX.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] DOTYPE
*> \verbatim
*>          DOTYPE is LOGICAL array, dimension (NTYPES)
*>          The matrix types to be used for testing.  Matrices of type j
*>          (for 1 <= j <= NTYPES) are used for testing if DOTYPE(j) =
*>          .TRUE.; if DOTYPE(j) = .FALSE., then type j is not used.
*> \endverbatim
*>
*> \param[in] NN
*> \verbatim
*>          NN is INTEGER
*>          The number of values of N contained in the vector NVAL.
*> \endverbatim
*>
*> \param[in] NVAL
*> \verbatim
*>          NVAL is INTEGER array, dimension (NN)
*>          The values of the matrix dimension N.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand side vectors to be generated for
*>          each linear system.
*> \endverbatim
*>
*> \param[in] THRESH
*> \verbatim
*>          THRESH is REAL
*>          The threshold value for the test ratios.  A result is
*>          included in the output file if RESULT >= THRESH.  To have
*>          every test ratio printed, use THRESH = 0.
*> \endverbatim
*>
*> \param[in] TSTERR
*> \verbatim
*>          TSTERR is LOGICAL
*>          Flag that indicates whether error exits are to be tested.
*> \endverbatim
*>
*> \param[in] NMAX
*> \verbatim
*>          NMAX is INTEGER
*>          The maximum value permitted for N, used in dimensioning the
*>          work arrays.
*> \endverbatim
*>
*> \param[out] A
*> \verbatim
*>          A is COMPLEX array, dimension (NMAX*(NMAX+1)/2)
*> \endverbatim
*>
*> \param[out] AFAC
*> \verbatim
*>          AFAC is COMPLEX array, dimension (NMAX*(NMAX+1)/2)
*> \endverbatim
*>
*> \param[out] ASAV
*> \verbatim
*>          ASAV is COMPLEX array, dimension (NMAX*(NMAX+1)/2)
*> \endverbatim
*>
*> \param[out] B
*> \verbatim
*>          B is COMPLEX array, dimension (NMAX*NRHS)
*> \endverbatim
*>
*> \param[out] BSAV
*> \verbatim
*>          BSAV is COMPLEX array, dimension (NMAX*NRHS)
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is COMPLEX array, dimension (NMAX*NRHS)
*> \endverbatim
*>
*> \param[out] XACT
*> \verbatim
*>          XACT is COMPLEX array, dimension (NMAX*NRHS)
*> \endverbatim
*>
*> \param[out] S
*> \verbatim
*>          S is REAL array, dimension (NMAX)
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension
*>                      (NMAX*max(3,NRHS))
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension (NMAX+2*NRHS)
*> \endverbatim
*>
*> \param[in] NOUT
*> \verbatim
*>          NOUT is INTEGER
*>          The unit number for output.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complex_lin
*
*  =====================================================================
      SUBROUTINE CDRVPP( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
     $                   A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK,
     $                   RWORK, NOUT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      INTEGER            NMAX, NN, NOUT, NRHS
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            NVAL( * )
      REAL               RWORK( * ), S( * )
      COMPLEX            A( * ), AFAC( * ), ASAV( * ), B( * ),
     $                   BSAV( * ), WORK( * ), X( * ), XACT( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      INTEGER            NTYPES
      PARAMETER          ( NTYPES = 9 )
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 6 )
*     ..
*     .. Local Scalars ..
      LOGICAL            EQUIL, NOFACT, PREFAC, ZEROT
      CHARACTER          DIST, EQUED, FACT, PACKIT, TYPE, UPLO, XTYPE
      CHARACTER*3        PATH
      INTEGER            I, IEQUED, IFACT, IMAT, IN, INFO, IOFF, IUPLO,
     $                   IZERO, K, K1, KL, KU, LDA, MODE, N, NERRS,
     $                   NFACT, NFAIL, NIMAT, NPP, NRUN, NT
      REAL               AINVNM, AMAX, ANORM, CNDNUM, RCOND, RCONDC,
     $                   ROLDC, SCOND
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
*     ..
*     .. Local Arrays ..
      CHARACTER          EQUEDS( 2 ), FACTS( 3 ), PACKS( 2 ), UPLOS( 2 )
      INTEGER            ISEED( 4 ), ISEEDY( 4 )
      REAL               RESULT( NTESTS )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               CLANHP, SGET06
      EXTERNAL           LSAME, CLANHP, SGET06
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALADHD, ALAERH, ALASVM, CCOPY, CERRVX, CGET04,
     $                   CLACPY, CLAIPD, CLAQHP, CLARHS, CLASET, CLATB4,
     $                   CLATMS, CPPEQU, CPPSV, CPPSVX, CPPT01, CPPT02,
     $                   CPPT05, CPPTRF, CPPTRI
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
      INTEGER            INFOT, NUNIT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MAX
*     ..
*     .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' / , FACTS / 'F', 'N', 'E' / ,
     $                   PACKS / 'C', 'R' / , EQUEDS / 'N', 'Y' /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'PP'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
*
*     Test the error exits
*
      IF( TSTERR )
     $   CALL CERRVX( PATH, NOUT )
      INFOT = 0
*
*     Do for each value of N in NVAL
*
      DO 140 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         NPP = N*( N+1 ) / 2
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 )
     $      NIMAT = 1
*
         DO 130 IMAT = 1, NIMAT
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) )
     $         GO TO 130
*
*           Skip types 3, 4, or 5 if the matrix size is too small.
*
            ZEROT = IMAT.GE.3 .AND. IMAT.LE.5
            IF( ZEROT .AND. N.LT.IMAT-2 )
     $         GO TO 130
*
*           Do first for UPLO = 'U', then for UPLO = 'L'
*
            DO 120 IUPLO = 1, 2
               UPLO = UPLOS( IUPLO )
               PACKIT = PACKS( IUPLO )
*
*              Set up parameters with CLATB4 and generate a test matrix
*              with CLATMS.
*
               CALL CLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE,
     $                      CNDNUM, DIST )
               RCONDC = ONE / CNDNUM
*
               SRNAMT = 'CLATMS'
               CALL CLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE,
     $                      CNDNUM, ANORM, KL, KU, PACKIT, A, LDA, WORK,
     $                      INFO )
*
*              Check error code from CLATMS.
*
               IF( INFO.NE.0 ) THEN
                  CALL ALAERH( PATH, 'CLATMS', INFO, 0, UPLO, N, N, -1,
     $                         -1, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 120
               END IF
*
*              For types 3-5, zero one row and column of the matrix to
*              test that INFO is returned correctly.
*
               IF( ZEROT ) THEN
                  IF( IMAT.EQ.3 ) THEN
                     IZERO = 1
                  ELSE IF( IMAT.EQ.4 ) THEN
                     IZERO = N
                  ELSE
                     IZERO = N / 2 + 1
                  END IF
*
*                 Set row and column IZERO of A to 0.
*
                  IF( IUPLO.EQ.1 ) THEN
                     IOFF = ( IZERO-1 )*IZERO / 2
                     DO 20 I = 1, IZERO - 1
                        A( IOFF+I ) = ZERO
   20                CONTINUE
                     IOFF = IOFF + IZERO
                     DO 30 I = IZERO, N
                        A( IOFF ) = ZERO
                        IOFF = IOFF + I
   30                CONTINUE
                  ELSE
                     IOFF = IZERO
                     DO 40 I = 1, IZERO - 1
                        A( IOFF ) = ZERO
                        IOFF = IOFF + N - I
   40                CONTINUE
                     IOFF = IOFF - IZERO
                     DO 50 I = IZERO, N
                        A( IOFF+I ) = ZERO
   50                CONTINUE
                  END IF
               ELSE
                  IZERO = 0
               END IF
*
*              Set the imaginary part of the diagonals.
*
               IF( IUPLO.EQ.1 ) THEN
                  CALL CLAIPD( N, A, 2, 1 )
               ELSE
                  CALL CLAIPD( N, A, N, -1 )
               END IF
*
*              Save a copy of the matrix A in ASAV.
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
               CALL CCOPY( NPP, A, 1, ASAV, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CCOPY : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
*
               DO 110 IEQUED = 1, 2
                  EQUED = EQUEDS( IEQUED )
                  IF( IEQUED.EQ.1 ) THEN
                     NFACT = 3
                  ELSE
                     NFACT = 1
                  END IF
*
                  DO 100 IFACT = 1, NFACT
                     FACT = FACTS( IFACT )
                     PREFAC = LSAME( FACT, 'F' )
                     NOFACT = LSAME( FACT, 'N' )
                     EQUIL = LSAME( FACT, 'E' )
*
                     IF( ZEROT ) THEN
                        IF( PREFAC )
     $                     GO TO 100
                        RCONDC = ZERO
*
                     ELSE IF( .NOT.LSAME( FACT, 'N' ) ) THEN
*
*                       Compute the condition number for comparison with
*                       the value returned by CPPSVX (FACT = 'N' reuses
*                       the condition number from the previous iteration
*                          with FACT = 'F').
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                        CALL CCOPY( NPP, ASAV, 1, AFAC, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CCOPY : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
                        IF( EQUIL .OR. IEQUED.GT.1 ) THEN
*
*                          Compute row and column scale factors to
*                          equilibrate the matrix A.
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                           CALL CPPEQU( UPLO, N, AFAC, S, SCOND, AMAX,
     $                                  INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CPPEQU : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
                           IF( INFO.EQ.0 .AND. N.GT.0 ) THEN
                              IF( IEQUED.GT.1 )
     $                           SCOND = ZERO
*
*                             Equilibrate the matrix.
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                              CALL CLAQHP( UPLO, N, AFAC, S, SCOND,
     $                                     AMAX, EQUED )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CLAQHP : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
                           END IF
                        END IF
*
*                       Save the condition number of the
*                       non-equilibrated system for use in CGET04.
*
                        IF( EQUIL )
     $                     ROLDC = RCONDC
*
*                       Compute the 1-norm of A.
*
                        ANORM = CLANHP( '1', UPLO, N, AFAC, RWORK )
*
*                       Factor the matrix A.
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                        CALL CPPTRF( UPLO, N, AFAC, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CPPTRF : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
*
*                       Form the inverse of A.
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                        CALL CCOPY( NPP, AFAC, 1, A, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CCOPY : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                        CALL CPPTRI( UPLO, N, A, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CPPTRI : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
*
*                       Compute the 1-norm condition number of A.
*
                        AINVNM = CLANHP( '1', UPLO, N, A, RWORK )
                        IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
                           RCONDC = ONE
                        ELSE
                           RCONDC = ( ONE / ANORM ) / AINVNM
                        END IF
                     END IF
*
*                    Restore the matrix A.
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                     CALL CCOPY( NPP, ASAV, 1, A, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CCOPY : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
*
*                    Form an exact solution and set the right hand side.
*
                     SRNAMT = 'CLARHS'
                     CALL CLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU,
     $                            NRHS, A, LDA, XACT, LDA, B, LDA,
     $                            ISEED, INFO )
                     XTYPE = 'C'
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                     CALL CLACPY( 'Full', N, NRHS, B, LDA, BSAV, LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
*
                     IF( NOFACT ) THEN
*
*                       --- Test CPPSV  ---
*
*                       Compute the L*L' or U'*U factorization of the
*                       matrix and solve the system.
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                        CALL CCOPY( NPP, A, 1, AFAC, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CCOPY : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                        CALL CLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
*
                        SRNAMT = 'CPPSV '
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                        CALL CPPSV( UPLO, N, NRHS, AFAC, X, LDA, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CPPSV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
*
*                       Check error code from CPPSV .
*
                        IF( INFO.NE.IZERO ) THEN
                           CALL ALAERH( PATH, 'CPPSV ', INFO, IZERO,
     $                                  UPLO, N, N, -1, -1, NRHS, IMAT,
     $                                  NFAIL, NERRS, NOUT )
                           GO TO 70
                        ELSE IF( INFO.NE.0 ) THEN
                           GO TO 70
                        END IF
*
*                       Reconstruct matrix from factors and compute
*                       residual.
*
                        CALL CPPT01( UPLO, N, A, AFAC, RWORK,
     $                               RESULT( 1 ) )
*
*                       Compute residual of the computed solution.
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                        CALL CLACPY( 'Full', N, NRHS, B, LDA, WORK,
     $                               LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
                        CALL CPPT02( UPLO, N, NRHS, A, X, LDA, WORK,
     $                               LDA, RWORK, RESULT( 2 ) )
*
*                       Check solution from generated exact solution.
*
                        CALL CGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC,
     $                               RESULT( 3 ) )
                        NT = 3
*
*                       Print information about the tests that did not
*                       pass the threshold.
*
                        DO 60 K = 1, NT
                           IF( RESULT( K ).GE.THRESH ) THEN
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                           CALL ALADHD( NOUT, PATH )
                              WRITE( NOUT, FMT = 9999 )'CPPSV ', UPLO,
     $                           N, IMAT, K, RESULT( K )
                              NFAIL = NFAIL + 1
                           END IF
   60                   CONTINUE
                        NRUN = NRUN + NT
   70                   CONTINUE
                     END IF
*
*                    --- Test CPPSVX ---
*
                     IF( .NOT.PREFAC .AND. NPP.GT.0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                        CALL CLASET( 'Full', NPP, 1, CMPLX( ZERO ),
     $                               CMPLX( ZERO ), AFAC, NPP )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CLASET : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
       ENDIF
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                     CALL CLASET( 'Full', N, NRHS, CMPLX( ZERO ),
     $                            CMPLX( ZERO ), X, LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CLASET : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
                     IF( IEQUED.GT.1 .AND. N.GT.0 ) THEN
*
*                       Equilibrate the matrix if FACT='F' and
*                       EQUED='Y'.
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                        CALL CLAQHP( UPLO, N, A, S, SCOND, AMAX, EQUED )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CLAQHP : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
                     END IF
*
*                    Solve the system and compute the condition number
*                    and error bounds using CPPSVX.
*
                     SRNAMT = 'CPPSVX'
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                     CALL CPPSVX( FACT, UPLO, N, NRHS, A, AFAC, EQUED,
     $                            S, B, LDA, X, LDA, RCOND, RWORK,
     $                            RWORK( NRHS+1 ), WORK,
     $                            RWORK( 2*NRHS+1 ), INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CPPSVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
*
*                    Check the error code from CPPSVX.
*
                     IF( INFO.NE.IZERO ) THEN
                        CALL ALAERH( PATH, 'CPPSVX', INFO, IZERO,
     $                               FACT // UPLO, N, N, -1, -1, NRHS,
     $                               IMAT, NFAIL, NERRS, NOUT )
                        GO TO 90
                     END IF
*
                     IF( INFO.EQ.0 ) THEN
                        IF( .NOT.PREFAC ) THEN
*
*                          Reconstruct matrix from factors and compute
*                          residual.
*
                           CALL CPPT01( UPLO, N, A, AFAC,
     $                                  RWORK( 2*NRHS+1 ), RESULT( 1 ) )
                           K1 = 1
                        ELSE
                           K1 = 2
                        END IF
*
*                       Compute residual of the computed solution.
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
                        CALL CLACPY( 'Full', N, NRHS, BSAV, LDA, WORK,
     $                               LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
                        CALL CPPT02( UPLO, N, NRHS, ASAV, X, LDA, WORK,
     $                               LDA, RWORK( 2*NRHS+1 ),
     $                               RESULT( 2 ) )
*
*                       Check solution from generated exact solution.
*
                        IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED,
     $                      'N' ) ) ) THEN
                           CALL CGET04( N, NRHS, X, LDA, XACT, LDA,
     $                                  RCONDC, RESULT( 3 ) )
                        ELSE
                           CALL CGET04( N, NRHS, X, LDA, XACT, LDA,
     $                                  ROLDC, RESULT( 3 ) )
                        END IF
*
*                       Check the error bounds from iterative
*                       refinement.
*
                        CALL CPPT05( UPLO, N, NRHS, ASAV, B, LDA, X,
     $                               LDA, XACT, LDA, RWORK,
     $                               RWORK( NRHS+1 ), RESULT( 4 ) )
                     ELSE
                        K1 = 6
                     END IF
*
*                    Compare RCOND from CPPSVX with the computed value
*                    in RCONDC.
*
                     RESULT( 6 ) = SGET06( RCOND, RCONDC )
*
*                    Print information about the tests that did not pass
*                    the threshold.
*
                     DO 80 K = K1, 6
                        IF( RESULT( K ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                        CALL ALADHD( NOUT, PATH )
                           IF( PREFAC ) THEN
                              WRITE( NOUT, FMT = 9997 )'CPPSVX', FACT,
     $                           UPLO, N, EQUED, IMAT, K, RESULT( K )
                           ELSE
                              WRITE( NOUT, FMT = 9998 )'CPPSVX', FACT,
     $                           UPLO, N, IMAT, K, RESULT( K )
                           END IF
                           NFAIL = NFAIL + 1
                        END IF
   80                CONTINUE
                     NRUN = NRUN + 7 - K1
   90                CONTINUE
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I1,
     $      ', test(', I1, ')=', G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5,
     $      ', type ', I1, ', test(', I1, ')=', G12.5 )
 9997 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5,
     $      ', EQUED=''', A1, ''', type ', I1, ', test(', I1, ')=',
     $      G12.5 )
      RETURN
*
*     End of CDRVPP
*
      END

