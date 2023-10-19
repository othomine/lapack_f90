!> \brief \b DDRVAC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DDRVAC( DOTYPE, NM, MVAL, NNS, NSVAL, THRESH, NMAX,
!                          A, AFAC, B, X, WORK,
!                          RWORK, SWORK, NOUT )
!
!       .. Scalar Arguments ..
!       INTEGER            NMAX, NM, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            MVAL( * ), NSVAL( * )
!       REAL               SWORK(*)
!       DOUBLE PRECISION   A( * ), AFAC( * ), B( * ),
!      $                   RWORK( * ), WORK( * ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DDRVAC tests DSPOSV.
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
!>          The number of values of N contained in the vector MVAL.
!> \endverbatim
!>
!> \param[in] MVAL
!> \verbatim
!>          MVAL is INTEGER array, dimension (NM)
!>          The values of the matrix dimension N.
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
!>          THRESH is DOUBLE PRECISION
!>          The threshold value for the test ratios.  A result is
!>          included in the output file if RESULT >= THRESH.  To have
!>          every test ratio printed, use THRESH = 0.
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
!>          A is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension
!>                      (NMAX*max(3,NSMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension
!>                      (max(2*NMAX,2*NSMAX+NWORK))
!> \endverbatim
!>
!> \param[out] SWORK
!> \verbatim
!>          SWORK is REAL array, dimension
!>                      (NMAX*(NSMAX+NMAX))
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
!> \ingroup double_lin
!
!  =====================================================================
   SUBROUTINE DDRVAC( DOTYPE, NM, MVAL, NNS, NSVAL, THRESH, NMAX, &
                      A, AFAC, B, X, WORK, &
                      RWORK, SWORK, NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            NMAX, NM, NNS, NOUT
   DOUBLE PRECISION   THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            MVAL( * ), NSVAL( * )
   REAL               SWORK(*)
   DOUBLE PRECISION   A( * ), AFAC( * ), B( * ), &
                      RWORK( * ), WORK( * ), X( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO
   PARAMETER          ( ZERO = 0.0D+0 )
   INTEGER            NTYPES
   PARAMETER          ( NTYPES = 9 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 1 )
!     ..
!     .. Local Scalars ..
   LOGICAL            ZEROT
   CHARACTER          DIST, TYPE, UPLO, XTYPE
   CHARACTER*3        PATH
   INTEGER            I, IM, IMAT, INFO, IOFF, IRHS, IUPLO, &
                      IZERO, KL, KU, LDA, MODE, N, &
                      NERRS, NFAIL, NIMAT, NRHS, NRUN
   DOUBLE PRECISION   ANORM, CNDNUM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          UPLOS( 2 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   DOUBLE PRECISION   RESULT( NTESTS )
!     ..
!     .. Local Variables ..
   INTEGER            ITER, KASE
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAERH, DLACPY, &
                      DLARHS, DLASET, DLATB4, DLATMS, &
                      DPOT06, DSPOSV
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          DBLE, MAX, SQRT
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
   KASE = 0
   PATH( 1: 1 ) = 'Double precision'
   PATH( 2: 3 ) = 'PO'
   NRUN = 0
   NFAIL = 0
   NERRS = 0
   DO I = 1, 4
      ISEED( I ) = ISEEDY( I )
   ENDDO
!
   INFOT = 0
!
!     Do for each value of N in MVAL
!
   DO IM = 1, NM
      N = MVAL( IM )
      LDA = MAX( N, 1 )
      NIMAT = NTYPES
      IF( N <= 0 ) &
         NIMAT = 1
!
      DO IMAT = 1, NIMAT
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
         IF( .NOT.DOTYPE( IMAT ) ) &
            GO TO 110
!
!           Skip types 3, 4, or 5 if the matrix size is too small.
!
         ZEROT = IMAT >= 3 .AND. IMAT <= 5
         IF( ZEROT .AND. N < IMAT-2 ) &
            GO TO 110
!
!           Do first for UPLO = 'U', then for UPLO = 'L'
!
         DO IUPLO = 1, 2
            UPLO = UPLOS( IUPLO )
!
!              Set up parameters with DLATB4 and generate a test matrix
!              with DLATMS.
!
            CALL DLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, &
                         CNDNUM, DIST )
!
            SRNAMT = 'DLATMS'
            CALL DLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, &
                         CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, &
                         INFO )
!
!              Check error code from DLATMS.
!
            IF( INFO /= 0 ) THEN
               CALL ALAERH( PATH, 'DLATMS', INFO, 0, UPLO, N, N, -1, &
                            -1, -1, IMAT, NFAIL, NERRS, NOUT )
               GO TO 100
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
            DO IRHS = 1, NNS
               NRHS = NSVAL( IRHS )
               XTYPE = 'N'
!
!                 Form an exact solution and set the right hand side.
!
               SRNAMT = 'DLARHS'
               CALL DLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU, &
                            NRHS, A, LDA, X, LDA, B, LDA, &
                            ISEED, INFO )
!
!                 Compute the L*L' or U'*U factorization of the
!                 matrix and solve the system.
!
               SRNAMT = 'DSPOSV '
               KASE = KASE + 1
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DLACPY( 'All', N, N, A, LDA, AFAC, LDA)
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DSPOSV( UPLO, N, NRHS, AFAC, LDA, B, LDA, X, LDA, &
                            WORK, SWORK, ITER, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DSPOSV : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif

               IF (ITER < 0) THEN
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL DLACPY( 'All', N, N, A, LDA, AFAC, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
               ENDIF
!
!                 Check error code from DSPOSV .
!
               IF( INFO /= IZERO ) THEN
!
                  IF( NFAIL == 0 .AND. NERRS == 0 ) &
                     CALL ALAHD( NOUT, PATH )
                  NERRS = NERRS + 1
!
                  IF( INFO /= IZERO .AND. IZERO /= 0 ) THEN
                     WRITE( NOUT, FMT = 9988 )'DSPOSV',INFO,IZERO,N, &
                        IMAT
                  ELSE
                     WRITE( NOUT, FMT = 9975 )'DSPOSV',INFO,N,IMAT
                  END IF
               END IF
!
!                 Skip the remaining test if the matrix is singular.
!
               IF( INFO /= 0 ) &
                  GO TO 110
!
!                 Check the quality of the solution
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DLACPY( 'All', N, NRHS, B, LDA, WORK, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
               CALL DPOT06( UPLO, N, NRHS, A, LDA, X, LDA, WORK, &
                  LDA, RWORK, RESULT( 1 ) )
!
!                 Check if the test passes the testing.
!                 Print information about the tests that did not
!                 pass the testing.
!
!                 If iterative refinement has been used and claimed to
!                 be successful (ITER>0), we want
!                 NORM1(B - A*X)/(NORM1(A)*NORM1(X)*EPS*SRQT(N)) < 1
!
!                 If double precision has been used (ITER<0), we want
!                 NORM1(B - A*X)/(NORM1(A)*NORM1(X)*EPS) < THRES
!                 (Cf. the linear solver testing routines)
!
               IF ((THRESH <= 0.0E+00) &
                  .OR.((ITER >= 0).AND.(N > 0) &
                  .AND.(RESULT(1) >= SQRT(DBLE(N)))) &
                  .OR.((ITER < 0).AND.(RESULT(1) >= THRESH))) THEN
!
                  IF( NFAIL == 0 .AND. NERRS == 0 ) THEN
                     WRITE( NOUT, FMT = 8999 )'DPO'
                     WRITE( NOUT, FMT = '( '' Matrix types:'' )' )
                     WRITE( NOUT, FMT = 8979 )
                     WRITE( NOUT, FMT = '( '' Test ratios:'' )' )
                     WRITE( NOUT, FMT = 8960 )1
                     WRITE( NOUT, FMT = '( '' Messages:'' )' )
                  END IF
!
                  WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, IMAT, 1, &
                     RESULT( 1 )
!
                  NFAIL = NFAIL + 1
!
               END IF
!
               NRUN = NRUN + 1
!
            ENDDO
  100       CONTINUE
            ENDDO
  110    CONTINUE
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   IF( NFAIL > 0 ) THEN
      WRITE( NOUT, FMT = 9996 )'DSPOSV', NFAIL, NRUN
   ELSE
      WRITE( NOUT, FMT = 9995 )'DSPOSV', NRUN
   END IF
   IF( NERRS > 0 ) THEN
      WRITE( NOUT, FMT = 9994 )NERRS
   END IF
!
 9998 FORMAT( ' UPLO=''', A1, ''', N =', I5, ', NRHS=', I3, ', type ', &
         I2, ', test(', I2, ') =', G12.5 )
 9996 FORMAT( 1X, A6, ': ', I6, ' out of ', I6, &
         ' tests failed to pass the threshold' )
 9995 FORMAT( /1X, 'All tests for ', A6, &
         ' routines passed the threshold ( ', I6, ' tests run)' )
 9994 FORMAT( 6X, I6, ' error messages recorded' )
!
!     SUBNAM, INFO, INFOE, N, IMAT
!
 9988 FORMAT( ' *** ', A6, ' returned with INFO =', I5, ' instead of ', &
         I5, / ' ==> N =', I5, ', type ', &
         I2 )
!
!     SUBNAM, INFO, N, IMAT
!
 9975 FORMAT( ' *** Error code from ', A6, '=', I5, ' for M=', I5, &
         ', type ', I2 )
 8999 FORMAT( / 1X, A3, ':  positive definite dense matrices' )
 8979 FORMAT( 4X, '1. Diagonal', 24X, '7. Last n/2 columns zero', / 4X, &
         '2. Upper triangular', 16X, &
         '8. Random, CNDNUM = sqrt(0.1/EPS)', / 4X, &
         '3. Lower triangular', 16X, '9. Random, CNDNUM = 0.1/EPS', &
         / 4X, '4. Random, CNDNUM = 2', 13X, &
         '10. Scaled near underflow', / 4X, '5. First column zero', &
         14X, '11. Scaled near overflow', / 4X, &
         '6. Last column zero' )
 8960 FORMAT( 3X, I2, ': norm_1( B - A * X )  / ', &
         '( norm_1(A) * norm_1(X) * EPS * SQRT(N) ) > 1 if ITERREF', &
         / 4x, 'or norm_1( B - A * X )  / ', &
         '( norm_1(A) * norm_1(X) * EPS ) > THRES if DPOTRF' )

   RETURN
!
!     End of DDRVAC
!
END
                                                                                                                                                                                                                                                                                                            




