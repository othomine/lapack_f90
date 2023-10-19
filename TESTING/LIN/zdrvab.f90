!> \brief \b ZDRVAB
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZDRVAB( DOTYPE, NM, MVAL, NNS,
!                          NSVAL, THRESH, NMAX, A, AFAC, B,
!                          X, WORK, RWORK, SWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       INTEGER            NM, NMAX, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            MVAL( * ), NSVAL( * ), IWORK( * )
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX            SWORK( * )
!       COMPLEX*16         A( * ), AFAC( * ), B( * ),
!      $                   WORK( * ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZDRVAB tests ZCGESV
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
!>          The maximum value permitted for M or N, used in dimensioning
!>          the work arrays.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (NMAX*NSMAX)
!>          where NSMAX is the largest entry in NSVAL.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension
!>                      (NMAX*max(3,NSMAX*2))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension
!>                      NMAX
!> \endverbatim
!>
!> \param[out] SWORK
!> \verbatim
!>          SWORK is COMPLEX array, dimension
!>                      (NMAX*(NSMAX+NMAX))
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension
!>                      NMAX
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
!> \ingroup complex16_lin
!
!  =====================================================================
   SUBROUTINE ZDRVAB( DOTYPE, NM, MVAL, NNS, &
                      NSVAL, THRESH, NMAX, A, AFAC, B, &
                      X, WORK, RWORK, SWORK, IWORK, NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            NM, NMAX, NNS, NOUT
   DOUBLE PRECISION   THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            MVAL( * ), NSVAL( * ), IWORK( * )
   DOUBLE PRECISION   RWORK( * )
   COMPLEX            SWORK( * )
   COMPLEX*16         A( * ), AFAC( * ), B( * ), &
                      WORK( * ), X( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO
   PARAMETER          ( ZERO = 0.0D+0 )
   INTEGER            NTYPES
   PARAMETER          ( NTYPES = 11 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 1 )
!     ..
!     .. Local Scalars ..
   LOGICAL            ZEROT
   CHARACTER          DIST, TRANS, TYPE, XTYPE
   CHARACTER*3        PATH
   INTEGER            I, IM, IMAT, INFO, IOFF, IRHS, &
                      IZERO, KL, KU, LDA, M, MODE, N, &
                      NERRS, NFAIL, NIMAT, NRHS, NRUN
   DOUBLE PRECISION   ANORM, CNDNUM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   DOUBLE PRECISION   RESULT( NTESTS )
!     ..
!     .. Local Variables ..
   INTEGER            ITER, KASE
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAERH, ALAHD, ZGET08, ZLACPY, ZLARHS, ZLASET, &
                      ZLATB4, ZLATMS
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          DCMPLX, DBLE, MAX, MIN, SQRT
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
   DATA               ISEEDY / 2006, 2007, 2008, 2009 /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   KASE = 0
   PATH( 1: 1 ) = 'Zomplex precision'
   PATH( 2: 3 ) = 'GE'
   NRUN = 0
   NFAIL = 0
   NERRS = 0
   DO I = 1, 4
      ISEED( I ) = ISEEDY( I )
   ENDDO
!
   INFOT = 0
!
!     Do for each value of M in MVAL
!
   DO IM = 1, NM
      M = MVAL( IM )
      LDA = MAX( 1, M )
!
      N = M
      NIMAT = NTYPES
      IF( M <= 0 .OR. N <= 0 ) &
         NIMAT = 1
!
      DO IMAT = 1, NIMAT
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
         IF( .NOT.DOTYPE( IMAT ) ) &
            GO TO 100
!
!           Skip types 5, 6, or 7 if the matrix size is too small.
!
         ZEROT = IMAT >= 5 .AND. IMAT <= 7
         IF( ZEROT .AND. N < IMAT-4 ) &
            GO TO 100
!
!           Set up parameters with ZLATB4 and generate a test matrix
!           with ZLATMS.
!
         CALL ZLATB4( PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE, &
                      CNDNUM, DIST )
!
         SRNAMT = 'ZLATMS'
         CALL ZLATMS( M, N, DIST, ISEED, TYPE, RWORK, MODE, &
                      CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, &
                      WORK, INFO )
!
!           Check error code from ZLATMS.
!
         IF( INFO /= 0 ) THEN
            CALL ALAERH( PATH, 'ZLATMS', INFO, 0, ' ', M, N, -1, &
                         -1, -1, IMAT, NFAIL, NERRS, NOUT )
            GO TO 100
         END IF
!
!           For types 5-7, zero one or more columns of the matrix to
!           test that INFO is returned correctly.
!
         IF( ZEROT ) THEN
            IF( IMAT == 5 ) THEN
               IZERO = 1
            ELSE IF( IMAT == 6 ) THEN
               IZERO = MIN( M, N )
            ELSE
               IZERO = MIN( M, N ) / 2 + 1
            END IF
            IOFF = ( IZERO-1 )*LDA
            IF( IMAT < 7 ) THEN
               DO I = 1, M
                  A( IOFF+I ) = ZERO
               ENDDO
            ELSE
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLASET( 'Full', M, N-IZERO+1, DCMPLX(ZERO), &
                            DCMPLX(ZERO), A( IOFF+1 ), LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
            END IF
         ELSE
            IZERO = 0
         END IF
!
         DO IRHS = 1, NNS
            NRHS = NSVAL( IRHS )
            XTYPE = 'N'
            TRANS = 'N'
!
            SRNAMT = 'ZLARHS'
            CALL ZLARHS( PATH, XTYPE, ' ', TRANS, N, N, KL, &
                         KU, NRHS, A, LDA, X, LDA, B, &
                         LDA, ISEED, INFO )
!
            SRNAMT = 'ZCGESV'
!
            KASE = KASE + 1
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZLACPY( 'Full', M, N, A, LDA, AFAC, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZCGESV( N, NRHS, A, LDA, IWORK, B, LDA, X, LDA, &
                         WORK, SWORK, RWORK, ITER, INFO)
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZCGESV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
            IF (ITER < 0) THEN
#ifdef _TIMER
                call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                CALL ZLACPY( 'Full', M, N, AFAC, LDA, A, LDA )
#ifdef _TIMER
                call system_clock(count_rate=nb_periods_sec,count=S2_time)
                open(file='results.out', unit=10, position = 'append')
                write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                      real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                close(10)
#endif
            ENDIF
!
!              Check error code from ZCGESV. This should be the same as
!              the one of DGETRF.
!
            IF( INFO /= IZERO ) THEN
!
               IF( NFAIL == 0 .AND. NERRS == 0 ) &
                  CALL ALAHD( NOUT, PATH )
               NERRS = NERRS + 1
!
               IF( INFO /= IZERO .AND. IZERO /= 0 ) THEN
                  WRITE( NOUT, FMT = 9988 )'ZCGESV',INFO, &
                            IZERO,M,IMAT
               ELSE
                  WRITE( NOUT, FMT = 9975 )'ZCGESV',INFO, &
                            M, IMAT
               END IF
            END IF
!
!              Skip the remaining test if the matrix is singular.
!
            IF( INFO /= 0 ) &
               GO TO 100
!
!              Check the quality of the solution
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
            CALL ZGET08( TRANS, N, N, NRHS, A, LDA, X, LDA, WORK, &
                         LDA, RWORK, RESULT( 1 ) )
!
!              Check if the test passes the testing.
!              Print information about the tests that did not
!              pass the testing.
!
!              If iterative refinement has been used and claimed to
!              be successful (ITER>0), we want
!                NORMI(B - A*X)/(NORMI(A)*NORMI(X)*EPS*SRQT(N)) < 1
!
!              If double precision has been used (ITER<0), we want
!                NORMI(B - A*X)/(NORMI(A)*NORMI(X)*EPS) < THRES
!              (Cf. the linear solver testing routines)
!
            IF ((THRESH <= 0.0E+00) &
               .OR.((ITER >= 0).AND.(N > 0) &
                    .AND.(RESULT(1) >= SQRT(DBLE(N)))) &
               .OR.((ITER < 0).AND.(RESULT(1) >= THRESH))) THEN
!
               IF( NFAIL == 0 .AND. NERRS == 0 ) THEN
                  WRITE( NOUT, FMT = 8999 )'DGE'
                  WRITE( NOUT, FMT = '( '' Matrix types:'' )' )
                  WRITE( NOUT, FMT = 8979 )
                  WRITE( NOUT, FMT = '( '' Test ratios:'' )' )
                  WRITE( NOUT, FMT = 8960 )1
                  WRITE( NOUT, FMT = '( '' Messages:'' )' )
               END IF
!
               WRITE( NOUT, FMT = 9998 )TRANS, N, NRHS, &
                  IMAT, 1, RESULT( 1 )
               NFAIL = NFAIL + 1
            END IF
            NRUN = NRUN + 1
         ENDDO
  100    CONTINUE
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   IF( NFAIL > 0 ) THEN
      WRITE( NOUT, FMT = 9996 )'ZCGESV', NFAIL, NRUN
   ELSE
      WRITE( NOUT, FMT = 9995 )'ZCGESV', NRUN
   END IF
   IF( NERRS > 0 ) THEN
      WRITE( NOUT, FMT = 9994 )NERRS
   END IF
!
 9998 FORMAT( ' TRANS=''', A1, ''', N =', I5, ', NRHS=', I3, ', type ', &
         I2, ', test(', I2, ') =', G12.5 )
 9996 FORMAT( 1X, A6, ': ', I6, ' out of ', I6, &
         ' tests failed to pass the threshold' )
 9995 FORMAT( /1X, 'All tests for ', A6, &
         ' routines passed the threshold ( ', I6, ' tests run)' )
 9994 FORMAT( 6X, I6, ' error messages recorded' )
!
!     SUBNAM, INFO, INFOE, M, IMAT
!
 9988 FORMAT( ' *** ', A6, ' returned with INFO =', I5, ' instead of ', &
         I5, / ' ==> M =', I5, ', type ', &
         I2 )
!
!     SUBNAM, INFO, M, IMAT
!
 9975 FORMAT( ' *** Error code from ', A6, '=', I5, ' for M=', I5, &
         ', type ', I2 )
 8999 FORMAT( / 1X, A3, ':  General dense matrices' )
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
         '( norm_1(A) * norm_1(X) * EPS ) > THRES if DGETRF' )
   RETURN
!
!     End of ZDRVAB
!
END
                                                                                                                                                                                                                                                                                                            




