!> \brief \b ZCHKAA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       PROGRAM ZCHKAA
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZCHKAA is the main test program for the COMPLEX*16 linear equation
!> routines.
!>
!> The program must be driven by a short data file. The first 15 records
!> (not including the first comment  line) specify problem dimensions
!> and program options using list-directed input. The remaining lines
!> specify the LAPACK test paths and the number of matrix types to use
!> in testing.  An annotated example of a data file can be obtained by
!> deleting the first 3 characters from the following 42 lines:
!> Data file for testing COMPLEX*16 LAPACK linear equation routines
!> 7                      Number of values of M
!> 0 1 2 3 5 10 16        Values of M (row dimension)
!> 7                      Number of values of N
!> 0 1 2 3 5 10 16        Values of N (column dimension)
!> 1                      Number of values of NRHS
!> 2                      Values of NRHS (number of right hand sides)
!> 5                      Number of values of NB
!> 1 3 3 3 20             Values of NB (the blocksize)
!> 1 0 5 9 1              Values of NX (crossover point)
!> 3                      Number of values of RANK
!> 30 50 90               Values of rank (as a % of N)
!> 30.0                   Threshold value of test ratio
!> T                      Put T to test the LAPACK routines
!> T                      Put T to test the driver routines
!> T                      Put T to test the error exits
!> ZGE   11               List types on next line if 0 < NTYPES < 11
!> ZGB    8               List types on next line if 0 < NTYPES <  8
!> ZGT   12               List types on next line if 0 < NTYPES < 12
!> ZPO    9               List types on next line if 0 < NTYPES <  9
!> ZPS    9               List types on next line if 0 < NTYPES <  9
!> ZPP    9               List types on next line if 0 < NTYPES <  9
!> ZPB    8               List types on next line if 0 < NTYPES <  8
!> ZPT   12               List types on next line if 0 < NTYPES < 12
!> ZHE   10               List types on next line if 0 < NTYPES < 10
!> ZHR   10               List types on next line if 0 < NTYPES < 10
!> ZHK   10               List types on next line if 0 < NTYPES < 10
!> ZHA   10               List types on next line if 0 < NTYPES < 10
!> ZH2   10               List types on next line if 0 < NTYPES < 10
!> ZSA   11               List types on next line if 0 < NTYPES < 10
!> ZS2   11               List types on next line if 0 < NTYPES < 10
!> ZHP   10               List types on next line if 0 < NTYPES < 10
!> ZSY   11               List types on next line if 0 < NTYPES < 11
!> ZSR   11               List types on next line if 0 < NTYPES < 11
!> ZSK   11               List types on next line if 0 < NTYPES < 11
!> ZSP   11               List types on next line if 0 < NTYPES < 11
!> ZTR   18               List types on next line if 0 < NTYPES < 18
!> ZTP   18               List types on next line if 0 < NTYPES < 18
!> ZTB   17               List types on next line if 0 < NTYPES < 17
!> ZQR    8               List types on next line if 0 < NTYPES <  8
!> ZRQ    8               List types on next line if 0 < NTYPES <  8
!> ZLQ    8               List types on next line if 0 < NTYPES <  8
!> ZQL    8               List types on next line if 0 < NTYPES <  8
!> ZQP    6               List types on next line if 0 < NTYPES <  6
!> ZTZ    3               List types on next line if 0 < NTYPES <  3
!> ZLS    6               List types on next line if 0 < NTYPES <  6
!> ZEQ
!> ZQT
!> ZQX
!> ZTS
!> ZHH
!> \endverbatim
!
!  Parameters:
!  ==========
!
!> \verbatim
!>  NMAX    INTEGER
!>          The maximum allowable value for M and N.
!>
!>  MAXIN   INTEGER
!>          The number of different values that can be used for each of
!>          M, N, NRHS, NB, NX and RANK
!>
!>  MAXRHS  INTEGER
!>          The maximum number of right hand sides
!>
!>  MATMAX  INTEGER
!>          The maximum number of matrix types to use for testing
!>
!>  NIN     INTEGER
!>          The unit number for input
!>
!>  NOUT    INTEGER
!>          The unit number for output
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
!> \ingroup complex16_lin
!
!  =====================================================================
   PROGRAM ZCHKAA
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            NMAX
   PARAMETER          ( NMAX = 132 )
   INTEGER            MAXIN
   PARAMETER          ( MAXIN = 12 )
   INTEGER            MAXRHS
   PARAMETER          ( MAXRHS = 16 )
   INTEGER            MATMAX
   PARAMETER          ( MATMAX = 30 )
   INTEGER            NIN, NOUT
   PARAMETER          ( NIN = 5, NOUT = 6 )
   INTEGER            KDMAX
   PARAMETER          ( KDMAX = NMAX+( NMAX+1 ) / 4 )
!     ..
!     .. Local Scalars ..
   LOGICAL            FATAL, TSTCHK, TSTDRV, TSTERR
   CHARACTER          C1
   CHARACTER*2        C2
   CHARACTER*3        PATH
   CHARACTER*10       INTSTR
   CHARACTER*72       ALINE
   INTEGER            I, IC, J, K, LA, LAFAC, LDA, NB, NM, NMATS, NN, &
                      NNB, NNB2, NNS, NRHS, NTYPES, NRANK, &
                      VERS_MAJOR, VERS_MINOR, VERS_PATCH
   DOUBLE PRECISION   EPS, THREQ, THRESH
   INTEGER(8)         nb_periods_sec, S1, S2, S1T, S2T
   REAL               STOT
!     ..
!     .. Local Arrays ..
   LOGICAL            DOTYPE( MATMAX )
   INTEGER            IWORK( 25*NMAX ), MVAL( MAXIN ), &
                      NBVAL( MAXIN ), NBVAL2( MAXIN ), &
                      NSVAL( MAXIN ), NVAL( MAXIN ), NXVAL( MAXIN ), &
                      RANKVAL( MAXIN ), PIV( NMAX )
   DOUBLE PRECISION   S( 2*NMAX )
   COMPLEX*16         E( NMAX )
!
!    .. Allocatable Arrays ..
   INTEGER AllocateStatus
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::   RWORK
   COMPLEX*16, DIMENSION(:,:), ALLOCATABLE::  A, B, WORK
!     ..
!     .. External Functions ..
   LOGICAL            LSAME, LSAMEN
   DOUBLE PRECISION   DLAMCH, DSECND
   EXTERNAL           LSAME, LSAMEN, DLAMCH, DSECND
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAREQ, ZCHKEQ, ZCHKGB, ZCHKGE, ZCHKGT, ZCHKHE, &
                      ZCHKHE_ROOK, ZCHKHE_RK, ZCHKHE_AA, ZCHKHP, &
                      ZCHKLQ, ZCHKUNHR_COL, ZCHKPB, ZCHKPO, ZCHKPS, &
                      ZCHKPP, ZCHKPT, ZCHKQ3, ZCHKQL, ZCHKQR, ZCHKRQ, &
                      ZCHKSP, ZCHKSY, ZCHKSY_ROOK, ZCHKSY_RK, &
                      ZCHKSY_AA, ZCHKTB, ZCHKTP, ZCHKTR, ZCHKTZ, &
                      ZDRVGB, ZDRVGE, ZDRVGT, ZDRVHE, ZDRVHE_ROOK, &
                      ZDRVHE_RK, ZDRVHE_AA, ZDRVHE_AA_2STAGE, ZDRVHP, &
                      ZDRVLS, ZDRVPB,  ZDRVPO, ZDRVPP, ZDRVPT, &
                      ZDRVSP, ZDRVSY, ZDRVSY_ROOK, ZDRVSY_RK, &
                      ZDRVSY_AA, ZDRVSY_AA_2STAGE, ILAVER, ZCHKQRT, &
                      ZCHKQRTP, ZCHKLQT, ZCHKLQTP, ZCHKTSQR
!     ..
!     .. Scalars in Common ..
   LOGICAL            LERR, OK
   CHARACTER*32       SRNAMT
   INTEGER            INFOT, NUNIT
!     ..
!     .. Arrays in Common ..
   INTEGER            IPARMS( 100 )
!     ..
!     .. Common blocks ..
   COMMON             / INFOC / INFOT, NUNIT, OK, LERR
   COMMON             / SRNAMC / SRNAMT
   COMMON             / CLAENV / IPARMS
!     ..
!     .. Data statements ..
   DATA               THREQ / 2.0D0 / , INTSTR / '0123456789' /
!
!     .. Allocate memory dynamically ..
   ALLOCATE (RWORK( 150*NMAX+2*MAXRHS ), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
   ALLOCATE (A ((KDMAX+1) * NMAX, 7), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
   ALLOCATE (B (NMAX * MAXRHS, 4), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
   ALLOCATE (WORK (NMAX, NMAX+MAXRHS+10), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!     ..
!     .. Executable Statements ..
!
   call system_clock(count_rate=nb_periods_sec,count=S1T)
   LDA = NMAX
   FATAL = .FALSE.
!
!     Read a dummy line.
!
   READ( NIN, FMT = * )
!
!     Report values of parameters.
!
   CALL ILAVER( VERS_MAJOR, VERS_MINOR, VERS_PATCH )
   WRITE( NOUT, FMT = 9994 ) VERS_MAJOR, VERS_MINOR, VERS_PATCH
!
!     Read the values of M
!
   READ( NIN, FMT = * )NM
   IF( NM < 1 ) THEN
      WRITE( NOUT, FMT = 9996 )' NM ', NM, 1
      NM = 0
      FATAL = .TRUE.
   ELSE IF( NM > MAXIN ) THEN
      WRITE( NOUT, FMT = 9995 )' NM ', NM, MAXIN
      NM = 0
      FATAL = .TRUE.
   END IF
   READ( NIN, FMT = * )( MVAL( I ), I = 1, NM )
   DO I = 1, NM
      IF( MVAL( I ) < 0 ) THEN
         WRITE( NOUT, FMT = 9996 )' M  ', MVAL( I ), 0
         FATAL = .TRUE.
      ELSE IF( MVAL( I ) > NMAX ) THEN
         WRITE( NOUT, FMT = 9995 )' M  ', MVAL( I ), NMAX
         FATAL = .TRUE.
      END IF
   ENDDO
   IF( NM > 0 ) &
      WRITE( NOUT, FMT = 9993 )'M   ', ( MVAL( I ), I = 1, NM )
!
!     Read the values of N
!
   READ( NIN, FMT = * )NN
   IF( NN < 1 ) THEN
      WRITE( NOUT, FMT = 9996 )' NN ', NN, 1
      NN = 0
      FATAL = .TRUE.
   ELSE IF( NN > MAXIN ) THEN
      WRITE( NOUT, FMT = 9995 )' NN ', NN, MAXIN
      NN = 0
      FATAL = .TRUE.
   END IF
   READ( NIN, FMT = * )( NVAL( I ), I = 1, NN )
   DO I = 1, NN
      IF( NVAL( I ) < 0 ) THEN
         WRITE( NOUT, FMT = 9996 )' N  ', NVAL( I ), 0
         FATAL = .TRUE.
      ELSE IF( NVAL( I ) > NMAX ) THEN
         WRITE( NOUT, FMT = 9995 )' N  ', NVAL( I ), NMAX
         FATAL = .TRUE.
      END IF
   ENDDO
   IF( NN > 0 ) &
      WRITE( NOUT, FMT = 9993 )'N   ', ( NVAL( I ), I = 1, NN )
!
!     Read the values of NRHS
!
   READ( NIN, FMT = * )NNS
   IF( NNS < 1 ) THEN
      WRITE( NOUT, FMT = 9996 )' NNS', NNS, 1
      NNS = 0
      FATAL = .TRUE.
   ELSE IF( NNS > MAXIN ) THEN
      WRITE( NOUT, FMT = 9995 )' NNS', NNS, MAXIN
      NNS = 0
      FATAL = .TRUE.
   END IF
   READ( NIN, FMT = * )( NSVAL( I ), I = 1, NNS )
   DO I = 1, NNS
      IF( NSVAL( I ) < 0 ) THEN
         WRITE( NOUT, FMT = 9996 )'NRHS', NSVAL( I ), 0
         FATAL = .TRUE.
      ELSE IF( NSVAL( I ) > MAXRHS ) THEN
         WRITE( NOUT, FMT = 9995 )'NRHS', NSVAL( I ), MAXRHS
         FATAL = .TRUE.
      END IF
   ENDDO
   IF( NNS > 0 ) &
      WRITE( NOUT, FMT = 9993 )'NRHS', ( NSVAL( I ), I = 1, NNS )
!
!     Read the values of NB
!
   READ( NIN, FMT = * )NNB
   IF( NNB < 1 ) THEN
      WRITE( NOUT, FMT = 9996 )'NNB ', NNB, 1
      NNB = 0
      FATAL = .TRUE.
   ELSE IF( NNB > MAXIN ) THEN
      WRITE( NOUT, FMT = 9995 )'NNB ', NNB, MAXIN
      NNB = 0
      FATAL = .TRUE.
   END IF
   READ( NIN, FMT = * )( NBVAL( I ), I = 1, NNB )
   DO I = 1, NNB
      IF( NBVAL( I ) < 0 ) THEN
         WRITE( NOUT, FMT = 9996 )' NB ', NBVAL( I ), 0
         FATAL = .TRUE.
      END IF
   ENDDO
   IF( NNB > 0 ) &
      WRITE( NOUT, FMT = 9993 )'NB  ', ( NBVAL( I ), I = 1, NNB )
!
!     Set NBVAL2 to be the set of unique values of NB
!
   NNB2 = 0
   DO I = 1, NNB
      NB = NBVAL( I )
      DO J = 1, NNB2
         IF( NB == NBVAL2( J ) ) &
            GO TO 60
      ENDDO
      NNB2 = NNB2 + 1
      NBVAL2( NNB2 ) = NB
60 CONTINUE
   ENDDO
!
!     Read the values of NX
!
   READ( NIN, FMT = * )( NXVAL( I ), I = 1, NNB )
   DO I = 1, NNB
      IF( NXVAL( I ) < 0 ) THEN
         WRITE( NOUT, FMT = 9996 )' NX ', NXVAL( I ), 0
         FATAL = .TRUE.
      END IF
   ENDDO
   IF( NNB > 0 ) &
      WRITE( NOUT, FMT = 9993 )'NX  ', ( NXVAL( I ), I = 1, NNB )
!
!     Read the values of RANKVAL
!
   READ( NIN, FMT = * )NRANK
   IF( NN < 1 ) THEN
      WRITE( NOUT, FMT = 9996 )' NRANK ', NRANK, 1
      NRANK = 0
      FATAL = .TRUE.
   ELSE IF( NN > MAXIN ) THEN
      WRITE( NOUT, FMT = 9995 )' NRANK ', NRANK, MAXIN
      NRANK = 0
      FATAL = .TRUE.
   END IF
   READ( NIN, FMT = * )( RANKVAL( I ), I = 1, NRANK )
   DO I = 1, NRANK
      IF( RANKVAL( I ) < 0 ) THEN
         WRITE( NOUT, FMT = 9996 )' RANK  ', RANKVAL( I ), 0
         FATAL = .TRUE.
      ELSE IF( RANKVAL( I ) > 100 ) THEN
         WRITE( NOUT, FMT = 9995 )' RANK  ', RANKVAL( I ), 100
         FATAL = .TRUE.
      END IF
   END DO
   IF( NRANK > 0 ) &
      WRITE( NOUT, FMT = 9993 )'RANK % OF N', &
      ( RANKVAL( I ), I = 1, NRANK )
!
!     Read the threshold value for the test ratios.
!
   READ( NIN, FMT = * )THRESH
   WRITE( NOUT, FMT = 9992 )THRESH
!
!     Read the flag that indicates whether to test the LAPACK routines.
!
   READ( NIN, FMT = * )TSTCHK
!
!     Read the flag that indicates whether to test the driver routines.
!
   READ( NIN, FMT = * )TSTDRV
!
!     Read the flag that indicates whether to test the error exits.
!
   READ( NIN, FMT = * )TSTERR
!
   IF( FATAL ) THEN
      WRITE( NOUT, FMT = 9999 )
      STOP
   END IF
!
!     Calculate and print the machine dependent constants.
!
   EPS = DLAMCH( 'Underflow threshold' )
   WRITE( NOUT, FMT = 9991 )'underflow', EPS
   EPS = DLAMCH( 'Overflow threshold' )
   WRITE( NOUT, FMT = 9991 )'overflow ', EPS
   EPS = DLAMCH( 'Epsilon' )
   WRITE( NOUT, FMT = 9991 )'precision', EPS
   WRITE( NOUT, FMT = * )
   NRHS = NSVAL( 1 )
!
80 CONTINUE
!
!     Read a test path and the number of matrix types to use.
!
   READ( NIN, FMT = '(A72)', END = 140 )ALINE
   PATH = ALINE( 1: 3 )
   NMATS = MATMAX
   I = 3
90 CONTINUE
   I = I + 1
   IF( I > 72 ) &
      GO TO 130
   IF( ALINE( I: I ) == ' ' ) &
      GO TO 90
   NMATS = 0
  100 CONTINUE
   C1 = ALINE( I: I )
   DO K = 1, 10
      IF( C1 == INTSTR( K: K ) ) THEN
         IC = K - 1
         GO TO 120
      END IF
      ENDDO
   GO TO 130
  120 CONTINUE
   NMATS = NMATS*10 + IC
   I = I + 1
   IF( I > 72 ) &
      GO TO 130
   GO TO 100
  130 CONTINUE
   C1 = PATH( 1: 1 )
   C2 = PATH( 2: 3 )
!
!     Check first character for correct precision.
!
   IF( .NOT.LSAME( C1, 'Zomplex precision' ) ) THEN
      WRITE( NOUT, FMT = 9990 )PATH
!
   ELSE IF( NMATS <= 0 ) THEN
!
!        Check for a positive number of tests requested.
!
      WRITE( NOUT, FMT = 9989 )PATH
!
   ELSE IF( LSAMEN( 2, C2, 'GE' ) ) THEN
!
!        GE:  general matrices
!
      NTYPES = 11
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKGE( DOTYPE, NM, MVAL, NN, NVAL, NNB2, NBVAL2, NNS, &
                      NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), &
                      A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), &
                      B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKGE : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVGE( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, &
                      A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, &
                      RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVGE : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'GB' ) ) THEN
!
!        GB:  general banded matrices
!
      LA = ( 2*KDMAX+1 )*NMAX
      LAFAC = ( 3*KDMAX+1 )*NMAX
      NTYPES = 8
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKGB( DOTYPE, NM, MVAL, NN, NVAL, NNB2, NBVAL2, NNS, &
                      NSVAL, THRESH, TSTERR, A( 1, 1 ), LA, &
                      A( 1, 3 ), LAFAC, B( 1, 1 ), B( 1, 2 ), &
                      B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKGB : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVGB( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                      A( 1, 1 ), LA, A( 1, 3 ), LAFAC, A( 1, 6 ), &
                      B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, &
                      WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVGB : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'GT' ) ) THEN
!
!        GT:  general tridiagonal matrices
!
      NTYPES = 12
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKGT( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), &
                      B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKGT : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVGT( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                      A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), &
                      B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVGT : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'PO' ) ) THEN
!
!        PO:  positive definite matrices
!
      NTYPES = 9
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKPO( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                      THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                      A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                      WORK, RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKPO : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVPO( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, &
                      A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, &
                      RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVPO : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'PS' ) ) THEN
!
!        PS:  positive semi-definite matrices
!
      NTYPES = 9
!
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKPS( DOTYPE, NN, NVAL, NNB2, NBVAL2, NRANK, &
                      RANKVAL, THRESH, TSTERR, LDA, A( 1, 1 ), &
                      A( 1, 2 ), A( 1, 3 ), PIV, WORK, RWORK, &
                      NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKPS : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'PP' ) ) THEN
!
!        PP:  positive definite packed matrices
!
      NTYPES = 9
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKPP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                      B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, &
                      NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKPP : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVPP( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, &
                      A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, &
                      RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVPP : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'PB' ) ) THEN
!
!        PB:  positive definite banded matrices
!
      NTYPES = 8
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKPB( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                      THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                      A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                      WORK, RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKPB : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVPB( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, &
                      A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, &
                      RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVPB : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'PT' ) ) THEN
!
!        PT:  positive definite tridiagonal matrices
!
      NTYPES = 12
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKPT( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      A( 1, 1 ), S, A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), &
                      B( 1, 3 ), WORK, RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKPT : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVPT( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                      A( 1, 1 ), S, A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), &
                      B( 1, 3 ), WORK, RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVPT : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'HE' ) ) THEN
!
!        HE:  Hermitian indefinite matrices
!
      NTYPES = 10
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKHE( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                      THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                      A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                      WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKHE : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVHE( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, &
                      A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, &
                      NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVHE : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF

   ELSE IF( LSAMEN( 2, C2, 'HR' ) ) THEN
!
!        HR:  Hermitian indefinite matrices,
!             with bounded Bunch-Kaufman (rook) pivoting algorithm,
!
      NTYPES = 10
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKHE_ROOK(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                          THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                          A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                          WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKHE_ROOK : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVHE_ROOK( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                           LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                           B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, &
                           RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVHE_ROOK : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'HK' ) ) THEN
!
!        HK:  Hermitian indefinite matrices,
!             with bounded Bunch-Kaufman (rook) pivoting algorithm,
!             different matrix storage format than HR path version.
!
      NTYPES = 10
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKHE_RK ( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                          THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                          E, A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), &
                          B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKHE_RK : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVHE_RK( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                         LDA, A( 1, 1 ), A( 1, 2 ), E, A( 1, 3 ), &
                         B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, &
                         RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVHE_RK : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'HA' ) ) THEN
!
!        HA:  Hermitian matrices,
!             Aasen Algorithm
!
      NTYPES = 10
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKHE_AA( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, &
                            NSVAL, THRESH, TSTERR, LDA, &
                            A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                            B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                            WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKHE_AA : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVHE_AA( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                            LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                                 B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                            WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVHE_AA : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'H2' ) ) THEN
!
!        H2:  Hermitian matrices,
!             with partial (Aasen's) pivoting algorithm
!
      NTYPES = 10
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKHE_AA_2STAGE( DOTYPE, NN, NVAL, NNB2, NBVAL2, &
                            NNS, NSVAL, THRESH, TSTERR, LDA, &
                            A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                            B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                            WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKHE_AA_2STAGE : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVHE_AA_2STAGE( &
                            DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                            LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                                 B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                            WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVHE_AA_2STAGE : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
!
   ELSE IF( LSAMEN( 2, C2, 'HP' ) ) THEN
!
!        HP:  Hermitian indefinite packed matrices
!
      NTYPES = 10
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKHP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                      B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, &
                      IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKHP : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVHP( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, &
                      A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, &
                      NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVHP : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'SY' ) ) THEN
!
!        SY:  symmetric indefinite matrices,
!             with partial (Bunch-Kaufman) pivoting algorithm
!
      NTYPES = 11
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKSY( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                      THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                      A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                      WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKSY : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVSY( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, &
                      A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, &
                      NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVSY : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'SR' ) ) THEN
!
!        SR:  symmetric indefinite matrices,
!             with bounded Bunch-Kaufman (rook) pivoting algorithm
!
      NTYPES = 11
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKSY_ROOK(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                          THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                          A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                          WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKSY_ROOK : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVSY_ROOK( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                           LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                           B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, &
                           RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVSY_ROOK : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'SK' ) ) THEN
!
!        SK:  symmetric indefinite matrices,
!             with bounded Bunch-Kaufman (rook) pivoting algorithm,
!             different matrix storage format than SR path version.
!
      NTYPES = 11
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKSY_RK( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                         THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                         E, A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), &
                         B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKSY_RK : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVSY_RK( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                         LDA, A( 1, 1 ), A( 1, 2 ), E, A( 1, 3 ), &
                         B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, &
                         RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVSY_RK : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'SA' ) ) THEN
!
!        SA:  symmetric indefinite matrices with Aasen's algorithm,
!
      NTYPES = 11
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKSY_AA( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                         THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                         A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), &
                         B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKSY_AA : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVSY_AA( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                         LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                         B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, &
                         RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVSY_AA : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'S2' ) ) THEN
!
!        S2:  symmetric indefinite matrices with Aasen's algorithm
!             2 stage
!
      NTYPES = 11
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKSY_AA_2STAGE( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, &
                         NSVAL, THRESH, TSTERR, LDA, &
                         A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                         B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                         WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKSY_AA_2STAGE : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVSY_AA_2STAGE( &
                         DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                         LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                         B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, &
                         RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVSY_AA_2STAGE : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'SP' ) ) THEN
!
!        SP:  symmetric indefinite packed matrices,
!             with partial (Bunch-Kaufman) pivoting algorithm
!
      NTYPES = 11
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKSP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                      B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, &
                      IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKSP : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVSP( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, &
                      A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, &
                      NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVSP : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'TR' ) ) THEN
!
!        TR:  triangular matrices
!
      NTYPES = 18
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKTR( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                      THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                      B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, &
                      NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKTR : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'TP' ) ) THEN
!
!        TP:  triangular packed matrices
!
      NTYPES = 18
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKTP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKTP : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'TB' ) ) THEN
!
!        TB:  triangular banded matrices
!
      NTYPES = 17
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKTB( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKTB : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'QR' ) ) THEN
!
!        QR:  QR factorization
!
      NTYPES = 8
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKQR( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, &
                      NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), &
                      A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), &
                      B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), &
                      WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKQR : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'LQ' ) ) THEN
!
!        LQ:  LQ factorization
!
      NTYPES = 8
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKLQ( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, &
                      NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), &
                      A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), &
                      B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), &
                      WORK, RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKLQ : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'QL' ) ) THEN
!
!        QL:  QL factorization
!
      NTYPES = 8
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKQL( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, &
                      NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), &
                      A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), &
                      B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), &
                      WORK, RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKQL : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'RQ' ) ) THEN
!
!        RQ:  RQ factorization
!
      NTYPES = 8
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKRQ( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, &
                      NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), &
                      A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), &
                      B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), &
                      WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKRQ : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'EQ' ) ) THEN
!
!        EQ:  Equilibration routines for general and positive definite
!             matrices (THREQ should be between 2 and 10)
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKEQ( THREQ, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKEQ : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'TZ' ) ) THEN
!
!        TZ:  Trapezoidal matrix
!
      NTYPES = 3
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKTZ( DOTYPE, NM, MVAL, NN, NVAL, THRESH, TSTERR, &
                      A( 1, 1 ), A( 1, 2 ), S( 1 ), &
                      B( 1, 1 ), WORK, RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKTZ : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'QP' ) ) THEN
!
!        QP:  QR factorization with pivoting
!
      NTYPES = 6
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKQ3( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, &
                      THRESH, A( 1, 1 ), A( 1, 2 ), S( 1 ), &
                      B( 1, 1 ), WORK, RWORK, IWORK, &
                      NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKQ3 : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'LS' ) ) THEN
!
!        LS:  Least squares drivers
!
      NTYPES = 6
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVLS( DOTYPE, NM, MVAL, NN, NVAL, NNS, NSVAL, NNB, &
                      NBVAL, NXVAL, THRESH, TSTERR, A( 1, 1 ), &
                      A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), &
                      S( 1 ), S( NMAX+1 ), NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVLS : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
!
   ELSE IF( LSAMEN( 2, C2, 'QT' ) ) THEN
!
!        QT:  QRT routines for general matrices
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKQRT( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, &
                       NBVAL, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKQRT : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'QX' ) ) THEN
!
!        QX:  QRT routines for triangular-pentagonal matrices
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKQRTP( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, &
                        NBVAL, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKQRTP : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'TQ' ) ) THEN
!
!        TQ:  LQT routines for general matrices
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKLQT( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, &
                       NBVAL, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKLQT : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'XQ' ) ) THEN
!
!        XQ:  LQT routines for triangular-pentagonal matrices
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKLQTP( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, &
                        NBVAL, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKLQTP : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'TS' ) ) THEN
!
!        TS:  QR routines for tall-skinny matrices
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKTSQR( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, &
                        NBVAL, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKTSQR : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'TQ' ) ) THEN
!
!        TQ:  LQT routines for general matrices
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKLQT( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, &
                       NBVAL, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKLQT : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'XQ' ) ) THEN
!
!        XQ:  LQT routines for triangular-pentagonal matrices
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKLQTP( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, &
                        NBVAL, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKLQTP : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'TS' ) ) THEN
!
!        TS:  QR routines for tall-skinny matrices
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKTSQR( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, &
                        NBVAL, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKTSQR : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'HH' ) ) THEN
!
!        HH:  Householder reconstruction for tall-skinny matrices
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKUNHR_COL( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, &
                            NBVAL, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKUNHR_COL : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 ) PATH
      END IF
!
   ELSE
!
      WRITE( NOUT, FMT = 9990 )PATH
   END IF
!
!     Go back to get another input line.
!
   GO TO 80
!
!     Branch to this line when the last record is read.
!
  140 CONTINUE
   CLOSE ( NIN )
   call system_clock(count_rate=nb_periods_sec,count=S2T)
   WRITE( NOUT, FMT = 9998 )
   WRITE( NOUT, FMT = 9997 ) real(S2T - S1T)/real(nb_periods_sec)
!
   DEALLOCATE (A, STAT = AllocateStatus)
   DEALLOCATE (B, STAT = AllocateStatus)
   DEALLOCATE (RWORK, STAT = AllocateStatus)
   DEALLOCATE (WORK,  STAT = AllocateStatus)
!
 9999 FORMAT( / ' Execution not attempted due to input errors' )
 9998 FORMAT( / ' End of tests' )
 9997 FORMAT( ' Total time used = ', F16.8, ' seconds', / )
 9996 FORMAT( ' Invalid input value: ', A4, '=', I6, '; must be >=', &
         I6 )
 9995 FORMAT( ' Invalid input value: ', A4, '=', I6, '; must be <=', &
         I6 )
 9994 FORMAT( ' Tests of the COMPLEX*16 LAPACK routines ', &
         / ' LAPACK VERSION ', I1, '.', I1, '.', I1, &
         / / ' The following parameter values will be used:' )
 9993 FORMAT( 4X, A4, ':  ', 10I6, / 11X, 10I6 )
 9992 FORMAT( / ' Routines pass computational tests if test ratio is ', &
         'less than', F8.2, / )
 9991 FORMAT( ' Relative machine ', A, ' is taken to be', D16.6 )
 9990 FORMAT( / 1X, A3, ':  Unrecognized path name' )
 9989 FORMAT( / 1X, A3, ' routines were not tested' )
 9988 FORMAT( / 1X, A3, ' driver routines were not tested' )
!
!     End of ZCHKAA
!
   END
