!> \brief \b CCHKAA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       PROGRAM CCHKAA
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CCHKAA is the main test program for the COMPLEX linear equation
!> routines.
!>
!> The program must be driven by a short data file. The first 15 records
!> (not including the first comment  line) specify problem dimensions
!> and program options using list-directed input. The remaining lines
!> specify the LAPACK test paths and the number of matrix types to use
!> in testing.  An annotated example of a data file can be obtained by
!> deleting the first 3 characters from the following 42 lines:
!> Data file for testing COMPLEX LAPACK linear equation routines
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
!> CGE   11               List types on next line if 0 < NTYPES < 11
!> CGB    8               List types on next line if 0 < NTYPES <  8
!> CGT   12               List types on next line if 0 < NTYPES < 12
!> CPO    9               List types on next line if 0 < NTYPES <  9
!> CPO    9               List types on next line if 0 < NTYPES <  9
!> CPP    9               List types on next line if 0 < NTYPES <  9
!> CPB    8               List types on next line if 0 < NTYPES <  8
!> CPT   12               List types on next line if 0 < NTYPES < 12
!> CHE   10               List types on next line if 0 < NTYPES < 10
!> CHR   10               List types on next line if 0 < NTYPES < 10
!> CHK   10               List types on next line if 0 < NTYPES < 10
!> CHA   10               List types on next line if 0 < NTYPES < 10
!> CH2   10               List types on next line if 0 < NTYPES < 10
!> CSA   11               List types on next line if 0 < NTYPES < 10
!> CS2   11               List types on next line if 0 < NTYPES < 10
!> CHP   10               List types on next line if 0 < NTYPES < 10
!> CSY   11               List types on next line if 0 < NTYPES < 11
!> CSK   11               List types on next line if 0 < NTYPES < 11
!> CSR   11               List types on next line if 0 < NTYPES < 11
!> CSP   11               List types on next line if 0 < NTYPES < 11
!> CTR   18               List types on next line if 0 < NTYPES < 18
!> CTP   18               List types on next line if 0 < NTYPES < 18
!> CTB   17               List types on next line if 0 < NTYPES < 17
!> CQR    8               List types on next line if 0 < NTYPES <  8
!> CRQ    8               List types on next line if 0 < NTYPES <  8
!> CLQ    8               List types on next line if 0 < NTYPES <  8
!> CQL    8               List types on next line if 0 < NTYPES <  8
!> CQP    6               List types on next line if 0 < NTYPES <  6
!> CTZ    3               List types on next line if 0 < NTYPES <  3
!> CLS    6               List types on next line if 0 < NTYPES <  6
!> CEQ
!> CQT
!> CQX
!> CTS
!> CHH
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
!> \ingroup complex_lin
!
!  =====================================================================
   PROGRAM CCHKAA
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
   REAL               EPS, THREQ, THRESH
   INTEGER(8)         nb_periods_sec, S1, S2, S1T, S2T
   REAL               STOT
!     ..
!     .. Local Arrays ..
   LOGICAL            DOTYPE( MATMAX )
   INTEGER            IWORK( 25*NMAX ), MVAL( MAXIN ), &
                      NBVAL( MAXIN ), NBVAL2( MAXIN ), &
                      NSVAL( MAXIN ), NVAL( MAXIN ), NXVAL( MAXIN ), &
                      RANKVAL( MAXIN ), PIV( NMAX )
   REAL               S( 2*NMAX )
   COMPLEX            E( NMAX )
!     ..
!     .. Allocatable Arrays ..
   INTEGER AllocateStatus
   REAL, DIMENSION(:), ALLOCATABLE :: RWORK
   COMPLEX, DIMENSION(:,:), ALLOCATABLE :: A, B, WORK
!     ..
!     .. External Functions ..
   LOGICAL            LSAME, LSAMEN
   REAL               SECOND, SLAMCH
   EXTERNAL           LSAME, LSAMEN, SECOND, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAREQ, CCHKEQ, CCHKGB, CCHKGE, CCHKGT, CCHKHE, &
                      CCHKHE_ROOK, CCHKHE_RK, CCHKHE_AA, CCHKHP, &
                      CCHKLQ, CCHKUNHR_COL, CCHKPB, CCHKPO, CCHKPS, &
                      CCHKPP, CCHKPT, CCHKQ3, CCHKQL, CCHKQR, CCHKRQ, &
                      CCHKSP, CCHKSY, CCHKSY_ROOK, CCHKSY_RK, &
                      CCHKSY_AA, CCHKTB,  CCHKTP, CCHKTR, CCHKTZ, &
                      CDRVGB, CDRVGE, CDRVGT, CDRVHE, CDRVHE_ROOK, &
                      CDRVHE_RK, CDRVHE_AA, CDRVHP, CDRVLS, CDRVPB, &
                      CDRVPO, CDRVPP, CDRVPT, CDRVSP, CDRVSY, &
                      CDRVSY_ROOK, CDRVSY_RK, CDRVSY_AA, ILAVER, &
                      CCHKQRT, CCHKQRTP
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
   COMMON             / CLAENV / IPARMS
   COMMON             / INFOC / INFOT, NUNIT, OK, LERR
   COMMON             / SRNAMC / SRNAMT
!     ..
!     .. Data statements ..
   DATA               THREQ / 2.0 / , INTSTR / '0123456789' /
!     ..
!     .. Allocate memory dynamically ..
!
   ALLOCATE ( A( ( KDMAX+1 )*NMAX, 7 ), STAT = AllocateStatus )
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
   ALLOCATE ( B( NMAX*MAXRHS, 4 ), STAT = AllocateStatus )
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
   ALLOCATE ( WORK( NMAX, NMAX+MAXRHS+10 ), STAT = AllocateStatus )
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
   ALLOCATE ( RWORK( 150*NMAX+2*MAXRHS ), STAT = AllocateStatus )
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
   EPS = SLAMCH( 'Underflow threshold' )
   WRITE( NOUT, FMT = 9991 )'underflow', EPS
   EPS = SLAMCH( 'Overflow threshold' )
   WRITE( NOUT, FMT = 9991 )'overflow ', EPS
   EPS = SLAMCH( 'Epsilon' )
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
   IF( .NOT.LSAME( C1, 'Complex precision' ) ) THEN
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
         CALL CCHKGE( DOTYPE, NM, MVAL, NN, NVAL, NNB2, NBVAL2, NNS, &
                      NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), &
                      A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), &
                      B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKGE : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVGE( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, &
                      A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, &
                      RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVGE : ', &
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
         CALL CCHKGB( DOTYPE, NM, MVAL, NN, NVAL, NNB2, NBVAL2, NNS, &
                      NSVAL, THRESH, TSTERR, A( 1, 1 ), LA, &
                      A( 1, 3 ), LAFAC, B( 1, 1 ), B( 1, 2 ), &
                      B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKGB : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVGB( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                      A( 1, 1 ), LA, A( 1, 3 ), LAFAC, A( 1, 6 ), &
                      B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, &
                      WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVGB : ', &
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
         CALL CCHKGT( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), &
                      B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKGT : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVGT( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                      A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), &
                      B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVGT : ', &
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
         CALL CCHKPO( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                      THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                      A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                      WORK, RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKPO : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVPO( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, &
                      A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, &
                      RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVPO : ', &
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
         CALL CCHKPS( DOTYPE, NN, NVAL, NNB2, NBVAL2, NRANK, &
                      RANKVAL, THRESH, TSTERR, LDA, A( 1, 1 ), &
                      A( 1, 2 ), A( 1, 3 ), PIV, WORK, RWORK, &
                      NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKPS : ', &
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
         CALL CCHKPP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                      B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, &
                      NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKPP : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVPP( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, &
                      A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, &
                      RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVPP : ', &
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
         CALL CCHKPB( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                      THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                      A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                      WORK, RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKPB : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVPB( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, &
                      A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, &
                      RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVPB : ', &
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
         CALL CCHKPT( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      A( 1, 1 ), S, A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), &
                      B( 1, 3 ), WORK, RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKPT : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVPT( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                      A( 1, 1 ), S, A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), &
                      B( 1, 3 ), WORK, RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVPT : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'HE' ) ) THEN
!
!        HE:  Hermitian indefinite matrices,
!             with partial (Bunch-Kaufman) pivoting algorithm
!
      NTYPES = 10
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CCHKHE( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                      THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                      A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                      WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKHE : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVHE( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, &
                      A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, &
                      NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVHE : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'HR' ) ) THEN
!
!        HR:  Hermitian indefinite matrices,
!             with bounded Bunch-Kaufman (rook) pivoting algorithm
!
      NTYPES = 10
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CCHKHE_ROOK(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                          THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                          A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                          WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKHE_ROOK : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVHE_ROOK( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                           LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                           B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, &
                           RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVHE_ROOK : ', &
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
         CALL CCHKHE_RK( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                         THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                         E, A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), &
                         B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKHE_RK : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVHE_RK( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                         LDA, A( 1, 1 ), A( 1, 2 ), E, A( 1, 3 ), &
                         B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, &
                         RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVHE_RK : ', &
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
         CALL CCHKHE_AA( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, &
                         NSVAL, THRESH, TSTERR, LDA, &
                         A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                         B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                         WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKHE_AA : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVHE_AA( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                         LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                              B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                         WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVHE_AA : ', &
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
         CALL CCHKHE_AA_2STAGE( DOTYPE, NN, NVAL, NNB2, NBVAL2, &
                            NNS, NSVAL, THRESH, TSTERR, LDA, &
                            A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                            B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                            WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKHE_AA_2STAGE : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVHE_AA_2STAGE( &
                            DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                            LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                                 B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                            WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVHE_AA_2STAGE : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9988 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'HP' ) ) THEN
!
!        HP:  Hermitian indefinite packed matrices,
!             with partial (Bunch-Kaufman) pivoting algorithm
!
      NTYPES = 10
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CCHKHP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                      B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, &
                      IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKHP : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVHP( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, &
                      A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, &
                      NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVHP : ', &
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
         CALL CCHKSY( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                      THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                      A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                      WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKSY : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVSY( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, &
                      A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, &
                      NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVSY : ', &
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
         CALL CCHKSY_ROOK(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                          THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                          A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                          WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKSY_ROOK : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVSY_ROOK( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                           LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                           B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, &
                           RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVSY_ROOK : ', &
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
         CALL CCHKSY_RK( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                         THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                         E, A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), &
                         B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKSY_RK : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVSY_RK( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                         LDA, A( 1, 1 ), A( 1, 2 ), E, A( 1, 3 ), &
                         B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, &
                         RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVSY_RK : ', &
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
         CALL CCHKSY_AA( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                         THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                         A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), &
                         B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKSY_AA : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVSY_AA( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                         LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                         B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, &
                         RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVSY_AA : ', &
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
         CALL CCHKSY_AA_2STAGE( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, &
                         NSVAL, THRESH, TSTERR, LDA, &
                         A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                         B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                         WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKSY_AA_2STAGE : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVSY_AA_2STAGE( &
                         DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, &
                         LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                         B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, &
                         RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVSY_AA_2STAGE : ', &
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
         CALL CCHKSP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                      B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, &
                      IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKSP : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
      IF( TSTDRV ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CDRVSP( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, &
                      A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, &
                      NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVSP : ', &
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
         CALL CCHKTR( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, &
                      THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), &
                      B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, &
                      NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKTR : ', &
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
         CALL CCHKTP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKTP : ', &
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
         CALL CCHKTB( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), &
                      B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKTB : ', &
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
         CALL CCHKQR( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, &
                      NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), &
                      A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), &
                      B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), &
                      WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKQR : ', &
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
         CALL CCHKLQ( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, &
                      NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), &
                      A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), &
                      B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), &
                      WORK, RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKLQ : ', &
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
         CALL CCHKQL( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, &
                      NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), &
                      A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), &
                      B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), &
                      WORK, RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKQL : ', &
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
         CALL CCHKRQ( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, &
                      NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), &
                      A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), &
                      B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), &
                      WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKRQ : ', &
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
         CALL CCHKEQ( THREQ, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKEQ : ', &
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
         CALL CCHKTZ( DOTYPE, NM, MVAL, NN, NVAL, THRESH, TSTERR, &
                      A( 1, 1 ), A( 1, 2 ), S( 1 ), &
                      B( 1, 1 ), WORK, RWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKTZ : ', &
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
         CALL CCHKQ3( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, &
                      THRESH, A( 1, 1 ), A( 1, 2 ), S( 1 ), &
                      B( 1, 1 ), WORK, RWORK, IWORK, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKQ3 : ', &
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
         CALL CDRVLS( DOTYPE, NM, MVAL, NN, NVAL, NNS, NSVAL, NNB, &
                      NBVAL, NXVAL, THRESH, TSTERR, A( 1, 1 ), &
                      A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), &
                      S( 1 ), S( NMAX+1 ), NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CDRVLS : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
      ELSE
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
!
   ELSE IF( LSAMEN( 2, C2, 'QT' ) ) THEN
!
!        QT:  QRT routines for general matrices
!
      IF( TSTCHK ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL CCHKQRT( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, &
                       NBVAL, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKQRT : ', &
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
         CALL CCHKQRTP( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, &
                        NBVAL, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKQRTP : ', &
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
         CALL CCHKLQT( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, &
                       NBVAL, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKLQT : ', &
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
         CALL CCHKLQTP( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, &
                        NBVAL, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKLQTP : ', &
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
         CALL CCHKTSQR( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, &
                        NBVAL, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKTSQR : ', &
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
         CALL CCHKUNHR_COL( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, &
                            NBVAL, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CCHKUNHR_COL : ', &
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
   DEALLOCATE (WORK, STAT = AllocateStatus)
   DEALLOCATE (RWORK,  STAT = AllocateStatus)
!
 9999 FORMAT( / ' Execution not attempted due to input errors' )
 9998 FORMAT( / ' End of tests' )
 9997 FORMAT( ' Total time used = ', F16.8, ' seconds', / )
 9996 FORMAT( ' Invalid input value: ', A4, '=', I6, '; must be >=', &
         I6 )
 9995 FORMAT( ' Invalid input value: ', A4, '=', I6, '; must be <=', &
         I6 )
 9994 FORMAT( ' Tests of the COMPLEX LAPACK routines ', &
         / ' LAPACK VERSION ', I1, '.', I1, '.', I1, &
         / / ' The following parameter values will be used:' )
 9993 FORMAT( 4X, A4, ':  ', 10I6, / 11X, 10I6 )
 9992 FORMAT( / ' Routines pass computational tests if test ratio is ', &
         'less than', F8.2, / )
 9991 FORMAT( ' Relative machine ', A, ' is taken to be', E16.6 )
 9990 FORMAT( / 1X, A3, ':  Unrecognized path name' )
 9989 FORMAT( / 1X, A3, ' routines were not tested' )
 9988 FORMAT( / 1X, A3, ' driver routines were not tested' )
!
!     End of CCHKAA
!
   END
