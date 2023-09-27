!> \brief \b CCHKRFP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       PROGRAM CCHKRFP
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CCHKRFP is the main test program for the COMPLEX linear equation
!> routines with RFP storage format
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \verbatim
!>  MAXIN   INTEGER
!>          The number of different values that can be used for each of
!>          M, N, or NB
!>
!>  MAXRHS  INTEGER
!>          The maximum number of right hand sides
!>
!>  NTYPES  INTEGER
!>
!>  NMAX    INTEGER
!>          The maximum allowable value for N.
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
   PROGRAM CCHKRFP
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            MAXIN
   PARAMETER          ( MAXIN = 12 )
   INTEGER            NMAX
   PARAMETER          ( NMAX =  50 )
   INTEGER            MAXRHS
   PARAMETER          ( MAXRHS = 16 )
   INTEGER            NTYPES
   PARAMETER          ( NTYPES = 9 )
   INTEGER            NIN, NOUT
   PARAMETER          ( NIN = 5, NOUT = 6 )
!     ..
!     .. Local Scalars ..
   LOGICAL            FATAL, TSTERR
   INTEGER            VERS_MAJOR, VERS_MINOR, VERS_PATCH
   INTEGER            I, NN, NNS, NNT
   REAL               EPS, THRESH

   INTEGER(8)         nb_periods_sec, S1, S2, S1T, S2T
   REAL               STOT
!     ..
!     .. Local Arrays ..
   INTEGER            NVAL( MAXIN ), NSVAL( MAXIN ), NTVAL( NTYPES )
   COMPLEX            WORKA( NMAX, NMAX )
   COMPLEX            WORKASAV( NMAX, NMAX )
   COMPLEX            WORKB( NMAX, MAXRHS )
   COMPLEX            WORKXACT( NMAX, MAXRHS )
   COMPLEX            WORKBSAV( NMAX, MAXRHS )
   COMPLEX            WORKX( NMAX, MAXRHS )
   COMPLEX            WORKAFAC( NMAX, NMAX )
   COMPLEX            WORKAINV( NMAX, NMAX )
   COMPLEX            WORKARF( (NMAX*(NMAX+1))/2 )
   COMPLEX            WORKAP( (NMAX*(NMAX+1))/2 )
   COMPLEX            WORKARFINV( (NMAX*(NMAX+1))/2 )
   COMPLEX            C_WORK_CLATMS( 3 * NMAX )
   COMPLEX            C_WORK_CPOT02( NMAX, MAXRHS )
   COMPLEX            C_WORK_CPOT03( NMAX, NMAX )
   REAL               S_WORK_CLATMS( NMAX )
   REAL               S_WORK_CLANHE( NMAX )
   REAL               S_WORK_CPOT01( NMAX )
   REAL               S_WORK_CPOT02( NMAX )
   REAL               S_WORK_CPOT03( NMAX )
!     ..
!     .. External Functions ..
   REAL               SLAMCH, SECOND
   EXTERNAL           SLAMCH, SECOND
!     ..
!     .. External Subroutines ..
   EXTERNAL           ILAVER, CDRVRFP, CDRVRF1, CDRVRF2, CDRVRF3, &
                      CDRVRF4
!     ..
!     .. Executable Statements ..
!
   call system_clock(count_rate=nb_periods_sec,count=S1T)
   FATAL = .FALSE.
!
!     Read a dummy line.
!
   READ( NIN, FMT = * )
!
!     Report LAPACK version tag (e.g. LAPACK-3.2.0)
!
   CALL ILAVER( VERS_MAJOR, VERS_MINOR, VERS_PATCH )
   WRITE( NOUT, FMT = 9994 ) VERS_MAJOR, VERS_MINOR, VERS_PATCH
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
         WRITE( NOUT, FMT = 9996 )' M  ', NVAL( I ), 0
         FATAL = .TRUE.
      ELSE IF( NVAL( I ) > NMAX ) THEN
         WRITE( NOUT, FMT = 9995 )' M  ', NVAL( I ), NMAX
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
!     Read the matrix types
!
   READ( NIN, FMT = * )NNT
   IF( NNT < 1 ) THEN
      WRITE( NOUT, FMT = 9996 )' NMA', NNT, 1
      NNT = 0
      FATAL = .TRUE.
   ELSE IF( NNT > NTYPES ) THEN
      WRITE( NOUT, FMT = 9995 )' NMA', NNT, NTYPES
      NNT = 0
      FATAL = .TRUE.
   END IF
   READ( NIN, FMT = * )( NTVAL( I ), I = 1, NNT )
   DO I = 1, NNT
      IF( NTVAL( I ) < 0 ) THEN
         WRITE( NOUT, FMT = 9996 )'TYPE', NTVAL( I ), 0
         FATAL = .TRUE.
      ELSE IF( NTVAL( I ) > NTYPES ) THEN
         WRITE( NOUT, FMT = 9995 )'TYPE', NTVAL( I ), NTYPES
         FATAL = .TRUE.
      END IF
      ENDDO
   IF( NNT > 0 ) &
      WRITE( NOUT, FMT = 9993 )'TYPE', ( NTVAL( I ), I = 1, NNT )
!
!     Read the threshold value for the test ratios.
!
   READ( NIN, FMT = * )THRESH
   WRITE( NOUT, FMT = 9992 )THRESH
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
!
!     Test the error exit of:
!
   IF( TSTERR ) &
      CALL CERRRFP( NOUT )
!
!    Test the routines: cpftrf, cpftri, cpftrs (as in CDRVPO).
!    This also tests the routines: ctfsm, ctftri, ctfttr, ctrttf.
!
   call system_clock(count_rate=nb_periods_sec,count=S1)
   CALL CDRVRFP( NOUT, NN, NVAL, NNS, NSVAL, NNT, NTVAL, THRESH, &
                 WORKA, WORKASAV, WORKAFAC, WORKAINV, WORKB, &
                 WORKBSAV, WORKXACT, WORKX, WORKARF, WORKARFINV, &
                 C_WORK_CLATMS, C_WORK_CPOT02, &
                 C_WORK_CPOT03, S_WORK_CLATMS, S_WORK_CLANHE, &
                 S_WORK_CPOT01, S_WORK_CPOT02, S_WORK_CPOT03 )
   call system_clock(count_rate=nb_periods_sec,count=S2)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CDRVRFP : ', &
         real(S2-S1)/real(nb_periods_sec), ' s'
   close(10)
!
!    Test the routine: clanhf
!
   call system_clock(count_rate=nb_periods_sec,count=S1)
   CALL CDRVRF1( NOUT, NN, NVAL, THRESH, WORKA, NMAX, WORKARF, &
                 S_WORK_CLANHE )
   call system_clock(count_rate=nb_periods_sec,count=S2)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CDRVRF1 : ', &
         real(S2-S1)/real(nb_periods_sec), ' s'
   close(10)
!
!    Test the conversion routines:
!       chfttp, ctpthf, ctfttr, ctrttf, ctrttp and ctpttr.
!
   call system_clock(count_rate=nb_periods_sec,count=S1)
   CALL CDRVRF2( NOUT, NN, NVAL, WORKA, NMAX, WORKARF, &
                 WORKAP, WORKASAV )
   call system_clock(count_rate=nb_periods_sec,count=S2)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CDRVRF2 : ', &
         real(S2-S1)/real(nb_periods_sec), ' s'
   close(10)
!
!    Test the routine: ctfsm
!
   call system_clock(count_rate=nb_periods_sec,count=S1)
   CALL CDRVRF3( NOUT, NN, NVAL, THRESH, WORKA, NMAX, WORKARF, &
                 WORKAINV, WORKAFAC, S_WORK_CLANHE, &
                 C_WORK_CPOT03, C_WORK_CPOT02 )
   call system_clock(count_rate=nb_periods_sec,count=S2)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CDRVRF3 : ', &
         real(S2-S1)/real(nb_periods_sec), ' s'
   close(10)
!
!
!    Test the routine: chfrk
!
   call system_clock(count_rate=nb_periods_sec,count=S1)
   CALL CDRVRF4( NOUT, NN, NVAL, THRESH, WORKA, WORKAFAC, NMAX, &
                 WORKARF, WORKAINV, NMAX, S_WORK_CLANHE)
   call system_clock(count_rate=nb_periods_sec,count=S2)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CDRVRF4 : ', &
         real(S2-S1)/real(nb_periods_sec), ' s'
   close(10)
!
   CLOSE ( NIN )
   call system_clock(count_rate=nb_periods_sec,count=S2T)
   WRITE( NOUT, FMT = 9998 )
   WRITE( NOUT, FMT = 9997 ) real(S2T - S1T)/real(nb_periods_sec)
!
 9999 FORMAT( / ' Execution not attempted due to input errors' )
 9998 FORMAT( / ' End of tests' )
 9997 FORMAT( ' Total time used = ', F16.8, ' seconds', / )
 9996 FORMAT( ' !! Invalid input value: ', A4, '=', I6, '; must be >=', &
         I6 )
 9995 FORMAT( ' !! Invalid input value: ', A4, '=', I6, '; must be <=', &
         I6 )
 9994 FORMAT( /  ' Tests of the COMPLEX LAPACK RFP routines ', &
         / ' LAPACK VERSION ', I1, '.', I1, '.', I1, &
         / / ' The following parameter values will be used:' )
 9993 FORMAT( 4X, A4, ':  ', 10I6, / 11X, 10I6 )
 9992 FORMAT( / ' Routines pass computational tests if test ratio is ', &
         'less than', F8.2, / )
 9991 FORMAT( ' Relative machine ', A, ' is taken to be', D16.6 )
!
!     End of CCHKRFP
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
