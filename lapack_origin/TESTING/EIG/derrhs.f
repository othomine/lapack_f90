*> \brief \b DERRHS
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DERRHS( PATH, NUNIT )
*
*       .. Scalar Arguments ..
*       CHARACTER*3        PATH
*       INTEGER            NUNIT
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DERRHS tests the error exits for DGEBAK, DGEBAL, DGEHRD, DGEHD2,
*> DORGHR, DORMHR, DHSEQR, DHSEIN, DTREVC, and DTREVC3.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] PATH
*> \verbatim
*>          PATH is CHARACTER*3
*>          The LAPACK path name for the routines to be tested.
*> \endverbatim
*>
*> \param[in] NUNIT
*> \verbatim
*>          NUNIT is INTEGER
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
*> \ingroup double_eig
*
*  =====================================================================
      SUBROUTINE DERRHS( PATH, NUNIT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      INTEGER            NUNIT
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NMAX, LW
      PARAMETER          ( NMAX = 3, LW = ( NMAX+2 )*( NMAX+2 )+NMAX )
*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      INTEGER            I, IHI, ILO, INFO, J, M, NT
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
*     ..
*     .. Local Arrays ..
      LOGICAL            SEL( NMAX )
      INTEGER            IFAILL( NMAX ), IFAILR( NMAX )
      DOUBLE PRECISION   A( NMAX, NMAX ), C( NMAX, NMAX ), S( NMAX ),
     $                   TAU( NMAX ), VL( NMAX, NMAX ),
     $                   VR( NMAX, NMAX ), W( LW ), WI( NMAX ),
     $                   WR( NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHKXER, DGEBAK, DGEBAL, DGEHRD, DHSEIN, DHSEQR,
     $                   DORGHR, DORMHR, DTREVC, DTREVC3, DGEHD2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
      INTEGER            INFOT, NOUT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )
*
*     Set the variables to innocuous values.
*
      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = 1.D0 / DBLE( I+J )
   10    CONTINUE
         WI( J ) = DBLE( J )
         SEL( J ) = .TRUE.
   20 CONTINUE
      OK = .TRUE.
      NT = 0
*
*     Test error exits of the nonsymmetric eigenvalue routines.
*
      IF( LSAMEN( 2, C2, 'HS' ) ) THEN
*
*        DGEBAL
*
         SRNAMT = 'DGEBAL'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEBAL( '/', 0, A, 1, ILO, IHI, S, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEBAL : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEBAL', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEBAL( 'N', -1, A, 1, ILO, IHI, S, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEBAL : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEBAL', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEBAL( 'N', 2, A, 1, ILO, IHI, S, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEBAL : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEBAL', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        DGEBAK
*
         SRNAMT = 'DGEBAK'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEBAK( '/', 'R', 0, 1, 0, S, 0, A, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEBAK : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEBAK', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEBAK( 'N', '/', 0, 1, 0, S, 0, A, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEBAK : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEBAK', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEBAK( 'N', 'R', -1, 1, 0, S, 0, A, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEBAK : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEBAK', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEBAK( 'N', 'R', 0, 0, 0, S, 0, A, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEBAK : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEBAK', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEBAK( 'N', 'R', 0, 2, 0, S, 0, A, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEBAK : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEBAK', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEBAK( 'N', 'R', 2, 2, 1, S, 0, A, 2, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEBAK : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEBAK', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEBAK( 'N', 'R', 0, 1, 1, S, 0, A, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEBAK : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEBAK', INFOT, NOUT, LERR, OK )
         INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEBAK( 'N', 'R', 0, 1, 0, S, -1, A, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEBAK : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEBAK', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEBAK( 'N', 'R', 2, 1, 2, S, 0, A, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEBAK : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEBAK', INFOT, NOUT, LERR, OK )
         NT = NT + 9
*
*        DGEHRD
*
         SRNAMT = 'DGEHRD'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEHRD( -1, 1, 1, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEHRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEHRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEHRD( 0, 0, 0, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEHRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEHRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEHRD( 0, 2, 0, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEHRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEHRD', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEHRD( 1, 1, 0, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEHRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEHRD', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEHRD( 0, 1, 1, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEHRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEHRD', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEHRD( 2, 1, 1, A, 1, TAU, W, 2, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEHRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEHRD', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEHRD( 2, 1, 2, A, 2, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEHRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEHRD', INFOT, NOUT, LERR, OK )
         NT = NT + 7
*
*        DGEHD2
*
         SRNAMT = 'DGEHD2'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEHD2( -1, 1, 1, A, 1, TAU, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEHD2 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEHD2', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEHD2( 0, 0, 0, A, 1, TAU, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEHD2 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEHD2', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEHD2( 0, 2, 0, A, 1, TAU, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEHD2 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEHD2', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEHD2( 1, 1, 0, A, 1, TAU, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEHD2 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEHD2', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEHD2( 0, 1, 1, A, 1, TAU, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEHD2 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEHD2', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DGEHD2( 2, 1, 1, A, 1, TAU, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DGEHD2 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DGEHD2', INFOT, NOUT, LERR, OK )
         NT = NT + 6
*
*        DORGHR
*
         SRNAMT = 'DORGHR'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORGHR( -1, 1, 1, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORGHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORGHR', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORGHR( 0, 0, 0, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORGHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORGHR', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORGHR( 0, 2, 0, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORGHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORGHR', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORGHR( 1, 1, 0, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORGHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORGHR', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORGHR( 0, 1, 1, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORGHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORGHR', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORGHR( 2, 1, 1, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORGHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORGHR', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORGHR( 3, 1, 3, A, 3, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORGHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORGHR', INFOT, NOUT, LERR, OK )
         NT = NT + 7
*
*        DORMHR
*
         SRNAMT = 'DORMHR'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORMHR( '/', 'N', 0, 0, 1, 0, A, 1, TAU, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORMHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORMHR', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORMHR( 'L', '/', 0, 0, 1, 0, A, 1, TAU, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORMHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORMHR', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORMHR( 'L', 'N', -1, 0, 1, 0, A, 1, TAU, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORMHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORMHR', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORMHR( 'L', 'N', 0, -1, 1, 0, A, 1, TAU, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORMHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORMHR', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORMHR( 'L', 'N', 0, 0, 0, 0, A, 1, TAU, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORMHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORMHR', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORMHR( 'L', 'N', 0, 0, 2, 0, A, 1, TAU, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORMHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORMHR', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORMHR( 'L', 'N', 1, 2, 2, 1, A, 1, TAU, C, 1, W, 2,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORMHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORMHR', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORMHR( 'R', 'N', 2, 1, 2, 1, A, 1, TAU, C, 2, W, 2,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORMHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORMHR', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORMHR( 'L', 'N', 1, 1, 1, 0, A, 1, TAU, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORMHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORMHR', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORMHR( 'L', 'N', 0, 1, 1, 1, A, 1, TAU, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORMHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORMHR', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORMHR( 'R', 'N', 1, 0, 1, 1, A, 1, TAU, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORMHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORMHR', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORMHR( 'L', 'N', 2, 1, 1, 1, A, 1, TAU, C, 2, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORMHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORMHR', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORMHR( 'R', 'N', 1, 2, 1, 1, A, 1, TAU, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORMHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORMHR', INFOT, NOUT, LERR, OK )
         INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORMHR( 'L', 'N', 2, 1, 1, 1, A, 2, TAU, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORMHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORMHR', INFOT, NOUT, LERR, OK )
         INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORMHR( 'L', 'N', 1, 2, 1, 1, A, 1, TAU, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORMHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORMHR', INFOT, NOUT, LERR, OK )
         INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DORMHR( 'R', 'N', 2, 1, 1, 1, A, 1, TAU, C, 2, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DORMHR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DORMHR', INFOT, NOUT, LERR, OK )
         NT = NT + 16
*
*        DHSEQR
*
         SRNAMT = 'DHSEQR'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DHSEQR( '/', 'N', 0, 1, 0, A, 1, WR, WI, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEQR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DHSEQR', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DHSEQR( 'E', '/', 0, 1, 0, A, 1, WR, WI, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEQR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DHSEQR', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DHSEQR( 'E', 'N', -1, 1, 0, A, 1, WR, WI, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEQR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DHSEQR', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DHSEQR( 'E', 'N', 0, 0, 0, A, 1, WR, WI, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEQR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DHSEQR', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DHSEQR( 'E', 'N', 0, 2, 0, A, 1, WR, WI, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEQR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DHSEQR', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DHSEQR( 'E', 'N', 1, 1, 0, A, 1, WR, WI, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEQR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DHSEQR', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DHSEQR( 'E', 'N', 1, 1, 2, A, 1, WR, WI, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEQR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DHSEQR', INFOT, NOUT, LERR, OK )
         INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DHSEQR( 'E', 'N', 2, 1, 2, A, 1, WR, WI, C, 2, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEQR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DHSEQR', INFOT, NOUT, LERR, OK )
         INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DHSEQR( 'E', 'V', 2, 1, 2, A, 2, WR, WI, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEQR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DHSEQR', INFOT, NOUT, LERR, OK )
         INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DHSEQR( 'E', 'N', 2, 1, 2, A, 2, WR, WI, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEQR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DHSEQR', INFOT, NOUT, LERR, OK )
         NT = NT + 10
*
*        DHSEIN
*
         SRNAMT = 'DHSEIN'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DHSEIN( '/', 'N', 'N', SEL, 0, A, 1, WR, WI, VL, 1, VR, 1,
     $                0, M, W, IFAILL, IFAILR, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEIN : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DHSEIN', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DHSEIN( 'R', '/', 'N', SEL, 0, A, 1, WR, WI, VL, 1, VR, 1,
     $                0, M, W, IFAILL, IFAILR, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEIN : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DHSEIN', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DHSEIN( 'R', 'N', '/', SEL, 0, A, 1, WR, WI, VL, 1, VR, 1,
     $                0, M, W, IFAILL, IFAILR, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEIN : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DHSEIN', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DHSEIN( 'R', 'N', 'N', SEL, -1, A, 1, WR, WI, VL, 1, VR,
     $                1, 0, M, W, IFAILL, IFAILR, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEIN : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DHSEIN', INFOT, NOUT, LERR, OK )
         INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DHSEIN( 'R', 'N', 'N', SEL, 2, A, 1, WR, WI, VL, 1, VR, 2,
     $                4, M, W, IFAILL, IFAILR, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEIN : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DHSEIN', INFOT, NOUT, LERR, OK )
         INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DHSEIN( 'L', 'N', 'N', SEL, 2, A, 2, WR, WI, VL, 1, VR, 1,
     $                4, M, W, IFAILL, IFAILR, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEIN : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DHSEIN', INFOT, NOUT, LERR, OK )
         INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DHSEIN( 'R', 'N', 'N', SEL, 2, A, 2, WR, WI, VL, 1, VR, 1,
     $                4, M, W, IFAILL, IFAILR, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEIN : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DHSEIN', INFOT, NOUT, LERR, OK )
         INFOT = 14
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DHSEIN( 'R', 'N', 'N', SEL, 2, A, 2, WR, WI, VL, 1, VR, 2,
     $                1, M, W, IFAILL, IFAILR, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEIN : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DHSEIN', INFOT, NOUT, LERR, OK )
         NT = NT + 8
*
*        DTREVC
*
         SRNAMT = 'DTREVC'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DTREVC( '/', 'A', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DTREVC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DTREVC', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DTREVC( 'L', '/', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DTREVC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DTREVC', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DTREVC( 'L', 'A', SEL, -1, A, 1, VL, 1, VR, 1, 0, M, W,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DTREVC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DTREVC', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DTREVC( 'L', 'A', SEL, 2, A, 1, VL, 2, VR, 1, 4, M, W,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DTREVC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DTREVC', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DTREVC( 'L', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DTREVC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DTREVC', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DTREVC( 'R', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DTREVC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DTREVC', INFOT, NOUT, LERR, OK )
         INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DTREVC( 'L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 1, M, W,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DTREVC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DTREVC', INFOT, NOUT, LERR, OK )
         NT = NT + 7
*
*        DTREVC3
*
         SRNAMT = 'DTREVC3'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DTREVC3( '/', 'A', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W,
     $                 LW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DTREVC3 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DTREVC3', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DTREVC3( 'L', '/', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W,
     $                 LW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DTREVC3 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DTREVC3', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DTREVC3( 'L', 'A', SEL, -1, A, 1, VL, 1, VR, 1, 0, M, W,
     $                 LW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DTREVC3 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DTREVC3', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DTREVC3( 'L', 'A', SEL, 2, A, 1, VL, 2, VR, 1, 4, M, W,
     $                 LW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DTREVC3 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DTREVC3', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DTREVC3( 'L', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W,
     $                 LW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DTREVC3 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DTREVC3', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DTREVC3( 'R', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W,
     $                 LW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DTREVC3 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DTREVC3', INFOT, NOUT, LERR, OK )
         INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DTREVC3( 'L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 1, M, W,
     $                 LW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DTREVC3 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DTREVC3', INFOT, NOUT, LERR, OK )
         INFOT = 14
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL DTREVC3( 'L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 2, M, W,
     $                 2, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : DTREVC3 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'DTREVC3', INFOT, NOUT, LERR, OK )
         NT = NT + 8
      END IF
*
*     Print a summary line.
*
      IF( OK ) THEN
         WRITE( NOUT, FMT = 9999 )PATH, NT
      ELSE
         WRITE( NOUT, FMT = 9998 )PATH
      END IF
*
 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits',
     $      ' (', I3, ' tests done)' )
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ',
     $      'exits ***' )
*
      RETURN
*
*     End of DERRHS
*
      END

