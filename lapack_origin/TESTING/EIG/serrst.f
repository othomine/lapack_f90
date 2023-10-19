*> \brief \b SERRST
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SERRST( PATH, NUNIT )
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
*> SERRST tests the error exits for SSYTRD, SSYTD2, SORGTR, SORMTR, SSPTRD,
*> SOPGTR, SOPMTR, SSTEQR, SSTERF, SSTEBZ, SSTEIN, SPTEQR, SSBTRD,
*> SSYEV, SSYEVX, SSYEVD, SSBEV, SSBEVX, SSBEVD,
*> SSPEV, SSPEVX, SSPEVD, SSTEV, SSTEVX, SSTEVD, and SSTEDC.
*> SSYEVD_2STAGE, SSYEVR_2STAGE, SSYEVX_2STAGE,
*> SSYEV_2STAGE, SSBEV_2STAGE, SSBEVD_2STAGE,
*> SSBEVX_2STAGE, SSYTRD_2STAGE, SSYTRD_SY2SB,
*> SSYTRD_SB2ST
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
*> \ingroup single_eig
*
*  =====================================================================
      SUBROUTINE SERRST( PATH, NUNIT )
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
*     NMAX has to be at least 3 or LIW may be too small
*     .. Parameters ..
      INTEGER            NMAX, LIW, LW
      PARAMETER          ( NMAX = 3, LIW = 12*NMAX, LW = 20*NMAX )
*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      INTEGER            I, INFO, J, M, N, NSPLIT, NT
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
*     ..
*     .. Local Arrays ..
      INTEGER            I1( NMAX ), I2( NMAX ), I3( NMAX ), IW( LIW )
      REAL               A( NMAX, NMAX ), C( NMAX, NMAX ), D( NMAX ),
     $                   E( NMAX ), Q( NMAX, NMAX ), R( NMAX ),
     $                   TAU( NMAX ), W( LW ), X( NMAX ),
     $                   Z( NMAX, NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHKXER, SOPGTR, SOPMTR, SORGTR, SORMTR, SPTEQR,
     $                   SSBEV, SSBEVD, SSBEVX, SSBTRD, SSPEV, SSPEVD,
     $                   SSPEVX, SSPTRD, SSTEBZ, SSTEDC, SSTEIN, SSTEQR,
     $                   SSTERF, SSTEV, SSTEVD, SSTEVR, SSTEVX, SSYEV,
     $                   SSYEVD, SSYEVR, SSYEVX, SSYTRD, SSYTD2,
     $                   SSYEVD_2STAGE, SSYEVR_2STAGE, SSYEVX_2STAGE,
     $                   SSYEV_2STAGE, SSBEV_2STAGE, SSBEVD_2STAGE,
     $                   SSBEVX_2STAGE, SSYTRD_2STAGE, SSYTRD_SY2SB,
     $                   SSYTRD_SB2ST
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
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
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
            A( I, J ) = 1. / REAL( I+J )
   10    CONTINUE
   20 CONTINUE
      DO 30 J = 1, NMAX
         D( J ) = REAL( J )
         E( J ) = 0.0
         I1( J ) = J
         I2( J ) = J
         TAU( J ) = 1.
   30 CONTINUE
      OK = .TRUE.
      NT = 0
*
*     Test error exits for the ST path.
*
      IF( LSAMEN( 2, C2, 'ST' ) ) THEN
*
*        SSYTRD
*
         SRNAMT = 'SSYTRD'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD( '/', 0, A, 1, D, E, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD( 'U', -1, A, 1, D, E, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD( 'U', 2, A, 1, D, E, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD( 'U', 0, A, 1, D, E, TAU, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD', INFOT, NOUT, LERR, OK )
         NT = NT + 4
*
*        SSYTD2
*
         SRNAMT = 'SSYTD2'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTD2( '/', 0, A, 1, D, E, TAU, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTD2 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTD2', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTD2( 'U', -1, A, 1, D, E, TAU, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTD2 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTD2', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTD2( 'U', 2, A, 1, D, E, TAU, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTD2 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTD2', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        SSYTRD_2STAGE
*
         SRNAMT = 'SSYTRD_2STAGE'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_2STAGE( '/', 'U', 0, A, 1, D, E, TAU,
     $                                  C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_2STAGE( 'H', 'U', 0, A, 1, D, E, TAU,
     $                                  C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_2STAGE( 'N', '/', 0, A, 1, D, E, TAU,
     $                                  C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_2STAGE( 'N', 'U', -1, A, 1, D, E, TAU,
     $                                  C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_2STAGE( 'N', 'U', 2, A, 1, D, E, TAU,
     $                                  C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_2STAGE( 'N', 'U', 0, A, 1, D, E, TAU,
     $                                  C, 0, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_2STAGE( 'N', 'U', 0, A, 1, D, E, TAU,
     $                                  C, 1, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_2STAGE', INFOT, NOUT, LERR, OK )
         NT = NT + 7
*
*        SSYTRD_SY2SB
*
         SRNAMT = 'SSYTRD_SY2SB'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SY2SB( '/', 0, 0, A, 1, C, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SY2SB : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SY2SB', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SY2SB( 'U', -1, 0, A, 1, C, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SY2SB : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SY2SB', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SY2SB( 'U', 0, -1, A, 1, C, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SY2SB : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SY2SB', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SY2SB( 'U', 2, 0, A, 1, C, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SY2SB : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SY2SB', INFOT, NOUT, LERR, OK )
         INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SY2SB( 'U', 0, 2, A, 1, C, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SY2SB : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SY2SB', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SY2SB( 'U', 0, 0, A, 1, C, 1, TAU, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SY2SB : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SY2SB', INFOT, NOUT, LERR, OK )
         NT = NT + 6
*
*        SSYTRD_SB2ST
*
         SRNAMT = 'SSYTRD_SB2ST'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SB2ST( '/', 'N', 'U', 0, 0, A, 1, D, E,
     $                                    C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SB2ST : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SB2ST', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SB2ST( 'Y', '/', 'U', 0, 0, A, 1, D, E,
     $                                    C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SB2ST : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SB2ST', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SB2ST( 'Y', 'H', 'U', 0, 0, A, 1, D, E,
     $                                    C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SB2ST : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SB2ST', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SB2ST( 'Y', 'N', '/', 0, 0, A, 1, D, E,
     $                                    C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SB2ST : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SB2ST', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SB2ST( 'Y', 'N', 'U', -1, 0, A, 1, D, E,
     $                                    C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SB2ST : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SB2ST', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SB2ST( 'Y', 'N', 'U', 0, -1, A, 1, D, E,
     $                                    C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SB2ST : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SB2ST', INFOT, NOUT, LERR, OK )
         INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SB2ST( 'Y', 'N', 'U', 0, 1, A, 1, D, E,
     $                                    C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SB2ST : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SB2ST', INFOT, NOUT, LERR, OK )
         INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SB2ST( 'Y', 'N', 'U', 0, 0, A, 1, D, E,
     $                                    C, 0, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SB2ST : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SB2ST', INFOT, NOUT, LERR, OK )
         INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SB2ST( 'Y', 'N', 'U', 0, 0, A, 1, D, E,
     $                                    C, 1, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SB2ST : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SB2ST', INFOT, NOUT, LERR, OK )
         NT = NT + 9
*
*        SORGTR
*
         SRNAMT = 'SORGTR'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SORGTR( '/', 0, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SORGTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SORGTR', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SORGTR( 'U', -1, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SORGTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SORGTR', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SORGTR( 'U', 2, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SORGTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SORGTR', INFOT, NOUT, LERR, OK )
         INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SORGTR( 'U', 3, A, 3, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SORGTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SORGTR', INFOT, NOUT, LERR, OK )
         NT = NT + 4
*
*        SORMTR
*
         SRNAMT = 'SORMTR'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SORMTR( '/', 'U', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SORMTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SORMTR( 'L', '/', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SORMTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SORMTR( 'L', 'U', '/', 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SORMTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SORMTR( 'L', 'U', 'N', -1, 0, A, 1, TAU, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SORMTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SORMTR( 'L', 'U', 'N', 0, -1, A, 1, TAU, C, 1, W, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SORMTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SORMTR( 'L', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SORMTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SORMTR( 'R', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SORMTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SORMTR( 'L', 'U', 'N', 2, 0, A, 2, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SORMTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SORMTR( 'L', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SORMTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SORMTR( 'R', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SORMTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         NT = NT + 10
*
*        SSPTRD
*
         SRNAMT = 'SSPTRD'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPTRD( '/', 0, A, D, E, TAU, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPTRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPTRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPTRD( 'U', -1, A, D, E, TAU, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPTRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPTRD', INFOT, NOUT, LERR, OK )
         NT = NT + 2
*
*        SOPGTR
*
         SRNAMT = 'SOPGTR'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SOPGTR( '/', 0, A, TAU, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SOPGTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SOPGTR', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SOPGTR( 'U', -1, A, TAU, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SOPGTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SOPGTR', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SOPGTR( 'U', 2, A, TAU, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SOPGTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SOPGTR', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        SOPMTR
*
         SRNAMT = 'SOPMTR'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SOPMTR( '/', 'U', 'N', 0, 0, A, TAU, C, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SOPMTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SOPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SOPMTR( 'L', '/', 'N', 0, 0, A, TAU, C, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SOPMTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SOPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SOPMTR( 'L', 'U', '/', 0, 0, A, TAU, C, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SOPMTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SOPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SOPMTR( 'L', 'U', 'N', -1, 0, A, TAU, C, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SOPMTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SOPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SOPMTR( 'L', 'U', 'N', 0, -1, A, TAU, C, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SOPMTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SOPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SOPMTR( 'L', 'U', 'N', 2, 0, A, TAU, C, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SOPMTR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SOPMTR', INFOT, NOUT, LERR, OK )
         NT = NT + 6
*
*        SPTEQR
*
         SRNAMT = 'SPTEQR'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SPTEQR( '/', 0, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SPTEQR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SPTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SPTEQR( 'N', -1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SPTEQR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SPTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SPTEQR( 'V', 2, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SPTEQR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SPTEQR', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        SSTEBZ
*
         SRNAMT = 'SSTEBZ'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEBZ( '/', 'E', 0, 0.0, 1.0, 1, 0, 0.0, D, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEBZ : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEBZ( 'A', '/', 0, 0.0, 0.0, 0, 0, 0.0, D, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEBZ : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEBZ( 'A', 'E', -1, 0.0, 0.0, 0, 0, 0.0, D, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEBZ : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEBZ( 'V', 'E', 0, 0.0, 0.0, 0, 0, 0.0, D, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEBZ : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEBZ( 'I', 'E', 0, 0.0, 0.0, 0, 0, 0.0, D, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEBZ : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEBZ( 'I', 'E', 1, 0.0, 0.0, 2, 1, 0.0, D, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEBZ : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEBZ( 'I', 'E', 1, 0.0, 0.0, 1, 0, 0.0, D, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEBZ : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEBZ( 'I', 'E', 1, 0.0, 0.0, 1, 2, 0.0, D, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEBZ : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEBZ', INFOT, NOUT, LERR, OK )
         NT = NT + 8
*
*        SSTEIN
*
         SRNAMT = 'SSTEIN'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEIN( -1, D, E, 0, X, I1, I2, Z, 1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEIN : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEIN( 0, D, E, -1, X, I1, I2, Z, 1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEIN : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEIN( 0, D, E, 1, X, I1, I2, Z, 1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEIN : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEIN( 2, D, E, 0, X, I1, I2, Z, 1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEIN : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEIN', INFOT, NOUT, LERR, OK )
         NT = NT + 4
*
*        SSTEQR
*
         SRNAMT = 'SSTEQR'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEQR( '/', 0, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEQR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEQR( 'N', -1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEQR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEQR( 'V', 2, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEQR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEQR', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        SSTERF
*
         SRNAMT = 'SSTERF'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTERF( -1, D, E, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTERF : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTERF', INFOT, NOUT, LERR, OK )
         NT = NT + 1
*
*        SSTEDC
*
         SRNAMT = 'SSTEDC'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEDC( '/', 0, D, E, Z, 1, W, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEDC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEDC( 'N', -1, D, E, Z, 1, W, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEDC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEDC( 'V', 2, D, E, Z, 1, W, 23, IW, 28, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEDC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEDC( 'N', 1, D, E, Z, 1, W, 0, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEDC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEDC( 'I', 2, D, E, Z, 2, W, 0, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEDC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEDC( 'V', 2, D, E, Z, 2, W, 0, IW, 28, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEDC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEDC( 'N', 1, D, E, Z, 1, W, 1, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEDC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEDC( 'I', 2, D, E, Z, 2, W, 19, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEDC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEDC( 'V', 2, D, E, Z, 2, W, 23, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEDC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEDC', INFOT, NOUT, LERR, OK )
         NT = NT + 9
*
*        SSTEVD
*
         SRNAMT = 'SSTEVD'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVD( '/', 0, D, E, Z, 1, W, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVD', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVD( 'N', -1, D, E, Z, 1, W, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVD', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVD( 'V', 2, D, E, Z, 1, W, 19, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVD', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVD( 'N', 1, D, E, Z, 1, W, 0, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVD', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVD( 'V', 2, D, E, Z, 2, W, 12, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVD', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVD( 'N', 0, D, E, Z, 1, W, 1, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVD', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVD( 'V', 2, D, E, Z, 2, W, 19, IW, 11, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVD', INFOT, NOUT, LERR, OK )
         NT = NT + 7
*
*        SSTEV
*
         SRNAMT = 'SSTEV '
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEV( '/', 0, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEV( 'N', -1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEV ', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEV( 'V', 2, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEV ', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        SSTEVX
*
         SRNAMT = 'SSTEVX'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVX( '/', 'A', 0, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVX( 'N', '/', 0, D, E, 0.0, 1.0, 1, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVX( 'N', 'A', -1, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVX( 'N', 'V', 1, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVX( 'N', 'I', 1, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVX( 'N', 'I', 1, D, E, 0.0, 0.0, 2, 1, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVX( 'N', 'I', 2, D, E, 0.0, 0.0, 2, 1, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVX( 'N', 'I', 1, D, E, 0.0, 0.0, 1, 2, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 14
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVX( 'V', 'A', 2, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVX', INFOT, NOUT, LERR, OK )
         NT = NT + 9
*
*        SSTEVR
*
         N = 1
         SRNAMT = 'SSTEVR'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVR( '/', 'A', 0, D, E, 0.0, 0.0, 1, 1, 0.0, M, R, Z,
     $                1, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVR', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVR( 'V', '/', 0, D, E, 0.0, 0.0, 1, 1, 0.0, M, R, Z,
     $                1, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVR', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVR( 'V', 'A', -1, D, E, 0.0, 0.0, 1, 1, 0.0, M, R, Z,
     $                1, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVR', INFOT, NOUT, LERR, OK )
         INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVR( 'V', 'V', 1, D, E, 0.0, 0.0, 1, 1, 0.0, M, R, Z,
     $                1, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVR', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVR( 'V', 'I', 1, D, E, 0.0, 0.0, 0, 1, 0.0, M, W, Z,
     $                1, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVR', INFOT, NOUT, LERR, OK )
         INFOT = 9
         N = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVR( 'V', 'I', 2, D, E, 0.0, 0.0, 2, 1, 0.0, M, W, Z,
     $                1, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVR', INFOT, NOUT, LERR, OK )
         INFOT = 14
         N = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVR( 'V', 'I', 1, D, E, 0.0, 0.0, 1, 1, 0.0, M, W, Z,
     $                0, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVR', INFOT, NOUT, LERR, OK )
         INFOT = 17
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVR( 'V', 'I', 1, D, E, 0.0, 0.0, 1, 1, 0.0, M, W, Z,
     $                1, IW, X, 20*N-1, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVR', INFOT, NOUT, LERR, OK )
         INFOT = 19
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSTEVR( 'V', 'I', 1, D, E, 0.0, 0.0, 1, 1, 0.0, M, W, Z,
     $                1, IW, X, 20*N, IW( 2*N+1 ), 10*N-1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSTEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSTEVR', INFOT, NOUT, LERR, OK )
         NT = NT + 9
*
*        SSYEVD
*
         SRNAMT = 'SSYEVD'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD( '/', 'U', 0, A, 1, X, W, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD( 'N', '/', 0, A, 1, X, W, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD( 'N', 'U', -1, A, 1, X, W, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD( 'N', 'U', 2, A, 1, X, W, 3, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD( 'N', 'U', 1, A, 1, X, W, 0, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD( 'N', 'U', 2, A, 2, X, W, 4, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD( 'V', 'U', 2, A, 2, X, W, 20, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD( 'N', 'U', 1, A, 1, X, W, 1, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD( 'N', 'U', 2, A, 2, X, W, 5, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD( 'V', 'U', 2, A, 2, X, W, 27, IW, 11, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD', INFOT, NOUT, LERR, OK )
         NT = NT + 10
*
*        SSYEVD_2STAGE
*
         SRNAMT = 'SSYEVD_2STAGE'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD_2STAGE( '/', 'U', 0, A, 1, X, W, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD_2STAGE( 'V', 'U', 0, A, 1, X, W, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD_2STAGE( 'N', '/', 0, A, 1, X, W, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD_2STAGE( 'N', 'U', -1, A, 1, X, W, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD_2STAGE( 'N', 'U', 2, A, 1, X, W, 3, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD_2STAGE( 'N', 'U', 1, A, 1, X, W, 0, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD_2STAGE( 'N', 'U', 2, A, 2, X, W, 4, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
*         INFOT = 8
*         CALL SSYEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 20, IW, 12, INFO )
*         CALL CHKXER( 'SSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD_2STAGE( 'N', 'U', 1, A, 1, X, W, 1, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVD_2STAGE( 'N', 'U', 2, A, 2, X, W, 25, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
*         INFOT = 10
*         CALL SSYEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 27, IW, 11, INFO )
*         CALL CHKXER( 'SSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
         NT = NT + 9
*
*        SSYEVR
*
         SRNAMT = 'SSYEVR'
         N = 1
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR( '/', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R,
     $                Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR( 'V', '/', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R,
     $                Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR( 'V', 'A', '/', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M,
     $                R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR( 'V', 'A', 'U', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M,
     $                R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR( 'V', 'A', 'U', 2, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R,
     $                Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR( 'V', 'V', 'U', 1, A, 1, 0.0E0, 0.0E0, 1, 1, 0.0,
     $                M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR( 'V', 'I', 'U', 1, A, 1, 0.0E0, 0.0E0, 0, 1, 0.0,
     $                M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR', INFOT, NOUT, LERR, OK )
         INFOT = 10
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR( 'V', 'I', 'U', 2, A, 2, 0.0E0, 0.0E0, 2, 1, 0.0,
     $                M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR', INFOT, NOUT, LERR, OK )
         INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR( 'V', 'I', 'U', 1, A, 1, 0.0E0, 0.0E0, 1, 1, 0.0,
     $                M, R, Z, 0, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR', INFOT, NOUT, LERR, OK )
         INFOT = 18
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR( 'V', 'I', 'U', 1, A, 1, 0.0E0, 0.0E0, 1, 1, 0.0,
     $                M, R, Z, 1, IW, Q, 26*N-1, IW( 2*N+1 ), 10*N,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR', INFOT, NOUT, LERR, OK )
         INFOT = 20
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR( 'V', 'I', 'U', 1, A, 1, 0.0E0, 0.0E0, 1, 1, 0.0,
     $                M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N-1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR', INFOT, NOUT, LERR, OK )
         NT = NT + 11
*
*        SSYEVR_2STAGE
*
         SRNAMT = 'SSYEVR_2STAGE'
         N = 1
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR_2STAGE( '/', 'A', 'U', 0, A, 1,
     $                0.0E0, 0.0E0, 1, 1, 0.0E0,
     $                M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR_2STAGE( 'V', 'A', 'U', 0, A, 1,
     $                0.0E0, 0.0E0, 1, 1, 0.0E0,
     $                M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR_2STAGE( 'N', '/', 'U', 0, A, 1,
     $                0.0E0, 0.0E0, 1, 1, 0.0E0,
     $                M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR_2STAGE( 'N', 'A', '/', -1, A, 1,
     $                0.0E0, 0.0E0, 1, 1, 0.0E0,
     $                M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR_2STAGE( 'N', 'A', 'U', -1, A, 1,
     $                0.0E0, 0.0E0, 1, 1, 0.0E0,
     $                M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR_2STAGE( 'N', 'A', 'U', 2, A, 1,
     $                0.0E0, 0.0E0, 1, 1, 0.0E0,
     $                M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR_2STAGE( 'N', 'V', 'U', 1, A, 1,
     $                0.0E0, 0.0E0, 1, 1, 0.0E0,
     $                M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR_2STAGE( 'N', 'I', 'U', 1, A, 1,
     $                0.0E0, 0.0E0, 0, 1, 0.0E0,
     $                M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR_2STAGE( 'N', 'I', 'U', 2, A, 2,
     $                0.0E0, 0.0E0, 2, 1, 0.0E0,
     $                M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR_2STAGE( 'N', 'I', 'U', 1, A, 1,
     $                0.0E0, 0.0E0, 1, 1, 0.0E0,
     $                M, R, Z, 0, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 18
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR_2STAGE( 'N', 'I', 'U', 1, A, 1,
     $                0.0E0, 0.0E0, 1, 1, 0.0E0,
     $                M, R, Z, 1, IW, Q, 0, IW( 2*N+1 ), 10*N,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 20
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVR_2STAGE( 'N', 'I', 'U', 1, A, 1,
     $                0.0E0, 0.0E0, 1, 1, 0.0E0,
     $                M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 0,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVR_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVR_2STAGE', INFOT, NOUT, LERR, OK )
         NT = NT + 12
*
*        SSYEV
*
         SRNAMT = 'SSYEV '
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEV( '/', 'U', 0, A, 1, X, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEV( 'N', '/', 0, A, 1, X, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEV( 'N', 'U', -1, A, 1, X, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEV( 'N', 'U', 2, A, 1, X, W, 3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEV( 'N', 'U', 1, A, 1, X, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEV ', INFOT, NOUT, LERR, OK )
         NT = NT + 5
*
*        SSYEV_2STAGE
*
         SRNAMT = 'SSYEV_2STAGE '
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEV_2STAGE( '/', 'U', 0, A, 1, X, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEV_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEV_2STAGE ', INFOT, NOUT, LERR, OK )
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEV_2STAGE( 'V', 'U', 0, A, 1, X, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEV_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEV_2STAGE ', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEV_2STAGE( 'N', '/', 0, A, 1, X, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEV_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEV_2STAGE ', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEV_2STAGE( 'N', 'U', -1, A, 1, X, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEV_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEV_2STAGE ', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEV_2STAGE( 'N', 'U', 2, A, 1, X, W, 3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEV_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEV_2STAGE ', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEV_2STAGE( 'N', 'U', 1, A, 1, X, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEV_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEV_2STAGE ', INFOT, NOUT, LERR, OK )
         NT = NT + 6
*
*        SSYEVX
*
         SRNAMT = 'SSYEVX'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX( '/', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 1, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX( 'N', '/', 'U', 0, A, 1, 0.0, 1.0, 1, 0, 0.0, M, X,
     $                Z, 1, W, 1, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX( 'N', 'A', '/', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 1, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX( 'N', 'A', 'U', -1, A, 1, 0.0, 0.0, 0, 0, 0.0, M,
     $                X, Z, 1, W, 1, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX( 'N', 'A', 'U', 2, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 16, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX( 'N', 'V', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 8, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX( 'N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 8, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX( 'N', 'I', 'U', 1, A, 1, 0.0, 0.0, 2, 1, 0.0, M, X,
     $                Z, 1, W, 8, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX( 'N', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, X,
     $                Z, 1, W, 16, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX( 'N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 2, 0.0, M, X,
     $                Z, 1, W, 8, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX( 'V', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 16, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 17
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX( 'V', 'A', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 0, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX', INFOT, NOUT, LERR, OK )
         NT = NT + 12
*
*        SSYEVX_2STAGE
*
         SRNAMT = 'SSYEVX_2STAGE'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX_2STAGE( '/', 'A', 'U', 0, A, 1,
     $                 0.0E0, 0.0E0, 0, 0, 0.0E0,
     $                M, X, Z, 1, W, 1, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX_2STAGE( 'V', 'A', 'U', 0, A, 1,
     $                 0.0E0, 0.0E0, 0, 0, 0.0E0,
     $                M, X, Z, 1, W, 1, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX_2STAGE( 'N', '/', 'U', 0, A, 1,
     $                0.0E0, 1.0E0, 1, 0, 0.0E0,
     $                M, X, Z, 1, W, 1, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX_2STAGE( 'N', 'A', '/', 0, A, 1,
     $                0.0E0, 0.0E0, 0, 0, 0.0E0,
     $                M, X, Z, 1, W, 1, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX_2STAGE( 'N', 'A', 'U', -1, A, 1,
     $                0.0E0, 0.0E0, 0, 0, 0.0E0,
     $                M, X, Z, 1, W, 1, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX_2STAGE( 'N', 'A', 'U', 2, A, 1,
     $                0.0E0, 0.0E0, 0, 0, 0.0E0,
     $                M, X, Z, 1, W, 16, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX_2STAGE( 'N', 'V', 'U', 1, A, 1,
     $                0.0E0, 0.0E0, 0, 0, 0.0E0,
     $                M, X, Z, 1, W, 8, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX_2STAGE( 'N', 'I', 'U', 1, A, 1,
     $                0.0E0, 0.0E0, 0, 0, 0.0E0,
     $                M, X, Z, 1, W, 8, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX_2STAGE( 'N', 'I', 'U', 1, A, 1,
     $                0.0E0, 0.0E0, 2, 1, 0.0E0,
     $                M, X, Z, 1, W, 8, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX_2STAGE( 'N', 'I', 'U', 2, A, 2,
     $                0.0E0, 0.0E0, 2, 1, 0.0E0,
     $                M, X, Z, 1, W, 16, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX_2STAGE( 'N', 'I', 'U', 1, A, 1,
     $                0.0E0, 0.0E0, 1, 2, 0.0E0,
     $                M, X, Z, 1, W, 8, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX_2STAGE( 'N', 'A', 'U', 2, A, 2,
     $                0.0E0, 0.0E0, 0, 0, 0.0E0,
     $                M, X, Z, 0, W, 16, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 17
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYEVX_2STAGE( 'N', 'A', 'U', 1, A, 1,
     $                0.0E0, 0.0E0, 0, 0, 0.0E0,
     $                M, X, Z, 1, W, 0, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYEVX_2STAGE', INFOT, NOUT, LERR, OK )
         NT = NT + 13
*
*        SSPEVD
*
         SRNAMT = 'SSPEVD'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVD( '/', 'U', 0, A, X, Z, 1, W, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVD', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVD( 'N', '/', 0, A, X, Z, 1, W, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVD', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVD( 'N', 'U', -1, A, X, Z, 1, W, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVD', INFOT, NOUT, LERR, OK )
         INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVD( 'V', 'U', 2, A, X, Z, 1, W, 23, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVD', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVD( 'N', 'U', 1, A, X, Z, 1, W, 0, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVD', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVD( 'N', 'U', 2, A, X, Z, 1, W, 3, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVD', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVD( 'V', 'U', 2, A, X, Z, 2, W, 16, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVD', INFOT, NOUT, LERR, OK )
         INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVD( 'N', 'U', 1, A, X, Z, 1, W, 1, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVD', INFOT, NOUT, LERR, OK )
         INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVD( 'N', 'U', 2, A, X, Z, 1, W, 4, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVD', INFOT, NOUT, LERR, OK )
         INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVD( 'V', 'U', 2, A, X, Z, 2, W, 23, IW, 11, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVD', INFOT, NOUT, LERR, OK )
         NT = NT + 10
*
*        SSPEV
*
         SRNAMT = 'SSPEV '
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEV( '/', 'U', 0, A, W, Z, 1, X, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEV( 'N', '/', 0, A, W, Z, 1, X, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEV ', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEV( 'N', 'U', -1, A, W, Z, 1, X, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEV ', INFOT, NOUT, LERR, OK )
         INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEV( 'V', 'U', 2, A, W, Z, 1, X, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEV ', INFOT, NOUT, LERR, OK )
         NT = NT + 4
*
*        SSPEVX
*
         SRNAMT = 'SSPEVX'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVX( '/', 'A', 'U', 0, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVX( 'N', '/', 'U', 0, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVX( 'N', 'A', '/', 0, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVX( 'N', 'A', 'U', -1, A, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVX', INFOT, NOUT, LERR, OK )
         INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVX( 'N', 'V', 'U', 1, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVX( 'N', 'I', 'U', 1, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVX( 'N', 'I', 'U', 1, A, 0.0, 0.0, 2, 1, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVX( 'N', 'I', 'U', 2, A, 0.0, 0.0, 2, 1, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVX( 'N', 'I', 'U', 1, A, 0.0, 0.0, 1, 2, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVX', INFOT, NOUT, LERR, OK )
         INFOT = 14
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSPEVX( 'V', 'A', 'U', 2, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSPEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSPEVX', INFOT, NOUT, LERR, OK )
         NT = NT + 10
*
*     Test error exits for the SB path.
*
      ELSE IF( LSAMEN( 2, C2, 'SB' ) ) THEN
*
*        SSBTRD
*
         SRNAMT = 'SSBTRD'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBTRD( '/', 'U', 0, 0, A, 1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBTRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBTRD( 'N', '/', 0, 0, A, 1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBTRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBTRD( 'N', 'U', -1, 0, A, 1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBTRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBTRD( 'N', 'U', 0, -1, A, 1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBTRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBTRD( 'N', 'U', 1, 1, A, 1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBTRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBTRD( 'V', 'U', 2, 0, A, 1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBTRD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBTRD', INFOT, NOUT, LERR, OK )
         NT = NT + 6
*
*        SSYTRD_SB2ST
*
         SRNAMT = 'SSYTRD_SB2ST'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SB2ST( '/', 'N', 'U', 0, 0, A, 1, D, E,
     $                                    C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SB2ST : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SB2ST', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SB2ST( 'N', '/', 'U', 0, 0, A, 1, D, E,
     $                                    C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SB2ST : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SB2ST', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SB2ST( 'N', 'H', 'U', 0, 0, A, 1, D, E,
     $                                    C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SB2ST : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SB2ST', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SB2ST( 'N', 'N', '/', 0, 0, A, 1, D, E,
     $                                    C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SB2ST : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SB2ST', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SB2ST( 'N', 'N', 'U', -1, 0, A, 1, D, E,
     $                                    C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SB2ST : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SB2ST', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SB2ST( 'N', 'N', 'U', 0, -1, A, 1, D, E,
     $                                    C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SB2ST : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SB2ST', INFOT, NOUT, LERR, OK )
         INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SB2ST( 'N', 'N', 'U', 0, 1, A, 1, D, E,
     $                                    C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SB2ST : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SB2ST', INFOT, NOUT, LERR, OK )
         INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SB2ST( 'N', 'N', 'U', 0, 0, A, 1, D, E,
     $                                    C, 0, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SB2ST : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SB2ST', INFOT, NOUT, LERR, OK )
         INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSYTRD_SB2ST( 'N', 'N', 'U', 0, 0, A, 1, D, E,
     $                                    C, 1, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSYTRD_SB2ST : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSYTRD_SB2ST', INFOT, NOUT, LERR, OK )
         NT = NT + 9
*
*        SSBEVD
*
         SRNAMT = 'SSBEVD'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD( '/', 'U', 0, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD( 'N', '/', 0, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD( 'N', 'U', -1, 0, A, 1, X, Z, 1, W, 1, IW, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD( 'N', 'U', 0, -1, A, 1, X, Z, 1, W, 1, IW, 1,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD( 'N', 'U', 2, 1, A, 1, X, Z, 1, W, 4, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD( 'V', 'U', 2, 1, A, 2, X, Z, 1, W, 25, IW, 12,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD', INFOT, NOUT, LERR, OK )
         INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD( 'N', 'U', 1, 0, A, 1, X, Z, 1, W, 0, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD', INFOT, NOUT, LERR, OK )
         INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD( 'N', 'U', 2, 0, A, 1, X, Z, 1, W, 3, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD', INFOT, NOUT, LERR, OK )
         INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD( 'V', 'U', 2, 0, A, 1, X, Z, 2, W, 18, IW, 12,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD', INFOT, NOUT, LERR, OK )
         INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD( 'N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD', INFOT, NOUT, LERR, OK )
         INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD( 'V', 'U', 2, 0, A, 1, X, Z, 2, W, 25, IW, 11,
     $                INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD', INFOT, NOUT, LERR, OK )
         NT = NT + 11
*
*        SSBEVD_2STAGE
*
         SRNAMT = 'SSBEVD_2STAGE'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD_2STAGE( '/', 'U', 0, 0, A, 1, X, Z, 1, W,
     $                                        1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD_2STAGE( 'V', 'U', 0, 0, A, 1, X, Z, 1, W,
     $                                        1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD_2STAGE( 'N', '/', 0, 0, A, 1, X, Z, 1, W,
     $                                        1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD_2STAGE( 'N', 'U', -1, 0, A, 1, X, Z, 1, W,
     $                                         1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD_2STAGE( 'N', 'U', 0, -1, A, 1, X, Z, 1, W,
     $                                         1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD_2STAGE( 'N', 'U', 2, 1, A, 1, X, Z, 1, W,
     $                                        4, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
*         INFOT = 9
*         CALL SSBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 1, W,
*     $                                      25, IW, 12, INFO )
*         CALL CHKXER( 'SSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD_2STAGE( 'N', 'U', 1, 0, A, 1, X, Z, 1, W,
     $                                        0, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD_2STAGE( 'N', 'U', 2, 0, A, 1, X, Z, 1, W,
     $                                        3, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
*         INFOT = 11
*         CALL SSBEVD_2STAGE( 'V', 'U', 2, 0, A, 1, X, Z, 2, W,
*     $                                      18, IW, 12, INFO )
*         CALL CHKXER( 'SSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVD_2STAGE( 'N', 'U', 1, 0, A, 1, X, Z, 1, W,
     $                                        1, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVD_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
*         INFOT = 13
*         CALL SSBEVD_2STAGE( 'V', 'U', 2, 0, A, 1, X, Z, 2, W,
*     $                                      25, IW, 11, INFO )
*         CALL CHKXER( 'SSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
*         NT = NT + 12
         NT = NT + 9
*
*        SSBEV
*
         SRNAMT = 'SSBEV '
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEV( '/', 'U', 0, 0, A, 1, X, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEV( 'N', '/', 0, 0, A, 1, X, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEV ', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEV( 'N', 'U', -1, 0, A, 1, X, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEV ', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEV( 'N', 'U', 0, -1, A, 1, X, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEV ', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEV( 'N', 'U', 2, 1, A, 1, X, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEV ', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEV( 'V', 'U', 2, 0, A, 1, X, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEV : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEV ', INFOT, NOUT, LERR, OK )
         NT = NT + 6
*
*        SSBEV_2STAGE
*
         SRNAMT = 'SSBEV_2STAGE '
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEV_2STAGE( '/', 'U', 0, 0, A, 1, X, Z, 1, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEV_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEV_2STAGE ', INFOT, NOUT, LERR, OK )
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEV_2STAGE( 'V', 'U', 0, 0, A, 1, X, Z, 1, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEV_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEV_2STAGE ', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEV_2STAGE( 'N', '/', 0, 0, A, 1, X, Z, 1, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEV_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEV_2STAGE ', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEV_2STAGE( 'N', 'U', -1, 0, A, 1, X, Z, 1, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEV_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEV_2STAGE ', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEV_2STAGE( 'N', 'U', 0, -1, A, 1, X, Z, 1, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEV_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEV_2STAGE ', INFOT, NOUT, LERR, OK )
         INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEV_2STAGE( 'N', 'U', 2, 1, A, 1, X, Z, 1, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEV_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEV_2STAGE ', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEV_2STAGE( 'N', 'U', 2, 0, A, 1, X, Z, 0, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEV_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEV_2STAGE ', INFOT, NOUT, LERR, OK )
         INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEV_2STAGE( 'N', 'U', 0, 0, A, 1, X, Z, 1, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEV_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEV_2STAGE ', INFOT, NOUT, LERR, OK )
         NT = NT + 8
*
*        SSBEVX
*
         SRNAMT = 'SSBEVX'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX( '/', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX( 'N', '/', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX( 'N', 'A', '/', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX( 'N', 'A', 'U', -1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX( 'N', 'A', 'U', 0, -1, A, 1, Q, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX', INFOT, NOUT, LERR, OK )
         INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX( 'N', 'A', 'U', 2, 1, A, 1, Q, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX( 'V', 'A', 'U', 2, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 2, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX', INFOT, NOUT, LERR, OK )
         INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX( 'N', 'V', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX', INFOT, NOUT, LERR, OK )
         INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX( 'N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX', INFOT, NOUT, LERR, OK )
         INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX( 'N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 2, 1,
     $                0.0, M, X, Z, 1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX', INFOT, NOUT, LERR, OK )
         INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX( 'N', 'I', 'U', 2, 0, A, 1, Q, 1, 0.0, 0.0, 2, 1,
     $                0.0, M, X, Z, 1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX', INFOT, NOUT, LERR, OK )
         INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX( 'N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 1, 2,
     $                0.0, M, X, Z, 1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX', INFOT, NOUT, LERR, OK )
         INFOT = 18
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX( 'V', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX', INFOT, NOUT, LERR, OK )
         NT = NT + 13
*
*        SSBEVX_2STAGE
*
         SRNAMT = 'SSBEVX_2STAGE'
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX_2STAGE( '/', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0E0,
     $          0.0E0, 0, 0, 0.0E0, M, X, Z, 1, W, 0, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX_2STAGE( 'V', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0E0,
     $          0.0E0, 0, 0, 0.0E0, M, X, Z, 1, W, 0, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX_2STAGE( 'N', '/', 'U', 0, 0, A, 1, Q, 1, 0.0E0,
     $          0.0E0, 0, 0, 0.0E0, M, X, Z, 1, W, 0, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX_2STAGE( 'N', 'A', '/', 0, 0, A, 1, Q, 1, 0.0E0,
     $          0.0E0, 0, 0, 0.0E0, M, X, Z, 1, W, 0, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX_2STAGE( 'N', 'A', 'U', -1, 0, A, 1, Q, 1, 0.0E0,
     $           0.0E0, 0, 0, 0.0E0, M, X, Z, 1, W, 0, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX_2STAGE( 'N', 'A', 'U', 0, -1, A, 1, Q, 1, 0.0E0,
     $           0.0E0, 0, 0, 0.0E0, M, X, Z, 1, W, 0, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX_2STAGE( 'N', 'A', 'U', 2, 1, A, 1, Q, 1, 0.0E0,
     $          0.0E0, 0, 0, 0.0E0, M, X, Z, 1, W, 0, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
*         INFOT = 9
*         CALL SSBEVX_2STAGE( 'V', 'A', 'U', 2, 0, A, 1, Q, 1, 0.0E0,
*     $          0.0E0, 0, 0, 0.0E0, M, X, Z, 2, W, 0, IW, I3, INFO )
*         CALL CHKXER( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX_2STAGE( 'N', 'V', 'U', 1, 0, A, 1, Q, 1, 0.0E0,
     $          0.0E0, 0, 0, 0.0E0, M, X, Z, 1, W, 0, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX_2STAGE( 'N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0E0,
     $          0.0E0, 0, 0, 0.0E0, M, X, Z, 1, W, 0, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX_2STAGE( 'N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0E0,
     $          0.0E0, 2, 1, 0.0E0, M, X, Z, 1, W, 0, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX_2STAGE( 'N', 'I', 'U', 2, 0, A, 1, Q, 1, 0.0E0,
     $          0.0E0, 2, 1, 0.0E0, M, X, Z, 1, W, 0, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX_2STAGE( 'N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0E0,
     $          0.0E0, 1, 2, 0.0E0, M, X, Z, 1, W, 0, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
*         INFOT = 18
*         CALL SSBEVX_2STAGE( 'V', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0E0,
*     $          0.0E0, 0, 0, 0.0E0, M, X, Z, 1, W, 0, IW, I3, INFO )
*         CALL CHKXER( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 20
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
         CALL SSBEVX_2STAGE( 'N', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0E0,
     $           0.0E0, 0, 0, 0.0E0, M, X, Z, 1, W, 0, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SSBEVX_2STAGE : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
         CALL CHKXER( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
*         NT = NT + 15
         NT = NT + 13
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
*     End of SERRST
*
      END

