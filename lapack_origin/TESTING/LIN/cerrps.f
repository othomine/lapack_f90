*> \brief \b CERRPS
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CERRPS( PATH, NUNIT )
*
*       .. Scalar Arguments ..
*       INTEGER            NUNIT
*       CHARACTER*3        PATH
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CERRPS tests the error exits for the COMPLEX routines
*> for CPSTRF..
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
*> \ingroup complex_lin
*
*  =====================================================================
      SUBROUTINE CERRPS( PATH, NUNIT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            NUNIT
      CHARACTER*3        PATH
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NMAX
      PARAMETER          ( NMAX = 4 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J, RANK
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
*     ..
*     .. Local Arrays ..
      COMPLEX            A( NMAX, NMAX )
      REAL               RWORK( 2*NMAX )
      INTEGER            PIV( NMAX )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, CPSTF2, CPSTRF
*     ..
*     .. Scalars in Common ..
      INTEGER            INFOT, NOUT
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
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
*
*     Set the variables to innocuous values.
*
      DO 110 J = 1, NMAX
         DO 100 I = 1, NMAX
            A( I, J ) = 1.0 / REAL( I+J )
*
  100    CONTINUE
         PIV( J ) = J
         RWORK( J ) = 0.
         RWORK( NMAX+J ) = 0.
*
  110 CONTINUE
      OK = .TRUE.
*
*
*        Test error exits of the routines that use the Cholesky
*        decomposition of an Hermitian positive semidefinite matrix.
*
*        CPSTRF
*
      SRNAMT = 'CPSTRF'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
      CALL CPSTRF( '/', 0, A, 1, PIV, RANK, -1.0, RWORK, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CPSTRF : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CPSTRF', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
      CALL CPSTRF( 'U', -1, A, 1, PIV, RANK, -1.0, RWORK, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CPSTRF : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CPSTRF', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
      CALL CPSTRF( 'U', 2, A, 1, PIV, RANK, -1.0, RWORK, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CPSTRF : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CPSTRF', INFOT, NOUT, LERR, OK )
*
*        CPSTF2
*
      SRNAMT = 'CPSTF2'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
      CALL CPSTF2( '/', 0, A, 1, PIV, RANK, -1.0, RWORK, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CPSTF2 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CPSTF2', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
      CALL CPSTF2( 'U', -1, A, 1, PIV, RANK, -1.0, RWORK, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CPSTF2 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CPSTF2', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
      CALL CPSTF2( 'U', 2, A, 1, PIV, RANK, -1.0, RWORK, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : CPSTF2 : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CPSTF2', INFOT, NOUT, LERR, OK )
*
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of CERRPS
*
      END

