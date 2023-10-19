!> \brief \b ZERRTSQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZERRTSQR( PATH, NUNIT )
!
!       .. Scalar Arguments ..
!       CHARACTER*3        PATH
!       INTEGER            NUNIT
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZERRTSQR tests the error exits for the ZOUBLE PRECISION routines
!> that use the TSQR decomposition of a general matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] PATH
!> \verbatim
!>          PATH is CHARACTER*3
!>          The LAPACK path name for the routines to be tested.
!> \endverbatim
!>
!> \param[in] NUNIT
!> \verbatim
!>          NUNIT is INTEGER
!>          The unit number for output.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Zenver
!> \author NAG Ltd.
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
!
!> \ingroup double_lin
!
!  =====================================================================
   SUBROUTINE ZERRTSQR( PATH, NUNIT )
   IMPLICIT NONE
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER*3        PATH
   INTEGER            NUNIT
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            NMAX
   PARAMETER          ( NMAX = 2 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, INFO, J, MB, NB
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   COMPLEX*16         A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), &
                      C( NMAX, NMAX ), TAU(NMAX)
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAESM, CHKXER, ZGEQR, &
                      ZGEMQR, ZGELQ, ZGEMLQ
!     ..
!     .. Scalars in Common ..
   LOGICAL            LERR, OK
   CHARACTER*32       SRNAMT
   INTEGER            INFOT, NOUT
!     ..
!     .. Common blocks ..
   COMMON             / INFOC / INFOT, NOUT, OK, LERR
   COMMON             / SRNAMC / SRNAMT
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          DBLE
!     ..
!     .. Executable Statements ..
!
   NOUT = NUNIT
   WRITE( NOUT, FMT = * )
!
!     Set the variables to innocuous values.
!
   DO J = 1, NMAX
      DO I = 1, NMAX
         A( I, J ) = 1.D0 / DBLE( I+J )
         C( I, J ) = 1.D0 / DBLE( I+J )
         T( I, J ) = 1.D0 / DBLE( I+J )
      END DO
      W( J ) = 0.D0
   END DO
   OK = .TRUE.
!
!     Error exits for TS factorization
!
!     ZGEQR
!
   SRNAMT = 'ZGEQR'
   INFOT = 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEQR( -1, 0, A, 1, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEQR', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEQR( 0, -1, A, 1, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEQR', INFOT, NOUT, LERR, OK )
   INFOT = 4
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEQR( 1, 1, A, 0, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEQR', INFOT, NOUT, LERR, OK )
   INFOT = 6
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEQR( 3, 2, A, 3, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEQR', INFOT, NOUT, LERR, OK )
   INFOT = 8
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEQR( 3, 2, A, 3, TAU, 8, W, 0, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEQR', INFOT, NOUT, LERR, OK )
!
!     ZLATSQR
!
   MB = 1
   NB = 1
   SRNAMT = 'ZLATSQR'
   INFOT = 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLATSQR( -1, 0, MB, NB, A, 1, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLATSQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZLATSQR', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLATSQR( 1, 2, MB, NB, A, 1, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLATSQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZLATSQR', INFOT, NOUT, LERR, OK )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLATSQR( 0, -1, MB, NB, A, 1, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLATSQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZLATSQR', INFOT, NOUT, LERR, OK )
   INFOT = 3
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLATSQR( 2, 1, -1, NB, A, 2, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLATSQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZLATSQR', INFOT, NOUT, LERR, OK )
   INFOT = 4
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLATSQR( 2, 1, MB, 2, A, 2, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLATSQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZLATSQR', INFOT, NOUT, LERR, OK )
   INFOT = 6
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLATSQR( 2, 1, MB, NB, A, 1, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLATSQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZLATSQR', INFOT, NOUT, LERR, OK )
   INFOT = 8
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLATSQR( 2, 1, MB, NB, A, 2, TAU, 0, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLATSQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZLATSQR', INFOT, NOUT, LERR, OK )
   INFOT = 10
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLATSQR( 2, 1, MB, NB, A, 2, TAU, 2, W, 0, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLATSQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZLATSQR', INFOT, NOUT, LERR, OK )
!
!     ZGEMQR
!
   TAU(1)=1
   TAU(2)=1
   SRNAMT = 'ZGEMQR'
   NB=1
   INFOT = 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( '/', 'N', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'L', '/', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
   INFOT = 3
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'L', 'N', -1, 0, 0, A, 1, TAU, 1, C, 1, W,1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
   INFOT = 4
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'L', 'N', 0, -1, 0, A, 1, TAU, 1, C, 1, W,1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'L', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'R', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
   INFOT = 7
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'L', 'N', 2, 1, 0, A, 0, TAU, 1, C, 1, W, 1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
   INFOT = 9
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'R', 'N', 2, 2, 1, A, 2, TAU, 0, C, 1, W, 1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
   INFOT = 9
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'L', 'N', 2, 2, 1, A, 2, TAU, 0, C, 1, W, 1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
   INFOT = 11
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'L', 'N', 2, 1, 1, A, 2, TAU, 6, C, 0, W, 1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
   INFOT = 13
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'L', 'N', 2, 2, 1, A, 2, TAU, 6, C, 2, W, 0,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
!
!     ZGELQ
!
   SRNAMT = 'ZGELQ'
   INFOT = 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGELQ( -1, 0, A, 1, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGELQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGELQ', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGELQ( 0, -1, A, 1, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGELQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGELQ', INFOT, NOUT, LERR, OK )
   INFOT = 4
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGELQ( 1, 1, A, 0, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGELQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGELQ', INFOT, NOUT, LERR, OK )
   INFOT = 6
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGELQ( 2, 3, A, 3, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGELQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGELQ', INFOT, NOUT, LERR, OK )
   INFOT = 8
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGELQ( 2, 3, A, 3, TAU, 8, W, 0, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGELQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGELQ', INFOT, NOUT, LERR, OK )
!
!     ZLASWLQ
!
   MB = 1
   NB = 1
   SRNAMT = 'ZLASWLQ'
   INFOT = 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASWLQ( -1, 0, MB, NB, A, 1, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASWLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZLASWLQ', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASWLQ( 2, 1, MB, NB, A, 1, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASWLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZLASWLQ', INFOT, NOUT, LERR, OK )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASWLQ( 0, -1, MB, NB, A, 1, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASWLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZLASWLQ', INFOT, NOUT, LERR, OK )
   INFOT = 3
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASWLQ( 1, 2, -1, NB, A, 1, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASWLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZLASWLQ', INFOT, NOUT, LERR, OK )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASWLQ( 1, 1, 2, NB, A, 1, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASWLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZLASWLQ', INFOT, NOUT, LERR, OK )
   INFOT = 4
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASWLQ( 1, 2, MB, -1, A, 1, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASWLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZLASWLQ', INFOT, NOUT, LERR, OK )
   INFOT = 6
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASWLQ( 1, 2, MB, NB, A, 0, TAU, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASWLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZLASWLQ', INFOT, NOUT, LERR, OK )
   INFOT = 8
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASWLQ( 1, 2, MB, NB, A, 1, TAU, 0, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASWLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZLASWLQ', INFOT, NOUT, LERR, OK )
   INFOT = 10
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASWLQ( 1, 2, MB, NB, A, 1, TAU, 1, W, 0, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASWLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZLASWLQ', INFOT, NOUT, LERR, OK )
!
!     ZGEMLQ
!
   TAU(1)=1
   TAU(2)=1
   SRNAMT = 'ZGEMLQ'
   NB=1
   INFOT = 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( '/', 'N', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'L', '/', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 3
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'L', 'N', -1, 0, 0, A, 1, TAU, 1, C, 1, W,1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 4
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'L', 'N', 0, -1, 0, A, 1, TAU, 1, C, 1, W,1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'L', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'R', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 7
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'L', 'N', 1, 2, 0, A, 0, TAU, 1, C, 1, W, 1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 9
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'R', 'N', 2, 2, 1, A, 1, TAU, 0, C, 1, W, 1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 9
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'L', 'N', 2, 2, 1, A, 1, TAU, 0, C, 1, W, 1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 11
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'L', 'N', 1, 2, 1, A, 1, TAU, 6, C, 0, W, 1,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 13
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'L', 'N', 2, 2, 1, A, 2, TAU, 6, C, 2, W, 0,INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
!
!     Print a summary line.
!
   CALL ALAESM( PATH, OK, NOUT )
!
   RETURN
!
!     End of ZERRTSQR
!
END
                                                                                                                                                                                                                                                                                                            




