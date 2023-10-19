!> \brief \b ZERRLQ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZERRLQ( PATH, NUNIT )
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
!> ZERRLQ tests the error exits for the COMPLEX*16 routines
!> that use the LQ decomposition of a general matrix.
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
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
!
!> \ingroup complex16_lin
!
!  =====================================================================
   SUBROUTINE ZERRLQ( PATH, NUNIT )
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
   INTEGER            I, INFO, J
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   COMPLEX*16         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), &
                      W( NMAX ), X( NMAX )
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAESM, CHKXER, ZGELQ2, ZGELQF, ZUNGL2, &
                      ZUNGLQ, ZUNML2, ZUNMLQ
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
   INTRINSIC          DBLE, DCMPLX
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
         A( I, J ) = DCMPLX( 1.D0 / DBLE( I+J ), &
                     -1.D0 / DBLE( I+J ) )
         AF( I, J ) = DCMPLX( 1.D0 / DBLE( I+J ), &
                      -1.D0 / DBLE( I+J ) )
      ENDDO
      B( J ) = 0.D0
      W( J ) = 0.D0
      X( J ) = 0.D0
   ENDDO
   OK = .TRUE.
!
!     Error exits for LQ factorization
!
!     ZGELQF
!
   SRNAMT = 'ZGELQF'
   INFOT = 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGELQF( -1, 0, A, 1, B, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGELQF : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGELQF', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGELQF( 0, -1, A, 1, B, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGELQF : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGELQF', INFOT, NOUT, LERR, OK )
   INFOT = 4
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGELQF( 2, 1, A, 1, B, W, 2, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGELQF : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGELQF', INFOT, NOUT, LERR, OK )
   INFOT = 7
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGELQF( 2, 1, A, 2, B, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGELQF : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGELQF', INFOT, NOUT, LERR, OK )
!
!     ZGELQ2
!
   SRNAMT = 'ZGELQ2'
   INFOT = 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGELQ2( -1, 0, A, 1, B, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGELQ2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGELQ2', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGELQ2( 0, -1, A, 1, B, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGELQ2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGELQ2', INFOT, NOUT, LERR, OK )
   INFOT = 4
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGELQ2( 2, 1, A, 1, B, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGELQ2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGELQ2', INFOT, NOUT, LERR, OK )
!
!     ZUNGLQ
!
   SRNAMT = 'ZUNGLQ'
   INFOT = 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGLQ( -1, 0, 0, A, 1, X, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGLQ', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGLQ( 0, -1, 0, A, 1, X, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGLQ', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGLQ( 2, 1, 0, A, 2, X, W, 2, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGLQ', INFOT, NOUT, LERR, OK )
   INFOT = 3
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGLQ( 0, 0, -1, A, 1, X, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGLQ', INFOT, NOUT, LERR, OK )
   INFOT = 3
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGLQ( 1, 1, 2, A, 1, X, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGLQ', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGLQ( 2, 2, 0, A, 1, X, W, 2, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGLQ', INFOT, NOUT, LERR, OK )
   INFOT = 8
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGLQ( 2, 2, 0, A, 2, X, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGLQ', INFOT, NOUT, LERR, OK )
!
!     ZUNGL2
!
   SRNAMT = 'ZUNGL2'
   INFOT = 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGL2( -1, 0, 0, A, 1, X, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGL2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGL2', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGL2( 0, -1, 0, A, 1, X, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGL2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGL2', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGL2( 2, 1, 0, A, 2, X, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGL2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGL2', INFOT, NOUT, LERR, OK )
   INFOT = 3
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGL2( 0, 0, -1, A, 1, X, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGL2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGL2', INFOT, NOUT, LERR, OK )
   INFOT = 3
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGL2( 1, 1, 2, A, 1, X, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGL2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGL2', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGL2( 2, 2, 0, A, 1, X, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGL2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGL2', INFOT, NOUT, LERR, OK )
!
!     ZUNMLQ
!
   SRNAMT = 'ZUNMLQ'
   INFOT = 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMLQ( '/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMLQ( 'L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 3
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMLQ( 'L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 4
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMLQ( 'L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMLQ( 'L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMLQ( 'L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMLQ( 'R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 7
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMLQ( 'L', 'N', 2, 0, 2, A, 1, X, AF, 2, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 7
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMLQ( 'R', 'N', 0, 2, 2, A, 1, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 10
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMLQ( 'L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 12
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMLQ( 'L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMLQ', INFOT, NOUT, LERR, OK )
   INFOT = 12
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMLQ( 'R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMLQ', INFOT, NOUT, LERR, OK )
!
!     ZUNML2
!
   SRNAMT = 'ZUNML2'
   INFOT = 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNML2( '/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNML2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNML2', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNML2( 'L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNML2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNML2', INFOT, NOUT, LERR, OK )
   INFOT = 3
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNML2( 'L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNML2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNML2', INFOT, NOUT, LERR, OK )
   INFOT = 4
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNML2( 'L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNML2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNML2', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNML2( 'L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNML2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNML2', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNML2( 'L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNML2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNML2', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNML2( 'R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNML2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNML2', INFOT, NOUT, LERR, OK )
   INFOT = 7
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNML2( 'L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNML2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNML2', INFOT, NOUT, LERR, OK )
   INFOT = 7
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNML2( 'R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNML2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNML2', INFOT, NOUT, LERR, OK )
   INFOT = 10
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNML2( 'L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNML2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNML2', INFOT, NOUT, LERR, OK )
!
!     Print a summary line.
!
   CALL ALAESM( PATH, OK, NOUT )
!
   RETURN
!
!     End of ZERRLQ
!
END
                                                                                                                                                                                                                                                                                                            




