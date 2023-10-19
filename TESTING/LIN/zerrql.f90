!> \brief \b ZERRQL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZERRQL( PATH, NUNIT )
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
!> ZERRQL tests the error exits for the COMPLEX*16 routines
!> that use the QL decomposition of a general matrix.
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
   SUBROUTINE ZERRQL( PATH, NUNIT )
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
   EXTERNAL           ALAESM, CHKXER, ZGEQL2, ZGEQLF, ZGEQLS, ZUNG2L, &
                      ZUNGQL, ZUNM2L, ZUNMQL
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
!     Error exits for QL factorization
!
!     ZGEQLF
!
   SRNAMT = 'ZGEQLF'
   INFOT = 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEQLF( -1, 0, A, 1, B, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEQLF : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEQLF', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEQLF( 0, -1, A, 1, B, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEQLF : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEQLF', INFOT, NOUT, LERR, OK )
   INFOT = 4
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEQLF( 2, 1, A, 1, B, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEQLF : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEQLF', INFOT, NOUT, LERR, OK )
   INFOT = 7
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEQLF( 1, 2, A, 1, B, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEQLF : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEQLF', INFOT, NOUT, LERR, OK )
!
!     ZGEQL2
!
   SRNAMT = 'ZGEQL2'
   INFOT = 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEQL2( -1, 0, A, 1, B, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEQL2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEQL2', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEQL2( 0, -1, A, 1, B, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEQL2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEQL2', INFOT, NOUT, LERR, OK )
   INFOT = 4
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEQL2( 2, 1, A, 1, B, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEQL2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZGEQL2', INFOT, NOUT, LERR, OK )
!
!     ZGEQLS
!
   SRNAMT = 'ZGEQLS'
   INFOT = 1
   CALL ZGEQLS( -1, 0, 0, A, 1, X, B, 1, W, 1, INFO )
   CALL CHKXER( 'ZGEQLS', INFOT, NOUT, LERR, OK )
   INFOT = 2
   CALL ZGEQLS( 0, -1, 0, A, 1, X, B, 1, W, 1, INFO )
   CALL CHKXER( 'ZGEQLS', INFOT, NOUT, LERR, OK )
   INFOT = 2
   CALL ZGEQLS( 1, 2, 0, A, 1, X, B, 1, W, 1, INFO )
   CALL CHKXER( 'ZGEQLS', INFOT, NOUT, LERR, OK )
   INFOT = 3
   CALL ZGEQLS( 0, 0, -1, A, 1, X, B, 1, W, 1, INFO )
   CALL CHKXER( 'ZGEQLS', INFOT, NOUT, LERR, OK )
   INFOT = 5
   CALL ZGEQLS( 2, 1, 0, A, 1, X, B, 2, W, 1, INFO )
   CALL CHKXER( 'ZGEQLS', INFOT, NOUT, LERR, OK )
   INFOT = 8
   CALL ZGEQLS( 2, 1, 0, A, 2, X, B, 1, W, 1, INFO )
   CALL CHKXER( 'ZGEQLS', INFOT, NOUT, LERR, OK )
   INFOT = 10
   CALL ZGEQLS( 1, 1, 2, A, 1, X, B, 1, W, 1, INFO )
   CALL CHKXER( 'ZGEQLS', INFOT, NOUT, LERR, OK )
!
!     ZUNGQL
!
   SRNAMT = 'ZUNGQL'
   INFOT = 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGQL( -1, 0, 0, A, 1, X, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGQL', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGQL( 0, -1, 0, A, 1, X, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGQL', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGQL( 1, 2, 0, A, 1, X, W, 2, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGQL', INFOT, NOUT, LERR, OK )
   INFOT = 3
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGQL( 0, 0, -1, A, 1, X, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGQL', INFOT, NOUT, LERR, OK )
   INFOT = 3
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGQL( 1, 1, 2, A, 1, X, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGQL', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGQL( 2, 1, 0, A, 1, X, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGQL', INFOT, NOUT, LERR, OK )
   INFOT = 8
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGQL( 2, 2, 0, A, 2, X, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNGQL', INFOT, NOUT, LERR, OK )
!
!     ZUNG2L
!
   SRNAMT = 'ZUNG2L'
   INFOT = 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNG2L( -1, 0, 0, A, 1, X, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNG2L : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNG2L', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNG2L( 0, -1, 0, A, 1, X, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNG2L : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNG2L', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNG2L( 1, 2, 0, A, 1, X, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNG2L : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNG2L', INFOT, NOUT, LERR, OK )
   INFOT = 3
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNG2L( 0, 0, -1, A, 1, X, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNG2L : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNG2L', INFOT, NOUT, LERR, OK )
   INFOT = 3
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNG2L( 2, 1, 2, A, 2, X, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNG2L : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNG2L', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNG2L( 2, 1, 0, A, 1, X, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNG2L : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNG2L', INFOT, NOUT, LERR, OK )
!
!     ZUNMQL
!
   SRNAMT = 'ZUNMQL'
   INFOT = 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMQL( '/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMQL', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMQL( 'L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMQL', INFOT, NOUT, LERR, OK )
   INFOT = 3
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMQL( 'L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMQL', INFOT, NOUT, LERR, OK )
   INFOT = 4
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMQL( 'L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMQL', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMQL( 'L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMQL', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMQL( 'L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMQL', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMQL( 'R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMQL', INFOT, NOUT, LERR, OK )
   INFOT = 7
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMQL( 'L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMQL', INFOT, NOUT, LERR, OK )
   INFOT = 7
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMQL( 'R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMQL', INFOT, NOUT, LERR, OK )
   INFOT = 10
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMQL( 'L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMQL', INFOT, NOUT, LERR, OK )
   INFOT = 12
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMQL( 'L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMQL', INFOT, NOUT, LERR, OK )
   INFOT = 12
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMQL( 'R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNMQL', INFOT, NOUT, LERR, OK )
!
!     ZUNM2L
!
   SRNAMT = 'ZUNM2L'
   INFOT = 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNM2L( '/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNM2L : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNM2L', INFOT, NOUT, LERR, OK )
   INFOT = 2
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNM2L( 'L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNM2L : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNM2L', INFOT, NOUT, LERR, OK )
   INFOT = 3
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNM2L( 'L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNM2L : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNM2L', INFOT, NOUT, LERR, OK )
   INFOT = 4
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNM2L( 'L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNM2L : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNM2L', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNM2L( 'L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNM2L : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNM2L', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNM2L( 'L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNM2L : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNM2L', INFOT, NOUT, LERR, OK )
   INFOT = 5
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNM2L( 'R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNM2L : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNM2L', INFOT, NOUT, LERR, OK )
   INFOT = 7
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNM2L( 'L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNM2L : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNM2L', INFOT, NOUT, LERR, OK )
   INFOT = 7
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNM2L( 'R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNM2L : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNM2L', INFOT, NOUT, LERR, OK )
   INFOT = 10
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNM2L( 'L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNM2L : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL CHKXER( 'ZUNM2L', INFOT, NOUT, LERR, OK )
!
!     Print a summary line.
!
   CALL ALAESM( PATH, OK, NOUT )
!
   RETURN
!
!     End of ZERRQL
!
END
                                                                                                                                                                                                                                                                                                            




