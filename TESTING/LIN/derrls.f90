!> \brief \b DERRLS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DERRLS( PATH, NUNIT )
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
!> DERRLS tests the error exits for the DOUBLE PRECISION least squares
!> driver routines (DGELS, DGELST, DGETSLS, SGELSS, SGELSY, SGELSD).
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
!> \ingroup double_lin
!
!  =====================================================================
   SUBROUTINE DERRLS( PATH, NUNIT )
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
   CHARACTER*2        C2
   INTEGER            INFO, IRNK
   DOUBLE PRECISION   RCOND
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   INTEGER            IP( NMAX )
   DOUBLE PRECISION   A( NMAX, NMAX ), B( NMAX, NMAX ), S( NMAX ), &
                      W( NMAX )
!     ..
!     .. External Functions ..
   LOGICAL            LSAMEN
   EXTERNAL           LSAMEN
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAESM, CHKXER, DGELS, DGELSD, DGELSS, DGELST, &
                      DGELSY, DGETSLS
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
!     .. Executable Statements ..
!
   NOUT = NUNIT
   WRITE( NOUT, FMT = * )
   C2 = PATH( 2: 3 )
   A( 1, 1 ) = 1.0D+0
   A( 1, 2 ) = 2.0D+0
   A( 2, 2 ) = 3.0D+0
   A( 2, 1 ) = 4.0D+0
   OK = .TRUE.
!
   IF( LSAMEN( 2, C2, 'LS' ) ) THEN
!
!        Test error exits for the least squares driver routines.
!
!        DGELS
!
      SRNAMT = 'DGELS '
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELS( '/', 0, 0, 0, A, 1, B, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELS( 'N', -1, 0, 0, A, 1, B, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELS( 'N', 0, -1, 0, A, 1, B, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELS( 'N', 0, 0, -1, A, 1, B, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELS( 'N', 2, 0, 0, A, 1, B, 2, W, 2, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELS( 'N', 2, 0, 0, A, 2, B, 1, W, 2, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELS( 'N', 0, 2, 0, A, 1, B, 1, W, 2, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELS', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELS( 'N', 1, 1, 0, A, 1, B, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
!
!        DGELST
!
      SRNAMT = 'DGELST'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELST( '/', 0, 0, 0, A, 1, B, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELST', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELST( 'N', -1, 0, 0, A, 1, B, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELST', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELST( 'N', 0, -1, 0, A, 1, B, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELST', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELST( 'N', 0, 0, -1, A, 1, B, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELST', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELST( 'N', 2, 0, 0, A, 1, B, 2, W, 2, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELST', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELST( 'N', 2, 0, 0, A, 2, B, 1, W, 2, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELST', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELST( 'N', 0, 2, 0, A, 1, B, 1, W, 2, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELST', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELST( 'N', 1, 1, 0, A, 1, B, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELST', INFOT, NOUT, LERR, OK )
!
!        DGETSLS
!
      SRNAMT = 'DGETSLS'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGETSLS( '/', 0, 0, 0, A, 1, B, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGETSLS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGETSLS', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGETSLS( 'N', -1, 0, 0, A, 1, B, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGETSLS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGETSLS', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGETSLS( 'N', 0, -1, 0, A, 1, B, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGETSLS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGETSLS', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGETSLS( 'N', 0, 0, -1, A, 1, B, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGETSLS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGETSLS', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGETSLS( 'N', 2, 0, 0, A, 1, B, 2, W, 2, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGETSLS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGETSLS', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGETSLS( 'N', 2, 0, 0, A, 2, B, 1, W, 2, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGETSLS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGETSLS', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGETSLS( 'N', 0, 2, 0, A, 1, B, 1, W, 2, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGETSLS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGETSLS', INFOT, NOUT, LERR, OK )
!
!        DGELSS
!
      SRNAMT = 'DGELSS'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELSS( -1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELSS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELSS', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELSS( 0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELSS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELSS', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELSS( 0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELSS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELSS', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELSS( 2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 2, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELSS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELSS', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELSS( 2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 2, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELSS : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELSS', INFOT, NOUT, LERR, OK )
!
!        DGELSY
!
      SRNAMT = 'DGELSY'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELSY( -1, 0, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELSY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELSY', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELSY( 0, -1, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELSY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELSY', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELSY( 0, 0, -1, A, 1, B, 1, IP, RCOND, IRNK, W, 10, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELSY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELSY', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELSY( 2, 0, 0, A, 1, B, 2, IP, RCOND, IRNK, W, 10, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELSY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELSY', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELSY( 2, 0, 0, A, 2, B, 1, IP, RCOND, IRNK, W, 10, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELSY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELSY', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELSY( 2, 2, 1, A, 2, B, 2, IP, RCOND, IRNK, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELSY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELSY', INFOT, NOUT, LERR, OK )
!
!        DGELSD
!
      SRNAMT = 'DGELSD'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELSD( -1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, IP, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELSD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELSD', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELSD( 0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, IP, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELSD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELSD', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELSD( 0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 10, IP, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELSD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELSD', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELSD( 2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 10, IP, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELSD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELSD', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELSD( 2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 10, IP, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELSD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELSD', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGELSD( 2, 2, 1, A, 2, B, 2, S, RCOND, IRNK, W, 1, IP, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGELSD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'DGELSD', INFOT, NOUT, LERR, OK )
   END IF
!
!     Print a summary line.
!
   CALL ALAESM( PATH, OK, NOUT )
!
   RETURN
!
!     End of DERRLS
!
END
                                                                                                                                                                                                                                                                                                            




