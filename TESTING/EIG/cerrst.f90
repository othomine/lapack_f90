!> \brief \b CERRST
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CERRST( PATH, NUNIT )
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
!> CERRST tests the error exits for CHETRD, CHETD2, CUNGTR, CUNMTR, CHPTRD,
!> CUNGTR, CUPMTR, CSTEQR, CSTEIN, CPTEQR, CHBTRD,
!> CHEEV, CHEEVX, CHEEVD, CHBEV, CHBEVX, CHBEVD,
!> CHPEV, CHPEVX, CHPEVD, and CSTEDC.
!> CHEEVD_2STAGE, CHEEVR_2STAGE, CHEEVX_2STAGE,
!> CHEEV_2STAGE, CHBEV_2STAGE, CHBEVD_2STAGE,
!> CHBEVX_2STAGE, CHETRD_2STAGE, CHETRD_HE2HB,
!> CHETRD_HB2ST
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
!> \ingroup complex_eig
!
!  =====================================================================
   SUBROUTINE CERRST( PATH, NUNIT )
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
   INTEGER            NMAX, LIW, LW
   PARAMETER          ( NMAX = 3, LIW = 12*NMAX, LW = 20*NMAX )
!     ..
!     .. Local Scalars ..
   CHARACTER*2        C2
   INTEGER            I, INFO, J, M, N, NT
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   INTEGER            I1( NMAX ), I2( NMAX ), I3( NMAX ), IW( LIW )
   REAL               D( NMAX ), E( NMAX ), R( LW ), RW( LW ), &
                      X( NMAX )
   COMPLEX            A( NMAX, NMAX ), C( NMAX, NMAX ), &
                      Q( NMAX, NMAX ), TAU( NMAX ), W( LW ), &
                      Z( NMAX, NMAX )
!     ..
!     .. External Functions ..
   LOGICAL            LSAMEN
   EXTERNAL           LSAMEN
!     ..
!     .. External Subroutines ..
   EXTERNAL           CHBEV, CHBEVD, CHBEVX, CHBTRD, CHEEV, CHEEVD, &
                      CHEEVR, CHEEVX, CHETRD, CHKXER, CHPEV, CHPEVD, &
                      CHPEVX, CHPTRD, CPTEQR, CSTEDC, CSTEIN, CSTEQR, &
                      CUNGTR, CUNMTR, CUPGTR, CUPMTR, CHETD2, &
                      CHEEVD_2STAGE, CHEEVR_2STAGE, CHEEVX_2STAGE, &
                      CHEEV_2STAGE, CHBEV_2STAGE, CHBEVD_2STAGE, &
                      CHBEVX_2STAGE, CHETRD_2STAGE, CHETRD_HE2HB, &
                      CHETRD_HB2ST
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
   INTRINSIC          REAL
!     ..
!     .. Executable Statements ..
!
   NOUT = NUNIT
   WRITE( NOUT, FMT = * )
   C2 = PATH( 2: 3 )
!
!     Set the variables to innocuous values.
!
   DO J = 1, NMAX
      DO I = 1, NMAX
         A( I, J ) = 1. / REAL( I+J )
      ENDDO
   ENDDO
   FORALL (J = 1:NMAX)
      D( J ) = REAL( J )
      I1( J ) = J
      I2( J ) = J
   ENDFORALL
   E(1:NMAX) = 0.0
   TAU(1:NMAX) = 1.
   OK = .TRUE.
   NT = 0
!
!     Test error exits for the ST path.
!
   IF( LSAMEN( 2, C2, 'ST' ) ) THEN
!
!        CHETRD
!
      SRNAMT = 'CHETRD'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD( '/', 0, A, 1, D, E, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD( 'U', -1, A, 1, D, E, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD( 'U', 2, A, 1, D, E, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD( 'U', 0, A, 1, D, E, TAU, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD', INFOT, NOUT, LERR, OK )
      NT = NT + 4
!
!        CHETD2
!
      SRNAMT = 'CHETD2'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETD2( '/', 0, A, 1, D, E, TAU, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETD2 : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETD2', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETD2( 'U', -1, A, 1, D, E, TAU, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETD2 : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETD2', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETD2( 'U', 2, A, 1, D, E, TAU, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETD2 : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETD2', INFOT, NOUT, LERR, OK )
      NT = NT + 3
!
!        CHETRD_2STAGE
!
      SRNAMT = 'CHETRD_2STAGE'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_2STAGE( '/', 'U', 0, A, 1, D, E, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_2STAGE( 'H', 'U', 0, A, 1, D, E, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_2STAGE( 'N', '/', 0, A, 1, D, E, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_2STAGE( 'N', 'U', -1, A, 1, D, E, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_2STAGE( 'N', 'U', 2, A, 1, D, E, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_2STAGE( 'N', 'U', 0, A, 1, D, E, TAU, C, 0, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_2STAGE( 'N', 'U', 0, A, 1, D, E, TAU, C, 1, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_2STAGE', INFOT, NOUT, LERR, OK )
      NT = NT + 7
!
!        CHETRD_HE2HB
!
      SRNAMT = 'CHETRD_HE2HB'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HE2HB( '/', 0, 0, A, 1, C, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HE2HB : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HE2HB', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HE2HB( 'U', -1, 0, A, 1, C, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HE2HB : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HE2HB', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HE2HB( 'U', 0, -1, A, 1, C, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HE2HB : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HE2HB', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HE2HB( 'U', 2, 0, A, 1, C, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HE2HB : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HE2HB', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HE2HB( 'U', 0, 2, A, 1, C, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HE2HB : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HE2HB', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HE2HB( 'U', 0, 0, A, 1, C, 1, TAU, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HE2HB : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HE2HB', INFOT, NOUT, LERR, OK )
      NT = NT + 6
!
!        CHETRD_HB2ST
!
      SRNAMT = 'CHETRD_HB2ST'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HB2ST( '/', 'N', 'U', 0, 0, A, 1, D, E, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HB2ST( 'Y', '/', 'U', 0, 0, A, 1, D, E, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HB2ST( 'Y', 'H', 'U', 0, 0, A, 1, D, E, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HB2ST( 'Y', 'N', '/', 0, 0, A, 1, D, E, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HB2ST( 'Y', 'N', 'U', -1, 0, A, 1, D, E, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HB2ST( 'Y', 'N', 'U', 0, -1, A, 1, D, E, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HB2ST( 'Y', 'N', 'U', 0, 1, A, 1, D, E, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HB2ST( 'Y', 'N', 'U', 0, 0, A, 1, D, E, C, 0, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HB2ST( 'Y', 'N', 'U', 0, 0, A, 1, D, E, C, 1, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      NT = NT + 9
!
!        CUNGTR
!
      SRNAMT = 'CUNGTR'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUNGTR( '/', 0, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUNGTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUNGTR', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUNGTR( 'U', -1, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUNGTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUNGTR', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUNGTR( 'U', 2, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUNGTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUNGTR', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUNGTR( 'U', 3, A, 3, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUNGTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUNGTR', INFOT, NOUT, LERR, OK )
      NT = NT + 4
!
!        CUNMTR
!
      SRNAMT = 'CUNMTR'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUNMTR( '/', 'U', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUNMTR( 'L', '/', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUNMTR( 'L', 'U', '/', 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUNMTR( 'L', 'U', 'N', -1, 0, A, 1, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUNMTR( 'L', 'U', 'N', 0, -1, A, 1, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUNMTR( 'L', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUNMTR( 'R', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUNMTR( 'L', 'U', 'N', 2, 0, A, 2, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUNMTR( 'L', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUNMTR( 'R', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
      NT = NT + 10
!
!        CHPTRD
!
      SRNAMT = 'CHPTRD'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPTRD( '/', 0, A, D, E, TAU, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPTRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPTRD', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPTRD( 'U', -1, A, D, E, TAU, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPTRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPTRD', INFOT, NOUT, LERR, OK )
      NT = NT + 2
!
!        CUPGTR
!
      SRNAMT = 'CUPGTR'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUPGTR( '/', 0, A, TAU, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUPGTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUPGTR', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUPGTR( 'U', -1, A, TAU, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUPGTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUPGTR', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUPGTR( 'U', 2, A, TAU, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUPGTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUPGTR', INFOT, NOUT, LERR, OK )
      NT = NT + 3
!
!        CUPMTR
!
      SRNAMT = 'CUPMTR'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUPMTR( '/', 'U', 'N', 0, 0, A, TAU, C, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUPMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUPMTR', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUPMTR( 'L', '/', 'N', 0, 0, A, TAU, C, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUPMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUPMTR', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUPMTR( 'L', 'U', '/', 0, 0, A, TAU, C, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUPMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUPMTR', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUPMTR( 'L', 'U', 'N', -1, 0, A, TAU, C, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUPMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUPMTR', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUPMTR( 'L', 'U', 'N', 0, -1, A, TAU, C, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUPMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUPMTR', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CUPMTR( 'L', 'U', 'N', 2, 0, A, TAU, C, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CUPMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CUPMTR', INFOT, NOUT, LERR, OK )
      NT = NT + 6
!
!        CPTEQR
!
      SRNAMT = 'CPTEQR'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CPTEQR( '/', 0, D, E, Z, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CPTEQR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CPTEQR', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CPTEQR( 'N', -1, D, E, Z, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CPTEQR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CPTEQR', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CPTEQR( 'V', 2, D, E, Z, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CPTEQR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CPTEQR', INFOT, NOUT, LERR, OK )
      NT = NT + 3
!
!        CSTEIN
!
      SRNAMT = 'CSTEIN'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSTEIN( -1, D, E, 0, X, I1, I2, Z, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSTEIN : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CSTEIN', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSTEIN( 0, D, E, -1, X, I1, I2, Z, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSTEIN : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CSTEIN', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSTEIN( 0, D, E, 1, X, I1, I2, Z, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSTEIN : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CSTEIN', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSTEIN( 2, D, E, 0, X, I1, I2, Z, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSTEIN : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CSTEIN', INFOT, NOUT, LERR, OK )
      NT = NT + 4
!
!        CSTEQR
!
      SRNAMT = 'CSTEQR'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSTEQR( '/', 0, D, E, Z, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSTEQR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CSTEQR', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSTEQR( 'N', -1, D, E, Z, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSTEQR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CSTEQR', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSTEQR( 'V', 2, D, E, Z, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSTEQR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CSTEQR', INFOT, NOUT, LERR, OK )
      NT = NT + 3
!
!        CSTEDC
!
      SRNAMT = 'CSTEDC'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSTEDC( '/', 0, D, E, Z, 1, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSTEDC( 'N', -1, D, E, Z, 1, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSTEDC( 'V', 2, D, E, Z, 1, W, 4, RW, 23, IW, 28, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSTEDC( 'N', 2, D, E, Z, 1, W, 0, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSTEDC( 'V', 2, D, E, Z, 2, W, 0, RW, 23, IW, 28, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSTEDC( 'N', 2, D, E, Z, 1, W, 1, RW, 0, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSTEDC( 'I', 2, D, E, Z, 2, W, 1, RW, 1, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSTEDC( 'V', 2, D, E, Z, 2, W, 4, RW, 1, IW, 28, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSTEDC( 'N', 2, D, E, Z, 1, W, 1, RW, 1, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSTEDC( 'I', 2, D, E, Z, 2, W, 1, RW, 23, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CSTEDC( 'V', 2, D, E, Z, 2, W, 4, RW, 23, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CSTEDC', INFOT, NOUT, LERR, OK )
      NT = NT + 11
!
!        CHEEVD
!
      SRNAMT = 'CHEEVD'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD( '/', 'U', 0, A, 1, X, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD( 'N', '/', 0, A, 1, X, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD( 'N', 'U', -1, A, 1, X, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD( 'N', 'U', 2, A, 1, X, W, 3, RW, 2, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD( 'N', 'U', 1, A, 1, X, W, 0, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD( 'N', 'U', 2, A, 2, X, W, 2, RW, 2, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD( 'V', 'U', 2, A, 2, X, W, 3, RW, 25, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD( 'N', 'U', 1, A, 1, X, W, 1, RW, 0, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD( 'N', 'U', 2, A, 2, X, W, 3, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD( 'V', 'U', 2, A, 2, X, W, 8, RW, 18, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD( 'N', 'U', 1, A, 1, X, W, 1, RW, 1, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD( 'V', 'U', 2, A, 2, X, W, 8, RW, 25, IW, 11, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD', INFOT, NOUT, LERR, OK )
      NT = NT + 12
!
!        CHEEVD_2STAGE
!
      SRNAMT = 'CHEEVD_2STAGE'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD_2STAGE( '/', 'U', 0, A, 1, X, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD_2STAGE( 'V', 'U', 0, A, 1, X, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD_2STAGE( 'N', '/', 0, A, 1, X, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD_2STAGE( 'N', 'U', -1, A, 1, X, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD_2STAGE( 'N', 'U', 2, A, 1, X, W, 3, RW, 2, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD_2STAGE( 'N', 'U', 1, A, 1, X, W, 0, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD_2STAGE( 'N', 'U', 2, A, 2, X, W, 2, RW, 2, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
!         INFOT = 8
!         CALL CHEEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 3,
!     $                            RW, 25, IW, 12, INFO )
!         CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD_2STAGE( 'N', 'U', 1, A, 1, X, W, 1, RW, 0, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD_2STAGE( 'N', 'U', 2, A, 2, X, W, 25, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
!         INFOT = 10
!         CALL CHEEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 8,
!     $                            RW, 18, IW, 12, INFO )
!         CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVD_2STAGE( 'N', 'U', 1, A, 1, X, W, 1, RW, 1, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 12
!         CALL CHEEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 8,
!     $                            RW, 25, IW, 11, INFO )
!         CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      NT = NT + 10
!
!        CHEEV
!
      SRNAMT = 'CHEEV '
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEV( '/', 'U', 0, A, 1, X, W, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEV ', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEV( 'N', '/', 0, A, 1, X, W, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEV ', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEV( 'N', 'U', -1, A, 1, X, W, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEV ', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEV( 'N', 'U', 2, A, 1, X, W, 3, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEV ', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEV( 'N', 'U', 2, A, 2, X, W, 2, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEV ', INFOT, NOUT, LERR, OK )
      NT = NT + 5
!
!        CHEEV_2STAGE
!
      SRNAMT = 'CHEEV_2STAGE '
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEV_2STAGE( '/', 'U', 0, A, 1, X, W, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEV_2STAGE( 'V', 'U', 0, A, 1, X, W, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEV_2STAGE( 'N', '/', 0, A, 1, X, W, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEV_2STAGE( 'N', 'U', -1, A, 1, X, W, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEV_2STAGE( 'N', 'U', 2, A, 1, X, W, 3, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEV_2STAGE( 'N', 'U', 2, A, 2, X, W, 2, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEV_2STAGE ', INFOT, NOUT, LERR, OK )
      NT = NT + 6
!
!        CHEEVX
!
      SRNAMT = 'CHEEVX'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX( '/', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX( 'V', '/', 'U', 0, A, 1, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX( 'V', 'A', '/', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX( 'V', 'A', 'U', -1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX( 'V', 'A', 'U', 2, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, 3, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX( 'V', 'V', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX( 'V', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX( 'V', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 2, W, 3, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX', INFOT, NOUT, LERR, OK )
      INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX( 'V', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 3, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX', INFOT, NOUT, LERR, OK )
      INFOT = 17
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX( 'V', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, 2, RW, IW, I1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX', INFOT, NOUT, LERR, OK )
      NT = NT + 10
!
!        CHEEVX_2STAGE
!
      SRNAMT = 'CHEEVX_2STAGE'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX_2STAGE( '/', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX_2STAGE( 'V', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX_2STAGE( 'N', '/', 'U', 0, A, 1, 0.0, 1.0, 1, 0, 0.0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX_2STAGE( 'N', 'A', '/', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX_2STAGE( 'N', 'A', 'U', -1, A, 1, 0.0, 0.0, 0, 0, 0.0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX_2STAGE( 'N', 'A', 'U', 2, A, 1, 0.0, 0.0, 0, 0, 0.0, &
                   M, X, Z, 2, W, 3, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX_2STAGE( 'N', 'V', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX_2STAGE( 'N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX_2STAGE( 'N', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, &
                   M, X, Z, 2, W, 3, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX_2STAGE( 'N', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, &
                   M, X, Z, 0, W, 3, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 17
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVX_2STAGE( 'N', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, &
                   M, X, Z, 2, W, 0, RW, IW, I1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      NT = NT + 11
!
!        CHEEVR
!
      SRNAMT = 'CHEEVR'
      N = 1
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR( '/', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, &
                   Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR( 'V', '/', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, &
                   Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR( 'V', 'A', '/', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, &
                   R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR( 'V', 'A', 'U', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, &
                   R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR( 'V', 'A', 'U', 2, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, &
                   Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR( 'V', 'V', 'U', 1, A, 1, 0.0E0, 0.0E0, 1, 1, 0.0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR( 'V', 'I', 'U', 1, A, 1, 0.0E0, 0.0E0, 0, 1, 0.0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 10
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR( 'V', 'I', 'U', 2, A, 2, 0.0E0, 0.0E0, 2, 1, 0.0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR( 'V', 'I', 'U', 1, A, 1, 0.0E0, 0.0E0, 1, 1, 0.0, &
                   M, R, Z, 0, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 18
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR( 'V', 'I', 'U', 1, A, 1, 0.0E0, 0.0E0, 1, 1, 0.0, &
                   M, R, Z, 1, IW, Q, 2*N-1, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 20
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR( 'V', 'I', 'U', 1, A, 1, 0.0E0, 0.0E0, 1, 1, 0.0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N-1, IW( 2*N-1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 22
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR( 'V', 'I', 'U', 1, A, 1, 0.0E0, 0.0E0, 1, 1, 0.0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW, 10*N-1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR', INFOT, NOUT, LERR, OK )
      NT = NT + 12
!
!        CHEEVR_2STAGE
!
      SRNAMT = 'CHEEVR_2STAGE'
      N = 1
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR_2STAGE( '/', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR_2STAGE( 'V', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR_2STAGE( 'N', '/', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR_2STAGE( 'N', 'A', '/', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR_2STAGE( 'N', 'A', 'U', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR_2STAGE( 'N', 'A', 'U', 2, A, 1, 0.0, 0.0, 1, 1, 0.0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR_2STAGE( 'N', 'V', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR_2STAGE( 'N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 1, 0.0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR_2STAGE( 'N', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR_2STAGE( 'N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, &
                   M, R, Z, 0, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 18
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR_2STAGE( 'N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, &
                   M, R, Z, 1, IW, Q, 2*N-1, RW, 24*N, IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 20
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR_2STAGE( 'N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, &
                   M, R, Z, 1, IW, Q, 26*N, RW, 24*N-1, IW( 2*N-1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 22
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEEVR_2STAGE( 'N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, &
                   M, R, Z, 1, IW, Q, 26*N, RW, 24*N, IW, 10*N-1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      NT = NT + 13
!
!        CHPEVD
!
      SRNAMT = 'CHPEVD'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVD( '/', 'U', 0, A, X, Z, 1, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVD( 'N', '/', 0, A, X, Z, 1, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVD( 'N', 'U', -1, A, X, Z, 1, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVD( 'V', 'U', 2, A, X, Z, 1, W, 4, RW, 25, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVD( 'N', 'U', 1, A, X, Z, 1, W, 0, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVD( 'N', 'U', 2, A, X, Z, 2, W, 1, RW, 2, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVD( 'V', 'U', 2, A, X, Z, 2, W, 2, RW, 25, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVD( 'N', 'U', 1, A, X, Z, 1, W, 1, RW, 0, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVD( 'N', 'U', 2, A, X, Z, 2, W, 2, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVD( 'V', 'U', 2, A, X, Z, 2, W, 4, RW, 18, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVD( 'N', 'U', 1, A, X, Z, 1, W, 1, RW, 1, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVD( 'N', 'U', 2, A, X, Z, 2, W, 2, RW, 2, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVD( 'V', 'U', 2, A, X, Z, 2, W, 4, RW, 25, IW, 2, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVD', INFOT, NOUT, LERR, OK )
      NT = NT + 13
!
!        CHPEV
!
      SRNAMT = 'CHPEV '
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEV( '/', 'U', 0, A, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEV ', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEV( 'N', '/', 0, A, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEV ', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEV( 'N', 'U', -1, A, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEV ', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEV( 'V', 'U', 2, A, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEV ', INFOT, NOUT, LERR, OK )
      NT = NT + 4
!
!        CHPEVX
!
      SRNAMT = 'CHPEVX'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVX( '/', 'A', 'U', 0, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVX', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVX( 'V', '/', 'U', 0, A, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVX', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVX( 'V', 'A', '/', 0, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVX', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVX( 'V', 'A', 'U', -1, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVX', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVX( 'V', 'V', 'U', 1, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVX', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVX( 'V', 'I', 'U', 1, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVX', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVX( 'V', 'I', 'U', 2, A, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 2, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVX', INFOT, NOUT, LERR, OK )
      INFOT = 14
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHPEVX( 'V', 'A', 'U', 2, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHPEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHPEVX', INFOT, NOUT, LERR, OK )
      NT = NT + 8
!
!     Test error exits for the HB path.
!
   ELSE IF( LSAMEN( 2, C2, 'HB' ) ) THEN
!
!        CHBTRD
!
      SRNAMT = 'CHBTRD'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBTRD( '/', 'U', 0, 0, A, 1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBTRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBTRD', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBTRD( 'N', '/', 0, 0, A, 1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBTRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBTRD', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBTRD( 'N', 'U', -1, 0, A, 1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBTRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBTRD', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBTRD( 'N', 'U', 0, -1, A, 1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBTRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBTRD', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBTRD( 'N', 'U', 1, 1, A, 1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBTRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBTRD', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBTRD( 'V', 'U', 2, 0, A, 1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBTRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBTRD', INFOT, NOUT, LERR, OK )
      NT = NT + 6
!
!        CHETRD_HB2ST
!
      SRNAMT = 'CHETRD_HB2ST'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HB2ST( '/', 'N', 'U', 0, 0, A, 1, D, E, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HB2ST( 'N', '/', 'U', 0, 0, A, 1, D, E, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HB2ST( 'N', 'H', 'U', 0, 0, A, 1, D, E, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HB2ST( 'N', 'N', '/', 0, 0, A, 1, D, E, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HB2ST( 'N', 'N', 'U', -1, 0, A, 1, D, E, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HB2ST( 'N', 'N', 'U', 0, -1, A, 1, D, E, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HB2ST( 'N', 'N', 'U', 0, 1, A, 1, D, E, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HB2ST( 'N', 'N', 'U', 0, 0, A, 1, D, E, C, 0, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHETRD_HB2ST( 'N', 'N', 'U', 0, 0, A, 1, D, E, C, 1, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      NT = NT + 9
!
!        CHBEVD
!
      SRNAMT = 'CHBEVD'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD( '/', 'U', 0, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD( 'N', '/', 0, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD( 'N', 'U', -1, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD( 'N', 'U', 0, -1, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD( 'N', 'U', 2, 1, A, 1, X, Z, 1, W, 2, RW, 2, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD( 'V', 'U', 2, 1, A, 2, X, Z, 1, W, 8, RW, 25, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD( 'N', 'U', 1, 0, A, 1, X, Z, 1, W, 0, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD( 'N', 'U', 2, 1, A, 2, X, Z, 2, W, 1, RW, 2, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD( 'V', 'U', 2, 1, A, 2, X, Z, 2, W, 2, RW, 25, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD( 'N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, RW, 0, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD( 'N', 'U', 2, 1, A, 2, X, Z, 2, W, 2, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD( 'V', 'U', 2, 1, A, 2, X, Z, 2, W, 8, RW, 2, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD( 'N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD( 'N', 'U', 2, 1, A, 2, X, Z, 2, W, 2, RW, 2, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD( 'V', 'U', 2, 1, A, 2, X, Z, 2, W, 8, RW, 25, IW, 2, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD', INFOT, NOUT, LERR, OK )
      NT = NT + 15
!
!        CHBEVD_2STAGE
!
      SRNAMT = 'CHBEVD_2STAGE'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD_2STAGE( '/', 'U', 0, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD_2STAGE( 'V', 'U', 0, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD_2STAGE( 'N', '/', 0, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD_2STAGE( 'N', 'U', -1, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD_2STAGE( 'N', 'U', 0, -1, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD_2STAGE( 'N', 'U', 2, 1, A, 1, X, Z, 1, W, 2, RW, 2, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD_2STAGE( 'N', 'U', 2, 1, A, 2, X, Z, 0, W, 8, RW, 25, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD_2STAGE( 'N', 'U', 1, 0, A, 1, X, Z, 1, W, 0, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD_2STAGE( 'N', 'U', 2, 1, A, 2, X, Z, 2, W, 1, RW, 2, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
!         INFOT = 11
!         CALL CHBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 2,
!     $                         W, 2, RW, 25, IW, 12, INFO )
!         CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD_2STAGE( 'N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, RW, 0, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD_2STAGE( 'N', 'U', 2, 1, A, 2, X, Z, 2, W, 25, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
!         INFOT = 13
!         CALL CHBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 2,
!     $                          W, 25, RW, 2, IW, 12, INFO )
!         CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD_2STAGE( 'N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVD_2STAGE( 'N', 'U', 2, 1, A, 2, X, Z, 2, W, 25, RW, 2, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
!         INFOT = 15
!         CALL CHBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 2,
!     $                          W, 25, RW, 25, IW, 2, INFO )
!         CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      NT = NT + 13
!
!        CHBEV
!
      SRNAMT = 'CHBEV '
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEV( '/', 'U', 0, 0, A, 1, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEV ', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEV( 'N', '/', 0, 0, A, 1, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEV ', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEV( 'N', 'U', -1, 0, A, 1, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEV ', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEV( 'N', 'U', 0, -1, A, 1, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEV ', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEV( 'N', 'U', 2, 1, A, 1, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEV ', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEV( 'V', 'U', 2, 0, A, 1, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEV ', INFOT, NOUT, LERR, OK )
      NT = NT + 6
!
!        CHBEV_2STAGE
!
      SRNAMT = 'CHBEV_2STAGE '
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEV_2STAGE( '/', 'U', 0, 0, A, 1, X, Z, 1, W, 0, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEV_2STAGE( 'V', 'U', 0, 0, A, 1, X, Z, 1, W, 0, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEV_2STAGE( 'N', '/', 0, 0, A, 1, X, Z, 1, W, 0, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEV_2STAGE( 'N', 'U', -1, 0, A, 1, X, Z, 1, W, 0, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEV_2STAGE( 'N', 'U', 0, -1, A, 1, X, Z, 1, W, 0, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEV_2STAGE( 'N', 'U', 2, 1, A, 1, X, Z, 1, W, 0, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEV_2STAGE( 'N', 'U', 2, 0, A, 1, X, Z, 0, W, 0, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEV_2STAGE( 'N', 'U', 2, 0, A, 1, X, Z, 1, W, 0, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEV_2STAGE ', INFOT, NOUT, LERR, OK )
      NT = NT + 8
!
!        CHBEVX
!
      SRNAMT = 'CHBEVX'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX( '/', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX( 'V', '/', 'U', 0, 0, A, 1, Q, 1, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX( 'V', 'A', '/', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX( 'V', 'A', 'U', -1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX( 'V', 'A', 'U', 0, -1, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX( 'V', 'A', 'U', 2, 1, A, 1, Q, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX( 'V', 'A', 'U', 2, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX( 'V', 'V', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX( 'V', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX( 'V', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX', INFOT, NOUT, LERR, OK )
      INFOT = 18
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX( 'V', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX', INFOT, NOUT, LERR, OK )
      NT = NT + 11
!
!        CHBEVX_2STAGE
!
      SRNAMT = 'CHBEVX_2STAGE'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX_2STAGE( '/', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX_2STAGE( 'V', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX_2STAGE( 'N', '/', 'U', 0, 0, A, 1, Q, 1, 0.0, 1.0, 1, 0, 0.0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX_2STAGE( 'N', 'A', '/', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX_2STAGE( 'N', 'A', 'U', -1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX_2STAGE( 'N', 'A', 'U', 0, -1, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX_2STAGE( 'N', 'A', 'U', 2, 1, A, 1, Q, 2, 0.0, 0.0, 0, 0, 0.0, &
                          M, X, Z, 2, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
!         INFOT = 9
!         CALL CHBEVX_2STAGE( 'V', 'A', 'U', 2, 0, A, 1, Q, 1,
!     $                       0.0, 0.0, 0, 0, 0.0,
!     $                       M, X, Z, 2, W, 0, RW, IW, I3, INFO )
!         CALL CHKXER( 'CHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX_2STAGE( 'N', 'V', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX_2STAGE( 'N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX_2STAGE( 'N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 1, 2, 0.0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 18
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX_2STAGE( 'N', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0, 0.0, 0, 0, 0.0, &
                          M, X, Z, 0, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 20
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHBEVX_2STAGE( 'N', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0, 0.0, 0, 0, 0.0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'CHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      NT = NT + 12
   END IF
!
!     Print a summary line.
!
   IF( OK ) THEN
      WRITE( NOUT, FMT = 9999 )PATH, NT
   ELSE
      WRITE( NOUT, FMT = 9998 )PATH
   END IF
!
 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits', ' (', I3, ' tests done)' )
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ', 'exits ***' )
!
   RETURN
!
!     End of CERRST
!
END




