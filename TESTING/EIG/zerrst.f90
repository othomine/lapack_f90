!> \brief \b ZERRST
!
!  @precisions fortran z -> c
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZERRST( PATH, NUNIT )
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
!> ZERRST tests the error exits for ZHETRD, ZHETD2, ZUNGTR, CUNMTR, ZHPTRD,
!> ZUNGTR, ZUPMTR, ZSTEQR, CSTEIN, ZPTEQR, ZHBTRD,
!> ZHEEV, CHEEVX, CHEEVD, ZHBEV, CHBEVX, CHBEVD,
!> ZHPEV, CHPEVX, CHPEVD, and ZSTEDC.
!> ZHEEVD_2STAGE, ZHEEVR_2STAGE, ZHEEVX_2STAGE,
!> ZHEEV_2STAGE, ZHBEV_2STAGE, ZHBEVD_2STAGE,
!> ZHBEVX_2STAGE, ZHETRD_2STAGE
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
!> \ingroup complex16_eig
!
!  =====================================================================
   SUBROUTINE ZERRST( PATH, NUNIT )
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
   DOUBLE PRECISION   D( NMAX ), E( NMAX ), R( LW ), RW( LW ), &
                      X( NMAX )
   COMPLEX*16         A( NMAX, NMAX ), C( NMAX, NMAX ), &
                      Q( NMAX, NMAX ), TAU( NMAX ), W( LW ), &
                      Z( NMAX, NMAX )
!     ..
!     .. External Functions ..
   LOGICAL            LSAMEN
   EXTERNAL           LSAMEN
!     ..
!     .. External Subroutines ..
   EXTERNAL           CHKXER, ZHBEV, ZHBEVD, ZHBEVX, ZHBTRD, ZHEEV, &
                      ZHEEVD, ZHEEVR, ZHEEVX, ZHETRD, ZHPEV, ZHPEVD, &
                      ZHPEVX, ZHPTRD, ZPTEQR, ZSTEDC, ZSTEIN, ZSTEQR, &
                      ZUNGTR, ZUNMTR, ZUPGTR, ZUPMTR, ZHETD2, &
                      ZHEEVD_2STAGE, ZHEEVR_2STAGE, ZHEEVX_2STAGE, &
                      ZHEEV_2STAGE, ZHBEV_2STAGE, ZHBEVD_2STAGE, &
                      ZHBEVX_2STAGE, ZHETRD_2STAGE
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
!
!     Set the variables to innocuous values.
!
   DO J = 1, NMAX
      DO I = 1, NMAX
         A( I, J ) = 1.D0 / DBLE( I+J )
      ENDDO
   ENDDO
   FORALL (J = 1:NMAX)
      D( J ) = DBLE( J )
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
!        ZHETRD
!
      SRNAMT = 'ZHETRD'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD( '/', 0, A, 1, D, E, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD( 'U', -1, A, 1, D, E, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD( 'U', 2, A, 1, D, E, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD( 'U', 0, A, 1, D, E, TAU, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD', INFOT, NOUT, LERR, OK )
      NT = NT + 4
!
!        ZHETD2
!
      SRNAMT = 'ZHETD2'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETD2( '/', 0, A, 1, D, E, TAU, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETD2 : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETD2', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETD2( 'U', -1, A, 1, D, E, TAU, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETD2 : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETD2', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETD2( 'U', 2, A, 1, D, E, TAU, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETD2 : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETD2', INFOT, NOUT, LERR, OK )
      NT = NT + 3
!
!        ZHETRD_2STAGE
!
      SRNAMT = 'ZHETRD_2STAGE'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_2STAGE( '/', 'U', 0, A, 1, D, E, TAU, &
                                     C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_2STAGE( 'H', 'U', 0, A, 1, D, E, TAU, &
                                     C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_2STAGE( 'N', '/', 0, A, 1, D, E, TAU, &
                                     C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_2STAGE( 'N', 'U', -1, A, 1, D, E, TAU, &
                                     C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_2STAGE( 'N', 'U', 2, A, 1, D, E, TAU, &
                                     C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_2STAGE( 'N', 'U', 0, A, 1, D, E, TAU, &
                                     C, 0, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_2STAGE( 'N', 'U', 0, A, 1, D, E, TAU, &
                                     C, 1, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_2STAGE', INFOT, NOUT, LERR, OK )
      NT = NT + 7
!
!        ZHETRD_HE2HB
!
      SRNAMT = 'ZHETRD_HE2HB'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HE2HB( '/', 0, 0, A, 1, C, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HE2HB : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HE2HB', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HE2HB( 'U', -1, 0, A, 1, C, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HE2HB : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HE2HB', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HE2HB( 'U', 0, -1, A, 1, C, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HE2HB : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HE2HB', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HE2HB( 'U', 2, 0, A, 1, C, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HE2HB : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HE2HB', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HE2HB( 'U', 0, 2, A, 1, C, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HE2HB : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HE2HB', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HE2HB( 'U', 0, 0, A, 1, C, 1, TAU, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HE2HB : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HE2HB', INFOT, NOUT, LERR, OK )
      NT = NT + 6
!
!        ZHETRD_HB2ST
!
      SRNAMT = 'ZHETRD_HB2ST'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HB2ST( '/', 'N', 'U', 0, 0, A, 1, D, E, &
                                       C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HB2ST( 'Y', '/', 'U', 0, 0, A, 1, D, E, &
                                       C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HB2ST( 'Y', 'H', 'U', 0, 0, A, 1, D, E, &
                                       C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HB2ST( 'Y', 'N', '/', 0, 0, A, 1, D, E, &
                                       C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HB2ST( 'Y', 'N', 'U', -1, 0, A, 1, D, E, &
                                       C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HB2ST( 'Y', 'N', 'U', 0, -1, A, 1, D, E, &
                                       C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HB2ST( 'Y', 'N', 'U', 0, 1, A, 1, D, E, &
                                       C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HB2ST( 'Y', 'N', 'U', 0, 0, A, 1, D, E, &
                                       C, 0, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HB2ST( 'Y', 'N', 'U', 0, 0, A, 1, D, E, &
                                       C, 1, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      NT = NT + 9
!
!        ZUNGTR
!
      SRNAMT = 'ZUNGTR'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUNGTR( '/', 0, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUNGTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUNGTR', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUNGTR( 'U', -1, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUNGTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUNGTR', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUNGTR( 'U', 2, A, 1, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUNGTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUNGTR', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUNGTR( 'U', 3, A, 3, TAU, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUNGTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUNGTR', INFOT, NOUT, LERR, OK )
      NT = NT + 4
!
!        ZUNMTR
!
      SRNAMT = 'ZUNMTR'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUNMTR( '/', 'U', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUNMTR', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUNMTR( 'L', '/', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUNMTR', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUNMTR( 'L', 'U', '/', 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUNMTR', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUNMTR( 'L', 'U', 'N', -1, 0, A, 1, TAU, C, 1, W, 1, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUNMTR', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUNMTR( 'L', 'U', 'N', 0, -1, A, 1, TAU, C, 1, W, 1, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUNMTR', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUNMTR( 'L', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUNMTR', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUNMTR( 'R', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUNMTR', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUNMTR( 'L', 'U', 'N', 2, 0, A, 2, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUNMTR', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUNMTR( 'L', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUNMTR', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUNMTR( 'R', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUNMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUNMTR', INFOT, NOUT, LERR, OK )
      NT = NT + 10
!
!        ZHPTRD
!
      SRNAMT = 'ZHPTRD'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPTRD( '/', 0, A, D, E, TAU, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPTRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPTRD', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPTRD( 'U', -1, A, D, E, TAU, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPTRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPTRD', INFOT, NOUT, LERR, OK )
      NT = NT + 2
!
!        ZUPGTR
!
      SRNAMT = 'ZUPGTR'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUPGTR( '/', 0, A, TAU, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUPGTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUPGTR', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUPGTR( 'U', -1, A, TAU, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUPGTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUPGTR', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUPGTR( 'U', 2, A, TAU, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUPGTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUPGTR', INFOT, NOUT, LERR, OK )
      NT = NT + 3
!
!        ZUPMTR
!
      SRNAMT = 'ZUPMTR'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUPMTR( '/', 'U', 'N', 0, 0, A, TAU, C, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUPMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUPMTR', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUPMTR( 'L', '/', 'N', 0, 0, A, TAU, C, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUPMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUPMTR', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUPMTR( 'L', 'U', '/', 0, 0, A, TAU, C, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUPMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUPMTR', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUPMTR( 'L', 'U', 'N', -1, 0, A, TAU, C, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUPMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUPMTR', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUPMTR( 'L', 'U', 'N', 0, -1, A, TAU, C, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUPMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUPMTR', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZUPMTR( 'L', 'U', 'N', 2, 0, A, TAU, C, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZUPMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZUPMTR', INFOT, NOUT, LERR, OK )
      NT = NT + 6
!
!        ZPTEQR
!
      SRNAMT = 'ZPTEQR'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZPTEQR( '/', 0, D, E, Z, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZPTEQR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZPTEQR', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZPTEQR( 'N', -1, D, E, Z, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZPTEQR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZPTEQR', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZPTEQR( 'V', 2, D, E, Z, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZPTEQR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZPTEQR', INFOT, NOUT, LERR, OK )
      NT = NT + 3
!
!        ZSTEIN
!
      SRNAMT = 'ZSTEIN'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZSTEIN( -1, D, E, 0, X, I1, I2, Z, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZSTEIN : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZSTEIN', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZSTEIN( 0, D, E, -1, X, I1, I2, Z, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZSTEIN : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZSTEIN', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZSTEIN( 0, D, E, 1, X, I1, I2, Z, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZSTEIN : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZSTEIN', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZSTEIN( 2, D, E, 0, X, I1, I2, Z, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZSTEIN : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZSTEIN', INFOT, NOUT, LERR, OK )
      NT = NT + 4
!
!        ZSTEQR
!
      SRNAMT = 'ZSTEQR'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZSTEQR( '/', 0, D, E, Z, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZSTEQR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZSTEQR', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZSTEQR( 'N', -1, D, E, Z, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZSTEQR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZSTEQR', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZSTEQR( 'V', 2, D, E, Z, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZSTEQR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZSTEQR', INFOT, NOUT, LERR, OK )
      NT = NT + 3
!
!        ZSTEDC
!
      SRNAMT = 'ZSTEDC'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZSTEDC( '/', 0, D, E, Z, 1, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZSTEDC( 'N', -1, D, E, Z, 1, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZSTEDC( 'V', 2, D, E, Z, 1, W, 4, RW, 23, IW, 28, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZSTEDC( 'N', 2, D, E, Z, 1, W, 0, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZSTEDC( 'V', 2, D, E, Z, 2, W, 0, RW, 23, IW, 28, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZSTEDC( 'N', 2, D, E, Z, 1, W, 1, RW, 0, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZSTEDC( 'I', 2, D, E, Z, 2, W, 1, RW, 1, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZSTEDC( 'V', 2, D, E, Z, 2, W, 4, RW, 1, IW, 28, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZSTEDC( 'N', 2, D, E, Z, 1, W, 1, RW, 1, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZSTEDC( 'I', 2, D, E, Z, 2, W, 1, RW, 23, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZSTEDC', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZSTEDC( 'V', 2, D, E, Z, 2, W, 4, RW, 23, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZSTEDC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZSTEDC', INFOT, NOUT, LERR, OK )
      NT = NT + 11
!
!        ZHEEVD
!
      SRNAMT = 'ZHEEVD'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD( '/', 'U', 0, A, 1, X, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD( 'N', '/', 0, A, 1, X, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD( 'N', 'U', -1, A, 1, X, W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD( 'N', 'U', 2, A, 1, X, W, 3, RW, 2, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD( 'N', 'U', 1, A, 1, X, W, 0, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD( 'N', 'U', 2, A, 2, X, W, 2, RW, 2, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD( 'V', 'U', 2, A, 2, X, W, 3, RW, 25, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD( 'N', 'U', 1, A, 1, X, W, 1, RW, 0, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD( 'N', 'U', 2, A, 2, X, W, 3, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD( 'V', 'U', 2, A, 2, X, W, 8, RW, 18, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD( 'N', 'U', 1, A, 1, X, W, 1, RW, 1, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD( 'V', 'U', 2, A, 2, X, W, 8, RW, 25, IW, 11, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD', INFOT, NOUT, LERR, OK )
      NT = NT + 12
!
!        ZHEEVD_2STAGE
!
      SRNAMT = 'ZHEEVD_2STAGE'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD_2STAGE( '/', 'U', 0, A, 1, X, W, 1, &
                                 RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD_2STAGE( 'V', 'U', 0, A, 1, X, W, 1, &
                                 RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD_2STAGE( 'N', '/', 0, A, 1, X, W, 1, &
                                 RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD_2STAGE( 'N', 'U', -1, A, 1, X, W, 1, &
                                  RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD_2STAGE( 'N', 'U', 2, A, 1, X, W, 3, &
                                 RW, 2, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD_2STAGE( 'N', 'U', 1, A, 1, X, W, 0, &
                                 RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD_2STAGE( 'N', 'U', 2, A, 2, X, W, 2, &
                                 RW, 2, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
!         INFOT = 8
!         CALL ZHEEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 3,
!     $                            RW, 25, IW, 12, INFO )
!         CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD_2STAGE( 'N', 'U', 1, A, 1, X, W, 1, &
                                 RW, 0, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD_2STAGE( 'N', 'U', 2, A, 2, X, W, 25, &
                                 RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
!         INFOT = 10
!         CALL ZHEEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 8,
!     $                            RW, 18, IW, 12, INFO )
!         CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVD_2STAGE( 'N', 'U', 1, A, 1, X, W, 1, &
                                 RW, 1, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 12
!         CALL ZHEEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 8,
!     $                            RW, 25, IW, 11, INFO )
!         CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
      NT = NT + 10
!
!        ZHEEV
!
      SRNAMT = 'ZHEEV '
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEV( '/', 'U', 0, A, 1, X, W, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEV ', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEV( 'N', '/', 0, A, 1, X, W, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEV ', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEV( 'N', 'U', -1, A, 1, X, W, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEV ', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEV( 'N', 'U', 2, A, 1, X, W, 3, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEV ', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEV( 'N', 'U', 2, A, 2, X, W, 2, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEV ', INFOT, NOUT, LERR, OK )
      NT = NT + 5
!
!        ZHEEV_2STAGE
!
      SRNAMT = 'ZHEEV_2STAGE '
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEV_2STAGE( '/', 'U', 0, A, 1, X, W, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEV_2STAGE( 'V', 'U', 0, A, 1, X, W, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEV_2STAGE( 'N', '/', 0, A, 1, X, W, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEV_2STAGE( 'N', 'U', -1, A, 1, X, W, 1, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEV_2STAGE( 'N', 'U', 2, A, 1, X, W, 3, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEV_2STAGE( 'N', 'U', 2, A, 2, X, W, 2, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEV_2STAGE ', INFOT, NOUT, LERR, OK )
      NT = NT + 6
!
!        ZHEEVX
!
      SRNAMT = 'ZHEEVX'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX( '/', 'A', 'U', 0, A, 1, 0.0D0, 0.0D0, 0, 0, 0.0D0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX( 'V', '/', 'U', 0, A, 1, 0.0D0, 1.0D0, 1, 0, 0.0D0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX( 'V', 'A', '/', 0, A, 1, 0.0D0, 0.0D0, 0, 0, 0.0D0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX( 'V', 'A', 'U', -1, A, 1, 0.0D0, 0.0D0, 0, 0, &
                   0.0D0, M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX( 'V', 'A', 'U', 2, A, 1, 0.0D0, 0.0D0, 0, 0, 0.0D0, &
                   M, X, Z, 2, W, 3, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX( 'V', 'V', 'U', 1, A, 1, 0.0D0, 0.0D0, 0, 0, 0.0D0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX( 'V', 'I', 'U', 1, A, 1, 0.0D0, 0.0D0, 0, 0, 0.0D0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX( 'V', 'I', 'U', 2, A, 2, 0.0D0, 0.0D0, 2, 1, 0.0D0, &
                   M, X, Z, 2, W, 3, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX', INFOT, NOUT, LERR, OK )
      INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX( 'V', 'A', 'U', 2, A, 2, 0.0D0, 0.0D0, 0, 0, 0.0D0, &
                   M, X, Z, 1, W, 3, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX', INFOT, NOUT, LERR, OK )
      INFOT = 17
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX( 'V', 'A', 'U', 2, A, 2, 0.0D0, 0.0D0, 0, 0, 0.0D0, &
                   M, X, Z, 2, W, 2, RW, IW, I1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX', INFOT, NOUT, LERR, OK )
      NT = NT + 10
!
!        ZHEEVX_2STAGE
!
      SRNAMT = 'ZHEEVX_2STAGE'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX_2STAGE( '/', 'A', 'U', 0, A, 1, &
                   0.0D0, 0.0D0, 0, 0, 0.0D0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX_2STAGE( 'V', 'A', 'U', 0, A, 1, &
                   0.0D0, 0.0D0, 0, 0, 0.0D0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX_2STAGE( 'N', '/', 'U', 0, A, 1, &
                   0.0D0, 1.0D0, 1, 0, 0.0D0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX_2STAGE( 'N', 'A', '/', 0, A, 1, &
                   0.0D0, 0.0D0, 0, 0, 0.0D0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX_2STAGE( 'N', 'A', 'U', -1, A, 1, &
                   0.0D0, 0.0D0, 0, 0, 0.0D0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX_2STAGE( 'N', 'A', 'U', 2, A, 1, &
                   0.0D0, 0.0D0, 0, 0, 0.0D0, &
                   M, X, Z, 2, W, 3, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX_2STAGE( 'N', 'V', 'U', 1, A, 1, &
                   0.0D0, 0.0D0, 0, 0, 0.0D0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX_2STAGE( 'N', 'I', 'U', 1, A, 1, &
                   0.0D0, 0.0D0, 0, 0, 0.0D0, &
                   M, X, Z, 1, W, 1, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX_2STAGE( 'N', 'I', 'U', 2, A, 2, &
                   0.0D0, 0.0D0, 2, 1, 0.0D0, &
                   M, X, Z, 2, W, 3, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX_2STAGE( 'N', 'A', 'U', 2, A, 2, &
                   0.0D0, 0.0D0, 0, 0, 0.0D0, &
                   M, X, Z, 0, W, 3, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 17
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVX_2STAGE( 'N', 'A', 'U', 2, A, 2, &
                   0.0D0, 0.0D0, 0, 0, 0.0D0, &
                   M, X, Z, 2, W, 0, RW, IW, I1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK )
      NT = NT + 11
!
!        ZHEEVR
!
      SRNAMT = 'ZHEEVR'
      N = 1
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR( '/', 'A', 'U', 0, A, 1, 0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR( 'V', '/', 'U', 0, A, 1, 0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR( 'V', 'A', '/', -1, A, 1, 0.0D0, 0.0D0, 1, 1, &
                   0.0D0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, &
                   IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR( 'V', 'A', 'U', -1, A, 1, 0.0D0, 0.0D0, 1, 1, &
                   0.0D0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, &
                   IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR( 'V', 'A', 'U', 2, A, 1, 0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR( 'V', 'V', 'U', 1, A, 1, 0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR( 'V', 'I', 'U', 1, A, 1, 0.0D0, 0.0D0, 0, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 10
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR( 'V', 'I', 'U', 2, A, 2, 0.0D0, 0.0D0, 2, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR( 'V', 'I', 'U', 1, A, 1, 0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 0, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 18
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR( 'V', 'I', 'U', 1, A, 1, 0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N-1, RW, 24*N, IW( 2*N+1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 20
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR( 'V', 'I', 'U', 1, A, 1, 0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N-1, IW( 2*N-1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR', INFOT, NOUT, LERR, OK )
      INFOT = 22
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR( 'V', 'I', 'U', 1, A, 1, 0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW, 10*N-1, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR', INFOT, NOUT, LERR, OK )
      NT = NT + 12
!
!        ZHEEVR_2STAGE
!
      SRNAMT = 'ZHEEVR_2STAGE'
      N = 1
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR_2STAGE( '/', 'A', 'U', 0, A, 1, &
                   0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR_2STAGE( 'V', 'A', 'U', 0, A, 1, &
                   0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR_2STAGE( 'N', '/', 'U', 0, A, 1, &
                   0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR_2STAGE( 'N', 'A', '/', -1, A, 1, &
                   0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, &
                   IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR_2STAGE( 'N', 'A', 'U', -1, A, 1, &
                   0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, &
                   IW( 2*N+1 ), 10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR_2STAGE( 'N', 'A', 'U', 2, A, 1, &
                   0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR_2STAGE( 'N', 'V', 'U', 1, A, 1, &
                   0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR_2STAGE( 'N', 'I', 'U', 1, A, 1, &
                   0.0D0, 0.0D0, 0, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR_2STAGE( 'N', 'I', 'U', 2, A, 2, &
                   0.0D0, 0.0D0, 2, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR_2STAGE( 'N', 'I', 'U', 1, A, 1, &
                   0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 0, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 18
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR_2STAGE( 'N', 'I', 'U', 1, A, 1, &
                   0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 2*N-1, RW, 24*N, IW( 2*N+1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 20
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR_2STAGE( 'N', 'I', 'U', 1, A, 1, &
                   0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 26*N, RW, 24*N-1, IW( 2*N-1 ), &
                   10*N, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 22
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHEEVR_2STAGE( 'N', 'I', 'U', 1, A, 1, &
                   0.0D0, 0.0D0, 1, 1, 0.0D0, &
                   M, R, Z, 1, IW, Q, 26*N, RW, 24*N, IW, 10*N-1, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHEEVR_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK )
      NT = NT + 13
!
!        ZHPEVD
!
      SRNAMT = 'ZHPEVD'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVD( '/', 'U', 0, A, X, Z, 1, W, 1, RW, 1, IW, 1, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVD( 'N', '/', 0, A, X, Z, 1, W, 1, RW, 1, IW, 1, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVD( 'N', 'U', -1, A, X, Z, 1, W, 1, RW, 1, IW, 1, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVD( 'V', 'U', 2, A, X, Z, 1, W, 4, RW, 25, IW, 12, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVD( 'N', 'U', 1, A, X, Z, 1, W, 0, RW, 1, IW, 1, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVD( 'N', 'U', 2, A, X, Z, 2, W, 1, RW, 2, IW, 1, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVD( 'V', 'U', 2, A, X, Z, 2, W, 2, RW, 25, IW, 12, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVD( 'N', 'U', 1, A, X, Z, 1, W, 1, RW, 0, IW, 1, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVD( 'N', 'U', 2, A, X, Z, 2, W, 2, RW, 1, IW, 1, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVD( 'V', 'U', 2, A, X, Z, 2, W, 4, RW, 18, IW, 12, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVD( 'N', 'U', 1, A, X, Z, 1, W, 1, RW, 1, IW, 0, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVD( 'N', 'U', 2, A, X, Z, 2, W, 2, RW, 2, IW, 0, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVD', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVD( 'V', 'U', 2, A, X, Z, 2, W, 4, RW, 25, IW, 2, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVD', INFOT, NOUT, LERR, OK )
      NT = NT + 13
!
!        ZHPEV
!
      SRNAMT = 'ZHPEV '
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEV( '/', 'U', 0, A, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEV ', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEV( 'N', '/', 0, A, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEV ', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEV( 'N', 'U', -1, A, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEV ', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEV( 'V', 'U', 2, A, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEV ', INFOT, NOUT, LERR, OK )
      NT = NT + 4
!
!        ZHPEVX
!
      SRNAMT = 'ZHPEVX'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVX( '/', 'A', 'U', 0, A, 0.0D0, 0.0D0, 0, 0, 0.0D0, M, &
                   X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVX', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVX( 'V', '/', 'U', 0, A, 0.0D0, 1.0D0, 1, 0, 0.0D0, M, &
                   X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVX', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVX( 'V', 'A', '/', 0, A, 0.0D0, 0.0D0, 0, 0, 0.0D0, M, &
                   X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVX', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVX( 'V', 'A', 'U', -1, A, 0.0D0, 0.0D0, 0, 0, 0.0D0, &
                   M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVX', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVX( 'V', 'V', 'U', 1, A, 0.0D0, 0.0D0, 0, 0, 0.0D0, M, &
                   X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVX', INFOT, NOUT, LERR, OK )
      INFOT = 8
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVX( 'V', 'I', 'U', 1, A, 0.0D0, 0.0D0, 0, 0, 0.0D0, M, &
                   X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVX', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVX( 'V', 'I', 'U', 2, A, 0.0D0, 0.0D0, 2, 1, 0.0D0, M, &
                   X, Z, 2, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVX', INFOT, NOUT, LERR, OK )
      INFOT = 14
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPEVX( 'V', 'A', 'U', 2, A, 0.0D0, 0.0D0, 0, 0, 0.0D0, M, &
                   X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHPEVX', INFOT, NOUT, LERR, OK )
      NT = NT + 8
!
!     Test error exits for the HB path.
!
   ELSE IF( LSAMEN( 2, C2, 'HB' ) ) THEN
!
!        ZHBTRD
!
      SRNAMT = 'ZHBTRD'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBTRD( '/', 'U', 0, 0, A, 1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBTRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBTRD', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBTRD( 'N', '/', 0, 0, A, 1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBTRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBTRD', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBTRD( 'N', 'U', -1, 0, A, 1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBTRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBTRD', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBTRD( 'N', 'U', 0, -1, A, 1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBTRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBTRD', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBTRD( 'N', 'U', 1, 1, A, 1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBTRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBTRD', INFOT, NOUT, LERR, OK )
      INFOT = 10
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBTRD( 'V', 'U', 2, 0, A, 1, D, E, Z, 1, W, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBTRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBTRD', INFOT, NOUT, LERR, OK )
      NT = NT + 6
!
!        ZHETRD_HB2ST
!
      SRNAMT = 'ZHETRD_HB2ST'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HB2ST( '/', 'N', 'U', 0, 0, A, 1, D, E, &
                                       C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HB2ST( 'N', '/', 'U', 0, 0, A, 1, D, E, &
                                       C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HB2ST( 'N', 'H', 'U', 0, 0, A, 1, D, E, &
                                       C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HB2ST( 'N', 'N', '/', 0, 0, A, 1, D, E, &
                                       C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HB2ST( 'N', 'N', 'U', -1, 0, A, 1, D, E, &
                                       C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HB2ST( 'N', 'N', 'U', 0, -1, A, 1, D, E, &
                                       C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HB2ST( 'N', 'N', 'U', 0, 1, A, 1, D, E, &
                                       C, 1, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HB2ST( 'N', 'N', 'U', 0, 0, A, 1, D, E, &
                                       C, 0, W, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHETRD_HB2ST( 'N', 'N', 'U', 0, 0, A, 1, D, E, &
                                       C, 1, W, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHETRD_HB2ST : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHETRD_HB2ST', INFOT, NOUT, LERR, OK )
      NT = NT + 9
!
!        ZHBEVD
!
      SRNAMT = 'ZHBEVD'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD( '/', 'U', 0, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD( 'N', '/', 0, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD( 'N', 'U', -1, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, &
                   1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD( 'N', 'U', 0, -1, A, 1, X, Z, 1, W, 1, RW, 1, IW, &
                   1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD( 'N', 'U', 2, 1, A, 1, X, Z, 1, W, 2, RW, 2, IW, 1, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD( 'V', 'U', 2, 1, A, 2, X, Z, 1, W, 8, RW, 25, IW, &
                   12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD( 'N', 'U', 1, 0, A, 1, X, Z, 1, W, 0, RW, 1, IW, 1, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD( 'N', 'U', 2, 1, A, 2, X, Z, 2, W, 1, RW, 2, IW, 1, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD( 'V', 'U', 2, 1, A, 2, X, Z, 2, W, 2, RW, 25, IW, &
                   12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD( 'N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, RW, 0, IW, 1, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD( 'N', 'U', 2, 1, A, 2, X, Z, 2, W, 2, RW, 1, IW, 1, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD( 'V', 'U', 2, 1, A, 2, X, Z, 2, W, 8, RW, 2, IW, &
                   12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD( 'N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 0, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD( 'N', 'U', 2, 1, A, 2, X, Z, 2, W, 2, RW, 2, IW, 0, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD', INFOT, NOUT, LERR, OK )
      INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD( 'V', 'U', 2, 1, A, 2, X, Z, 2, W, 8, RW, 25, IW, &
                   2, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD', INFOT, NOUT, LERR, OK )
      NT = NT + 15
!
!        ZHBEVD_2STAGE
!
      SRNAMT = 'ZHBEVD_2STAGE'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD_2STAGE( '/', 'U', 0, 0, A, 1, X, Z, 1, &
                              W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD_2STAGE( 'V', 'U', 0, 0, A, 1, X, Z, 1, &
                              W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD_2STAGE( 'N', '/', 0, 0, A, 1, X, Z, 1, &
                              W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD_2STAGE( 'N', 'U', -1, 0, A, 1, X, Z, 1, &
                               W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD_2STAGE( 'N', 'U', 0, -1, A, 1, X, Z, 1, &
                               W, 1, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD_2STAGE( 'N', 'U', 2, 1, A, 1, X, Z, 1, &
                              W, 2, RW, 2, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD_2STAGE( 'N', 'U', 2, 1, A, 2, X, Z, 0, &
                            W, 8, RW, 25, IW, 12, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD_2STAGE( 'N', 'U', 1, 0, A, 1, X, Z, 1, &
                              W, 0, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD_2STAGE( 'N', 'U', 2, 1, A, 2, X, Z, 2, &
                              W, 1, RW, 2, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
!         INFOT = 11
!         CALL ZHBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 2,
!     $                         W, 2, RW, 25, IW, 12, INFO )
!         CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD_2STAGE( 'N', 'U', 1, 0, A, 1, X, Z, 1, &
                              W, 1, RW, 0, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD_2STAGE( 'N', 'U', 2, 1, A, 2, X, Z, 2, &
                              W, 25, RW, 1, IW, 1, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
!         INFOT = 13
!         CALL ZHBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 2,
!     $                          W, 25, RW, 2, IW, 12, INFO )
!         CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD_2STAGE( 'N', 'U', 1, 0, A, 1, X, Z, 1, &
                              W, 1, RW, 1, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 15
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVD_2STAGE( 'N', 'U', 2, 1, A, 2, X, Z, 2, &
                              W, 25, RW, 2, IW, 0, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVD_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
!         INFOT = 15
!         CALL ZHBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 2,
!     $                          W, 25, RW, 25, IW, 2, INFO )
!         CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
      NT = NT + 13
!
!        ZHBEV
!
      SRNAMT = 'ZHBEV '
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEV( '/', 'U', 0, 0, A, 1, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEV ', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEV( 'N', '/', 0, 0, A, 1, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEV ', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEV( 'N', 'U', -1, 0, A, 1, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEV ', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEV( 'N', 'U', 0, -1, A, 1, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEV ', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEV( 'N', 'U', 2, 1, A, 1, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEV ', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEV( 'V', 'U', 2, 0, A, 1, X, Z, 1, W, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEV ', INFOT, NOUT, LERR, OK )
      NT = NT + 6
!
!        ZHBEV_2STAGE
!
      SRNAMT = 'ZHBEV_2STAGE '
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEV_2STAGE( '/', 'U', 0, 0, A, 1, X, &
                           Z, 1, W, 0, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEV_2STAGE( 'V', 'U', 0, 0, A, 1, X, &
                           Z, 1, W, 0, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEV_2STAGE( 'N', '/', 0, 0, A, 1, X, &
                           Z, 1, W, 0, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEV_2STAGE( 'N', 'U', -1, 0, A, 1, X, &
                            Z, 1, W, 0, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEV_2STAGE( 'N', 'U', 0, -1, A, 1, X, &
                            Z, 1, W, 0, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 6
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEV_2STAGE( 'N', 'U', 2, 1, A, 1, X, &
                           Z, 1, W, 0, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEV_2STAGE( 'N', 'U', 2, 0, A, 1, X, &
                           Z, 0, W, 0, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEV_2STAGE ', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEV_2STAGE( 'N', 'U', 2, 0, A, 1, X, &
                           Z, 1, W, 0, RW, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEV_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEV_2STAGE ', INFOT, NOUT, LERR, OK )
      NT = NT + 8
!
!        ZHBEVX
!
      SRNAMT = 'ZHBEVX'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX( '/', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0D0, 0.0D0, 0, &
                   0, 0.0D0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX( 'V', '/', 'U', 0, 0, A, 1, Q, 1, 0.0D0, 1.0D0, 1, &
                   0, 0.0D0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX( 'V', 'A', '/', 0, 0, A, 1, Q, 1, 0.0D0, 0.0D0, 0, &
                   0, 0.0D0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX( 'V', 'A', 'U', -1, 0, A, 1, Q, 1, 0.0D0, 0.0D0, 0, &
                   0, 0.0D0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX( 'V', 'A', 'U', 0, -1, A, 1, Q, 1, 0.0D0, 0.0D0, 0, &
                   0, 0.0D0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX( 'V', 'A', 'U', 2, 1, A, 1, Q, 2, 0.0D0, 0.0D0, 0, &
                   0, 0.0D0, M, X, Z, 2, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX', INFOT, NOUT, LERR, OK )
      INFOT = 9
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX( 'V', 'A', 'U', 2, 0, A, 1, Q, 1, 0.0D0, 0.0D0, 0, &
                   0, 0.0D0, M, X, Z, 2, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX( 'V', 'V', 'U', 1, 0, A, 1, Q, 1, 0.0D0, 0.0D0, 0, &
                   0, 0.0D0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX( 'V', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0D0, 0.0D0, 0, &
                   0, 0.0D0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX( 'V', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0D0, 0.0D0, 1, &
                   2, 0.0D0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX', INFOT, NOUT, LERR, OK )
      INFOT = 18
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX( 'V', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0D0, 0.0D0, 0, &
                   0, 0.0D0, M, X, Z, 1, W, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX', INFOT, NOUT, LERR, OK )
      NT = NT + 11
!
!        ZHBEVX_2STAGE
!
      SRNAMT = 'ZHBEVX_2STAGE'
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX_2STAGE( '/', 'A', 'U', 0, 0, A, 1, Q, 1, &
                          0.0D0, 0.0D0, 0, 0, 0.0D0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      INFOT = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX_2STAGE( 'V', 'A', 'U', 0, 0, A, 1, Q, 1, &
                          0.0D0, 0.0D0, 0, 0, 0.0D0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 2
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX_2STAGE( 'N', '/', 'U', 0, 0, A, 1, Q, 1, &
                          0.0D0, 1.0D0, 1, 0, 0.0D0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 3
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX_2STAGE( 'N', 'A', '/', 0, 0, A, 1, Q, 1, &
                          0.0D0, 0.0D0, 0, 0, 0.0D0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      INFOT = 4
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX_2STAGE( 'N', 'A', 'U', -1, 0, A, 1, Q, 1, &
                          0.0D0, 0.0D0, 0, 0, 0.0D0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 5
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX_2STAGE( 'N', 'A', 'U', 0, -1, A, 1, Q, 1, &
                          0.0D0, 0.0D0, 0, 0, 0.0D0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 7
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX_2STAGE( 'N', 'A', 'U', 2, 1, A, 1, Q, 2, &
                          0.0D0, 0.0D0, 0, 0, 0.0D0, &
                          M, X, Z, 2, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
!         INFOT = 9
!         CALL ZHBEVX_2STAGE( 'V', 'A', 'U', 2, 0, A, 1, Q, 1,
!     $                       0.0D0, 0.0D0, 0, 0, 0.0D0,
!     $                       M, X, Z, 2, W, 0, RW, IW, I3, INFO )
!         CALL CHKXER( 'ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 11
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX_2STAGE( 'N', 'V', 'U', 1, 0, A, 1, Q, 1, &
                          0.0D0, 0.0D0, 0, 0, 0.0D0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 12
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX_2STAGE( 'N', 'I', 'U', 1, 0, A, 1, Q, 1, &
                          0.0D0, 0.0D0, 0, 0, 0.0D0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 13
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX_2STAGE( 'N', 'I', 'U', 1, 0, A, 1, Q, 1, &
                          0.0D0, 0.0D0, 1, 2, 0.0D0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 18
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX_2STAGE( 'N', 'A', 'U', 2, 0, A, 1, Q, 2, &
                          0.0D0, 0.0D0, 0, 0, 0.0D0, &
                          M, X, Z, 0, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
      INFOT = 20
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHBEVX_2STAGE( 'N', 'A', 'U', 2, 0, A, 1, Q, 2, &
                          0.0D0, 0.0D0, 0, 0, 0.0D0, &
                          M, X, Z, 1, W, 0, RW, IW, I3, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHBEVX_2STAGE : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      CALL CHKXER( 'ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
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
!     End of ZERRST
!
END




