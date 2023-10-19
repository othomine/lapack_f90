!> \brief \b ZTSQR01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTSQR01(TSSW, M,N, MB, NB, RESULT)
!
!       .. Scalar Arguments ..
!       INTEGER M, N, MB
!       .. Return values ..
!       DOUBLE PRECISION RESULT(6)
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTSQR01 tests ZGEQR , ZGELQ, ZGEMLQ and ZGEMQR.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TSSW
!> \verbatim
!>          TSSW is CHARACTER
!>          'TS' for testing tall skinny QR
!>               and anything else for testing short wide LQ
!> \endverbatim
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          Number of rows in test matrix.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          Number of columns in test matrix.
!> \endverbatim
!> \param[in] MB
!> \verbatim
!>          MB is INTEGER
!>          Number of row in row block in test matrix.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          Number of columns in column block test matrix.
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (6)
!>          Results of each of the six tests below.
!>
!>          RESULT(1) = | A - Q R | or | A - L Q |
!>          RESULT(2) = | I - Q^H Q | or | I - Q Q^H |
!>          RESULT(3) = | Q C - Q C |
!>          RESULT(4) = | Q^H C - Q^H C |
!>          RESULT(5) = | C Q - C Q |
!>          RESULT(6) = | C Q^H - C Q^H |
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
!  =====================================================================
   SUBROUTINE ZTSQR01(TSSW, M, N, MB, NB, RESULT)
   IMPLICIT NONE
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER         TSSW
   INTEGER           M, N, MB, NB
!     .. Return values ..
   DOUBLE PRECISION  RESULT(6)
!
!  =====================================================================
!
!     ..
!     .. Local allocatable arrays
   COMPLEX*16, ALLOCATABLE :: AF(:,:), Q(:,:), &
     R(:,:), WORK( : ), T(:), &
     CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:), LQ(:,:)
   DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
!
!     .. Parameters ..
   DOUBLE PRECISION ZERO
   COMPLEX*16 ONE, CZERO
   PARAMETER( ZERO = 0.0, ONE = (1.0,0.0), CZERO=(0.0,0.0) )
!     ..
!     .. Local Scalars ..
   LOGICAL TESTZEROS, TS
   INTEGER INFO, J, K, L, LWORK, TSIZE, MNB
   DOUBLE PRECISION   ANORM, EPS, RESID, CNORM, DNORM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   INTEGER            ISEED( 4 )
   COMPLEX*16         TQUERY( 5 ), WORKQUERY( 1 )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION DLAMCH, ZLANGE, ZLANSY
   LOGICAL  LSAME
   INTEGER ILAENV
   EXTERNAL DLAMCH, ZLANGE, ZLANSY, LSAME, ILAENV
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC  MAX, MIN
!     .. Scalars in Common ..
   CHARACTER*32       srnamt
!     ..
!     .. Common blocks ..
   COMMON             / srnamc / srnamt
!     ..
!     .. Data statements ..
   DATA ISEED / 1988, 1989, 1990, 1991 /
!
!     TEST TALL SKINNY OR SHORT WIDE
!
   TS = LSAME(TSSW, 'TS')
!
!     TEST MATRICES WITH HALF OF MATRIX BEING ZEROS
!
   TESTZEROS = .FALSE.
!
   EPS = DLAMCH( 'Epsilon' )
   K = MIN(M,N)
   L = MAX(M,N,1)
   MNB = MAX ( MB, NB)
   LWORK = MAX(3,L)*MNB
!
!     Dynamically allocate local arrays
!
   ALLOCATE ( A(M,N), AF(M,N), Q(L,L), R(M,L), RWORK(L), &
              C(M,N), CF(M,N), &
              D(N,M), DF(N,M), LQ(L,N) )
!
!     Put random numbers into A and copy to AF
!
   DO J=1,N
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLARNV( 2, ISEED, M, A( 1, J ) )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   END DO
   IF (TESTZEROS) THEN
      IF (M >= 4) THEN
         DO J=1,N
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZLARNV( 2, ISEED, M/2, A( M/4, J ) )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZLARNV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         END DO
      END IF
   END IF
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'Full', M, N, A, M, AF, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   IF (TS) THEN
!
!     Factor the matrix A in the array AF.
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEQR( M, N, AF, M, TQUERY, -1, WORKQUERY, -1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   TSIZE = INT( TQUERY( 1 ) )
   LWORK = INT( WORKQUERY( 1 ) )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'L', 'N', M, M, K, AF, M, TQUERY, TSIZE, CF, M, &
                WORKQUERY, -1, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'L', 'N', M, N, K, AF, M, TQUERY, TSIZE, CF, M, &
                WORKQUERY, -1, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'L', 'C', M, N, K, AF, M, TQUERY, TSIZE, CF, M, &
                WORKQUERY, -1, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'R', 'N', N, M, K, AF, M, TQUERY, TSIZE, DF, N, &
                WORKQUERY, -1, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'R', 'C', N, M, K, AF, M, TQUERY, TSIZE, DF, N, &
                WORKQUERY, -1, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) )
   ALLOCATE ( T( TSIZE ) )
   ALLOCATE ( WORK( LWORK ) )
   srnamt = 'ZGEQR'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEQR( M, N, AF, M, T, TSIZE, WORK, LWORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Generate the m-by-m matrix Q
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', M, M, CZERO, ONE, Q, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   srnamt = 'ZGEMQR'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'L', 'N', M, M, K, AF, M, T, TSIZE, Q, M, &
                 WORK, LWORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Copy R
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', M, N, CZERO, CZERO, R, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'Upper', M, N, AF, M, R, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute |R - Q'*A| / |A| and store in RESULT(1)
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'C', 'N', M, N, M, -ONE, Q, M, A, M, ONE, R, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   ANORM = ZLANGE( '1', M, N, A, M, RWORK )
   RESID = ZLANGE( '1', M, N, R, M, RWORK )
   IF( ANORM > ZERO ) THEN
      RESULT( 1 ) = RESID / (EPS*MAX(1,M)*ANORM)
   ELSE
      RESULT( 1 ) = ZERO
   END IF
!
!     Compute |I - Q'*Q| and store in RESULT(2)
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', M, M, CZERO, ONE, R, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZHERK( 'U', 'C', M, M, DREAL(-ONE), Q, M, DREAL(ONE), R, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZHERK : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = ZLANSY( '1', 'Upper', M, R, M, RWORK )
   RESULT( 2 ) = RESID / (EPS*MAX(1,M))
!
!     Generate random m-by-n matrix C and a copy CF
!
   DO J=1,N
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLARNV( 2, ISEED, M, C( 1, J ) )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   END DO
   CNORM = ZLANGE( '1', M, N, C, M, RWORK)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'Full', M, N, C, M, CF, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Apply Q to C as Q*C
!
   srnamt = 'ZGEMQR'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'L', 'N', M, N, K, AF, M, T, TSIZE, CF, M, &
                WORK, LWORK, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute |Q*C - Q*C| / |C|
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'N', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = ZLANGE( '1', M, N, CF, M, RWORK )
   IF( CNORM > ZERO ) THEN
      RESULT( 3 ) = RESID / (EPS*MAX(1,M)*CNORM)
   ELSE
      RESULT( 3 ) = ZERO
   END IF
!
!     Copy C into CF again
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'Full', M, N, C, M, CF, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Apply Q to C as QT*C
!
   srnamt = 'ZGEMQR'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'L', 'C', M, N, K, AF, M, T, TSIZE, CF, M, &
                WORK, LWORK, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute |QT*C - QT*C| / |C|
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'C', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = ZLANGE( '1', M, N, CF, M, RWORK )
   IF( CNORM > ZERO ) THEN
      RESULT( 4 ) = RESID / (EPS*MAX(1,M)*CNORM)
   ELSE
      RESULT( 4 ) = ZERO
   END IF
!
!     Generate random n-by-m matrix D and a copy DF
!
   DO J=1,M
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLARNV( 2, ISEED, N, D( 1, J ) )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   END DO
   DNORM = ZLANGE( '1', N, M, D, N, RWORK)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'Full', N, M, D, N, DF, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Apply Q to D as D*Q
!
   srnamt = 'ZGEMQR'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'R', 'N', N, M, K, AF, M, T, TSIZE, DF, N, &
                WORK, LWORK, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute |D*Q - D*Q| / |D|
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'N', 'N', N, M, M, -ONE, D, N, Q, M, ONE, DF, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = ZLANGE( '1', N, M, DF, N, RWORK )
   IF( DNORM > ZERO ) THEN
      RESULT( 5 ) = RESID / (EPS*MAX(1,M)*DNORM)
   ELSE
      RESULT( 5 ) = ZERO
   END IF
!
!     Copy D into DF again
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'Full', N, M, D, N, DF, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Apply Q to D as D*QT
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQR( 'R', 'C', N, M, K, AF, M, T, TSIZE, DF, N, &
                WORK, LWORK, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute |D*QT - D*QT| / |D|
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'N', 'C', N, M, M, -ONE, D, N, Q, M, ONE, DF, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = ZLANGE( '1', N, M, DF, N, RWORK )
   IF( CNORM > ZERO ) THEN
      RESULT( 6 ) = RESID / (EPS*MAX(1,M)*DNORM)
   ELSE
      RESULT( 6 ) = ZERO
   END IF
!
!     Short and wide
!
   ELSE
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGELQ( M, N, AF, M, TQUERY, -1, WORKQUERY, -1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGELQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   TSIZE = INT( TQUERY( 1 ) )
   LWORK = INT( WORKQUERY( 1 ) )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'R', 'N', N, N, K, AF, M, TQUERY, TSIZE, Q, N, &
                 WORKQUERY, -1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'L', 'N', N, M, K, AF, M, TQUERY, TSIZE, DF, N, &
                WORKQUERY, -1, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'L', 'C', N, M, K, AF, M, TQUERY, TSIZE, DF, N, &
                WORKQUERY, -1, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'R', 'N', M, N, K, AF, M, TQUERY, TSIZE, CF, M, &
                WORKQUERY, -1, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'R', 'C', M, N, K, AF, M, TQUERY, TSIZE, CF, M, &
                WORKQUERY, -1, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) )
   ALLOCATE ( T( TSIZE ) )
   ALLOCATE ( WORK( LWORK ) )
   srnamt = 'ZGELQ'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGELQ( M, N, AF, M, T, TSIZE, WORK, LWORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGELQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!
!     Generate the n-by-n matrix Q
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', N, N, CZERO, ONE, Q, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   srnamt = 'ZGEMLQ'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'R', 'N', N, N, K, AF, M, T, TSIZE, Q, N, &
                 WORK, LWORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Copy R
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', M, N, CZERO, CZERO, LQ, L )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'Lower', M, N, AF, M, LQ, L )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute |L - A*Q'| / |A| and store in RESULT(1)
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'N', 'C', M, N, N, -ONE, A, M, Q, N, ONE, LQ, L )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   ANORM = ZLANGE( '1', M, N, A, M, RWORK )
   RESID = ZLANGE( '1', M, N, LQ, L, RWORK )
   IF( ANORM > ZERO ) THEN
      RESULT( 1 ) = RESID / (EPS*MAX(1,N)*ANORM)
   ELSE
      RESULT( 1 ) = ZERO
   END IF
!
!     Compute |I - Q'*Q| and store in RESULT(2)
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', N, N, CZERO, ONE, LQ, L )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZHERK( 'U', 'C', N, N, DREAL(-ONE), Q, N, DREAL(ONE), LQ, L)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZHERK : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = ZLANSY( '1', 'Upper', N, LQ, L, RWORK )
   RESULT( 2 ) = RESID / (EPS*MAX(1,N))
!
!     Generate random m-by-n matrix C and a copy CF
!
   DO J=1,M
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLARNV( 2, ISEED, N, D( 1, J ) )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   END DO
   DNORM = ZLANGE( '1', N, M, D, N, RWORK)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'Full', N, M, D, N, DF, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Apply Q to C as Q*C
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'L', 'N', N, M, K, AF, M, T, TSIZE, DF, N, &
                WORK, LWORK, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute |Q*D - Q*D| / |D|
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'N', 'N', N, M, N, -ONE, Q, N, D, N, ONE, DF, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = ZLANGE( '1', N, M, DF, N, RWORK )
   IF( DNORM > ZERO ) THEN
      RESULT( 3 ) = RESID / (EPS*MAX(1,N)*DNORM)
   ELSE
      RESULT( 3 ) = ZERO
   END IF
!
!     Copy D into DF again
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'Full', N, M, D, N, DF, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Apply Q to D as QT*D
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'L', 'C', N, M, K, AF, M, T, TSIZE, DF, N, &
                WORK, LWORK, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute |QT*D - QT*D| / |D|
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'C', 'N', N, M, N, -ONE, Q, N, D, N, ONE, DF, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = ZLANGE( '1', N, M, DF, N, RWORK )
   IF( DNORM > ZERO ) THEN
      RESULT( 4 ) = RESID / (EPS*MAX(1,N)*DNORM)
   ELSE
      RESULT( 4 ) = ZERO
   END IF
!
!     Generate random n-by-m matrix D and a copy DF
!
   DO J=1,N
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLARNV( 2, ISEED, M, C( 1, J ) )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   END DO
   CNORM = ZLANGE( '1', M, N, C, M, RWORK)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'Full', M, N, C, M, CF, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Apply Q to C as C*Q
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'R', 'N', M, N, K, AF, M, T, TSIZE, CF, M, &
                WORK, LWORK, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute |C*Q - C*Q| / |C|
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'N', 'N', M, N, N, -ONE, C, M, Q, N, ONE, CF, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = ZLANGE( '1', N, M, DF, N, RWORK )
   IF( CNORM > ZERO ) THEN
      RESULT( 5 ) = RESID / (EPS*MAX(1,N)*CNORM)
   ELSE
      RESULT( 5 ) = ZERO
   END IF
!
!     Copy C into CF again
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'Full', M, N, C, M, CF, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Apply Q to D as D*QT
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMLQ( 'R', 'C', M, N, K, AF, M, T, TSIZE, CF, M, &
                WORK, LWORK, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute |C*QT - C*QT| / |C|
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'N', 'C', M, N, N, -ONE, C, M, Q, N, ONE, CF, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = ZLANGE( '1', M, N, CF, M, RWORK )
   IF( CNORM > ZERO ) THEN
      RESULT( 6 ) = RESID / (EPS*MAX(1,N)*CNORM)
   ELSE
      RESULT( 6 ) = ZERO
   END IF
!
   END IF
!
!     Deallocate all arrays
!
   DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T, C, D, CF, DF)
!
   RETURN
END
                                                                                                                                                                                                                                                                                                            




