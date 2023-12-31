!> \brief \b DQRT04
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DQRT04(M,N,NB,RESULT)
!
!       .. Scalar Arguments ..
!       INTEGER M, N, NB, LDT
!       .. Return values ..
!       DOUBLE PRECISION RESULT(6)
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DQRT04 tests DGEQRT and DGEMQRT.
!> \endverbatim
!
!  Arguments:
!  ==========
!
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
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          Block size of test matrix.  NB <= Min(M,N).
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (6)
!>          Results of each of the six tests below.
!>
!>          RESULT(1) = | A - Q R |
!>          RESULT(2) = | I - Q^H Q |
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
!> \ingroup double_lin
!
!  =====================================================================
   SUBROUTINE DQRT04(M,N,NB,RESULT)
   IMPLICIT NONE
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER M, N, NB, LDT
!     .. Return values ..
   DOUBLE PRECISION RESULT(6)
!
!  =====================================================================
!
!     ..
!     .. Local allocatable arrays
   DOUBLE PRECISION, ALLOCATABLE :: AF(:,:), Q(:,:), &
     R(:,:), RWORK(:), WORK( : ), T(:,:), &
     CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:)
!
!     .. Parameters ..
   DOUBLE PRECISION ONE, ZERO
   PARAMETER( ZERO = 0.0, ONE = 1.0 )
!     ..
!     .. Local Scalars ..
   INTEGER INFO, J, K, L, LWORK
   DOUBLE PRECISION   ANORM, EPS, RESID, CNORM, DNORM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   INTEGER            ISEED( 4 )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION DLAMCH, DLANGE, DLANSY
   LOGICAL  LSAME
   EXTERNAL DLAMCH, DLANGE, DLANSY, LSAME
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC  MAX, MIN
!     ..
!     .. Data statements ..
   DATA ISEED / 1988, 1989, 1990, 1991 /
!
   EPS = DLAMCH( 'Epsilon' )
   K = MIN(M,N)
   L = MAX(M,N)
   LWORK = MAX(2,L)*MAX(2,L)*NB
!
!     Dynamically allocate local arrays
!
   ALLOCATE ( A(M,N), AF(M,N), Q(M,M), R(M,L), RWORK(L), &
              WORK(LWORK), T(NB,N), C(M,N), CF(M,N), &
              D(N,M), DF(N,M) )
!
!     Put random numbers into A and copy to AF
!
   LDT=NB
   DO J=1,N
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLARNV( 2, ISEED, M, A( 1, J ) )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   END DO
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLACPY( 'Full', M, N, A, M, AF, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Factor the matrix A in the array AF.
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DGEQRT( M, N, NB, AF, M, T, LDT, WORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DGEQRT : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Generate the m-by-m matrix Q
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLASET( 'Full', M, M, ZERO, ONE, Q, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DGEMQRT( 'R', 'N', M, M, K, NB, AF, M, T, LDT, Q, M, &
                 WORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DGEMQRT : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Copy R
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLASET( 'Full', M, N, ZERO, ZERO, R, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLACPY( 'Upper', M, N, AF, M, R, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute |R - Q'*A| / |A| and store in RESULT(1)
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DGEMM( 'T', 'N', M, N, M, -ONE, Q, M, A, M, ONE, R, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   ANORM = DLANGE( '1', M, N, A, M, RWORK )
   RESID = DLANGE( '1', M, N, R, M, RWORK )
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
   CALL DLASET( 'Full', M, M, ZERO, ONE, R, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DSYRK( 'U', 'C', M, M, -ONE, Q, M, ONE, R, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DSYRK : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = DLANSY( '1', 'Upper', M, R, M, RWORK )
   RESULT( 2 ) = RESID / (EPS*MAX(1,M))
!
!     Generate random m-by-n matrix C and a copy CF
!
   DO J=1,N
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLARNV( 2, ISEED, M, C( 1, J ) )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   END DO
   CNORM = DLANGE( '1', M, N, C, M, RWORK)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLACPY( 'Full', M, N, C, M, CF, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Apply Q to C as Q*C
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DGEMQRT( 'L', 'N', M, N, K, NB, AF, M, T, NB, CF, M, &
                WORK, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DGEMQRT : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute |Q*C - Q*C| / |C|
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DGEMM( 'N', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = DLANGE( '1', M, N, CF, M, RWORK )
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
   CALL DLACPY( 'Full', M, N, C, M, CF, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Apply Q to C as QT*C
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DGEMQRT( 'L', 'T', M, N, K, NB, AF, M, T, NB, CF, M, &
                WORK, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DGEMQRT : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute |QT*C - QT*C| / |C|
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DGEMM( 'T', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = DLANGE( '1', M, N, CF, M, RWORK )
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
      CALL DLARNV( 2, ISEED, N, D( 1, J ) )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   END DO
   DNORM = DLANGE( '1', N, M, D, N, RWORK)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLACPY( 'Full', N, M, D, N, DF, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Apply Q to D as D*Q
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DGEMQRT( 'R', 'N', N, M, K, NB, AF, M, T, NB, DF, N, &
                WORK, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DGEMQRT : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute |D*Q - D*Q| / |D|
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DGEMM( 'N', 'N', N, M, M, -ONE, D, N, Q, M, ONE, DF, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = DLANGE( '1', N, M, DF, N, RWORK )
   IF( CNORM > ZERO ) THEN
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
   CALL DLACPY( 'Full', N, M, D, N, DF, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Apply Q to D as D*QT
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DGEMQRT( 'R', 'T', N, M, K, NB, AF, M, T, NB, DF, N, &
                WORK, INFO)
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DGEMQRT : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute |D*QT - D*QT| / |D|
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DGEMM( 'N', 'T', N, M, M, -ONE, D, N, Q, M, ONE, DF, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = DLANGE( '1', N, M, DF, N, RWORK )
   IF( CNORM > ZERO ) THEN
      RESULT( 6 ) = RESID / (EPS*MAX(1,M)*DNORM)
   ELSE
      RESULT( 6 ) = ZERO
   END IF
!
!     Deallocate all arrays
!
   DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T, C, D, CF, DF)
!
   RETURN
   END

                                                                                                                                                                                                                                                                                                            




