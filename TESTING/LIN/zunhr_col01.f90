!> \brief \b ZUNHR_COL01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNHR_COL01( M, N, MB1, NB1, NB2, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER           M, N, MB1, NB1, NB2
!       .. Return values ..
!       DOUBLE PRECISION  RESULT(6)
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNHR_COL01 tests ZUNGTSQR and ZUNHR_COL using ZLATSQR, ZGEMQRT.
!> Therefore, ZLATSQR (part of ZGEQR), ZGEMQRT (part of ZGEMQR)
!> have to be tested before this test.
!>
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
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          Number of columns in test matrix.
!> \endverbatim
!> \param[in] MB1
!> \verbatim
!>          MB1 is INTEGER
!>          Number of row in row block in an input test matrix.
!> \endverbatim
!>
!> \param[in] NB1
!> \verbatim
!>          NB1 is INTEGER
!>          Number of columns in column block an input test matrix.
!> \endverbatim
!>
!> \param[in] NB2
!> \verbatim
!>          NB2 is INTEGER
!>          Number of columns in column block in an output test matrix.
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (6)
!>          Results of each of the six tests below.
!>
!>            A is a m-by-n test input matrix to be factored.
!>            so that A = Q_gr * ( R )
!>                               ( 0 ),
!>
!>            Q_qr is an implicit m-by-m unitary Q matrix, the result
!>            of factorization in blocked WY-representation,
!>            stored in ZGEQRT output format.
!>
!>            R is a n-by-n upper-triangular matrix,
!>
!>            0 is a (m-n)-by-n zero matrix,
!>
!>            Q is an explicit m-by-m unitary matrix Q = Q_gr * I
!>
!>            C is an m-by-n random matrix,
!>
!>            D is an n-by-m random matrix.
!>
!>          The six tests are:
!>
!>          RESULT(1) = |R - (Q**H) * A| / ( eps * m * |A| )
!>            is equivalent to test for | A - Q * R | / (eps * m * |A|),
!>
!>          RESULT(2) = |I - (Q**H) * Q| / ( eps * m ),
!>
!>          RESULT(3) = | Q_qr * C - Q * C | / (eps * m * |C|),
!>
!>          RESULT(4) = | (Q_gr**H) * C - (Q**H) * C | / (eps * m * |C|)
!>
!>          RESULT(5) = | D * Q_qr - D * Q | / (eps * m * |D|)
!>
!>          RESULT(6) = | D * (Q_qr**H) - D * (Q**H) | / (eps * m * |D|),
!>
!>          where:
!>            Q_qr * C, (Q_gr**H) * C, D * Q_qr, D * (Q_qr**H) are
!>            computed using ZGEMQRT,
!>
!>            Q * C, (Q**H) * C, D * Q, D * (Q**H)  are
!>            computed using ZGEMM.
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
   SUBROUTINE ZUNHR_COL01( M, N, MB1, NB1, NB2, RESULT )
   IMPLICIT NONE
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER           M, N, MB1, NB1, NB2
!     .. Return values ..
   DOUBLE PRECISION  RESULT(6)
!
!  =====================================================================
!
!     ..
!     .. Local allocatable arrays
   COMPLEX*16      , ALLOCATABLE ::  A(:,:), AF(:,:), Q(:,:), R(:,:), &
                      WORK( : ), T1(:,:), T2(:,:), DIAG(:), &
                      C(:,:), CF(:,:), D(:,:), DF(:,:)
   DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO
   PARAMETER          ( ZERO = 0.0D+0 )
   COMPLEX*16         CONE, CZERO
   PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ), &
                        CZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
   LOGICAL            TESTZEROS
   INTEGER            INFO, I, J, K, L, LWORK, NB1_UB, NB2_UB, NRB
   DOUBLE PRECISION   ANORM, EPS, RESID, CNORM, DNORM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   INTEGER            ISEED( 4 )
   COMPLEX*16         WORKQUERY( 1 )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH, ZLANGE, ZLANSY
   EXTERNAL           DLAMCH, ZLANGE, ZLANSY
!     ..
!     .. External Subroutines ..
   EXTERNAL           ZLACPY, ZLARNV, ZLASET, ZLATSQR, ZUNHR_COL, &
                      ZUNGTSQR, ZSCAL, ZGEMM, ZGEMQRT, ZHERK
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          CEILING, DBLE, MAX, MIN
!     ..
!     .. Scalars in Common ..
   CHARACTER(LEN=32)  SRNAMT
!     ..
!     .. Common blocks ..
   COMMON             / SRMNAMC / SRNAMT
!     ..
!     .. Data statements ..
   DATA ISEED / 1988, 1989, 1990, 1991 /
!
!     TEST MATRICES WITH HALF OF MATRIX BEING ZEROS
!
   TESTZEROS = .FALSE.
!
   EPS = DLAMCH( 'Epsilon' )
   K = MIN( M, N )
   L = MAX( M, N, 1)
!
!     Dynamically allocate local arrays
!
   ALLOCATE ( A(M,N), AF(M,N), Q(L,L), R(M,L), RWORK(L), &
              C(M,N), CF(M,N), &
              D(N,M), DF(N,M) )
!
!     Put random numbers into A and copy to AF
!
   DO J = 1, N
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
   IF( TESTZEROS ) THEN
      IF( M >= 4 ) THEN
         DO J = 1, N
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
!     Number of row blocks in ZLATSQR
!
   NRB = MAX( 1, CEILING( DBLE( M - N ) / DBLE( MB1 - N ) ) )
!
   ALLOCATE ( T1( NB1, N * NRB ) )
   ALLOCATE ( T2( NB2, N ) )
   ALLOCATE ( DIAG( N ) )
!
!     Begin determine LWORK for the array WORK and allocate memory.
!
!     ZLATSQR requires NB1 to be bounded by N.
!
   NB1_UB = MIN( NB1, N)
!
!     ZGEMQRT requires NB2 to be bounded by N.
!
   NB2_UB = MIN( NB2, N)
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLATSQR( M, N, MB1, NB1_UB, AF, M, T1, NB1, &
                 WORKQUERY, -1, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLATSQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   LWORK = INT( WORKQUERY( 1 ) )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGTSQR( M, N, MB1, NB1, AF, M, T1, NB1, WORKQUERY, -1, &
                  INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGTSQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif

   LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) )
!
!     In ZGEMQRT, WORK is N*NB2_UB if SIDE = 'L',
!                or  M*NB2_UB if SIDE = 'R'.
!
   LWORK = MAX( LWORK, NB2_UB * N, NB2_UB * M )
!
   ALLOCATE ( WORK( LWORK ) )
!
!     End allocate memory for WORK.
!
!
!     Begin Householder reconstruction routines
!
!     Factor the matrix A in the array AF.
!
   SRNAMT = 'ZLATSQR'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLATSQR( M, N, MB1, NB1_UB, AF, M, T1, NB1, WORK, LWORK, &
                 INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLATSQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Copy the factor R into the array R.
!
   SRNAMT = 'ZLACPY'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'U', N, N, AF, M, R, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Reconstruct the orthogonal matrix Q.
!
   SRNAMT = 'ZUNGTSQR'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGTSQR( M, N, MB1, NB1, AF, M, T1, NB1, WORK, LWORK, &
                  INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGTSQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Perform the Householder reconstruction, the result is stored
!     the arrays AF and T2.
!
   SRNAMT = 'ZUNHR_COL'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNHR_COL( M, N, NB2, AF, M, T2, NB2, DIAG, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNHR_COL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute the factor R_hr corresponding to the Householder
!     reconstructed Q_hr and place it in the upper triangle of AF to
!     match the Q storage format in ZGEQRT. R_hr = R_tsqr * S,
!     this means changing the sign of I-th row of the matrix R_tsqr
!     according to sign of of I-th diagonal element DIAG(I) of the
!     matrix S.
!
   SRNAMT = 'ZLACPY'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'U', N, N, R, M, AF, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   DO I = 1, N
      IF( DIAG( I ) == -CONE ) THEN
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZSCAL( N+1-I, -CONE, AF( I, I ), M )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZSCAL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      END IF
   END DO
!
!     End Householder reconstruction routines.
!
!
!     Generate the m-by-m matrix Q
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', M, M, CZERO, CONE, Q, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   SRNAMT = 'ZGEMQRT'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQRT( 'L', 'N', M, M, K, NB2_UB, AF, M, T2, NB2, Q, M, &
                 WORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQRT : ',&
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
!
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
!     TEST 1
!     Compute |R - (Q**H)*A| / ( eps * m * |A| ) and store in RESULT(1)
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'C', 'N', M, N, M, -CONE, Q, M, A, M, CONE, R, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   ANORM = ZLANGE( '1', M, N, A, M, RWORK )
   RESID = ZLANGE( '1', M, N, R, M, RWORK )
   IF( ANORM > ZERO ) THEN
      RESULT( 1 ) = RESID / ( EPS * MAX( 1, M ) * ANORM )
   ELSE
      RESULT( 1 ) = ZERO
   END IF
!
!     TEST 2
!     Compute |I - (Q**H)*Q| / ( eps * m ) and store in RESULT(2)
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', M, M, CZERO, CONE, R, M )
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
   CALL ZHERK( 'U', 'C', M, M, REAL(-CONE), Q, M, REAL(CONE), R, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZHERK : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = ZLANSY( '1', 'Upper', M, R, M, RWORK )
   RESULT( 2 ) = RESID / ( EPS * MAX( 1, M ) )
!
!     Generate random m-by-n matrix C
!
   DO J = 1, N
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
   CNORM = ZLANGE( '1', M, N, C, M, RWORK )
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
!     Apply Q to C as Q*C = CF
!
   SRNAMT = 'ZGEMQRT'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQRT( 'L', 'N', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, &
                  WORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQRT : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     TEST 3
!     Compute |CF - Q*C| / ( eps *  m * |C| )
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'N', 'N', M, N, M, -CONE, Q, M, C, M, CONE, CF, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = ZLANGE( '1', M, N, CF, M, RWORK )
   IF( CNORM > ZERO ) THEN
      RESULT( 3 ) = RESID / ( EPS * MAX( 1, M ) * CNORM )
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
!     Apply Q to C as (Q**H)*C = CF
!
   SRNAMT = 'ZGEMQRT'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQRT( 'L', 'C', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, &
                  WORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQRT : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     TEST 4
!     Compute |CF - (Q**H)*C| / ( eps * m * |C|)
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'C', 'N', M, N, M, -CONE, Q, M, C, M, CONE, CF, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = ZLANGE( '1', M, N, CF, M, RWORK )
   IF( CNORM > ZERO ) THEN
      RESULT( 4 ) = RESID / ( EPS * MAX( 1, M ) * CNORM )
   ELSE
      RESULT( 4 ) = ZERO
   END IF
!
!     Generate random n-by-m matrix D and a copy DF
!
   DO J = 1, M
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
   DNORM = ZLANGE( '1', N, M, D, N, RWORK )
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
!     Apply Q to D as D*Q = DF
!
   SRNAMT = 'ZGEMQRT'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQRT( 'R', 'N', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, &
                  WORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQRT : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     TEST 5
!     Compute |DF - D*Q| / ( eps * m * |D| )
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'N', 'N', N, M, M, -CONE, D, N, Q, M, CONE, DF, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = ZLANGE( '1', N, M, DF, N, RWORK )
   IF( DNORM > ZERO ) THEN
      RESULT( 5 ) = RESID / ( EPS * MAX( 1, M ) * DNORM )
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
!     Apply Q to D as D*QT = DF
!
   SRNAMT = 'ZGEMQRT'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMQRT( 'R', 'C', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, &
                  WORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMQRT : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     TEST 6
!     Compute |DF - D*(Q**H)| / ( eps * m * |D| )
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'N', 'C', N, M, M, -CONE, D, N, Q, M, CONE, DF, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RESID = ZLANGE( '1', N, M, DF, N, RWORK )
   IF( DNORM > ZERO ) THEN
      RESULT( 6 ) = RESID / ( EPS * MAX( 1, M ) * DNORM )
   ELSE
      RESULT( 6 ) = ZERO
   END IF
!
!     Deallocate all arrays
!
   DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T1, T2, DIAG, &
                C, D, CF, DF )
!
   RETURN
!
!     End of ZUNHR_COL01
!
END
                                                                                                                                                                                                                                                                                                            




