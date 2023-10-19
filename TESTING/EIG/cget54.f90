!> \brief \b CGET54
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGET54( N, A, LDA, B, LDB, S, LDS, T, LDT, U, LDU, V,
!                          LDV, WORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LDS, LDT, LDU, LDV, N
!       REAL               RESULT
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), B( LDB, * ), S( LDS, * ),
!      $                   T( LDT, * ), U( LDU, * ), V( LDV, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGET54 checks a generalized decomposition of the form
!>
!>          A = U*S*V'  and B = U*T* V'
!>
!> where ' means conjugate transpose and U and V are unitary.
!>
!> Specifically,
!>
!>   RESULT = ||( A - U*S*V', B - U*T*V' )|| / (||( A, B )||*n*ulp )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The size of the matrix.  If it is zero, SGET54 does nothing.
!>          It must be at least zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA, N)
!>          The original (unfactored) matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  It must be at least 1
!>          and at least N.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB, N)
!>          The original (unfactored) matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of B.  It must be at least 1
!>          and at least N.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is COMPLEX array, dimension (LDS, N)
!>          The factored matrix S.
!> \endverbatim
!>
!> \param[in] LDS
!> \verbatim
!>          LDS is INTEGER
!>          The leading dimension of S.  It must be at least 1
!>          and at least N.
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is COMPLEX array, dimension (LDT, N)
!>          The factored matrix T.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of T.  It must be at least 1
!>          and at least N.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is COMPLEX array, dimension (LDU, N)
!>          The orthogonal matrix on the left-hand side in the
!>          decomposition.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U.  LDU must be at least N and
!>          at least 1.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is COMPLEX array, dimension (LDV, N)
!>          The orthogonal matrix on the left-hand side in the
!>          decomposition.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of V.  LDV must be at least N and
!>          at least 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (3*N**2)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL
!>          The value RESULT, It is currently limited to 1/ulp, to
!>          avoid overflow. Errors are flagged by RESULT=10/ulp.
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
   SUBROUTINE CGET54( N, A, LDA, B, LDB, S, LDS, T, LDT, U, LDU, V, LDV, WORK, RESULT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, LDB, LDS, LDT, LDU, LDV, N
   REAL               RESULT
!     ..
!     .. Array Arguments ..
   COMPLEX            A( LDA, * ), B( LDB, * ), S( LDS, * ), &
                      T( LDT, * ), U( LDU, * ), V( LDV, * ), &
                      WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   REAL               ABNORM, ULP, UNFL, WNORM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   REAL               DUM( 1 )
!     ..
!     .. External Functions ..
   REAL               CLANGE, SLAMCH
   EXTERNAL           CLANGE, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMM, CLACPY
!     ..
!     .. Executable Statements ..
!
   RESULT = 0.0E+0
   IF( N <= 0 ) RETURN
!
!     Constants
!
   UNFL = SLAMCH( 'Safe minimum' )
   ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
!
!     compute the norm of (A,B)
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CLACPY( 'Full', N, N, A, LDA, WORK, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CLACPY( 'Full', N, N, B, LDB, WORK( N*N+1 ), N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   ABNORM = MAX( CLANGE( '1', N, 2*N, WORK, N, DUM ), UNFL )
!
!     Compute W1 = A - U*S*V', and put in the array WORK(1:N*N)
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CLACPY( ' ', N, N, A, LDA, WORK, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CGEMM( 'N', 'N', N, N, N, (1.0E+0,0.0E+0), U, LDU, S, LDS, (0.0E+0,0.0E+0), &
               WORK( N*N+1 ), N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CGEMM( 'N', 'C', N, N, N, -(1.0E+0,0.0E+0), WORK( N*N+1 ), N, V, LDV, &
               (1.0E+0,0.0E+0), WORK, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute W2 = B - U*T*V', and put in the workarray W(N*N+1:2*N*N)
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CLACPY( ' ', N, N, B, LDB, WORK( N*N+1 ), N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CGEMM( 'N', 'N', N, N, N, (1.0E+0,0.0E+0), U, LDU, T, LDT, (0.0E+0,0.0E+0), &
               WORK( 2*N*N+1 ), N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CGEMM( 'N', 'C', N, N, N, -(1.0E+0,0.0E+0), WORK( 2*N*N+1 ), N, V, LDV, &
               (1.0E+0,0.0E+0), WORK( N*N+1 ), N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute norm(W)/ ( ulp*norm((A,B)) )
!
   WNORM = CLANGE( '1', N, 2*N, WORK, N, DUM )
!
   IF( ABNORM > WNORM ) THEN
      RESULT = ( WNORM / ABNORM ) / ( 2*N*ULP )
   ELSE
      IF( ABNORM < 1.0E+0 ) THEN
         RESULT = ( MIN( WNORM, 2*N*ABNORM ) / ABNORM ) / ( 2*N*ULP )
      ELSE
         RESULT = MIN( WNORM / ABNORM, REAL( 2*N ) ) / ( 2*N*ULP )
      END IF
   END IF
!
   RETURN
!
!     End of CGET54
!
END




