!> \brief \b SLSETS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLSETS( M, P, N, A, AF, LDA, B, BF, LDB, C, CF,
!                          D, DF, X, WORK, LWORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LWORK, M, P, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), AF( LDA, * ), B( LDB, * ),
!      $                   BF( LDB, * ), RESULT( 2 ), RWORK( * ),
!      $                   C( * ), D( * ), CF( * ), DF( * ),
!      $                   WORK( LWORK ), X( * )
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLSETS tests SGGLSE - a subroutine for solving linear equality
!> constrained least square problem (LSE).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          P is INTEGER
!>          The number of rows of the matrix B.  P >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrices A and B.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The M-by-N matrix A.
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is REAL array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A, AF, Q and R.
!>          LDA >= max(M,N).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (LDB,N)
!>          The P-by-N matrix A.
!> \endverbatim
!>
!> \param[out] BF
!> \verbatim
!>          BF is REAL array, dimension (LDB,N)
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the arrays B, BF, V and S.
!>          LDB >= max(P,N).
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is REAL array, dimension( M )
!>          the vector C in the LSE problem.
!> \endverbatim
!>
!> \param[out] CF
!> \verbatim
!>          CF is REAL array, dimension( M )
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension( P )
!>          the vector D in the LSE problem.
!> \endverbatim
!>
!> \param[out] DF
!> \verbatim
!>          DF is REAL array, dimension( P )
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is REAL array, dimension( N )
!>          solution vector X in the LSE problem.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (2)
!>          The test ratios:
!>            RESULT(1) = norm( A*x - c )/ norm(A)*norm(X)*EPS
!>            RESULT(2) = norm( B*x - d )/ norm(B)*norm(X)*EPS
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
!> \ingroup single_eig
!
!  =====================================================================
   SUBROUTINE SLSETS( M, P, N, A, AF, LDA, B, BF, LDB, C, CF, &
                      D, DF, X, WORK, LWORK, RWORK, RESULT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, LDB, LWORK, M, P, N
!     ..
!     .. Array Arguments ..
   REAL               A( LDA, * ), AF( LDA, * ), B( LDB, * ), &
                      BF( LDB, * ), RESULT( 2 ), RWORK( * ), &
                      C( * ), D( * ), CF( * ), DF( * ), &
                      WORK( LWORK ), X( * )
!
!  ====================================================================
!
!     ..
!     .. Local Scalars ..
   INTEGER            INFO
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Subroutines ..
   EXTERNAL           SGGLSE, SLACPY, SGET02
!     ..
!     .. Executable Statements ..
!
!     Copy the matrices A and B to the arrays AF and BF,
!     and the vectors C and D to the arrays CF and DF,
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SLACPY( 'Full', M, N, A, LDA, AF, LDA )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SLACPY( 'Full', P, N, B, LDB, BF, LDB )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SCOPY( M, C, 1, CF, 1 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SCOPY( P, D, 1, DF, 1 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Solve LSE problem
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SGGLSE( M, N, P, AF, LDA, BF, LDB, CF, DF, X, &
                WORK, LWORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SGGLSE : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Test the residual for the solution of LSE
!
!     Compute RESULT(1) = norm( A*x - c ) / norm(A)*norm(X)*EPS
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SCOPY( M, C, 1, CF, 1 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SCOPY( P, D, 1, DF, 1 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CALL SGET02( 'No transpose', M, N, 1, A, LDA, X, N, CF, M, &
                RWORK, RESULT( 1 ) )
!
!     Compute result(2) = norm( B*x - d ) / norm(B)*norm(X)*EPS
!
   CALL SGET02( 'No transpose', P, N, 1, B, LDB, X, N, DF, P, &
                RWORK, RESULT( 2 ) )
!
   RETURN
!
!     End of SLSETS
!
END




