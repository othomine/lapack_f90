!> \brief \b ZGLMTS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGLMTS( N, M, P, A, AF, LDA, B, BF, LDB, D, DF, X, U,
!                          WORK, LWORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LWORK, M, N, P
!       DOUBLE PRECISION   RESULT
!       ..
!       .. Array Arguments ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGLMTS tests ZGGGLM - a subroutine for solving the generalized
!> linear model problem.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows of the matrices A and B.  N >= 0.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of columns of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          P is INTEGER
!>          The number of columns of the matrix B.  P >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,M)
!>          The N-by-M matrix A.
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is COMPLEX*16 array, dimension (LDA,M)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A, AF. LDA >= max(M,N).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,P)
!>          The N-by-P matrix A.
!> \endverbatim
!>
!> \param[out] BF
!> \verbatim
!>          BF is COMPLEX*16 array, dimension (LDB,P)
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the arrays B, BF. LDB >= max(P,N).
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is COMPLEX*16 array, dimension( N )
!>          On input, the left hand side of the GLM.
!> \endverbatim
!>
!> \param[out] DF
!> \verbatim
!>          DF is COMPLEX*16 array, dimension( N )
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension( M )
!>          solution vector X in the GLM problem.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX*16 array, dimension( P )
!>          solution vector U in the GLM problem.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
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
!>          RWORK is DOUBLE PRECISION array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION
!>          The test ratio:
!>                           norm( d - A*x - B*u )
!>            RESULT = -----------------------------------------
!>                     (norm(A)+norm(B))*(norm(x)+norm(u))*EPS
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
   SUBROUTINE ZGLMTS( N, M, P, A, AF, LDA, B, BF, LDB, D, DF, X, U, &
                      WORK, LWORK, RWORK, RESULT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, LDB, LWORK, M, N, P
   DOUBLE PRECISION   RESULT
!     ..
!     .. Array Arguments ..
!
!  ====================================================================
!
   DOUBLE PRECISION   RWORK( * )
   COMPLEX*16         A( LDA, * ), AF( LDA, * ), B( LDB, * ), &
                      BF( LDB, * ), D( * ), DF( * ), U( * ), &
                      WORK( LWORK ), X( * )
!     ..
!     .. Local Scalars ..
   INTEGER            INFO
   DOUBLE PRECISION   ANORM, BNORM, DNORM, EPS, UNFL, XNORM, YNORM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH, DZASUM, ZLANGE
   EXTERNAL           DLAMCH, DZASUM, ZLANGE
!     ..
!     .. External Subroutines ..
!
   EXTERNAL           ZCOPY, ZGEMV, ZGGGLM, ZLACPY
!     ..
!     .. Executable Statements ..
!
   EPS = DLAMCH( 'Epsilon' )
   UNFL = DLAMCH( 'Safe minimum' )
   ANORM = MAX( ZLANGE( '1', N, M, A, LDA, RWORK ), UNFL )
   BNORM = MAX( ZLANGE( '1', N, P, B, LDB, RWORK ), UNFL )
!
!     Copy the matrices A and B to the arrays AF and BF,
!     and the vector D the array DF.
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'Full', N, M, A, LDA, AF, LDA )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'Full', N, P, B, LDB, BF, LDB )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZCOPY( N, D, 1, DF, 1 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZCOPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Solve GLM problem
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGGGLM( N, M, P, AF, LDA, BF, LDB, DF, X, U, WORK, LWORK, &
                INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGGGLM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Test the residual for the solution of LSE
!
!                       norm( d - A*x - B*u )
!       RESULT = -----------------------------------------
!                (norm(A)+norm(B))*(norm(x)+norm(u))*EPS
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZCOPY( N, D, 1, DF, 1 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZCOPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMV( 'No transpose', N, M, -(1.0D0,0.0D0), A, LDA, X, 1, (1.0D0,0.0D0), DF, &
               1 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMV : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMV( 'No transpose', N, P, -(1.0D0,0.0D0), B, LDB, U, 1, (1.0D0,0.0D0), DF, &
               1 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMV : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   DNORM = DZASUM( N, DF, 1 )
   XNORM = DZASUM( M, X, 1 ) + DZASUM( P, U, 1 )
   YNORM = ANORM + BNORM
!
   IF( XNORM <= 0.0D0 ) THEN
      RESULT = 0.0D0
   ELSE
      RESULT = ( ( DNORM / YNORM ) / XNORM ) / EPS
   END IF
!
   RETURN
!
!     End of ZGLMTS
!
END




