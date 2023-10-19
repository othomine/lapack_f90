!> \brief \b DLQT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLQT01( M, N, A, AF, Q, L, LDA, TAU, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), AF( LDA, * ), L( LDA, * ),
!      $                   Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLQT01 tests DGELQF, which computes the LQ factorization of an m-by-n
!> matrix A, and partially tests DORGLQ which forms the n-by-n
!> orthogonal matrix Q.
!>
!> DLQT01 compares L with A*Q', and checks that Q is orthogonal.
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
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m-by-n matrix A.
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is DOUBLE PRECISION array, dimension (LDA,N)
!>          Details of the LQ factorization of A, as returned by DGELQF.
!>          See DGELQF for further details.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDA,N)
!>          The n-by-n orthogonal matrix Q.
!> \endverbatim
!>
!> \param[out] L
!> \verbatim
!>          L is DOUBLE PRECISION array, dimension (LDA,max(M,N))
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A, AF, Q and L.
!>          LDA >= max(M,N).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors, as returned
!>          by DGELQF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
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
!>          RWORK is DOUBLE PRECISION array, dimension (max(M,N))
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (2)
!>          The test ratios:
!>          RESULT(1) = norm( L - A*Q' ) / ( N * norm(A) * EPS )
!>          RESULT(2) = norm( I - Q*Q' ) / ( N * EPS )
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
   SUBROUTINE DLQT01( M, N, A, AF, Q, L, LDA, TAU, WORK, LWORK, &
                      RWORK, RESULT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   A( LDA, * ), AF( LDA, * ), L( LDA, * ), &
                      Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), &
                      WORK( LWORK )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
   DOUBLE PRECISION   ROGUE
   PARAMETER          ( ROGUE = -1.0D+10 )
!     ..
!     .. Local Scalars ..
   INTEGER            INFO, MINMN
   DOUBLE PRECISION   ANORM, EPS, RESID
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH, DLANGE, DLANSY
   EXTERNAL           DLAMCH, DLANGE, DLANSY
!     ..
!     .. External Subroutines ..
   EXTERNAL           DGELQF, DGEMM, DLACPY, DLASET, DORGLQ, DSYRK
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          DBLE, MAX, MIN
!     ..
!     .. Scalars in Common ..
   CHARACTER*32       SRNAMT
!     ..
!     .. Common blocks ..
   COMMON             / SRNAMC / SRNAMT
!     ..
!     .. Executable Statements ..
!
   MINMN = MIN( M, N )
   EPS = DLAMCH( 'Epsilon' )
!
!     Copy the matrix A to the array AF.
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLACPY( 'Full', M, N, A, LDA, AF, LDA )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Factorize the matrix A in the array AF.
!
   SRNAMT = 'DGELQF'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DGELQF( M, N, AF, LDA, TAU, WORK, LWORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DGELQF : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Copy details of Q
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLASET( 'Full', N, N, ROGUE, ROGUE, Q, LDA )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   IF( N > 1 )  THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLACPY( 'Upper', M, N-1, AF( 1, 2 ), LDA, Q( 1, 2 ), LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   ENDIF
!
!     Generate the n-by-n matrix Q
!
   SRNAMT = 'DORGLQ'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DORGLQ( N, N, MINMN, Q, LDA, TAU, WORK, LWORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DORGLQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Copy L
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLASET( 'Full', M, N, ZERO, ZERO, L, LDA )
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
   CALL DLACPY( 'Lower', M, N, AF, LDA, L, LDA )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute L - A*Q'
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DGEMM( 'No transpose', 'Transpose', M, N, N, -ONE, A, LDA, Q, &
               LDA, ONE, L, LDA )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute norm( L - Q'*A ) / ( N * norm(A) * EPS ) .
!
   ANORM = DLANGE( '1', M, N, A, LDA, RWORK )
   RESID = DLANGE( '1', M, N, L, LDA, RWORK )
   IF( ANORM > ZERO ) THEN
      RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, N ) ) ) / ANORM ) / EPS
   ELSE
      RESULT( 1 ) = ZERO
   END IF
!
!     Compute I - Q*Q'
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLASET( 'Full', N, N, ZERO, ONE, L, LDA )
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
   CALL DSYRK( 'Upper', 'No transpose', N, N, -ONE, Q, LDA, ONE, L, &
               LDA )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DSYRK : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute norm( I - Q*Q' ) / ( N * EPS ) .
!
   RESID = DLANSY( '1', 'Upper', N, L, LDA, RWORK )
!
   RESULT( 2 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / EPS
!
   RETURN
!
!     End of DLQT01
!
END
                                                                                                                                                                                                                                                                                                            




