!> \brief \b ZGRQTS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGRQTS( M, P, N, A, AF, Q, R, LDA, TAUA, B, BF, Z, T,
!                          BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LWORK, M, N, P
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RESULT( 4 ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), AF( LDA, * ), B( LDB, * ),
!      $                   BF( LDB, * ), BWK( LDB, * ), Q( LDA, * ),
!      $                   R( LDA, * ), T( LDB, * ), TAUA( * ), TAUB( * ),
!      $                   WORK( LWORK ), Z( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGRQTS tests ZGGRQF, which computes the GRQ factorization of an
!> M-by-N matrix A and a P-by-N matrix B: A = R*Q and B = Z*T*Q.
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The M-by-N matrix A.
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is COMPLEX*16 array, dimension (LDA,N)
!>          Details of the GRQ factorization of A and B, as returned
!>          by ZGGRQF, see CGGRQF for further details.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDA,N)
!>          The N-by-N unitary matrix Q.
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is COMPLEX*16 array, dimension (LDA,MAX(M,N))
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A, AF, R and Q.
!>          LDA >= max(M,N).
!> \endverbatim
!>
!> \param[out] TAUA
!> \verbatim
!>          TAUA is COMPLEX*16 array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors, as returned
!>          by DGGQRC.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,N)
!>          On entry, the P-by-N matrix A.
!> \endverbatim
!>
!> \param[out] BF
!> \verbatim
!>          BF is COMPLEX*16 array, dimension (LDB,N)
!>          Details of the GQR factorization of A and B, as returned
!>          by ZGGRQF, see CGGRQF for further details.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDB,P)
!>          The P-by-P unitary matrix Z.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDB,max(P,N))
!> \endverbatim
!>
!> \param[out] BWK
!> \verbatim
!>          BWK is COMPLEX*16 array, dimension (LDB,N)
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the arrays B, BF, Z and T.
!>          LDB >= max(P,N).
!> \endverbatim
!>
!> \param[out] TAUB
!> \verbatim
!>          TAUB is COMPLEX*16 array, dimension (min(P,N))
!>          The scalar factors of the elementary reflectors, as returned
!>          by DGGRQF.
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
!>          The dimension of the array WORK, LWORK >= max(M,P,N)**2.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (4)
!>          The test ratios:
!>            RESULT(1) = norm( R - A*Q' ) / ( MAX(M,N)*norm(A)*ULP)
!>            RESULT(2) = norm( T*Q - Z'*B ) / (MAX(P,N)*norm(B)*ULP)
!>            RESULT(3) = norm( I - Q'*Q ) / ( N*ULP )
!>            RESULT(4) = norm( I - Z'*Z ) / ( P*ULP )
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
   SUBROUTINE ZGRQTS( M, P, N, A, AF, Q, R, LDA, TAUA, B, BF, Z, T, &
                      BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, LDB, LWORK, M, N, P
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   RESULT( 4 ), RWORK( * )
   COMPLEX*16         A( LDA, * ), AF( LDA, * ), B( LDB, * ), &
                      BF( LDB, * ), BWK( LDB, * ), Q( LDA, * ), &
                      R( LDA, * ), T( LDB, * ), TAUA( * ), TAUB( * ), &
                      WORK( LWORK ), Z( LDB, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   COMPLEX*16         CROGUE
   PARAMETER          ( CROGUE = ( -1.0D+10, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
   INTEGER            INFO
   DOUBLE PRECISION   ANORM, BNORM, RESID, ULP, UNFL
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH, ZLANGE, ZLANHE
   EXTERNAL           DLAMCH, ZLANGE, ZLANHE
!     ..
!     .. External Subroutines ..
   EXTERNAL           ZGEMM, ZGGRQF, ZHERK, ZLACPY, ZLASET, ZUNGQR, &
                      ZUNGRQ
!     ..
!     .. Executable Statements ..
!
   ULP = DLAMCH( 'Precision' )
   UNFL = DLAMCH( 'Safe minimum' )
!
!     Copy the matrix A to the array AF.
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'Full', M, N, A, LDA, AF, LDA )
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
   CALL ZLACPY( 'Full', P, N, B, LDB, BF, LDB )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   ANORM = MAX( ZLANGE( '1', M, N, A, LDA, RWORK ), UNFL )
   BNORM = MAX( ZLANGE( '1', P, N, B, LDB, RWORK ), UNFL )
!
!     Factorize the matrices A and B in the arrays AF and BF.
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGGRQF( M, P, N, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, &
                INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGGRQF : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Generate the N-by-N matrix Q
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', N, N, CROGUE, CROGUE, Q, LDA )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   IF( M <= N ) THEN
      IF( M > 0 .AND. M < N )  THEN
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZLACPY( 'Full', M, N-M, AF, LDA, Q( N-M+1, 1 ), LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDIF
      IF( M > 1 )  THEN
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZLACPY( 'Lower', M-1, M-1, AF( 2, N-M+1 ), LDA, &
                      Q( N-M+2, N-M+1 ), LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDIF
   ELSE
      IF( N > 1 )  THEN
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZLACPY( 'Lower', N-1, N-1, AF( M-N+2, 1 ), LDA, &
                      Q( 2, 1 ), LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDIF
   END IF
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGRQ( N, N, MIN( M, N ), Q, LDA, TAUA, WORK, LWORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGRQ : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Generate the P-by-P matrix Z
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', P, P, CROGUE, CROGUE, Z, LDB )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   IF( P > 1 )  THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLACPY( 'Lower', P-1, N, BF( 2, 1 ), LDB, Z( 2, 1 ), LDB )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   ENDIF
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNGQR( P, P, MIN( P, N ), Z, LDB, TAUB, WORK, LWORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNGQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Copy R
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', M, N, (0.0D+0,0.0D+0), (0.0D+0,0.0D+0), R, LDA )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   IF( M <= N ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLACPY( 'Upper', M, M, AF( 1, N-M+1 ), LDA, R( 1, N-M+1 ), &
                   LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   ELSE
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLACPY( 'Full', M-N, N, AF, LDA, R, LDA )
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
      CALL ZLACPY( 'Upper', N, N, AF( M-N+1, 1 ), LDA, R( M-N+1, 1 ), &
                   LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   END IF
!
!     Copy T
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', P, N, (0.0D+0,0.0D+0), (0.0D+0,0.0D+0), T, LDB )
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
   CALL ZLACPY( 'Upper', P, N, BF, LDB, T, LDB )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute R - A*Q'
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'No transpose', 'Conjugate transpose', M, N, N, -(1.0D0,0.0D0), &
               A, LDA, Q, LDA, (1.0D0,0.0D0), R, LDA )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute norm( R - A*Q' ) / ( MAX(M,N)*norm(A)*ULP ) .
!
   RESID = ZLANGE( '1', M, N, R, LDA, RWORK )
   IF( ANORM > 0.0D0 ) THEN
      RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M, N ) ) ) / ANORM ) / &
                    ULP
   ELSE
      RESULT( 1 ) = 0.0D0
   END IF
!
!     Compute T*Q - Z'*B
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'Conjugate transpose', 'No transpose', P, N, P, (1.0D0,0.0D0), &
               Z, LDB, B, LDB, (0.0D+0,0.0D+0), BWK, LDB )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'No transpose', 'No transpose', P, N, N, (1.0D0,0.0D0), T, LDB, &
               Q, LDA, -(1.0D0,0.0D0), BWK, LDB )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute norm( T*Q - Z'*B ) / ( MAX(P,N)*norm(A)*ULP ) .
!
   RESID = ZLANGE( '1', P, N, BWK, LDB, RWORK )
   IF( BNORM > 0.0D0 ) THEN
      RESULT( 2 ) = ( ( RESID / DBLE( MAX( 1, P, M ) ) ) / BNORM ) / &
                    ULP
   ELSE
      RESULT( 2 ) = 0.0D0
   END IF
!
!     Compute I - Q*Q'
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', N, N, (0.0D+0,0.0D+0), (1.0D0,0.0D0), R, LDA )
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
   CALL ZHERK( 'Upper', 'No Transpose', N, N, -1.0D0, Q, LDA, 1.0D0, R, &
               LDA )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZHERK : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute norm( I - Q'*Q ) / ( N * ULP ) .
!
   RESID = ZLANHE( '1', 'Upper', N, R, LDA, RWORK )
   RESULT( 3 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / ULP
!
!     Compute I - Z'*Z
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', P, P, (0.0D+0,0.0D+0), (1.0D0,0.0D0), T, LDB )
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
   CALL ZHERK( 'Upper', 'Conjugate transpose', P, P, -1.0D0, Z, LDB, &
               1.0D0, T, LDB )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZHERK : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute norm( I - Z'*Z ) / ( P*ULP ) .
!
   RESID = ZLANHE( '1', 'Upper', P, T, LDB, RWORK )
   RESULT( 4 ) = ( RESID / DBLE( MAX( 1, P ) ) ) / ULP
!
   RETURN
!
!     End of ZGRQTS
!
END




