!> \brief \b ZGSVTS3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGSVTS3( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V,
!                           LDV, Q, LDQ, ALPHA, BETA, R, LDR, IWORK, WORK,
!                           LWORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LDQ, LDR, LDU, LDV, LWORK, M, N, P
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   ALPHA( * ), BETA( * ), RESULT( 6 ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), AF( LDA, * ), B( LDB, * ),
!      $                   BF( LDB, * ), Q( LDQ, * ), R( LDR, * ),
!      $                   U( LDU, * ), V( LDV, * ), WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGSVTS3 tests ZGGSVD3, which computes the GSVD of an M-by-N matrix A
!> and a P-by-N matrix B:
!>              U'*A*Q = D1*R and V'*B*Q = D2*R.
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
!>          A is COMPLEX*16 array, dimension (LDA,M)
!>          The M-by-N matrix A.
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is COMPLEX*16 array, dimension (LDA,N)
!>          Details of the GSVD of A and B, as returned by ZGGSVD3,
!>          see ZGGSVD3 for further details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A and AF.
!>          LDA >= max( 1,M ).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,P)
!>          On entry, the P-by-N matrix B.
!> \endverbatim
!>
!> \param[out] BF
!> \verbatim
!>          BF is COMPLEX*16 array, dimension (LDB,N)
!>          Details of the GSVD of A and B, as returned by ZGGSVD3,
!>          see ZGGSVD3 for further details.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the arrays B and BF.
!>          LDB >= max(1,P).
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX*16 array, dimension(LDU,M)
!>          The M by M unitary matrix U.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U. LDU >= max(1,M).
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension(LDV,M)
!>          The P by P unitary matrix V.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V. LDV >= max(1,P).
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension(LDQ,N)
!>          The N by N unitary matrix Q.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q. LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[out] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION array, dimension (N)
!>
!>          The generalized singular value pairs of A and B, the
!>          ``diagonal'' matrices D1 and D2 are constructed from
!>          ALPHA and BETA, see subroutine ZGGSVD3 for details.
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is COMPLEX*16 array, dimension(LDQ,N)
!>          The upper triangular matrix R.
!> \endverbatim
!>
!> \param[in] LDR
!> \verbatim
!>          LDR is INTEGER
!>          The leading dimension of the array R. LDR >= max(1,N).
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N)
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
!>          The dimension of the array WORK,
!>          LWORK >= max(M,P,N)*max(M,P,N).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (max(M,P,N))
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (6)
!>          The test ratios:
!>          RESULT(1) = norm( U'*A*Q - D1*R ) / ( MAX(M,N)*norm(A)*ULP)
!>          RESULT(2) = norm( V'*B*Q - D2*R ) / ( MAX(P,N)*norm(B)*ULP)
!>          RESULT(3) = norm( I - U'*U ) / ( M*ULP )
!>          RESULT(4) = norm( I - V'*V ) / ( P*ULP )
!>          RESULT(5) = norm( I - Q'*Q ) / ( N*ULP )
!>          RESULT(6) = 0        if ALPHA is in decreasing order;
!>                    = ULPINV   otherwise.
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
   SUBROUTINE ZGSVTS3( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V, &
                       LDV, Q, LDQ, ALPHA, BETA, R, LDR, IWORK, WORK, &
                       LWORK, RWORK, RESULT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, LDB, LDQ, LDR, LDU, LDV, LWORK, M, N, P
!     ..
!     .. Array Arguments ..
   INTEGER            IWORK( * )
   DOUBLE PRECISION   ALPHA( * ), BETA( * ), RESULT( 6 ), RWORK( * )
   COMPLEX*16         A( LDA, * ), AF( LDA, * ), B( LDB, * ), &
                      BF( LDB, * ), Q( LDQ, * ), R( LDR, * ), &
                      U( LDU, * ), V( LDV, * ), WORK( LWORK )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, INFO, J, K, L
   DOUBLE PRECISION   ANORM, BNORM, RESID, TEMP, ULP, ULPINV, UNFL
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH, ZLANGE, ZLANHE
   EXTERNAL           DLAMCH, ZLANGE, ZLANHE
!     ..
!     .. External Subroutines ..
   EXTERNAL           DCOPY, ZGEMM, ZGGSVD3, ZHERK, ZLACPY, ZLASET
!     ..
!     .. Executable Statements ..
!
   ULP = DLAMCH( 'Precision' )
   ULPINV = 1.0D0 / ULP
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
   CALL ZGGSVD3( 'U', 'V', 'Q', M, N, P, K, L, AF, LDA, BF, LDB, &
                 ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, LWORK, &
                 RWORK, IWORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGGSVD3 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Copy R
!
   DO I = 1, MIN( K+L, M )
      R(I,I:K+L) = AF(I,N-K-L+I:N)
   ENDDO
!
   IF( M-K-L < 0 ) THEN
      DO I = M + 1, K + L
         R(I,I:K+L) = BF(I-K,N-K-L+I:N)
      ENDDO
   END IF
!
!     Compute A:= U'*A*Q - D1*R
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'No transpose', 'No transpose', M, N, N, (1.0D0,0.0D0), A, LDA, &
               Q, LDQ, (0.0D+0,0.0D+0), WORK, LDA )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'Conjugate transpose', 'No transpose', M, N, M, (1.0D0,0.0D0), &
               U, LDU, WORK, LDA, (0.0D+0,0.0D+0), A, LDA )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   DO I = 1, K
      A(I,N-K-L+I:N) = A(I,N-K-L+I:N) - R(I,I:K+L)
   ENDDO
!
   DO I = K + 1, MIN( K+L, M )
      A(I,N-K-L+I:N-K-L+K+L) = A(I,N-K-L+I:N)-ALPHA(I)*R(I,I:K+L)
   ENDDO
!
!     Compute norm( U'*A*Q - D1*R ) / ( MAX(1,M,N)*norm(A)*ULP ) .
!
   RESID = ZLANGE( '1', M, N, A, LDA, RWORK )
   IF( ANORM > 0.0D0 ) THEN
      RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M, N ) ) ) / ANORM ) / &
                    ULP
   ELSE
      RESULT( 1 ) = 0.0D0
   END IF
!
!     Compute B := V'*B*Q - D2*R
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'No transpose', 'No transpose', P, N, N, (1.0D0,0.0D0), B, LDB, &
               Q, LDQ, (0.0D+0,0.0D+0), WORK, LDB )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'Conjugate transpose', 'No transpose', P, N, P, (1.0D0,0.0D0), &
               V, LDV, WORK, LDB, (0.0D+0,0.0D+0), B, LDB )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   DO I = 1, L
      B(I,N-L+I:N) = B(I,N-L+I:N) - BETA(K+I)*R(K+I,K+I:K+L)
   ENDDO
!
!     Compute norm( V'*B*Q - D2*R ) / ( MAX(P,N)*norm(B)*ULP ) .
!
   RESID = ZLANGE( '1', P, N, B, LDB, RWORK )
   IF( BNORM > 0.0D0 ) THEN
      RESULT( 2 ) = ( ( RESID / DBLE( MAX( 1, P, N ) ) ) / BNORM ) / &
                    ULP
   ELSE
      RESULT( 2 ) = 0.0D0
   END IF
!
!     Compute I - U'*U
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', M, M, (0.0D+0,0.0D+0), (1.0D0,0.0D0), WORK, LDQ )
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
   CALL ZHERK( 'Upper', 'Conjugate transpose', M, M, -1.0D0, U, LDU, &
               1.0D0, WORK, LDU )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZHERK : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute norm( I - U'*U ) / ( M * ULP ) .
!
   RESID = ZLANHE( '1', 'Upper', M, WORK, LDU, RWORK )
   RESULT( 3 ) = ( RESID / DBLE( MAX( 1, M ) ) ) / ULP
!
!     Compute I - V'*V
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', P, P, (0.0D+0,0.0D+0), (1.0D0,0.0D0), WORK, LDV )
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
   CALL ZHERK( 'Upper', 'Conjugate transpose', P, P, -1.0D0, V, LDV, &
               1.0D0, WORK, LDV )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZHERK : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute norm( I - V'*V ) / ( P * ULP ) .
!
   RESID = ZLANHE( '1', 'Upper', P, WORK, LDV, RWORK )
   RESULT( 4 ) = ( RESID / DBLE( MAX( 1, P ) ) ) / ULP
!
!     Compute I - Q'*Q
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', N, N, (0.0D+0,0.0D+0), (1.0D0,0.0D0), WORK, LDQ )
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
   CALL ZHERK( 'Upper', 'Conjugate transpose', N, N, -1.0D0, Q, LDQ, &
               1.0D0, WORK, LDQ )
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
   RESID = ZLANHE( '1', 'Upper', N, WORK, LDQ, RWORK )
   RESULT( 5 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / ULP
!
!     Check sorting
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DCOPY( N, ALPHA, 1, RWORK, 1 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   DO I = K + 1, MIN( K+L, M )
      J = IWORK( I )
      IF( I /= J ) THEN
         TEMP = RWORK( I )
         RWORK( I ) = RWORK( J )
         RWORK( J ) = TEMP
      END IF
   ENDDO
!
   RESULT( 6 ) = 0.0D0
   DO I = K + 1, MIN( K+L, M ) - 1
      IF( RWORK( I ) < RWORK( I+1 ) ) RESULT( 6 ) = ULPINV
   ENDDO
!
   RETURN
!
!     End of ZGSVTS3
!
END




