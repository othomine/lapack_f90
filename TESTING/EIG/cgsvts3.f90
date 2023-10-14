!> \brief \b CGSVTS3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGSVTS3( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V,
!                           LDV, Q, LDQ, ALPHA, BETA, R, LDR, IWORK, WORK,
!                           LWORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LDQ, LDR, LDU, LDV, LWORK, M, N, P
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               ALPHA( * ), BETA( * ), RESULT( 6 ), RWORK( * )
!       COMPLEX            A( LDA, * ), AF( LDA, * ), B( LDB, * ),
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
!> CGSVTS3 tests CGGSVD3, which computes the GSVD of an M-by-N matrix A
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
!>          A is COMPLEX array, dimension (LDA,M)
!>          The M-by-N matrix A.
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is COMPLEX array, dimension (LDA,N)
!>          Details of the GSVD of A and B, as returned by CGGSVD3,
!>          see CGGSVD3 for further details.
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
!>          B is COMPLEX array, dimension (LDB,P)
!>          On entry, the P-by-N matrix B.
!> \endverbatim
!>
!> \param[out] BF
!> \verbatim
!>          BF is COMPLEX array, dimension (LDB,N)
!>          Details of the GSVD of A and B, as returned by CGGSVD3,
!>          see CGGSVD3 for further details.
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
!>          U is COMPLEX array, dimension(LDU,M)
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
!>          V is COMPLEX array, dimension(LDV,M)
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
!>          Q is COMPLEX array, dimension(LDQ,N)
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
!>          ALPHA is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is REAL array, dimension (N)
!>
!>          The generalized singular value pairs of A and B, the
!>          ``diagonal'' matrices D1 and D2 are constructed from
!>          ALPHA and BETA, see subroutine CGGSVD3 for details.
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is COMPLEX array, dimension(LDQ,N)
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
!>          WORK is COMPLEX array, dimension (LWORK)
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
!>          RWORK is REAL array, dimension (max(M,P,N))
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (6)
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
!
!> \ingroup complex_eig
!
!  =====================================================================
   SUBROUTINE CGSVTS3( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V, &
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
   REAL               ALPHA( * ), BETA( * ), RESULT( 6 ), RWORK( * )
   COMPLEX            A( LDA, * ), AF( LDA, * ), B( LDB, * ), &
                      BF( LDB, * ), Q( LDQ, * ), R( LDR, * ), &
                      U( LDU, * ), V( LDV, * ), WORK( LWORK )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, INFO, J, K, L
   REAL               ANORM, BNORM, RESID, TEMP, ULP, ULPINV, UNFL
!     ..
!     .. External Functions ..
   REAL               CLANGE, CLANHE, SLAMCH
   EXTERNAL           CLANGE, CLANHE, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMM, CGGSVD3, CHERK, CLACPY, CLASET, SCOPY
!     ..
!     .. Executable Statements ..
!
   ULP = SLAMCH( 'Precision' )
   ULPINV = 1.0E+0 / ULP
   UNFL = SLAMCH( 'Safe minimum' )
!
!     Copy the matrix A to the array AF.
!
   CALL CLACPY( 'Full', M, N, A, LDA, AF, LDA )
   CALL CLACPY( 'Full', P, N, B, LDB, BF, LDB )
!
   ANORM = MAX( CLANGE( '1', M, N, A, LDA, RWORK ), UNFL )
   BNORM = MAX( CLANGE( '1', P, N, B, LDB, RWORK ), UNFL )
!
!     Factorize the matrices A and B in the arrays AF and BF.
!
   CALL CGGSVD3( 'U', 'V', 'Q', M, N, P, K, L, AF, LDA, BF, LDB, &
                 ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, LWORK, &
                 RWORK, IWORK, INFO )
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
   CALL CGEMM( 'No transpose', 'No transpose', M, N, N, (1.0E+0,0.0E+0), A, LDA, &
               Q, LDQ, (0.0E+0,0.0E+0), WORK, LDA )
!
   CALL CGEMM( 'Conjugate transpose', 'No transpose', M, N, M, (1.0E+0,0.0E+0), &
               U, LDU, WORK, LDA, (0.0E+0,0.0E+0), A, LDA )
!
   DO I = 1, K
      DO J = I, K + L
         A( I, N-K-L+J ) = A( I, N-K-L+J ) - R( I, J )
      ENDDO
   ENDDO
!
   DO I = K + 1, MIN( K+L, M )
      A(I,N-K-L+I:N-K-L+K+L) = A(I,N-K-L+I:N)-ALPHA(I)*R(I,I:K+L)
   ENDDO
!
!     Compute norm( U'*A*Q - D1*R ) / ( MAX(1,M,N)*norm(A)*ULP ) .
!
   RESID = CLANGE( '1', M, N, A, LDA, RWORK )
   IF( ANORM > 0.0E+0 ) THEN
      RESULT( 1 ) = ( ( RESID / REAL( MAX( 1, M, N ) ) ) / ANORM ) / ULP
   ELSE
      RESULT( 1 ) = 0.0E+0
   END IF
!
!     Compute B := V'*B*Q - D2*R
!
   CALL CGEMM( 'No transpose', 'No transpose', P, N, N, (1.0E+0,0.0E+0), B, LDB, &
               Q, LDQ, (0.0E+0,0.0E+0), WORK, LDB )
!
   CALL CGEMM( 'Conjugate transpose', 'No transpose', P, N, P, (1.0E+0,0.0E+0), &
               V, LDV, WORK, LDB, (0.0E+0,0.0E+0), B, LDB )
!
   DO I = 1, L
      B(I,N-L+I:N) = B(I,N-L+I:N) - BETA(K+I)*R(K+I,K+I:K+L)
   ENDDO
!
!     Compute norm( V'*B*Q - D2*R ) / ( MAX(P,N)*norm(B)*ULP ) .
!
   RESID = CLANGE( '1', P, N, B, LDB, RWORK )
   IF( BNORM > 0.0E+0 ) THEN
      RESULT( 2 ) = ( ( RESID / REAL( MAX( 1, P, N ) ) ) / BNORM ) / ULP
   ELSE
      RESULT( 2 ) = 0.0E+0
   END IF
!
!     Compute I - U'*U
!
   CALL CLASET( 'Full', M, M, (0.0E+0,0.0E+0), (1.0E+0,0.0E+0), WORK, LDQ )
   CALL CHERK( 'Upper', 'Conjugate transpose', M, M, -1.0E+0, U, LDU, 1.0E+0, WORK, LDU )
!
!     Compute norm( I - U'*U ) / ( M * ULP ) .
!
   RESID = CLANHE( '1', 'Upper', M, WORK, LDU, RWORK )
   RESULT( 3 ) = ( RESID / REAL( MAX( 1, M ) ) ) / ULP
!
!     Compute I - V'*V
!
   CALL CLASET( 'Full', P, P, (0.0E+0,0.0E+0), (1.0E+0,0.0E+0), WORK, LDV )
   CALL CHERK( 'Upper', 'Conjugate transpose', P, P, -1.0E+0, V, LDV, 1.0E+0, WORK, LDV )
!
!     Compute norm( I - V'*V ) / ( P * ULP ) .
!
   RESID = CLANHE( '1', 'Upper', P, WORK, LDV, RWORK )
   RESULT( 4 ) = ( RESID / REAL( MAX( 1, P ) ) ) / ULP
!
!     Compute I - Q'*Q
!
   CALL CLASET( 'Full', N, N, (0.0E+0,0.0E+0), (1.0E+0,0.0E+0), WORK, LDQ )
   CALL CHERK( 'Upper', 'Conjugate transpose', N, N, -1.0E+0, Q, LDQ, 1.0E+0, WORK, LDQ )
!
!     Compute norm( I - Q'*Q ) / ( N * ULP ) .
!
   RESID = CLANHE( '1', 'Upper', N, WORK, LDQ, RWORK )
   RESULT( 5 ) = ( RESID / REAL( MAX( 1, N ) ) ) / ULP
!
!     Check sorting
!
   CALL SCOPY( N, ALPHA, 1, RWORK, 1 )
   DO I = K + 1, MIN( K+L, M )
      J = IWORK( I )
      IF( I /= J ) THEN
         TEMP = RWORK( I )
         RWORK( I ) = RWORK( J )
         RWORK( J ) = TEMP
      END IF
   ENDDO
!
   RESULT( 6 ) = 0.0E+0
   DO I = K + 1, MIN( K+L, M ) - 1
      IF( RWORK( I ) < RWORK( I+1 ) ) RESULT( 6 ) = ULPINV
   ENDDO
!
   RETURN
!
!     End of CGSVTS3
!
END

