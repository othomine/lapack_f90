!> \brief \b CBDT03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CBDT03( UPLO, N, KD, D, E, U, LDU, S, VT, LDVT, WORK,
!                          RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            KD, LDU, LDVT, N
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * ), S( * )
!       COMPLEX            U( LDU, * ), VT( LDVT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CBDT03 reconstructs a bidiagonal matrix B from its SVD:
!>    S = U' * B * V
!> where U and V are orthogonal matrices and S is diagonal.
!>
!> The test ratio to test the singular value decomposition is
!>    RESID = norm( B - U * S * VT ) / ( n * norm(B) * EPS )
!> where VT = V' and EPS is the machine precision.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix B is upper or lower bidiagonal.
!>          = 'U':  Upper bidiagonal
!>          = 'L':  Lower bidiagonal
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix B.
!> \endverbatim
!>
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The bandwidth of the bidiagonal matrix B.  If KD = 1, the
!>          matrix B is bidiagonal, and if KD = 0, B is diagonal and E is
!>          not referenced.  If KD is greater than 1, it is assumed to be
!>          1, and if KD is less than 0, it is assumed to be 0.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The n diagonal elements of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          The (n-1) superdiagonal elements of the bidiagonal matrix B
!>          if UPLO = 'U', or the (n-1) subdiagonal elements of B if
!>          UPLO = 'L'.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is COMPLEX array, dimension (LDU,N)
!>          The n by n orthogonal matrix U in the reduction B = U'*A*P.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= max(1,N)
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is REAL array, dimension (N)
!>          The singular values from the SVD of B, sorted in decreasing
!>          order.
!> \endverbatim
!>
!> \param[in] VT
!> \verbatim
!>          VT is COMPLEX array, dimension (LDVT,N)
!>          The n by n orthogonal matrix V' in the reduction
!>          B = U * S * V'.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>          The leading dimension of the array VT.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
!>          The test ratio:  norm(B - U * S * V') / ( n * norm(A) * EPS )
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
   SUBROUTINE CBDT03( UPLO, N, KD, D, E, U, LDU, S, VT, LDVT, WORK, &
                      RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            KD, LDU, LDVT, N
   REAL               RESID
!     ..
!     .. Array Arguments ..
   REAL               D( * ), E( * ), S( * )
   COMPLEX            U( LDU, * ), VT( LDVT, * ), WORK( * )
!     ..
!
! ======================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, J
   REAL               BNORM, EPS
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ISAMAX
   REAL               SCASUM, SLAMCH
   EXTERNAL           LSAME, ISAMAX, SCASUM, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMV
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
   RESID = 0.0E+0
   IF( N <= 0 ) RETURN
!
!     Compute B - U * S * V' one column at a time.
!
   BNORM = 0.0E+0
   IF( KD >= 1 ) THEN
!
!        B is bidiagonal.
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!           B is upper bidiagonal.
!
         DO J = 1, N
            WORK(N+1:2*N) = S(1:N)*VT(1:N,J)
            CALL CGEMV( 'No transpose', N, N, -CMPLX( 1.0E+0 ), U, LDU, &
                        WORK( N+1 ), 1, CMPLX( 0.0E+0 ), WORK, 1 )
            WORK( J ) = WORK( J ) + D( J )
            IF( J > 1 ) THEN
               WORK( J-1 ) = WORK( J-1 ) + E( J-1 )
               BNORM = MAX( BNORM, ABS( D( J ) )+ABS( E( J-1 ) ) )
            ELSE
               BNORM = MAX( BNORM, ABS( D( J ) ) )
            END IF
            RESID = MAX( RESID, SCASUM( N, WORK, 1 ) )
         ENDDO
      ELSE
!
!           B is lower bidiagonal.
!
         DO J = 1, N
            WORK(N+1:2*N) = S(1:N)*VT(1:N,J)
            CALL CGEMV( 'No transpose', N, N, -CMPLX( 1.0E+0 ), U, LDU, &
                        WORK( N+1 ), 1, CMPLX( 0.0E+0 ), WORK, 1 )
            WORK( J ) = WORK( J ) + D( J )
            IF( J < N ) THEN
               WORK( J+1 ) = WORK( J+1 ) + E( J )
               BNORM = MAX( BNORM, ABS( D( J ) )+ABS( E( J ) ) )
            ELSE
               BNORM = MAX( BNORM, ABS( D( J ) ) )
            END IF
            RESID = MAX( RESID, SCASUM( N, WORK, 1 ) )
         ENDDO
      END IF
   ELSE
!
!        B is diagonal.
!
      DO J = 1, N
         WORK(N+1:2*N) = S(1:N)*VT(1:N,J)
         CALL CGEMV( 'No transpose', N, N, -CMPLX( 1.0E+0 ), U, LDU, &
                     WORK( N+1 ), 1, CMPLX( 0.0E+0 ), WORK, 1 )
         WORK( J ) = WORK( J ) + D( J )
         RESID = MAX( RESID, SCASUM( N, WORK, 1 ) )
      ENDDO
      J = ISAMAX( N, D, 1 )
      BNORM = ABS( D( J ) )
   END IF
!
!     Compute norm(B - U * S * V') / ( n * norm(B) * EPS )
!
   EPS = SLAMCH( 'Precision' )
!
   IF( BNORM <= 0.0E+0 ) THEN
      IF( RESID /= 0.0E+0 ) RESID = 1.0E+0 / EPS
   ELSE
      IF( BNORM >= RESID ) THEN
         RESID = ( RESID / BNORM ) / ( REAL( N )*EPS )
      ELSE
         IF( BNORM < 1.0E+0 ) THEN
            RESID = ( MIN( RESID, REAL( N )*BNORM ) / BNORM ) / ( REAL( N )*EPS )
         ELSE
            RESID = MIN( RESID / BNORM, REAL( N ) ) / ( REAL( N )*EPS )
         END IF
      END IF
   END IF
!
   RETURN
!
!     End of CBDT03
!
END
