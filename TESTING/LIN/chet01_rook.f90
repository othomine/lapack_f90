!> \brief \b CHET01_ROOK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHET01_ROOK( UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C, LDC,
!                               RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDAFAC, LDC, N
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               RWORK( * )
!       COMPLEX            A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHET01_ROOK reconstructs a complex Hermitian indefinite matrix A from its
!> block L*D*L' or U*D*U' factorization and computes the residual
!>    norm( C - A ) / ( N * norm(A) * EPS ),
!> where C is the reconstructed matrix, EPS is the machine epsilon,
!> L' is the transpose of L, and U' is the transpose of U.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          complex Hermitian matrix A is stored:
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows and columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The original complex Hermitian matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N)
!> \endverbatim
!>
!> \param[in] AFAC
!> \verbatim
!>          AFAC is COMPLEX array, dimension (LDAFAC,N)
!>          The factored form of the matrix A.  AFAC contains the block
!>          diagonal matrix D and the multipliers used to obtain the
!>          factor L or U from the block L*D*L' or U*D*U' factorization
!>          as computed by CSYTRF_ROOK.
!> \endverbatim
!>
!> \param[in] LDAFAC
!> \verbatim
!>          LDAFAC is INTEGER
!>          The leading dimension of the array AFAC.  LDAFAC >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices from CSYTRF_ROOK.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC,N)
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C.  LDC >= max(1,N).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
!>          If UPLO = 'L', norm(L*D*L' - A) / ( N * norm(A) * EPS )
!>          If UPLO = 'U', norm(U*D*U' - A) / ( N * norm(A) * EPS )
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
!> \ingroup complex_lin
!
!  =====================================================================
   SUBROUTINE CHET01_ROOK( UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C, &
                           LDC, RWORK, RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            LDA, LDAFAC, LDC, N
   REAL               RESID
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * )
   REAL               RWORK( * )
   COMPLEX            A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO, ONE
   PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
   COMPLEX            CZERO, CONE
   PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ), &
                      CONE = ( 1.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
   INTEGER            I, INFO, J
   REAL               ANORM, EPS
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               CLANHE, SLAMCH
   EXTERNAL           LSAME, CLANHE, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           CLASET, CLAVHE_ROOK
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          AIMAG, REAL
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0.
!
   IF( N <= 0 ) THEN
      RESID = ZERO
      RETURN
   END IF
!
!     Determine EPS and the norm of A.
!
   EPS = SLAMCH( 'Epsilon' )
   ANORM = CLANHE( '1', UPLO, N, A, LDA, RWORK )
!
!     Check the imaginary parts of the diagonal elements and return with
!     an error code if any are nonzero.
!
   DO J = 1, N
      IF( AIMAG( AFAC( J, J ) ) /= ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
   ENDDO
!
!     Initialize C to the identity matrix.
!
   CALL CLASET( 'Full', N, N, CZERO, CONE, C, LDC )
!
!     Call CLAVHE_ROOK to form the product D * U' (or D * L' ).
!
   CALL CLAVHE_ROOK( UPLO, 'Conjugate', 'Non-unit', N, N, AFAC, &
                     LDAFAC, IPIV, C, LDC, INFO )
!
!     Call CLAVHE_ROOK again to multiply by U (or L ).
!
   CALL CLAVHE_ROOK( UPLO, 'No transpose', 'Unit', N, N, AFAC, &
                     LDAFAC, IPIV, C, LDC, INFO )
!
!     Compute the difference  C - A .
!
   IF( LSAME( UPLO, 'U' ) ) THEN
      DO J = 1, N
         DO I = 1, J - 1
            C( I, J ) = C( I, J ) - A( I, J )
         ENDDO
         C( J, J ) = C( J, J ) - REAL( A( J, J ) )
      ENDDO
   ELSE
      DO J = 1, N
         C( J, J ) = C( J, J ) - REAL( A( J, J ) )
         DO I = J + 1, N
            C( I, J ) = C( I, J ) - A( I, J )
         ENDDO
      ENDDO
   END IF
!
!     Compute norm( C - A ) / ( N * norm(A) * EPS )
!
   RESID = CLANHE( '1', UPLO, N, C, LDC, RWORK )
!
   IF( ANORM <= ZERO ) THEN
      IF( RESID /= ZERO ) &
         RESID = ONE / EPS
   ELSE
      RESID = ( ( RESID/REAL( N ) )/ANORM ) / EPS
   END IF
!
   RETURN
!
!     End of CHET01_ROOK
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
