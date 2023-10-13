!> \brief \b CBDT02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CBDT02( M, N, B, LDB, C, LDC, U, LDU, WORK, RWORK,
!                          RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            LDB, LDC, LDU, M, N
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            B( LDB, * ), C( LDC, * ), U( LDU, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CBDT02 tests the change of basis C = U**H * B by computing the
!> residual
!>
!>    RESID = norm(B - U * C) / ( max(m,n) * norm(B) * EPS ),
!>
!> where B and C are M by N matrices, U is an M by M orthogonal matrix,
!> and EPS is the machine precision.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrices B and C and the order of
!>          the matrix Q.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrices B and C.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,N)
!>          The m by n matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,M).
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC,N)
!>          The m by n matrix C, assumed to contain U**H * B.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C.  LDC >= max(1,M).
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is COMPLEX array, dimension (LDU,M)
!>          The m by m orthogonal matrix U.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (M)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
!>          RESID = norm(B - U * C) / ( max(m,n) * norm(B) * EPS ),
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
   SUBROUTINE CBDT02( M, N, B, LDB, C, LDC, U, LDU, WORK, RWORK, &
                      RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDB, LDC, LDU, M, N
   REAL               RESID
!     ..
!     .. Array Arguments ..
   REAL               RWORK( * )
   COMPLEX            B( LDB, * ), C( LDC, * ), U( LDU, * ), &
                      WORK( * )
!     ..
!
! ======================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            J
   REAL               BNORM, EPS, REALMN
!     ..
!     .. External Functions ..
   REAL               CLANGE, SCASUM, SLAMCH
   EXTERNAL           CLANGE, SCASUM, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           CCOPY, CGEMV
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
   RESID = 0.0E+0
   IF( M <= 0 .OR. N <= 0 ) &
      RETURN
   REALMN = REAL( MAX( M, N ) )
   EPS = SLAMCH( 'Precision' )
!
!     Compute norm(B - U * C)
!
   DO J = 1, N
      CALL CCOPY( M, B( 1, J ), 1, WORK, 1 )
      CALL CGEMV( 'No transpose', M, M, -CMPLX( 1.0E+0 ), U, LDU, &
                  C( 1, J ), 1, CMPLX( 1.0E+0 ), WORK, 1 )
      RESID = MAX( RESID, SCASUM( M, WORK, 1 ) )
   ENDDO
!
!     Compute norm of B.
!
   BNORM = CLANGE( '1', M, N, B, LDB, RWORK )
!
   IF( BNORM <= 0.0E+0 ) THEN
      IF( RESID /= 0.0E+0 ) &
         RESID = 1.0E+0 / EPS
   ELSE
      IF( BNORM >= RESID ) THEN
         RESID = ( RESID / BNORM ) / ( REALMN*EPS )
      ELSE
         IF( BNORM < 1.0E+0 ) THEN
            RESID = ( MIN( RESID, REALMN*BNORM ) / BNORM ) / &
                    ( REALMN*EPS )
         ELSE
            RESID = MIN( RESID / BNORM, REALMN ) / ( REALMN*EPS )
         END IF
      END IF
   END IF
   RETURN
!
!     End of CBDT02
!
END
