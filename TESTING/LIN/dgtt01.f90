!> \brief \b DGTT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGTT01( N, DL, D, DU, DLF, DF, DUF, DU2, IPIV, WORK,
!                          LDWORK, RWORK, RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            LDWORK, N
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   D( * ), DF( * ), DL( * ), DLF( * ), DU( * ),
!      $                   DU2( * ), DUF( * ), RWORK( * ),
!      $                   WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGTT01 reconstructs a tridiagonal matrix A from its LU factorization
!> and computes the residual
!>    norm(L*U - A) / ( norm(A) * EPS ),
!> where EPS is the machine epsilon.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] DL
!> \verbatim
!>          DL is DOUBLE PRECISION array, dimension (N-1)
!>          The (n-1) sub-diagonal elements of A.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The diagonal elements of A.
!> \endverbatim
!>
!> \param[in] DU
!> \verbatim
!>          DU is DOUBLE PRECISION array, dimension (N-1)
!>          The (n-1) super-diagonal elements of A.
!> \endverbatim
!>
!> \param[in] DLF
!> \verbatim
!>          DLF is DOUBLE PRECISION array, dimension (N-1)
!>          The (n-1) multipliers that define the matrix L from the
!>          LU factorization of A.
!> \endverbatim
!>
!> \param[in] DF
!> \verbatim
!>          DF is DOUBLE PRECISION array, dimension (N)
!>          The n diagonal elements of the upper triangular matrix U from
!>          the LU factorization of A.
!> \endverbatim
!>
!> \param[in] DUF
!> \verbatim
!>          DUF is DOUBLE PRECISION array, dimension (N-1)
!>          The (n-1) elements of the first super-diagonal of U.
!> \endverbatim
!>
!> \param[in] DU2
!> \verbatim
!>          DU2 is DOUBLE PRECISION array, dimension (N-2)
!>          The (n-2) elements of the second super-diagonal of U.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices; for 1 <= i <= n, row i of the matrix was
!>          interchanged with row IPIV(i).  IPIV(i) will always be either
!>          i or i+1; IPIV(i) = i indicates a row interchange was not
!>          required.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LDWORK,N)
!> \endverbatim
!>
!> \param[in] LDWORK
!> \verbatim
!>          LDWORK is INTEGER
!>          The leading dimension of the array WORK.  LDWORK >= max(1,N).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          The scaled residual:  norm(L*U - A) / (norm(A) * EPS)
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
!> \ingroup double_lin
!
!  =====================================================================
   SUBROUTINE DGTT01( N, DL, D, DU, DLF, DF, DUF, DU2, IPIV, WORK, &
                      LDWORK, RWORK, RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDWORK, N
   DOUBLE PRECISION   RESID
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * )
   DOUBLE PRECISION   D( * ), DF( * ), DL( * ), DLF( * ), DU( * ), &
                      DU2( * ), DUF( * ), RWORK( * ), &
                      WORK( LDWORK, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ONE, ZERO
   PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, IP, J, LASTJ
   DOUBLE PRECISION   ANORM, EPS, LI
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH, DLANGT, DLANHS
   EXTERNAL           DLAMCH, DLANGT, DLANHS
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MIN
!     ..
!     .. External Subroutines ..
   EXTERNAL           DAXPY, DSWAP
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
   IF( N <= 0 ) THEN
      RESID = ZERO
      RETURN
   END IF
!
   EPS = DLAMCH( 'Epsilon' )
!
!     Copy the matrix U to WORK.
!
   DO J = 1, N
      DO I = 1, N
         WORK( I, J ) = ZERO
      ENDDO
   ENDDO
   DO I = 1, N
      IF( I == 1 ) THEN
         WORK( I, I ) = DF( I )
         IF( N >= 2 ) &
            WORK( I, I+1 ) = DUF( I )
         IF( N >= 3 ) &
            WORK( I, I+2 ) = DU2( I )
      ELSE IF( I == N ) THEN
         WORK( I, I ) = DF( I )
      ELSE
         WORK( I, I ) = DF( I )
         WORK( I, I+1 ) = DUF( I )
         IF( I < N-1 ) &
            WORK( I, I+2 ) = DU2( I )
      END IF
   ENDDO
!
!     Multiply on the left by L.
!
   LASTJ = N
   DO I = N - 1, 1, -1
      LI = DLF( I )
      CALL DAXPY( LASTJ-I+1, LI, WORK( I, I ), LDWORK, &
                  WORK( I+1, I ), LDWORK )
      IP = IPIV( I )
      IF( IP == I ) THEN
         LASTJ = MIN( I+2, N )
      ELSE
         CALL DSWAP( LASTJ-I+1, WORK( I, I ), LDWORK, WORK( I+1, I ), &
                     LDWORK )
      END IF
   ENDDO
!
!     Subtract the matrix A.
!
   WORK( 1, 1 ) = WORK( 1, 1 ) - D( 1 )
   IF( N > 1 ) THEN
      WORK( 1, 2 ) = WORK( 1, 2 ) - DU( 1 )
      WORK( N, N-1 ) = WORK( N, N-1 ) - DL( N-1 )
      WORK( N, N ) = WORK( N, N ) - D( N )
      DO I = 2, N - 1
         WORK( I, I-1 ) = WORK( I, I-1 ) - DL( I-1 )
         WORK( I, I ) = WORK( I, I ) - D( I )
         WORK( I, I+1 ) = WORK( I, I+1 ) - DU( I )
      ENDDO
   END IF
!
!     Compute the 1-norm of the tridiagonal matrix A.
!
   ANORM = DLANGT( '1', N, DL, D, DU )
!
!     Compute the 1-norm of WORK, which is only guaranteed to be
!     upper Hessenberg.
!
   RESID = DLANHS( '1', N, WORK, LDWORK, RWORK )
!
!     Compute norm(L*U - A) / (norm(A) * EPS)
!
   IF( ANORM <= ZERO ) THEN
      IF( RESID /= ZERO ) &
         RESID = ONE / EPS
   ELSE
      RESID = ( RESID / ANORM ) / EPS
   END IF
!
   RETURN
!
!     End of DGTT01
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
