!> \brief \b ZLAPTM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAPTM( UPLO, N, NRHS, ALPHA, D, E, X, LDX, BETA, B,
!                          LDB )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDB, LDX, N, NRHS
!       DOUBLE PRECISION   ALPHA, BETA
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * )
!       COMPLEX*16         B( LDB, * ), E( * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAPTM multiplies an N by NRHS matrix X by a Hermitian tridiagonal
!> matrix A and stores the result in a matrix B.  The operation has the
!> form
!>
!>    B := alpha * A * X + beta * B
!>
!> where alpha may be either 1. or -1. and beta may be 0., 1., or -1.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER
!>          Specifies whether the superdiagonal or the subdiagonal of the
!>          tridiagonal matrix A is stored.
!>          = 'U':  Upper, E is the superdiagonal of A.
!>          = 'L':  Lower, E is the subdiagonal of A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrices X and B.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION
!>          The scalar alpha.  ALPHA must be 1. or -1.; otherwise,
!>          it is assumed to be 0.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The n diagonal elements of the tridiagonal matrix A.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is COMPLEX*16 array, dimension (N-1)
!>          The (n-1) subdiagonal or superdiagonal elements of A.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (LDX,NRHS)
!>          The N by NRHS matrix X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(N,1).
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION
!>          The scalar beta.  BETA must be 0., 1., or -1.; otherwise,
!>          it is assumed to be 1.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>          On entry, the N by NRHS matrix B.
!>          On exit, B is overwritten by the matrix expression
!>          B := alpha * A * X + beta * B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(N,1).
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
!> \ingroup complex16_lin
!
!  =====================================================================
   SUBROUTINE ZLAPTM( UPLO, N, NRHS, ALPHA, D, E, X, LDX, BETA, B, &
                      LDB )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            LDB, LDX, N, NRHS
   DOUBLE PRECISION   ALPHA, BETA
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   D( * )
   COMPLEX*16         B( LDB, * ), E( * ), X( LDX, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ONE, ZERO
   PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, J
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          DCONJG
!     ..
!     .. Executable Statements ..
!
   IF( N == 0 ) &
      RETURN
!
   IF( BETA == ZERO ) THEN
      DO J = 1, NRHS
         DO I = 1, N
            B( I, J ) = ZERO
         ENDDO
      ENDDO
   ELSE IF( BETA == -ONE ) THEN
      DO J = 1, NRHS
         DO I = 1, N
            B( I, J ) = -B( I, J )
         ENDDO
      ENDDO
   END IF
!
   IF( ALPHA == ONE ) THEN
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!           Compute B := B + A*X, where E is the superdiagonal of A.
!
         DO J = 1, NRHS
            IF( N == 1 ) THEN
               B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J )
            ELSE
               B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J ) + &
                           E( 1 )*X( 2, J )
               B( N, J ) = B( N, J ) + DCONJG( E( N-1 ) )* &
                           X( N-1, J ) + D( N )*X( N, J )
               DO I = 2, N - 1
                  B( I, J ) = B( I, J ) + DCONJG( E( I-1 ) )* &
                              X( I-1, J ) + D( I )*X( I, J ) + &
                              E( I )*X( I+1, J )
               ENDDO
            END IF
         ENDDO
      ELSE
!
!           Compute B := B + A*X, where E is the subdiagonal of A.
!
         DO J = 1, NRHS
            IF( N == 1 ) THEN
               B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J )
            ELSE
               B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J ) + &
                           DCONJG( E( 1 ) )*X( 2, J )
               B( N, J ) = B( N, J ) + E( N-1 )*X( N-1, J ) + &
                           D( N )*X( N, J )
               DO I = 2, N - 1
                  B( I, J ) = B( I, J ) + E( I-1 )*X( I-1, J ) + &
                              D( I )*X( I, J ) + &
                              DCONJG( E( I ) )*X( I+1, J )
               ENDDO
            END IF
         ENDDO
      END IF
   ELSE IF( ALPHA == -ONE ) THEN
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!           Compute B := B - A*X, where E is the superdiagonal of A.
!
         DO J = 1, NRHS
            IF( N == 1 ) THEN
               B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J )
            ELSE
               B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J ) - &
                           E( 1 )*X( 2, J )
               B( N, J ) = B( N, J ) - DCONJG( E( N-1 ) )* &
                           X( N-1, J ) - D( N )*X( N, J )
               DO I = 2, N - 1
                  B( I, J ) = B( I, J ) - DCONJG( E( I-1 ) )* &
                              X( I-1, J ) - D( I )*X( I, J ) - &
                              E( I )*X( I+1, J )
               ENDDO
            END IF
            ENDDO
      ELSE
!
!           Compute B := B - A*X, where E is the subdiagonal of A.
!
         DO J = 1, NRHS
            IF( N == 1 ) THEN
               B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J )
            ELSE
               B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J ) - &
                           DCONJG( E( 1 ) )*X( 2, J )
               B( N, J ) = B( N, J ) - E( N-1 )*X( N-1, J ) - &
                           D( N )*X( N, J )
               DO I = 2, N - 1
                  B( I, J ) = B( I, J ) - E( I-1 )*X( I-1, J ) - &
                              D( I )*X( I, J ) - &
                              DCONJG( E( I ) )*X( I+1, J )
                  ENDDO
            END IF
            ENDDO
      END IF
   END IF
   RETURN
!
!     End of ZLAPTM
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        



