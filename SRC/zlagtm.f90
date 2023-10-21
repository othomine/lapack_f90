!> \brief \b ZLAGTM performs a matrix-matrix product of the form C = αAB+βC, where A is a tridiagonal matrix, B and C are rectangular matrices, and α and β are scalars, which may be 0, 1, or -1.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAGTM + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlagtm.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlagtm.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlagtm.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAGTM( TRANS, N, NRHS, ALPHA, DL, D, DU, X, LDX, BETA,
!                          B, LDB )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            LDB, LDX, N, NRHS
!       DOUBLE PRECISION   ALPHA, BETA
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         B( LDB, * ), D( * ), DL( * ), DU( * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAGTM performs a matrix-matrix product of the form
!>
!>    B := alpha * A * X + beta * B
!>
!> where A is a tridiagonal matrix of order N, B and X are N by NRHS
!> matrices, and alpha and beta are real scalars, each of which may be
!> 0., 1., or -1.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the operation applied to A.
!>          = 'N':  No transpose, B := alpha * A * X + beta * B
!>          = 'T':  Transpose,    B := alpha * A**T * X + beta * B
!>          = 'C':  Conjugate transpose, B := alpha * A**H * X + beta * B
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
!>          The scalar alpha.  ALPHA must be 0., 1., or -1.; otherwise,
!>          it is assumed to be 0.
!> \endverbatim
!>
!> \param[in] DL
!> \verbatim
!>          DL is COMPLEX*16 array, dimension (N-1)
!>          The (n-1) sub-diagonal elements of T.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is COMPLEX*16 array, dimension (N)
!>          The diagonal elements of T.
!> \endverbatim
!>
!> \param[in] DU
!> \verbatim
!>          DU is COMPLEX*16 array, dimension (N-1)
!>          The (n-1) super-diagonal elements of T.
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
!> \ingroup lagtm
!
!  =====================================================================
   SUBROUTINE ZLAGTM( TRANS, N, NRHS, ALPHA, DL, D, DU, X, LDX, BETA, &
                      B, LDB )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          TRANS
   INTEGER            LDB, LDX, N, NRHS
   DOUBLE PRECISION   ALPHA, BETA
!     ..
!     .. Array Arguments ..
   COMPLEX*16         B( LDB, * ), D( * ), DL( * ), DU( * ), &
                      X( LDX, * )
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
!     Multiply B by BETA if BETA /= 1.
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
      IF( LSAME( TRANS, 'N' ) ) THEN
!
!           Compute B := B + A*X
!
         DO J = 1, NRHS
            IF( N == 1 ) THEN
               B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J )
            ELSE
               B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J ) + &
                           DU( 1 )*X( 2, J )
               B( N, J ) = B( N, J ) + DL( N-1 )*X( N-1, J ) + &
                           D( N )*X( N, J )
               DO I = 2, N - 1
                  B( I, J ) = B( I, J ) + DL( I-1 )*X( I-1, J ) + &
                              D( I )*X( I, J ) + DU( I )*X( I+1, J )
               ENDDO
            END IF
         ENDDO
      ELSE IF( LSAME( TRANS, 'T' ) ) THEN
!
!           Compute B := B + A**T * X
!
         DO J = 1, NRHS
            IF( N == 1 ) THEN
               B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J )
            ELSE
               B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J ) + &
                           DL( 1 )*X( 2, J )
               B( N, J ) = B( N, J ) + DU( N-1 )*X( N-1, J ) + &
                           D( N )*X( N, J )
               DO I = 2, N - 1
                  B( I, J ) = B( I, J ) + DU( I-1 )*X( I-1, J ) + &
                              D( I )*X( I, J ) + DL( I )*X( I+1, J )
               ENDDO
            END IF
         ENDDO
      ELSE IF( LSAME( TRANS, 'C' ) ) THEN
!
!           Compute B := B + A**H * X
!
         DO J = 1, NRHS
            IF( N == 1 ) THEN
               B( 1, J ) = B( 1, J ) + DCONJG( D( 1 ) )*X( 1, J )
            ELSE
               B( 1, J ) = B( 1, J ) + DCONJG( D( 1 ) )*X( 1, J ) + &
                           DCONJG( DL( 1 ) )*X( 2, J )
               B( N, J ) = B( N, J ) + DCONJG( DU( N-1 ) )* &
                           X( N-1, J ) + DCONJG( D( N ) )*X( N, J )
               DO I = 2, N - 1
                  B( I, J ) = B( I, J ) + DCONJG( DU( I-1 ) )* &
                              X( I-1, J ) + DCONJG( D( I ) )* &
                              X( I, J ) + DCONJG( DL( I ) )* &
                              X( I+1, J )
               ENDDO
            END IF
            ENDDO
      END IF
   ELSE IF( ALPHA == -ONE ) THEN
      IF( LSAME( TRANS, 'N' ) ) THEN
!
!           Compute B := B - A*X
!
         DO J = 1, NRHS
            IF( N == 1 ) THEN
               B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J )
            ELSE
               B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J ) - &
                           DU( 1 )*X( 2, J )
               B( N, J ) = B( N, J ) - DL( N-1 )*X( N-1, J ) - &
                           D( N )*X( N, J )
               DO I = 2, N - 1
                  B( I, J ) = B( I, J ) - DL( I-1 )*X( I-1, J ) - &
                              D( I )*X( I, J ) - DU( I )*X( I+1, J )
                  ENDDO
            END IF
            ENDDO
      ELSE IF( LSAME( TRANS, 'T' ) ) THEN
!
!           Compute B := B - A**T *X
!
         DO J = 1, NRHS
            IF( N == 1 ) THEN
               B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J )
            ELSE
               B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J ) - &
                           DL( 1 )*X( 2, J )
               B( N, J ) = B( N, J ) - DU( N-1 )*X( N-1, J ) - &
                           D( N )*X( N, J )
               DO I = 2, N - 1
                  B( I, J ) = B( I, J ) - DU( I-1 )*X( I-1, J ) - &
                              D( I )*X( I, J ) - DL( I )*X( I+1, J )
                  ENDDO
            END IF
            ENDDO
      ELSE IF( LSAME( TRANS, 'C' ) ) THEN
!
!           Compute B := B - A**H *X
!
         DO J = 1, NRHS
            IF( N == 1 ) THEN
               B( 1, J ) = B( 1, J ) - DCONJG( D( 1 ) )*X( 1, J )
            ELSE
               B( 1, J ) = B( 1, J ) - DCONJG( D( 1 ) )*X( 1, J ) - &
                           DCONJG( DL( 1 ) )*X( 2, J )
               B( N, J ) = B( N, J ) - DCONJG( DU( N-1 ) )* &
                           X( N-1, J ) - DCONJG( D( N ) )*X( N, J )
               DO I = 2, N - 1
                  B( I, J ) = B( I, J ) - DCONJG( DU( I-1 ) )* &
                              X( I-1, J ) - DCONJG( D( I ) )* &
                              X( I, J ) - DCONJG( DL( I ) )* &
                              X( I+1, J )
                  ENDDO
            END IF
            ENDDO
      END IF
   END IF
   RETURN
!
!     End of ZLAGTM
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

