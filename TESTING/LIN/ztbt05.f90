!> \brief \b ZTBT05
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTBT05( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B,
!                          LDB, X, LDX, XACT, LDXACT, FERR, BERR, RESLTS )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, TRANS, UPLO
!       INTEGER            KD, LDAB, LDB, LDX, LDXACT, N, NRHS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   BERR( * ), FERR( * ), RESLTS( * )
!       COMPLEX*16         AB( LDAB, * ), B( LDB, * ), X( LDX, * ),
!      $                   XACT( LDXACT, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTBT05 tests the error bounds from iterative refinement for the
!> computed solution to a system of equations A*X = B, where A is a
!> triangular band matrix.
!>
!> RESLTS(1) = test of the error bound
!>           = norm(X - XACT) / ( norm(X) * FERR )
!>
!> A large value is returned if this ratio is not less than one.
!>
!> RESLTS(2) = residual from the iterative refinement routine
!>           = the maximum of BERR / ( NZ*EPS + (*) ), where
!>             (*) = NZ*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )
!>             and NZ = max. number of nonzeros in any row of A, plus 1
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A is upper or lower triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the form of the system of equations.
!>          = 'N':  A * X = B  (No transpose)
!>          = 'T':  A'* X = B  (Transpose)
!>          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A is unit triangular.
!>          = 'N':  Non-unit triangular
!>          = 'U':  Unit triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows of the matrices X, B, and XACT, and the
!>          order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of super-diagonals of the matrix A if UPLO = 'U',
!>          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns of the matrices X, B, and XACT.
!>          NRHS >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is COMPLEX*16 array, dimension (LDAB,N)
!>          The upper or lower triangular band matrix A, stored in the
!>          first kd+1 rows of the array. The j-th column of A is stored
!>          in the j-th column of the array AB as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>          If DIAG = 'U', the diagonal elements of A are not referenced
!>          and are assumed to be 1.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>          The right hand side vectors for the system of linear
!>          equations.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (LDX,NRHS)
!>          The computed solution vectors.  Each vector is stored as a
!>          column of the matrix X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,N).
!> \endverbatim
!>
!> \param[in] XACT
!> \verbatim
!>          XACT is COMPLEX*16 array, dimension (LDX,NRHS)
!>          The exact solution vectors.  Each vector is stored as a
!>          column of the matrix XACT.
!> \endverbatim
!>
!> \param[in] LDXACT
!> \verbatim
!>          LDXACT is INTEGER
!>          The leading dimension of the array XACT.  LDXACT >= max(1,N).
!> \endverbatim
!>
!> \param[in] FERR
!> \verbatim
!>          FERR is DOUBLE PRECISION array, dimension (NRHS)
!>          The estimated forward error bounds for each solution vector
!>          X.  If XTRUE is the true solution, FERR bounds the magnitude
!>          of the largest entry in (X - XTRUE) divided by the magnitude
!>          of the largest entry in X.
!> \endverbatim
!>
!> \param[in] BERR
!> \verbatim
!>          BERR is DOUBLE PRECISION array, dimension (NRHS)
!>          The componentwise relative backward error of each solution
!>          vector (i.e., the smallest relative change in any entry of A
!>          or B that makes X an exact solution).
!> \endverbatim
!>
!> \param[out] RESLTS
!> \verbatim
!>          RESLTS is DOUBLE PRECISION array, dimension (2)
!>          The maximum over the NRHS solution vectors of the ratios:
!>          RESLTS(1) = norm(X - XACT) / ( norm(X) * FERR )
!>          RESLTS(2) = BERR / ( NZ*EPS + (*) )
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
!> \ingroup complex16_lin
!
!  =====================================================================
   SUBROUTINE ZTBT05( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, &
                      LDB, X, LDX, XACT, LDXACT, FERR, BERR, RESLTS )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          DIAG, TRANS, UPLO
   INTEGER            KD, LDAB, LDB, LDX, LDXACT, N, NRHS
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   BERR( * ), FERR( * ), RESLTS( * )
   COMPLEX*16         AB( LDAB, * ), B( LDB, * ), X( LDX, * ), &
                      XACT( LDXACT, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            NOTRAN, UNIT, UPPER
   INTEGER            I, IFU, IMAX, J, K, NZ
   DOUBLE PRECISION   AXBI, DIFF, EPS, ERRBND, OVFL, TMP, UNFL, XNORM
   COMPLEX*16         ZDUM
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            IZAMAX
   DOUBLE PRECISION   DLAMCH
   EXTERNAL           LSAME, IZAMAX, DLAMCH
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, DBLE, DIMAG, MAX, MIN
!     ..
!     .. Statement Functions ..
   DOUBLE PRECISION   CABS1
!     ..
!     .. Statement Function definitions ..
   CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0 or NRHS = 0.
!
   IF( N <= 0 .OR. NRHS <= 0 ) THEN
      RESLTS( 1 ) = ZERO
      RESLTS( 2 ) = ZERO
      RETURN
   END IF
!
   EPS = DLAMCH( 'Epsilon' )
   UNFL = DLAMCH( 'Safe minimum' )
   OVFL = ONE / UNFL
   UPPER = LSAME( UPLO, 'U' )
   NOTRAN = LSAME( TRANS, 'N' )
   UNIT = LSAME( DIAG, 'U' )
   NZ = MIN( KD, N-1 ) + 1
!
!     Test 1:  Compute the maximum of
!        norm(X - XACT) / ( norm(X) * FERR )
!     over all the vectors X and XACT using the infinity-norm.
!
   ERRBND = ZERO
   DO J = 1, NRHS
      IMAX = IZAMAX( N, X( 1, J ), 1 )
      XNORM = MAX( CABS1( X( IMAX, J ) ), UNFL )
      DIFF = ZERO
      DO I = 1, N
         DIFF = MAX( DIFF, CABS1( X( I, J )-XACT( I, J ) ) )
      ENDDO
!
      IF( XNORM > ONE ) THEN
         GO TO 20
      ELSE IF( DIFF <= OVFL*XNORM ) THEN
         GO TO 20
      ELSE
         ERRBND = ONE / EPS
         GO TO 30
      END IF
!
20    CONTINUE
      IF( DIFF / XNORM <= FERR( J ) ) THEN
         ERRBND = MAX( ERRBND, ( DIFF / XNORM ) / FERR( J ) )
      ELSE
         ERRBND = ONE / EPS
      END IF
30 CONTINUE
   ENDDO
   RESLTS( 1 ) = ERRBND
!
!     Test 2:  Compute the maximum of BERR / ( NZ*EPS + (*) ), where
!     (*) = NZ*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )
!
   IFU = 0
   IF( UNIT ) &
      IFU = 1
   DO K = 1, NRHS
      DO I = 1, N
         TMP = CABS1( B( I, K ) )
         IF( UPPER ) THEN
            IF( .NOT.NOTRAN ) THEN
               DO J = MAX( I-KD, 1 ), I - IFU
                  TMP = TMP + CABS1( AB( KD+1-I+J, I ) )* &
                        CABS1( X( J, K ) )
               ENDDO
               IF( UNIT ) &
                  TMP = TMP + CABS1( X( I, K ) )
            ELSE
               IF( UNIT ) &
                  TMP = TMP + CABS1( X( I, K ) )
               DO J = I + IFU, MIN( I+KD, N )
                  TMP = TMP + CABS1( AB( KD+1+I-J, J ) )* &
                        CABS1( X( J, K ) )
               ENDDO
            END IF
         ELSE
            IF( NOTRAN ) THEN
               DO J = MAX( I-KD, 1 ), I - IFU
                  TMP = TMP + CABS1( AB( 1+I-J, J ) )* &
                        CABS1( X( J, K ) )
               ENDDO
               IF( UNIT ) &
                  TMP = TMP + CABS1( X( I, K ) )
            ELSE
               IF( UNIT ) &
                  TMP = TMP + CABS1( X( I, K ) )
               DO J = I + IFU, MIN( I+KD, N )
                  TMP = TMP + CABS1( AB( 1+J-I, I ) )* &
                        CABS1( X( J, K ) )
               ENDDO
            END IF
         END IF
         IF( I == 1 ) THEN
            AXBI = TMP
         ELSE
            AXBI = MIN( AXBI, TMP )
         END IF
      ENDDO
      TMP = BERR( K ) / ( NZ*EPS+NZ*UNFL / MAX( AXBI, NZ*UNFL ) )
      IF( K == 1 ) THEN
         RESLTS( 2 ) = TMP
      ELSE
         RESLTS( 2 ) = MAX( RESLTS( 2 ), TMP )
      END IF
   ENDDO
!
   RETURN
!
!     End of ZTBT05
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
