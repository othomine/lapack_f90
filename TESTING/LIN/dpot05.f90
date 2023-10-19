!> \brief \b DPOT05
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DPOT05( UPLO, N, NRHS, A, LDA, B, LDB, X, LDX, XACT,
!                          LDXACT, FERR, BERR, RESLTS )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDB, LDX, LDXACT, N, NRHS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), BERR( * ), FERR( * ),
!      $                   RESLTS( * ), X( LDX, * ), XACT( LDXACT, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DPOT05 tests the error bounds from iterative refinement for the
!> computed solution to a system of equations A*X = B, where A is a
!> symmetric n by n matrix.
!>
!> RESLTS(1) = test of the error bound
!>           = norm(X - XACT) / ( norm(X) * FERR )
!>
!> A large value is returned if this ratio is not less than one.
!>
!> RESLTS(2) = residual from the iterative refinement routine
!>           = the maximum of BERR / ( (n+1)*EPS + (*) ), where
!>             (*) = (n+1)*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix A is stored.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows of the matrices X, B, and XACT, and the
!>          order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns of the matrices X, B, and XACT.
!>          NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The symmetric matrix A.  If UPLO = 'U', the leading n by n
!>          upper triangular part of A contains the upper triangular part
!>          of the matrix A, and the strictly lower triangular part of A
!>          is not referenced.  If UPLO = 'L', the leading n by n lower
!>          triangular part of A contains the lower triangular part of
!>          the matrix A, and the strictly upper triangular part of A is
!>          not referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
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
!>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)
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
!>          XACT is DOUBLE PRECISION array, dimension (LDX,NRHS)
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
!>          RESLTS(2) = BERR / ( (n+1)*EPS + (*) )
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
!> \ingroup double_lin
!
!  =====================================================================
   SUBROUTINE DPOT05( UPLO, N, NRHS, A, LDA, B, LDB, X, LDX, XACT, &
                      LDXACT, FERR, BERR, RESLTS )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            LDA, LDB, LDX, LDXACT, N, NRHS
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), BERR( * ), FERR( * ), &
                      RESLTS( * ), X( LDX, * ), XACT( LDXACT, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            UPPER
   INTEGER            I, IMAX, J, K
   DOUBLE PRECISION   AXBI, DIFF, EPS, ERRBND, OVFL, TMP, UNFL, XNORM
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            IDAMAX
   DOUBLE PRECISION   DLAMCH
   EXTERNAL           LSAME, IDAMAX, DLAMCH
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, MAX, MIN
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
!
!     Test 1:  Compute the maximum of
!        norm(X - XACT) / ( norm(X) * FERR )
!     over all the vectors X and XACT using the infinity-norm.
!
   ERRBND = ZERO
   DO J = 1, NRHS
      IMAX = IDAMAX( N, X( 1, J ), 1 )
      XNORM = MAX( ABS( X( IMAX, J ) ), UNFL )
      DIFF = ZERO
      DO I = 1, N
         DIFF = MAX( DIFF, ABS( X( I, J )-XACT( I, J ) ) )
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
!     Test 2:  Compute the maximum of BERR / ( (n+1)*EPS + (*) ), where
!     (*) = (n+1)*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )
!
   DO K = 1, NRHS
      DO I = 1, N
         TMP = ABS( B( I, K ) )
         IF( UPPER ) THEN
            DO J = 1, I
               TMP = TMP + ABS( A( J, I ) )*ABS( X( J, K ) )
            ENDDO
            DO J = I + 1, N
               TMP = TMP + ABS( A( I, J ) )*ABS( X( J, K ) )
            ENDDO
         ELSE
            DO J = 1, I - 1
               TMP = TMP + ABS( A( I, J ) )*ABS( X( J, K ) )
            ENDDO
            DO J = I, N
               TMP = TMP + ABS( A( J, I ) )*ABS( X( J, K ) )
            ENDDO
         END IF
         IF( I == 1 ) THEN
            AXBI = TMP
         ELSE
            AXBI = MIN( AXBI, TMP )
         END IF
      ENDDO
      TMP = BERR( K ) / ( ( N+1 )*EPS+( N+1 )*UNFL / &
            MAX( AXBI, ( N+1 )*UNFL ) )
      IF( K == 1 ) THEN
         RESLTS( 2 ) = TMP
      ELSE
         RESLTS( 2 ) = MAX( RESLTS( 2 ), TMP )
      END IF
   ENDDO
!
   RETURN
!
!     End of DPOT05
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        



