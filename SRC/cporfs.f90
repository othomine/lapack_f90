!> \brief \b CPORFS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CPORFS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cporfs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cporfs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cporfs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X,
!                          LDX, FERR, BERR, WORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
!       ..
!       .. Array Arguments ..
!       REAL               BERR( * ), FERR( * ), RWORK( * )
!       COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
!      $                   WORK( * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CPORFS improves the computed solution to a system of linear
!> equations when the coefficient matrix is Hermitian positive definite,
!> and provides error bounds and backward error estimates for the
!> solution.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
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
!>          of the matrices B and X.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The Hermitian matrix A.  If UPLO = 'U', the leading N-by-N
!>          upper triangular part of A contains the upper triangular part
!>          of the matrix A, and the strictly lower triangular part of A
!>          is not referenced.  If UPLO = 'L', the leading N-by-N lower
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
!> \param[in] AF
!> \verbatim
!>          AF is COMPLEX array, dimension (LDAF,N)
!>          The triangular factor U or L from the Cholesky factorization
!>          A = U**H*U or A = L*L**H, as computed by CPOTRF.
!> \endverbatim
!>
!> \param[in] LDAF
!> \verbatim
!>          LDAF is INTEGER
!>          The leading dimension of the array AF.  LDAF >= max(1,N).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,NRHS)
!>          The right hand side matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX array, dimension (LDX,NRHS)
!>          On entry, the solution matrix X, as computed by CPOTRS.
!>          On exit, the improved solution matrix X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,N).
!> \endverbatim
!>
!> \param[out] FERR
!> \verbatim
!>          FERR is REAL array, dimension (NRHS)
!>          The estimated forward error bound for each solution vector
!>          X(j) (the j-th column of the solution matrix X).
!>          If XTRUE is the true solution corresponding to X(j), FERR(j)
!>          is an estimated upper bound for the magnitude of the largest
!>          element in (X(j) - XTRUE) divided by the magnitude of the
!>          largest element in X(j).  The estimate is as reliable as
!>          the estimate for RCOND, and is almost always a slight
!>          overestimate of the true error.
!> \endverbatim
!>
!> \param[out] BERR
!> \verbatim
!>          BERR is REAL array, dimension (NRHS)
!>          The componentwise relative backward error of each solution
!>          vector X(j) (i.e., the smallest relative change in
!>          any element of A or B that makes X(j) an exact solution).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!> \par Internal Parameters:
!  =========================
!>
!> \verbatim
!>  ITMAX is the maximum number of steps of iterative refinement.
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
!> \ingroup porfs
!
!  =====================================================================
   SUBROUTINE CPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X, &
                      LDX, FERR, BERR, WORK, RWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
!     ..
!     .. Array Arguments ..
   REAL               BERR( * ), FERR( * ), RWORK( * )
   COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ), &
                      WORK( * ), X( LDX, * )
!     ..
!
!  ====================================================================
!
!     .. Parameters ..
   INTEGER            ITMAX
   PARAMETER          ( ITMAX = 5 )
!     ..
!     .. Local Scalars ..
   LOGICAL            UPPER
   INTEGER            COUNT, I, J, K, KASE, NZ
   REAL               EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK
   COMPLEX            ZDUM
!     ..
!     .. Local Arrays ..
   INTEGER            ISAVE( 3 )
!     ..
!     .. External Subroutines ..
   EXTERNAL           CHEMV, CLACN2, CPOTRS, XERBLA
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               SLAMCH
   EXTERNAL           LSAME, SLAMCH
!     ..
!     .. Statement Functions ..
   REAL               CABS1
!     ..
!     .. Statement Function definitions ..
   CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   UPPER = LSAME( UPLO, 'U' )
   IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( NRHS < 0 ) THEN
      INFO = -3
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -5
   ELSE IF( LDAF < MAX( 1, N ) ) THEN
      INFO = -7
   ELSE IF( LDB < MAX( 1, N ) ) THEN
      INFO = -9
   ELSE IF( LDX < MAX( 1, N ) ) THEN
      INFO = -11
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CPORFS', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 .OR. NRHS == 0 ) THEN
      FERR(1:NRHS) = 0.0E+0
      BERR(1:NRHS) = 0.0E+0
      RETURN
   END IF
!
!     NZ = maximum number of nonzero elements in each row of A, plus 1
!
   NZ = N + 1
   EPS = SLAMCH( 'Epsilon' )
   SAFMIN = SLAMCH( 'Safe minimum' )
   SAFE1 = NZ*SAFMIN
   SAFE2 = SAFE1 / EPS
!
!     Do for each right hand side
!
   DO J = 1, NRHS
!
      COUNT = 1
      LSTRES = 3.0E+0
20    CONTINUE
!
!        Loop until stopping criterion is satisfied.
!
!        Compute residual R = B - A * X
!
      WORK(1:N) = B(1:N,J)
      CALL CHEMV( UPLO, N, -(1.0E+0,0.0E+0), A, LDA, X( 1, J ), 1, (1.0E+0,0.0E+0), WORK, 1 )
!
!        Compute componentwise relative backward error from formula
!
!        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )
!
!        where abs(Z) is the componentwise absolute value of the matrix
!        or vector Z.  If the i-th component of the denominator is less
!        than SAFE2, then SAFE1 is added to the i-th components of the
!        numerator and denominator before dividing.
!
      DO I = 1, N
         RWORK( I ) = CABS1( B( I, J ) )
      ENDDO
!
!        Compute abs(A)*abs(X) + abs(B).
!
      IF( UPPER ) THEN
         DO K = 1, N
            S = 0.0E+0
            XK = CABS1( X( K, J ) )
            DO I = 1, K - 1
               RWORK( I ) = RWORK( I ) + CABS1( A( I, K ) )*XK
               S = S + CABS1( A( I, K ) )*CABS1( X( I, J ) )
            ENDDO
            RWORK( K ) = RWORK( K ) + ABS( REAL( A( K, K ) ) )*XK + S
         ENDDO
      ELSE
         DO K = 1, N
            S = 0.0E+0
            XK = CABS1( X( K, J ) )
            RWORK( K ) = RWORK( K ) + ABS( REAL( A( K, K ) ) )*XK
            DO I = K + 1, N
               RWORK( I ) = RWORK( I ) + CABS1( A( I, K ) )*XK
               S = S + CABS1( A( I, K ) )*CABS1( X( I, J ) )
            ENDDO
            RWORK( K ) = RWORK( K ) + S
         ENDDO
      END IF
      S = 0.0E+0
      DO I = 1, N
         IF( RWORK( I ) > SAFE2 ) THEN
            S = MAX( S, CABS1( WORK( I ) ) / RWORK( I ) )
         ELSE
            S = MAX( S, ( CABS1( WORK( I ) )+SAFE1 ) / &
                ( RWORK( I )+SAFE1 ) )
         END IF
      ENDDO
      BERR( J ) = S
!
!        Test stopping criterion. Continue iterating if
!           1) The residual BERR(J) is larger than machine epsilon, and
!           2) BERR(J) decreased by at least a factor of 2 during the
!              last iteration, and
!           3) At most ITMAX iterations tried.
!
      IF( BERR( J ) > EPS .AND. 2.0E+0*BERR( J ) <= LSTRES .AND. &
          COUNT <= ITMAX ) THEN
!
!           Update solution and try again.
!
         CALL CPOTRS( UPLO, N, 1, AF, LDAF, WORK, N, INFO )
         X(1:N,J) = X(1:N,J) + WORK(1:N)
         LSTRES = BERR( J )
         COUNT = COUNT + 1
         GO TO 20
      END IF
!
!        Bound error from formula
!
!        norm(X - XTRUE) / norm(X) .le. FERR =
!        norm( abs(inv(A))*
!           ( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) / norm(X)
!
!        where
!          norm(Z) is the magnitude of the largest component of Z
!          inv(A) is the inverse of A
!          abs(Z) is the componentwise absolute value of the matrix or
!             vector Z
!          NZ is the maximum number of nonzeros in any row of A, plus 1
!          EPS is machine epsilon
!
!        The i-th component of abs(R)+NZ*EPS*(abs(A)*abs(X)+abs(B))
!        is incremented by SAFE1 if the i-th component of
!        abs(A)*abs(X) + abs(B) is less than SAFE2.
!
!        Use CLACN2 to estimate the infinity-norm of the matrix
!           inv(A) * diag(W),
!        where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) )))
!
      DO I = 1, N
         IF( RWORK( I ) > SAFE2 ) THEN
            RWORK( I ) = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I )
         ELSE
            RWORK( I ) = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I ) + &
                         SAFE1
         END IF
      ENDDO
!
      KASE = 0
  100    CONTINUE
      CALL CLACN2( N, WORK( N+1 ), WORK, FERR( J ), KASE, ISAVE )
      IF( KASE /= 0 ) THEN
         IF( KASE == 1 ) THEN
!
!              Multiply by diag(W)*inv(A**H).
!
            CALL CPOTRS( UPLO, N, 1, AF, LDAF, WORK, N, INFO )
            WORK(1:N) = RWORK(1:N)*WORK(1:N)
         ELSE IF( KASE == 2 ) THEN
!
!              Multiply by inv(A)*diag(W).
!
            WORK(1:N) = RWORK(1:N)*WORK(1:N)
            CALL CPOTRS( UPLO, N, 1, AF, LDAF, WORK, N, INFO )
         END IF
         GO TO 100
      END IF
!
!        Normalize error.
!
      LSTRES = 0.0E+0
      DO I = 1, N
         LSTRES = MAX( LSTRES, CABS1( X( I, J ) ) )
      ENDDO
      IF( LSTRES /= 0.0E+0 ) FERR( J ) = FERR( J ) / LSTRES
!
      ENDDO
!
   RETURN
!
!     End of CPORFS
!
END
