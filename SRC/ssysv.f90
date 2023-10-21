!> \brief <b> SSYSV computes the solution to system of linear equations A * X = B for SY matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSYSV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssysv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssysv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssysv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
!                         LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSYSV computes the solution to a real system of linear equations
!>    A * X = B,
!> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
!> matrices.
!>
!> The diagonal pivoting method is used to factor A as
!>    A = U * D * U**T,  if UPLO = 'U', or
!>    A = L * D * L**T,  if UPLO = 'L',
!> where U (or L) is a product of permutation and unit upper (lower)
!> triangular matrices, and D is symmetric and block diagonal with
!> 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
!> used to solve the system of equations A * X = B.
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
!>          The number of linear equations, i.e., the order of the
!>          matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, if INFO = 0, the block diagonal matrix D and the
!>          multipliers used to obtain the factor U or L from the
!>          factorization A = U*D*U**T or A = L*D*L**T as computed by
!>          SSYTRF.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D, as
!>          determined by SSYTRF.  If IPIV(k) > 0, then rows and columns
!>          k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1
!>          diagonal block.  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0,
!>          then rows and columns k-1 and -IPIV(k) were interchanged and
!>          D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If UPLO = 'L' and
!>          IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and
!>          -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2
!>          diagonal block.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,NRHS)
!>          On entry, the N-by-NRHS right hand side matrix B.
!>          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of WORK.  LWORK >= 1, and for best performance
!>          LWORK >= max(1,N*NB), where NB is the optimal blocksize for
!>          SSYTRF.
!>          for LWORK < N, TRS will be done with Level BLAS 2
!>          for LWORK >= N, TRS will be done with Level BLAS 3
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization
!>               has been completed, but the block diagonal matrix D is
!>               exactly singular, so the solution could not be computed.
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
!> \ingroup hesv
!
!  =====================================================================
   SUBROUTINE SSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, &
                     LWORK, INFO )
!
!  -- LAPACK driver routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * )
   REAL               A( LDA, * ), B( LDB, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
   LOGICAL            LQUERY
   INTEGER            LWKOPT
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA, SSYTRF, SSYTRS, SSYTRS2
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   LQUERY = ( LWORK == -1 )
   IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( NRHS < 0 ) THEN
      INFO = -3
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -5
   ELSE IF( LDB < MAX( 1, N ) ) THEN
      INFO = -8
   ELSE IF( LWORK < 1 .AND. .NOT.LQUERY ) THEN
      INFO = -10
   END IF
!
   IF( INFO == 0 ) THEN
      IF( N == 0 ) THEN
         LWKOPT = 1
      ELSE
         CALL SSYTRF( UPLO, N, A, LDA, IPIV, WORK, -1, INFO )
         LWKOPT = INT( WORK( 1 ) )
      END IF
      WORK( 1 ) = LWKOPT
   END IF
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'SSYSV ', -INFO )
      RETURN
   ELSE IF( LQUERY ) THEN
      RETURN
   END IF
!
!     Compute the factorization A = U*D*U**T or A = L*D*L**T.
!
   CALL SSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
   IF( INFO == 0 ) THEN
!
!        Solve the system A*X = B, overwriting B with X.
!
      IF ( LWORK < N ) THEN
!
!        Solve with TRS ( Use Level BLAS 2)
!
         CALL SSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
      ELSE
!
!        Solve with TRS2 ( Use Level BLAS 3)
!
         CALL SSYTRS2( UPLO,N,NRHS,A,LDA,IPIV,B,LDB,WORK,INFO )
!
      END IF
!
   END IF
!
   WORK( 1 ) = LWKOPT
!
   RETURN
!
!     End of SSYSV
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

