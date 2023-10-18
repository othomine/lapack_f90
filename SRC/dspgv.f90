!> \brief \b DSPGV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSPGV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dspgv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dspgv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dspgv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK,
!                         INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, ITYPE, LDZ, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AP( * ), BP( * ), W( * ), WORK( * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSPGV computes all the eigenvalues and, optionally, the eigenvectors
!> of a real generalized symmetric-definite eigenproblem, of the form
!> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
!> Here A and B are assumed to be symmetric, stored in packed format,
!> and B is also positive definite.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          Specifies the problem type to be solved:
!>          = 1:  A*x = (lambda)*B*x
!>          = 2:  A*B*x = (lambda)*x
!>          = 3:  B*A*x = (lambda)*x
!> \endverbatim
!>
!> \param[in] JOBZ
!> \verbatim
!>          JOBZ is CHARACTER*1
!>          = 'N':  Compute eigenvalues only;
!>          = 'V':  Compute eigenvalues and eigenvectors.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangles of A and B are stored;
!>          = 'L':  Lower triangles of A and B are stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A and B.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the symmetric matrix
!>          A, packed columnwise in a linear array.  The j-th column of A
!>          is stored in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!>
!>          On exit, the contents of AP are destroyed.
!> \endverbatim
!>
!> \param[in,out] BP
!> \verbatim
!>          BP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the symmetric matrix
!>          B, packed columnwise in a linear array.  The j-th column of B
!>          is stored in the array BP as follows:
!>          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;
!>          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.
!>
!>          On exit, the triangular factor U or L from the Cholesky
!>          factorization B = U**T*U or B = L*L**T, in the same storage
!>          format as B.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (N)
!>          If INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
!>          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
!>          eigenvectors.  The eigenvectors are normalized as follows:
!>          if ITYPE = 1 or 2, Z**T*B*Z = I;
!>          if ITYPE = 3, Z**T*inv(B)*Z = I.
!>          If JOBZ = 'N', then Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1, and if
!>          JOBZ = 'V', LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (3*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  DPPTRF or DSPEV returned an error code:
!>             <= N:  if INFO = i, DSPEV failed to converge;
!>                    i off-diagonal elements of an intermediate
!>                    tridiagonal form did not converge to zero.
!>             > N:   if INFO = n + i, for 1 <= i <= n, then the leading
!>                    principal minor of order i of B is not positive.
!>                    The factorization of B could not be completed and
!>                    no eigenvalues or eigenvectors were computed.
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
!> \ingroup hpgv
!
!  =====================================================================
   SUBROUTINE DSPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, &
                     INFO )
!
!  -- LAPACK driver routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          JOBZ, UPLO
   INTEGER            INFO, ITYPE, LDZ, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   AP( * ), BP( * ), W( * ), WORK( * ), &
                      Z( LDZ, * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
   LOGICAL            UPPER, WANTZ
   CHARACTER          TRANS
   INTEGER            J, NEIG
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           DPPTRF, DSPEV, DSPGST, DTPMV, DTPSV, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   WANTZ = LSAME( JOBZ, 'V' )
   UPPER = LSAME( UPLO, 'U' )
!
   INFO = 0
   IF( ITYPE < 1 .OR. ITYPE > 3 ) THEN
      INFO = -1
   ELSE IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
      INFO = -2
   ELSE IF( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
      INFO = -3
   ELSE IF( N < 0 ) THEN
      INFO = -4
   ELSE IF( LDZ < 1 .OR. ( WANTZ .AND. LDZ < N ) ) THEN
      INFO = -9
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'DSPGV ', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) &
      RETURN
!
!     Form a Cholesky factorization of B.
!
   CALL DPPTRF( UPLO, N, BP, INFO )
   IF( INFO /= 0 ) THEN
      INFO = N + INFO
      RETURN
   END IF
!
!     Transform problem to standard eigenvalue problem and solve.
!
   CALL DSPGST( ITYPE, UPLO, N, AP, BP, INFO )
   CALL DSPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO )
!
   IF( WANTZ ) THEN
!
!        Backtransform eigenvectors to the original problem.
!
      NEIG = N
      IF( INFO > 0 ) &
         NEIG = INFO - 1
      IF( ITYPE == 1 .OR. ITYPE == 2 ) THEN
!
!           For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
!           backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y
!
         IF( UPPER ) THEN
            TRANS = 'N'
         ELSE
            TRANS = 'T'
         END IF
!
         DO J = 1, NEIG
            CALL DTPSV( UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), &
                        1 )
         ENDDO
!
      ELSE IF( ITYPE == 3 ) THEN
!
!           For B*A*x=(lambda)*x;
!           backtransform eigenvectors: x = L*y or U**T*y
!
         IF( UPPER ) THEN
            TRANS = 'T'
         ELSE
            TRANS = 'N'
         END IF
!
         DO J = 1, NEIG
            CALL DTPMV( UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), &
                        1 )
         ENDDO
      END IF
   END IF
   RETURN
!
!     End of DSPGV
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
