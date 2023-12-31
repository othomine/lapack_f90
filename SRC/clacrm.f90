!> \brief \b CLACRM multiplies a complex matrix by a square real matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLACRM + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacrm.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacrm.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacrm.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLACRM( M, N, A, LDA, B, LDB, C, LDC, RWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LDC, M, N
!       ..
!       .. Array Arguments ..
!       REAL               B( LDB, * ), RWORK( * )
!       COMPLEX            A( LDA, * ), C( LDC, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLACRM performs a very simple matrix-matrix multiplication:
!>          C := A * B,
!> where A is M by N and complex; B is N by N and real;
!> C is M by N and complex.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A and of the matrix C.
!>          M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns and rows of the matrix B and
!>          the number of columns of the matrix C.
!>          N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA, N)
!>          On entry, A contains the M by N matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >=max(1,M).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (LDB, N)
!>          On entry, B contains the N by N matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >=max(1,N).
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC, N)
!>          On exit, C contains the M by N matrix C.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >=max(1,N).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (2*M*N)
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
!> \ingroup lacrm
!
!  =====================================================================
   SUBROUTINE CLACRM( M, N, A, LDA, B, LDB, C, LDC, RWORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, LDB, LDC, M, N
!     ..
!     .. Array Arguments ..
   REAL               B( LDB, * ), RWORK( * )
   COMPLEX            A( LDA, * ), C( LDC, * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, J, L
!     ..
!     .. External Subroutines ..
   EXTERNAL           SGEMM
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible.
!
   IF( ( M == 0 ) .OR. ( N == 0 ) ) RETURN
!
   DO J = 1, N
      RWORK((J-1)*M+1:J*M) = REAL(A(1:M,J))
   ENDDO
!
   L = M*N + 1
   CALL SGEMM( 'N', 'N', M, N, N, 1.0E+0, RWORK, M, B, LDB, 0.0E+0, &
               RWORK( L ), M )
   DO J = 1, N
      C(1:M,J) = RWORK(L+(J-1)*M:L+J*M-1)
   ENDDO
!
   DO J = 1, N
      RWORK((J-1)*M+1:J*M) = AIMAG(A(1:M,J))
   ENDDO
   CALL SGEMM( 'N', 'N', M, N, N, 1.0E+0, RWORK, M, B, LDB, 0.0E+0, &
               RWORK( L ), M )
   DO J = 1, N
      C(1:M,J) = CMPLX(REAL(C(1:M,J)),RWORK(L+(J-1)*M:L+J*M-1))
   ENDDO
!
   RETURN
!
!     End of CLACRM
!
END
