!> \brief \b DGESC2 solves a system of linear equations using the LU factorization with complete pivoting computed by sgetc2.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGESC2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesc2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesc2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesc2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, N
!       DOUBLE PRECISION   SCALE
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), JPIV( * )
!       DOUBLE PRECISION   A( LDA, * ), RHS( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGESC2 solves a system of linear equations
!>
!>           A * X = scale* RHS
!>
!> with a general N-by-N matrix A using the LU factorization with
!> complete pivoting computed by DGETC2.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the  LU part of the factorization of the n-by-n
!>          matrix A computed by DGETC2:  A = P * L * U * Q
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1, N).
!> \endverbatim
!>
!> \param[in,out] RHS
!> \verbatim
!>          RHS is DOUBLE PRECISION array, dimension (N).
!>          On entry, the right hand side vector b.
!>          On exit, the solution vector X.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N).
!>          The pivot indices; for 1 <= i <= N, row i of the
!>          matrix has been interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[in] JPIV
!> \verbatim
!>          JPIV is INTEGER array, dimension (N).
!>          The pivot indices; for 1 <= j <= N, column j of the
!>          matrix has been interchanged with column JPIV(j).
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          On exit, SCALE contains the scale factor. SCALE is chosen
!>          0 <= SCALE <= 1 to prevent overflow in the solution.
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
!> \ingroup gesc2
!
!> \par Contributors:
!  ==================
!>
!>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
!>     Umea University, S-901 87 Umea, Sweden.
!
!  =====================================================================
   SUBROUTINE DGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, N
   DOUBLE PRECISION   SCALE
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * ), JPIV( * )
   DOUBLE PRECISION   A( LDA, * ), RHS( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, J
   DOUBLE PRECISION   BIGNUM, EPS, SMLNUM, TEMP
!     ..
!     .. External Subroutines ..
   EXTERNAL           DLASWP, DSCAL
!     ..
!     .. External Functions ..
   INTEGER            IDAMAX
   DOUBLE PRECISION   DLAMCH
   EXTERNAL           IDAMAX, DLAMCH
!     ..
!     .. Executable Statements ..
!
!      Set constant to control overflow
!
   EPS = DLAMCH( 'P' )
   SMLNUM = DLAMCH( 'S' ) / EPS
   BIGNUM = 1.0D0 / SMLNUM
!
!     Apply permutations IPIV to RHS
!
   CALL DLASWP( 1, RHS, LDA, 1, N-1, IPIV, 1 )
!
!     Solve for L part
!
   DO I = 1, N - 1
      RHS(I+1:N) = RHS(I+1:N) - A(I+1:N,I)*RHS(I)
   ENDDO
!
!     Solve for U part
!
   SCALE = 1.0D0
!
!     Check for scaling
!
   I = IDAMAX( N, RHS, 1 )
   IF( 2.0D0*SMLNUM*ABS( RHS( I ) ) > ABS( A( N, N ) ) ) THEN
      TEMP = ( 1.0D0 / 2.0D0 ) / ABS( RHS( I ) )
      RHS(1:N) = TEMP*RHS(1:N)
      SCALE = SCALE*TEMP
   END IF
!
   DO I = N, 1, -1
      TEMP = 1.0D0 / A( I, I )
      RHS( I ) = RHS( I )*TEMP
      RHS( I ) = RHS( I ) - SUM(RHS(I+1:N)*( A( I,I+1:N)*TEMP ))
   ENDDO
!
!     Apply permutations JPIV to the solution (RHS)
!
   CALL DLASWP( 1, RHS, LDA, 1, N-1, JPIV, -1 )
   RETURN
!
!     End of DGESC2
!
END
