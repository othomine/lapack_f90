!> \brief \b CGESC2 solves a system of linear equations using the LU factorization with complete pivoting computed by sgetc2.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGESC2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgesc2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgesc2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgesc2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, N
!       REAL               SCALE
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), JPIV( * )
!       COMPLEX            A( LDA, * ), RHS( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGESC2 solves a system of linear equations
!>
!>           A * X = scale* RHS
!>
!> with a general N-by-N matrix A using the LU factorization with
!> complete pivoting computed by CGETC2.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA, N)
!>          On entry, the  LU part of the factorization of the n-by-n
!>          matrix A computed by CGETC2:  A = P * L * U * Q
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
!>          RHS is COMPLEX array, dimension N.
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
!>          SCALE is REAL
!>           On exit, SCALE contains the scale factor. SCALE is chosen
!>           0 <= SCALE <= 1 to prevent overflow in the solution.
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
   SUBROUTINE CGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, N
   REAL               SCALE
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * ), JPIV( * )
   COMPLEX            A( LDA, * ), RHS( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, J
   REAL               BIGNUM, EPS, SMLNUM
   COMPLEX            TEMP
!     ..
!     .. External Subroutines ..
   EXTERNAL           CLASWP, CSCAL
!     ..
!     .. External Functions ..
   INTEGER            ICAMAX
   REAL               SLAMCH
   EXTERNAL           ICAMAX, SLAMCH
!     ..
!     .. Executable Statements ..
!
!     Set constant to control overflow
!
   EPS = SLAMCH( 'P' )
   SMLNUM = SLAMCH( 'S' ) / EPS
   BIGNUM = 1.0E+0 / SMLNUM
!
!     Apply permutations IPIV to RHS
!
   CALL CLASWP( 1, RHS, LDA, 1, N-1, IPIV, 1 )
!
!     Solve for L part
!
   DO I = 1, N - 1
      RHS(I+1:N) = RHS(I+1:N) - A(I+1:N,I)*RHS(I)
   ENDDO
!
!     Solve for U part
!
   SCALE = 1.0E+0
!
!     Check for scaling
!
   I = ICAMAX( N, RHS, 1 )
   IF( 2.0E+0*SMLNUM*ABS( RHS( I ) ) > ABS( A( N, N ) ) ) THEN
      TEMP = CMPLX(0.5E+0,0.0E+0)/ABS(RHS( I ) )
      RHS(1:N) = TEMP*RHS(1:N)
      SCALE = SCALE*REAL( TEMP )
   END IF
   DO I = N, 1, -1
      RHS( I ) = (RHS( I ) - SUM(RHS(I+1:N)*A(I,I+1:N)))/A(I,I)
   ENDDO
!
!     Apply permutations JPIV to the solution (RHS)
!
   CALL CLASWP( 1, RHS, LDA, 1, N-1, JPIV, -1 )
   RETURN
!
!     End of CGESC2
!
END
