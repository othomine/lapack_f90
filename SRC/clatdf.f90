!> \brief \b CLATDF uses the LU factorization of the n-by-n matrix computed by sgetc2 and computes a contribution to the reciprocal Dif-estimate.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLATDF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clatdf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clatdf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clatdf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLATDF( IJOB, N, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV,
!                          JPIV )
!
!       .. Scalar Arguments ..
!       INTEGER            IJOB, LDZ, N
!       REAL               RDSCAL, RDSUM
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), JPIV( * )
!       COMPLEX            RHS( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLATDF computes the contribution to the reciprocal Dif-estimate
!> by solving for x in Z * x = b, where b is chosen such that the norm
!> of x is as large as possible. It is assumed that LU decomposition
!> of Z has been computed by CGETC2. On entry RHS = f holds the
!> contribution from earlier solved sub-systems, and on return RHS = x.
!>
!> The factorization of Z returned by CGETC2 has the form
!> Z = P * L * U * Q, where P and Q are permutation matrices. L is lower
!> triangular with unit diagonal elements and U is upper triangular.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] IJOB
!> \verbatim
!>          IJOB is INTEGER
!>          IJOB = 2: First compute an approximative null-vector e
!>              of Z using CGECON, e is normalized and solve for
!>              Zx = +-e - f with the sign giving the greater value of
!>              2-norm(x).  About 5 times as expensive as Default.
!>          IJOB .ne. 2: Local look ahead strategy where
!>              all entries of the r.h.s. b is chosen as either +1 or
!>              -1.  Default.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix Z.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ, N)
!>          On entry, the LU part of the factorization of the n-by-n
!>          matrix Z computed by CGETC2:  Z = P * L * U * Q
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDA >= max(1, N).
!> \endverbatim
!>
!> \param[in,out] RHS
!> \verbatim
!>          RHS is COMPLEX array, dimension (N).
!>          On entry, RHS contains contributions from other subsystems.
!>          On exit, RHS contains the solution of the subsystem with
!>          entries according to the value of IJOB (see above).
!> \endverbatim
!>
!> \param[in,out] RDSUM
!> \verbatim
!>          RDSUM is REAL
!>          On entry, the sum of squares of computed contributions to
!>          the Dif-estimate under computation by CTGSYL, where the
!>          scaling factor RDSCAL (see below) has been factored out.
!>          On exit, the corresponding sum of squares updated with the
!>          contributions from the current sub-system.
!>          If TRANS = 'T' RDSUM is not touched.
!>          NOTE: RDSUM only makes sense when CTGSY2 is called by CTGSYL.
!> \endverbatim
!>
!> \param[in,out] RDSCAL
!> \verbatim
!>          RDSCAL is REAL
!>          On entry, scaling factor used to prevent overflow in RDSUM.
!>          On exit, RDSCAL is updated w.r.t. the current contributions
!>          in RDSUM.
!>          If TRANS = 'T', RDSCAL is not touched.
!>          NOTE: RDSCAL only makes sense when CTGSY2 is called by
!>          CTGSYL.
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
!> \ingroup latdf
!
!> \par Further Details:
!  =====================
!>
!>  This routine is a further developed implementation of algorithm
!>  BSOLVE in [1] using complete pivoting in the LU factorization.
!
!> \par Contributors:
!  ==================
!>
!>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
!>     Umea University, S-901 87 Umea, Sweden.
!
!> \par References:
!  ================
!>
!>   [1]   Bo Kagstrom and Lars Westin,
!>         Generalized Schur Methods with Condition Estimators for
!>         Solving the Generalized Sylvester Equation, IEEE Transactions
!>         on Automatic Control, Vol. 34, No. 7, July 1989, pp 745-751.
!>
!>   [2]   Peter Poromaa,
!>         On Efficient and Robust Estimators for the Separation
!>         between two Regular Matrix Pairs with Applications in
!>         Condition Estimation. Report UMINF-95.05, Department of
!>         Computing Science, Umea University, S-901 87 Umea, Sweden,
!>         1995.
!
!  =====================================================================
   SUBROUTINE CLATDF( IJOB, N, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, &
                      JPIV )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            IJOB, LDZ, N
   REAL               RDSCAL, RDSUM
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * ), JPIV( * )
   COMPLEX            RHS( * ), Z( LDZ, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            MAXDIM
   PARAMETER          ( MAXDIM = 2 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, INFO, J, K
   REAL               RTEMP, SCALE, SMINU, SPLUS
   COMPLEX            BM, BP, PMONE, TEMP
!     ..
!     .. Local Arrays ..
   REAL               RWORK( MAXDIM )
   COMPLEX            WORK( 4*MAXDIM ), XM( MAXDIM ), XP( MAXDIM )
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGECON, CGESC2, CLASSQ, CLASWP
!     ..
!     .. External Functions ..
   REAL               SCASUM
   COMPLEX            CDOTC
   EXTERNAL           SCASUM, CDOTC
!     ..
!     .. Executable Statements ..
!
   IF( IJOB /= 2 ) THEN
!
!        Apply permutations IPIV to RHS
!
      CALL CLASWP( 1, RHS, LDZ, 1, N-1, IPIV, 1 )
!
!        Solve for L-part choosing RHS either to +1 or -1.
!
      PMONE = -(1.0E+0,0.0E+0)
      DO J = 1, N - 1
         BP = RHS( J ) + (1.0E+0,0.0E+0)
         BM = RHS( J ) - (1.0E+0,0.0E+0)
         SPLUS = 1.0E+0
!
!           Look-ahead for L- part RHS(1:N-1) = +-1
!           SPLUS and SMIN computed more efficiently than in BSOLVE[1].
!
         SPLUS = SPLUS + REAL( CDOTC( N-J, Z( J+1, J ), 1, Z( J+1, &
                 J ), 1 ) )
         SMINU = REAL( CDOTC( N-J, Z( J+1, J ), 1, RHS( J+1 ), 1 ) )
         SPLUS = SPLUS*REAL( RHS( J ) )
         IF( SPLUS > SMINU ) THEN
            RHS( J ) = BP
         ELSE IF( SMINU > SPLUS ) THEN
            RHS( J ) = BM
         ELSE
!
!              In this case the updating sums are equal and we can
!              choose RHS(J) +1 or -1. The first time this happens we
!              choose -1, thereafter +1. This is a simple way to get
!              good estimates of matrices like Byers well-known example
!              (see [1]). (Not done in BSOLVE.)
!
            RHS( J ) = RHS( J ) + PMONE
            PMONE = (1.0E+0,0.0E+0)
         END IF
!
!           Compute the remaining r.h.s.
!
         RHS(J+1:N) = RHS(J+1:N)-RHS(J)*Z(J+1:N,J)
      ENDDO
!
!        Solve for U- part, lockahead for RHS(N) = +-1. This is not done
!        In BSOLVE and will hopefully give us a better estimate because
!        any ill-conditioning of the original matrix is transferred to U
!        and not to L. U(N, N) is an approximation to sigma_min(LU).
!
      WORK(1:N-1) = RHS(1:N-1)
      WORK( N ) = RHS( N ) + (1.0E+0,0.0E+0)
      RHS( N ) = RHS( N ) - (1.0E+0,0.0E+0)
      SPLUS = 0.0E+0
      SMINU = 0.0E+0
      DO I = N, 1, -1
         TEMP = (1.0E+0,0.0E+0) / Z( I, I )
         WORK(I) = TEMP*(WORK(I) - SUM(WORK(I+1:N)*Z(I,I+1:N)))
         RHS(I) = TEMP*(RHS(I) - SUM(RHS(I+1:N)*Z(I,I+1:N)))
         SPLUS = SPLUS + ABS( WORK( I ) )
         SMINU = SMINU + ABS( RHS( I ) )
      ENDDO
      IF( SPLUS > SMINU ) RHS(1:N) = WORK(1:N)
!
!        Apply the permutations JPIV to the computed solution (RHS)
!
      CALL CLASWP( 1, RHS, LDZ, 1, N-1, JPIV, -1 )
!
!        Compute the sum of squares
!
      CALL CLASSQ( N, RHS, 1, RDSCAL, RDSUM )
      RETURN
   END IF
!
!     ENTRY IJOB = 2
!
!     Compute approximate nullvector XM of Z
!
   CALL CGECON( 'I', N, Z, LDZ, 1.0E+0, RTEMP, WORK, RWORK, INFO )
   XM(1:N) = WORK(N+1:2*N)
!
!     Compute RHS
!
   CALL CLASWP( 1, XM, LDZ, 1, N-1, IPIV, -1 )
   TEMP = (1.0E+0,0.0E+0) / SQRT( CDOTC( N, XM, 1, XM, 1 ) )
   XM(1:N) = TEMP*XM(1:N)
   XP(1:N) = XM(1:N)
   XP(1:N) = XP(1:N) + RHS(1:N)
   RHS(1:N) = RHS(1:N) - XM(1:N)
   CALL CGESC2( N, Z, LDZ, RHS, IPIV, JPIV, SCALE )
   CALL CGESC2( N, Z, LDZ, XP, IPIV, JPIV, SCALE )
   IF( SCASUM( N, XP, 1 ) > SCASUM( N, RHS, 1 ) ) RHS(1:N) = XP(1:N)
!
!     Compute the sum of squares
!
   CALL CLASSQ( N, RHS, 1, RDSCAL, RDSUM )
   RETURN
!
!     End of CLATDF
!
END
