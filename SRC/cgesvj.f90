!> \brief <b> CGESVJ </b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGESVJ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgesvj.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgesvj.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgesvj.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGESVJ( JOBA, JOBU, JOBV, M, N, A, LDA, SVA, MV, V,
!                          LDV, CWORK, LWORK, RWORK, LRWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDV, LWORK, LRWORK, M, MV, N
!       CHARACTER*1        JOBA, JOBU, JOBV
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ),  V( LDV, * ), CWORK( LWORK )
!       REAL               RWORK( LRWORK ),  SVA( N )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGESVJ computes the singular value decomposition (SVD) of a complex
!> M-by-N matrix A, where M >= N. The SVD of A is written as
!>                                    [++]   [xx]   [x0]   [xx]
!>              A = U * SIGMA * V^*,  [++] = [xx] * [ox] * [xx]
!>                                    [++]   [xx]
!> where SIGMA is an N-by-N diagonal matrix, U is an M-by-N orthonormal
!> matrix, and V is an N-by-N unitary matrix. The diagonal elements
!> of SIGMA are the singular values of A. The columns of U and V are the
!> left and the right singular vectors of A, respectively.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBA
!> \verbatim
!>          JOBA is CHARACTER*1
!>          Specifies the structure of A.
!>          = 'L': The input matrix A is lower triangular;
!>          = 'U': The input matrix A is upper triangular;
!>          = 'G': The input matrix A is general M-by-N matrix, M >= N.
!> \endverbatim
!>
!> \param[in] JOBU
!> \verbatim
!>          JOBU is CHARACTER*1
!>          Specifies whether to compute the left singular vectors
!>          (columns of U):
!>          = 'U' or 'F': The left singular vectors corresponding to the nonzero
!>                 singular values are computed and returned in the leading
!>                 columns of A. See more details in the description of A.
!>                 The default numerical orthogonality threshold is set to
!>                 approximately TOL=CTOL*EPS, CTOL=SQRT(M), EPS=SLAMCH('E').
!>          = 'C': Analogous to JOBU='U', except that user can control the
!>                 level of numerical orthogonality of the computed left
!>                 singular vectors. TOL can be set to TOL = CTOL*EPS, where
!>                 CTOL is given on input in the array WORK.
!>                 No CTOL smaller than 1.0E+0 is allowed. CTOL greater
!>                 than 1 / EPS is meaningless. The option 'C'
!>                 can be used if M*EPS is satisfactory orthogonality
!>                 of the computed left singular vectors, so CTOL=M could
!>                 save few sweeps of Jacobi rotations.
!>                 See the descriptions of A and WORK(1).
!>          = 'N': The matrix U is not computed. However, see the
!>                 description of A.
!> \endverbatim
!>
!> \param[in] JOBV
!> \verbatim
!>          JOBV is CHARACTER*1
!>          Specifies whether to compute the right singular vectors, that
!>          is, the matrix V:
!>          = 'V' or 'J': the matrix V is computed and returned in the array V
!>          = 'A':  the Jacobi rotations are applied to the MV-by-N
!>                  array V. In other words, the right singular vector
!>                  matrix V is not computed explicitly; instead it is
!>                  applied to an MV-by-N matrix initially stored in the
!>                  first MV rows of V.
!>          = 'N':  the matrix V is not computed and the array V is not
!>                  referenced
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the input matrix A. 1/SLAMCH('E') > M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the input matrix A.
!>          M >= N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit,
!>          If JOBU = 'U' .OR. JOBU = 'C':
!>                 If INFO = 0 :
!>                 RANKA orthonormal columns of U are returned in the
!>                 leading RANKA columns of the array A. Here RANKA <= N
!>                 is the number of computed singular values of A that are
!>                 above the underflow threshold SLAMCH('S'). The singular
!>                 vectors corresponding to underflowed or zero singular
!>                 values are not computed. The value of RANKA is returned
!>                 in the array RWORK as RANKA=NINT(RWORK(2)). Also see the
!>                 descriptions of SVA and RWORK. The computed columns of U
!>                 are mutually numerically orthogonal up to approximately
!>                 TOL=SQRT(M)*EPS (default); or TOL=CTOL*EPS (JOBU = 'C'),
!>                 see the description of JOBU.
!>                 If INFO > 0,
!>                 the procedure CGESVJ did not converge in the given number
!>                 of iterations (sweeps). In that case, the computed
!>                 columns of U may not be orthogonal up to TOL. The output
!>                 U (stored in A), SIGMA (given by the computed singular
!>                 values in SVA(1:N)) and V is still a decomposition of the
!>                 input matrix A in the sense that the residual
!>                 || A - SCALE * U * SIGMA * V^* ||_2 / ||A||_2 is small.
!>          If JOBU = 'N':
!>                 If INFO = 0 :
!>                 Note that the left singular vectors are 'for free' in the
!>                 one-sided Jacobi SVD algorithm. However, if only the
!>                 singular values are needed, the level of numerical
!>                 orthogonality of U is not an issue and iterations are
!>                 stopped when the columns of the iterated matrix are
!>                 numerically orthogonal up to approximately M*EPS. Thus,
!>                 on exit, A contains the columns of U scaled with the
!>                 corresponding singular values.
!>                 If INFO > 0 :
!>                 the procedure CGESVJ did not converge in the given number
!>                 of iterations (sweeps).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] SVA
!> \verbatim
!>          SVA is REAL array, dimension (N)
!>          On exit,
!>          If INFO = 0 :
!>          depending on the value SCALE = RWORK(1), we have:
!>                 If SCALE = 1.0E+0:
!>                 SVA(1:N) contains the computed singular values of A.
!>                 During the computation SVA contains the Euclidean column
!>                 norms of the iterated matrices in the array A.
!>                 If SCALE  /=  1.0E+0:
!>                 The singular values of A are SCALE*SVA(1:N), and this
!>                 factored representation is due to the fact that some of the
!>                 singular values of A might underflow or overflow.
!>
!>          If INFO > 0 :
!>          the procedure CGESVJ did not converge in the given number of
!>          iterations (sweeps) and SCALE*SVA(1:N) may not be accurate.
!> \endverbatim
!>
!> \param[in] MV
!> \verbatim
!>          MV is INTEGER
!>          If JOBV = 'A', then the product of Jacobi rotations in CGESVJ
!>          is applied to the first MV rows of V. See the description of JOBV.
!> \endverbatim
!>
!> \param[in,out] V
!> \verbatim
!>          V is COMPLEX array, dimension (LDV,N)
!>          If JOBV = 'V', then V contains on exit the N-by-N matrix of
!>                         the right singular vectors;
!>          If JOBV = 'A', then V contains the product of the computed right
!>                         singular vector matrix and the initial matrix in
!>                         the array V.
!>          If JOBV = 'N', then V is not referenced.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V, LDV >= 1.
!>          If JOBV = 'V', then LDV >= max(1,N).
!>          If JOBV = 'A', then LDV >= max(1,MV) .
!> \endverbatim
!>
!> \param[in,out] CWORK
!> \verbatim
!>          CWORK is COMPLEX array, dimension (max(1,LWORK))
!>          Used as workspace.
!>          If on entry LWORK = -1, then a workspace query is assumed and
!>          no computation is done; CWORK(1) is set to the minial (and optimal)
!>          length of CWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER.
!>          Length of CWORK, LWORK >= M+N.
!> \endverbatim
!>
!> \param[in,out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (max(6,LRWORK))
!>          On entry,
!>          If JOBU = 'C' :
!>          RWORK(1) = CTOL, where CTOL defines the threshold for convergence.
!>                    The process stops if all columns of A are mutually
!>                    orthogonal up to CTOL*EPS, EPS=SLAMCH('E').
!>                    It is required that CTOL >= 1.0E+0, i.e. it is not
!>                    allowed to force the routine to obtain orthogonality
!>                    below EPSILON.
!>          On exit,
!>          RWORK(1) = SCALE is the scaling factor such that SCALE*SVA(1:N)
!>                    are the computed singular values of A.
!>                    (See description of SVA().)
!>          RWORK(2) = NINT(RWORK(2)) is the number of the computed nonzero
!>                    singular values.
!>          RWORK(3) = NINT(RWORK(3)) is the number of the computed singular
!>                    values that are larger than the underflow threshold.
!>          RWORK(4) = NINT(RWORK(4)) is the number of sweeps of Jacobi
!>                    rotations needed for numerical convergence.
!>          RWORK(5) = max_{i /= j} |COS(A(:,i),A(:,j))| in the last sweep.
!>                    This is useful information in cases when CGESVJ did
!>                    not converge, as it can be used to estimate whether
!>                    the output is still useful and for post festum analysis.
!>          RWORK(6) = the largest absolute value over all sines of the
!>                    Jacobi rotation angles in the last sweep. It can be
!>                    useful for a post festum analysis.
!>         If on entry LRWORK = -1, then a workspace query is assumed and
!>         no computation is done; RWORK(1) is set to the minial (and optimal)
!>         length of RWORK.
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>         LRWORK is INTEGER
!>         Length of RWORK, LRWORK >= MAX(6,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, then the i-th argument had an illegal value
!>          > 0:  CGESVJ did not converge in the maximal allowed number
!>                (NSWEEP=30) of sweeps. The output may still be useful.
!>                See the description of RWORK.
!> \endverbatim
!>
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
!
!> \ingroup gesvj
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!> The orthogonal N-by-N matrix V is obtained as a product of Jacobi plane
!> rotations. In the case of underflow of the tangent of the Jacobi angle, a
!> modified Jacobi transformation of Drmac [3] is used. Pivot strategy uses
!> column interchanges of de Rijk [1]. The relative accuracy of the computed
!> singular values and the accuracy of the computed singular vectors (in
!> angle metric) is as guaranteed by the theory of Demmel and Veselic [2].
!> The condition number that determines the accuracy in the full rank case
!> is essentially min_{D=diag} kappa(A*D), where kappa(.) is the
!> spectral condition number. The best performance of this Jacobi SVD
!> procedure is achieved if used in an  accelerated version of Drmac and
!> Veselic [4,5], and it is the kernel routine in the SIGMA library [6].
!> Some tuning parameters (marked with [TP]) are available for the
!> implementer.
!> The computational range for the nonzero singular values is the  machine
!> number interval ( UNDERFLOW , OVERFLOW ). In extreme cases, even
!> denormalized singular values can be computed with the corresponding
!> gradual loss of accurate digits.
!> \endverbatim
!
!> \par Contributor:
!  ==================
!>
!> \verbatim
!>
!>  ============
!>
!>  Zlatko Drmac (Zagreb, Croatia)
!>
!> \endverbatim
!
!> \par References:
!  ================
!>
!> \verbatim
!>
!> [1] P. P. M. De Rijk: A one-sided Jacobi algorithm for computing the
!>    singular value decomposition on a vector computer.
!>    SIAM J. Sci. Stat. Comp., Vol. 10 (1998), pp. 359-371.
!> [2] J. Demmel and K. Veselic: Jacobi method is more accurate than QR.
!> [3] Z. Drmac: Implementation of Jacobi rotations for accurate singular
!>    value computation in floating point arithmetic.
!>    SIAM J. Sci. Comp., Vol. 18 (1997), pp. 1200-1222.
!> [4] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm I.
!>    SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1322-1342.
!>    LAPACK Working note 169.
!> [5] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm II.
!>    SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1343-1362.
!>    LAPACK Working note 170.
!> [6] Z. Drmac: SIGMA - mathematical software library for accurate SVD, PSV,
!>    QSVD, (H,K)-SVD computations.
!>    Department of Mathematics, University of Zagreb, 2008, 2015.
!> \endverbatim
!
!> \par Bugs, examples and comments:
!  =================================
!>
!> \verbatim
!>  ===========================
!>  Please report all bugs and send interesting test examples and comments to
!>  drmac@math.hr. Thank you.
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE CGESVJ( JOBA, JOBU, JOBV, M, N, A, LDA, SVA, MV, V, &
                      LDV, CWORK, LWORK, RWORK, LRWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
   IMPLICIT NONE
!     .. Scalar Arguments ..
   INTEGER            INFO, LDA, LDV, LWORK, LRWORK, M, MV, N
   CHARACTER*1        JOBA, JOBU, JOBV
!     ..
!     .. Array Arguments ..
   COMPLEX            A( LDA, * ),  V( LDV, * ), CWORK( LWORK )
   REAL               RWORK( LRWORK ), SVA( N )
!     ..
!
!  =====================================================================
!
!     .. Local Parameters ..
   INTEGER      NSWEEP
   PARAMETER  ( NSWEEP = 30 )
!     ..
!     .. Local Scalars ..
   COMPLEX A_TMP( LDA ), V_TMP( LDV )
   COMPLEX AAPQ, OMPQ
   REAL    AAPP, AAPP0, AAPQ1, AAQQ, APOAQ, AQOAP, BIG, &
           BIGTHETA, CS, CTOL, EPSLN, MXAAPQ, &
           MXSINJ, ROOTBIG, ROOTEPS, ROOTSFMIN, ROOTTOL, &
           SKL, SFMIN, SMALL, SN, T, TEMP1, THETA, THSIGN, TOL
   INTEGER BLSKIP, EMPTSW, i, ibr, IERR, igl, IJBLSK, ir1, &
           ISWROT, jbc, jgl, KBL, LKAHEAD, MVL, N2, N34, &
           N4, NBL, NOTROT, p, PSKIPPED, q, ROWSKIP, SWBAND
   LOGICAL APPLV, GOSCALE, LOWER, LQUERY, LSVEC, NOSCALE, ROTOK, &
           RSVEC, UCTOL, UPPER
!     ..
!     ..
!     .. External Functions ..
!     ..
!     from BLAS
   REAL               SCNRM2
   COMPLEX            CDOTC
   EXTERNAL           CDOTC, SCNRM2
   INTEGER            ISAMAX
   EXTERNAL           ISAMAX
!     from LAPACK
   REAL               SLAMCH
   EXTERNAL           SLAMCH
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!     ..
!     from BLAS
   EXTERNAL           CCOPY, CROT, CSSCAL
!     from LAPACK
   EXTERNAL           CLASCL, CLASET, CLASSQ, SLASCL, XERBLA
   EXTERNAL           CGSVJ0, CGSVJ1
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
   LSVEC = LSAME( JOBU, 'U' ) .OR. LSAME( JOBU, 'F' )
   UCTOL = LSAME( JOBU, 'C' )
   RSVEC = LSAME( JOBV, 'V' ) .OR. LSAME( JOBV, 'J' )
   APPLV = LSAME( JOBV, 'A' )
   UPPER = LSAME( JOBA, 'U' )
   LOWER = LSAME( JOBA, 'L' )
!
   LQUERY = ( LWORK  ==  -1 ) .OR. ( LRWORK  ==  -1 )
   IF( .NOT.( UPPER .OR. LOWER .OR. LSAME( JOBA, 'G' ) ) ) THEN
      INFO = -1
   ELSE IF( .NOT.( LSVEC .OR. UCTOL .OR. LSAME( JOBU, 'N' ) ) ) THEN
      INFO = -2
   ELSE IF( .NOT.( RSVEC .OR. APPLV .OR. LSAME( JOBV, 'N' ) ) ) THEN
      INFO = -3
   ELSE IF( M < 0 ) THEN
      INFO = -4
   ELSE IF( ( N < 0 ) .OR. ( N > M ) ) THEN
      INFO = -5
   ELSE IF( LDA < M ) THEN
      INFO = -7
   ELSE IF( MV < 0 ) THEN
      INFO = -9
   ELSE IF( ( RSVEC .AND. ( LDV < N ) ) .OR. &
             ( APPLV .AND. ( LDV < MV ) ) ) THEN
      INFO = -11
   ELSE IF( UCTOL .AND. ( RWORK( 1 ) <= 1.0E+0 ) ) THEN
      INFO = -12
   ELSE IF( LWORK < ( M+N ) .AND. ( .NOT.LQUERY ) ) THEN
      INFO = -13
   ELSE IF( LRWORK < MAX( N, 6 ) .AND. ( .NOT.LQUERY ) ) THEN
      INFO = -15
   ELSE
      INFO = 0
   END IF
!
!     #:(
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CGESVJ', -INFO )
      RETURN
   ELSE IF ( LQUERY ) THEN
      CWORK(1) = M + N
      RWORK(1) = MAX( N, 6 )
      RETURN
   END IF
!
! #:) Quick return for void matrix
!
   IF( ( M == 0 ) .OR. ( N == 0 ) )RETURN
!
!     Set numerical parameters
!     The stopping criterion for Jacobi rotations is
!
!     max_{i<>j}|A(:,i)^* * A(:,j)| / (||A(:,i)||*||A(:,j)||) < CTOL*EPS
!
!     where EPS is the round-off and CTOL is defined as follows:
!
   IF( UCTOL ) THEN
!        ... user controlled
      CTOL = RWORK( 1 )
   ELSE
!        ... default
      IF( LSVEC .OR. RSVEC .OR. APPLV ) THEN
         CTOL = SQRT( REAL( M ) )
      ELSE
         CTOL = REAL( M )
      END IF
   END IF
!     ... and the machine dependent parameters are
![!]  (Make sure that SLAMCH() works properly on the target machine.)
!
   EPSLN = SLAMCH( 'Epsilon' )
   ROOTEPS = SQRT( EPSLN )
   SFMIN = SLAMCH( 'SafeMinimum' )
   ROOTSFMIN = SQRT( SFMIN )
   SMALL = SFMIN / EPSLN
!      BIG = SLAMCH( 'Overflow' )
   BIG     = 1.0E+0  / SFMIN
   ROOTBIG = 1.0E+0 / ROOTSFMIN
!     LARGE = BIG / SQRT( REAL( M*N ) )
   BIGTHETA = 1.0E+0 / ROOTEPS
!
   TOL = CTOL*EPSLN
   ROOTTOL = SQRT( TOL )
!
   IF( REAL( M )*EPSLN >= 1.0E+0 ) THEN
      INFO = -4
      CALL XERBLA( 'CGESVJ', -INFO )
      RETURN
   END IF
!
!     Initialize the right singular vector matrix.
!
   IF( RSVEC ) THEN
      MVL = N
      CALL CLASET( 'A', MVL, N, (0.0E+0,0.0E+0), (1.0E+0,0.0E+0), V, LDV )
   ELSE IF( APPLV ) THEN
      MVL = MV
   END IF
   RSVEC = RSVEC .OR. APPLV
!
!     Initialize SVA( 1:N ) = ( ||A e_i||_2, i = 1:N )
!(!)  If necessary, scale A to protect the largest singular value
!     from overflow. It is possible that saving the largest singular
!     value destroys the information about the small ones.
!     This initial scaling is almost minimal in the sense that the
!     goal is to make sure that no column norm overflows, and that
!     SQRT(N)*max_i SVA(i) does not overflow. If INFinite entries
!     in A are detected, the procedure returns with INFO=-6.
!
   SKL = 1.0E+0 / SQRT( REAL( M )*REAL( N ) )
   NOSCALE = .TRUE.
   GOSCALE = .TRUE.
!
   IF( LOWER ) THEN
!        the input matrix is M-by-N lower triangular (trapezoidal)
      DO p = 1, N
         AAPP = 0.0E+0
         AAQQ = 1.0E+0
         CALL CLASSQ( M-p+1, A( p, p ), 1, AAPP, AAQQ )
         IF( AAPP > BIG ) THEN
            INFO = -6
            CALL XERBLA( 'CGESVJ', -INFO )
            RETURN
         END IF
         AAQQ = SQRT( AAQQ )
         IF( ( AAPP < ( BIG / AAQQ ) ) .AND. NOSCALE ) THEN
            SVA( p ) = AAPP*AAQQ
         ELSE
            NOSCALE = .FALSE.
            SVA( p ) = AAPP*( AAQQ*SKL )
            IF( GOSCALE ) THEN
               GOSCALE = .FALSE.
               SVA(1:p-1) = SVA(1:p-1)*SKL
            END IF
         END IF
      ENDDO
   ELSE IF( UPPER ) THEN
!        the input matrix is M-by-N upper triangular (trapezoidal)
      DO p = 1, N
         AAPP = 0.0E+0
         AAQQ = 1.0E+0
         CALL CLASSQ( p, A( 1, p ), 1, AAPP, AAQQ )
         IF( AAPP > BIG ) THEN
            INFO = -6
            CALL XERBLA( 'CGESVJ', -INFO )
            RETURN
         END IF
         AAQQ = SQRT( AAQQ )
         IF( ( AAPP < ( BIG / AAQQ ) ) .AND. NOSCALE ) THEN
            SVA( p ) = AAPP*AAQQ
         ELSE
            NOSCALE = .FALSE.
            SVA( p ) = AAPP*( AAQQ*SKL )
            IF( GOSCALE ) THEN
               GOSCALE = .FALSE.
               SVA(1:p-1) = SVA(1:p-1)*SKL
            END IF
         END IF
      ENDDO
   ELSE
!        the input matrix is M-by-N general dense
      DO p = 1, N
         AAPP = 0.0E+0
         AAQQ = 1.0E+0
         CALL CLASSQ( M, A( 1, p ), 1, AAPP, AAQQ )
         IF( AAPP > BIG ) THEN
            INFO = -6
            CALL XERBLA( 'CGESVJ', -INFO )
            RETURN
         END IF
         AAQQ = SQRT( AAQQ )
         IF( ( AAPP < ( BIG / AAQQ ) ) .AND. NOSCALE ) THEN
            SVA( p ) = AAPP*AAQQ
         ELSE
            NOSCALE = .FALSE.
            SVA( p ) = AAPP*( AAQQ*SKL )
            IF( GOSCALE ) THEN
               GOSCALE = .FALSE.
               SVA(1:p-1) = SVA(1:p-1)*SKL
            END IF
         END IF
         ENDDO
   END IF
!
   IF( NOSCALE )SKL = 1.0E+0
!
!     Move the smaller part of the spectrum from the underflow threshold
!(!)  Start by determining the position of the nonzero entries of the
!     array SVA() relative to ( SFMIN, BIG ).
!
   AAPP = 0.0E+0
   AAQQ = BIG
   DO p = 1, N
      IF( SVA( p ) /= 0.0E+0 ) AAQQ = MIN( AAQQ, SVA( p ) )
      AAPP = MAX( AAPP, SVA( p ) )
   ENDDO
!
! #:) Quick return for zero matrix
!
   IF( AAPP == 0.0E+0 ) THEN
      IF( LSVEC )CALL CLASET( 'G', M, N, (0.0E+0,0.0E+0), (1.0E+0,0.0E+0), A, LDA )
      RWORK( 1 ) = 1.0E+0
      RWORK( 2:6 ) = 0.0E+0
      RETURN
   END IF
!
! #:) Quick return for one-column matrix
!
   IF( N == 1 ) THEN
      IF( LSVEC )CALL CLASCL( 'G', 0, 0, SVA( 1 ), SKL, M, 1, &
                              A( 1, 1 ), LDA, IERR )
      RWORK( 1 ) = 1.0E+0 / SKL
      IF( SVA( 1 ) >= SFMIN ) THEN
         RWORK( 2 ) = 1.0E+0
      ELSE
         RWORK( 2 ) = 0.0E+0
      END IF
      RWORK( 3:6 ) = 0.0E+0
      RETURN
   END IF
!
!     Protect small singular values from underflow, and try to
!     avoid underflows/overflows in computing Jacobi rotations.
!
   SN = SQRT( SFMIN / EPSLN )
   TEMP1 = SQRT( BIG / REAL( N ) )
   IF( ( AAPP <= SN ) .OR. ( AAQQ >= TEMP1 ) .OR. &
       ( ( SN <= AAQQ ) .AND. ( AAPP <= TEMP1 ) ) ) THEN
      TEMP1 = MIN( BIG, TEMP1 / AAPP )
!         AAQQ  = AAQQ*TEMP1
!         AAPP  = AAPP*TEMP1
   ELSE IF( ( AAQQ <= SN ) .AND. ( AAPP <= TEMP1 ) ) THEN
      TEMP1 = MIN( SN / AAQQ, BIG / ( AAPP*SQRT( REAL( N ) ) ) )
!         AAQQ  = AAQQ*TEMP1
!         AAPP  = AAPP*TEMP1
   ELSE IF( ( AAQQ >= SN ) .AND. ( AAPP >= TEMP1 ) ) THEN
      TEMP1 = MAX( SN / AAQQ, TEMP1 / AAPP )
!         AAQQ  = AAQQ*TEMP1
!         AAPP  = AAPP*TEMP1
   ELSE IF( ( AAQQ <= SN ) .AND. ( AAPP >= TEMP1 ) ) THEN
      TEMP1 = MIN( SN / AAQQ, BIG / ( SQRT( REAL( N ) )*AAPP ) )
!         AAQQ  = AAQQ*TEMP1
!         AAPP  = AAPP*TEMP1
   ELSE
      TEMP1 = 1.0E+0
   END IF
!
!     Scale, if necessary
!
   IF( TEMP1 /= 1.0E+0 ) THEN
      CALL SLASCL( 'G', 0, 0, 1.0E+0, TEMP1, N, 1, SVA, N, IERR )
   END IF
   SKL = TEMP1*SKL
   IF( SKL /= 1.0E+0 ) THEN
      CALL CLASCL( JOBA, 0, 0, 1.0E+0, SKL, M, N, A, LDA, IERR )
      SKL = 1.0E+0 / SKL
   END IF
!
!     Row-cyclic Jacobi SVD algorithm with column pivoting
!
   EMPTSW = ( N*( N-1 ) ) / 2
   NOTROT = 0

   CWORK(1:N) = (1.0E+0,0.0E+0)
!
!
!
   SWBAND = 3
![TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective
!     if CGESVJ is used as a computational routine in the preconditioned
!     Jacobi SVD algorithm CGEJSV. For sweeps i=1:SWBAND the procedure
!     works on pivots inside a band-like region around the diagonal.
!     The boundaries are determined dynamically, based on the number of
!     pivots above a threshold.
!
   KBL = MIN( 8, N )
![TP] KBL is a tuning parameter that defines the tile size in the
!     tiling of the p-q loops of pivot pairs. In general, an optimal
!     value of KBL depends on the matrix dimensions and on the
!     parameters of the computer's memory.
!
   NBL = N / KBL
   IF( ( NBL*KBL ) /= N )NBL = NBL + 1
!
   BLSKIP = KBL**2
![TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.
!
   ROWSKIP = MIN( 5, KBL )
![TP] ROWSKIP is a tuning parameter.
!
   LKAHEAD = 1
![TP] LKAHEAD is a tuning parameter.
!
!     Quasi block transformations, using the lower (upper) triangular
!     structure of the input matrix. The quasi-block-cycling usually
!     invokes cubic convergence. Big part of this cycle is done inside
!     canonical subspaces of dimensions less than M.
!
   IF( ( LOWER .OR. UPPER ) .AND. ( N > MAX( 64, 4*KBL ) ) ) THEN
![TP] The number of partition levels and the actual partition are
!     tuning parameters.
      N4 = N / 4
      N2 = N / 2
      N34 = 3*N4
      IF( APPLV ) THEN
         q = 0
      ELSE
         q = 1
      END IF
!
      IF( LOWER ) THEN
!
!     This works very well on lower triangular matrices, in particular
!     in the framework of the preconditioned Jacobi SVD (xGEJSV).
!     The idea is simple:
!     [+ 0 0 0]   Note that Jacobi transformations of [0 0]
!     [+ + 0 0]                                       [0 0]
!     [+ + x 0]   actually work on [x 0]              [x 0]
!     [+ + x x]                    [x x].             [x x]
!
         CALL CGSVJ0( JOBV, M-N34, N-N34, A( N34+1, N34+1 ), LDA, &
                      CWORK( N34+1 ), SVA( N34+1 ), MVL, &
                      V( N34*q+1, N34+1 ), LDV, EPSLN, SFMIN, TOL, &
                      2, CWORK( N+1 ), LWORK-N, IERR )

         CALL CGSVJ0( JOBV, M-N2, N34-N2, A( N2+1, N2+1 ), LDA, &
                      CWORK( N2+1 ), SVA( N2+1 ), MVL, &
                      V( N2*q+1, N2+1 ), LDV, EPSLN, SFMIN, TOL, 2, &
                      CWORK( N+1 ), LWORK-N, IERR )

         CALL CGSVJ1( JOBV, M-N2, N-N2, N4, A( N2+1, N2+1 ), LDA, &
                      CWORK( N2+1 ), SVA( N2+1 ), MVL, &
                      V( N2*q+1, N2+1 ), LDV, EPSLN, SFMIN, TOL, 1, &
                      CWORK( N+1 ), LWORK-N, IERR )
!
         CALL CGSVJ0( JOBV, M-N4, N2-N4, A( N4+1, N4+1 ), LDA, &
                      CWORK( N4+1 ), SVA( N4+1 ), MVL, &
                      V( N4*q+1, N4+1 ), LDV, EPSLN, SFMIN, TOL, 1, &
                      CWORK( N+1 ), LWORK-N, IERR )
!
         CALL CGSVJ0( JOBV, M, N4, A, LDA, CWORK, SVA, MVL, V, LDV, &
                      EPSLN, SFMIN, TOL, 1, CWORK( N+1 ), LWORK-N, &
                      IERR )
!
         CALL CGSVJ1( JOBV, M, N2, N4, A, LDA, CWORK, SVA, MVL, V, &
                      LDV, EPSLN, SFMIN, TOL, 1, CWORK( N+1 ), &
                      LWORK-N, IERR )
!
!
      ELSE IF( UPPER ) THEN
!
!
         CALL CGSVJ0( JOBV, N4, N4, A, LDA, CWORK, SVA, MVL, V, LDV, &
                      EPSLN, SFMIN, TOL, 2, CWORK( N+1 ), LWORK-N, &
                      IERR )
!
         CALL CGSVJ0( JOBV, N2, N4, A( 1, N4+1 ), LDA, CWORK( N4+1 ), &
                      SVA( N4+1 ), MVL, V( N4*q+1, N4+1 ), LDV, &
                      EPSLN, SFMIN, TOL, 1, CWORK( N+1 ), LWORK-N, &
                      IERR )
!
         CALL CGSVJ1( JOBV, N2, N2, N4, A, LDA, CWORK, SVA, MVL, V, &
                      LDV, EPSLN, SFMIN, TOL, 1, CWORK( N+1 ), &
                      LWORK-N, IERR )
!
         CALL CGSVJ0( JOBV, N2+N4, N4, A( 1, N2+1 ), LDA, &
                      CWORK( N2+1 ), SVA( N2+1 ), MVL, &
                      V( N2*q+1, N2+1 ), LDV, EPSLN, SFMIN, TOL, 1, &
                      CWORK( N+1 ), LWORK-N, IERR )

      END IF
!
   END IF
!
!     .. Row-cyclic pivot strategy with de Rijk's pivoting ..
!
   DO i = 1, NSWEEP
!
!     .. go go go ...
!
      MXAAPQ = 0.0E+0
      MXSINJ = 0.0E+0
      ISWROT = 0
!
      NOTROT = 0
      PSKIPPED = 0
!
!     Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs
!     1 <= p < q <= N. This is the first step toward a blocked implementation
!     of the rotations. New implementation, based on block transformations,
!     is under development.
!
      DO ibr = 1, NBL
!
         igl = ( ibr-1 )*KBL + 1
!
         DO ir1 = 0, MIN( LKAHEAD, NBL-ibr )
!
            igl = igl + ir1*KBL
!
            DO p = igl, MIN( igl+KBL-1, N-1 )
!
!     .. de Rijk's pivoting
!
               q = ISAMAX( N-p+1, SVA( p ), 1 ) + p - 1
               IF( p /= q ) THEN
                  A_TMP(1:M) = A(1:M,p)
                  A(1:M,p) = A(1:M,q)
                  A(1:M,q) = A_TMP(1:M)
                  IF( RSVEC ) THEN
                     V_TMP(1:MVL) = V(1:MVL,p)
                     V(1:MVL,p) = V(1:MVL,q)
                     V(1:MVL,q) = V_TMP(1:MVL)
                  ENDIF
                  TEMP1 = SVA( p )
                  SVA( p ) = SVA( q )
                  SVA( q ) = TEMP1
                  AAPQ = CWORK(p)
                  CWORK(p) = CWORK(q)
                  CWORK(q) = AAPQ
               END IF
!
               IF( ir1 == 0 ) THEN
!
!        Column norms are periodically updated by explicit
!        norm computation.
![!]     Caveat:
!        Unfortunately, some BLAS implementations compute SCNRM2(M,A(1,p),1)
!        as SQRT(S=CDOTC(M,A(1,p),1,A(1,p),1)), which may cause the result to
!        overflow for ||A(:,p)||_2 > SQRT(overflow_threshold), and to
!        underflow for ||A(:,p)||_2 < SQRT(underflow_threshold).
!        Hence, SCNRM2 cannot be trusted, not even in the case when
!        the true norm is far from the under(over)flow boundaries.
!        If properly implemented SCNRM2 is available, the IF-THEN-ELSE-END IF
!        below should be replaced with "AAPP = SCNRM2( M, A(1,p), 1 )".
!
                  IF( ( SVA( p ) < ROOTBIG ) .AND. &
                       ( SVA( p ) > ROOTSFMIN ) ) THEN
                     SVA( p ) = SCNRM2( M, A( 1, p ), 1 )
                  ELSE
                     TEMP1 = 0.0E+0
                     AAPP = 1.0E+0
                     CALL CLASSQ( M, A( 1, p ), 1, TEMP1, AAPP )
                     SVA( p ) = TEMP1*SQRT( AAPP )
                  END IF
                  AAPP = SVA( p )
               ELSE
                  AAPP = SVA( p )
               END IF
!
               IF( AAPP > 0.0E+0 ) THEN
!
                  PSKIPPED = 0
!
                  DO q = p + 1, MIN( igl+KBL-1, N )
!
                     AAQQ = SVA( q )
!
                     IF( AAQQ > 0.0E+0 ) THEN
!
                        AAPP0 = AAPP
                        IF( AAQQ >= 1.0E+0 ) THEN
                           ROTOK = ( SMALL*AAPP ) <= AAQQ
                           IF( AAPP < ( BIG / AAQQ ) ) THEN
                              AAPQ = ( CDOTC( M, A( 1, p ), 1, &
                                      A( 1, q ), 1 ) / AAQQ ) / AAPP
                           ELSE
                              CWORK(N+1:N+M) = A(1:M,p)
                              CALL CLASCL( 'G', 0, 0, AAPP, 1.0E+0, &
                                   M, 1, CWORK(N+1), LDA, IERR )
                              AAPQ = CDOTC( M, CWORK(N+1), 1, &
                                      A( 1, q ), 1 ) / AAQQ
                           END IF
                        ELSE
                           ROTOK = AAPP <= ( AAQQ / SMALL )
                           IF( AAPP > ( SMALL / AAQQ ) ) THEN
                              AAPQ = ( CDOTC( M, A( 1, p ), 1, &
                                       A( 1, q ), 1 ) / AAPP ) / AAQQ
                           ELSE
                              CWORK(N+1:N+M) = A(1:M,q)
                              CALL CLASCL( 'G', 0, 0, AAQQ, &
                                            1.0E+0, M, 1, &
                                            CWORK(N+1), LDA, IERR )
                              AAPQ = CDOTC( M, A(1, p ), 1, &
                                      CWORK(N+1), 1 ) / AAPP
                           END IF
                        END IF
!
!                           AAPQ = AAPQ * CONJG( CWORK(p) ) * CWORK(q)
                        AAPQ1  = -ABS(AAPQ)
                        MXAAPQ = MAX( MXAAPQ, -AAPQ1 )
!
!        TO rotate or NOT to rotate, THAT is the question ...
!
                        IF( ABS( AAPQ1 ) > TOL ) THEN
                            OMPQ = AAPQ / ABS(AAPQ)
!
!           .. rotate
![RTD]      ROTATED = ROTATED + 1.0E+0
!
                           IF( ir1 == 0 ) THEN
                              NOTROT = 0
                              PSKIPPED = 0
                              ISWROT = ISWROT + 1
                           END IF
!
                           IF( ROTOK ) THEN
!
                              AQOAP = AAQQ / AAPP
                              APOAQ = AAPP / AAQQ
                              THETA = -0.5E+0*ABS( AQOAP-APOAQ )/AAPQ1
!
                              IF( ABS( THETA ) > BIGTHETA ) THEN
!
                                 T  = 0.5E+0 / THETA
                                 CS = 1.0E+0

                                 CALL CROT( M, A(1,p), 1, A(1,q), 1, &
                                             CS, CONJG(OMPQ)*T )
                                 IF ( RSVEC ) THEN
                                     CALL CROT( MVL, V(1,p), 1, &
                                     V(1,q), 1, CS, CONJG(OMPQ)*T )
                                 END IF

                                 SVA( q ) = AAQQ*SQRT( MAX( 0.0E+0, &
                                             1.0E+0+T*APOAQ*AAPQ1 ) )
                                 AAPP = AAPP*SQRT( MAX( 0.0E+0, &
                                             1.0E+0-T*AQOAP*AAPQ1 ) )
                                 MXSINJ = MAX( MXSINJ, ABS( T ) )
!
                              ELSE
!
!                 .. choose correct signum for THETA and rotate
!
                                 THSIGN = -SIGN( 1.0E+0, AAPQ1 )
                                 T = 1.0E+0 / ( THETA+THSIGN* &
                                      SQRT( 1.0E+0+THETA*THETA ) )
                                 CS = SQRT( 1.0E+0 / ( 1.0E+0+T*T ) )
                                 SN = T*CS
!
                                 MXSINJ = MAX( MXSINJ, ABS( SN ) )
                                 SVA( q ) = AAQQ*SQRT( MAX( 0.0E+0, &
                                             1.0E+0+T*APOAQ*AAPQ1 ) )
                                 AAPP = AAPP*SQRT( MAX( 0.0E+0, &
                                         1.0E+0-T*AQOAP*AAPQ1 ) )
!
                                 CALL CROT( M, A(1,p), 1, A(1,q), 1, &
                                             CS, CONJG(OMPQ)*SN )
                                 IF ( RSVEC ) THEN
                                     CALL CROT( MVL, V(1,p), 1, &
                                     V(1,q), 1, CS, CONJG(OMPQ)*SN )
                                 END IF
                              END IF
                              CWORK(p) = -CWORK(q) * OMPQ
!
                              ELSE
!              .. have to use modified Gram-Schmidt like transformation
                              CWORK(N+1:N+M) = A(1:M,p)
                              CALL CLASCL( 'G', 0, 0, AAPP, 1.0E+0, M, &
                                           1, CWORK(N+1), LDA, &
                                           IERR )
                              CALL CLASCL( 'G', 0, 0, AAQQ, 1.0E+0, M, &
                                           1, A( 1, q ), LDA, IERR )
                              A(1:M,q) = A(1:M,q) - AAPQ*CWORK(N+1:N+M)
                              CALL CLASCL( 'G', 0, 0, 1.0E+0, AAQQ, M, &
                                           1, A( 1, q ), LDA, IERR )
                              SVA( q ) = AAQQ*SQRT( MAX( 0.0E+0, &
                                         1.0E+0-AAPQ1*AAPQ1 ) )
                              MXSINJ = MAX( MXSINJ, SFMIN )
                           END IF
!           END IF ROTOK THEN ... ELSE
!
!           In the case of cancellation in updating SVA(q), SVA(p)
!           recompute SVA(q), SVA(p).
!
                           IF( ( SVA( q ) / AAQQ )**2 <= ROOTEPS ) THEN
                              IF( ( AAQQ < ROOTBIG ) .AND. &
                                  ( AAQQ > ROOTSFMIN ) ) THEN
                                 SVA( q ) = SCNRM2( M, A( 1, q ), 1 )
                              ELSE
                                 T = 0.0E+0
                                 AAQQ = 1.0E+0
                                 CALL CLASSQ( M, A( 1, q ), 1, T, AAQQ )
                                 SVA( q ) = T*SQRT( AAQQ )
                              END IF
                           END IF
                           IF( ( AAPP / AAPP0 ) <= ROOTEPS ) THEN
                              IF( ( AAPP < ROOTBIG ) .AND. &
                                  ( AAPP > ROOTSFMIN ) ) THEN
                                 AAPP = SCNRM2( M, A( 1, p ), 1 )
                              ELSE
                                 T = 0.0E+0
                                 AAPP = 1.0E+0
                                 CALL CLASSQ( M, A( 1, p ), 1, T, AAPP )
                                 AAPP = T*SQRT( AAPP )
                              END IF
                              SVA( p ) = AAPP
                           END IF
!
                        ELSE
!                             A(:,p) and A(:,q) already numerically orthogonal
                           IF( ir1 == 0 )NOTROT = NOTROT + 1
![RTD]      SKIPPED  = SKIPPED + 1
                           PSKIPPED = PSKIPPED + 1
                        END IF
                     ELSE
!                          A(:,q) is zero column
                        IF( ir1 == 0 )NOTROT = NOTROT + 1
                        PSKIPPED = PSKIPPED + 1
                     END IF
!
                     IF( ( i <= SWBAND ) .AND. &
                         ( PSKIPPED > ROWSKIP ) ) THEN
                        IF( ir1 == 0 )AAPP = -AAPP
                        NOTROT = 0
                        GO TO 2103
                     END IF
!
                     ENDDO
!     END q-LOOP
!
 2103                CONTINUE
!     bailed out of q-loop
!
                  SVA( p ) = AAPP
!
               ELSE
                  SVA( p ) = AAPP
                  IF( ( ir1 == 0 ) .AND. ( AAPP == 0.0E+0 ) ) &
                      NOTROT = NOTROT + MIN( igl+KBL-1, N ) - p
               END IF
!
               ENDDO
!     end of the p-loop
!     end of doing the block ( ibr, ibr )
            ENDDO
!     end of ir1-loop
!
! ... go to the off diagonal blocks
!
         igl = ( ibr-1 )*KBL + 1
!
         DO jbc = ibr + 1, NBL
!
            jgl = ( jbc-1 )*KBL + 1
!
!        doing the block at ( ibr, jbc )
!
            IJBLSK = 0
            DO p = igl, MIN( igl+KBL-1, N )
!
               AAPP = SVA( p )
               IF( AAPP > 0.0E+0 ) THEN
!
                  PSKIPPED = 0
!
                  DO q = jgl, MIN( jgl+KBL-1, N )
!
                     AAQQ = SVA( q )
                     IF( AAQQ > 0.0E+0 ) THEN
                        AAPP0 = AAPP
!
!     .. M x 2 Jacobi SVD ..
!
!        Safe Gram matrix computation
!
                        IF( AAQQ >= 1.0E+0 ) THEN
                           IF( AAPP >= AAQQ ) THEN
                              ROTOK = ( SMALL*AAPP ) <= AAQQ
                           ELSE
                              ROTOK = ( SMALL*AAQQ ) <= AAPP
                           END IF
                           IF( AAPP < ( BIG / AAQQ ) ) THEN
                              AAPQ = ( CDOTC( M, A( 1, p ), 1, &
                                     A( 1, q ), 1 ) / AAQQ ) / AAPP
                           ELSE
                              CWORK(N+1:N+M) = A(1:M,p)
                              CALL CLASCL( 'G', 0, 0, AAPP, &
                                           1.0E+0, M, 1, &
                                           CWORK(N+1), LDA, IERR )
                              AAPQ = CDOTC( M, CWORK(N+1), 1, &
                                     A( 1, q ), 1 ) / AAQQ
                           END IF
                        ELSE
                           IF( AAPP >= AAQQ ) THEN
                              ROTOK = AAPP <= ( AAQQ / SMALL )
                           ELSE
                              ROTOK = AAQQ <= ( AAPP / SMALL )
                           END IF
                           IF( AAPP > ( SMALL / AAQQ ) ) THEN
                              AAPQ = ( CDOTC( M, A( 1, p ), 1, &
                                    A( 1, q ), 1 ) / MAX(AAQQ,AAPP) ) &
                                                   / MIN(AAQQ,AAPP)
                           ELSE
                              CWORK(N+1:N+M) = A(1:M,q)
                              CALL CLASCL( 'G', 0, 0, AAQQ, &
                                           1.0E+0, M, 1, &
                                           CWORK(N+1), LDA, IERR )
                              AAPQ = CDOTC( M, A( 1, p ), 1, &
                                     CWORK(N+1),  1 ) / AAPP
                           END IF
                        END IF
!
!                           AAPQ = AAPQ * CONJG(CWORK(p))*CWORK(q)
                        AAPQ1  = -ABS(AAPQ)
                        MXAAPQ = MAX( MXAAPQ, -AAPQ1 )
!
!        TO rotate or NOT to rotate, THAT is the question ...
!
                        IF( ABS( AAPQ1 ) > TOL ) THEN
                           OMPQ = AAPQ / ABS(AAPQ)
                           NOTROT = 0
![RTD]      ROTATED  = ROTATED + 1
                           PSKIPPED = 0
                           ISWROT = ISWROT + 1
!
                           IF( ROTOK ) THEN
!
                              AQOAP = AAQQ / AAPP
                              APOAQ = AAPP / AAQQ
                              THETA = -0.5E+0*ABS( AQOAP-APOAQ )/ AAPQ1
                              IF( AAQQ > AAPP0 )THETA = -THETA
!
                              IF( ABS( THETA ) > BIGTHETA ) THEN
                                 T  = 0.5E+0 / THETA
                                 CS = 1.0E+0
                                 CALL CROT( M, A(1,p), 1, A(1,q), 1, &
                                             CS, CONJG(OMPQ)*T )
                                 IF( RSVEC ) THEN
                                     CALL CROT( MVL, V(1,p), 1, &
                                     V(1,q), 1, CS, CONJG(OMPQ)*T )
                                 END IF
                                 SVA( q ) = AAQQ*SQRT( MAX( 0.0E+0, &
                                            1.0E+0+T*APOAQ*AAPQ1 ) )
                                 AAPP = AAPP*SQRT( MAX( 0.0E+0, &
                                        1.0E+0-T*AQOAP*AAPQ1 ) )
                                 MXSINJ = MAX( MXSINJ, ABS( T ) )
                              ELSE
!
!                 .. choose correct signum for THETA and rotate
!
                                 THSIGN = -SIGN( 1.0E+0, AAPQ1 )
                                 IF( AAQQ > AAPP0 )THSIGN = -THSIGN
                                 T = 1.0E+0 / ( THETA+THSIGN* &
                                     SQRT( 1.0E+0+THETA*THETA ) )
                                 CS = SQRT( 1.0E+0 / ( 1.0E+0+T*T ) )
                                 SN = T*CS
                                 MXSINJ = MAX( MXSINJ, ABS( SN ) )
                                 SVA( q ) = AAQQ*SQRT( MAX( 0.0E+0, &
                                            1.0E+0+T*APOAQ*AAPQ1 ) )
                                 AAPP = AAPP*SQRT( MAX( 0.0E+0, &
                                            1.0E+0-T*AQOAP*AAPQ1 ) )
!
                                 CALL CROT( M, A(1,p), 1, A(1,q), 1, &
                                             CS, CONJG(OMPQ)*SN )
                                 IF( RSVEC ) THEN
                                     CALL CROT( MVL, V(1,p), 1, &
                                     V(1,q), 1, CS, CONJG(OMPQ)*SN )
                                 END IF
                              END IF
                              CWORK(p) = -CWORK(q) * OMPQ
!
                           ELSE
!              .. have to use modified Gram-Schmidt like transformation
                            IF( AAPP > AAQQ ) THEN
                                 CWORK(N+1:N+M) = A(1:M,p)
                                 CALL CLASCL( 'G', 0, 0, AAPP, 1.0E+0, &
                                              M, 1, CWORK(N+1),LDA, &
                                              IERR )
                                 CALL CLASCL( 'G', 0, 0, AAQQ, 1.0E+0, &
                                              M, 1, A( 1, q ), LDA, &
                                              IERR )
                                 A(1:M,q) = A(1:M,q) - AAPQ*CWORK(N+1:N+M)
                                 CALL CLASCL( 'G', 0, 0, 1.0E+0, AAQQ, &
                                              M, 1, A( 1, q ), LDA, &
                                              IERR )
                                 SVA( q ) = AAQQ*SQRT( MAX( 0.0E+0, &
                                            1.0E+0-AAPQ1*AAPQ1 ) )
                                 MXSINJ = MAX( MXSINJ, SFMIN )
                            ELSE
                                 CWORK(N+1:N+M) = A(1:M,q)
                                 CALL CLASCL( 'G', 0, 0, AAQQ, 1.0E+0, &
                                              M, 1, CWORK(N+1),LDA, &
                                              IERR )
                                 CALL CLASCL( 'G', 0, 0, AAPP, 1.0E+0, &
                                              M, 1, A( 1, p ), LDA, &
                                              IERR )
                                 A(1:M,p) = A(1:M,p)-CONJG(AAPQ)*CWORK(N+1:N+M)
                                 CALL CLASCL( 'G', 0, 0, 1.0E+0, AAPP, &
                                              M, 1, A( 1, p ), LDA, &
                                              IERR )
                                 SVA( p ) = AAPP*SQRT( MAX( 0.0E+0, &
                                            1.0E+0-AAPQ1*AAPQ1 ) )
                                 MXSINJ = MAX( MXSINJ, SFMIN )
                            END IF
                           END IF
!           END IF ROTOK THEN ... ELSE
!
!           In the case of cancellation in updating SVA(q), SVA(p)
!           .. recompute SVA(q), SVA(p)
                           IF( ( SVA( q ) / AAQQ )**2 <= ROOTEPS ) &
                               THEN
                              IF( ( AAQQ < ROOTBIG ) .AND. &
                                  ( AAQQ > ROOTSFMIN ) ) THEN
                                 SVA( q ) = SCNRM2( M, A( 1, q ), 1)
                               ELSE
                                 T = 0.0E+0
                                 AAQQ = 1.0E+0
                                 CALL CLASSQ( M, A( 1, q ), 1, T, &
                                              AAQQ )
                                 SVA( q ) = T*SQRT( AAQQ )
                              END IF
                           END IF
                           IF( ( AAPP / AAPP0 )**2 <= ROOTEPS ) THEN
                              IF( ( AAPP < ROOTBIG ) .AND. &
                                  ( AAPP > ROOTSFMIN ) ) THEN
                                 AAPP = SCNRM2( M, A( 1, p ), 1 )
                              ELSE
                                 T = 0.0E+0
                                 AAPP = 1.0E+0
                                 CALL CLASSQ( M, A( 1, p ), 1, T, &
                                              AAPP )
                                 AAPP = T*SQRT( AAPP )
                              END IF
                              SVA( p ) = AAPP
                           END IF
!              end of OK rotation
                        ELSE
                           NOTROT = NOTROT + 1
![RTD]      SKIPPED  = SKIPPED  + 1
                           PSKIPPED = PSKIPPED + 1
                           IJBLSK = IJBLSK + 1
                        END IF
                     ELSE
                        NOTROT = NOTROT + 1
                        PSKIPPED = PSKIPPED + 1
                        IJBLSK = IJBLSK + 1
                     END IF
!
                     IF( ( i <= SWBAND ) .AND. ( IJBLSK >= BLSKIP ) ) &
                         THEN
                        SVA( p ) = AAPP
                        NOTROT = 0
                        GO TO 2011
                     END IF
                     IF( ( i <= SWBAND ) .AND. &
                         ( PSKIPPED > ROWSKIP ) ) THEN
                        AAPP = -AAPP
                        NOTROT = 0
                        EXIT
                     END IF
!
                  ENDDO
!        end of the q-loop
!
                  SVA( p ) = AAPP
!
               ELSE
!
                  IF( AAPP == 0.0E+0 )NOTROT = NOTROT + &
                      MIN( jgl+KBL-1, N ) - jgl + 1
                  IF( AAPP < 0.0E+0 )NOTROT = 0
!
               END IF
!
            ENDDO
!     end of the p-loop
         ENDDO
!     end of the jbc-loop
 2011       CONTINUE
!2011 bailed out of the jbc-loop
         DO p = igl, MIN( igl+KBL-1, N )
            SVA( p ) = ABS( SVA( p ) )
         ENDDO
!**
      ENDDO
!2000 :: end of the ibr-loop
!
!     .. update SVA(N)
      IF( ( SVA( N ) < ROOTBIG ) .AND. ( SVA( N ) > ROOTSFMIN ) ) THEN
         SVA( N ) = SCNRM2( M, A( 1, N ), 1 )
      ELSE
         T = 0.0E+0
         AAPP = 1.0E+0
         CALL CLASSQ( M, A( 1, N ), 1, T, AAPP )
         SVA( N ) = T*SQRT( AAPP )
      END IF
!
!     Additional steering devices
!
      IF( ( i < SWBAND ) .AND. ( ( MXAAPQ <= ROOTTOL ) .OR. &
          ( ISWROT <= N ) ) )SWBAND = i
!
      IF( ( i > SWBAND+1 ) .AND. ( MXAAPQ < SQRT( REAL( N ) )* &
          TOL ) .AND. ( REAL( N )*MXAAPQ*MXSINJ < TOL ) ) THEN
         GO TO 1994
      END IF
!
      IF( NOTROT >= EMPTSW )GO TO 1994
!
   ENDDO
!     end i=1:NSWEEP loop
!
! #:( Reaching this point means that the procedure has not converged.
   INFO = NSWEEP - 1
   GO TO 1995
!
 1994 CONTINUE
! #:) Reaching this point means numerical convergence after the i-th
!     sweep.
!
   INFO = 0
! #:) INFO = 0 confirms successful iterations.
 1995 CONTINUE
!
!     Sort the singular values and find how many are above
!     the underflow threshold.
!
   N2 = 0
   N4 = 0
   DO p = 1, N - 1
      q = ISAMAX( N-p+1, SVA( p ), 1 ) + p - 1
      IF( p /= q ) THEN
         TEMP1 = SVA( p )
         SVA( p ) = SVA( q )
         SVA( q ) = TEMP1
         A_TMP(1:M) = A(1:M,p)
         A(1:M,p) = A(1:M,q)
         A(1:M,q) = A_TMP(1:M)
         IF( RSVEC ) THEN
            V_TMP(1:MVL) = V(1:MVL,p)
            V(1:MVL,p) = V(1:MVL,q)
            V(1:MVL,q) = V_TMP(1:MVL)
         ENDIF
      END IF
      IF( SVA( p ) /= 0.0E+0 ) THEN
         N4 = N4 + 1
         IF( SVA( p )*SKL > SFMIN )N2 = N2 + 1
      END IF
   ENDDO
   IF( SVA( N ) /= 0.0E+0 ) THEN
      N4 = N4 + 1
      IF( SVA( N )*SKL > SFMIN )N2 = N2 + 1
   END IF
!
!     Normalize the left singular vectors.
!
   IF( LSVEC .OR. UCTOL ) THEN
      DO p = 1, N4
!           CALL CSSCAL( M, 1.0E+0 / SVA( p ), A( 1, p ), 1 )
         CALL CLASCL( 'G',0,0, SVA(p), 1.0E+0, M, 1, A(1,p), M, IERR )
      ENDDO
   END IF
!
!     Scale the product of Jacobi rotations.
!
   IF( RSVEC ) THEN
         DO p = 1, N
            V(1:MVL,p) = V(1:MVL,p) / SCNRM2( MVL, V( 1, p ), 1 )
         ENDDO
   END IF
!
!     Undo scaling, if necessary (and possible).
   IF( ( ( SKL > 1.0E+0 ) .AND. ( SVA( 1 ) < ( BIG / SKL ) ) ) &
       .OR. ( ( SKL < 1.0E+0 ) .AND. ( SVA( MAX( N2, 1 ) )  > &
       ( SFMIN / SKL ) ) ) ) THEN
      SVA(1:N) = SKL*SVA(1:N)
      SKL = 1.0E+0
   END IF
!
   RWORK( 1 ) = SKL
!     The singular values of A are SKL*SVA(1:N). If SKL /= 1.0E+0
!     then some of the singular values may overflow or underflow and
!     the spectrum is given in this factored representation.
!
   RWORK( 2 ) = REAL( N4 )
!     N4 is the number of computed nonzero singular values of A.
!
   RWORK( 3 ) = REAL( N2 )
!     N2 is the number of singular values of A greater than SFMIN.
!     If N2<N, SVA(N2:N) contains 0.0E+0S and/or denormalized numbers
!     that may carry some information.
!
   RWORK( 4 ) = REAL( i )
!     i is the index of the last sweep before declaring convergence.
!
   RWORK( 5 ) = MXAAPQ
!     MXAAPQ is the largest absolute value of scaled pivots in the
!     last sweep
!
   RWORK( 6 ) = MXSINJ
!     MXSINJ is the largest absolute value of the sines of Jacobi angles
!     in the last sweep
!
   RETURN
!     ..
!     .. END OF CGESVJ
!     ..
   END
!
