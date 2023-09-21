!> \brief <b> CHBEVX_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
!
!  @generated from zhbevx_2stage.f, fortran z -> c, Sat Nov  5 23:18:22 2016
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHBEVX_2STAGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbevx_2stage.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbevx_2stage.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbevx_2stage.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHBEVX_2STAGE( JOBZ, RANGE, UPLO, N, KD, AB, LDAB,
!                                 Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W,
!                                 Z, LDZ, WORK, LWORK, RWORK, IWORK,
!                                 IFAIL, INFO )
!
!       IMPLICIT NONE
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, RANGE, UPLO
!       INTEGER            IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N, LWORK
!       REAL               ABSTOL, VL, VU
!       ..
!       .. Array Arguments ..
!       INTEGER            IFAIL( * ), IWORK( * )
!       REAL               RWORK( * ), W( * )
!       COMPLEX            AB( LDAB, * ), Q( LDQ, * ), WORK( * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHBEVX_2STAGE computes selected eigenvalues and, optionally, eigenvectors
!> of a complex Hermitian band matrix A using the 2stage technique for
!> the reduction to tridiagonal.  Eigenvalues and eigenvectors
!> can be selected by specifying either a range of values or a range of
!> indices for the desired eigenvalues.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBZ
!> \verbatim
!>          JOBZ is CHARACTER*1
!>          = 'N':  Compute eigenvalues only;
!>          = 'V':  Compute eigenvalues and eigenvectors.
!>                  Not available in this release.
!> \endverbatim
!>
!> \param[in] RANGE
!> \verbatim
!>          RANGE is CHARACTER*1
!>          = 'A': all eigenvalues will be found;
!>          = 'V': all eigenvalues in the half-open interval (VL,VU]
!>                 will be found;
!>          = 'I': the IL-th through IU-th eigenvalues will be found.
!> \endverbatim
!>
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
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of superdiagonals of the matrix A if UPLO = 'U',
!>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is COMPLEX array, dimension (LDAB, N)
!>          On entry, the upper or lower triangle of the Hermitian band
!>          matrix A, stored in the first KD+1 rows of the array.  The
!>          j-th column of A is stored in the j-th column of the array AB
!>          as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>
!>          On exit, AB is overwritten by values generated during the
!>          reduction to tridiagonal form.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD + 1.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX array, dimension (LDQ, N)
!>          If JOBZ = 'V', the N-by-N unitary matrix used in the
!>                          reduction to tridiagonal form.
!>          If JOBZ = 'N', the array Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  If JOBZ = 'V', then
!>          LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[in] VL
!> \verbatim
!>          VL is REAL
!>          If RANGE='V', the lower bound of the interval to
!>          be searched for eigenvalues. VL < VU.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] VU
!> \verbatim
!>          VU is REAL
!>          If RANGE='V', the upper bound of the interval to
!>          be searched for eigenvalues. VL < VU.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] IL
!> \verbatim
!>          IL is INTEGER
!>          If RANGE='I', the index of the
!>          smallest eigenvalue to be returned.
!>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!>          Not referenced if RANGE = 'A' or 'V'.
!> \endverbatim
!>
!> \param[in] IU
!> \verbatim
!>          IU is INTEGER
!>          If RANGE='I', the index of the
!>          largest eigenvalue to be returned.
!>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!>          Not referenced if RANGE = 'A' or 'V'.
!> \endverbatim
!>
!> \param[in] ABSTOL
!> \verbatim
!>          ABSTOL is REAL
!>          The absolute error tolerance for the eigenvalues.
!>          An approximate eigenvalue is accepted as converged
!>          when it is determined to lie in an interval [a,b]
!>          of width less than or equal to
!>
!>                  ABSTOL + EPS *   max( |a|,|b| ) ,
!>
!>          where EPS is the machine precision.  If ABSTOL is less than
!>          or equal to zero, then  EPS*|T|  will be used in its place,
!>          where |T| is the 1-norm of the tridiagonal matrix obtained
!>          by reducing AB to tridiagonal form.
!>
!>          Eigenvalues will be computed most accurately when ABSTOL is
!>          set to twice the underflow threshold 2*SLAMCH('S'), not zero.
!>          If this routine returns with INFO>0, indicating that some
!>          eigenvectors did not converge, try setting ABSTOL to
!>          2*SLAMCH('S').
!>
!>          See "Computing Small Singular Values of Bidiagonal Matrices
!>          with Guaranteed High Relative Accuracy," by Demmel and
!>          Kahan, LAPACK Working Note #3.
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The total number of eigenvalues found.  0 <= M <= N.
!>          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is REAL array, dimension (N)
!>          The first M elements contain the selected eigenvalues in
!>          ascending order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ, max(1,M))
!>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
!>          contain the orthonormal eigenvectors of the matrix A
!>          corresponding to the selected eigenvalues, with the i-th
!>          column of Z holding the eigenvector associated with W(i).
!>          If an eigenvector fails to converge, then that column of Z
!>          contains the latest approximation to the eigenvector, and the
!>          index of the eigenvector is returned in IFAIL.
!>          If JOBZ = 'N', then Z is not referenced.
!>          Note: the user must ensure that at least max(1,M) columns are
!>          supplied in the array Z; if RANGE = 'V', the exact value of M
!>          is not known in advance and an upper bound must be used.
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
!>          WORK is COMPLEX array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK. LWORK >= 1, when N <= 1;
!>          otherwise
!>          If JOBZ = 'N' and N > 1, LWORK must be queried.
!>                                   LWORK = MAX(1, dimension) where
!>                                   dimension = (2KD+1)*N + KD*NTHREADS
!>                                   where KD is the size of the band.
!>                                   NTHREADS is the number of threads used when
!>                                   openMP compilation is enabled, otherwise =1.
!>          If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal sizes of the WORK, RWORK and
!>          IWORK arrays, returns these values as the first entries of
!>          the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (7*N)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (5*N)
!> \endverbatim
!>
!> \param[out] IFAIL
!> \verbatim
!>          IFAIL is INTEGER array, dimension (N)
!>          If JOBZ = 'V', then if INFO = 0, the first M elements of
!>          IFAIL are zero.  If INFO > 0, then IFAIL contains the
!>          indices of the eigenvectors that failed to converge.
!>          If JOBZ = 'N', then IFAIL is not referenced.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, then i eigenvectors failed to converge.
!>                Their indices are stored in array IFAIL.
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
!> \ingroup hbevx_2stage
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  All details about the 2stage techniques are available in:
!>
!>  Azzam Haidar, Hatem Ltaief, and Jack Dongarra.
!>  Parallel reduction to condensed forms for symmetric eigenvalue problems
!>  using aggregated fine-grained and memory-aware kernels. In Proceedings
!>  of 2011 International Conference for High Performance Computing,
!>  Networking, Storage and Analysis (SC '11), New York, NY, USA,
!>  Article 8 , 11 pages.
!>  http://doi.acm.org/10.1145/2063384.2063394
!>
!>  A. Haidar, J. Kurzak, P. Luszczek, 2013.
!>  An improved parallel singular value algorithm and its implementation
!>  for multicore hardware, In Proceedings of 2013 International Conference
!>  for High Performance Computing, Networking, Storage and Analysis (SC '13).
!>  Denver, Colorado, USA, 2013.
!>  Article 90, 12 pages.
!>  http://doi.acm.org/10.1145/2503210.2503292
!>
!>  A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra.
!>  A novel hybrid CPU-GPU generalized eigensolver for electronic structure
!>  calculations based on fine-grained memory aware tasks.
!>  International Journal of High Performance Computing Applications.
!>  Volume 28 Issue 2, Pages 196-209, May 2014.
!>  http://hpc.sagepub.com/content/28/2/196
!>
!> \endverbatim
!
!  =====================================================================
   SUBROUTINE CHBEVX_2STAGE( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, &
                             Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, &
                             Z, LDZ, WORK, LWORK, RWORK, IWORK, &
                             IFAIL, INFO )
!
   IMPLICIT NONE
!
!  -- LAPACK driver routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          JOBZ, RANGE, UPLO
   INTEGER            IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N, LWORK
   REAL               ABSTOL, VL, VU
!     ..
!     .. Array Arguments ..
   INTEGER            IFAIL( * ), IWORK( * )
   REAL               RWORK( * ), W( * )
   COMPLEX            AB( LDAB, * ), Q( LDQ, * ), WORK( * ), &
                      Z( LDZ, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO, ONE
   PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
   COMPLEX            CZERO, CONE
   PARAMETER          ( CZERO = ( 0.0E0, 0.0E0 ), &
                      CONE = ( 1.0E0, 0.0E0 ) )
!     ..
!     .. Local Scalars ..
   LOGICAL            ALLEIG, INDEIG, LOWER, TEST, VALEIG, WANTZ, &
                      LQUERY
   CHARACTER          ORDER
   INTEGER            I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL, &
                      INDISP, INDIWK, INDRWK, INDWRK, ISCALE, ITMP1, &
                      LLWORK, LWMIN, LHTRD, LWTRD, IB, INDHOUS, &
                      J, JJ, NSPLIT
   REAL               ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, &
                      SIGMA, SMLNUM, TMP1, VLL, VUU
   COMPLEX            CTMP1
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ILAENV2STAGE
   REAL               SLAMCH, CLANHB
   EXTERNAL           LSAME, SLAMCH, CLANHB, ILAENV2STAGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           SCOPY, SSCAL, SSTEBZ, SSTERF, XERBLA, CCOPY, &
                      CGEMV, CLACPY, CLASCL, CSTEIN, CSTEQR, &
                      CSWAP, CHETRD_HB2ST
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          REAL, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   WANTZ = LSAME( JOBZ, 'V' )
   ALLEIG = LSAME( RANGE, 'A' )
   VALEIG = LSAME( RANGE, 'V' )
   INDEIG = LSAME( RANGE, 'I' )
   LOWER = LSAME( UPLO, 'L' )
   LQUERY = ( LWORK == -1 )
!
   INFO = 0
   IF( .NOT.( LSAME( JOBZ, 'N' ) ) ) THEN
      INFO = -1
   ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
      INFO = -2
   ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
      INFO = -3
   ELSE IF( N < 0 ) THEN
      INFO = -4
   ELSE IF( KD < 0 ) THEN
      INFO = -5
   ELSE IF( LDAB < KD+1 ) THEN
      INFO = -7
   ELSE IF( WANTZ .AND. LDQ < MAX( 1, N ) ) THEN
      INFO = -9
   ELSE
      IF( VALEIG ) THEN
         IF( N > 0 .AND. VU <= VL ) &
            INFO = -11
      ELSE IF( INDEIG ) THEN
         IF( IL < 1 .OR. IL > MAX( 1, N ) ) THEN
            INFO = -12
         ELSE IF( IU < MIN( N, IL ) .OR. IU > N ) THEN
            INFO = -13
         END IF
      END IF
   END IF
   IF( INFO == 0 ) THEN
      IF( LDZ < 1 .OR. ( WANTZ .AND. LDZ < N ) ) &
         INFO = -18
   END IF
!
   IF( INFO == 0 ) THEN
      IF( N <= 1 ) THEN
         LWMIN = 1
         WORK( 1 ) = LWMIN
      ELSE
         IB    = ILAENV2STAGE( 2, 'CHETRD_HB2ST', JOBZ, &
                               N, KD, -1, -1 )
         LHTRD = ILAENV2STAGE( 3, 'CHETRD_HB2ST', JOBZ, &
                               N, KD, IB, -1 )
         LWTRD = ILAENV2STAGE( 4, 'CHETRD_HB2ST', JOBZ, &
                               N, KD, IB, -1 )
         LWMIN = LHTRD + LWTRD
         WORK( 1 )  = LWMIN
      ENDIF
!
      IF( LWORK < LWMIN .AND. .NOT.LQUERY ) &
         INFO = -20
   END IF
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CHBEVX_2STAGE', -INFO )
      RETURN
   ELSE IF( LQUERY ) THEN
      RETURN
   END IF
!
!     Quick return if possible
!
   M = 0
   IF( N == 0 ) &
      RETURN
!
   IF( N == 1 ) THEN
      M = 1
      IF( LOWER ) THEN
         CTMP1 = AB( 1, 1 )
      ELSE
         CTMP1 = AB( KD+1, 1 )
      END IF
      TMP1 = REAL( CTMP1 )
      IF( VALEIG ) THEN
         IF( .NOT.( VL < TMP1 .AND. VU >= TMP1 ) ) &
            M = 0
      END IF
      IF( M == 1 ) THEN
         W( 1 ) = REAL( CTMP1 )
         IF( WANTZ ) &
            Z( 1, 1 ) = CONE
      END IF
      RETURN
   END IF
!
!     Get machine constants.
!
   SAFMIN = SLAMCH( 'Safe minimum' )
   EPS    = SLAMCH( 'Precision' )
   SMLNUM = SAFMIN / EPS
   BIGNUM = ONE / SMLNUM
   RMIN   = SQRT( SMLNUM )
   RMAX   = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
!
!     Scale matrix to allowable range, if necessary.
!
   ISCALE = 0
   ABSTLL = ABSTOL
   IF( VALEIG ) THEN
      VLL = VL
      VUU = VU
   ELSE
      VLL = ZERO
      VUU = ZERO
   END IF
   ANRM = CLANHB( 'M', UPLO, N, KD, AB, LDAB, RWORK )
   IF( ANRM > ZERO .AND. ANRM < RMIN ) THEN
      ISCALE = 1
      SIGMA = RMIN / ANRM
   ELSE IF( ANRM > RMAX ) THEN
      ISCALE = 1
      SIGMA = RMAX / ANRM
   END IF
   IF( ISCALE == 1 ) THEN
      IF( LOWER ) THEN
         CALL CLASCL( 'B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO )
      ELSE
         CALL CLASCL( 'Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO )
      END IF
      IF( ABSTOL > 0 ) &
         ABSTLL = ABSTOL*SIGMA
      IF( VALEIG ) THEN
         VLL = VL*SIGMA
         VUU = VU*SIGMA
      END IF
   END IF
!
!     Call CHBTRD_HB2ST to reduce Hermitian band matrix to tridiagonal form.
!
   INDD = 1
   INDE = INDD + N
   INDRWK = INDE + N
!
   INDHOUS = 1
   INDWRK  = INDHOUS + LHTRD
   LLWORK  = LWORK - INDWRK + 1
!
   CALL CHETRD_HB2ST( 'N', JOBZ, UPLO, N, KD, AB, LDAB, &
                       RWORK( INDD ), RWORK( INDE ), WORK( INDHOUS ), &
                       LHTRD, WORK( INDWRK ), LLWORK, IINFO )
!
!     If all eigenvalues are desired and ABSTOL is less than or equal
!     to zero, then call SSTERF or CSTEQR.  If this fails for some
!     eigenvalue, then try SSTEBZ.
!
   TEST = .FALSE.
   IF (INDEIG) THEN
      IF (IL == 1 .AND. IU == N) THEN
         TEST = .TRUE.
      END IF
   END IF
   IF ((ALLEIG .OR. TEST) .AND. (ABSTOL <= ZERO)) THEN
      CALL SCOPY( N, RWORK( INDD ), 1, W, 1 )
      INDEE = INDRWK + 2*N
      IF( .NOT.WANTZ ) THEN
         CALL SCOPY( N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 )
         CALL SSTERF( N, W, RWORK( INDEE ), INFO )
      ELSE
         CALL CLACPY( 'A', N, N, Q, LDQ, Z, LDZ )
         CALL SCOPY( N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 )
         CALL CSTEQR( JOBZ, N, W, RWORK( INDEE ), Z, LDZ, &
                      RWORK( INDRWK ), INFO )
         IF( INFO == 0 ) THEN
            DO I = 1, N
               IFAIL( I ) = 0
            ENDDO
         END IF
      END IF
      IF( INFO == 0 ) THEN
         M = N
         GO TO 30
      END IF
      INFO = 0
   END IF
!
!     Otherwise, call SSTEBZ and, if eigenvectors are desired, CSTEIN.
!
   IF( WANTZ ) THEN
      ORDER = 'B'
   ELSE
      ORDER = 'E'
   END IF
   INDIBL = 1
   INDISP = INDIBL + N
   INDIWK = INDISP + N
   CALL SSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, &
                RWORK( INDD ), RWORK( INDE ), M, NSPLIT, W, &
                IWORK( INDIBL ), IWORK( INDISP ), RWORK( INDRWK ), &
                IWORK( INDIWK ), INFO )
!
   IF( WANTZ ) THEN
      CALL CSTEIN( N, RWORK( INDD ), RWORK( INDE ), M, W, &
                   IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, &
                   RWORK( INDRWK ), IWORK( INDIWK ), IFAIL, INFO )
!
!        Apply unitary matrix used in reduction to tridiagonal
!        form to eigenvectors returned by CSTEIN.
!
      DO J = 1, M
         CALL CCOPY( N, Z( 1, J ), 1, WORK( 1 ), 1 )
         CALL CGEMV( 'N', N, N, CONE, Q, LDQ, WORK, 1, CZERO, &
                     Z( 1, J ), 1 )
      ENDDO
   END IF
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
30 CONTINUE
   IF( ISCALE == 1 ) THEN
      IF( INFO == 0 ) THEN
         IMAX = M
      ELSE
         IMAX = INFO - 1
      END IF
      CALL SSCAL( IMAX, ONE / SIGMA, W, 1 )
   END IF
!
!     If eigenvalues are not in order, then sort them, along with
!     eigenvectors.
!
   IF( WANTZ ) THEN
      DO J = 1, M - 1
         I = 0
         TMP1 = W( J )
         DO JJ = J + 1, M
            IF( W( JJ ) < TMP1 ) THEN
               I = JJ
               TMP1 = W( JJ )
            END IF
         ENDDO
!
         IF( I /= 0 ) THEN
            ITMP1 = IWORK( INDIBL+I-1 )
            W( I ) = W( J )
            IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
            W( J ) = TMP1
            IWORK( INDIBL+J-1 ) = ITMP1
            CALL CSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
            IF( INFO /= 0 ) THEN
               ITMP1 = IFAIL( I )
               IFAIL( I ) = IFAIL( J )
               IFAIL( J ) = ITMP1
            END IF
         END IF
      ENDDO
   END IF
!
!     Set WORK(1) to optimal workspace size.
!
   WORK( 1 ) = LWMIN
!
   RETURN
!
!     End of CHBEVX_2STAGE
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        