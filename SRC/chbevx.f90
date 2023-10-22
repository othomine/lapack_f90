!> \brief <b> CHBEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHBEVX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbevx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbevx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbevx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHBEVX( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, LDQ, VL,
!                          VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK,
!                          IWORK, IFAIL, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, RANGE, UPLO
!       INTEGER            IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N
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
!> CHBEVX computes selected eigenvalues and, optionally, eigenvectors
!> of a complex Hermitian band matrix A.  Eigenvalues and eigenvectors
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
!>          WORK is COMPLEX array, dimension (N)
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
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
!
!> \ingroup hbevx
!
!  =====================================================================
   SUBROUTINE CHBEVX( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, LDQ, VL, &
                      VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK, &
                      IWORK, IFAIL, INFO )
!
!  -- LAPACK driver routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          JOBZ, RANGE, UPLO
   INTEGER            IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N
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
!     ..
!     .. Local Scalars ..
   LOGICAL            ALLEIG, INDEIG, LOWER, TEST, VALEIG, WANTZ
   CHARACTER          ORDER
   INTEGER            I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL, &
                      INDISP, INDIWK, INDRWK, INDWRK, ISCALE, ITMP1, &
                      J, JJ, NSPLIT
   REAL               ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, &
                      SIGMA, SMLNUM, TMP1, VLL, VUU
   COMPLEX            CTMP1, Z_TMP( LDZ )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               CLANHB, SLAMCH
   EXTERNAL           LSAME, CLANHB, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMV, CHBTRD, CLACPY, CLASCL, CSTEIN, &
                      CSTEQR, SSTEBZ, SSTERF, XERBLA
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
!
   INFO = 0
   IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
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
      IF( LDZ < 1 .OR. ( WANTZ .AND. LDZ < N ) ) INFO = -18
   END IF
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CHBEVX', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   M = 0
   IF( N == 0 ) RETURN
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
         IF( .NOT.( VL < TMP1 .AND. VU >= TMP1 ) ) M = 0
      END IF
      IF( M == 1 ) THEN
         W( 1 ) = REAL( CTMP1 )
         IF( WANTZ ) Z( 1, 1 ) = (1.0E+0,0.0E+0)
      END IF
      RETURN
   END IF
!
!     Get machine constants.
!
   SAFMIN = SLAMCH( 'Safe minimum' )
   EPS = SLAMCH( 'Precision' )
   SMLNUM = SAFMIN / EPS
   BIGNUM = 1.0E+0 / SMLNUM
   RMIN = SQRT( SMLNUM )
   RMAX = MIN( SQRT( BIGNUM ), 1.0E+0 / SQRT( SQRT( SAFMIN ) ) )
!
!     Scale matrix to allowable range, if necessary.
!
   ISCALE = 0
   ABSTLL = ABSTOL
   IF ( VALEIG ) THEN
      VLL = VL
      VUU = VU
   ELSE
      VLL = 0.0E+0
      VUU = 0.0E+0
   ENDIF
   ANRM = CLANHB( 'M', UPLO, N, KD, AB, LDAB, RWORK )
   IF( ANRM > 0.0E+0 .AND. ANRM < RMIN ) THEN
      ISCALE = 1
      SIGMA = RMIN / ANRM
   ELSE IF( ANRM > RMAX ) THEN
      ISCALE = 1
      SIGMA = RMAX / ANRM
   END IF
   IF( ISCALE == 1 ) THEN
      IF( LOWER ) THEN
         CALL CLASCL( 'B', KD, KD, 1.0E+0, SIGMA, N, N, AB, LDAB, INFO )
      ELSE
         CALL CLASCL( 'Q', KD, KD, 1.0E+0, SIGMA, N, N, AB, LDAB, INFO )
      END IF
      IF( ABSTOL > 0 ) &
         ABSTLL = ABSTOL*SIGMA
      IF( VALEIG ) THEN
         VLL = VL*SIGMA
         VUU = VU*SIGMA
      END IF
   END IF
!
!     Call CHBTRD to reduce Hermitian band matrix to tridiagonal form.
!
   INDD = 1
   INDE = INDD + N
   INDRWK = INDE + N
   INDWRK = 1
   CALL CHBTRD( JOBZ, UPLO, N, KD, AB, LDAB, RWORK( INDD ), &
                RWORK( INDE ), Q, LDQ, WORK( INDWRK ), IINFO )
!
!     If all eigenvalues are desired and ABSTOL is less than or equal
!     to zero, then call SSTERF or CSTEQR.  If this fails for some
!     eigenvalue, then try SSTEBZ.
!
   TEST = .FALSE.
   IF (INDEIG) THEN
      IF (IL == 1 .AND. IU == N) TEST = .TRUE.
   END IF
   IF ((ALLEIG .OR. TEST) .AND. (ABSTOL <= 0.0E+0)) THEN
      W(1:N) = RWORK( INDD:INDD+N-1 )
      INDEE = INDRWK + 2*N
      RWORK( INDEE:INDEE+N-2 ) = RWORK( INDE:INDE+N-2 )
      IF( .NOT.WANTZ ) THEN
         CALL SSTERF( N, W, RWORK( INDEE ), INFO )
      ELSE
         Z(1:N,1:N) = Q(1:N,1:N)
         CALL CSTEQR( JOBZ, N, W, RWORK( INDEE ), Z, LDZ, &
                      RWORK( INDRWK ), INFO )
         IF( INFO == 0 ) IFAIL(1:N) = 0
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
         WORK(1:N) = Z(1:N,J)
         CALL CGEMV( 'N', N, N, (1.0E+0,0.0E+0), Q, LDQ, WORK, 1, (0.0E+0,0.0E+0), &
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
      W(1:IMAX) = W(1:IMAX)/SIGMA
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
            Z_TMP(1:N) = Z(1:N,I)
            Z(1:N,I) = Z(1:N,J)
            Z(1:N,J) = Z_TMP(1:N)
            IF( INFO /= 0 ) THEN
               ITMP1 = IFAIL( I )
               IFAIL( I ) = IFAIL( J )
               IFAIL( J ) = ITMP1
            END IF
         END IF
      ENDDO
   END IF
!
   RETURN
!
!     End of CHBEVX
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

