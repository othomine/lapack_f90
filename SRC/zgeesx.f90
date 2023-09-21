!> \brief <b> ZGEESX computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGEESX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeesx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeesx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeesx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGEESX( JOBVS, SORT, SELECT, SENSE, N, A, LDA, SDIM, W,
!                          VS, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK,
!                          BWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBVS, SENSE, SORT
!       INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM
!       DOUBLE PRECISION   RCONDE, RCONDV
!       ..
!       .. Array Arguments ..
!       LOGICAL            BWORK( * )
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * )
!       ..
!       .. Function Arguments ..
!       LOGICAL            SELECT
!       EXTERNAL           SELECT
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGEESX computes for an N-by-N complex nonsymmetric matrix A, the
!> eigenvalues, the Schur form T, and, optionally, the matrix of Schur
!> vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).
!>
!> Optionally, it also orders the eigenvalues on the diagonal of the
!> Schur form so that selected eigenvalues are at the top left;
!> computes a reciprocal condition number for the average of the
!> selected eigenvalues (RCONDE); and computes a reciprocal condition
!> number for the right invariant subspace corresponding to the
!> selected eigenvalues (RCONDV).  The leading columns of Z form an
!> orthonormal basis for this invariant subspace.
!>
!> For further explanation of the reciprocal condition numbers RCONDE
!> and RCONDV, see Section 4.10 of the LAPACK Users' Guide (where
!> these quantities are called s and sep respectively).
!>
!> A complex matrix is in Schur form if it is upper triangular.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBVS
!> \verbatim
!>          JOBVS is CHARACTER*1
!>          = 'N': Schur vectors are not computed;
!>          = 'V': Schur vectors are computed.
!> \endverbatim
!>
!> \param[in] SORT
!> \verbatim
!>          SORT is CHARACTER*1
!>          Specifies whether or not to order the eigenvalues on the
!>          diagonal of the Schur form.
!>          = 'N': Eigenvalues are not ordered;
!>          = 'S': Eigenvalues are ordered (see SELECT).
!> \endverbatim
!>
!> \param[in] SELECT
!> \verbatim
!>          SELECT is a LOGICAL FUNCTION of one COMPLEX*16 argument
!>          SELECT must be declared EXTERNAL in the calling subroutine.
!>          If SORT = 'S', SELECT is used to select eigenvalues to order
!>          to the top left of the Schur form.
!>          If SORT = 'N', SELECT is not referenced.
!>          An eigenvalue W(j) is selected if SELECT(W(j)) is true.
!> \endverbatim
!>
!> \param[in] SENSE
!> \verbatim
!>          SENSE is CHARACTER*1
!>          Determines which reciprocal condition numbers are computed.
!>          = 'N': None are computed;
!>          = 'E': Computed for average of selected eigenvalues only;
!>          = 'V': Computed for selected right invariant subspace only;
!>          = 'B': Computed for both.
!>          If SENSE = 'E', 'V' or 'B', SORT must equal 'S'.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A. N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA, N)
!>          On entry, the N-by-N matrix A.
!>          On exit, A is overwritten by its Schur form T.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] SDIM
!> \verbatim
!>          SDIM is INTEGER
!>          If SORT = 'N', SDIM = 0.
!>          If SORT = 'S', SDIM = number of eigenvalues for which
!>                         SELECT is true.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX*16 array, dimension (N)
!>          W contains the computed eigenvalues, in the same order
!>          that they appear on the diagonal of the output Schur form T.
!> \endverbatim
!>
!> \param[out] VS
!> \verbatim
!>          VS is COMPLEX*16 array, dimension (LDVS,N)
!>          If JOBVS = 'V', VS contains the unitary matrix Z of Schur
!>          vectors.
!>          If JOBVS = 'N', VS is not referenced.
!> \endverbatim
!>
!> \param[in] LDVS
!> \verbatim
!>          LDVS is INTEGER
!>          The leading dimension of the array VS.  LDVS >= 1, and if
!>          JOBVS = 'V', LDVS >= N.
!> \endverbatim
!>
!> \param[out] RCONDE
!> \verbatim
!>          RCONDE is DOUBLE PRECISION
!>          If SENSE = 'E' or 'B', RCONDE contains the reciprocal
!>          condition number for the average of the selected eigenvalues.
!>          Not referenced if SENSE = 'N' or 'V'.
!> \endverbatim
!>
!> \param[out] RCONDV
!> \verbatim
!>          RCONDV is DOUBLE PRECISION
!>          If SENSE = 'V' or 'B', RCONDV contains the reciprocal
!>          condition number for the selected right invariant subspace.
!>          Not referenced if SENSE = 'N' or 'E'.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= max(1,2*N).
!>          Also, if SENSE = 'E' or 'V' or 'B', LWORK >= 2*SDIM*(N-SDIM),
!>          where SDIM is the number of selected eigenvalues computed by
!>          this routine.  Note that 2*SDIM*(N-SDIM) <= N*N/2. Note also
!>          that an error is only returned if LWORK < max(1,2*N), but if
!>          SENSE = 'E' or 'V' or 'B' this may not be large enough.
!>          For good performance, LWORK must generally be larger.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates upper bound on the optimal size of the
!>          array WORK, returns this value as the first entry of the WORK
!>          array, and no error message related to LWORK is issued by
!>          XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] BWORK
!> \verbatim
!>          BWORK is LOGICAL array, dimension (N)
!>          Not referenced if SORT = 'N'.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value.
!>          > 0: if INFO = i, and i is
!>             <= N: the QR algorithm failed to compute all the
!>                   eigenvalues; elements 1:ILO-1 and i+1:N of W
!>                   contain those eigenvalues which have converged; if
!>                   JOBVS = 'V', VS contains the transformation which
!>                   reduces A to its partially converged Schur form.
!>             = N+1: the eigenvalues could not be reordered because some
!>                   eigenvalues were too close to separate (the problem
!>                   is very ill-conditioned);
!>             = N+2: after reordering, roundoff changed values of some
!>                   complex eigenvalues so that leading eigenvalues in
!>                   the Schur form no longer satisfy SELECT=.TRUE.  This
!>                   could also be caused by underflow due to scaling.
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
!> \ingroup geesx
!
!  =====================================================================
   SUBROUTINE ZGEESX( JOBVS, SORT, SELECT, SENSE, N, A, LDA, SDIM, W, &
                      VS, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, &
                      BWORK, INFO )
!
!  -- LAPACK driver routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          JOBVS, SENSE, SORT
   INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM
   DOUBLE PRECISION   RCONDE, RCONDV
!     ..
!     .. Array Arguments ..
   LOGICAL            BWORK( * )
   DOUBLE PRECISION   RWORK( * )
   COMPLEX*16         A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * )
!     ..
!     .. Function Arguments ..
   LOGICAL            SELECT
   EXTERNAL           SELECT
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            LQUERY, SCALEA, WANTSB, WANTSE, WANTSN, WANTST, &
                      WANTSV, WANTVS
   INTEGER            HSWORK, I, IBAL, ICOND, IERR, IEVAL, IHI, ILO, &
                      ITAU, IWRK, LWRK, MAXWRK, MINWRK
   DOUBLE PRECISION   ANRM, BIGNUM, CSCALE, EPS, SMLNUM
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   DUM( 1 )
!     ..
!     .. External Subroutines ..
   EXTERNAL           DLASCL, XERBLA, ZCOPY, ZGEBAK, ZGEBAL, ZGEHRD, &
                      ZHSEQR, ZLACPY, ZLASCL, ZTRSEN, ZUNGHR
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ILAENV
   DOUBLE PRECISION   DLAMCH, ZLANGE
   EXTERNAL           LSAME, ILAENV, DLAMCH, ZLANGE
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
   INFO = 0
   WANTVS = LSAME( JOBVS, 'V' )
   WANTST = LSAME( SORT, 'S' )
   WANTSN = LSAME( SENSE, 'N' )
   WANTSE = LSAME( SENSE, 'E' )
   WANTSV = LSAME( SENSE, 'V' )
   WANTSB = LSAME( SENSE, 'B' )
   LQUERY = ( LWORK == -1 )
!
   IF( ( .NOT.WANTVS ) .AND. ( .NOT.LSAME( JOBVS, 'N' ) ) ) THEN
      INFO = -1
   ELSE IF( ( .NOT.WANTST ) .AND. ( .NOT.LSAME( SORT, 'N' ) ) ) THEN
      INFO = -2
   ELSE IF( .NOT.( WANTSN .OR. WANTSE .OR. WANTSV .OR. WANTSB ) .OR. &
            ( .NOT.WANTST .AND. .NOT.WANTSN ) ) THEN
      INFO = -4
   ELSE IF( N < 0 ) THEN
      INFO = -5
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -7
   ELSE IF( LDVS < 1 .OR. ( WANTVS .AND. LDVS < N ) ) THEN
      INFO = -11
   END IF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of real workspace needed at that point in the
!       code, as well as the preferred amount for good performance.
!       CWorkspace refers to complex workspace, and RWorkspace to real
!       workspace. NB refers to the optimal block size for the
!       immediately following subroutine, as returned by ILAENV.
!       HSWORK refers to the workspace preferred by ZHSEQR, as
!       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
!       the worst case.
!       If SENSE = 'E', 'V' or 'B', then the amount of workspace needed
!       depends on SDIM, which is computed by the routine ZTRSEN later
!       in the code.)
!
   IF( INFO == 0 ) THEN
      IF( N == 0 ) THEN
         MINWRK = 1
         LWRK = 1
      ELSE
         MAXWRK = N + N*ILAENV( 1, 'ZGEHRD', ' ', N, 1, N, 0 )
         MINWRK = 2*N
!
         CALL ZHSEQR( 'S', JOBVS, N, 1, N, A, LDA, W, VS, LDVS, &
                WORK, -1, IEVAL )
         HSWORK = INT( WORK( 1 ) )
!
         IF( .NOT.WANTVS ) THEN
            MAXWRK = MAX( MAXWRK, HSWORK )
         ELSE
            MAXWRK = MAX( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'ZUNGHR', &
                          ' ', N, 1, N, -1 ) )
            MAXWRK = MAX( MAXWRK, HSWORK )
         END IF
         LWRK = MAXWRK
         IF( .NOT.WANTSN ) &
            LWRK = MAX( LWRK, ( N*N )/2 )
      END IF
      WORK( 1 ) = LWRK
!
      IF( LWORK < MINWRK .AND. .NOT.LQUERY ) THEN
         INFO = -15
      END IF
   END IF
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'ZGEESX', -INFO )
      RETURN
   ELSE IF( LQUERY ) THEN
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) THEN
      SDIM = 0
      RETURN
   END IF
!
!     Get machine constants
!
   EPS = DLAMCH( 'P' )
   SMLNUM = DLAMCH( 'S' )
   BIGNUM = ONE / SMLNUM
   SMLNUM = SQRT( SMLNUM ) / EPS
   BIGNUM = ONE / SMLNUM
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
   ANRM = ZLANGE( 'M', N, N, A, LDA, DUM )
   SCALEA = .FALSE.
   IF( ANRM > ZERO .AND. ANRM < SMLNUM ) THEN
      SCALEA = .TRUE.
      CSCALE = SMLNUM
   ELSE IF( ANRM > BIGNUM ) THEN
      SCALEA = .TRUE.
      CSCALE = BIGNUM
   END IF
   IF( SCALEA ) &
      CALL ZLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
!
!
!     Permute the matrix to make it more nearly triangular
!     (CWorkspace: none)
!     (RWorkspace: need N)
!
   IBAL = 1
   CALL ZGEBAL( 'P', N, A, LDA, ILO, IHI, RWORK( IBAL ), IERR )
!
!     Reduce to upper Hessenberg form
!     (CWorkspace: need 2*N, prefer N+N*NB)
!     (RWorkspace: none)
!
   ITAU = 1
   IWRK = N + ITAU
   CALL ZGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), &
                LWORK-IWRK+1, IERR )
!
   IF( WANTVS ) THEN
!
!        Copy Householder vectors to VS
!
      CALL ZLACPY( 'L', N, N, A, LDA, VS, LDVS )
!
!        Generate unitary matrix in VS
!        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
!        (RWorkspace: none)
!
      CALL ZUNGHR( N, ILO, IHI, VS, LDVS, WORK( ITAU ), WORK( IWRK ), &
                   LWORK-IWRK+1, IERR )
   END IF
!
   SDIM = 0
!
!     Perform QR iteration, accumulating Schur vectors in VS if desired
!     (CWorkspace: need 1, prefer HSWORK (see comments) )
!     (RWorkspace: none)
!
   IWRK = ITAU
   CALL ZHSEQR( 'S', JOBVS, N, ILO, IHI, A, LDA, W, VS, LDVS, &
                WORK( IWRK ), LWORK-IWRK+1, IEVAL )
   IF( IEVAL > 0 ) &
      INFO = IEVAL
!
!     Sort eigenvalues if desired
!
   IF( WANTST .AND. INFO == 0 ) THEN
      IF( SCALEA ) &
         CALL ZLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, W, N, IERR )
      DO I = 1, N
         BWORK( I ) = SELECT( W( I ) )
      ENDDO
!
!        Reorder eigenvalues, transform Schur vectors, and compute
!        reciprocal condition numbers
!        (CWorkspace: if SENSE is not 'N', need 2*SDIM*(N-SDIM)
!                     otherwise, need none )
!        (RWorkspace: none)
!
      CALL ZTRSEN( SENSE, JOBVS, BWORK, N, A, LDA, VS, LDVS, W, SDIM, &
                   RCONDE, RCONDV, WORK( IWRK ), LWORK-IWRK+1, &
                   ICOND )
      IF( .NOT.WANTSN ) &
         MAXWRK = MAX( MAXWRK, 2*SDIM*( N-SDIM ) )
      IF( ICOND == -14 ) THEN
!
!           Not enough complex workspace
!
         INFO = -15
      END IF
   END IF
!
   IF( WANTVS ) THEN
!
!        Undo balancing
!        (CWorkspace: none)
!        (RWorkspace: need N)
!
      CALL ZGEBAK( 'P', 'R', N, ILO, IHI, RWORK( IBAL ), N, VS, LDVS, &
                   IERR )
   END IF
!
   IF( SCALEA ) THEN
!
!        Undo scaling for the Schur form of A
!
      CALL ZLASCL( 'U', 0, 0, CSCALE, ANRM, N, N, A, LDA, IERR )
      CALL ZCOPY( N, A, LDA+1, W, 1 )
      IF( ( WANTSV .OR. WANTSB ) .AND. INFO == 0 ) THEN
         DUM( 1 ) = RCONDV
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, 1, 1, DUM, 1, IERR )
         RCONDV = DUM( 1 )
      END IF
   END IF
!
   WORK( 1 ) = MAXWRK
   RETURN
!
!     End of ZGEESX
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        