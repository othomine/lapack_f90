!> \brief \b DGSVJ1 pre-processor for the routine dgesvj, applies Jacobi rotations targeting only particular pivots.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGSVJ1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgsvj1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgsvj1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgsvj1.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGSVJ1( JOBV, M, N, N1, A, LDA, D, SVA, MV, V, LDV,
!                          EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   EPS, SFMIN, TOL
!       INTEGER            INFO, LDA, LDV, LWORK, M, MV, N, N1, NSWEEP
!       CHARACTER*1        JOBV
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), D( N ), SVA( N ), V( LDV, * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGSVJ1 is called from DGESVJ as a pre-processor and that is its main
!> purpose. It applies Jacobi rotations in the same way as DGESVJ does, but
!> it targets only particular pivots and it does not check convergence
!> (stopping criterion). Few tuning parameters (marked by [TP]) are
!> available for the implementer.
!>
!> Further Details
!> ~~~~~~~~~~~~~~~
!> DGSVJ1 applies few sweeps of Jacobi rotations in the column space of
!> the input M-by-N matrix A. The pivot pairs are taken from the (1,2)
!> off-diagonal block in the corresponding N-by-N Gram matrix A^T * A. The
!> block-entries (tiles) of the (1,2) off-diagonal block are marked by the
!> [x]'s in the following scheme:
!>
!>    | *  *  * [x] [x] [x]|
!>    | *  *  * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.
!>    | *  *  * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
!>    |[x] [x] [x] *  *  * |
!>    |[x] [x] [x] *  *  * |
!>    |[x] [x] [x] *  *  * |
!>
!> In terms of the columns of A, the first N1 columns are rotated 'against'
!> the remaining N-N1 columns, trying to increase the angle between the
!> corresponding subspaces. The off-diagonal block is N1-by(N-N1) and it is
!> tiled using quadratic tiles of side KBL. Here, KBL is a tuning parameter.
!> The number of sweeps is given in NSWEEP and the orthogonality threshold
!> is given in TOL.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBV
!> \verbatim
!>          JOBV is CHARACTER*1
!>          Specifies whether the output from this procedure is used
!>          to compute the matrix V:
!>          = 'V': the product of the Jacobi rotations is accumulated
!>                 by postmultiplying the N-by-N array V.
!>                (See the description of V.)
!>          = 'A': the product of the Jacobi rotations is accumulated
!>                 by postmultiplying the MV-by-N array V.
!>                (See the descriptions of MV and V.)
!>          = 'N': the Jacobi rotations are not accumulated.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the input matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the input matrix A.
!>          M >= N >= 0.
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!>          N1 specifies the 2 x 2 block partition, the first N1 columns are
!>          rotated 'against' the remaining N-N1 columns of A.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, M-by-N matrix A, such that A*diag(D) represents
!>          the input matrix.
!>          On exit,
!>          A_onexit * D_onexit represents the input matrix A*diag(D)
!>          post-multiplied by a sequence of Jacobi rotations, where the
!>          rotation threshold and the total number of sweeps are given in
!>          TOL and NSWEEP, respectively.
!>          (See the descriptions of N1, D, TOL and NSWEEP.)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The array D accumulates the scaling factors from the fast scaled
!>          Jacobi rotations.
!>          On entry, A*diag(D) represents the input matrix.
!>          On exit, A_onexit*diag(D_onexit) represents the input matrix
!>          post-multiplied by a sequence of Jacobi rotations, where the
!>          rotation threshold and the total number of sweeps are given in
!>          TOL and NSWEEP, respectively.
!>          (See the descriptions of N1, A, TOL and NSWEEP.)
!> \endverbatim
!>
!> \param[in,out] SVA
!> \verbatim
!>          SVA is DOUBLE PRECISION array, dimension (N)
!>          On entry, SVA contains the Euclidean norms of the columns of
!>          the matrix A*diag(D).
!>          On exit, SVA contains the Euclidean norms of the columns of
!>          the matrix onexit*diag(D_onexit).
!> \endverbatim
!>
!> \param[in] MV
!> \verbatim
!>          MV is INTEGER
!>          If JOBV = 'A', then MV rows of V are post-multiplied by a
!>                         sequence of Jacobi rotations.
!>          If JOBV = 'N', then MV is not referenced.
!> \endverbatim
!>
!> \param[in,out] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension (LDV,N)
!>          If JOBV = 'V', then N rows of V are post-multiplied by a
!>                         sequence of Jacobi rotations.
!>          If JOBV = 'A', then MV rows of V are post-multiplied by a
!>                         sequence of Jacobi rotations.
!>          If JOBV = 'N', then V is not referenced.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V,  LDV >= 1.
!>          If JOBV = 'V', LDV >= N.
!>          If JOBV = 'A', LDV >= MV.
!> \endverbatim
!>
!> \param[in] EPS
!> \verbatim
!>          EPS is DOUBLE PRECISION
!>          EPS = DLAMCH('Epsilon')
!> \endverbatim
!>
!> \param[in] SFMIN
!> \verbatim
!>          SFMIN is DOUBLE PRECISION
!>          SFMIN = DLAMCH('Safe Minimum')
!> \endverbatim
!>
!> \param[in] TOL
!> \verbatim
!>          TOL is DOUBLE PRECISION
!>          TOL is the threshold for Jacobi rotations. For a pair
!>          A(:,p), A(:,q) of pivot columns, the Jacobi rotation is
!>          applied only if DABS(COS(angle(A(:,p),A(:,q)))) > TOL.
!> \endverbatim
!>
!> \param[in] NSWEEP
!> \verbatim
!>          NSWEEP is INTEGER
!>          NSWEEP is the number of sweeps of Jacobi rotations to be
!>          performed.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          LWORK is the dimension of WORK. LWORK >= M.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, then the i-th argument had an illegal value
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
!> \ingroup gsvj1
!
!> \par Contributors:
!  ==================
!>
!> Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany)
!
!  =====================================================================
   SUBROUTINE DGSVJ1( JOBV, M, N, N1, A, LDA, D, SVA, MV, V, LDV, &
                      EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   DOUBLE PRECISION   EPS, SFMIN, TOL
   INTEGER            INFO, LDA, LDV, LWORK, M, MV, N, N1, NSWEEP
   CHARACTER*1        JOBV
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   A( LDA, * ), D( N ), SVA( N ), V( LDV, * ), &
                      WORK( LWORK )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   DOUBLE PRECISION   AAPP, AAPP0, AAPQ, AAQQ, APOAQ, AQOAP, BIG, &
                      BIGTHETA, CS, LARGE, MXAAPQ, MXSINJ, ROOTBIG, &
                      ROOTEPS, ROOTSFMIN, ROOTTOL, SMALL, SN, T, &
                      TEMP1, THETA, THSIGN
   INTEGER            BLSKIP, EMPTSW, i, ibr, igl, IERR, IJBLSK, &
                      ISWROT, jbc, jgl, KBL, MVL, NOTROT, nblc, nblr, &
                      p, PSKIPPED, q, ROWSKIP, SWBAND
   LOGICAL            APPLV, ROTOK, RSVEC
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   FASTR( 5 )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DNRM2
   INTEGER            IDAMAX
   LOGICAL            LSAME
   EXTERNAL           IDAMAX, LSAME, DNRM2
!     ..
!     .. External Subroutines ..
   EXTERNAL           DLASCL, DLASSQ, DROTM, DSWAP, &
                      XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   APPLV = LSAME( JOBV, 'A' )
   RSVEC = LSAME( JOBV, 'V' )
   IF( .NOT.( RSVEC .OR. APPLV .OR. LSAME( JOBV, 'N' ) ) ) THEN
      INFO = -1
   ELSE IF( M < 0 ) THEN
      INFO = -2
   ELSE IF( ( N < 0 ) .OR. ( N > M ) ) THEN
      INFO = -3
   ELSE IF( N1 < 0 ) THEN
      INFO = -4
   ELSE IF( LDA < M ) THEN
      INFO = -6
   ELSE IF( ( RSVEC.OR.APPLV ) .AND. ( MV < 0 ) ) THEN
      INFO = -9
   ELSE IF( ( RSVEC.AND.( LDV < N ) ).OR. &
            ( APPLV.AND.( LDV < MV ) )  ) THEN
      INFO = -11
   ELSE IF( TOL <= EPS ) THEN
      INFO = -14
   ELSE IF( NSWEEP < 0 ) THEN
      INFO = -15
   ELSE IF( LWORK < M ) THEN
      INFO = -17
   ELSE
      INFO = 0
   END IF
!
!     #:(
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'DGSVJ1', -INFO )
      RETURN
   END IF
!
   IF( RSVEC ) THEN
      MVL = N
   ELSE IF( APPLV ) THEN
      MVL = MV
   END IF
   RSVEC = RSVEC .OR. APPLV

   ROOTEPS = DSQRT( EPS )
   ROOTSFMIN = DSQRT( SFMIN )
   SMALL = SFMIN / EPS
   BIG = 1.0D0 / SFMIN
   ROOTBIG = 1.0D0 / ROOTSFMIN
   LARGE = BIG / DSQRT( DBLE( M*N ) )
   BIGTHETA = 1.0D0 / ROOTEPS
   ROOTTOL = DSQRT( TOL )
!
!     .. Initialize the right singular vector matrix ..
!
!     RSVEC = LSAME( JOBV, 'Y' )
!
   EMPTSW = N1*( N-N1 )
   NOTROT = 0
   FASTR( 1 ) = 0.0D0
!
!     .. Row-cyclic pivot strategy with de Rijk's pivoting ..
!
   KBL = MIN( 8, N )
   NBLR = N1 / KBL
   IF( ( NBLR*KBL ) /= N1 )NBLR = NBLR + 1

!     .. the tiling is nblr-by-nblc [tiles]

   NBLC = ( N-N1 ) / KBL
   IF( ( NBLC*KBL ) /= ( N-N1 ) )NBLC = NBLC + 1
   BLSKIP = ( KBL**2 ) + 1
![TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.

   ROWSKIP = MIN( 5, KBL )
![TP] ROWSKIP is a tuning parameter.
   SWBAND = 0
![TP] SWBAND is a tuning parameter. It is meaningful and effective
!     if SGESVJ is used as a computational routine in the preconditioned
!     Jacobi SVD algorithm SGESVJ.
!
!
!     | *   *   * [x] [x] [x]|
!     | *   *   * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.
!     | *   *   * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
!     |[x] [x] [x] *   *   * |
!     |[x] [x] [x] *   *   * |
!     |[x] [x] [x] *   *   * |
!
!
   DO i = 1, NSWEEP
!     .. go go go ...
!
      MXAAPQ = 0.0D0
      MXSINJ = 0.0D0
      ISWROT = 0
!
      NOTROT = 0
      PSKIPPED = 0
!
      DO ibr = 1, NBLR

         igl = ( ibr-1 )*KBL + 1
!
!
!........................................................
! ... go to the off diagonal blocks

         igl = ( ibr-1 )*KBL + 1

         DO jbc = 1, NBLC

            jgl = N1 + ( jbc-1 )*KBL + 1

!        doing the block at ( ibr, jbc )

            IJBLSK = 0
            DO p = igl, MIN( igl+KBL-1, N1 )

               AAPP = SVA( p )

               IF( AAPP > 0.0D0 ) THEN

                  PSKIPPED = 0

                  DO q = jgl, MIN( jgl+KBL-1, N )
!
                     AAQQ = SVA( q )

                     IF( AAQQ > 0.0D0 ) THEN
                        AAPP0 = AAPP
!
!     .. M x 2 Jacobi SVD ..
!
!        .. Safe Gram matrix computation ..
!
                        IF( AAQQ >= 1.0D0 ) THEN
                           IF( AAPP >= AAQQ ) THEN
                              ROTOK = ( SMALL*AAPP ) <= AAQQ
                           ELSE
                              ROTOK = ( SMALL*AAQQ ) <= AAPP
                           END IF
                           IF( AAPP < ( BIG / AAQQ ) ) THEN
                              AAPQ = SUM(A(1:M,p)*A(1:M,q))*D( p )*D( q ) / AAQQ / AAPP
                           ELSE
                              WORK(1:M) = A(1:M,p)
                              CALL DLASCL( 'G', 0, 0, AAPP, D( p ), &
                                           M, 1, WORK, LDA, IERR )
                              AAPQ = SUM(WORK(1:M)*A(1:M,q))*D( q ) / AAQQ
                           END IF
                        ELSE
                           IF( AAPP >= AAQQ ) THEN
                              ROTOK = AAPP <= ( AAQQ / SMALL )
                           ELSE
                              ROTOK = AAQQ <= ( AAPP / SMALL )
                           END IF
                           IF( AAPP > ( SMALL / AAQQ ) ) THEN
                              AAPQ = SUM(A(1:M,p)*A(1:M,q))*D( p )*D( q ) / AAQQ / AAPP
                           ELSE
                              WORK(1:M) = A(1:M,q)
                              CALL DLASCL( 'G', 0, 0, AAQQ, D( q ), &
                                           M, 1, WORK, LDA, IERR )
                              AAPQ = SUM(WORK(1:M)*A(1:M,p))*D( p ) / AAPP
                           END IF
                        END IF

                        MXAAPQ = MAX( MXAAPQ, DABS( AAPQ ) )

!        TO rotate or NOT to rotate, THAT is the question ...
!
                        IF( DABS( AAPQ ) > TOL ) THEN
                           NOTROT = 0
!           ROTATED  = ROTATED + 1
                           PSKIPPED = 0
                           ISWROT = ISWROT + 1
!
                           IF( ROTOK ) THEN
!
                              AQOAP = AAQQ / AAPP
                              APOAQ = AAPP / AAQQ
                              THETA = -0.5D0*DABS(AQOAP-APOAQ) / AAPQ
                              IF( AAQQ > AAPP0 )THETA = -THETA

                              IF( DABS( THETA ) > BIGTHETA ) THEN
                                 T = 0.5D0 / THETA
                                 FASTR( 3 ) = T*D( p ) / D( q )
                                 FASTR( 4 ) = -T*D( q ) / D( p )
                                 CALL DROTM( M, A( 1, p ), 1, &
                                             A( 1, q ), 1, FASTR )
                                 IF( RSVEC )CALL DROTM( MVL, &
                                                 V( 1, p ), 1, &
                                                 V( 1, q ), 1, &
                                                 FASTR )
                                 SVA( q ) = AAQQ*DSQRT( MAX( 0.0D0, &
                                            1.0D0+T*APOAQ*AAPQ ) )
                                 AAPP = AAPP*DSQRT( MAX( 0.0D0, &
                                        1.0D0-T*AQOAP*AAPQ ) )
                                 MXSINJ = MAX( MXSINJ, DABS( T ) )
                              ELSE
!
!                 .. choose correct signum for THETA and rotate
!
                                 THSIGN = -DSIGN( 1.0D0, AAPQ )
                                 IF( AAQQ > AAPP0 )THSIGN = -THSIGN
                                 T = 1.0D0 / ( THETA+THSIGN* &
                                     DSQRT( 1.0D0+THETA*THETA ) )
                                 CS = DSQRT( 1.0D0 / ( 1.0D0+T*T ) )
                                 SN = T*CS
                                 MXSINJ = MAX( MXSINJ, DABS( SN ) )
                                 SVA( q ) = AAQQ*DSQRT( MAX( 0.0D0, &
                                            1.0D0+T*APOAQ*AAPQ ) )
                                 AAPP = AAPP*DSQRT( MAX( 0.0D0, &
                                       1.0D0-T*AQOAP*AAPQ ) )

                                 APOAQ = D( p ) / D( q )
                                 AQOAP = D( q ) / D( p )
                                 IF( D( p ) >= 1.0D0 ) THEN
!
                                    IF( D( q ) >= 1.0D0 ) THEN
                                       FASTR( 3 ) = T*APOAQ
                                       FASTR( 4 ) = -T*AQOAP
                                       D( p ) = D( p )*CS
                                       D( q ) = D( q )*CS
                                       CALL DROTM( M, A( 1, p ), 1, &
                                                   A( 1, q ), 1, &
                                                   FASTR )
                                       IF( RSVEC )CALL DROTM( MVL, &
                                           V( 1, p ), 1, V( 1, q ), &
                                           1, FASTR )
                                    ELSE
                                       A(1:M,p) = A(1:M,p)-T*AQOAP*A(1:M,q)
                                       A(1:M,q) = A(1:M,q)+CS*SN*APOAQ*A(1:M,p)
                                       IF( RSVEC ) THEN
                                          V(1:MVL,p) = V(1:MVL,p)-T*AQOAP*V(1:MVL,q)
                                          V(1:MVL,q) = V(1:MVL,q)+CS*SN*APOAQ*V(1:MVL,p)
                                       END IF
                                       D( p ) = D( p )*CS
                                       D( q ) = D( q ) / CS
                                    END IF
                                 ELSE
                                    IF( D( q ) >= 1.0D0 ) THEN
                                       A(1:M,q) = A(1:M,q)+T*APOAQ*A(1:M,p)
                                       A(1:M,p) = A(1:M,p)-CS*SN*AQOAP*A(1:M,q)
                                       IF( RSVEC ) THEN
                                          V(1:MVL,q) = V(1:MVL,q)+T*APOAQ*V(1:MVL,p)
                                          V(1:MVL,p) = V(1:MVL,p)-CS*SN*AQOAP*V(1:MVL,q)
                                       END IF
                                       D( p ) = D( p ) / CS
                                       D( q ) = D( q )*CS
                                    ELSE
                                       IF( D( p ) >= D( q ) ) THEN
                                          A(1:M,p) = A(1:M,p)-T*AQOAP*A(1:M,q)
                                          A(1:M,q) = A(1:M,q)+CS*SN*APOAQ*A(1:M,p)
                                          D( p ) = D( p )*CS
                                          D( q ) = D( q ) / CS
                                          IF( RSVEC ) THEN
                                             V(1:MVL,p) = V(1:MVL,p)-T*AQOAP*V(1:MVL,q)
                                             V(1:MVL,q) = V(1:MVL,q)+CS*SN*APOAQ*V(1:MVL,p)
                                          END IF
                                       ELSE
                                          A(1:M,q) = A(1:M,q)+T*APOAQ*A(1:M,p)
                                          A(1:M,p) = A(1:M,p)-CS*SN*AQOAP*A(1:M,q)
                                          D( p ) = D( p ) / CS
                                          D( q ) = D( q )*CS
                                          IF( RSVEC ) THEN
                                             V(1:MVL,q) = V(1:MVL,q)+T*APOAQ*V(1:MVL,p)
                                             V(1:MVL,p) = V(1:MVL,p)-CS*SN*AQOAP*V(1:MVL,q)
                                          END IF
                                       END IF
                                    END IF
                                 END IF
                              END IF

                           ELSE
                              IF( AAPP > AAQQ ) THEN
                                 WORK(1:M) = A(1:M,p)
                                 CALL DLASCL( 'G', 0, 0, AAPP, 1.0D0, &
                                              M, 1, WORK, LDA, IERR )
                                 CALL DLASCL( 'G', 0, 0, AAQQ, 1.0D0, &
                                              M, 1, A( 1, q ), LDA, &
                                              IERR )
                                 TEMP1 = -AAPQ*D( p ) / D( q )
                                 A(1:M,q) = A(1:M,q) + TEMP1*WORK(1:M)
                                 CALL DLASCL( 'G', 0, 0, 1.0D0, AAQQ, &
                                              M, 1, A( 1, q ), LDA, &
                                              IERR )
                                 SVA( q ) = AAQQ*DSQRT( MAX( 0.0D0, &
                                            1.0D0-AAPQ*AAPQ ) )
                                 MXSINJ = MAX( MXSINJ, SFMIN )
                              ELSE
                                 WORK(1:M) = A(1:M,q)
                                 CALL DLASCL( 'G', 0, 0, AAQQ, 1.0D0, &
                                              M, 1, WORK, LDA, IERR )
                                 CALL DLASCL( 'G', 0, 0, AAPP, 1.0D0, &
                                              M, 1, A( 1, p ), LDA, &
                                              IERR )
                                 TEMP1 = -AAPQ*D( q ) / D( p )
                                 A(1:M,p) = A(1:M,p) + TEMP1*WORK(1:M)
                                 CALL DLASCL( 'G', 0, 0, 1.0D0, AAPP, &
                                              M, 1, A( 1, p ), LDA, &
                                              IERR )
                                 SVA( p ) = AAPP*DSQRT( MAX( 0.0D0, &
                                            1.0D0-AAPQ*AAPQ ) )
                                 MXSINJ = MAX( MXSINJ, SFMIN )
                              END IF
                           END IF
!           END IF ROTOK THEN ... ELSE
!
!           In the case of cancellation in updating SVA(q)
!           .. recompute SVA(q)
                           IF( ( SVA( q ) / AAQQ )**2 <= ROOTEPS ) &
                               THEN
                              IF( ( AAQQ < ROOTBIG ) .AND. &
                                  ( AAQQ > ROOTSFMIN ) ) THEN
                                 SVA( q ) = DNRM2( M, A( 1, q ), 1 )* &
                                            D( q )
                              ELSE
                                 T = 0.0D0
                                 AAQQ = 1.0D0
                                 CALL DLASSQ( M, A( 1, q ), 1, T, &
                                              AAQQ )
                                 SVA( q ) = T*DSQRT( AAQQ )*D( q )
                              END IF
                           END IF
                           IF( ( AAPP / AAPP0 )**2 <= ROOTEPS ) THEN
                              IF( ( AAPP < ROOTBIG ) .AND. &
                                  ( AAPP > ROOTSFMIN ) ) THEN
                                 AAPP = DNRM2( M, A( 1, p ), 1 )* &
                                        D( p )
                              ELSE
                                 T = 0.0D0
                                 AAPP = 1.0D0
                                 CALL DLASSQ( M, A( 1, p ), 1, T, &
                                              AAPP )
                                 AAPP = T*DSQRT( AAPP )*D( p )
                              END IF
                              SVA( p ) = AAPP
                           END IF
!              end of OK rotation
                        ELSE
                           NOTROT = NOTROT + 1
!           SKIPPED  = SKIPPED  + 1
                           PSKIPPED = PSKIPPED + 1
                           IJBLSK = IJBLSK + 1
                        END IF
                     ELSE
                        NOTROT = NOTROT + 1
                        PSKIPPED = PSKIPPED + 1
                        IJBLSK = IJBLSK + 1
                     END IF

!      IF ( NOTROT  >=  EMPTSW )  GO TO 2011
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
                        GO TO 2203
                     END IF

!
                     ENDDO
!        end of the q-loop
 2203                CONTINUE

                  SVA( p ) = AAPP
!
               ELSE
                  IF( AAPP == 0.0D0 )NOTROT = NOTROT + &
                      MIN( jgl+KBL-1, N ) - jgl + 1
                  IF( AAPP < 0.0D0 )NOTROT = 0
!**      IF ( NOTROT  >=  EMPTSW )  GO TO 2011
               END IF

               ENDDO
!     end of the p-loop
            ENDDO
!     end of the jbc-loop
 2011       CONTINUE
!2011 bailed out of the jbc-loop
         DO p = igl, MIN( igl+KBL-1, N )
            SVA( p ) = DABS( SVA( p ) )
         ENDDO
!**   IF ( NOTROT  >=  EMPTSW ) GO TO 1994
      ENDDO
!2000 :: end of the ibr-loop
!
!     .. update SVA(N)
      IF( ( SVA( N ) < ROOTBIG ) .AND. ( SVA( N ) > ROOTSFMIN ) ) &
          THEN
         SVA( N ) = DNRM2( M, A( 1, N ), 1 )*D( N )
      ELSE
         T = 0.0D0
         AAPP = 1.0D0
         CALL DLASSQ( M, A( 1, N ), 1, T, AAPP )
         SVA( N ) = T*DSQRT( AAPP )*D( N )
      END IF
!
!     Additional steering devices
!
      IF( ( i < SWBAND ) .AND. ( ( MXAAPQ <= ROOTTOL ) .OR. &
          ( ISWROT <= N ) ) )SWBAND = i

      IF( ( i > SWBAND+1 ) .AND. ( MXAAPQ < DBLE( N )*TOL ) .AND. &
          ( DBLE( N )*MXAAPQ*MXSINJ < TOL ) ) THEN
         GO TO 1994
      END IF

!
      IF( NOTROT >= EMPTSW )GO TO 1994

      ENDDO
!     end i=1:NSWEEP loop
! #:) Reaching this point means that the procedure has completed the given
!     number of sweeps.
   INFO = NSWEEP - 1
   GO TO 1995
 1994 CONTINUE
! #:) Reaching this point means that during the i-th sweep all pivots were
!     below the given threshold, causing early exit.

   INFO = 0
! #:) INFO = 0 confirms successful iterations.
 1995 CONTINUE
!
!     Sort the vector D
!
   DO p = 1, N - 1
      q = IDAMAX( N-p+1, SVA( p ), 1 ) + p - 1
      IF( p /= q ) THEN
         TEMP1 = SVA( p )
         SVA( p ) = SVA( q )
         SVA( q ) = TEMP1
         TEMP1 = D( p )
         D( p ) = D( q )
         D( q ) = TEMP1
         CALL DSWAP( M, A( 1, p ), 1, A( 1, q ), 1 )
         IF( RSVEC )CALL DSWAP( MVL, V( 1, p ), 1, V( 1, q ), 1 )
      END IF
   ENDDO
!
   RETURN
!     ..
!     .. END OF DGSVJ1
!     ..
END
