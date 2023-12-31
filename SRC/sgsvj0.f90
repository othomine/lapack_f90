!> \brief \b SGSVJ0 pre-processor for the routine sgesvj.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGSVJ0 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgsvj0.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgsvj0.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgsvj0.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGSVJ0( JOBV, M, N, A, LDA, D, SVA, MV, V, LDV, EPS,
!                          SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDV, LWORK, M, MV, N, NSWEEP
!       REAL               EPS, SFMIN, TOL
!       CHARACTER*1        JOBV
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), SVA( N ), D( N ), V( LDV, * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGSVJ0 is called from SGESVJ as a pre-processor and that is its main
!> purpose. It applies Jacobi rotations in the same way as SGESVJ does, but
!> it does not check convergence (stopping criterion). Few tuning
!> parameters (marked by [TP]) are available for the implementer.
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
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, M-by-N matrix A, such that A*diag(D) represents
!>          the input matrix.
!>          On exit,
!>          A_onexit * D_onexit represents the input matrix A*diag(D)
!>          post-multiplied by a sequence of Jacobi rotations, where the
!>          rotation threshold and the total number of sweeps are given in
!>          TOL and NSWEEP, respectively.
!>          (See the descriptions of D, TOL and NSWEEP.)
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
!>          D is REAL array, dimension (N)
!>          The array D accumulates the scaling factors from the fast scaled
!>          Jacobi rotations.
!>          On entry, A*diag(D) represents the input matrix.
!>          On exit, A_onexit*diag(D_onexit) represents the input matrix
!>          post-multiplied by a sequence of Jacobi rotations, where the
!>          rotation threshold and the total number of sweeps are given in
!>          TOL and NSWEEP, respectively.
!>          (See the descriptions of A, TOL and NSWEEP.)
!> \endverbatim
!>
!> \param[in,out] SVA
!> \verbatim
!>          SVA is REAL array, dimension (N)
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
!>                           sequence of Jacobi rotations.
!>          If JOBV = 'N',   then MV is not referenced.
!> \endverbatim
!>
!> \param[in,out] V
!> \verbatim
!>          V is REAL array, dimension (LDV,N)
!>          If JOBV = 'V' then N rows of V are post-multiplied by a
!>                           sequence of Jacobi rotations.
!>          If JOBV = 'A' then MV rows of V are post-multiplied by a
!>                           sequence of Jacobi rotations.
!>          If JOBV = 'N',   then V is not referenced.
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
!>          EPS is REAL
!>          EPS = SLAMCH('Epsilon')
!> \endverbatim
!>
!> \param[in] SFMIN
!> \verbatim
!>          SFMIN is REAL
!>          SFMIN = SLAMCH('Safe Minimum')
!> \endverbatim
!>
!> \param[in] TOL
!> \verbatim
!>          TOL is REAL
!>          TOL is the threshold for Jacobi rotations. For a pair
!>          A(:,p), A(:,q) of pivot columns, the Jacobi rotation is
!>          applied only if ABS(COS(angle(A(:,p),A(:,q)))) > TOL.
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
!>          WORK is REAL array, dimension (LWORK)
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
!> \ingroup gsvj0
!
!> \par Further Details:
!  =====================
!>
!> SGSVJ0 is used just to enable SGESVJ to call a simplified version of
!> itself to work on a submatrix of the original matrix.
!>
!> \par Contributors:
!  ==================
!>
!> Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany)
!>
!> \par Bugs, Examples and Comments:
!  =================================
!>
!> Please report all bugs and send interesting test examples and comments to
!> drmac@math.hr. Thank you.
!
!  =====================================================================
   SUBROUTINE SGSVJ0( JOBV, M, N, A, LDA, D, SVA, MV, V, LDV, EPS, &
                      SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, LDA, LDV, LWORK, M, MV, N, NSWEEP
   REAL               EPS, SFMIN, TOL
   CHARACTER*1        JOBV
!     ..
!     .. Array Arguments ..
   REAL               A( LDA, * ), SVA( N ), D( N ), V( LDV, * ), &
                      WORK( LWORK )
!     ..
!
!  =====================================================================
!
!     .. Local Parameters ..
   REAL               ZERO, HALF, ONE
   PARAMETER          ( ZERO = 0.0E0, HALF = 0.5E0, ONE = 1.0E0)
!     ..
!     .. Local Scalars ..
   REAL               AAPP, AAPP0, AAPQ, AAQQ, APOAQ, AQOAP, BIG, &
                      BIGTHETA, CS, MXAAPQ, MXSINJ, ROOTBIG, ROOTEPS, &
                      ROOTSFMIN, ROOTTOL, SMALL, SN, T, TEMP1, THETA, &
                      THSIGN
   INTEGER            BLSKIP, EMPTSW, i, ibr, IERR, igl, IJBLSK, ir1, &
                      ISWROT, jbc, jgl, KBL, LKAHEAD, MVL, NBL, &
                      NOTROT, p, PSKIPPED, q, ROWSKIP, SWBAND
   LOGICAL            APPLV, ROTOK, RSVEC
!     ..
!     .. Local Arrays ..
   REAL               FASTR( 5 )
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, MAX, FLOAT, MIN, SIGN, SQRT
!     ..
!     .. External Functions ..
   REAL               SDOT, SNRM2
   INTEGER            ISAMAX
   LOGICAL            LSAME
   EXTERNAL           ISAMAX, LSAME, SDOT, SNRM2
!     ..
!     .. External Subroutines ..
   EXTERNAL           SAXPY, SCOPY, SLASCL, SLASSQ, SROTM, SSWAP, &
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
   ELSE IF( LDA < M ) THEN
      INFO = -5
   ELSE IF( ( RSVEC.OR.APPLV ) .AND. ( MV < 0 ) ) THEN
      INFO = -8
   ELSE IF( ( RSVEC.AND.( LDV < N ) ).OR. &
            ( APPLV.AND.( LDV < MV ) ) ) THEN
      INFO = -10
   ELSE IF( TOL <= EPS ) THEN
      INFO = -13
   ELSE IF( NSWEEP < 0 ) THEN
      INFO = -14
   ELSE IF( LWORK < M ) THEN
      INFO = -16
   ELSE
      INFO = 0
   END IF
!
!     #:(
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'SGSVJ0', -INFO )
      RETURN
   END IF
!
   IF( RSVEC ) THEN
      MVL = N
   ELSE IF( APPLV ) THEN
      MVL = MV
   END IF
   RSVEC = RSVEC .OR. APPLV

   ROOTEPS = SQRT( EPS )
   ROOTSFMIN = SQRT( SFMIN )
   SMALL = SFMIN / EPS
   BIG = ONE / SFMIN
   ROOTBIG = ONE / ROOTSFMIN
   BIGTHETA = ONE / ROOTEPS
   ROOTTOL = SQRT( TOL )
!
!     .. Row-cyclic Jacobi SVD algorithm with column pivoting ..
!
   EMPTSW = ( N*( N-1 ) ) / 2
   NOTROT = 0
   FASTR( 1 ) = ZERO
!
!     .. Row-cyclic pivot strategy with de Rijk's pivoting ..
!

   SWBAND = 0
![TP] SWBAND is a tuning parameter. It is meaningful and effective
!     if SGESVJ is used as a computational routine in the preconditioned
!     Jacobi SVD algorithm SGESVJ. For sweeps i=1:SWBAND the procedure
!     ......

   KBL = MIN( 8, N )
![TP] KBL is a tuning parameter that defines the tile size in the
!     tiling of the p-q loops of pivot pairs. In general, an optimal
!     value of KBL depends on the matrix dimensions and on the
!     parameters of the computer's memory.
!
   NBL = N / KBL
   IF( ( NBL*KBL ) /= N )NBL = NBL + 1

   BLSKIP = ( KBL**2 ) + 1
![TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.

   ROWSKIP = MIN( 5, KBL )
![TP] ROWSKIP is a tuning parameter.

   LKAHEAD = 1
![TP] LKAHEAD is a tuning parameter.
   SWBAND = 0
   PSKIPPED = 0
!
   DO i = 1, NSWEEP
!     .. go go go ...
!
      MXAAPQ = ZERO
      MXSINJ = ZERO
      ISWROT = 0
!
      NOTROT = 0
      PSKIPPED = 0
!
      DO ibr = 1, NBL

         igl = ( ibr-1 )*KBL + 1
!
         DO ir1 = 0, MIN( LKAHEAD, NBL-ibr )
!
            igl = igl + ir1*KBL
!
            DO p = igl, MIN( igl+KBL-1, N-1 )

!     .. de Rijk's pivoting
               q = ISAMAX( N-p+1, SVA( p ), 1 ) + p - 1
               IF( p /= q ) THEN
                  CALL SSWAP( M, A( 1, p ), 1, A( 1, q ), 1 )
                  IF( RSVEC )CALL SSWAP( MVL, V( 1, p ), 1, &
                                         V( 1, q ), 1 )
                  TEMP1 = SVA( p )
                  SVA( p ) = SVA( q )
                  SVA( q ) = TEMP1
                  TEMP1 = D( p )
                  D( p ) = D( q )
                  D( q ) = TEMP1
               END IF
!
               IF( ir1 == 0 ) THEN
!
!        Column norms are periodically updated by explicit
!        norm computation.
!        Caveat:
!        Some BLAS implementations compute SNRM2(M,A(1,p),1)
!        as SQRT(SDOT(M,A(1,p),1,A(1,p),1)), which may result in
!        overflow for ||A(:,p)||_2 > SQRT(overflow_threshold), and
!        underflow for ||A(:,p)||_2 < SQRT(underflow_threshold).
!        Hence, SNRM2 cannot be trusted, not even in the case when
!        the true norm is far from the under(over)flow boundaries.
!        If properly implemented SNRM2 is available, the IF-THEN-ELSE
!        below should read "AAPP = SNRM2( M, A(1,p), 1 ) * D(p)".
!
                  IF( ( SVA( p ) < ROOTBIG ) .AND. &
                      ( SVA( p ) > ROOTSFMIN ) ) THEN
                     SVA( p ) = SNRM2( M, A( 1, p ), 1 )*D( p )
                  ELSE
                     TEMP1 = ZERO
                     AAPP = ONE
                     CALL SLASSQ( M, A( 1, p ), 1, TEMP1, AAPP )
                     SVA( p ) = TEMP1*SQRT( AAPP )*D( p )
                  END IF
                  AAPP = SVA( p )
               ELSE
                  AAPP = SVA( p )
               END IF

!
               IF( AAPP > ZERO ) THEN
!
                  PSKIPPED = 0
!
                  DO q = p + 1, MIN( igl+KBL-1, N )
!
                     AAQQ = SVA( q )

                     IF( AAQQ > ZERO ) THEN
!
                        AAPP0 = AAPP
                        IF( AAQQ >= ONE ) THEN
                           ROTOK = ( SMALL*AAPP ) <= AAQQ
                           IF( AAPP < ( BIG / AAQQ ) ) THEN
                              AAPQ = ( SDOT( M, A( 1, p ), 1, A( 1, &
                                     q ), 1 )*D( p )*D( q ) / AAQQ ) &
                                     / AAPP
                           ELSE
                              CALL SCOPY( M, A( 1, p ), 1, WORK, 1 )
                              CALL SLASCL( 'G', 0, 0, AAPP, D( p ), &
                                           M, 1, WORK, LDA, IERR )
                              AAPQ = SDOT( M, WORK, 1, A( 1, q ), &
                                     1 )*D( q ) / AAQQ
                           END IF
                        ELSE
                           ROTOK = AAPP <= ( AAQQ / SMALL )
                           IF( AAPP > ( SMALL / AAQQ ) ) THEN
                              AAPQ = ( SDOT( M, A( 1, p ), 1, A( 1, &
                                     q ), 1 )*D( p )*D( q ) / AAQQ ) &
                                     / AAPP
                           ELSE
                              CALL SCOPY( M, A( 1, q ), 1, WORK, 1 )
                              CALL SLASCL( 'G', 0, 0, AAQQ, D( q ), &
                                           M, 1, WORK, LDA, IERR )
                              AAPQ = SDOT( M, WORK, 1, A( 1, p ), &
                                     1 )*D( p ) / AAPP
                           END IF
                        END IF
!
                        MXAAPQ = MAX( MXAAPQ, ABS( AAPQ ) )
!
!        TO rotate or NOT to rotate, THAT is the question ...
!
                        IF( ABS( AAPQ ) > TOL ) THEN
!
!           .. rotate
!           ROTATED = ROTATED + ONE
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
                              THETA = -HALF*ABS( AQOAP-APOAQ ) / AAPQ
!
                              IF( ABS( THETA ) > BIGTHETA ) THEN
!
                                 T = HALF / THETA
                                 FASTR( 3 ) = T*D( p ) / D( q )
                                 FASTR( 4 ) = -T*D( q ) / D( p )
                                 CALL SROTM( M, A( 1, p ), 1, &
                                             A( 1, q ), 1, FASTR )
                                 IF( RSVEC )CALL SROTM( MVL, &
                                                 V( 1, p ), 1, &
                                                 V( 1, q ), 1, &
                                                 FASTR )
                                 SVA( q ) = AAQQ*SQRT( MAX( ZERO, &
                                            ONE+T*APOAQ*AAPQ ) )
                                 AAPP = AAPP*SQRT( MAX( ZERO, &
                                            ONE-T*AQOAP*AAPQ ) )
                                 MXSINJ = MAX( MXSINJ, ABS( T ) )
!
                              ELSE
!
!                 .. choose correct signum for THETA and rotate
!
                                 THSIGN = -SIGN( ONE, AAPQ )
                                 T = ONE / ( THETA+THSIGN* &
                                     SQRT( ONE+THETA*THETA ) )
                                 CS = SQRT( ONE / ( ONE+T*T ) )
                                 SN = T*CS
!
                                 MXSINJ = MAX( MXSINJ, ABS( SN ) )
                                 SVA( q ) = AAQQ*SQRT( MAX( ZERO, &
                                            ONE+T*APOAQ*AAPQ ) )
                                 AAPP = AAPP*SQRT( MAX( ZERO, &
                                        ONE-T*AQOAP*AAPQ ) )
!
                                 APOAQ = D( p ) / D( q )
                                 AQOAP = D( q ) / D( p )
                                 IF( D( p ) >= ONE ) THEN
                                    IF( D( q ) >= ONE ) THEN
                                       FASTR( 3 ) = T*APOAQ
                                       FASTR( 4 ) = -T*AQOAP
                                       D( p ) = D( p )*CS
                                       D( q ) = D( q )*CS
                                       CALL SROTM( M, A( 1, p ), 1, &
                                                   A( 1, q ), 1, &
                                                   FASTR )
                                       IF( RSVEC )CALL SROTM( MVL, &
                                           V( 1, p ), 1, V( 1, q ), &
                                           1, FASTR )
                                    ELSE
                                       CALL SAXPY( M, -T*AQOAP, &
                                                   A( 1, q ), 1, &
                                                   A( 1, p ), 1 )
                                       CALL SAXPY( M, CS*SN*APOAQ, &
                                                   A( 1, p ), 1, &
                                                   A( 1, q ), 1 )
                                       D( p ) = D( p )*CS
                                       D( q ) = D( q ) / CS
                                       IF( RSVEC ) THEN
                                          CALL SAXPY( MVL, -T*AQOAP, &
                                                      V( 1, q ), 1, &
                                                      V( 1, p ), 1 )
                                          CALL SAXPY( MVL, &
                                                      CS*SN*APOAQ, &
                                                      V( 1, p ), 1, &
                                                      V( 1, q ), 1 )
                                       END IF
                                    END IF
                                 ELSE
                                    IF( D( q ) >= ONE ) THEN
                                       CALL SAXPY( M, T*APOAQ, &
                                                   A( 1, p ), 1, &
                                                   A( 1, q ), 1 )
                                       CALL SAXPY( M, -CS*SN*AQOAP, &
                                                   A( 1, q ), 1, &
                                                   A( 1, p ), 1 )
                                       D( p ) = D( p ) / CS
                                       D( q ) = D( q )*CS
                                       IF( RSVEC ) THEN
                                          CALL SAXPY( MVL, T*APOAQ, &
                                                      V( 1, p ), 1, &
                                                      V( 1, q ), 1 )
                                          CALL SAXPY( MVL, &
                                                      -CS*SN*AQOAP, &
                                                      V( 1, q ), 1, &
                                                      V( 1, p ), 1 )
                                       END IF
                                    ELSE
                                       IF( D( p ) >= D( q ) ) THEN
                                          CALL SAXPY( M, -T*AQOAP, &
                                                      A( 1, q ), 1, &
                                                      A( 1, p ), 1 )
                                          CALL SAXPY( M, CS*SN*APOAQ, &
                                                      A( 1, p ), 1, &
                                                      A( 1, q ), 1 )
                                          D( p ) = D( p )*CS
                                          D( q ) = D( q ) / CS
                                          IF( RSVEC ) THEN
                                             CALL SAXPY( MVL, &
                                                  -T*AQOAP, &
                                                  V( 1, q ), 1, &
                                                  V( 1, p ), 1 )
                                             CALL SAXPY( MVL, &
                                                  CS*SN*APOAQ, &
                                                  V( 1, p ), 1, &
                                                  V( 1, q ), 1 )
                                          END IF
                                       ELSE
                                          CALL SAXPY( M, T*APOAQ, &
                                                      A( 1, p ), 1, &
                                                      A( 1, q ), 1 )
                                          CALL SAXPY( M, &
                                                      -CS*SN*AQOAP, &
                                                      A( 1, q ), 1, &
                                                      A( 1, p ), 1 )
                                          D( p ) = D( p ) / CS
                                          D( q ) = D( q )*CS
                                          IF( RSVEC ) THEN
                                             CALL SAXPY( MVL, &
                                                  T*APOAQ, V( 1, p ), &
                                                  1, V( 1, q ), 1 )
                                             CALL SAXPY( MVL, &
                                                  -CS*SN*AQOAP, &
                                                  V( 1, q ), 1, &
                                                  V( 1, p ), 1 )
                                          END IF
                                       END IF
                                    END IF
                                 END IF
                              END IF
!
                           ELSE
!              .. have to use modified Gram-Schmidt like transformation
                              CALL SCOPY( M, A( 1, p ), 1, WORK, 1 )
                              CALL SLASCL( 'G', 0, 0, AAPP, ONE, M, &
                                           1, WORK, LDA, IERR )
                              CALL SLASCL( 'G', 0, 0, AAQQ, ONE, M, &
                                           1, A( 1, q ), LDA, IERR )
                              TEMP1 = -AAPQ*D( p ) / D( q )
                              CALL SAXPY( M, TEMP1, WORK, 1, &
                                          A( 1, q ), 1 )
                              CALL SLASCL( 'G', 0, 0, ONE, AAQQ, M, &
                                           1, A( 1, q ), LDA, IERR )
                              SVA( q ) = AAQQ*SQRT( MAX( ZERO, &
                                         ONE-AAPQ*AAPQ ) )
                              MXSINJ = MAX( MXSINJ, SFMIN )
                           END IF
!           END IF ROTOK THEN ... ELSE
!
!           In the case of cancellation in updating SVA(q), SVA(p)
!           recompute SVA(q), SVA(p).
                           IF( ( SVA( q ) / AAQQ )**2 <= ROOTEPS ) &
                               THEN
                              IF( ( AAQQ < ROOTBIG ) .AND. &
                                  ( AAQQ > ROOTSFMIN ) ) THEN
                                 SVA( q ) = SNRM2( M, A( 1, q ), 1 )* &
                                            D( q )
                              ELSE
                                 T = ZERO
                                 AAQQ = ONE
                                 CALL SLASSQ( M, A( 1, q ), 1, T, &
                                              AAQQ )
                                 SVA( q ) = T*SQRT( AAQQ )*D( q )
                              END IF
                           END IF
                           IF( ( AAPP / AAPP0 ) <= ROOTEPS ) THEN
                              IF( ( AAPP < ROOTBIG ) .AND. &
                                  ( AAPP > ROOTSFMIN ) ) THEN
                                 AAPP = SNRM2( M, A( 1, p ), 1 )* &
                                        D( p )
                              ELSE
                                 T = ZERO
                                 AAPP = ONE
                                 CALL SLASSQ( M, A( 1, p ), 1, T, &
                                              AAPP )
                                 AAPP = T*SQRT( AAPP )*D( p )
                              END IF
                              SVA( p ) = AAPP
                           END IF
!
                        ELSE
!        A(:,p) and A(:,q) already numerically orthogonal
                           IF( ir1 == 0 )NOTROT = NOTROT + 1
                           PSKIPPED = PSKIPPED + 1
                        END IF
                     ELSE
!        A(:,q) is zero column
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

                  SVA( p ) = AAPP

               ELSE
                  SVA( p ) = AAPP
                  IF( ( ir1 == 0 ) .AND. ( AAPP == ZERO ) ) &
                      NOTROT = NOTROT + MIN( igl+KBL-1, N ) - p
               END IF
!
               ENDDO
!     end of the p-loop
!     end of doing the block ( ibr, ibr )
            ENDDO
!     end of ir1-loop
!
!........................................................
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
!
               IF( AAPP > ZERO ) THEN
!
                  PSKIPPED = 0
!
                  DO q = jgl, MIN( jgl+KBL-1, N )
!
                     AAQQ = SVA( q )
!
                     IF( AAQQ > ZERO ) THEN
                        AAPP0 = AAPP
!
!     .. M x 2 Jacobi SVD ..
!
!        .. Safe Gram matrix computation ..
!
                        IF( AAQQ >= ONE ) THEN
                           IF( AAPP >= AAQQ ) THEN
                              ROTOK = ( SMALL*AAPP ) <= AAQQ
                           ELSE
                              ROTOK = ( SMALL*AAQQ ) <= AAPP
                           END IF
                           IF( AAPP < ( BIG / AAQQ ) ) THEN
                              AAPQ = ( SDOT( M, A( 1, p ), 1, A( 1, &
                                     q ), 1 )*D( p )*D( q ) / AAQQ ) &
                                     / AAPP
                           ELSE
                              CALL SCOPY( M, A( 1, p ), 1, WORK, 1 )
                              CALL SLASCL( 'G', 0, 0, AAPP, D( p ), &
                                           M, 1, WORK, LDA, IERR )
                              AAPQ = SDOT( M, WORK, 1, A( 1, q ), &
                                     1 )*D( q ) / AAQQ
                           END IF
                        ELSE
                           IF( AAPP >= AAQQ ) THEN
                              ROTOK = AAPP <= ( AAQQ / SMALL )
                           ELSE
                              ROTOK = AAQQ <= ( AAPP / SMALL )
                           END IF
                           IF( AAPP > ( SMALL / AAQQ ) ) THEN
                              AAPQ = ( SDOT( M, A( 1, p ), 1, A( 1, &
                                     q ), 1 )*D( p )*D( q ) / AAQQ ) &
                                     / AAPP
                           ELSE
                              CALL SCOPY( M, A( 1, q ), 1, WORK, 1 )
                              CALL SLASCL( 'G', 0, 0, AAQQ, D( q ), &
                                           M, 1, WORK, LDA, IERR )
                              AAPQ = SDOT( M, WORK, 1, A( 1, p ), &
                                     1 )*D( p ) / AAPP
                           END IF
                        END IF
!
                        MXAAPQ = MAX( MXAAPQ, ABS( AAPQ ) )
!
!        TO rotate or NOT to rotate, THAT is the question ...
!
                        IF( ABS( AAPQ ) > TOL ) THEN
                           NOTROT = 0
!           ROTATED  = ROTATED + 1
                           PSKIPPED = 0
                           ISWROT = ISWROT + 1
!
                           IF( ROTOK ) THEN
!
                              AQOAP = AAQQ / AAPP
                              APOAQ = AAPP / AAQQ
                              THETA = -HALF*ABS( AQOAP-APOAQ ) / AAPQ
                              IF( AAQQ > AAPP0 )THETA = -THETA
!
                              IF( ABS( THETA ) > BIGTHETA ) THEN
                                 T = HALF / THETA
                                 FASTR( 3 ) = T*D( p ) / D( q )
                                 FASTR( 4 ) = -T*D( q ) / D( p )
                                 CALL SROTM( M, A( 1, p ), 1, &
                                             A( 1, q ), 1, FASTR )
                                 IF( RSVEC )CALL SROTM( MVL, &
                                                 V( 1, p ), 1, &
                                                 V( 1, q ), 1, &
                                                 FASTR )
                                 SVA( q ) = AAQQ*SQRT( MAX( ZERO, &
                                            ONE+T*APOAQ*AAPQ ) )
                                 AAPP = AAPP*SQRT( MAX( ZERO, &
                                        ONE-T*AQOAP*AAPQ ) )
                                 MXSINJ = MAX( MXSINJ, ABS( T ) )
                              ELSE
!
!                 .. choose correct signum for THETA and rotate
!
                                 THSIGN = -SIGN( ONE, AAPQ )
                                 IF( AAQQ > AAPP0 )THSIGN = -THSIGN
                                 T = ONE / ( THETA+THSIGN* &
                                     SQRT( ONE+THETA*THETA ) )
                                 CS = SQRT( ONE / ( ONE+T*T ) )
                                 SN = T*CS
                                 MXSINJ = MAX( MXSINJ, ABS( SN ) )
                                 SVA( q ) = AAQQ*SQRT( MAX( ZERO, &
                                            ONE+T*APOAQ*AAPQ ) )
                                 AAPP = AAPP*SQRT( MAX( ZERO, &
                                            ONE-T*AQOAP*AAPQ ) )
!
                                 APOAQ = D( p ) / D( q )
                                 AQOAP = D( q ) / D( p )
                                 IF( D( p ) >= ONE ) THEN
!
                                    IF( D( q ) >= ONE ) THEN
                                       FASTR( 3 ) = T*APOAQ
                                       FASTR( 4 ) = -T*AQOAP
                                       D( p ) = D( p )*CS
                                       D( q ) = D( q )*CS
                                       CALL SROTM( M, A( 1, p ), 1, &
                                                   A( 1, q ), 1, &
                                                   FASTR )
                                       IF( RSVEC )CALL SROTM( MVL, &
                                           V( 1, p ), 1, V( 1, q ), &
                                           1, FASTR )
                                    ELSE
                                       CALL SAXPY( M, -T*AQOAP, &
                                                   A( 1, q ), 1, &
                                                   A( 1, p ), 1 )
                                       CALL SAXPY( M, CS*SN*APOAQ, &
                                                   A( 1, p ), 1, &
                                                   A( 1, q ), 1 )
                                       IF( RSVEC ) THEN
                                          CALL SAXPY( MVL, -T*AQOAP, &
                                                      V( 1, q ), 1, &
                                                      V( 1, p ), 1 )
                                          CALL SAXPY( MVL, &
                                                      CS*SN*APOAQ, &
                                                      V( 1, p ), 1, &
                                                      V( 1, q ), 1 )
                                       END IF
                                       D( p ) = D( p )*CS
                                       D( q ) = D( q ) / CS
                                    END IF
                                 ELSE
                                    IF( D( q ) >= ONE ) THEN
                                       CALL SAXPY( M, T*APOAQ, &
                                                   A( 1, p ), 1, &
                                                   A( 1, q ), 1 )
                                       CALL SAXPY( M, -CS*SN*AQOAP, &
                                                   A( 1, q ), 1, &
                                                   A( 1, p ), 1 )
                                       IF( RSVEC ) THEN
                                          CALL SAXPY( MVL, T*APOAQ, &
                                                      V( 1, p ), 1, &
                                                      V( 1, q ), 1 )
                                          CALL SAXPY( MVL, &
                                                      -CS*SN*AQOAP, &
                                                      V( 1, q ), 1, &
                                                      V( 1, p ), 1 )
                                       END IF
                                       D( p ) = D( p ) / CS
                                       D( q ) = D( q )*CS
                                    ELSE
                                       IF( D( p ) >= D( q ) ) THEN
                                          CALL SAXPY( M, -T*AQOAP, &
                                                      A( 1, q ), 1, &
                                                      A( 1, p ), 1 )
                                          CALL SAXPY( M, CS*SN*APOAQ, &
                                                      A( 1, p ), 1, &
                                                      A( 1, q ), 1 )
                                          D( p ) = D( p )*CS
                                          D( q ) = D( q ) / CS
                                          IF( RSVEC ) THEN
                                             CALL SAXPY( MVL, &
                                                  -T*AQOAP, &
                                                  V( 1, q ), 1, &
                                                  V( 1, p ), 1 )
                                             CALL SAXPY( MVL, &
                                                  CS*SN*APOAQ, &
                                                  V( 1, p ), 1, &
                                                  V( 1, q ), 1 )
                                          END IF
                                       ELSE
                                          CALL SAXPY( M, T*APOAQ, &
                                                      A( 1, p ), 1, &
                                                      A( 1, q ), 1 )
                                          CALL SAXPY( M, &
                                                      -CS*SN*AQOAP, &
                                                      A( 1, q ), 1, &
                                                      A( 1, p ), 1 )
                                          D( p ) = D( p ) / CS
                                          D( q ) = D( q )*CS
                                          IF( RSVEC ) THEN
                                             CALL SAXPY( MVL, &
                                                  T*APOAQ, V( 1, p ), &
                                                  1, V( 1, q ), 1 )
                                             CALL SAXPY( MVL, &
                                                  -CS*SN*AQOAP, &
                                                  V( 1, q ), 1, &
                                                  V( 1, p ), 1 )
                                          END IF
                                       END IF
                                    END IF
                                 END IF
                              END IF
!
                           ELSE
                              IF( AAPP > AAQQ ) THEN
                                 CALL SCOPY( M, A( 1, p ), 1, WORK, &
                                             1 )
                                 CALL SLASCL( 'G', 0, 0, AAPP, ONE, &
                                              M, 1, WORK, LDA, IERR )
                                 CALL SLASCL( 'G', 0, 0, AAQQ, ONE, &
                                              M, 1, A( 1, q ), LDA, &
                                              IERR )
                                 TEMP1 = -AAPQ*D( p ) / D( q )
                                 CALL SAXPY( M, TEMP1, WORK, 1, &
                                             A( 1, q ), 1 )
                                 CALL SLASCL( 'G', 0, 0, ONE, AAQQ, &
                                              M, 1, A( 1, q ), LDA, &
                                              IERR )
                                 SVA( q ) = AAQQ*SQRT( MAX( ZERO, &
                                            ONE-AAPQ*AAPQ ) )
                                 MXSINJ = MAX( MXSINJ, SFMIN )
                              ELSE
                                 CALL SCOPY( M, A( 1, q ), 1, WORK, &
                                             1 )
                                 CALL SLASCL( 'G', 0, 0, AAQQ, ONE, &
                                              M, 1, WORK, LDA, IERR )
                                 CALL SLASCL( 'G', 0, 0, AAPP, ONE, &
                                              M, 1, A( 1, p ), LDA, &
                                              IERR )
                                 TEMP1 = -AAPQ*D( q ) / D( p )
                                 CALL SAXPY( M, TEMP1, WORK, 1, &
                                             A( 1, p ), 1 )
                                 CALL SLASCL( 'G', 0, 0, ONE, AAPP, &
                                              M, 1, A( 1, p ), LDA, &
                                              IERR )
                                 SVA( p ) = AAPP*SQRT( MAX( ZERO, &
                                            ONE-AAPQ*AAPQ ) )
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
                                 SVA( q ) = SNRM2( M, A( 1, q ), 1 )* &
                                            D( q )
                              ELSE
                                 T = ZERO
                                 AAQQ = ONE
                                 CALL SLASSQ( M, A( 1, q ), 1, T, &
                                              AAQQ )
                                 SVA( q ) = T*SQRT( AAQQ )*D( q )
                              END IF
                           END IF
                           IF( ( AAPP / AAPP0 )**2 <= ROOTEPS ) THEN
                              IF( ( AAPP < ROOTBIG ) .AND. &
                                  ( AAPP > ROOTSFMIN ) ) THEN
                                 AAPP = SNRM2( M, A( 1, p ), 1 )* &
                                        D( p )
                              ELSE
                                 T = ZERO
                                 AAPP = ONE
                                 CALL SLASSQ( M, A( 1, p ), 1, T, &
                                              AAPP )
                                 AAPP = T*SQRT( AAPP )*D( p )
                              END IF
                              SVA( p ) = AAPP
                           END IF
!              end of OK rotation
                        ELSE
                           NOTROT = NOTROT + 1
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
                        GO TO 2203
                     END IF
!
                     ENDDO
!        end of the q-loop
 2203                CONTINUE
!
                  SVA( p ) = AAPP
!
               ELSE
                  IF( AAPP == ZERO )NOTROT = NOTROT + &
                      MIN( jgl+KBL-1, N ) - jgl + 1
                  IF( AAPP < ZERO )NOTROT = 0
               END IF

               ENDDO
!     end of the p-loop
            ENDDO
!     end of the jbc-loop
 2011       CONTINUE
!2011 bailed out of the jbc-loop
         DO p = igl, MIN( igl+KBL-1, N )
            SVA( p ) = ABS( SVA( p ) )
            ENDDO
!
         ENDDO
!2000 :: end of the ibr-loop
!
!     .. update SVA(N)
      IF( ( SVA( N ) < ROOTBIG ) .AND. ( SVA( N ) > ROOTSFMIN ) ) &
          THEN
         SVA( N ) = SNRM2( M, A( 1, N ), 1 )*D( N )
      ELSE
         T = ZERO
         AAPP = ONE
         CALL SLASSQ( M, A( 1, N ), 1, T, AAPP )
         SVA( N ) = T*SQRT( AAPP )*D( N )
      END IF
!
!     Additional steering devices
!
      IF( ( i < SWBAND ) .AND. ( ( MXAAPQ <= ROOTTOL ) .OR. &
          ( ISWROT <= N ) ) )SWBAND = i
!
      IF( ( i > SWBAND+1 ) .AND. ( MXAAPQ < FLOAT( N )*TOL ) .AND. &
          ( FLOAT( N )*MXAAPQ*MXSINJ < TOL ) ) THEN
         GO TO 1994
      END IF
!
      IF( NOTROT >= EMPTSW )GO TO 1994

      ENDDO
!     end i=1:NSWEEP loop
! #:) Reaching this point means that the procedure has completed the given
!     number of iterations.
   INFO = NSWEEP - 1
   GO TO 1995
 1994 CONTINUE
! #:) Reaching this point means that during the i-th sweep all pivots were
!     below the given tolerance, causing early exit.
!
   INFO = 0
! #:) INFO = 0 confirms successful iterations.
 1995 CONTINUE
!
!     Sort the vector D.
   DO p = 1, N - 1
      q = ISAMAX( N-p+1, SVA( p ), 1 ) + p - 1
      IF( p /= q ) THEN
         TEMP1 = SVA( p )
         SVA( p ) = SVA( q )
         SVA( q ) = TEMP1
         TEMP1 = D( p )
         D( p ) = D( q )
         D( q ) = TEMP1
         CALL SSWAP( M, A( 1, p ), 1, A( 1, q ), 1 )
         IF( RSVEC )CALL SSWAP( MVL, V( 1, p ), 1, V( 1, q ), 1 )
      END IF
      ENDDO
!
   RETURN
!     ..
!     .. END OF SGSVJ0
!     ..
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

