!> \brief \b CSTEIN
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CSTEIN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cstein.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cstein.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cstein.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK,
!                          IWORK, IFAIL, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDZ, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IBLOCK( * ), IFAIL( * ), ISPLIT( * ),
!      $                   IWORK( * )
!       REAL               D( * ), E( * ), W( * ), WORK( * )
!       COMPLEX            Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSTEIN computes the eigenvectors of a real symmetric tridiagonal
!> matrix T corresponding to specified eigenvalues, using inverse
!> iteration.
!>
!> The maximum number of iterations allowed for each eigenvector is
!> specified by an internal parameter MAXITS (currently set to 5).
!>
!> Although the eigenvectors are real, they are stored in a complex
!> array, which may be passed to CUNMTR or CUPMTR for back
!> transformation to the eigenvectors of a complex Hermitian matrix
!> which was reduced to tridiagonal form.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The n diagonal elements of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          The (n-1) subdiagonal elements of the tridiagonal matrix
!>          T, stored in elements 1 to N-1.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of eigenvectors to be found.  0 <= M <= N.
!> \endverbatim
!>
!> \param[in] W
!> \verbatim
!>          W is REAL array, dimension (N)
!>          The first M elements of W contain the eigenvalues for
!>          which eigenvectors are to be computed.  The eigenvalues
!>          should be grouped by split-off block and ordered from
!>          smallest to largest within the block.  ( The output array
!>          W from SSTEBZ with ORDER = 'B' is expected here. )
!> \endverbatim
!>
!> \param[in] IBLOCK
!> \verbatim
!>          IBLOCK is INTEGER array, dimension (N)
!>          The submatrix indices associated with the corresponding
!>          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to
!>          the first submatrix from the top, =2 if W(i) belongs to
!>          the second submatrix, etc.  ( The output array IBLOCK
!>          from SSTEBZ is expected here. )
!> \endverbatim
!>
!> \param[in] ISPLIT
!> \verbatim
!>          ISPLIT is INTEGER array, dimension (N)
!>          The splitting points, at which T breaks up into submatrices.
!>          The first submatrix consists of rows/columns 1 to
!>          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1
!>          through ISPLIT( 2 ), etc.
!>          ( The output array ISPLIT from SSTEBZ is expected here. )
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ, M)
!>          The computed eigenvectors.  The eigenvector associated
!>          with the eigenvalue W(i) is stored in the i-th column of
!>          Z.  Any vector which fails to converge is set to its current
!>          iterate after MAXITS iterations.
!>          The imaginary parts of the eigenvectors are set to zero.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (5*N)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N)
!> \endverbatim
!>
!> \param[out] IFAIL
!> \verbatim
!>          IFAIL is INTEGER array, dimension (M)
!>          On normal exit, all elements of IFAIL are zero.
!>          If one or more eigenvectors fail to converge after
!>          MAXITS iterations, then their indices are stored in
!>          array IFAIL.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, then i eigenvectors failed to converge
!>               in MAXITS iterations.  Their indices are stored in
!>               array IFAIL.
!> \endverbatim
!
!> \par Internal Parameters:
!  =========================
!>
!> \verbatim
!>  MAXITS  INTEGER, default = 5
!>          The maximum number of iterations performed.
!>
!>  EXTRA   INTEGER, default = 2
!>          The number of iterations performed after norm growth
!>          criterion is satisfied, should be at least 1.
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
!> \ingroup stein
!
!  =====================================================================
   SUBROUTINE CSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK, &
                      IWORK, IFAIL, INFO )
   IMPLICIT NONE
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, LDZ, M, N
!     ..
!     .. Array Arguments ..
   INTEGER            IBLOCK( * ), IFAIL( * ), ISPLIT( * ), &
                      IWORK( * )
   REAL               D( * ), E( * ), W( * ), WORK( * )
   COMPLEX            Z( LDZ, * )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
   INTEGER            MAXITS, EXTRA
   PARAMETER          ( MAXITS = 5, EXTRA = 2 )
!     ..
!     .. Local Scalars ..
   INTEGER            B1, BLKSIZ, BN, GPIND, I, IINFO, INDRV1, &
                      INDRV2, INDRV3, INDRV4, INDRV5, ITS, J, J1, &
                      JBLK, JMAX, JR, NBLK, NRMCHK
   REAL               CTR, EPS, EPS1, NRM, ONENRM, ORTOL, PERTOL, &
                      SCL, SEP, STPCRT, TOL, XJ, XJM
!     ..
!     .. Local Arrays ..
   INTEGER            ISEED( 4 )
!     ..
!     .. External Functions ..
   INTEGER            ISAMAX
   REAL               SLAMCH, SNRM2
   EXTERNAL           ISAMAX, SLAMCH, SNRM2
!     ..
!     .. External Subroutines ..
   EXTERNAL           SCOPY, SLAGTF, SLAGTS, SLARNV, SSCAL, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   DO I = 1, M
      IFAIL( I ) = 0
   ENDDO
!
   IF( N < 0 ) THEN
      INFO = -1
   ELSE IF( M < 0 .OR. M > N ) THEN
      INFO = -4
   ELSE IF( LDZ < MAX( 1, N ) ) THEN
      INFO = -9
   ELSE
      DO J = 2, M
         IF( IBLOCK( J ) < IBLOCK( J-1 ) ) THEN
            INFO = -6
            GO TO 30
         END IF
         IF( IBLOCK( J ) == IBLOCK( J-1 ) .AND. W( J ) < W( J-1 ) ) THEN
            INFO = -5
            GO TO 30
         END IF
      ENDDO
30    CONTINUE
   END IF
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CSTEIN', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 .OR. M == 0 ) THEN
      RETURN
   ELSE IF( N == 1 ) THEN
      Z( 1, 1 ) = (1.0E+0,0.0E+0 )
      RETURN
   END IF
!
!     Get machine constants.
!
   EPS = SLAMCH( 'Precision' )
!
!     Initialize seed for random number generator SLARNV.
!
   DO I = 1, 4
      ISEED( I ) = 1
   ENDDO
!
!     Initialize pointers.
!
   INDRV1 = 0
   INDRV2 = INDRV1 + N
   INDRV3 = INDRV2 + N
   INDRV4 = INDRV3 + N
   INDRV5 = INDRV4 + N
!
!     Compute eigenvectors of matrix blocks.
!
   J1 = 1
   DO NBLK = 1, IBLOCK( M )
!
!        Find starting and ending indices of block nblk.
!
      IF( NBLK == 1 ) THEN
         B1 = 1
      ELSE
         B1 = ISPLIT( NBLK-1 ) + 1
      END IF
      BN = ISPLIT( NBLK )
      BLKSIZ = BN - B1 + 1
      IF( BLKSIZ == 1 ) GO TO 60
      GPIND = J1
!
!        Compute reorthogonalization criterion and stopping criterion.
!
      ONENRM = ABS( D( B1 ) ) + ABS( E( B1 ) )
      ONENRM = MAX( ONENRM, ABS( D( BN ) )+ABS( E( BN-1 ) ) )
      DO I = B1 + 1, BN - 1
         ONENRM = MAX( ONENRM, ABS( D( I ) )+ABS( E( I-1 ) )+ &
                  ABS( E( I ) ) )
      ENDDO
      ORTOL = 1.0E-3*ONENRM
!
      STPCRT = SQRT( 1.0E-1 / BLKSIZ )
!
!        Loop through eigenvalues of block nblk.
!
60    CONTINUE
      JBLK = 0
      DO J = J1, M
         IF( IBLOCK( J ) /= NBLK ) THEN
            J1 = J
            GO TO 180
         END IF
         JBLK = JBLK + 1
         XJ = W( J )
!
!           Skip all the work if the block size is one.
!
         IF( BLKSIZ == 1 ) THEN
            WORK( INDRV1+1 ) = 1.0E+0
            GO TO 140
         END IF
!
!           If eigenvalues j and j-1 are too close, add a relatively
!           small perturbation.
!
         IF( JBLK > 1 ) THEN
            EPS1 = ABS( EPS*XJ )
            PERTOL = 10.0E+0*EPS1
            SEP = XJ - XJM
            IF( SEP < PERTOL ) XJ = XJM + PERTOL
         END IF
!
         ITS = 0
         NRMCHK = 0
!
!           Get random starting vector.
!
         CALL SLARNV( 2, ISEED, BLKSIZ, WORK( INDRV1+1 ) )
!
!           Copy the matrix T so it won't be destroyed in factorization.
!
         WORK(INDRV4+1:INDRV4+BLKSIZ) = D(B1:B1+BLKSIZ-1)
         WORK(INDRV2+2:INDRV2+BLKSIZ) = E(B1:B1+BLKSIZ-2)
         WORK(INDRV3+1:INDRV3+BLKSIZ) = E(B1:B1+BLKSIZ-2)
!
!           Compute LU factors with partial pivoting  ( PT = LU )
!
         TOL = 0.0E+0
         CALL SLAGTF( BLKSIZ, WORK( INDRV4+1 ), XJ, WORK( INDRV2+2 ), &
                      WORK( INDRV3+1 ), TOL, WORK( INDRV5+1 ), IWORK, &
                      IINFO )
!
!           Update iteration count.
!
70       CONTINUE
         ITS = ITS + 1
         IF( ITS > MAXITS ) GO TO 120
!
!           Normalize and scale the righthand side vector Pb.
!
         JMAX = ISAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
         SCL = BLKSIZ*ONENRM*MAX( EPS, ABS( WORK( INDRV4+BLKSIZ ) ) ) / &
               ABS( WORK( INDRV1+JMAX ) )
         WORK(INDRV1+1:INDRV1+BLKSIZ) = SCL*WORK(INDRV1+1:INDRV1+BLKSIZ)
!
!           Solve the system LU = Pb.
!
         CALL SLAGTS( -1, BLKSIZ, WORK( INDRV4+1 ), WORK( INDRV2+2 ), &
                      WORK( INDRV3+1 ), WORK( INDRV5+1 ), IWORK, &
                      WORK( INDRV1+1 ), TOL, IINFO )
!
!           Reorthogonalize by modified Gram-Schmidt if eigenvalues are
!           close enough.
!
         IF( JBLK == 1 ) GO TO 110
         IF( ABS( XJ-XJM ) > ORTOL ) GPIND = J
         IF( GPIND /= J ) THEN
            DO I = GPIND, J - 1
               CTR = 0.0E+0
               DO JR = 1, BLKSIZ
                  CTR = CTR + WORK( INDRV1+JR )* REAL( Z( B1-1+JR, I ) )
               ENDDO
               DO JR = 1, BLKSIZ
                  WORK( INDRV1+JR ) = WORK( INDRV1+JR ) - CTR*REAL( Z( B1-1+JR, I ) )
               ENDDO
            ENDDO
         END IF
!
!           Check the infinity norm of the iterate.
!
  110       CONTINUE
         JMAX = ISAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
         NRM = ABS( WORK( INDRV1+JMAX ) )
!
!           Continue for additional iterations after norm reaches
!           stopping criterion.
!
         IF( NRM < STPCRT ) GO TO 70
         NRMCHK = NRMCHK + 1
         IF( NRMCHK < EXTRA+1 ) GO TO 70
!
         GO TO 130
!
!           If stopping criterion was not satisfied, update info and
!           store eigenvector number in array ifail.
!
  120       CONTINUE
         INFO = INFO + 1
         IFAIL( INFO ) = J
!
!           Accept iterate as jth eigenvector.
!
  130       CONTINUE
         SCL = 1.0E+0 / SNRM2( BLKSIZ, WORK( INDRV1+1 ), 1 )
         JMAX = ISAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
         IF( WORK( INDRV1+JMAX ) < 0.0E+0 ) &
            SCL = -SCL
         WORK(INDRV1+1:INDRV1+BLKSIZ) = SCL*WORK(INDRV1+1:INDRV1+BLKSIZ)
  140       CONTINUE
         Z(1:N,J) = (0.0E+0,0.0E+0 )
         DO I = 1, BLKSIZ
            Z( B1+I-1, J ) = CMPLX( WORK( INDRV1+I ), 0.0E+0 )
         ENDDO
!
!           Save the shift to check eigenvalue spacing at next
!           iteration.
!
         XJM = XJ
!
         ENDDO
  180 CONTINUE
      ENDDO
!
   RETURN
!
!     End of CSTEIN
!
END
