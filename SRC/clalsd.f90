!> \brief \b CLALSD uses the singular value decomposition of A to solve the least squares problem.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLALSD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clalsd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clalsd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clalsd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND,
!                          RANK, WORK, RWORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDB, N, NRHS, RANK, SMLSIZ
!       REAL               RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               D( * ), E( * ), RWORK( * )
!       COMPLEX            B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLALSD uses the singular value decomposition of A to solve the least
!> squares problem of finding X to minimize the Euclidean norm of each
!> column of A*X-B, where A is N-by-N upper bidiagonal, and X and B
!> are N-by-NRHS. The solution X overwrites B.
!>
!> The singular values of A smaller than RCOND times the largest
!> singular value are treated as zero in solving the least squares
!> problem; in this case a minimum norm solution is returned.
!> The actual singular values are returned in D in ascending order.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>         = 'U': D and E define an upper bidiagonal matrix.
!>         = 'L': D and E define a  lower bidiagonal matrix.
!> \endverbatim
!>
!> \param[in] SMLSIZ
!> \verbatim
!>          SMLSIZ is INTEGER
!>         The maximum size of the subproblems at the bottom of the
!>         computation tree.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The dimension of the  bidiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>         The number of columns of B. NRHS must be at least 1.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>         On entry D contains the main diagonal of the bidiagonal
!>         matrix. On exit, if INFO = 0, D contains its singular values.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>         Contains the super-diagonal entries of the bidiagonal matrix.
!>         On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,NRHS)
!>         On input, B contains the right hand sides of the least
!>         squares problem. On output, B contains the solution X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>         The leading dimension of B in the calling subprogram.
!>         LDB must be at least max(1,N).
!> \endverbatim
!>
!> \param[in] RCOND
!> \verbatim
!>          RCOND is REAL
!>         The singular values of A less than or equal to RCOND times
!>         the largest singular value are treated as zero in solving
!>         the least squares problem. If RCOND is negative,
!>         machine precision is used instead.
!>         For example, if diag(S)*X=B were the least squares problem,
!>         where diag(S) is a diagonal matrix of singular values, the
!>         solution would be X(i) = B(i) / S(i) if S(i) is greater than
!>         RCOND*max(S), and X(i) = 0 if S(i) is less than or equal to
!>         RCOND*max(S).
!> \endverbatim
!>
!> \param[out] RANK
!> \verbatim
!>          RANK is INTEGER
!>         The number of singular values of A greater than RCOND times
!>         the largest singular value.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (N * NRHS).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension at least
!>         (9*N + 2*N*SMLSIZ + 8*N*NLVL + 3*SMLSIZ*NRHS +
!>         MAX( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS ),
!>         where
!>         NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (3*N*NLVL + 11*N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>         = 0:  successful exit.
!>         < 0:  if INFO = -i, the i-th argument had an illegal value.
!>         > 0:  The algorithm failed to compute a singular value while
!>               working on the submatrix lying in rows and columns
!>               INFO/(N+1) through MOD(INFO,N+1).
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
!> \ingroup lalsd
!
!> \par Contributors:
!  ==================
!>
!>     Ming Gu and Ren-Cang Li, Computer Science Division, University of
!>       California at Berkeley, USA \n
!>     Osni Marques, LBNL/NERSC, USA \n
!
!  =====================================================================
   SUBROUTINE CLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND, &
                      RANK, WORK, RWORK, IWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INFO, LDB, N, NRHS, RANK, SMLSIZ
   REAL               RCOND
!     ..
!     .. Array Arguments ..
   INTEGER            IWORK( * )
   REAL               D( * ), E( * ), RWORK( * )
   COMPLEX            B( LDB, * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            BX, BXST, C, DIFL, DIFR, GIVCOL, GIVNUM, &
                      GIVPTR, I, ICMPQ1, ICMPQ2, IRWB, IRWIB, IRWRB, &
                      IRWU, IRWVT, IRWWRK, IWK, J, JCOL, JIMAG, &
                      JREAL, JROW, K, NLVL, NM1, NRWORK, NSIZE, NSUB, &
                      PERM, POLES, S, SIZEI, SMLSZP, SQRE, ST, ST1, &
                      U, VT, Z
   REAL               CS, EPS, ORGNRM, R, RCND, SN, TOL
!     ..
!     .. External Functions ..
   INTEGER            ISAMAX
   REAL               SLAMCH, SLANST
   EXTERNAL           ISAMAX, SLAMCH, SLANST
!     ..
!     .. External Subroutines ..
   EXTERNAL           CCOPY, CLACPY, CLALSA, CLASCL, CLASET, CSROT, &
                      SGEMM, SLARTG, SLASCL, SLASDA, SLASDQ, SLASET, &
                      SLASRT, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
!
   IF( N < 0 ) THEN
      INFO = -3
   ELSE IF( NRHS < 1 ) THEN
      INFO = -4
   ELSE IF( ( LDB < 1 ) .OR. ( LDB < N ) ) THEN
      INFO = -8
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CLALSD', -INFO )
      RETURN
   END IF
!
   EPS = SLAMCH( 'Epsilon' )
!
!     Set up the tolerance.
!
   IF( ( RCOND <= 0.0E+0 ) .OR. ( RCOND >= 1.0E+0 ) ) THEN
      RCND = EPS
   ELSE
      RCND = RCOND
   END IF
!
   RANK = 0
!
!     Quick return if possible.
!
   IF( N == 0 ) THEN
      RETURN
   ELSE IF( N == 1 ) THEN
      IF( D( 1 ) == 0.0E+0 ) THEN
         B(1,1:NRHS) = (0.0E+0,0.0E+0)
      ELSE
         RANK = 1
         CALL CLASCL( 'G', 0, 0, D( 1 ), 1.0E+0, 1, NRHS, B, LDB, INFO )
         D( 1 ) = ABS( D( 1 ) )
      END IF
      RETURN
   END IF
!
!     Rotate the matrix if it is lower bidiagonal.
!
   IF( UPLO == 'L' ) THEN
      DO I = 1, N - 1
         CALL SLARTG( D( I ), E( I ), CS, SN, R )
         D( I ) = R
         E( I ) = SN*D( I+1 )
         D( I+1 ) = CS*D( I+1 )
         IF( NRHS == 1 ) THEN
            CALL CSROT( 1, B( I, 1 ), 1, B( I+1, 1 ), 1, CS, SN )
         ELSE
            RWORK( I*2-1 ) = CS
            RWORK( I*2 ) = SN
         END IF
      ENDDO
      IF( NRHS > 1 ) THEN
         DO I = 1, NRHS
            DO J = 1, N - 1
               CS = RWORK( J*2-1 )
               SN = RWORK( J*2 )
               CALL CSROT( 1, B( J, I ), 1, B( J+1, I ), 1, CS, SN )
            ENDDO
         ENDDO
      END IF
   END IF
!
!     Scale.
!
   NM1 = N - 1
   ORGNRM = SLANST( 'M', N, D, E )
   IF( ORGNRM == 0.0E+0 ) THEN
      B(1:N,1:NRHS) = (0.0E+0,0.0E+0)
      RETURN
   END IF
!
   CALL SLASCL( 'G', 0, 0, ORGNRM, 1.0E+0, N, 1, D, N, INFO )
   CALL SLASCL( 'G', 0, 0, ORGNRM, 1.0E+0, NM1, 1, E, NM1, INFO )
!
!     If N is smaller than the minimum divide size SMLSIZ, then solve
!     the problem with another solver.
!
   IF( N <= SMLSIZ ) THEN
      IRWU = 1
      IRWVT = IRWU + N*N
      IRWWRK = IRWVT + N*N
      IRWRB = IRWWRK
      IRWIB = IRWRB + N*NRHS
      IRWB = IRWIB + N*NRHS
      CALL SLASET( 'A', N, N, 0.0E+0, 1.0E+0, RWORK( IRWU ), N )
      CALL SLASET( 'A', N, N, 0.0E+0, 1.0E+0, RWORK( IRWVT ), N )
      CALL SLASDQ( 'U', 0, N, N, N, 0, D, E, RWORK( IRWVT ), N, &
                   RWORK( IRWU ), N, RWORK( IRWWRK ), 1, &
                   RWORK( IRWWRK ), INFO )
      IF( INFO /= 0 ) RETURN
!
!        In the real version, B is passed to SLASDQ and multiplied
!        internally by Q**H. Here B is complex and that product is
!        computed below in two steps (real and imaginary parts).
!
      J = IRWB - 1
      DO JCOL = 1, NRHS
         DO JROW = 1, N
            J = J + 1
            RWORK( J ) = REAL( B( JROW, JCOL ) )
         ENDDO
      ENDDO
      CALL SGEMM( 'T', 'N', N, NRHS, N, 1.0E+0, RWORK( IRWU ), N, &
                  RWORK( IRWB ), N, 0.0E+0, RWORK( IRWRB ), N )
      J = IRWB - 1
      DO JCOL = 1, NRHS
         DO JROW = 1, N
            J = J + 1
            RWORK( J ) = AIMAG( B( JROW, JCOL ) )
         ENDDO
      ENDDO
      CALL SGEMM( 'T', 'N', N, NRHS, N, 1.0E+0, RWORK( IRWU ), N, &
                  RWORK( IRWB ), N, 0.0E+0, RWORK( IRWIB ), N )
      JREAL = IRWRB - 1
      JIMAG = IRWIB - 1
      DO JCOL = 1, NRHS
         DO JROW = 1, N
            JREAL = JREAL + 1
            JIMAG = JIMAG + 1
            B( JROW, JCOL ) = CMPLX( RWORK( JREAL ), RWORK( JIMAG ) )
         ENDDO
      ENDDO
!
      TOL = RCND*ABS( D( ISAMAX( N, D, 1 ) ) )
      DO I = 1, N
         IF( D( I ) <= TOL ) THEN
            B(1,1:NRHS) = (0.0E+0,0.0E+0)
         ELSE
            CALL CLASCL( 'G', 0, 0, D( I ), 1.0E+0, 1, NRHS, B( I, 1 ), LDB, INFO )
            RANK = RANK + 1
         END IF
      ENDDO
!
!        Since B is complex, the following call to SGEMM is performed
!        in two steps (real and imaginary parts). That is for V * B
!        (in the real version of the code V**H is stored in WORK).
!
!        CALL SGEMM( 'T', 'N', N, NRHS, N, 1.0E+0, WORK, N, B, LDB, 0.0E+0,
!    $               WORK( NWORK ), N )
!
      J = IRWB - 1
      DO JCOL = 1, NRHS
         DO JROW = 1, N
            J = J + 1
            RWORK( J ) = REAL( B( JROW, JCOL ) )
         ENDDO
      ENDDO
      CALL SGEMM( 'T', 'N', N, NRHS, N, 1.0E+0, RWORK( IRWVT ), N, &
                  RWORK( IRWB ), N, 0.0E+0, RWORK( IRWRB ), N )
      J = IRWB - 1
      DO JCOL = 1, NRHS
         DO JROW = 1, N
            J = J + 1
            RWORK( J ) = AIMAG( B( JROW, JCOL ) )
         ENDDO
      ENDDO
      CALL SGEMM( 'T', 'N', N, NRHS, N, 1.0E+0, RWORK( IRWVT ), N, &
                  RWORK( IRWB ), N, 0.0E+0, RWORK( IRWIB ), N )
      JREAL = IRWRB - 1
      JIMAG = IRWIB - 1
      DO JCOL = 1, NRHS
         DO JROW = 1, N
            JREAL = JREAL + 1
            JIMAG = JIMAG + 1
            B( JROW, JCOL ) = CMPLX( RWORK( JREAL ), RWORK( JIMAG ) )
         ENDDO
      ENDDO
!
!        Unscale.
!
      CALL SLASCL( 'G', 0, 0, 1.0E+0, ORGNRM, N, 1, D, N, INFO )
      CALL SLASRT( 'D', N, D, INFO )
      CALL CLASCL( 'G', 0, 0, ORGNRM, 1.0E+0, N, NRHS, B, LDB, INFO )
!
      RETURN
   END IF
!
!     Book-keeping and setting up some constants.
!
   NLVL = INT( LOG( REAL( N ) / REAL( SMLSIZ+1 ) ) / LOG( 2.0E+0 ) ) + 1
!
   SMLSZP = SMLSIZ + 1
!
   U = 1
   VT = 1 + SMLSIZ*N
   DIFL = VT + SMLSZP*N
   DIFR = DIFL + NLVL*N
   Z = DIFR + NLVL*N*2
   C = Z + NLVL*N
   S = C + N
   POLES = S + N
   GIVNUM = POLES + 2*NLVL*N
   NRWORK = GIVNUM + 2*NLVL*N
   BX = 1
!
   IRWRB = NRWORK
   IRWIB = IRWRB + SMLSIZ*NRHS
   IRWB = IRWIB + SMLSIZ*NRHS
!
   SIZEI = 1 + N
   K = SIZEI + N
   GIVPTR = K + N
   PERM = GIVPTR + N
   GIVCOL = PERM + NLVL*N
   IWK = GIVCOL + NLVL*N*2
!
   ST = 1
   SQRE = 0
   ICMPQ1 = 1
   ICMPQ2 = 0
   NSUB = 0
!
   WHERE (ABS( D( 1:N ) ) < EPS )
      D( 1:N ) = SIGN( EPS, D( 1:N ) )
   END WHERE
!
   DO I = 1, NM1
      IF( ( ABS( E( I ) ) < EPS ) .OR. ( I == NM1 ) ) THEN
         NSUB = NSUB + 1
         IWORK( NSUB ) = ST
!
!           Subproblem found. First determine its size and then
!           apply divide and conquer on it.
!
         IF( I < NM1 ) THEN
!
!              A subproblem with E(I) small for I < NM1.
!
            NSIZE = I - ST + 1
            IWORK( SIZEI+NSUB-1 ) = NSIZE
         ELSE IF( ABS( E( I ) ) >= EPS ) THEN
!
!              A subproblem with E(NM1) not too small but I = NM1.
!
            NSIZE = N - ST + 1
            IWORK( SIZEI+NSUB-1 ) = NSIZE
         ELSE
!
!              A subproblem with E(NM1) small. This implies an
!              1-by-1 subproblem at D(N), which is not solved
!              explicitly.
!
            NSIZE = I - ST + 1
            IWORK( SIZEI+NSUB-1 ) = NSIZE
            NSUB = NSUB + 1
            IWORK( NSUB ) = N
            IWORK( SIZEI+NSUB-1 ) = 1
            CALL CCOPY( NRHS, B( N, 1 ), LDB, WORK( BX+NM1 ), N )
         END IF
         ST1 = ST - 1
         IF( NSIZE == 1 ) THEN
!
!              This is a 1-by-1 subproblem and is not solved
!              explicitly.
!
            CALL CCOPY( NRHS, B( ST, 1 ), LDB, WORK( BX+ST1 ), N )
         ELSE IF( NSIZE <= SMLSIZ ) THEN
!
!              This is a small subproblem and is solved by SLASDQ.
!
            CALL SLASET( 'A', NSIZE, NSIZE, 0.0E+0, 1.0E+0, &
                         RWORK( VT+ST1 ), N )
            CALL SLASET( 'A', NSIZE, NSIZE, 0.0E+0, 1.0E+0, &
                         RWORK( U+ST1 ), N )
            CALL SLASDQ( 'U', 0, NSIZE, NSIZE, NSIZE, 0, D( ST ), &
                         E( ST ), RWORK( VT+ST1 ), N, RWORK( U+ST1 ), &
                         N, RWORK( NRWORK ), 1, RWORK( NRWORK ), &
                         INFO )
            IF( INFO /= 0 ) RETURN
!
!              In the real version, B is passed to SLASDQ and multiplied
!              internally by Q**H. Here B is complex and that product is
!              computed below in two steps (real and imaginary parts).
!
            J = IRWB - 1
            DO JCOL = 1, NRHS
               DO JROW = ST, ST + NSIZE - 1
                  J = J + 1
                  RWORK( J ) = REAL( B( JROW, JCOL ) )
               ENDDO
            ENDDO
            CALL SGEMM( 'T', 'N', NSIZE, NRHS, NSIZE, 1.0E+0, &
                        RWORK( U+ST1 ), N, RWORK( IRWB ), NSIZE, &
                        0.0E+0, RWORK( IRWRB ), NSIZE )
            J = IRWB - 1
            DO JCOL = 1, NRHS
               DO JROW = ST, ST + NSIZE - 1
                  J = J + 1
                  RWORK( J ) = AIMAG( B( JROW, JCOL ) )
               ENDDO
            ENDDO
            CALL SGEMM( 'T', 'N', NSIZE, NRHS, NSIZE, 1.0E+0, &
                        RWORK( U+ST1 ), N, RWORK( IRWB ), NSIZE, &
                        0.0E+0, RWORK( IRWIB ), NSIZE )
            JREAL = IRWRB - 1
            JIMAG = IRWIB - 1
            DO JCOL = 1, NRHS
               DO JROW = ST, ST + NSIZE - 1
                  JREAL = JREAL + 1
                  JIMAG = JIMAG + 1
                  B( JROW, JCOL ) = CMPLX( RWORK( JREAL ), &
                                    RWORK( JIMAG ) )
               ENDDO
            ENDDO
!
            CALL CLACPY( 'A', NSIZE, NRHS, B( ST, 1 ), LDB, &
                         WORK( BX+ST1 ), N )
         ELSE
!
!              A large problem. Solve it using divide and conquer.
!
            CALL SLASDA( ICMPQ1, SMLSIZ, NSIZE, SQRE, D( ST ), &
                         E( ST ), RWORK( U+ST1 ), N, RWORK( VT+ST1 ), &
                         IWORK( K+ST1 ), RWORK( DIFL+ST1 ), &
                         RWORK( DIFR+ST1 ), RWORK( Z+ST1 ), &
                         RWORK( POLES+ST1 ), IWORK( GIVPTR+ST1 ), &
                         IWORK( GIVCOL+ST1 ), N, IWORK( PERM+ST1 ), &
                         RWORK( GIVNUM+ST1 ), RWORK( C+ST1 ), &
                         RWORK( S+ST1 ), RWORK( NRWORK ), &
                         IWORK( IWK ), INFO )
            IF( INFO /= 0 ) RETURN
            BXST = BX + ST1
            CALL CLALSA( ICMPQ2, SMLSIZ, NSIZE, NRHS, B( ST, 1 ), &
                         LDB, WORK( BXST ), N, RWORK( U+ST1 ), N, &
                         RWORK( VT+ST1 ), IWORK( K+ST1 ), &
                         RWORK( DIFL+ST1 ), RWORK( DIFR+ST1 ), &
                         RWORK( Z+ST1 ), RWORK( POLES+ST1 ), &
                         IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N, &
                         IWORK( PERM+ST1 ), RWORK( GIVNUM+ST1 ), &
                         RWORK( C+ST1 ), RWORK( S+ST1 ), &
                         RWORK( NRWORK ), IWORK( IWK ), INFO )
            IF( INFO /= 0 ) RETURN
         END IF
         ST = I + 1
      END IF
   ENDDO
!
!     Apply the singular values and treat the tiny ones as zero.
!
   TOL = RCND*ABS( D( ISAMAX( N, D, 1 ) ) )
!
   DO I = 1, N
!
!        Some of the elements in D can be negative because 1-by-1
!        subproblems were not solved explicitly.
!
      IF( ABS( D( I ) ) <= TOL ) THEN
         CALL CLASET( 'A', 1, NRHS, (0.0E+0,0.0E+0), (0.0E+0,0.0E+0), WORK( BX+I-1 ), N )
      ELSE
         RANK = RANK + 1
         CALL CLASCL( 'G', 0, 0, D( I ), 1.0E+0, 1, NRHS, WORK( BX+I-1 ), N, INFO )
      END IF
      D( I ) = ABS( D( I ) )
   ENDDO
!
!     Now apply back the right singular vectors.
!
   ICMPQ2 = 1
   DO I = 1, NSUB
      ST = IWORK( I )
      ST1 = ST - 1
      NSIZE = IWORK( SIZEI+I-1 )
      BXST = BX + ST1
      IF( NSIZE == 1 ) THEN
         CALL CCOPY( NRHS, WORK( BXST ), N, B( ST, 1 ), LDB )
      ELSE IF( NSIZE <= SMLSIZ ) THEN
!
!           Since B and BX are complex, the following call to SGEMM
!           is performed in two steps (real and imaginary parts).
!
!           CALL SGEMM( 'T', 'N', NSIZE, NRHS, NSIZE, 1.0E+0,
!    $                  RWORK( VT+ST1 ), N, RWORK( BXST ), N, 0.0E+0,
!    $                  B( ST, 1 ), LDB )
!
         J = BXST - N - 1
         JREAL = IRWB - 1
         DO JCOL = 1, NRHS
            J = J + N
            DO JROW = 1, NSIZE
               JREAL = JREAL + 1
               RWORK( JREAL ) = REAL( WORK( J+JROW ) )
            ENDDO
         ENDDO
         CALL SGEMM( 'T', 'N', NSIZE, NRHS, NSIZE, 1.0E+0, &
                     RWORK( VT+ST1 ), N, RWORK( IRWB ), NSIZE, 0.0E+0, &
                     RWORK( IRWRB ), NSIZE )
         J = BXST - N - 1
         JIMAG = IRWB - 1
         DO JCOL = 1, NRHS
            J = J + N
            DO JROW = 1, NSIZE
               JIMAG = JIMAG + 1
               RWORK( JIMAG ) = AIMAG( WORK( J+JROW ) )
            ENDDO
         ENDDO
         CALL SGEMM( 'T', 'N', NSIZE, NRHS, NSIZE, 1.0E+0, &
                     RWORK( VT+ST1 ), N, RWORK( IRWB ), NSIZE, 0.0E+0, &
                     RWORK( IRWIB ), NSIZE )
         JREAL = IRWRB - 1
         JIMAG = IRWIB - 1
         DO JCOL = 1, NRHS
            DO JROW = ST, ST + NSIZE - 1
               JREAL = JREAL + 1
               JIMAG = JIMAG + 1
               B( JROW, JCOL ) = CMPLX( RWORK( JREAL ), &
                                 RWORK( JIMAG ) )
            ENDDO
         ENDDO
      ELSE
         CALL CLALSA( ICMPQ2, SMLSIZ, NSIZE, NRHS, WORK( BXST ), N, &
                      B( ST, 1 ), LDB, RWORK( U+ST1 ), N, &
                      RWORK( VT+ST1 ), IWORK( K+ST1 ), &
                      RWORK( DIFL+ST1 ), RWORK( DIFR+ST1 ), &
                      RWORK( Z+ST1 ), RWORK( POLES+ST1 ), &
                      IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N, &
                      IWORK( PERM+ST1 ), RWORK( GIVNUM+ST1 ), &
                      RWORK( C+ST1 ), RWORK( S+ST1 ), &
                      RWORK( NRWORK ), IWORK( IWK ), INFO )
         IF( INFO /= 0 ) RETURN
      END IF
   ENDDO
!
!     Unscale and sort the singular values.
!
   CALL SLASCL( 'G', 0, 0, 1.0E+0, ORGNRM, N, 1, D, N, INFO )
   CALL SLASRT( 'D', N, D, INFO )
   CALL CLASCL( 'G', 0, 0, ORGNRM, 1.0E+0, N, NRHS, B, LDB, INFO )
!
   RETURN
!
!     End of CLALSD
!
END
