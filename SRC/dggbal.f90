!> \brief \b DGGBAL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGGBAL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggbal.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggbal.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggbal.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE,
!                          RSCALE, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOB
!       INTEGER            IHI, ILO, INFO, LDA, LDB, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), LSCALE( * ),
!      $                   RSCALE( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGGBAL balances a pair of general real matrices (A,B).  This
!> involves, first, permuting A and B by similarity transformations to
!> isolate eigenvalues in the first 1 to ILO$-$1 and last IHI+1 to N
!> elements on the diagonal; and second, applying a diagonal similarity
!> transformation to rows and columns ILO to IHI to make the rows
!> and columns as close in norm as possible. Both steps are optional.
!>
!> Balancing may reduce the 1-norm of the matrices, and improve the
!> accuracy of the computed eigenvalues and/or eigenvectors in the
!> generalized eigenvalue problem A*x = lambda*B*x.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>          Specifies the operations to be performed on A and B:
!>          = 'N':  none:  simply set ILO = 1, IHI = N, LSCALE(I) = 1.0
!>                  and RSCALE(I) = 1.0 for i = 1,...,N.
!>          = 'P':  permute only;
!>          = 'S':  scale only;
!>          = 'B':  both permute and scale.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A and B.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the input matrix A.
!>          On exit,  A is overwritten by the balanced matrix.
!>          If JOB = 'N', A is not referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>          On entry, the input matrix B.
!>          On exit,  B is overwritten by the balanced matrix.
!>          If JOB = 'N', B is not referenced.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[out] IHI
!> \verbatim
!>          IHI is INTEGER
!>          ILO and IHI are set to integers such that on exit
!>          A(i,j) = 0 and B(i,j) = 0 if i > j and
!>          j = 1,...,ILO-1 or i = IHI+1,...,N.
!>          If JOB = 'N' or 'S', ILO = 1 and IHI = N.
!> \endverbatim
!>
!> \param[out] LSCALE
!> \verbatim
!>          LSCALE is DOUBLE PRECISION array, dimension (N)
!>          Details of the permutations and scaling factors applied
!>          to the left side of A and B.  If P(j) is the index of the
!>          row interchanged with row j, and D(j)
!>          is the scaling factor applied to row j, then
!>            LSCALE(j) = P(j)    for J = 1,...,ILO-1
!>                      = D(j)    for J = ILO,...,IHI
!>                      = P(j)    for J = IHI+1,...,N.
!>          The order in which the interchanges are made is N to IHI+1,
!>          then 1 to ILO-1.
!> \endverbatim
!>
!> \param[out] RSCALE
!> \verbatim
!>          RSCALE is DOUBLE PRECISION array, dimension (N)
!>          Details of the permutations and scaling factors applied
!>          to the right side of A and B.  If P(j) is the index of the
!>          column interchanged with column j, and D(j)
!>          is the scaling factor applied to column j, then
!>            LSCALE(j) = P(j)    for J = 1,...,ILO-1
!>                      = D(j)    for J = ILO,...,IHI
!>                      = P(j)    for J = IHI+1,...,N.
!>          The order in which the interchanges are made is N to IHI+1,
!>          then 1 to ILO-1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (lwork)
!>          lwork must be at least max(1,6*N) when JOB = 'S' or 'B', and
!>          at least 1 when JOB = 'N' or 'P'.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup ggbal
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  See R.C. WARD, Balancing the generalized eigenvalue problem,
!>                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE DGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE, &
                      RSCALE, WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          JOB
   INTEGER            IHI, ILO, INFO, LDA, LDB, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), LSCALE( * ), &
                      RSCALE( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, HALF, ONE
   PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
   DOUBLE PRECISION   THREE, SCLFAC
   PARAMETER          ( THREE = 3.0D+0, SCLFAC = 1.0D+1 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, ICAB, IFLOW, IP1, IR, IRAB, IT, J, JC, JP1, &
                      K, KOUNT, L, LCAB, LM1, LRAB, LSFMAX, LSFMIN, &
                      M, NR, NRP2
   DOUBLE PRECISION   ALPHA, BASL, BETA, CAB, CMAX, COEF, COEF2, &
                      COEF5, COR, EW, EWC, GAMMA, PGAMMA, RAB, SFMAX, &
                      SFMIN, SUM, T, TA, TB, TC
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            IDAMAX
   DOUBLE PRECISION   DDOT, DLAMCH
   EXTERNAL           LSAME, IDAMAX, DDOT, DLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           DAXPY, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, DBLE, INT, LOG10, MAX, MIN, SIGN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
   INFO = 0
   IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND. &
       .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -4
   ELSE IF( LDB < MAX( 1, N ) ) THEN
      INFO = -6
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'DGGBAL', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) THEN
      ILO = 1
      IHI = N
      RETURN
   END IF
!
   IF( N == 1 ) THEN
      ILO = 1
      IHI = N
      LSCALE( 1 ) = ONE
      RSCALE( 1 ) = ONE
      RETURN
   END IF
!
   IF( LSAME( JOB, 'N' ) ) THEN
      ILO = 1
      IHI = N
      DO I = 1, N
         LSCALE( I ) = ONE
         RSCALE( I ) = ONE
      ENDDO
      RETURN
   END IF
!
   K = 1
   L = N
   IF( LSAME( JOB, 'S' ) ) &
      GO TO 190
!
   GO TO 30
!
!     Permute the matrices A and B to isolate the eigenvalues.
!
!     Find row with one nonzero in columns 1 through L
!
20 CONTINUE
   L = LM1
   IF( L /= 1 ) &
      GO TO 30
!
   RSCALE( 1 ) = ONE
   LSCALE( 1 ) = ONE
   GO TO 190
!
30 CONTINUE
   LM1 = L - 1
   DO I = L, 1, -1
      DO J = 1, LM1
         JP1 = J + 1
         IF( A( I, J ) /= ZERO .OR. B( I, J ) /= ZERO ) &
            GO TO 50
      ENDDO
      J = L
      GO TO 70
!
50    CONTINUE
      DO J = JP1, L
         IF( A( I, J ) /= ZERO .OR. B( I, J ) /= ZERO ) &
            GO TO 80
      ENDDO
      J = JP1 - 1
!
70    CONTINUE
      M = L
      IFLOW = 1
      GO TO 160
80 CONTINUE
   ENDDO
   GO TO 100
!
!     Find column with one nonzero in rows K through N
!
90 CONTINUE
   K = K + 1
!
  100 CONTINUE
   DO J = K, L
      DO I = K, LM1
         IP1 = I + 1
         IF( A( I, J ) /= ZERO .OR. B( I, J ) /= ZERO ) &
            GO TO 120
         ENDDO
      I = L
      GO TO 140
  120    CONTINUE
      DO I = IP1, L
         IF( A( I, J ) /= ZERO .OR. B( I, J ) /= ZERO ) &
            GO TO 150
         ENDDO
      I = IP1 - 1
  140    CONTINUE
      M = K
      IFLOW = 2
      GO TO 160
  150 CONTINUE
      ENDDO
   GO TO 190
!
!     Permute rows M and I
!
  160 CONTINUE
   LSCALE( M ) = I
   IF( I == M ) &
      GO TO 170
   CALL DSWAP( N-K+1, A( I, K ), LDA, A( M, K ), LDA )
   CALL DSWAP( N-K+1, B( I, K ), LDB, B( M, K ), LDB )
!
!     Permute columns M and J
!
  170 CONTINUE
   RSCALE( M ) = J
   IF( J == M ) &
      GO TO 180
   CALL DSWAP( L, A( 1, J ), 1, A( 1, M ), 1 )
   CALL DSWAP( L, B( 1, J ), 1, B( 1, M ), 1 )
!
  180 CONTINUE
   GO TO ( 20, 90 )IFLOW
!
  190 CONTINUE
   ILO = K
   IHI = L
!
   IF( LSAME( JOB, 'P' ) ) THEN
      DO I = ILO, IHI
         LSCALE( I ) = ONE
         RSCALE( I ) = ONE
         ENDDO
      RETURN
   END IF
!
   IF( ILO == IHI ) &
      RETURN
!
!     Balance the submatrix in rows ILO to IHI.
!
   NR = IHI - ILO + 1
   DO I = ILO, IHI
      RSCALE( I ) = ZERO
      LSCALE( I ) = ZERO
!
      WORK( I ) = ZERO
      WORK( I+N ) = ZERO
      WORK( I+2*N ) = ZERO
      WORK( I+3*N ) = ZERO
      WORK( I+4*N ) = ZERO
      WORK( I+5*N ) = ZERO
      ENDDO
!
!     Compute right side vector in resulting linear equations
!
   BASL = LOG10( SCLFAC )
   DO I = ILO, IHI
      DO J = ILO, IHI
         TB = B( I, J )
         TA = A( I, J )
         IF( TA == ZERO ) &
            GO TO 210
         TA = LOG10( ABS( TA ) ) / BASL
  210       CONTINUE
         IF( TB == ZERO ) &
            GO TO 220
         TB = LOG10( ABS( TB ) ) / BASL
  220       CONTINUE
         WORK( I+4*N ) = WORK( I+4*N ) - TA - TB
         WORK( J+5*N ) = WORK( J+5*N ) - TA - TB
         ENDDO
      ENDDO
!
   COEF = ONE / DBLE( 2*NR )
   COEF2 = COEF*COEF
   COEF5 = HALF*COEF2
   NRP2 = NR + 2
   BETA = ZERO
   IT = 1
!
!     Start generalized conjugate gradient iteration
!
  250 CONTINUE
!
   GAMMA = DDOT( NR, WORK( ILO+4*N ), 1, WORK( ILO+4*N ), 1 ) + &
           DDOT( NR, WORK( ILO+5*N ), 1, WORK( ILO+5*N ), 1 )
!
   EW = ZERO
   EWC = ZERO
   DO I = ILO, IHI
      EW = EW + WORK( I+4*N )
      EWC = EWC + WORK( I+5*N )
      ENDDO
!
   GAMMA = COEF*GAMMA - COEF2*( EW**2+EWC**2 ) - COEF5*( EW-EWC )**2
   IF( GAMMA == ZERO ) &
      GO TO 350
   IF( IT /= 1 ) &
      BETA = GAMMA / PGAMMA
   T = COEF5*( EWC-THREE*EW )
   TC = COEF5*( EW-THREE*EWC )
!
   CALL DSCAL( NR, BETA, WORK( ILO ), 1 )
   CALL DSCAL( NR, BETA, WORK( ILO+N ), 1 )
!
   CALL DAXPY( NR, COEF, WORK( ILO+4*N ), 1, WORK( ILO+N ), 1 )
   CALL DAXPY( NR, COEF, WORK( ILO+5*N ), 1, WORK( ILO ), 1 )
!
   DO I = ILO, IHI
      WORK( I ) = WORK( I ) + TC
      WORK( I+N ) = WORK( I+N ) + T
      ENDDO
!
!     Apply matrix to vector
!
   DO I = ILO, IHI
      KOUNT = 0
      SUM = ZERO
      DO J = ILO, IHI
         IF( A( I, J ) == ZERO ) &
            GO TO 280
         KOUNT = KOUNT + 1
         SUM = SUM + WORK( J )
  280       CONTINUE
         IF( B( I, J ) == ZERO ) &
            GO TO 290
         KOUNT = KOUNT + 1
         SUM = SUM + WORK( J )
  290    CONTINUE
         ENDDO
      WORK( I+2*N ) = DBLE( KOUNT )*WORK( I+N ) + SUM
      ENDDO
!
   DO J = ILO, IHI
      KOUNT = 0
      SUM = ZERO
      DO I = ILO, IHI
         IF( A( I, J ) == ZERO ) &
            GO TO 310
         KOUNT = KOUNT + 1
         SUM = SUM + WORK( I+N )
  310       CONTINUE
         IF( B( I, J ) == ZERO ) &
            GO TO 320
         KOUNT = KOUNT + 1
         SUM = SUM + WORK( I+N )
  320    CONTINUE
         ENDDO
      WORK( J+3*N ) = DBLE( KOUNT )*WORK( J ) + SUM
      ENDDO
!
   SUM = DDOT( NR, WORK( ILO+N ), 1, WORK( ILO+2*N ), 1 ) + &
         DDOT( NR, WORK( ILO ), 1, WORK( ILO+3*N ), 1 )
   ALPHA = GAMMA / SUM
!
!     Determine correction to current iteration
!
   CMAX = ZERO
   DO I = ILO, IHI
      COR = ALPHA*WORK( I+N )
      IF( ABS( COR ) > CMAX ) &
         CMAX = ABS( COR )
      LSCALE( I ) = LSCALE( I ) + COR
      COR = ALPHA*WORK( I )
      IF( ABS( COR ) > CMAX ) &
         CMAX = ABS( COR )
      RSCALE( I ) = RSCALE( I ) + COR
      ENDDO
   IF( CMAX < HALF ) &
      GO TO 350
!
   CALL DAXPY( NR, -ALPHA, WORK( ILO+2*N ), 1, WORK( ILO+4*N ), 1 )
   CALL DAXPY( NR, -ALPHA, WORK( ILO+3*N ), 1, WORK( ILO+5*N ), 1 )
!
   PGAMMA = GAMMA
   IT = IT + 1
   IF( IT <= NRP2 ) &
      GO TO 250
!
!     End generalized conjugate gradient iteration
!
  350 CONTINUE
   SFMIN = DLAMCH( 'S' )
   SFMAX = ONE / SFMIN
   LSFMIN = INT( LOG10( SFMIN ) / BASL+ONE )
   LSFMAX = INT( LOG10( SFMAX ) / BASL )
   DO I = ILO, IHI
      IRAB = IDAMAX( N-ILO+1, A( I, ILO ), LDA )
      RAB = ABS( A( I, IRAB+ILO-1 ) )
      IRAB = IDAMAX( N-ILO+1, B( I, ILO ), LDB )
      RAB = MAX( RAB, ABS( B( I, IRAB+ILO-1 ) ) )
      LRAB = INT( LOG10( RAB+SFMIN ) / BASL+ONE )
      IR = INT(LSCALE( I ) + SIGN( HALF, LSCALE( I ) ))
      IR = MIN( MAX( IR, LSFMIN ), LSFMAX, LSFMAX-LRAB )
      LSCALE( I ) = SCLFAC**IR
      ICAB = IDAMAX( IHI, A( 1, I ), 1 )
      CAB = ABS( A( ICAB, I ) )
      ICAB = IDAMAX( IHI, B( 1, I ), 1 )
      CAB = MAX( CAB, ABS( B( ICAB, I ) ) )
      LCAB = INT( LOG10( CAB+SFMIN ) / BASL+ONE )
      JC = INT(RSCALE( I ) + SIGN( HALF, RSCALE( I ) ))
      JC = MIN( MAX( JC, LSFMIN ), LSFMAX, LSFMAX-LCAB )
      RSCALE( I ) = SCLFAC**JC
      ENDDO
!
!     Row scaling of matrices A and B
!
   DO I = ILO, IHI
      CALL DSCAL( N-ILO+1, LSCALE( I ), A( I, ILO ), LDA )
      CALL DSCAL( N-ILO+1, LSCALE( I ), B( I, ILO ), LDB )
      ENDDO
!
!     Column scaling of matrices A and B
!
   DO J = ILO, IHI
      CALL DSCAL( IHI, RSCALE( J ), A( 1, J ), 1 )
      CALL DSCAL( IHI, RSCALE( J ), B( 1, J ), 1 )
      ENDDO
!
   RETURN
!
!     End of DGGBAL
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

