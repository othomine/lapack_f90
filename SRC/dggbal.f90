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
   SUBROUTINE DGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, WORK, INFO )
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
   DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), LSCALE( * ), RSCALE( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   SCLFAC
   PARAMETER          ( SCLFAC = 1.0D+1 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, ICAB, IFLOW, IP1, IR, IRAB, IT, J, JC, JP1, &
                      K, KOUNT, L, LCAB, LM1, LRAB, LSFMAX, LSFMIN, &
                      M, NR, NRP2
   DOUBLE PRECISION   ALPHA, BASL, BETA, CAB, CMAX, COEF, COEF2, &
                      COEF5, COR, EW, EWC, GAMMA, PGAMMA, RAB, SFMAX, &
                      SFMIN, SOMME, T, TA, TB, TC
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
      LSCALE( 1 ) = 1.0D0
      RSCALE( 1 ) = 1.0D0
      RETURN
   END IF
!
   IF( LSAME( JOB, 'N' ) ) THEN
      ILO = 1
      IHI = N
      LSCALE(1:N) = 1.0D0
      RSCALE(1:N) = 1.0D0
      RETURN
   END IF
!
   K = 1
   L = N
   IF( LSAME( JOB, 'S' ) ) GO TO 190
!
   GO TO 30
!
!     Permute the matrices A and B to isolate the eigenvalues.
!
!     Find row with one nonzero in columns 1 through L
!
20 CONTINUE
   L = LM1
   IF( L /= 1 ) GO TO 30
!
   RSCALE( 1 ) = 1.0D0
   LSCALE( 1 ) = 1.0D0
   GO TO 190
!
30 CONTINUE
   LM1 = L - 1
   DO I = L, 1, -1
      DO J = 1, LM1
         JP1 = J + 1
         IF( A( I, J ) /= 0.0D0 .OR. B( I, J ) /= 0.0D0 ) GO TO 50
      ENDDO
      J = L
      GO TO 70
!
50    CONTINUE
      DO J = JP1, L
         IF( A( I, J ) /= 0.0D0 .OR. B( I, J ) /= 0.0D0 ) GO TO 80
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
         IF( A( I, J ) /= 0.0D0 .OR. B( I, J ) /= 0.0D0 ) GO TO 120
         ENDDO
      I = L
      GO TO 140
  120    CONTINUE
      DO I = IP1, L
         IF( A( I, J ) /= 0.0D0 .OR. B( I, J ) /= 0.0D0 ) GO TO 150
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
   IF( I /= M ) THEN
      CALL DSWAP( N-K+1, A( I, K ), LDA, A( M, K ), LDA )
      CALL DSWAP( N-K+1, B( I, K ), LDB, B( M, K ), LDB )
   ENDIF
!
!     Permute columns M and J
!
   RSCALE( M ) = J
   IF( J /= M ) THEN
      CALL DSWAP( L, A( 1, J ), 1, A( 1, M ), 1 )
      CALL DSWAP( L, B( 1, J ), 1, B( 1, M ), 1 )
   ENDIF
!
   GO TO ( 20, 90 )IFLOW
!
  190 CONTINUE
   ILO = K
   IHI = L
!
   IF( LSAME( JOB, 'P' ) ) THEN
      LSCALE(ILO:IHI) = 1.0D0
      RSCALE(ILO:IHI) = 1.0D0
      RETURN
   END IF
!
   IF( ILO == IHI ) RETURN
!
!     Balance the submatrix in rows ILO to IHI.
!
   NR = IHI - ILO + 1
   RSCALE(ILO:IHI) = 0.0D0
   LSCALE(ILO:IHI) = 0.0D0
!
   WORK(ILO:IHI) = 0.0D0
   WORK(ILO+N:IHI+N ) = 0.0D0
   WORK(ILO+2*N:IHI+2*N ) = 0.0D0
   WORK(ILO+3*N:IHI+3*N ) = 0.0D0
   WORK(ILO+4*N:IHI+4*N ) = 0.0D0
   WORK(ILO+5*N:IHI+5*N ) = 0.0D0
!
!     Compute right side vector in resulting linear equations
!
   BASL = LOG10( SCLFAC )
   DO I = ILO, IHI
      DO J = ILO, IHI
         TB = B( I, J )
         TA = A( I, J )
         IF( TA /= 0.0D0 ) TA = LOG10( ABS( TA ) ) / BASL
         IF( TB /= 0.0D0 ) TB = LOG10( ABS( TB ) ) / BASL
         WORK( I+4*N ) = WORK( I+4*N ) - TA - TB
         WORK( J+5*N ) = WORK( J+5*N ) - TA - TB
      ENDDO
   ENDDO
!
   COEF = 1.0D0 / DBLE( 2*NR )
   COEF2 = COEF*COEF
   COEF5 = 0.5D0*COEF2
   NRP2 = NR + 2
   BETA = 0.0D0
   IT = 1
!
!     Start generalized conjugate gradient iteration
!
  250 CONTINUE
!
   GAMMA = SUM(WORK(ILO+4*N:ILO+4*N+NR-1)**2) + SUM(WORK(ILO+5*N:ILO+5*N+NR-1)**2)
!
   EW = SUM(WORK(ILO+4*N:IHI+4*N))
   EWC = SUM(WORK(ILO+5*N:IHI+5*N))
!
   GAMMA = COEF*GAMMA - COEF2*( EW**2+EWC**2 ) - COEF5*( EW-EWC )**2
   IF( GAMMA == 0.0D0 ) GO TO 350
   IF( IT /= 1 ) BETA = GAMMA / PGAMMA
   T = COEF5*( EWC-3.0D0*EW )
   TC = COEF5*( EW-3.0D0*EWC )
!
   WORK(ILO:ILO+NR-1) = BETA*WORK(ILO:ILO+NR-1)
   WORK(ILO+N:ILO+N+NR-1) = BETA*WORK(ILO+N:ILO+N+NR-1)
!
   WORK(ILO+N:ILO+N+NR-1) = WORK(ILO+N:ILO+N+NR-1) + COEF*WORK(ILO+4*N:ILO+4*N+NR-1)
   WORK(ILO:ILO+NR-1) = WORK(ILO:ILO+NR-1) + COEF*WORK(ILO+5*N:ILO+5*N+NR-1)
!
   WORK(ILO:IHI) = WORK(ILO:IHI) + TC
   WORK(ILO+N:IHI+N) = WORK(ILO+N:IHI+N) + T
!
!     Apply matrix to vector
!
   DO I = ILO, IHI
      KOUNT = 0
      SOMME = 0.0D0
      DO J = ILO, IHI
         IF( A( I, J ) /= 0.0D0 ) THEN
            KOUNT = KOUNT + 1
            SOMME = SOMME + WORK( J )
         ENDIF
         IF( B( I, J ) /= 0.0D0 ) THEN
            KOUNT = KOUNT + 1
            SOMME = SOMME + WORK( J )
         ENDIF
      ENDDO
      WORK( I+2*N ) = DBLE( KOUNT )*WORK( I+N ) + SOMME
   ENDDO
!
   DO J = ILO, IHI
      KOUNT = 0
      SOMME = 0.0D0
      DO I = ILO, IHI
         IF( A( I, J ) /= 0.0D0 ) THEN
            KOUNT = KOUNT + 1
            SOMME = SOMME + WORK( I+N )
         ENDIF
         IF( B( I, J ) /= 0.0D0 ) THEN
            KOUNT = KOUNT + 1
            SOMME = SOMME + WORK( I+N )
         ENDIF
      ENDDO
      WORK( J+3*N ) = DBLE( KOUNT )*WORK( J ) + SOMME
   ENDDO
!
   SOMME = SUM(WORK(ILO+N:ILO+N+NR-1)*WORK(ILO+2*N:ILO+2*N+NR-1)) + &
           SUM(WORK(ILO:ILO+NR-1)*WORK(ILO+3*N:ILO+3*N+NR-1))
   ALPHA = GAMMA / SOMME
!
!     Determine correction to current iteration
!
   CMAX = 0.0D0
   DO I = ILO, IHI
      COR = ALPHA*WORK( I+N )
      IF( ABS( COR ) > CMAX ) CMAX = ABS( COR )
      LSCALE( I ) = LSCALE( I ) + COR
      COR = ALPHA*WORK( I )
      IF( ABS( COR ) > CMAX ) CMAX = ABS( COR )
      RSCALE( I ) = RSCALE( I ) + COR
   ENDDO
   IF( CMAX < 0.5D0 ) GO TO 350
!
   CALL DAXPY( NR, -ALPHA, WORK( ILO+2*N ), 1, WORK( ILO+4*N ), 1 )
   CALL DAXPY( NR, -ALPHA, WORK( ILO+3*N ), 1, WORK( ILO+5*N ), 1 )
!
   PGAMMA = GAMMA
   IT = IT + 1
   IF( IT <= NRP2 ) GO TO 250
!
!     End generalized conjugate gradient iteration
!
  350 CONTINUE
   SFMIN = DLAMCH( 'S' )
   SFMAX = 1.0D0 / SFMIN
   LSFMIN = INT( LOG10( SFMIN ) / BASL+1.0D0 )
   LSFMAX = INT( LOG10( SFMAX ) / BASL )
   DO I = ILO, IHI
      IRAB = IDAMAX( N-ILO+1, A( I, ILO ), LDA )
      RAB = ABS( A( I, IRAB+ILO-1 ) )
      IRAB = IDAMAX( N-ILO+1, B( I, ILO ), LDB )
      RAB = MAX( RAB, ABS( B( I, IRAB+ILO-1 ) ) )
      LRAB = INT( LOG10( RAB+SFMIN ) / BASL+1.0D0 )
      IR = INT(LSCALE( I ) + SIGN( 0.5D0, LSCALE( I ) ))
      IR = MIN( MAX( IR, LSFMIN ), LSFMAX, LSFMAX-LRAB )
      LSCALE( I ) = SCLFAC**IR
      ICAB = IDAMAX( IHI, A( 1, I ), 1 )
      CAB = ABS( A( ICAB, I ) )
      ICAB = IDAMAX( IHI, B( 1, I ), 1 )
      CAB = MAX( CAB, ABS( B( ICAB, I ) ) )
      LCAB = INT( LOG10( CAB+SFMIN ) / BASL+1.0D0 )
      JC = INT(RSCALE( I ) + SIGN( 0.5D0, RSCALE( I ) ))
      JC = MIN( MAX( JC, LSFMIN ), LSFMAX, LSFMAX-LCAB )
      RSCALE( I ) = SCLFAC**JC
   ENDDO
!
!     Row scaling of matrices A and B
!
   DO I = ILO, IHI
      A(I,ILO:N) = LSCALE(I)*A(I,ILO:N)
      B(I,ILO:N) = LSCALE(I)*B(I,ILO:N)
   ENDDO
!
!     Column scaling of matrices A and B
!
   DO J = ILO, IHI
      A(1:IHI,J) = RSCALE(J)*A(1:IHI,J)
      B(1:IHI,J) = RSCALE(J)*B(1:IHI,J)
   ENDDO
!
   RETURN
!
!     End of DGGBAL
!
END
