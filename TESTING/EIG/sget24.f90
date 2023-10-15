!> \brief \b SGET24
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGET24( COMP, JTYPE, THRESH, ISEED, NOUNIT, N, A, LDA,
!                          H, HT, WR, WI, WRT, WIT, WRTMP, WITMP, VS,
!                          LDVS, VS1, RCDEIN, RCDVIN, NSLCT, ISLCT,
!                          RESULT, WORK, LWORK, IWORK, BWORK, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            COMP
!       INTEGER            INFO, JTYPE, LDA, LDVS, LWORK, N, NOUNIT, NSLCT
!       REAL               RCDEIN, RCDVIN, THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            BWORK( * )
!       INTEGER            ISEED( 4 ), ISLCT( * ), IWORK( * )
!       REAL               A( LDA, * ), H( LDA, * ), HT( LDA, * ),
!      $                   RESULT( 17 ), VS( LDVS, * ), VS1( LDVS, * ),
!      $                   WI( * ), WIT( * ), WITMP( * ), WORK( * ),
!      $                   WR( * ), WRT( * ), WRTMP( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SGET24 checks the nonsymmetric eigenvalue (Schur form) problem
!>    expert driver SGEESX.
!>
!>    If COMP = .FALSE., the first 13 of the following tests will be
!>    be performed on the input matrix A, and also tests 14 and 15
!>    if LWORK is sufficiently large.
!>    If COMP = .TRUE., all 17 test will be performed.
!>
!>    (1)     0 if T is in Schur form, 1/ulp otherwise
!>           (no sorting of eigenvalues)
!>
!>    (2)     | A - VS T VS' | / ( n |A| ulp )
!>
!>      Here VS is the matrix of Schur eigenvectors, and T is in Schur
!>      form  (no sorting of eigenvalues).
!>
!>    (3)     | I - VS VS' | / ( n ulp ) (no sorting of eigenvalues).
!>
!>    (4)     0     if WR+sqrt(-1)*WI are eigenvalues of T
!>            1/ulp otherwise
!>            (no sorting of eigenvalues)
!>
!>    (5)     0     if T(with VS) = T(without VS),
!>            1/ulp otherwise
!>            (no sorting of eigenvalues)
!>
!>    (6)     0     if eigenvalues(with VS) = eigenvalues(without VS),
!>            1/ulp otherwise
!>            (no sorting of eigenvalues)
!>
!>    (7)     0 if T is in Schur form, 1/ulp otherwise
!>            (with sorting of eigenvalues)
!>
!>    (8)     | A - VS T VS' | / ( n |A| ulp )
!>
!>      Here VS is the matrix of Schur eigenvectors, and T is in Schur
!>      form  (with sorting of eigenvalues).
!>
!>    (9)     | I - VS VS' | / ( n ulp ) (with sorting of eigenvalues).
!>
!>    (10)    0     if WR+sqrt(-1)*WI are eigenvalues of T
!>            1/ulp otherwise
!>            If workspace sufficient, also compare WR, WI with and
!>            without reciprocal condition numbers
!>            (with sorting of eigenvalues)
!>
!>    (11)    0     if T(with VS) = T(without VS),
!>            1/ulp otherwise
!>            If workspace sufficient, also compare T with and without
!>            reciprocal condition numbers
!>            (with sorting of eigenvalues)
!>
!>    (12)    0     if eigenvalues(with VS) = eigenvalues(without VS),
!>            1/ulp otherwise
!>            If workspace sufficient, also compare VS with and without
!>            reciprocal condition numbers
!>            (with sorting of eigenvalues)
!>
!>    (13)    if sorting worked and SDIM is the number of
!>            eigenvalues which were SELECTed
!>            If workspace sufficient, also compare SDIM with and
!>            without reciprocal condition numbers
!>
!>    (14)    if RCONDE the same no matter if VS and/or RCONDV computed
!>
!>    (15)    if RCONDV the same no matter if VS and/or RCONDE computed
!>
!>    (16)  |RCONDE - RCDEIN| / cond(RCONDE)
!>
!>       RCONDE is the reciprocal average eigenvalue condition number
!>       computed by SGEESX and RCDEIN (the precomputed true value)
!>       is supplied as input.  cond(RCONDE) is the condition number
!>       of RCONDE, and takes errors in computing RCONDE into account,
!>       so that the resulting quantity should be O(ULP). cond(RCONDE)
!>       is essentially given by norm(A)/RCONDV.
!>
!>    (17)  |RCONDV - RCDVIN| / cond(RCONDV)
!>
!>       RCONDV is the reciprocal right invariant subspace condition
!>       number computed by SGEESX and RCDVIN (the precomputed true
!>       value) is supplied as input. cond(RCONDV) is the condition
!>       number of RCONDV, and takes errors in computing RCONDV into
!>       account, so that the resulting quantity should be O(ULP).
!>       cond(RCONDV) is essentially given by norm(A)/RCONDE.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] COMP
!> \verbatim
!>          COMP is LOGICAL
!>          COMP describes which input tests to perform:
!>            = .FALSE. if the computed condition numbers are not to
!>                      be tested against RCDVIN and RCDEIN
!>            = .TRUE.  if they are to be compared
!> \endverbatim
!>
!> \param[in] JTYPE
!> \verbatim
!>          JTYPE is INTEGER
!>          Type of input matrix. Used to label output if error occurs.
!> \endverbatim
!>
!> \param[in] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          If COMP = .FALSE., the random number generator seed
!>          used to produce matrix.
!>          If COMP = .TRUE., ISEED(1) = the number of the example.
!>          Used to label output if error occurs.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is REAL
!>          A test will count as "failed" if the "error", computed as
!>          described above, exceeds THRESH.  Note that the error
!>          is scaled to be O(1), so THRESH should be a reasonably
!>          small multiple of 1, e.g., 10 or 100.  In particular,
!>          it should not depend on the precision (single vs. double)
!>          or the size of the matrix.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] NOUNIT
!> \verbatim
!>          NOUNIT is INTEGER
!>          The FORTRAN unit number for printing out error messages
!>          (e.g., if a routine returns INFO not equal to 0.)
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The dimension of A. N must be at least 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA, N)
!>          Used to hold the matrix whose eigenvalues are to be
!>          computed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A, and H. LDA must be at
!>          least 1 and at least N.
!> \endverbatim
!>
!> \param[out] H
!> \verbatim
!>          H is REAL array, dimension (LDA, N)
!>          Another copy of the test matrix A, modified by SGEESX.
!> \endverbatim
!>
!> \param[out] HT
!> \verbatim
!>          HT is REAL array, dimension (LDA, N)
!>          Yet another copy of the test matrix A, modified by SGEESX.
!> \endverbatim
!>
!> \param[out] WR
!> \verbatim
!>          WR is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is REAL array, dimension (N)
!>
!>          The real and imaginary parts of the eigenvalues of A.
!>          On exit, WR + WI*i are the eigenvalues of the matrix in A.
!> \endverbatim
!>
!> \param[out] WRT
!> \verbatim
!>          WRT is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] WIT
!> \verbatim
!>          WIT is REAL array, dimension (N)
!>
!>          Like WR, WI, these arrays contain the eigenvalues of A,
!>          but those computed when SGEESX only computes a partial
!>          eigendecomposition, i.e. not Schur vectors
!> \endverbatim
!>
!> \param[out] WRTMP
!> \verbatim
!>          WRTMP is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] WITMP
!> \verbatim
!>          WITMP is REAL array, dimension (N)
!>
!>          Like WR, WI, these arrays contain the eigenvalues of A,
!>          but sorted by increasing real part.
!> \endverbatim
!>
!> \param[out] VS
!> \verbatim
!>          VS is REAL array, dimension (LDVS, N)
!>          VS holds the computed Schur vectors.
!> \endverbatim
!>
!> \param[in] LDVS
!> \verbatim
!>          LDVS is INTEGER
!>          Leading dimension of VS. Must be at least max(1, N).
!> \endverbatim
!>
!> \param[out] VS1
!> \verbatim
!>          VS1 is REAL array, dimension (LDVS, N)
!>          VS1 holds another copy of the computed Schur vectors.
!> \endverbatim
!>
!> \param[in] RCDEIN
!> \verbatim
!>          RCDEIN is REAL
!>          When COMP = .TRUE. RCDEIN holds the precomputed reciprocal
!>          condition number for the average of selected eigenvalues.
!> \endverbatim
!>
!> \param[in] RCDVIN
!> \verbatim
!>          RCDVIN is REAL
!>          When COMP = .TRUE. RCDVIN holds the precomputed reciprocal
!>          condition number for the selected right invariant subspace.
!> \endverbatim
!>
!> \param[in] NSLCT
!> \verbatim
!>          NSLCT is INTEGER
!>          When COMP = .TRUE. the number of selected eigenvalues
!>          corresponding to the precomputed values RCDEIN and RCDVIN.
!> \endverbatim
!>
!> \param[in] ISLCT
!> \verbatim
!>          ISLCT is INTEGER array, dimension (NSLCT)
!>          When COMP = .TRUE. ISLCT selects the eigenvalues of the
!>          input matrix corresponding to the precomputed values RCDEIN
!>          and RCDVIN. For I=1, ... ,NSLCT, if ISLCT(I) = J, then the
!>          eigenvalue with the J-th largest real part is selected.
!>          Not referenced if COMP = .FALSE.
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (17)
!>          The values computed by the 17 tests described above.
!>          The values are currently limited to 1/ulp, to avoid
!>          overflow.
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
!>          The number of entries in WORK to be passed to SGEESX. This
!>          must be at least 3*N, and N+N**2 if tests 14--16 are to
!>          be performed.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N*N)
!> \endverbatim
!>
!> \param[out] BWORK
!> \verbatim
!>          BWORK is LOGICAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          If 0,  successful exit.
!>          If <0, input parameter -INFO had an incorrect value.
!>          If >0, SGEESX returned an error code, the absolute
!>                 value of which is returned.
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
!> \ingroup single_eig
!
!  =====================================================================
   SUBROUTINE SGET24( COMP, JTYPE, THRESH, ISEED, NOUNIT, N, A, LDA, &
                      H, HT, WR, WI, WRT, WIT, WRTMP, WITMP, VS, &
                      LDVS, VS1, RCDEIN, RCDVIN, NSLCT, ISLCT, &
                      RESULT, WORK, LWORK, IWORK, BWORK, INFO )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            COMP
   INTEGER            INFO, JTYPE, LDA, LDVS, LWORK, N, NOUNIT, NSLCT
   REAL               RCDEIN, RCDVIN, THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            BWORK( * )
   INTEGER            ISEED( 4 ), ISLCT( * ), IWORK( * )
   REAL               A( LDA, * ), H( LDA, * ), HT( LDA, * ), &
                      RESULT( 17 ), VS( LDVS, * ), VS1( LDVS, * ), &
                      WI( * ), WIT( * ), WITMP( * ), WORK( * ), &
                      WR( * ), WRT( * ), WRTMP( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               EPSIN
   PARAMETER          ( EPSIN = 5.9605E-8 )
!     ..
!     .. Local Scalars ..
   CHARACTER          SORT
   INTEGER            I, IINFO, ISORT, ITMP, J, KMIN, KNTEIG, LIWORK, &
                      RSUB, SDIM, SDIM1
   REAL               ANORM, EPS, RCNDE1, RCNDV1, RCONDE, RCONDV, &
                      SMLNUM, TMP, TOL, TOLIN, ULP, ULPINV, V, VIMIN, &
                      VRMIN, WNORM
!     ..
!     .. Local Arrays ..
   INTEGER            IPNT( 20 )
!     ..
!     .. Arrays in Common ..
   LOGICAL            SELVAL( 20 )
   REAL               SELWI( 20 ), SELWR( 20 )
!     ..
!     .. Scalars in Common ..
   INTEGER            SELDIM, SELOPT
!     ..
!     .. Common blocks ..
   COMMON             / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
!     ..
!     .. External Functions ..
   LOGICAL            SSLECT
   REAL               SLAMCH, SLANGE
   EXTERNAL           SSLECT, SLAMCH, SLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           SCOPY, SGEESX, SGEMM, SLACPY, SORT01, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Check for errors
!
   INFO = 0
   IF( THRESH < 0.0E+0 ) THEN
      INFO = -3
   ELSE IF( NOUNIT <= 0 ) THEN
      INFO = -5
   ELSE IF( N < 0 ) THEN
      INFO = -6
   ELSE IF( LDA < 1 .OR. LDA < N ) THEN
      INFO = -8
   ELSE IF( LDVS < 1 .OR. LDVS < N ) THEN
      INFO = -18
   ELSE IF( LWORK < 3*N ) THEN
      INFO = -26
   END IF
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'SGET24', -INFO )
      RETURN
   END IF
!
!     Quick return if nothing to do
!
   RESULT(1:17) = -1.0E+0
!
   IF( N == 0 ) RETURN
!
!     Important constants
!
   SMLNUM = SLAMCH( 'Safe minimum' )
   ULP = SLAMCH( 'Precision' )
   ULPINV = 1.0E+0 / ULP
!
!     Perform tests (1)-(13)
!
   SELOPT = 0
   LIWORK = N*N
   DO ISORT = 0, 1
      IF( ISORT == 0 ) THEN
         SORT = 'N'
         RSUB = 0
      ELSE
         SORT = 'S'
         RSUB = 6
      END IF
!
!        Compute Schur form and Schur vectors, and test them
!
      CALL SLACPY( 'F', N, N, A, LDA, H, LDA )
      CALL SGEESX( 'V', SORT, SSLECT, 'N', N, H, LDA, SDIM, WR, WI, &
                   VS, LDVS, RCONDE, RCONDV, WORK, LWORK, IWORK, &
                   LIWORK, BWORK, IINFO )
      IF( IINFO /= 0 .AND. IINFO /= N+2 ) THEN
         RESULT( 1+RSUB ) = ULPINV
         IF( JTYPE /= 22 ) THEN
            WRITE( NOUNIT, FMT = 9998 )'SGEESX1', IINFO, N, JTYPE, &
               ISEED
         ELSE
            WRITE( NOUNIT, FMT = 9999 )'SGEESX1', IINFO, N, &
               ISEED( 1 )
         END IF
         INFO = ABS( IINFO )
         RETURN
      END IF
      IF( ISORT == 0 ) THEN
         CALL SCOPY( N, WR, 1, WRTMP, 1 )
         CALL SCOPY( N, WI, 1, WITMP, 1 )
      END IF
!
!        Do Test (1) or Test (7)
!
      RESULT( 1+RSUB ) = 0.0E+0
      DO J = 1, N - 2
         DO I = J + 2, N
            IF( H( I, J ) /= 0.0E+0 ) &
               RESULT( 1+RSUB ) = ULPINV
         ENDDO
      ENDDO
      DO I = 1, N - 2
         IF( H( I+1, I ) /= 0.0E+0 .AND. H( I+2, I+1 ) /= 0.0E+0 ) &
            RESULT( 1+RSUB ) = ULPINV
      ENDDO
      DO I = 1, N - 1
         IF( H( I+1, I ) /= 0.0E+0 ) THEN
            IF( H( I, I ) /= H( I+1, I+1 ) .OR. H( I, I+1 ) == &
                0.0E+0 .OR. SIGN( 1.0E+0, H( I+1, I ) ) == &
                SIGN( 1.0E+0, H( I, I+1 ) ) )RESULT( 1+RSUB ) = ULPINV
         END IF
      ENDDO
!
!        Test (2) or (8): Compute norm(A - Q*H*Q') / (norm(A) * N * ULP)
!
!        Copy A to VS1, used as workspace
!
      CALL SLACPY( ' ', N, N, A, LDA, VS1, LDVS )
!
!        Compute Q*H and store in HT.
!
      CALL SGEMM( 'No transpose', 'No transpose', N, N, N, 1.0E+0, VS, &
                  LDVS, H, LDA, 0.0E+0, HT, LDA )
!
!        Compute A - Q*H*Q'
!
      CALL SGEMM( 'No transpose', 'Transpose', N, N, N, -1.0E+0, HT, &
                  LDA, VS, LDVS, 1.0E+0, VS1, LDVS )
!
      ANORM = MAX( SLANGE( '1', N, N, A, LDA, WORK ), SMLNUM )
      WNORM = SLANGE( '1', N, N, VS1, LDVS, WORK )
!
      IF( ANORM > WNORM ) THEN
         RESULT( 2+RSUB ) = ( WNORM / ANORM ) / ( N*ULP )
      ELSE
         IF( ANORM < 1.0E+0 ) THEN
            RESULT( 2+RSUB ) = ( MIN( WNORM, N*ANORM ) / ANORM ) / &
                               ( N*ULP )
         ELSE
            RESULT( 2+RSUB ) = MIN( WNORM / ANORM, REAL( N ) ) / &
                               ( N*ULP )
         END IF
      END IF
!
!        Test (3) or (9):  Compute norm( I - Q'*Q ) / ( N * ULP )
!
      CALL SORT01( 'Columns', N, N, VS, LDVS, WORK, LWORK, &
                   RESULT( 3+RSUB ) )
!
!        Do Test (4) or Test (10)
!
      RESULT( 4+RSUB ) = 0.0E+0
      DO I = 1, N
         IF( H( I, I ) /= WR( I ) ) &
            RESULT( 4+RSUB ) = ULPINV
      ENDDO
      IF( N > 1 ) THEN
         IF( H( 2, 1 ) == 0.0E+0 .AND. WI( 1 ) /= 0.0E+0 ) &
            RESULT( 4+RSUB ) = ULPINV
         IF( H( N, N-1 ) == 0.0E+0 .AND. WI( N ) /= 0.0E+0 ) &
            RESULT( 4+RSUB ) = ULPINV
      END IF
      DO I = 1, N - 1
         IF( H( I+1, I ) /= 0.0E+0 ) THEN
            TMP = SQRT( ABS( H( I+1, I ) ) )* &
                  SQRT( ABS( H( I, I+1 ) ) )
            RESULT( 4+RSUB ) = MAX( RESULT( 4+RSUB ), &
                               ABS( WI( I )-TMP ) / &
                               MAX( ULP*TMP, SMLNUM ) )
            RESULT( 4+RSUB ) = MAX( RESULT( 4+RSUB ), &
                               ABS( WI( I+1 )+TMP ) / &
                               MAX( ULP*TMP, SMLNUM ) )
         ELSE IF( I > 1 ) THEN
            IF( H( I+1, I ) == 0.0E+0 .AND. H( I, I-1 ) == 0.0E+0 .AND. &
                WI( I ) /= 0.0E+0 )RESULT( 4+RSUB ) = ULPINV
         END IF
      ENDDO
!
!        Do Test (5) or Test (11)
!
      CALL SLACPY( 'F', N, N, A, LDA, HT, LDA )
      CALL SGEESX( 'N', SORT, SSLECT, 'N', N, HT, LDA, SDIM, WRT, &
                   WIT, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, IWORK, &
                   LIWORK, BWORK, IINFO )
      IF( IINFO /= 0 .AND. IINFO /= N+2 ) THEN
         RESULT( 5+RSUB ) = ULPINV
         IF( JTYPE /= 22 ) THEN
            WRITE( NOUNIT, FMT = 9998 )'SGEESX2', IINFO, N, JTYPE, &
               ISEED
         ELSE
            WRITE( NOUNIT, FMT = 9999 )'SGEESX2', IINFO, N, &
               ISEED( 1 )
         END IF
         INFO = ABS( IINFO )
         GO TO 250
      END IF
!
      RESULT( 5+RSUB ) = 0.0E+0
      DO J = 1, N
         DO I = 1, N
            IF( H( I, J ) /= HT( I, J ) ) RESULT( 5+RSUB ) = ULPINV
         ENDDO
      ENDDO
!
!        Do Test (6) or Test (12)
!
      RESULT( 6+RSUB ) = 0.0E+0
      DO I = 1, N
         IF( WR( I ) /= WRT( I ) .OR. WI( I ) /= WIT( I ) ) &
            RESULT( 6+RSUB ) = ULPINV
         ENDDO
!
!        Do Test (13)
!
      IF( ISORT == 1 ) THEN
         RESULT( 13 ) = 0.0E+0
         KNTEIG = 0
         DO I = 1, N
            IF( SSLECT( WR( I ), WI( I ) ) .OR. &
                SSLECT( WR( I ), -WI( I ) ) )KNTEIG = KNTEIG + 1
            IF( I < N ) THEN
               IF( ( SSLECT( WR( I+1 ), WI( I+1 ) ) .OR. &
                   SSLECT( WR( I+1 ), -WI( I+1 ) ) ) .AND. &
                   ( .NOT.( SSLECT( WR( I ), &
                   WI( I ) ) .OR. SSLECT( WR( I ), &
                   -WI( I ) ) ) ) .AND. IINFO /= N+2 )RESULT( 13 ) &
                   = ULPINV
            END IF
         ENDDO
         IF( SDIM /= KNTEIG ) RESULT( 13 ) = ULPINV
      END IF
!
   ENDDO
!
!     If there is enough workspace, perform tests (14) and (15)
!     as well as (10) through (13)
!
   IF( LWORK >= N+( N*N ) / 2 ) THEN
!
!        Compute both RCONDE and RCONDV with VS
!
      SORT = 'S'
      RESULT( 14 ) = 0.0E+0
      RESULT( 15 ) = 0.0E+0
      CALL SLACPY( 'F', N, N, A, LDA, HT, LDA )
      CALL SGEESX( 'V', SORT, SSLECT, 'B', N, HT, LDA, SDIM1, WRT, &
                   WIT, VS1, LDVS, RCONDE, RCONDV, WORK, LWORK, &
                   IWORK, LIWORK, BWORK, IINFO )
      IF( IINFO /= 0 .AND. IINFO /= N+2 ) THEN
         RESULT( 14 ) = ULPINV
         RESULT( 15 ) = ULPINV
         IF( JTYPE /= 22 ) THEN
            WRITE( NOUNIT, FMT = 9998 )'SGEESX3', IINFO, N, JTYPE, &
               ISEED
         ELSE
            WRITE( NOUNIT, FMT = 9999 )'SGEESX3', IINFO, N, &
               ISEED( 1 )
         END IF
         INFO = ABS( IINFO )
         GO TO 250
      END IF
!
!        Perform tests (10), (11), (12), and (13)
!
      IF (ANY(WR(1:N) /= WRT(1:N) .OR. WI(1:N) /= WIT(1:N))) RESULT( 10 ) = ULPINV
      IF (ANY(H(1:N,1:N) /= HT(1:N,1:N))) RESULT( 11 ) = ULPINV
      IF (ANY(VS(1:N,1:N) /= VS1(1:N,1:N))) RESULT( 12 ) = ULPINV
      IF( SDIM /= SDIM1 ) RESULT( 13 ) = ULPINV
!
!        Compute both RCONDE and RCONDV without VS, and compare
!
      CALL SLACPY( 'F', N, N, A, LDA, HT, LDA )
      CALL SGEESX( 'N', SORT, SSLECT, 'B', N, HT, LDA, SDIM1, WRT, &
                   WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, &
                   IWORK, LIWORK, BWORK, IINFO )
      IF( IINFO /= 0 .AND. IINFO /= N+2 ) THEN
         RESULT( 14 ) = ULPINV
         RESULT( 15 ) = ULPINV
         IF( JTYPE /= 22 ) THEN
            WRITE( NOUNIT, FMT = 9998 )'SGEESX4', IINFO, N, JTYPE, ISEED
         ELSE
            WRITE( NOUNIT, FMT = 9999 )'SGEESX4', IINFO, N, ISEED( 1 )
         END IF
         INFO = ABS( IINFO )
         GO TO 250
      END IF
!
!        Perform tests (14) and (15)
!
      IF( RCNDE1 /= RCONDE ) RESULT( 14 ) = ULPINV
      IF( RCNDV1 /= RCONDV ) RESULT( 15 ) = ULPINV
!
!        Perform tests (10), (11), (12), and (13)
!
      IF (ANY(WR(1:N) /= WRT(1:N) .OR. WI(1:N) /= WIT(1:N))) RESULT( 10 ) = ULPINV
      IF (ANY(H(1:N,1:N) /= HT(1:N,1:N))) RESULT( 11 ) = ULPINV
      IF (ANY(VS(1:N,1:N) /= VS1(1:N,1:N))) RESULT( 12 ) = ULPINV
      IF( SDIM /= SDIM1 ) RESULT( 13 ) = ULPINV
!
!        Compute RCONDE with VS, and compare
!
      CALL SLACPY( 'F', N, N, A, LDA, HT, LDA )
      CALL SGEESX( 'V', SORT, SSLECT, 'E', N, HT, LDA, SDIM1, WRT, &
                   WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, &
                   IWORK, LIWORK, BWORK, IINFO )
      IF( IINFO /= 0 .AND. IINFO /= N+2 ) THEN
         RESULT( 14 ) = ULPINV
         IF( JTYPE /= 22 ) THEN
            WRITE( NOUNIT, FMT = 9998 )'SGEESX5', IINFO, N, JTYPE, &
               ISEED
         ELSE
            WRITE( NOUNIT, FMT = 9999 )'SGEESX5', IINFO, N, &
               ISEED( 1 )
         END IF
         INFO = ABS( IINFO )
         GO TO 250
      END IF
!
!        Perform test (14)
!
      IF( RCNDE1 /= RCONDE ) &
         RESULT( 14 ) = ULPINV
!
!        Perform tests (10), (11), (12), and (13)
!
      IF (ANY(WR(1:N) /= WRT(1:N) .OR. WI(1:N) /= WIT(1:N))) RESULT( 10 ) = ULPINV
      IF (ANY(H(1:N,1:N) /= HT(1:N,1:N))) RESULT( 11 ) = ULPINV
      IF (ANY(VS(1:N,1:N) /= VS1(1:N,1:N))) RESULT( 12 ) = ULPINV
      IF( SDIM /= SDIM1 ) RESULT( 13 ) = ULPINV
!
!        Compute RCONDE without VS, and compare
!
      CALL SLACPY( 'F', N, N, A, LDA, HT, LDA )
      CALL SGEESX( 'N', SORT, SSLECT, 'E', N, HT, LDA, SDIM1, WRT, &
                   WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, &
                   IWORK, LIWORK, BWORK, IINFO )
      IF( IINFO /= 0 .AND. IINFO /= N+2 ) THEN
         RESULT( 14 ) = ULPINV
         IF( JTYPE /= 22 ) THEN
            WRITE( NOUNIT, FMT = 9998 )'SGEESX6', IINFO, N, JTYPE, &
               ISEED
         ELSE
            WRITE( NOUNIT, FMT = 9999 )'SGEESX6', IINFO, N, &
               ISEED( 1 )
         END IF
         INFO = ABS( IINFO )
         GO TO 250
      END IF
!
!        Perform test (14)
!
      IF( RCNDE1 /= RCONDE ) RESULT( 14 ) = ULPINV
!
!        Perform tests (10), (11), (12), and (13)
!
      IF (ANY(WR(1:N) /= WRT(1:N) .OR. WI(1:N) /= WIT(1:N))) RESULT( 10 ) = ULPINV
      IF (ANY(H(1:N,1:N) /= HT(1:N,1:N))) RESULT( 11 ) = ULPINV
      IF (ANY(VS(1:N,1:N) /= VS1(1:N,1:N))) RESULT( 12 ) = ULPINV
      IF( SDIM /= SDIM1 ) RESULT( 13 ) = ULPINV
!
!        Compute RCONDV with VS, and compare
!
      CALL SLACPY( 'F', N, N, A, LDA, HT, LDA )
      CALL SGEESX( 'V', SORT, SSLECT, 'V', N, HT, LDA, SDIM1, WRT, &
                   WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, &
                   IWORK, LIWORK, BWORK, IINFO )
      IF( IINFO /= 0 .AND. IINFO /= N+2 ) THEN
         RESULT( 15 ) = ULPINV
         IF( JTYPE /= 22 ) THEN
            WRITE( NOUNIT, FMT = 9998 )'SGEESX7', IINFO, N, JTYPE, &
               ISEED
         ELSE
            WRITE( NOUNIT, FMT = 9999 )'SGEESX7', IINFO, N, &
               ISEED( 1 )
         END IF
         INFO = ABS( IINFO )
         GO TO 250
      END IF
!
!        Perform test (15)
!
      IF( RCNDV1 /= RCONDV ) RESULT( 15 ) = ULPINV
!
!        Perform tests (10), (11), (12), and (13)
!
      IF (ANY(WR(1:N) /= WRT(1:N) .OR. WI(1:N) /= WIT(1:N))) RESULT( 10 ) = ULPINV
      IF (ANY(H(1:N,1:N) /= HT(1:N,1:N))) RESULT( 11 ) = ULPINV
      IF (ANY(VS(1:N,1:N) /= VS1(1:N,1:N))) RESULT( 12 ) = ULPINV
      IF( SDIM /= SDIM1 ) RESULT( 13 ) = ULPINV
!
!        Compute RCONDV without VS, and compare
!
      CALL SLACPY( 'F', N, N, A, LDA, HT, LDA )
      CALL SGEESX( 'N', SORT, SSLECT, 'V', N, HT, LDA, SDIM1, WRT, &
                   WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, &
                   IWORK, LIWORK, BWORK, IINFO )
      IF( IINFO /= 0 .AND. IINFO /= N+2 ) THEN
         RESULT( 15 ) = ULPINV
         IF( JTYPE /= 22 ) THEN
            WRITE( NOUNIT, FMT = 9998 )'SGEESX8', IINFO, N, JTYPE, &
               ISEED
         ELSE
            WRITE( NOUNIT, FMT = 9999 )'SGEESX8', IINFO, N, &
               ISEED( 1 )
         END IF
         INFO = ABS( IINFO )
         GO TO 250
      END IF
!
!        Perform test (15)
!
      IF( RCNDV1 /= RCONDV ) RESULT( 15 ) = ULPINV
!
!        Perform tests (10), (11), (12), and (13)
!
      IF (ANY(WR(1:N) /= WRT(1:N) .OR. WI(1:N) /= WIT(1:N))) RESULT( 10 ) = ULPINV
      IF (ANY(H(1:N,1:N) /= HT(1:N,1:N))) RESULT( 11 ) = ULPINV
      IF (ANY(VS(1:N,1:N) /= VS1(1:N,1:N))) RESULT( 12 ) = ULPINV
      IF( SDIM /= SDIM1 ) RESULT( 13 ) = ULPINV
!
   END IF
!
  250 CONTINUE
!
!     If there are precomputed reciprocal condition numbers, compare
!     computed values with them.
!
   IF( COMP ) THEN
!
!        First set up SELOPT, SELDIM, SELVAL, SELWR, and SELWI so that
!        the logical function SSLECT selects the eigenvalues specified
!        by NSLCT and ISLCT.
!
      SELDIM = N
      SELOPT = 1
      EPS = MAX( ULP, EPSIN )
      DO I = 1, N
         IPNT( I ) = I
      ENDDO
      SELVAL(1:N) = .FALSE.
      SELWR(1:N) = WRTMP(1:N)
      SELWI(1:N) = WITMP(1:N)
      DO I = 1, N - 1
         KMIN = I
         VRMIN = WRTMP( I )
         VIMIN = WITMP( I )
         DO J = I + 1, N
            IF( WRTMP( J ) < VRMIN ) THEN
               KMIN = J
               VRMIN = WRTMP( J )
               VIMIN = WITMP( J )
            END IF
            ENDDO
         WRTMP( KMIN ) = WRTMP( I )
         WITMP( KMIN ) = WITMP( I )
         WRTMP( I ) = VRMIN
         WITMP( I ) = VIMIN
         ITMP = IPNT( I )
         IPNT( I ) = IPNT( KMIN )
         IPNT( KMIN ) = ITMP
      ENDDO
      SELVAL( IPNT( ISLCT( 1:NSLCT ) ) ) = .TRUE.
!
!        Compute condition numbers
!
      CALL SLACPY( 'F', N, N, A, LDA, HT, LDA )
      CALL SGEESX( 'N', 'S', SSLECT, 'B', N, HT, LDA, SDIM1, WRT, &
                   WIT, VS1, LDVS, RCONDE, RCONDV, WORK, LWORK, &
                   IWORK, LIWORK, BWORK, IINFO )
      IF( IINFO /= 0 .AND. IINFO /= N+2 ) THEN
         RESULT( 16 ) = ULPINV
         RESULT( 17 ) = ULPINV
         WRITE( NOUNIT, FMT = 9999 )'SGEESX9', IINFO, N, ISEED( 1 )
         INFO = ABS( IINFO )
         GO TO 300
      END IF
!
!        Compare condition number for average of selected eigenvalues
!        taking its condition number into account
!
      ANORM = SLANGE( '1', N, N, A, LDA, WORK )
      V = MAX( REAL( N )*EPS*ANORM, SMLNUM )
      IF( ANORM == 0.0E+0 ) &
         V = 1.0E+0
      IF( V > RCONDV ) THEN
         TOL = 1.0E+0
      ELSE
         TOL = V / RCONDV
      END IF
      IF( V > RCDVIN ) THEN
         TOLIN = 1.0E+0
      ELSE
         TOLIN = V / RCDVIN
      END IF
      TOL = MAX( TOL, SMLNUM / EPS )
      TOLIN = MAX( TOLIN, SMLNUM / EPS )
      IF( EPS*( RCDEIN-TOLIN ) > RCONDE+TOL ) THEN
         RESULT( 16 ) = ULPINV
      ELSE IF( RCDEIN-TOLIN > RCONDE+TOL ) THEN
         RESULT( 16 ) = ( RCDEIN-TOLIN ) / ( RCONDE+TOL )
      ELSE IF( RCDEIN+TOLIN < EPS*( RCONDE-TOL ) ) THEN
         RESULT( 16 ) = ULPINV
      ELSE IF( RCDEIN+TOLIN < RCONDE-TOL ) THEN
         RESULT( 16 ) = ( RCONDE-TOL ) / ( RCDEIN+TOLIN )
      ELSE
         RESULT( 16 ) = 1.0E+0
      END IF
!
!        Compare condition numbers for right invariant subspace
!        taking its condition number into account
!
      IF( V > RCONDV*RCONDE ) THEN
         TOL = RCONDV
      ELSE
         TOL = V / RCONDE
      END IF
      IF( V > RCDVIN*RCDEIN ) THEN
         TOLIN = RCDVIN
      ELSE
         TOLIN = V / RCDEIN
      END IF
      TOL = MAX( TOL, SMLNUM / EPS )
      TOLIN = MAX( TOLIN, SMLNUM / EPS )
      IF( EPS*( RCDVIN-TOLIN ) > RCONDV+TOL ) THEN
         RESULT( 17 ) = ULPINV
      ELSE IF( RCDVIN-TOLIN > RCONDV+TOL ) THEN
         RESULT( 17 ) = ( RCDVIN-TOLIN ) / ( RCONDV+TOL )
      ELSE IF( RCDVIN+TOLIN < EPS*( RCONDV-TOL ) ) THEN
         RESULT( 17 ) = ULPINV
      ELSE IF( RCDVIN+TOLIN < RCONDV-TOL ) THEN
         RESULT( 17 ) = ( RCONDV-TOL ) / ( RCDVIN+TOLIN )
      ELSE
         RESULT( 17 ) = 1.0E+0
      END IF
!
  300    CONTINUE
!
   END IF
!
 9999 FORMAT( ' SGET24: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', &
         I6, ', INPUT EXAMPLE NUMBER = ', I4 )
 9998 FORMAT( ' SGET24: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', &
         I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
!
   RETURN
!
!     End of SGET24
!
END

