!> \brief \b ZGET24
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGET24( COMP, JTYPE, THRESH, ISEED, NOUNIT, N, A, LDA,
!                          H, HT, W, WT, WTMP, VS, LDVS, VS1, RCDEIN,
!                          RCDVIN, NSLCT, ISLCT, ISRT, RESULT, WORK,
!                          LWORK, RWORK, BWORK, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            COMP
!       INTEGER            INFO, ISRT, JTYPE, LDA, LDVS, LWORK, N, NOUNIT,
!      $                   NSLCT
!       DOUBLE PRECISION   RCDEIN, RCDVIN, THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            BWORK( * )
!       INTEGER            ISEED( 4 ), ISLCT( * )
!       DOUBLE PRECISION   RESULT( 17 ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), H( LDA, * ), HT( LDA, * ),
!      $                   VS( LDVS, * ), VS1( LDVS, * ), W( * ),
!      $                   WORK( * ), WT( * ), WTMP( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ZGET24 checks the nonsymmetric eigenvalue (Schur form) problem
!>    expert driver ZGEESX.
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
!>    (4)     0     if W are eigenvalues of T
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
!>    (10)    0     if W are eigenvalues of T
!>            1/ulp otherwise
!>            If workspace sufficient, also compare W with and
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
!>       computed by ZGEESX and RCDEIN (the precomputed true value)
!>       is supplied as input.  cond(RCONDE) is the condition number
!>       of RCONDE, and takes errors in computing RCONDE into account,
!>       so that the resulting quantity should be O(ULP). cond(RCONDE)
!>       is essentially given by norm(A)/RCONDV.
!>
!>    (17)  |RCONDV - RCDVIN| / cond(RCONDV)
!>
!>       RCONDV is the reciprocal right invariant subspace condition
!>       number computed by ZGEESX and RCDVIN (the precomputed true
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
!>          THRESH is DOUBLE PRECISION
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
!>          A is COMPLEX*16 array, dimension (LDA, N)
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
!>          H is COMPLEX*16 array, dimension (LDA, N)
!>          Another copy of the test matrix A, modified by ZGEESX.
!> \endverbatim
!>
!> \param[out] HT
!> \verbatim
!>          HT is COMPLEX*16 array, dimension (LDA, N)
!>          Yet another copy of the test matrix A, modified by ZGEESX.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX*16 array, dimension (N)
!>          The computed eigenvalues of A.
!> \endverbatim
!>
!> \param[out] WT
!> \verbatim
!>          WT is COMPLEX*16 array, dimension (N)
!>          Like W, this array contains the eigenvalues of A,
!>          but those computed when ZGEESX only computes a partial
!>          eigendecomposition, i.e. not Schur vectors
!> \endverbatim
!>
!> \param[out] WTMP
!> \verbatim
!>          WTMP is COMPLEX*16 array, dimension (N)
!>          Like W, this array contains the eigenvalues of A,
!>          but sorted by increasing real or imaginary part.
!> \endverbatim
!>
!> \param[out] VS
!> \verbatim
!>          VS is COMPLEX*16 array, dimension (LDVS, N)
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
!>          VS1 is COMPLEX*16 array, dimension (LDVS, N)
!>          VS1 holds another copy of the computed Schur vectors.
!> \endverbatim
!>
!> \param[in] RCDEIN
!> \verbatim
!>          RCDEIN is DOUBLE PRECISION
!>          When COMP = .TRUE. RCDEIN holds the precomputed reciprocal
!>          condition number for the average of selected eigenvalues.
!> \endverbatim
!>
!> \param[in] RCDVIN
!> \verbatim
!>          RCDVIN is DOUBLE PRECISION
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
!>          eigenvalue with the J-th largest real or imaginary part is
!>          selected. The real part is used if ISRT = 0, and the
!>          imaginary part if ISRT = 1.
!>          Not referenced if COMP = .FALSE.
!> \endverbatim
!>
!> \param[in] ISRT
!> \verbatim
!>          ISRT is INTEGER
!>          When COMP = .TRUE., ISRT describes how ISLCT is used to
!>          choose a subset of the spectrum.
!>          Not referenced if COMP = .FALSE.
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (17)
!>          The values computed by the 17 tests described above.
!>          The values are currently limited to 1/ulp, to avoid
!>          overflow.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (2*N*N)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The number of entries in WORK to be passed to ZGEESX. This
!>          must be at least 2*N, and N*(N+1)/2 if tests 14--16 are to
!>          be performed.
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
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          If 0,  successful exit.
!>          If <0, input parameter -INFO had an incorrect value.
!>          If >0, ZGEESX returned an error code, the absolute
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
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
!
!> \ingroup complex16_eig
!
!  =====================================================================
   SUBROUTINE ZGET24( COMP, JTYPE, THRESH, ISEED, NOUNIT, N, A, LDA, &
                      H, HT, W, WT, WTMP, VS, LDVS, VS1, RCDEIN, &
                      RCDVIN, NSLCT, ISLCT, ISRT, RESULT, WORK, &
                      LWORK, RWORK, BWORK, INFO )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            COMP
   INTEGER            INFO, ISRT, JTYPE, LDA, LDVS, LWORK, N, NOUNIT, &
                      NSLCT
   DOUBLE PRECISION   RCDEIN, RCDVIN, THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            BWORK( * )
   INTEGER            ISEED( 4 ), ISLCT( * )
   DOUBLE PRECISION   RESULT( 17 ), RWORK( * )
   COMPLEX*16         A( LDA, * ), H( LDA, * ), HT( LDA, * ), &
                      VS( LDVS, * ), VS1( LDVS, * ), W( * ), &
                      WORK( * ), WT( * ), WTMP( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   EPSIN
   PARAMETER          ( EPSIN = 5.9605D-8 )
!     ..
!     .. Local Scalars ..
   CHARACTER          SORT
   INTEGER            I, IINFO, ISORT, ITMP, J, KMIN, KNTEIG, RSUB, &
                      SDIM, SDIM1
   DOUBLE PRECISION   ANORM, EPS, RCNDE1, RCNDV1, RCONDE, RCONDV, &
                      SMLNUM, TOL, TOLIN, ULP, ULPINV, V, VRICMP, &
                      VRIMIN, WNORM
   COMPLEX*16         CTMP
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   INTEGER            IPNT( 20 )
!     ..
!     .. External Functions ..
   LOGICAL            ZSLECT
   DOUBLE PRECISION   DLAMCH, ZLANGE
   EXTERNAL           ZSLECT, DLAMCH, ZLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA, ZCOPY, ZGEESX, ZGEMM, ZLACPY, ZUNT01
!     ..
!     .. Arrays in Common ..
   LOGICAL            SELVAL( 20 )
   DOUBLE PRECISION   SELWI( 20 ), SELWR( 20 )
!     ..
!     .. Scalars in Common ..
   INTEGER            SELDIM, SELOPT
!     ..
!     .. Common blocks ..
   COMMON             / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
!     ..
!     .. Executable Statements ..
!
!     Check for errors
!
   INFO = 0
   IF( THRESH < 0.0D0 ) THEN
      INFO = -3
   ELSE IF( NOUNIT <= 0 ) THEN
      INFO = -5
   ELSE IF( N < 0 ) THEN
      INFO = -6
   ELSE IF( LDA < 1 .OR. LDA < N ) THEN
      INFO = -8
   ELSE IF( LDVS < 1 .OR. LDVS < N ) THEN
      INFO = -15
   ELSE IF( LWORK < 2*N ) THEN
      INFO = -24
   END IF
!
   IF( INFO /= 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'ZGET24', -INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : XERBLA : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      RETURN
   END IF
!
!     Quick return if nothing to do
!
   RESULT(1:17) = -1.0D+0
!
   IF( N == 0 ) RETURN
!
!     Important constants
!
   SMLNUM = DLAMCH( 'Safe minimum' )
   ULP = DLAMCH( 'Precision' )
   ULPINV = 1.0D0 / ULP
!
!     Perform tests (1)-(13)
!
   SELOPT = 0
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
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLACPY( 'F', N, N, A, LDA, H, LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZGEESX( 'V', SORT, ZSLECT, 'N', N, H, LDA, SDIM, W, VS, &
                   LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, &
                   IINFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZGEESX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( IINFO /= 0 ) THEN
         RESULT( 1+RSUB ) = ULPINV
         IF( JTYPE /= 22 ) THEN
            WRITE( NOUNIT, FMT = 9998 )'ZGEESX1', IINFO, N, JTYPE, &
               ISEED
         ELSE
            WRITE( NOUNIT, FMT = 9999 )'ZGEESX1', IINFO, N, &
               ISEED( 1 )
         END IF
         INFO = ABS( IINFO )
         RETURN
      END IF
      IF( ISORT == 0 ) THEN
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZCOPY( N, W, 1, WTMP, 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      END IF
!
!        Do Test (1) or Test (7)
!
      RESULT( 1+RSUB ) = 0.0D0
      DO J = 1, N - 1
         DO I = J + 1, N
            IF( H( I, J ) /= (0.0D+0,0.0D+0) ) RESULT( 1+RSUB ) = ULPINV
         ENDDO
      ENDDO
!
!        Test (2) or (8): Compute norm(A - Q*H*Q') / (norm(A) * N * ULP)
!
!        Copy A to VS1, used as workspace
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLACPY( ' ', N, N, A, LDA, VS1, LDVS )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        Compute Q*H and store in HT.
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZGEMM( 'No transpose', 'No transpose', N, N, N, (1.0D0,0.0D0), VS, &
                  LDVS, H, LDA, (0.0D+0,0.0D+0), HT, LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        Compute A - Q*H*Q'
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZGEMM( 'No transpose', 'Conjugate transpose', N, N, N, &
                  -(1.0D0,0.0D0), HT, LDA, VS, LDVS, (1.0D0,0.0D0), VS1, LDVS )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
      ANORM = MAX( ZLANGE( '1', N, N, A, LDA, RWORK ), SMLNUM )
      WNORM = ZLANGE( '1', N, N, VS1, LDVS, RWORK )
!
      IF( ANORM > WNORM ) THEN
         RESULT( 2+RSUB ) = ( WNORM / ANORM ) / ( N*ULP )
      ELSE
         IF( ANORM < 1.0D0 ) THEN
            RESULT( 2+RSUB ) = ( MIN( WNORM, N*ANORM ) / ANORM ) / &
                               ( N*ULP )
         ELSE
            RESULT( 2+RSUB ) = MIN( WNORM / ANORM, DBLE( N ) ) / &
                               ( N*ULP )
         END IF
      END IF
!
!        Test (3) or (9):  Compute norm( I - Q'*Q ) / ( N * ULP )
!
      CALL ZUNT01( 'Columns', N, N, VS, LDVS, WORK, LWORK, RWORK, &
                   RESULT( 3+RSUB ) )
!
!        Do Test (4) or Test (10)
!
      RESULT( 4+RSUB ) = 0.0D0
      DO I = 1, N
         IF( H( I, I ) /= W( I ) ) &
            RESULT( 4+RSUB ) = ULPINV
      ENDDO
!
!        Do Test (5) or Test (11)
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLACPY( 'F', N, N, A, LDA, HT, LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZGEESX( 'N', SORT, ZSLECT, 'N', N, HT, LDA, SDIM, WT, VS, &
                   LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, &
                   IINFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZGEESX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( IINFO /= 0 ) THEN
         RESULT( 5+RSUB ) = ULPINV
         IF( JTYPE /= 22 ) THEN
            WRITE( NOUNIT, FMT = 9998 )'ZGEESX2', IINFO, N, JTYPE, &
               ISEED
         ELSE
            WRITE( NOUNIT, FMT = 9999 )'ZGEESX2', IINFO, N, &
               ISEED( 1 )
         END IF
         INFO = ABS( IINFO )
         GO TO 220
      END IF
!
      RESULT( 5+RSUB ) = 0.0D0
      DO J = 1, N
         DO I = 1, N
            IF( H( I, J ) /= HT( I, J ) ) &
               RESULT( 5+RSUB ) = ULPINV
         ENDDO
      ENDDO
!
!        Do Test (6) or Test (12)
!
      RESULT( 6+RSUB ) = 0.0E+0
      IF (ANY(W(1:N) /= WT(1:N))) RESULT( 6+RSUB ) = ULPINV
!
!        Do Test (13)
!
      IF( ISORT == 1 ) THEN
         RESULT( 13 ) = 0.0D0
         KNTEIG = 0
         DO I = 1, N
            IF( ZSLECT( W( I ) ) ) KNTEIG = KNTEIG + 1
            IF( I < N ) THEN
               IF( ZSLECT( W( I+1 ) ) .AND. &
                   ( .NOT.ZSLECT( W( I ) ) ) )RESULT( 13 ) = ULPINV
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
   IF( LWORK >= ( N*( N+1 ) ) / 2 ) THEN
!
!        Compute both RCONDE and RCONDV with VS
!
      SORT = 'S'
      RESULT( 14:15 ) = 0.0D0
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLACPY( 'F', N, N, A, LDA, HT, LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZGEESX( 'V', SORT, ZSLECT, 'B', N, HT, LDA, SDIM1, WT, &
                   VS1, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, &
                   BWORK, IINFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZGEESX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( IINFO /= 0 ) THEN
         RESULT( 14 ) = ULPINV
         RESULT( 15 ) = ULPINV
         IF( JTYPE /= 22 ) THEN
            WRITE( NOUNIT, FMT = 9998 )'ZGEESX3', IINFO, N, JTYPE, &
               ISEED
         ELSE
            WRITE( NOUNIT, FMT = 9999 )'ZGEESX3', IINFO, N, &
               ISEED( 1 )
         END IF
         INFO = ABS( IINFO )
         GO TO 220
      END IF
!
!        Perform tests (10), (11), (12), and (13)
!
      IF (ANY(W(1:N) /= WT(1:N))) RESULT( 10 ) = ULPINV
      IF (ANY(H(1:N,1:N) /= HT(1:N,1:N))) RESULT( 11 ) = ULPINV
      IF (ANY(VS(1:N,1:N) /= VS1(1:N,1:N))) RESULT( 12 ) = ULPINV
      IF( SDIM /= SDIM1 ) RESULT( 13 ) = ULPINV
!
!        Compute both RCONDE and RCONDV without VS, and compare
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLACPY( 'F', N, N, A, LDA, HT, LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZGEESX( 'N', SORT, ZSLECT, 'B', N, HT, LDA, SDIM1, WT, &
                   VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, &
                   BWORK, IINFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZGEESX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( IINFO /= 0 ) THEN
         RESULT( 14 ) = ULPINV
         RESULT( 15 ) = ULPINV
         IF( JTYPE /= 22 ) THEN
            WRITE( NOUNIT, FMT = 9998 )'ZGEESX4', IINFO, N, JTYPE, ISEED
         ELSE
            WRITE( NOUNIT, FMT = 9999 )'ZGEESX4', IINFO, N, ISEED( 1 )
         END IF
         INFO = ABS( IINFO )
         GO TO 220
      END IF
!
!        Perform tests (14) and (15)
!
      IF( RCNDE1 /= RCONDE ) RESULT( 14 ) = ULPINV
      IF( RCNDV1 /= RCONDV ) RESULT( 15 ) = ULPINV
!
!        Perform tests (10), (11), (12), and (13)
!
      IF (ANY(W(1:N) /= WT(1:N))) RESULT( 10 ) = ULPINV
      IF (ANY(H(1:N,1:N) /= HT(1:N,1:N))) RESULT( 11 ) = ULPINV
      IF (ANY(VS(1:N,1:N) /= VS1(1:N,1:N))) RESULT( 12 ) = ULPINV
      IF( SDIM /= SDIM1 ) RESULT( 13 ) = ULPINV
!
!        Compute RCONDE with VS, and compare
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLACPY( 'F', N, N, A, LDA, HT, LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZGEESX( 'V', SORT, ZSLECT, 'E', N, HT, LDA, SDIM1, WT, &
                   VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, &
                   BWORK, IINFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZGEESX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( IINFO /= 0 ) THEN
         RESULT( 14 ) = ULPINV
         IF( JTYPE /= 22 ) THEN
            WRITE( NOUNIT, FMT = 9998 )'ZGEESX5', IINFO, N, JTYPE, &
               ISEED
         ELSE
            WRITE( NOUNIT, FMT = 9999 )'ZGEESX5', IINFO, N, &
               ISEED( 1 )
         END IF
         INFO = ABS( IINFO )
         GO TO 220
      END IF
!
!        Perform test (14)
!
      IF( RCNDE1 /= RCONDE ) RESULT( 14 ) = ULPINV
!
!        Perform tests (10), (11), (12), and (13)
!
      IF (ANY(W(1:N) /= WT(1:N))) RESULT( 10 ) = ULPINV
      IF (ANY(H(1:N,1:N) /= HT(1:N,1:N))) RESULT( 11 ) = ULPINV
      IF (ANY(VS(1:N,1:N) /= VS1(1:N,1:N))) RESULT( 12 ) = ULPINV
      IF( SDIM /= SDIM1 ) RESULT( 13 ) = ULPINV
!
!        Compute RCONDE without VS, and compare
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLACPY( 'F', N, N, A, LDA, HT, LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZGEESX( 'N', SORT, ZSLECT, 'E', N, HT, LDA, SDIM1, WT, &
                   VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, &
                   BWORK, IINFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZGEESX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( IINFO /= 0 ) THEN
         RESULT( 14 ) = ULPINV
         IF( JTYPE /= 22 ) THEN
            WRITE( NOUNIT, FMT = 9998 )'ZGEESX6', IINFO, N, JTYPE, ISEED
         ELSE
            WRITE( NOUNIT, FMT = 9999 )'ZGEESX6', IINFO, N, ISEED( 1 )
         END IF
         INFO = ABS( IINFO )
         GO TO 220
      END IF
!
!        Perform test (14)
!
      IF( RCNDE1 /= RCONDE ) RESULT( 14 ) = ULPINV
!
!        Perform tests (10), (11), (12), and (13)
!
      IF (ANY(W(1:N) /= WT(1:N))) RESULT( 10 ) = ULPINV
      IF (ANY(H(1:N,1:N) /= HT(1:N,1:N))) RESULT( 11 ) = ULPINV
      IF (ANY(VS(1:N,1:N) /= VS1(1:N,1:N))) RESULT( 12 ) = ULPINV
      IF( SDIM /= SDIM1 ) RESULT( 13 ) = ULPINV
!
!        Compute RCONDV with VS, and compare
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLACPY( 'F', N, N, A, LDA, HT, LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZGEESX( 'V', SORT, ZSLECT, 'V', N, HT, LDA, SDIM1, WT, &
                   VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, &
                   BWORK, IINFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZGEESX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( IINFO /= 0 ) THEN
         RESULT( 15 ) = ULPINV
         IF( JTYPE /= 22 ) THEN
            WRITE( NOUNIT, FMT = 9998 )'ZGEESX7', IINFO, N, JTYPE, ISEED
         ELSE
            WRITE( NOUNIT, FMT = 9999 )'ZGEESX7', IINFO, N, ISEED( 1 )
         END IF
         INFO = ABS( IINFO )
         GO TO 220
      END IF
!
!        Perform test (15)
!
      IF( RCNDV1 /= RCONDV ) RESULT( 15 ) = ULPINV
!
!        Perform tests (10), (11), (12), and (13)
!
      IF (ANY(W(1:N) /= WT(1:N))) RESULT( 10 ) = ULPINV
      IF (ANY(H(1:N,1:N) /= HT(1:N,1:N))) RESULT( 11 ) = ULPINV
      IF (ANY(VS(1:N,1:N) /= VS1(1:N,1:N))) RESULT( 12 ) = ULPINV
      IF( SDIM /= SDIM1 ) RESULT( 13 ) = ULPINV
!
!        Compute RCONDV without VS, and compare
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLACPY( 'F', N, N, A, LDA, HT, LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZGEESX( 'N', SORT, ZSLECT, 'V', N, HT, LDA, SDIM1, WT, &
                   VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, &
                   BWORK, IINFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZGEESX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( IINFO /= 0 ) THEN
         RESULT( 15 ) = ULPINV
         IF( JTYPE /= 22 ) THEN
            WRITE( NOUNIT, FMT = 9998 )'ZGEESX8', IINFO, N, JTYPE, ISEED
         ELSE
            WRITE( NOUNIT, FMT = 9999 )'ZGEESX8', IINFO, N, ISEED( 1 )
         END IF
         INFO = ABS( IINFO )
         GO TO 220
      END IF
!
!        Perform test (15)
!
      IF( RCNDV1 /= RCONDV ) RESULT( 15 ) = ULPINV
!
!        Perform tests (10), (11), (12), and (13)
!
      IF (ANY(W(1:N) /= WT(1:N))) RESULT( 10 ) = ULPINV
      IF (ANY(H(1:N,1:N) /= HT(1:N,1:N))) RESULT( 11 ) = ULPINV
      IF (ANY(VS(1:N,1:N) /= VS1(1:N,1:N))) RESULT( 12 ) = ULPINV
      IF( SDIM /= SDIM1 ) RESULT( 13 ) = ULPINV
!
   END IF
!
  220 CONTINUE
!
!     If there are precomputed reciprocal condition numbers, compare
!     computed values with them.
!
   IF( COMP ) THEN
!
!        First set up SELOPT, SELDIM, SELVAL, SELWR and SELWI so that
!        the logical function ZSLECT selects the eigenvalues specified
!        by NSLCT, ISLCT and ISRT.
!
      SELDIM = N
      SELOPT = 1
      EPS = MAX( ULP, EPSIN )
      DO I = 1, N
         IPNT( I ) = I
      ENDDO
      SELVAL(1:N) = .FALSE.
      SELWR(1:N) = DBLE( WTMP(1:N) )
      SELWI(1:N) = DIMAG( WTMP(1:N) )
      DO I = 1, N - 1
         KMIN = I
         IF( ISRT == 0 ) THEN
            VRIMIN = DBLE( WTMP( I ) )
         ELSE
            VRIMIN = DIMAG( WTMP( I ) )
         END IF
         DO J = I + 1, N
            IF( ISRT == 0 ) THEN
               VRICMP = DBLE( WTMP( J ) )
            ELSE
               VRICMP = DIMAG( WTMP( J ) )
            END IF
            IF( VRICMP < VRIMIN ) THEN
               KMIN = J
               VRIMIN = VRICMP
            END IF
            ENDDO
         CTMP = WTMP( KMIN )
         WTMP( KMIN ) = WTMP( I )
         WTMP( I ) = CTMP
         ITMP = IPNT( I )
         IPNT( I ) = IPNT( KMIN )
         IPNT( KMIN ) = ITMP
      ENDDO
      SELVAL( IPNT( ISLCT( 1:NSLCT ) ) ) = .TRUE.
!
!        Compute condition numbers
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLACPY( 'F', N, N, A, LDA, HT, LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZGEESX( 'N', 'S', ZSLECT, 'B', N, HT, LDA, SDIM1, WT, VS1, &
                   LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, &
                   IINFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZGEESX : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( IINFO /= 0 ) THEN
         RESULT( 16 ) = ULPINV
         RESULT( 17 ) = ULPINV
         WRITE( NOUNIT, FMT = 9999 )'ZGEESX9', IINFO, N, ISEED( 1 )
         INFO = ABS( IINFO )
         GO TO 270
      END IF
!
!        Compare condition number for average of selected eigenvalues
!        taking its condition number into account
!
      ANORM = ZLANGE( '1', N, N, A, LDA, RWORK )
      V = MAX( DBLE( N )*EPS*ANORM, SMLNUM )
      IF( ANORM == 0.0D0 ) &
         V = 1.0D0
      IF( V > RCONDV ) THEN
         TOL = 1.0D0
      ELSE
         TOL = V / RCONDV
      END IF
      IF( V > RCDVIN ) THEN
         TOLIN = 1.0D0
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
         RESULT( 16 ) = 1.0D0
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
         RESULT( 17 ) = 1.0D0
      END IF
!
  270    CONTINUE
!
   END IF
!
 9999 FORMAT( ' ZGET24: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', &
         I6, ', INPUT EXAMPLE NUMBER = ', I4 )
 9998 FORMAT( ' ZGET24: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', &
         I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
!
   RETURN
!
!     End of ZGET24
!
END




