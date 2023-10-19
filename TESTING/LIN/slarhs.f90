!> \brief \b SLARHS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLARHS( PATH, XTYPE, UPLO, TRANS, M, N, KL, KU, NRHS,
!                          A, LDA, X, LDX, B, LDB, ISEED, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS, UPLO, XTYPE
!       CHARACTER*3        PATH
!       INTEGER            INFO, KL, KU, LDA, LDB, LDX, M, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       REAL               A( LDA, * ), B( LDB, * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLARHS chooses a set of NRHS random solution vectors and sets
!> up the right hand sides for the linear system
!>    op(A) * X = B,
!> where op(A) = A or A**T, depending on TRANS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] PATH
!> \verbatim
!>          PATH is CHARACTER*3
!>          The type of the real matrix A.  PATH may be given in any
!>          combination of upper and lower case.  Valid types include
!>             xGE:  General m x n matrix
!>             xGB:  General banded matrix
!>             xPO:  Symmetric positive definite, 2-D storage
!>             xPP:  Symmetric positive definite packed
!>             xPB:  Symmetric positive definite banded
!>             xSY:  Symmetric indefinite, 2-D storage
!>             xSP:  Symmetric indefinite packed
!>             xSB:  Symmetric indefinite banded
!>             xTR:  Triangular
!>             xTP:  Triangular packed
!>             xTB:  Triangular banded
!>             xQR:  General m x n matrix
!>             xLQ:  General m x n matrix
!>             xQL:  General m x n matrix
!>             xRQ:  General m x n matrix
!>          where the leading character indicates the precision.
!> \endverbatim
!>
!> \param[in] XTYPE
!> \verbatim
!>          XTYPE is CHARACTER*1
!>          Specifies how the exact solution X will be determined:
!>          = 'N':  New solution; generate a random X.
!>          = 'C':  Computed; use value of X on entry.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          matrix A is stored, if A is symmetric.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Used only if A is nonsymmetric; specifies the operation
!>          applied to the matrix A.
!>          = 'N':  B := A    * X  (No transpose)
!>          = 'T':  B := A**T * X  (Transpose)
!>          = 'C':  B := A**H * X  (Conjugate transpose = Transpose)
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number or rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          Used only if A is a band matrix; specifies the number of
!>          subdiagonals of A if A is a general band matrix or if A is
!>          symmetric or triangular and UPLO = 'L'; specifies the number
!>          of superdiagonals of A if A is symmetric or triangular and
!>          UPLO = 'U'.  0 <= KL <= M-1.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          Used only if A is a general band matrix or if A is
!>          triangular.
!>
!>          If PATH = xGB, specifies the number of superdiagonals of A,
!>          and 0 <= KU <= N-1.
!>
!>          If PATH = xTR, xTP, or xTB, specifies whether or not the
!>          matrix has unit diagonal:
!>          = 1:  matrix has non-unit diagonal (default)
!>          = 2:  matrix has unit diagonal
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand side vectors in the system A*X = B.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The test matrix whose type is given by PATH.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          If PATH = xGB, LDA >= KL+KU+1.
!>          If PATH = xPB, xSB, xHB, or xTB, LDA >= KL+1.
!>          Otherwise, LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is or output) REAL array, dimension(LDX,NRHS)
!>          On entry, if XTYPE = 'C' (for 'Computed'), then X contains
!>          the exact solution to the system of linear equations.
!>          On exit, if XTYPE = 'N' (for 'New'), then X is initialized
!>          with random values.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  If TRANS = 'N',
!>          LDX >= max(1,N); if TRANS = 'T', LDX >= max(1,M).
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,NRHS)
!>          The right hand side vector(s) for the system of equations,
!>          computed from B = op(A) * X, where op(A) is determined by
!>          TRANS.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  If TRANS = 'N',
!>          LDB >= max(1,M); if TRANS = 'T', LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          The seed vector for the random number generator (used in
!>          SLATMS).  Modified on exit.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup single_lin
!
!  =====================================================================
   SUBROUTINE SLARHS( PATH, XTYPE, UPLO, TRANS, M, N, KL, KU, NRHS, &
                      A, LDA, X, LDX, B, LDB, ISEED, INFO )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          TRANS, UPLO, XTYPE
   CHARACTER*3        PATH
   INTEGER            INFO, KL, KU, LDA, LDB, LDX, M, N, NRHS
!     ..
!     .. Array Arguments ..
   INTEGER            ISEED( 4 )
   REAL               A( LDA, * ), B( LDB, * ), X( LDX, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ONE, ZERO
   PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            BAND, GEN, NOTRAN, QRS, SYM, TRAN, TRI
   CHARACTER          C1, DIAG
   CHARACTER*2        C2
   INTEGER            J, MB, NX
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME, LSAMEN
   EXTERNAL           LSAME, LSAMEN
!     ..
!     .. External Subroutines ..
   EXTERNAL           SGBMV, SGEMM, SLACPY, SLARNV, SSBMV, SSPMV, &
                      SSYMM, STBMV, STPMV, STRMM, XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   C1 = PATH( 1: 1 )
   C2 = PATH( 2: 3 )
   TRAN = LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' )
   NOTRAN = .NOT.TRAN
   GEN = LSAME( PATH( 2: 2 ), 'G' )
   QRS = LSAME( PATH( 2: 2 ), 'Q' ) .OR. LSAME( PATH( 3: 3 ), 'Q' )
   SYM = LSAME( PATH( 2: 2 ), 'P' ) .OR. LSAME( PATH( 2: 2 ), 'S' )
   TRI = LSAME( PATH( 2: 2 ), 'T' )
   BAND = LSAME( PATH( 3: 3 ), 'B' )
   IF( .NOT.LSAME( C1, 'Single precision' ) ) THEN
      INFO = -1
   ELSE IF( .NOT.( LSAME( XTYPE, 'N' ) .OR. LSAME( XTYPE, 'C' ) ) ) &
             THEN
      INFO = -2
   ELSE IF( ( SYM .OR. TRI ) .AND. .NOT. &
            ( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) ) ) THEN
      INFO = -3
   ELSE IF( ( GEN .OR. QRS ) .AND. .NOT. &
            ( TRAN .OR. LSAME( TRANS, 'N' ) ) ) THEN
      INFO = -4
   ELSE IF( M < 0 ) THEN
      INFO = -5
   ELSE IF( N < 0 ) THEN
      INFO = -6
   ELSE IF( BAND .AND. KL < 0 ) THEN
      INFO = -7
   ELSE IF( BAND .AND. KU < 0 ) THEN
      INFO = -8
   ELSE IF( NRHS < 0 ) THEN
      INFO = -9
   ELSE IF( ( .NOT.BAND .AND. LDA < MAX( 1, M ) ) .OR. &
            ( BAND .AND. ( SYM .OR. TRI ) .AND. LDA < KL+1 ) .OR. &
            ( BAND .AND. GEN .AND. LDA < KL+KU+1 ) ) THEN
      INFO = -11
   ELSE IF( ( NOTRAN .AND. LDX < MAX( 1, N ) ) .OR. &
            ( TRAN .AND. LDX < MAX( 1, M ) ) ) THEN
      INFO = -13
   ELSE IF( ( NOTRAN .AND. LDB < MAX( 1, M ) ) .OR. &
            ( TRAN .AND. LDB < MAX( 1, N ) ) ) THEN
      INFO = -15
   END IF
   IF( INFO /= 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'SLARHS', -INFO )
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
!     Initialize X to NRHS random vectors unless XTYPE = 'C'.
!
   IF( TRAN ) THEN
      NX = M
      MB = N
   ELSE
      NX = N
      MB = M
   END IF
   IF( .NOT.LSAME( XTYPE, 'C' ) ) THEN
      DO J = 1, NRHS
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SLARNV( 2, ISEED, N, X( 1, J ) )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SLARNV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
   END IF
!
!     Multiply X by op(A) using an appropriate
!     matrix multiply routine.
!
   IF( LSAMEN( 2, C2, 'GE' ) .OR. LSAMEN( 2, C2, 'QR' ) .OR. &
       LSAMEN( 2, C2, 'LQ' ) .OR. LSAMEN( 2, C2, 'QL' ) .OR. &
       LSAMEN( 2, C2, 'RQ' ) ) THEN
!
!        General matrix
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL SGEMM( TRANS, 'N', MB, NRHS, NX, ONE, A, LDA, X, LDX, &
                  ZERO, B, LDB )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
   ELSE IF( LSAMEN( 2, C2, 'PO' ) .OR. LSAMEN( 2, C2, 'SY' ) ) THEN
!
!        Symmetric matrix, 2-D storage
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL SSYMM( 'Left', UPLO, N, NRHS, ONE, A, LDA, X, LDX, ZERO, &
                  B, LDB )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : SSYMM : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
   ELSE IF( LSAMEN( 2, C2, 'GB' ) ) THEN
!
!        General matrix, band storage
!
      DO J = 1, NRHS
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SGBMV( TRANS, MB, NX, KL, KU, ONE, A, LDA, X( 1, J ), &
                     1, ZERO, B( 1, J ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SGBMV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
!
   ELSE IF( LSAMEN( 2, C2, 'PB' ) ) THEN
!
!        Symmetric matrix, band storage
!
      DO J = 1, NRHS
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SSBMV( UPLO, N, KL, ONE, A, LDA, X( 1, J ), 1, ZERO, &
                     B( 1, J ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SSBMV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
!
   ELSE IF( LSAMEN( 2, C2, 'PP' ) .OR. LSAMEN( 2, C2, 'SP' ) ) THEN
!
!        Symmetric matrix, packed storage
!
      DO J = 1, NRHS
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SSPMV( UPLO, N, ONE, A, X( 1, J ), 1, ZERO, B( 1, J ), &
                     1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SSPMV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
!
   ELSE IF( LSAMEN( 2, C2, 'TR' ) ) THEN
!
!        Triangular matrix.  Note that for triangular matrices,
!           KU = 1 => non-unit triangular
!           KU = 2 => unit triangular
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL SLACPY( 'Full', N, NRHS, X, LDX, B, LDB )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( KU == 2 ) THEN
         DIAG = 'U'
      ELSE
         DIAG = 'N'
      END IF
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL STRMM( 'Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B, &
                  LDB )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : STRMM : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
   ELSE IF( LSAMEN( 2, C2, 'TP' ) ) THEN
!
!        Triangular matrix, packed storage
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL SLACPY( 'Full', N, NRHS, X, LDX, B, LDB )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( KU == 2 ) THEN
         DIAG = 'U'
      ELSE
         DIAG = 'N'
      END IF
      DO J = 1, NRHS
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL STPMV( UPLO, TRANS, DIAG, N, A, B( 1, J ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : STPMV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
!
   ELSE IF( LSAMEN( 2, C2, 'TB' ) ) THEN
!
!        Triangular matrix, banded storage
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL SLACPY( 'Full', N, NRHS, X, LDX, B, LDB )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( KU == 2 ) THEN
         DIAG = 'U'
      ELSE
         DIAG = 'N'
      END IF
      DO J = 1, NRHS
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL STBMV( UPLO, TRANS, DIAG, N, KL, A, LDA, B( 1, J ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : STBMV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
!
   ELSE
!
!        If PATH is none of the above, return with an error code.
!
      INFO = -1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'SLARHS', -INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : XERBLA : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   END IF
!
   RETURN
!
!     End of SLARHS
!
END
                                                                                                                                                                                                                                                                                                            




