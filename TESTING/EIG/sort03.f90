!> \brief \b SORT03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SORT03( RC, MU, MV, N, K, U, LDU, V, LDV, WORK, LWORK,
!                          RESULT, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER*( * )    RC
!       INTEGER            INFO, K, LDU, LDV, LWORK, MU, MV, N
!       REAL               RESULT
!       ..
!       .. Array Arguments ..
!       REAL               U( LDU, * ), V( LDV, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SORT03 compares two orthogonal matrices U and V to see if their
!> corresponding rows or columns span the same spaces.  The rows are
!> checked if RC = 'R', and the columns are checked if RC = 'C'.
!>
!> RESULT is the maximum of
!>
!>    | V*V' - I | / ( MV ulp ), if RC = 'R', or
!>
!>    | V'*V - I | / ( MV ulp ), if RC = 'C',
!>
!> and the maximum over rows (or columns) 1 to K of
!>
!>    | U(i) - S*V(i) |/ ( N ulp )
!>
!> where S is +-1 (chosen to minimize the expression), U(i) is the i-th
!> row (column) of U, and V(i) is the i-th row (column) of V.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] RC
!> \verbatim
!>          RC is CHARACTER*1
!>          If RC = 'R' the rows of U and V are to be compared.
!>          If RC = 'C' the columns of U and V are to be compared.
!> \endverbatim
!>
!> \param[in] MU
!> \verbatim
!>          MU is INTEGER
!>          The number of rows of U if RC = 'R', and the number of
!>          columns if RC = 'C'.  If MU = 0 SORT03 does nothing.
!>          MU must be at least zero.
!> \endverbatim
!>
!> \param[in] MV
!> \verbatim
!>          MV is INTEGER
!>          The number of rows of V if RC = 'R', and the number of
!>          columns if RC = 'C'.  If MV = 0 SORT03 does nothing.
!>          MV must be at least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          If RC = 'R', the number of columns in the matrices U and V,
!>          and if RC = 'C', the number of rows in U and V.  If N = 0
!>          SORT03 does nothing.  N must be at least zero.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of rows or columns of U and V to compare.
!>          0 <= K <= max(MU,MV).
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is REAL array, dimension (LDU,N)
!>          The first matrix to compare.  If RC = 'R', U is MU by N, and
!>          if RC = 'C', U is N by MU.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U.  If RC = 'R', LDU >= max(1,MU),
!>          and if RC = 'C', LDU >= max(1,N).
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is REAL array, dimension (LDV,N)
!>          The second matrix to compare.  If RC = 'R', V is MV by N, and
!>          if RC = 'C', V is N by MV.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of V.  If RC = 'R', LDV >= max(1,MV),
!>          and if RC = 'C', LDV >= max(1,N).
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
!>          The length of the array WORK.  For best performance, LWORK
!>          should be at least N*N if RC = 'C' or M*M if RC = 'R', but
!>          the tests will be done even if LWORK is 0.
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL
!>          The value computed by the test described above.  RESULT is
!>          limited to 1/ulp to avoid overflow.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          0  indicates a successful exit
!>          -k indicates the k-th parameter had an illegal value
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
!> \ingroup single_eig
!
!  =====================================================================
   SUBROUTINE SORT03( RC, MU, MV, N, K, U, LDU, V, LDV, WORK, LWORK, &
                      RESULT, INFO )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER*( * )    RC
   INTEGER            INFO, K, LDU, LDV, LWORK, MU, MV, N
   REAL               RESULT
!     ..
!     .. Array Arguments ..
   REAL               U( LDU, * ), V( LDV, * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, IRC, J, LMX
   REAL               RES1, RES2, S, ULP
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ISAMAX
   REAL               SLAMCH
   EXTERNAL           LSAME, ISAMAX, SLAMCH
!     .. External Subroutines ..
   EXTERNAL           SORT01, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Check inputs
!
   INFO = 0
   IF( LSAME( RC, 'R' ) ) THEN
      IRC = 0
   ELSE IF( LSAME( RC, 'C' ) ) THEN
      IRC = 1
   ELSE
      IRC = -1
   END IF
   IF( IRC == -1 ) THEN
      INFO = -1
   ELSE IF( MU < 0 ) THEN
      INFO = -2
   ELSE IF( MV < 0 ) THEN
      INFO = -3
   ELSE IF( N < 0 ) THEN
      INFO = -4
   ELSE IF( K < 0 .OR. K > MAX( MU, MV ) ) THEN
      INFO = -5
   ELSE IF( ( IRC == 0 .AND. LDU < MAX( 1, MU ) ) .OR. &
            ( IRC == 1 .AND. LDU < MAX( 1, N ) ) ) THEN
      INFO = -7
   ELSE IF( ( IRC == 0 .AND. LDV < MAX( 1, MV ) ) .OR. &
            ( IRC == 1 .AND. LDV < MAX( 1, N ) ) ) THEN
      INFO = -9
   END IF
   IF( INFO /= 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'SORT03', -INFO )
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
!     Initialize result
!
   RESULT = 0.0E+0
   IF( MU == 0 .OR. MV == 0 .OR. N == 0 ) RETURN
!
!     Machine constants
!
   ULP = SLAMCH( 'Precision' )
!
   IF( IRC == 0 ) THEN
!
!        Compare rows
!
      RES1 = 0.0E+0
      DO I = 1, K
         LMX = ISAMAX( N, U( I, 1 ), LDU )
         S = SIGN( 1.0E+0, U( I, LMX ) )*SIGN( 1.0E+0, V( I, LMX ) )
         DO J = 1, N
            RES1 = MAX( RES1, ABS( U( I, J )-S*V( I, J ) ) )
         ENDDO
      ENDDO
      RES1 = RES1 / ( REAL( N )*ULP )
!
!        Compute orthogonality of rows of V.
!
      CALL SORT01( 'Rows', MV, N, V, LDV, WORK, LWORK, RES2 )
!
   ELSE
!
!        Compare columns
!
      RES1 = 0.0E+0
      DO I = 1, K
         LMX = ISAMAX( N, U( 1, I ), 1 )
         S = SIGN( 1.0E+0, U( LMX, I ) )*SIGN( 1.0E+0, V( LMX, I ) )
         DO J = 1, N
            RES1 = MAX( RES1, ABS( U( J, I )-S*V( J, I ) ) )
         ENDDO
      ENDDO
      RES1 = RES1 / ( REAL( N )*ULP )
!
!        Compute orthogonality of columns of V.
!
      CALL SORT01( 'Columns', N, MV, V, LDV, WORK, LWORK, RES2 )
   END IF
!
   RESULT = MIN( MAX( RES1, RES2 ), 1.0E+0 / ULP )
   RETURN
!
!     End of SORT03
!
END




