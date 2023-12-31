!> \brief \b DQRT17
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DQRT17( TRANS, IRESID, M, N, NRHS, A,
!                        LDA, X, LDX, B, LDB, C, WORK, LWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            IRESID, LDA, LDB, LDX, LWORK, M, N, NRHS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDB, * ),
!      $                   WORK( LWORK ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DQRT17 computes the ratio
!>
!>    norm(R**T * op(A)) / ( norm(A) * alpha * max(M,N,NRHS) * EPS ),
!>
!> where R = B - op(A)*X, op(A) is A or A**T, depending on TRANS, EPS
!> is the machine epsilon, and
!>
!>    alpha = norm(B) if IRESID = 1 (zero-residual problem)
!>    alpha = norm(R) if IRESID = 2 (otherwise).
!>
!> The norm used is the 1-norm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies whether or not the transpose of A is used.
!>          = 'N':  No transpose, op(A) = A.
!>          = 'T':  Transpose, op(A) = A**T.
!> \endverbatim
!>
!> \param[in] IRESID
!> \verbatim
!>          IRESID is INTEGER
!>          IRESID = 1 indicates zero-residual problem.
!>          IRESID = 2 indicates non-zero residual.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!>          If TRANS = 'N', the number of rows of the matrix B.
!>          If TRANS = 'T', the number of rows of the matrix X.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix  A.
!>          If TRANS = 'N', the number of rows of the matrix X.
!>          If TRANS = 'T', the number of rows of the matrix B.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns of the matrices X and B.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m-by-n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= M.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)
!>          If TRANS = 'N', the n-by-nrhs matrix X.
!>          If TRANS = 'T', the m-by-nrhs matrix X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.
!>          If TRANS = 'N', LDX >= N.
!>          If TRANS = 'T', LDX >= M.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!>          If TRANS = 'N', the m-by-nrhs matrix B.
!>          If TRANS = 'T', the n-by-nrhs matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.
!>          If TRANS = 'N', LDB >= M.
!>          If TRANS = 'T', LDB >= N.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDB,NRHS)
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
!>          The length of the array WORK.  LWORK >= NRHS*(M+N).
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
!> \ingroup double_lin
!
!  =====================================================================
   DOUBLE PRECISION FUNCTION DQRT17( TRANS, IRESID, M, N, NRHS, A, &
                    LDA, X, LDX, B, LDB, C, WORK, LWORK )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          TRANS
   INTEGER            IRESID, LDA, LDB, LDX, LWORK, M, N, NRHS
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDB, * ), &
                      WORK( LWORK ), X( LDX, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
   INTEGER            INFO, ISCL, NCOLS, NROWS
   DOUBLE PRECISION   ERR, NORMA, NORMB, NORMRS, SMLNUM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   RWORK( 1 )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   DOUBLE PRECISION   DLAMCH, DLANGE
   EXTERNAL           LSAME, DLAMCH, DLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           DGEMM, DLACPY, DLASCL, XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          DBLE, MAX
!     ..
!     .. Executable Statements ..
!
   DQRT17 = ZERO
!
   IF( LSAME( TRANS, 'N' ) ) THEN
      NROWS = M
      NCOLS = N
   ELSE IF( LSAME( TRANS, 'T' ) ) THEN
      NROWS = N
      NCOLS = M
   ELSE
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'DQRT17', 1 )
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
   IF( LWORK < NCOLS*NRHS ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'DQRT17', 13 )
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
   IF( M <= 0 .OR. N <= 0 .OR. NRHS <= 0 ) THEN
      RETURN
   END IF
!
   NORMA = DLANGE( 'One-norm', M, N, A, LDA, RWORK )
   SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
   ISCL = 0
!
!     compute residual and scale it
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLACPY( 'All', NROWS, NRHS, B, LDB, C, LDB )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DGEMM( TRANS, 'No transpose', NROWS, NRHS, NCOLS, -ONE, A, &
               LDA, X, LDX, ONE, C, LDB )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   NORMRS = DLANGE( 'Max', NROWS, NRHS, C, LDB, RWORK )
   IF( NORMRS > SMLNUM ) THEN
      ISCL = 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLASCL( 'General', 0, 0, NORMRS, ONE, NROWS, NRHS, C, LDB, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLASCL : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   END IF
!
!     compute R**T * op(A)
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DGEMM( 'Transpose', TRANS, NRHS, NCOLS, NROWS, ONE, C, LDB, &
               A, LDA, ZERO, WORK, NRHS )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     compute and properly scale error
!
   ERR = DLANGE( 'One-norm', NRHS, NCOLS, WORK, NRHS, RWORK )
   IF( NORMA /= ZERO ) &
      ERR = ERR / NORMA
!
   IF( ISCL == 1 ) &
      ERR = ERR*NORMRS
!
   IF( IRESID == 1 ) THEN
      NORMB = DLANGE( 'One-norm', NROWS, NRHS, B, LDB, RWORK )
      IF( NORMB /= ZERO ) &
         ERR = ERR / NORMB
   ELSE
      IF( NORMRS /= ZERO ) &
         ERR = ERR / NORMRS
   END IF
!
   DQRT17 = ERR / ( DLAMCH( 'Epsilon' )*DBLE( MAX( M, N, NRHS ) ) )
   RETURN
!
!     End of DQRT17
!
END
                                                                                                                                                                                                                                                                                                            




