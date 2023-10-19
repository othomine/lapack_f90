!> \brief \b ZBDT05
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZBDT05( M, N, A, LDA, S, NS, U, LDU,
!                          VT, LDVT, WORK, RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDU, LDVT, N, NS
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!      DOUBLE PRECISION   S( * )
!      COMPLEX*16         A( LDA, * ), U( * ), VT( LDVT, * ), WORK( * )
!       ..
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZBDT05 reconstructs a bidiagonal matrix B from its (partial) SVD:
!>    S = U' * B * V
!> where U and V are orthogonal matrices and S is diagonal.
!>
!> The test ratio to test the singular value decomposition is
!>    RESID = norm( S - U' * B * V ) / ( n * norm(B) * EPS )
!> where VT = V' and EPS is the machine precision.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrices A and U.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrices A and VT.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (NS)
!>          The singular values from the (partial) SVD of B, sorted in
!>          decreasing order.
!> \endverbatim
!>
!> \param[in] NS
!> \verbatim
!>          NS is INTEGER
!>          The number of singular values/vectors from the (partial)
!>          SVD of B.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is COMPLEX*16 array, dimension (LDU,NS)
!>          The n by ns orthogonal matrix U in S = U' * B * V.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= max(1,N)
!> \endverbatim
!>
!> \param[in] VT
!> \verbatim
!>          VT is COMPLEX*16 array, dimension (LDVT,N)
!>          The n by ns orthogonal matrix V in S = U' * B * V.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>          The leading dimension of the array VT.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (M,N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          The test ratio:  norm(S - U' * A * V) / ( n * norm(A) * EPS )
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
!> \ingroup double_eig
!
!  =====================================================================
   SUBROUTINE ZBDT05( M, N, A, LDA, S, NS, U, LDU, &
                       VT, LDVT, WORK, RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, LDU, LDVT, M, N, NS
   DOUBLE PRECISION   RESID
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   S( * )
   COMPLEX*16         A( LDA, * ), U( * ), VT( LDVT, * ), WORK( * )
!     ..
!
! ======================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, J
   DOUBLE PRECISION   ANORM, EPS
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   DUM( 1 )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            IDAMAX
   DOUBLE PRECISION   DASUM, DZASUM, DLAMCH, ZLANGE
   EXTERNAL           LSAME, IDAMAX, DASUM, DZASUM, DLAMCH, ZLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           ZGEMM
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible.
!
   RESID = 0.0D+0
   IF( MIN( M, N ) <= 0 .OR. NS <= 0 ) RETURN
!
   EPS = DLAMCH( 'Precision' )
   ANORM = ZLANGE( 'M', M, N, A, LDA, DUM )
!
!     Compute U' * A * V.
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'N', 'C', M, NS, N, (1.0D+0,0.0D+0), A, LDA, VT, &
               LDVT, (0.0D+0,0.0D+0), WORK( 1+NS*NS ), M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'C', 'N', NS, NS, M, -(1.0D+0,0.0D+0), U, LDU, WORK( 1+NS*NS ), &
               M, (0.0D+0,0.0D+0), WORK, NS )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     norm(S - U' * B * V)
!
   J = 0
   DO I = 1, NS
      WORK( J+I ) =  WORK( J+I ) + DCMPLX( S( I ), 0.0D+0 )
      RESID = MAX( RESID, DZASUM( NS, WORK( J+1 ), 1 ) )
      J = J + NS
   ENDDO
!
   IF( ANORM <= 0.0D+0 ) THEN
      IF( RESID /= 0.0D+0 ) &
         RESID = 1.0D+0 / EPS
   ELSE
      IF( ANORM >= RESID ) THEN
         RESID = ( RESID / ANORM ) / ( DBLE( N )*EPS )
      ELSE
         IF( ANORM < 1.0D+0 ) THEN
            RESID = ( MIN( RESID, DBLE( N )*ANORM ) / ANORM ) / &
                    ( DBLE( N )*EPS )
         ELSE
            RESID = MIN( RESID / ANORM, DBLE( N ) ) / &
                    ( DBLE( N )*EPS )
         END IF
      END IF
   END IF
!
   RETURN
!
!     End of ZBDT05
!
END




