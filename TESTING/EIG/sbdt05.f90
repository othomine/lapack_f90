!> \brief \b SBDT05
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SBDT05( M, N, A, LDA, S, NS, U, LDU,
!                          VT, LDVT, WORK, RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDU, LDVT, N, NS
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * ), S( * ), U( LDU, * ),
!      $                   VT( LDVT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SBDT05 reconstructs a bidiagonal matrix B from its (partial) SVD:
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
!>          A is REAL array, dimension (LDA,N)
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
!>          S is REAL array, dimension (NS)
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
!>          U is REAL array, dimension (LDU,NS)
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
!>          VT is REAL array, dimension (LDVT,N)
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
!>          WORK is REAL array, dimension (M,N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
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
   SUBROUTINE SBDT05( M, N, A, LDA, S, NS, U, LDU, &
                       VT, LDVT, WORK, RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, LDU, LDVT, M, N, NS
   REAL               RESID
!     ..
!     .. Array Arguments ..
   REAL               A( LDA, * ), S( * ), U( LDU, * ), &
                      VT( LDVT, * ), WORK( * )
!     ..
!
! ======================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, J
   REAL               ANORM, EPS
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ISAMAX
   REAL               SASUM, SLAMCH, SLANGE
   EXTERNAL           LSAME, ISAMAX, SASUM, SLAMCH, SLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           SGEMM
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible.
!
   RESID = 0.0E+0
   IF( MIN( M, N ) <= 0 .OR. NS <= 0 ) RETURN
!
   EPS = SLAMCH( 'Precision' )
   ANORM = SLANGE( 'M', M, N, A, LDA, WORK )
!
!     Compute U' * A * V.
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SGEMM( 'N', 'T', M, NS, N, 1.0E+0, A, LDA, VT, &
               LDVT, 0.0E+0, WORK( 1+NS*NS ), M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SGEMM( 'T', 'N', NS, NS, M, -1.0E+0, U, LDU, WORK( 1+NS*NS ), &
               M, 0.0E+0, WORK, NS )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     norm(S - U' * B * V)
!
   J = 0
   DO I = 1, NS
      WORK( J+I ) =  WORK( J+I ) + S( I )
      RESID = MAX( RESID, SASUM( NS, WORK( J+1 ), 1 ) )
      J = J + NS
   ENDDO
!
   IF( ANORM <= 0.0E+0 ) THEN
      IF( RESID /= 0.0E+0 ) &
         RESID = 1.0E+0 / EPS
   ELSE
      IF( ANORM >= RESID ) THEN
         RESID = ( RESID / ANORM ) / ( REAL( N )*EPS )
      ELSE
         IF( ANORM < 1.0E+0 ) THEN
            RESID = ( MIN( RESID, REAL( N )*ANORM ) / ANORM ) / &
                    ( REAL( N )*EPS )
         ELSE
            RESID = MIN( RESID / ANORM, REAL( N ) ) / &
                    ( REAL( N )*EPS )
         END IF
      END IF
   END IF
!
   RETURN
!
!     End of SBDT05
!
END



