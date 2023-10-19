!> \brief \b DBDT02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DBDT02( M, N, B, LDB, C, LDC, U, LDU, WORK, RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            LDB, LDC, LDU, M, N
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   B( LDB, * ), C( LDC, * ), U( LDU, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DBDT02 tests the change of basis C = U**H * B by computing the
!> residual
!>
!>    RESID = norm(B - U * C) / ( max(m,n) * norm(B) * EPS ),
!>
!> where B and C are M by N matrices, U is an M by M orthogonal matrix,
!> and EPS is the machine precision.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrices B and C and the order of
!>          the matrix Q.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrices B and C.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>          The m by n matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,M).
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!>          The m by n matrix C, assumed to contain U**H * B.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C.  LDC >= max(1,M).
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is DOUBLE PRECISION array, dimension (LDU,M)
!>          The m by m orthogonal matrix U.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          RESID = norm(B - U * C) / ( max(m,n) * norm(B) * EPS ),
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
   SUBROUTINE DBDT02( M, N, B, LDB, C, LDC, U, LDU, WORK, RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDB, LDC, LDU, M, N
   DOUBLE PRECISION   RESID
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   B( LDB, * ), C( LDC, * ), U( LDU, * ), &
                      WORK( * )
!     ..
!
! ======================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            J
   DOUBLE PRECISION   BNORM, EPS, REALMN
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DASUM, DLAMCH, DLANGE
   EXTERNAL           DASUM, DLAMCH, DLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           DCOPY, DGEMV
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
   RESID = 0.0D+0
   IF( M <= 0 .OR. N <= 0 ) RETURN
   REALMN = DBLE( MAX( M, N ) )
   EPS = DLAMCH( 'Precision' )
!
!     Compute norm(B - U * C)
!
   DO J = 1, N
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( M, B( 1, J ), 1, WORK, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGEMV( 'No transpose', M, M, -1.0D+0, U, LDU, C( 1, J ), 1, 1.0D+0, WORK, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGEMV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      RESID = MAX( RESID, DASUM( M, WORK, 1 ) )
   ENDDO
!
!     Compute norm of B.
!
   BNORM = DLANGE( '1', M, N, B, LDB, WORK )
!
   IF( BNORM <= 0.0D+0 ) THEN
      IF( RESID /= 0.0D+0 ) RESID = 1.0D+0 / EPS
   ELSE
      IF( BNORM >= RESID ) THEN
         RESID = ( RESID / BNORM ) / ( REALMN*EPS )
      ELSE
         IF( BNORM < 1.0D+0 ) THEN
            RESID = ( MIN( RESID, REALMN*BNORM ) / BNORM ) / ( REALMN*EPS )
         ELSE
            RESID = MIN( RESID / BNORM, REALMN ) / ( REALMN*EPS )
         END IF
      END IF
   END IF
   RETURN
!
!     End of DBDT02
!
END




