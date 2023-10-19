!> \brief \b ZUNT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNT01( ROWCOL, M, N, U, LDU, WORK, LWORK, RWORK,
!                          RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          ROWCOL
!       INTEGER            LDU, LWORK, M, N
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         U( LDU, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNT01 checks that the matrix U is unitary by computing the ratio
!>
!>    RESID = norm( I - U*U' ) / ( n * EPS ), if ROWCOL = 'R',
!> or
!>    RESID = norm( I - U'*U ) / ( m * EPS ), if ROWCOL = 'C'.
!>
!> Alternatively, if there isn't sufficient workspace to form
!> I - U*U' or I - U'*U, the ratio is computed as
!>
!>    RESID = abs( I - U*U' ) / ( n * EPS ), if ROWCOL = 'R',
!> or
!>    RESID = abs( I - U'*U ) / ( m * EPS ), if ROWCOL = 'C'.
!>
!> where EPS is the machine precision.  ROWCOL is used only if m = n;
!> if m > n, ROWCOL is assumed to be 'C', and if m < n, ROWCOL is
!> assumed to be 'R'.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ROWCOL
!> \verbatim
!>          ROWCOL is CHARACTER
!>          Specifies whether the rows or columns of U should be checked
!>          for orthogonality.  Used only if M = N.
!>          = 'R':  Check for orthogonal rows of U
!>          = 'C':  Check for orthogonal columns of U
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix U.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix U.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is COMPLEX*16 array, dimension (LDU,N)
!>          The unitary matrix U.  U is checked for orthogonal columns
!>          if m > n or if m = n and ROWCOL = 'C'.  U is checked for
!>          orthogonal rows if m < n or if m = n and ROWCOL = 'R'.
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
!>          WORK is COMPLEX*16 array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  For best performance, LWORK
!>          should be at least N*N if ROWCOL = 'C' or M*M if
!>          ROWCOL = 'R', but the test will be done even if LWORK is 0.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (min(M,N))
!>          Used only if LWORK is large enough to use the Level 3 BLAS
!>          code.
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          RESID = norm( I - U * U' ) / ( n * EPS ), if ROWCOL = 'R', or
!>          RESID = norm( I - U' * U ) / ( m * EPS ), if ROWCOL = 'C'.
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
   SUBROUTINE ZUNT01( ROWCOL, M, N, U, LDU, WORK, LWORK, RWORK, &
                      RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          ROWCOL
   INTEGER            LDU, LWORK, M, N
   DOUBLE PRECISION   RESID
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   RWORK( * )
   COMPLEX*16         U( LDU, * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   CHARACTER          TRANSU
   INTEGER            I, J, K, LDWORK, MNMIN
   DOUBLE PRECISION   EPS
   COMPLEX*16         TMP, ZDUM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   DOUBLE PRECISION   DLAMCH, ZLANSY
   COMPLEX*16         ZDOTC
   EXTERNAL           LSAME, DLAMCH, ZLANSY, ZDOTC
!     ..
!     .. External Subroutines ..
   EXTERNAL           ZHERK, ZLASET
!     ..
!     .. Statement Functions ..
   DOUBLE PRECISION   CABS1
!     ..
!     .. Statement Function definitions ..
   CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
!     ..
!     .. Executable Statements ..
!
   RESID = 0.0D0
!
!     Quick return if possible
!
   IF( M <= 0 .OR. N <= 0 ) RETURN
!
   EPS = DLAMCH( 'Precision' )
   IF( M < N .OR. ( M == N .AND. LSAME( ROWCOL, 'R' ) ) ) THEN
      TRANSU = 'N'
      K = N
   ELSE
      TRANSU = 'C'
      K = M
   END IF
   MNMIN = MIN( M, N )
!
   IF( ( MNMIN+1 )*MNMIN <= LWORK ) THEN
      LDWORK = MNMIN
   ELSE
      LDWORK = 0
   END IF
   IF( LDWORK > 0 ) THEN
!
!        Compute I - U*U' or I - U'*U.
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLASET( 'Upper', MNMIN, MNMIN, DCMPLX( 0.0D0 ), &
                   DCMPLX( 1.0D0 ), WORK, LDWORK )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHERK( 'Upper', TRANSU, MNMIN, K, -1.0D0, U, LDU, 1.0D0, WORK, &
                  LDWORK )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHERK : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        Compute norm( I - U*U' ) / ( K * EPS ) .
!
      RESID = ZLANSY( '1', 'Upper', MNMIN, WORK, LDWORK, RWORK )
      RESID = ( RESID / DBLE( K ) ) / EPS
   ELSE IF( TRANSU == 'C' ) THEN
!
!        Find the maximum element in abs( I - U'*U ) / ( m * EPS )
!
      DO J = 1, N
         DO I = 1, J
            IF( I /= J ) THEN
               TMP = 0.0D0
            ELSE
               TMP = 1.0D0
            END IF
            TMP = TMP - ZDOTC( M, U( 1, I ), 1, U( 1, J ), 1 )
            RESID = MAX( RESID, CABS1( TMP ) )
         ENDDO
      ENDDO
      RESID = ( RESID / DBLE( M ) ) / EPS
   ELSE
!
!        Find the maximum element in abs( I - U*U' ) / ( n * EPS )
!
      DO J = 1, M
         DO I = 1, J
            IF( I /= J ) THEN
               TMP = 0.0D0
            ELSE
               TMP = 1.0D0
            END IF
            TMP = TMP - ZDOTC( N, U( J, 1 ), LDU, U( I, 1 ), LDU )
            RESID = MAX( RESID, CABS1( TMP ) )
         ENDDO
      ENDDO
      RESID = ( RESID / DBLE( N ) ) / EPS
   END IF
   RETURN
!
!     End of ZUNT01
!
END




