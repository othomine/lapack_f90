!> \brief \b DSPT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSPT01( UPLO, N, A, AFAC, IPIV, C, LDC, RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDC, N
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( * ), AFAC( * ), C( LDC, * ), RWORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSPT01 reconstructs a symmetric indefinite packed matrix A from its
!> block L*D*L' or U*D*U' factorization and computes the residual
!>      norm( C - A ) / ( N * norm(A) * EPS ),
!> where C is the reconstructed matrix and EPS is the machine epsilon.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix A is stored:
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows and columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          The original symmetric matrix A, stored as a packed
!>          triangular matrix.
!> \endverbatim
!>
!> \param[in] AFAC
!> \verbatim
!>          AFAC is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          The factored form of the matrix A, stored as a packed
!>          triangular matrix.  AFAC contains the block diagonal matrix D
!>          and the multipliers used to obtain the factor L or U from the
!>          block L*D*L' or U*D*U' factorization as computed by DSPTRF.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices from DSPTRF.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C.  LDC >= max(1,N).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          If UPLO = 'L', norm(L*D*L' - A) / ( N * norm(A) * EPS )
!>          If UPLO = 'U', norm(U*D*U' - A) / ( N * norm(A) * EPS )
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
   SUBROUTINE DSPT01( UPLO, N, A, AFAC, IPIV, C, LDC, RWORK, RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            LDC, N
   DOUBLE PRECISION   RESID
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * )
   DOUBLE PRECISION   A( * ), AFAC( * ), C( LDC, * ), RWORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, INFO, J, JC
   DOUBLE PRECISION   ANORM, EPS
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   DOUBLE PRECISION   DLAMCH, DLANSP, DLANSY
   EXTERNAL           LSAME, DLAMCH, DLANSP, DLANSY
!     ..
!     .. External Subroutines ..
   EXTERNAL           DLASET, DLAVSP
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          DBLE
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0.
!
   IF( N <= 0 ) THEN
      RESID = ZERO
      RETURN
   END IF
!
!     Determine EPS and the norm of A.
!
   EPS = DLAMCH( 'Epsilon' )
   ANORM = DLANSP( '1', UPLO, N, A, RWORK )
!
!     Initialize C to the identity matrix.
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLASET( 'Full', N, N, ZERO, ONE, C, LDC )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Call DLAVSP to form the product D * U' (or D * L' ).
!
   CALL DLAVSP( UPLO, 'Transpose', 'Non-unit', N, N, AFAC, IPIV, C, &
                LDC, INFO )
!
!     Call DLAVSP again to multiply by U ( or L ).
!
   CALL DLAVSP( UPLO, 'No transpose', 'Unit', N, N, AFAC, IPIV, C, &
                LDC, INFO )
!
!     Compute the difference  C - A .
!
   IF( LSAME( UPLO, 'U' ) ) THEN
      JC = 0
      DO J = 1, N
         DO I = 1, J
            C( I, J ) = C( I, J ) - A( JC+I )
         ENDDO
         JC = JC + J
      ENDDO
   ELSE
      JC = 1
      DO J = 1, N
         DO I = J, N
            C( I, J ) = C( I, J ) - A( JC+I-J )
         ENDDO
         JC = JC + N - J + 1
      ENDDO
   END IF
!
!     Compute norm( C - A ) / ( N * norm(A) * EPS )
!
   RESID = DLANSY( '1', UPLO, N, C, LDC, RWORK )
!
   IF( ANORM <= ZERO ) THEN
      IF( RESID /= ZERO ) &
         RESID = ONE / EPS
   ELSE
      RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS
   END IF
!
   RETURN
!
!     End of DSPT01
!
END
                                                                                                                                                                                                                                                                                                            




