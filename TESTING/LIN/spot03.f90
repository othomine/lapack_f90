!> \brief \b SPOT03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SPOT03( UPLO, N, A, LDA, AINV, LDAINV, WORK, LDWORK,
!                          RWORK, RCOND, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDAINV, LDWORK, N
!       REAL               RCOND, RESID
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), AINV( LDAINV, * ), RWORK( * ),
!      $                   WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SPOT03 computes the residual for a symmetric matrix times its
!> inverse:
!>    norm( I - A*AINV ) / ( N * norm(A) * norm(AINV) * EPS ),
!> where EPS is the machine epsilon.
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
!>          A is REAL array, dimension (LDA,N)
!>          The original symmetric matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N)
!> \endverbatim
!>
!> \param[in,out] AINV
!> \verbatim
!>          AINV is REAL array, dimension (LDAINV,N)
!>          On entry, the inverse of the matrix A, stored as a symmetric
!>          matrix in the same format as A.
!>          In this version, AINV is expanded into a full matrix and
!>          multiplied by A, so the opposing triangle of AINV will be
!>          changed; i.e., if the upper triangular part of AINV is
!>          stored, the lower triangular part will be used as work space.
!> \endverbatim
!>
!> \param[in] LDAINV
!> \verbatim
!>          LDAINV is INTEGER
!>          The leading dimension of the array AINV.  LDAINV >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (LDWORK,N)
!> \endverbatim
!>
!> \param[in] LDWORK
!> \verbatim
!>          LDWORK is INTEGER
!>          The leading dimension of the array WORK.  LDWORK >= max(1,N).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is REAL
!>          The reciprocal of the condition number of A, computed as
!>          ( 1/norm(A) ) / norm(AINV).
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
!>          norm(I - A*AINV) / ( N * norm(A) * norm(AINV) * EPS )
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
   SUBROUTINE SPOT03( UPLO, N, A, LDA, AINV, LDAINV, WORK, LDWORK, &
                      RWORK, RCOND, RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            LDA, LDAINV, LDWORK, N
   REAL               RCOND, RESID
!     ..
!     .. Array Arguments ..
   REAL               A( LDA, * ), AINV( LDAINV, * ), RWORK( * ), &
                      WORK( LDWORK, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO, ONE
   PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, J
   REAL               AINVNM, ANORM, EPS
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               SLAMCH, SLANGE, SLANSY
   EXTERNAL           LSAME, SLAMCH, SLANGE, SLANSY
!     ..
!     .. External Subroutines ..
   EXTERNAL           SSYMM
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          REAL
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0.
!
   IF( N <= 0 ) THEN
      RCOND = ONE
      RESID = ZERO
      RETURN
   END IF
!
!     Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.
!
   EPS = SLAMCH( 'Epsilon' )
   ANORM = SLANSY( '1', UPLO, N, A, LDA, RWORK )
   AINVNM = SLANSY( '1', UPLO, N, AINV, LDAINV, RWORK )
   IF( ANORM <= ZERO .OR. AINVNM <= ZERO ) THEN
      RCOND = ZERO
      RESID = ONE / EPS
      RETURN
   END IF
   RCOND = ( ONE / ANORM ) / AINVNM
!
!     Expand AINV into a full matrix and call SSYMM to multiply
!     AINV on the left by A.
!
   IF( LSAME( UPLO, 'U' ) ) THEN
      DO J = 1, N
         DO I = 1, J - 1
            AINV( J, I ) = AINV( I, J )
         ENDDO
      ENDDO
   ELSE
      DO J = 1, N
         DO I = J + 1, N
            AINV( J, I ) = AINV( I, J )
         ENDDO
      ENDDO
   END IF
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SSYMM( 'Left', UPLO, N, N, -ONE, A, LDA, AINV, LDAINV, ZERO, &
               WORK, LDWORK )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SSYMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Add the identity matrix to WORK .
!
   DO I = 1, N
      WORK( I, I ) = WORK( I, I ) + ONE
   ENDDO
!
!     Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)
!
   RESID = SLANGE( '1', N, N, WORK, LDWORK, RWORK )
!
   RESID = ( ( RESID*RCOND ) / EPS ) / REAL( N )
!
   RETURN
!
!     End of SPOT03
!
END
                                                                                                                                                                                                                                                                                                            




