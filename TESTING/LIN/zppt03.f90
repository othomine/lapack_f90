!> \brief \b ZPPT03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZPPT03( UPLO, N, A, AINV, WORK, LDWORK, RWORK, RCOND,
!                          RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDWORK, N
!       DOUBLE PRECISION   RCOND, RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( * ), AINV( * ), WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZPPT03 computes the residual for a Hermitian packed matrix times its
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
!>          Hermitian matrix A is stored:
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
!>          A is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          The original Hermitian matrix A, stored as a packed
!>          triangular matrix.
!> \endverbatim
!>
!> \param[in] AINV
!> \verbatim
!>          AINV is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          The (Hermitian) inverse of the matrix A, stored as a packed
!>          triangular matrix.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LDWORK,N)
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
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is DOUBLE PRECISION
!>          The reciprocal of the condition number of A, computed as
!>          ( 1/norm(A) ) / norm(AINV).
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
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
!> \ingroup complex16_lin
!
!  =====================================================================
   SUBROUTINE ZPPT03( UPLO, N, A, AINV, WORK, LDWORK, RWORK, RCOND, &
                      RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            LDWORK, N
   DOUBLE PRECISION   RCOND, RESID
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   RWORK( * )
   COMPLEX*16         A( * ), AINV( * ), WORK( LDWORK, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
   COMPLEX*16         CZERO, CONE
   PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), &
                      CONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
   INTEGER            I, J, JJ
   DOUBLE PRECISION   AINVNM, ANORM, EPS
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   DOUBLE PRECISION   DLAMCH, ZLANGE, ZLANHP
   EXTERNAL           LSAME, DLAMCH, ZLANGE, ZLANHP
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          DBLE, DCONJG
!     ..
!     .. External Subroutines ..
   EXTERNAL           ZCOPY, ZHPMV
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
   EPS = DLAMCH( 'Epsilon' )
   ANORM = ZLANHP( '1', UPLO, N, A, RWORK )
   AINVNM = ZLANHP( '1', UPLO, N, AINV, RWORK )
   IF( ANORM <= ZERO .OR. AINVNM <= ZERO ) THEN
      RCOND = ZERO
      RESID = ONE / EPS
      RETURN
   END IF
   RCOND = ( ONE / ANORM ) / AINVNM
!
!     UPLO = 'U':
!     Copy the leading N-1 x N-1 submatrix of AINV to WORK(1:N,2:N) and
!     expand it to a full matrix, then multiply by A one column at a
!     time, moving the result one column to the left.
!
   IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Copy AINV
!
      JJ = 1
      DO J = 1, N - 1
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZCOPY( J, AINV( JJ ), 1, WORK( 1, J+1 ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         DO I = 1, J - 1
            WORK( J, I+1 ) = DCONJG( AINV( JJ+I-1 ) )
         ENDDO
         JJ = JJ + J
      ENDDO
      JJ = ( ( N-1 )*N ) / 2 + 1
      DO I = 1, N - 1
         WORK( N, I+1 ) = DCONJG( AINV( JJ+I-1 ) )
      ENDDO
!
!        Multiply by A
!
      DO J = 1, N - 1
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZHPMV( 'Upper', N, -CONE, A, WORK( 1, J+1 ), 1, CZERO, &
                     WORK( 1, J ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZHPMV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPMV( 'Upper', N, -CONE, A, AINV( JJ ), 1, CZERO, &
                  WORK( 1, N ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPMV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!     UPLO = 'L':
!     Copy the trailing N-1 x N-1 submatrix of AINV to WORK(1:N,1:N-1)
!     and multiply by A, moving each column to the right.
!
   ELSE
!
!        Copy AINV
!
      DO I = 1, N - 1
         WORK( 1, I ) = DCONJG( AINV( I+1 ) )
      ENDDO
      JJ = N + 1
      DO J = 2, N
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZCOPY( N-J+1, AINV( JJ ), 1, WORK( J, J-1 ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         DO I = 1, N - J
            WORK( J, J+I-1 ) = DCONJG( AINV( JJ+I ) )
         ENDDO
         JJ = JJ + N - J + 1
      ENDDO
!
!        Multiply by A
!
      DO J = N, 2, -1
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZHPMV( 'Lower', N, -CONE, A, WORK( 1, J-1 ), 1, CZERO, &
                     WORK( 1, J ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZHPMV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHPMV( 'Lower', N, -CONE, A, AINV( 1 ), 1, CZERO, &
                  WORK( 1, 1 ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHPMV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
   END IF
!
!     Add the identity matrix to WORK .
!
   DO I = 1, N
      WORK( I, I ) = WORK( I, I ) + CONE
   ENDDO
!
!     Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)
!
   RESID = ZLANGE( '1', N, N, WORK, LDWORK, RWORK )
!
   RESID = ( ( RESID*RCOND ) / EPS ) / DBLE( N )
!
   RETURN
!
!     End of ZPPT03
!
END
                                                                                                                                                                                                                                                                                                            




