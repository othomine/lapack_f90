!> \brief \b DSBT21
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSBT21( UPLO, N, KA, KS, A, LDA, D, E, U, LDU, WORK,
!                          RESULT )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            KA, KS, LDA, LDU, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), RESULT( 2 ),
!      $                   U( LDU, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSBT21  generally checks a decomposition of the form
!>
!>         A = U S U**T
!>
!> where **T means transpose, A is symmetric banded, U is
!> orthogonal, and S is diagonal (if KS=0) or symmetric
!> tridiagonal (if KS=1).
!>
!> Specifically:
!>
!>         RESULT(1) = | A - U S U**T | / ( |A| n ulp ) and
!>         RESULT(2) = | I - U U**T | / ( n ulp )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER
!>          If UPLO='U', the upper triangle of A and V will be used and
!>          the (strictly) lower triangle will not be referenced.
!>          If UPLO='L', the lower triangle of A and V will be used and
!>          the (strictly) upper triangle will not be referenced.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The size of the matrix.  If it is zero, DSBT21 does nothing.
!>          It must be at least zero.
!> \endverbatim
!>
!> \param[in] KA
!> \verbatim
!>          KA is INTEGER
!>          The bandwidth of the matrix A.  It must be at least zero.  If
!>          it is larger than N-1, then max( 0, N-1 ) will be used.
!> \endverbatim
!>
!> \param[in] KS
!> \verbatim
!>          KS is INTEGER
!>          The bandwidth of the matrix S.  It may only be zero or one.
!>          If zero, then S is diagonal, and E is not referenced.  If
!>          one, then S is symmetric tri-diagonal.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA, N)
!>          The original (unfactored) matrix.  It is assumed to be
!>          symmetric, and only the upper (UPLO='U') or only the lower
!>          (UPLO='L') will be referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  It must be at least 1
!>          and at least min( KA, N-1 ).
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The diagonal of the (symmetric tri-) diagonal matrix S.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          The off-diagonal of the (symmetric tri-) diagonal matrix S.
!>          E(1) is the (1,2) and (2,1) element, E(2) is the (2,3) and
!>          (3,2) element, etc.
!>          Not referenced if KS=0.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is DOUBLE PRECISION array, dimension (LDU, N)
!>          The orthogonal matrix in the decomposition, expressed as a
!>          dense matrix (i.e., not as a product of Householder
!>          transformations, Givens transformations, etc.)
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U.  LDU must be at least N and
!>          at least 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N**2+N)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (2)
!>          The values computed by the two tests described above.  The
!>          values are currently limited to 1/ulp, to avoid overflow.
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
   SUBROUTINE DSBT21( UPLO, N, KA, KS, A, LDA, D, E, U, LDU, WORK, &
                      RESULT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            KA, KS, LDA, LDU, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), RESULT( 2 ), &
                      U( LDU, * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            LOWER
   CHARACTER          CUPLO
   INTEGER            IKA, J, JC, JR, LW
   DOUBLE PRECISION   ANORM, ULP, UNFL, WNORM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   DOUBLE PRECISION   DLAMCH, DLANGE, DLANSB, DLANSP
   EXTERNAL           LSAME, DLAMCH, DLANGE, DLANSB, DLANSP
!     ..
!     .. External Subroutines ..
   EXTERNAL           DGEMM, DSPR, DSPR2
!     ..
!     .. Executable Statements ..
!
!     Constants
!
   RESULT( 1:2 ) = 0.0D+0
   IF( N <= 0 ) RETURN
!
   IKA = MAX( 0, MIN( N-1, KA ) )
   LW = ( N*( N+1 ) ) / 2
!
   IF( LSAME( UPLO, 'U' ) ) THEN
      LOWER = .FALSE.
      CUPLO = 'U'
   ELSE
      LOWER = .TRUE.
      CUPLO = 'L'
   END IF
!
   UNFL = DLAMCH( 'Safe minimum' )
   ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
!
!     Some Error Checks
!
!     Do Test 1
!
!     Norm of A:
!
   ANORM = MAX( DLANSB( '1', CUPLO, N, IKA, A, LDA, WORK ), UNFL )
!
!     Compute error matrix:    Error = A - U S U**T
!
!     Copy A from SB to SP storage format.
!
   J = 0
   DO JC = 1, N
      IF( LOWER ) THEN
         DO JR = 1, MIN( IKA+1, N+1-JC )
            J = J + 1
            WORK( J ) = A( JR, JC )
         ENDDO
         DO JR = IKA + 2, N + 1 - JC
            J = J + 1
            WORK( J ) = 0.0D+0
         ENDDO
      ELSE
         DO JR = IKA + 2, JC
            J = J + 1
            WORK( J ) = 0.0D+0
         ENDDO
         DO JR = MIN( IKA, JC-1 ), 0, -1
            J = J + 1
            WORK( J ) = A( IKA+1-JR, JC )
         ENDDO
      END IF
   ENDDO
!
   DO J = 1, N
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DSPR( CUPLO, N, -D( J ), U( 1, J ), 1, WORK )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DSPR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   ENDDO
!
   IF( N > 1 .AND. KS == 1 ) THEN
      DO J = 1, N - 1
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DSPR2( CUPLO, N, -E( J ), U( 1, J ), 1, U( 1, J+1 ), 1, &
                     WORK )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DSPR2 : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
   END IF
   WNORM = DLANSP( '1', CUPLO, N, WORK, WORK( LW+1 ) )
!
   IF( ANORM > WNORM ) THEN
      RESULT( 1 ) = ( WNORM / ANORM ) / ( N*ULP )
   ELSE
      IF( ANORM < 1.0D0 ) THEN
         RESULT( 1 ) = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
      ELSE
         RESULT( 1 ) = MIN( WNORM / ANORM, DBLE( N ) ) / ( N*ULP )
      END IF
   END IF
!
!     Do Test 2
!
!     Compute  U U**T - I
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DGEMM( 'N', 'C', N, N, N, 1.0D0, U, LDU, U, LDU, 0.0D+0, WORK, &
               N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   DO J = 1, N
      WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - 1.0D0
   ENDDO
!
   RESULT( 2 ) = MIN( DLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) ), &
                 DBLE( N ) ) / ( N*ULP )
!
   RETURN
!
!     End of DSBT21
!
END




