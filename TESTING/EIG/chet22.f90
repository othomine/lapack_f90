!> \brief \b CHET22
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHET22( ITYPE, UPLO, N, M, KBAND, A, LDA, D, E, U, LDU,
!                          V, LDV, TAU, WORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            ITYPE, KBAND, LDA, LDU, LDV, M, N
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * ), RESULT( 2 ), RWORK( * )
!       COMPLEX            A( LDA, * ), TAU( * ), U( LDU, * ),
!      $                   V( LDV, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      CHET22  generally checks a decomposition of the form
!>
!>              A U = U S
!>
!>      where A is complex Hermitian, the columns of U are orthonormal,
!>      and S is diagonal (if KBAND=0) or symmetric tridiagonal (if
!>      KBAND=1).  If ITYPE=1, then U is represented as a dense matrix,
!>      otherwise the U is expressed as a product of Householder
!>      transformations, whose vectors are stored in the array "V" and
!>      whose scaling constants are in "TAU"; we shall use the letter
!>      "V" to refer to the product of Householder transformations
!>      (which should be equal to U).
!>
!>      Specifically, if ITYPE=1, then:
!>
!>              RESULT(1) = | U**H A U - S | / ( |A| m ulp ) and
!>              RESULT(2) = | I - U**H U | / ( m ulp )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \verbatim
!>  ITYPE   INTEGER
!>          Specifies the type of tests to be performed.
!>          1: U expressed as a dense orthogonal matrix:
!>             RESULT(1) = | A - U S U**H | / ( |A| n ulp )  and
!>             RESULT(2) = | I - U U**H | / ( n ulp )
!>
!>  UPLO    CHARACTER
!>          If UPLO='U', the upper triangle of A will be used and the
!>          (strictly) lower triangle will not be referenced.  If
!>          UPLO='L', the lower triangle of A will be used and the
!>          (strictly) upper triangle will not be referenced.
!>          Not modified.
!>
!>  N       INTEGER
!>          The size of the matrix.  If it is zero, CHET22 does nothing.
!>          It must be at least zero.
!>          Not modified.
!>
!>  M       INTEGER
!>          The number of columns of U.  If it is zero, CHET22 does
!>          nothing.  It must be at least zero.
!>          Not modified.
!>
!>  KBAND   INTEGER
!>          The bandwidth of the matrix.  It may only be zero or one.
!>          If zero, then S is diagonal, and E is not referenced.  If
!>          one, then S is symmetric tri-diagonal.
!>          Not modified.
!>
!>  A       COMPLEX array, dimension (LDA , N)
!>          The original (unfactored) matrix.  It is assumed to be
!>          symmetric, and only the upper (UPLO='U') or only the lower
!>          (UPLO='L') will be referenced.
!>          Not modified.
!>
!>  LDA     INTEGER
!>          The leading dimension of A.  It must be at least 1
!>          and at least N.
!>          Not modified.
!>
!>  D       REAL array, dimension (N)
!>          The diagonal of the (symmetric tri-) diagonal matrix.
!>          Not modified.
!>
!>  E       REAL array, dimension (N)
!>          The off-diagonal of the (symmetric tri-) diagonal matrix.
!>          E(1) is ignored, E(2) is the (1,2) and (2,1) element, etc.
!>          Not referenced if KBAND=0.
!>          Not modified.
!>
!>  U       COMPLEX array, dimension (LDU, N)
!>          If ITYPE=1, this contains the orthogonal matrix in
!>          the decomposition, expressed as a dense matrix.
!>          Not modified.
!>
!>  LDU     INTEGER
!>          The leading dimension of U.  LDU must be at least N and
!>          at least 1.
!>          Not modified.
!>
!>  V       COMPLEX array, dimension (LDV, N)
!>          If ITYPE=2 or 3, the lower triangle of this array contains
!>          the Householder vectors used to describe the orthogonal
!>          matrix in the decomposition.  If ITYPE=1, then it is not
!>          referenced.
!>          Not modified.
!>
!>  LDV     INTEGER
!>          The leading dimension of V.  LDV must be at least N and
!>          at least 1.
!>          Not modified.
!>
!>  TAU     COMPLEX array, dimension (N)
!>          If ITYPE >= 2, then TAU(j) is the scalar factor of
!>          v(j) v(j)**H in the Householder transformation H(j) of
!>          the product  U = H(1)...H(n-2)
!>          If ITYPE < 2, then TAU is not referenced.
!>          Not modified.
!>
!>  WORK    COMPLEX array, dimension (2*N**2)
!>          Workspace.
!>          Modified.
!>
!>  RWORK   REAL array, dimension (N)
!>          Workspace.
!>          Modified.
!>
!>  RESULT  REAL array, dimension (2)
!>          The values computed by the two tests described above.  The
!>          values are currently limited to 1/ulp, to avoid overflow.
!>          RESULT(1) is always modified.  RESULT(2) is modified only
!>          if LDU is at least N.
!>          Modified.
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
!> \ingroup complex_eig
!
!  =====================================================================
   SUBROUTINE CHET22( ITYPE, UPLO, N, M, KBAND, A, LDA, D, E, U, LDU, &
                      V, LDV, TAU, WORK, RWORK, RESULT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            ITYPE, KBAND, LDA, LDU, LDV, M, N
!     ..
!     .. Array Arguments ..
   REAL               D( * ), E( * ), RESULT( 2 ), RWORK( * )
   COMPLEX            A( LDA, * ), TAU( * ), U( LDU, * ), &
                      V( LDV, * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            J, JJ, JJ1, JJ2, NN, NNP1
   REAL               ANORM, ULP, UNFL, WNORM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   REAL               CLANHE, SLAMCH
   EXTERNAL           CLANHE, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMM, CHEMM
!     ..
!     .. Executable Statements ..
!
   RESULT( 1:2 ) = 0.0E0
   IF( N <= 0 .OR. M <= 0 ) RETURN
!
   UNFL = SLAMCH( 'Safe minimum' )
   ULP = SLAMCH( 'Precision' )
!
!     Do Test 1
!
!     Norm of A:
!
   ANORM = MAX( CLANHE( '1', UPLO, N, A, LDA, RWORK ), UNFL )
!
!     Compute error matrix:
!
!     ITYPE=1: error = U**H A U - S
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CHEMM( 'L', UPLO, N, M, (1.0E0,0.0E0), A, LDA, U, LDU, (0.0E0,0.0E0), WORK, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CHEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   NN = N*N
   NNP1 = NN + 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CGEMM( 'C', 'N', M, M, N, (1.0E0,0.0E0), U, LDU, WORK, N, (0.0E0,0.0E0), WORK( NNP1 ), N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   DO J = 1, M
      JJ = NN + ( J-1 )*N + J
      WORK( JJ ) = WORK( JJ ) - D( J )
   ENDDO
   IF( KBAND == 1 .AND. N > 1 ) THEN
      DO J = 2, M
         JJ1 = NN + ( J-1 )*N + J - 1
         JJ2 = NN + ( J-2 )*N + J
         WORK( JJ1 ) = WORK( JJ1 ) - E( J-1 )
         WORK( JJ2 ) = WORK( JJ2 ) - E( J-1 )
      ENDDO
   END IF
   WNORM = CLANHE( '1', UPLO, M, WORK( NNP1 ), N, RWORK )
!
   IF( ANORM > WNORM ) THEN
      RESULT( 1 ) = ( WNORM / ANORM ) / ( M*ULP )
   ELSE
      IF( ANORM < 1.0E0 ) THEN
         RESULT( 1 ) = ( MIN( WNORM, M*ANORM ) / ANORM ) / ( M*ULP )
      ELSE
         RESULT( 1 ) = MIN( WNORM / ANORM, REAL( M ) ) / ( M*ULP )
      END IF
   END IF
!
!     Do Test 2
!
!     Compute  U**H U - I
!
   IF( ITYPE == 1 ) &
      CALL CUNT01( 'Columns', N, M, U, LDU, WORK, 2*N*N, RWORK, RESULT( 2 ) )
!
   RETURN
!
!     End of CHET22
!
END




