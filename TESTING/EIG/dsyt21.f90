!> \brief \b DSYT21
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYT21( ITYPE, UPLO, N, KBAND, A, LDA, D, E, U, LDU, V,
!                          LDV, TAU, WORK, RESULT )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            ITYPE, KBAND, LDA, LDU, LDV, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), RESULT( 2 ),
!      $                   TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSYT21 generally checks a decomposition of the form
!>
!>    A = U S U**T
!>
!> where **T means transpose, A is symmetric, U is orthogonal, and S is
!> diagonal (if KBAND=0) or symmetric tridiagonal (if KBAND=1).
!>
!> If ITYPE=1, then U is represented as a dense matrix; otherwise U is
!> expressed as a product of Householder transformations, whose vectors
!> are stored in the array "V" and whose scaling constants are in "TAU".
!> We shall use the letter "V" to refer to the product of Householder
!> transformations (which should be equal to U).
!>
!> Specifically, if ITYPE=1, then:
!>
!>    RESULT(1) = | A - U S U**T | / ( |A| n ulp ) and
!>    RESULT(2) = | I - U U**T | / ( n ulp )
!>
!> If ITYPE=2, then:
!>
!>    RESULT(1) = | A - V S V**T | / ( |A| n ulp )
!>
!> If ITYPE=3, then:
!>
!>    RESULT(1) = | I - V U**T | / ( n ulp )
!>
!> For ITYPE > 1, the transformation U is expressed as a product
!> V = H(1)...H(n-2),  where H(j) = I  -  tau(j) v(j) v(j)**T and each
!> vector v(j) has its first j elements 0 and the remaining n-j elements
!> stored in V(j+1:n,j).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          Specifies the type of tests to be performed.
!>          1: U expressed as a dense orthogonal matrix:
!>             RESULT(1) = | A - U S U**T | / ( |A| n ulp )  and
!>             RESULT(2) = | I - U U**T | / ( n ulp )
!>
!>          2: U expressed as a product V of Housholder transformations:
!>             RESULT(1) = | A - V S V**T | / ( |A| n ulp )
!>
!>          3: U expressed both as a dense orthogonal matrix and
!>             as a product of Housholder transformations:
!>             RESULT(1) = | I - V U**T | / ( n ulp )
!> \endverbatim
!>
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
!>          The size of the matrix.  If it is zero, DSYT21 does nothing.
!>          It must be at least zero.
!> \endverbatim
!>
!> \param[in] KBAND
!> \verbatim
!>          KBAND is INTEGER
!>          The bandwidth of the matrix.  It may only be zero or one.
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
!>          and at least N.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The diagonal of the (symmetric tri-) diagonal matrix.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          The off-diagonal of the (symmetric tri-) diagonal matrix.
!>          E(1) is the (1,2) and (2,1) element, E(2) is the (2,3) and
!>          (3,2) element, etc.
!>          Not referenced if KBAND=0.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is DOUBLE PRECISION array, dimension (LDU, N)
!>          If ITYPE=1 or 3, this contains the orthogonal matrix in
!>          the decomposition, expressed as a dense matrix.  If ITYPE=2,
!>          then it is not referenced.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U.  LDU must be at least N and
!>          at least 1.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension (LDV, N)
!>          If ITYPE=2 or 3, the columns of this array contain the
!>          Householder vectors used to describe the orthogonal matrix
!>          in the decomposition.  If UPLO='L', then the vectors are in
!>          the lower triangle, if UPLO='U', then in the upper
!>          triangle.
!>          *NOTE* If ITYPE=2 or 3, V is modified and restored.  The
!>          subdiagonal (if UPLO='L') or the superdiagonal (if UPLO='U')
!>          is set to one, and later reset to its original value, during
!>          the course of the calculation.
!>          If ITYPE=1, then it is neither referenced nor modified.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of V.  LDV must be at least N and
!>          at least 1.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (N)
!>          If ITYPE >= 2, then TAU(j) is the scalar factor of
!>          v(j) v(j)**T in the Householder transformation H(j) of
!>          the product  U = H(1)...H(n-2)
!>          If ITYPE < 2, then TAU is not referenced.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (2*N**2)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (2)
!>          The values computed by the two tests described above.  The
!>          values are currently limited to 1/ulp, to avoid overflow.
!>          RESULT(1) is always modified.  RESULT(2) is modified only
!>          if ITYPE=1.
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
   SUBROUTINE DSYT21( ITYPE, UPLO, N, KBAND, A, LDA, D, E, U, LDU, V, &
                      LDV, TAU, WORK, RESULT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            ITYPE, KBAND, LDA, LDU, LDV, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), RESULT( 2 ), &
                      TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            LOWER
   CHARACTER          CUPLO
   INTEGER            IINFO, J, JCOL, JR, JROW
   DOUBLE PRECISION   ANORM, ULP, UNFL, VSAVE, WNORM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   DOUBLE PRECISION   DLAMCH, DLANGE, DLANSY
   EXTERNAL           LSAME, DLAMCH, DLANGE, DLANSY
!     ..
!     .. External Subroutines ..
   EXTERNAL           DGEMM, DLACPY, DLARFY, DLASET, DORM2L, DORM2R, &
                      DSYR, DSYR2
!     ..
!     .. Executable Statements ..
!
   RESULT( 1 ) = 0.0D+0
   IF( ITYPE == 1 ) RESULT( 2 ) = 0.0D+0
   IF( N <= 0 ) RETURN
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
   IF( ITYPE < 1 .OR. ITYPE > 3 ) THEN
      RESULT( 1 ) = 10.0D0 / ULP
      RETURN
   END IF
!
!     Do Test 1
!
!     Norm of A:
!
   IF( ITYPE == 3 ) THEN
      ANORM = 1.0D0
   ELSE
      ANORM = MAX( DLANSY( '1', CUPLO, N, A, LDA, WORK ), UNFL )
   END IF
!
!     Compute error matrix:
!
   IF( ITYPE == 1 ) THEN
!
!        ITYPE=1: error = A - U S U**T
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLASET( 'Full', N, N, 0.0D+0, 0.0D+0, WORK, N )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLACPY( CUPLO, N, N, A, LDA, WORK, N )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
      DO J = 1, N
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DSYR( CUPLO, N, -D( J ), U( 1, J ), 1, WORK, N )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DSYR : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
!
      IF( N > 1 .AND. KBAND == 1 ) THEN
         DO J = 1, N - 1
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DSYR2( CUPLO, N, -E( J ), U( 1, J ), 1, U( 1, J+1 ), &
                        1, WORK, N )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DSYR2 : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ENDDO
      END IF
      WNORM = DLANSY( '1', CUPLO, N, WORK, N, WORK( N**2+1 ) )
!
   ELSE IF( ITYPE == 2 ) THEN
!
!        ITYPE=2: error = V S V**T - A
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLASET( 'Full', N, N, 0.0D+0, 0.0D+0, WORK, N )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
      IF( LOWER ) THEN
         WORK( N**2 ) = D( N )
         DO J = N - 1, 1, -1
            IF( KBAND == 1 ) THEN
               WORK( ( N+1 )*( J-1 )+2 ) = ( 1.0D0-TAU( J ) )*E( J )
               DO JR = J + 2, N
                  WORK( ( J-1 )*N+JR ) = -TAU( J )*E( J )*V( JR, J )
               ENDDO
            END IF
!
            VSAVE = V( J+1, J )
            V( J+1, J ) = 1.0D0
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLARFY( 'L', N-J, V( J+1, J ), 1, TAU( J ), &
                         WORK( ( N+1 )*J+1 ), N, WORK( N**2+1 ) )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLARFY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            V( J+1, J ) = VSAVE
            WORK( ( N+1 )*( J-1 )+1 ) = D( J )
         ENDDO
      ELSE
         WORK( 1 ) = D( 1 )
         DO J = 1, N - 1
            IF( KBAND == 1 ) THEN
               WORK( ( N+1 )*J ) = ( 1.0D0-TAU( J ) )*E( J )
               DO JR = 1, J - 1
                  WORK( J*N+JR ) = -TAU( J )*E( J )*V( JR, J+1 )
               ENDDO
            END IF
!
            VSAVE = V( J, J+1 )
            V( J, J+1 ) = 1.0D0
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLARFY( 'U', J, V( 1, J+1 ), 1, TAU( J ), WORK, N, &
                         WORK( N**2+1 ) )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLARFY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            V( J, J+1 ) = VSAVE
            WORK( ( N+1 )*J+1 ) = D( J+1 )
         ENDDO
      END IF
!
      DO JCOL = 1, N
         IF( LOWER ) THEN
            DO JROW = JCOL, N
               WORK( JROW+N*( JCOL-1 ) ) = WORK( JROW+N*( JCOL-1 ) ) &
                   - A( JROW, JCOL )
            ENDDO
         ELSE
            DO JROW = 1, JCOL
               WORK( JROW+N*( JCOL-1 ) ) = WORK( JROW+N*( JCOL-1 ) ) &
                   - A( JROW, JCOL )
            ENDDO
         END IF
      ENDDO
      WNORM = DLANSY( '1', CUPLO, N, WORK, N, WORK( N**2+1 ) )
!
   ELSE IF( ITYPE == 3 ) THEN
!
!        ITYPE=3: error = U V**T - I
!
      IF( N < 2 ) &
         RETURN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLACPY( ' ', N, N, U, LDU, WORK, N )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( LOWER ) THEN
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DORM2R( 'R', 'T', N, N-1, N-1, V( 2, 1 ), LDV, TAU, &
                      WORK( N+1 ), N, WORK( N**2+1 ), IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DORM2R : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ELSE
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DORM2L( 'R', 'T', N, N-1, N-1, V( 1, 2 ), LDV, TAU, &
                      WORK, N, WORK( N**2+1 ), IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DORM2L : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      END IF
      IF( IINFO /= 0 ) THEN
         RESULT( 1 ) = 10.0D0 / ULP
         RETURN
      END IF
!
      DO J = 1, N
         WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - 1.0D0
         ENDDO
!
      WNORM = DLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) )
   END IF
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
   IF( ITYPE == 1 ) THEN
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
      RESULT( 2 ) = MIN( DLANGE( '1', N, N, WORK, N, &
                    WORK( N**2+1 ) ), DBLE( N ) ) / ( N*ULP )
   END IF
!
   RETURN
!
!     End of DSYT21
!
END




