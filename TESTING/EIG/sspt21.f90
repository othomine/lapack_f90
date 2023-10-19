!> \brief \b SSPT21
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSPT21( ITYPE, UPLO, N, KBAND, AP, D, E, U, LDU, VP,
!                          TAU, WORK, RESULT )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            ITYPE, KBAND, LDU, N
!       ..
!       .. Array Arguments ..
!       REAL               AP( * ), D( * ), E( * ), RESULT( 2 ), TAU( * ),
!      $                   U( LDU, * ), VP( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSPT21  generally checks a decomposition of the form
!>
!>         A = U S U**T
!>
!> where **T means transpose, A is symmetric (stored in packed format), U
!> is orthogonal, and S is diagonal (if KBAND=0) or symmetric
!> tridiagonal (if KBAND=1).  If ITYPE=1, then U is represented as a
!> dense matrix, otherwise the U is expressed as a product of
!> Householder transformations, whose vectors are stored in the array
!> "V" and whose scaling constants are in "TAU"; we shall use the
!> letter "V" to refer to the product of Householder transformations
!> (which should be equal to U).
!>
!> Specifically, if ITYPE=1, then:
!>
!>         RESULT(1) = | A - U S U**T | / ( |A| n ulp ) and
!>         RESULT(2) = | I - U U**T | / ( n ulp )
!>
!> If ITYPE=2, then:
!>
!>         RESULT(1) = | A - V S V**T | / ( |A| n ulp )
!>
!> If ITYPE=3, then:
!>
!>         RESULT(1) = | I - V U**T | / ( n ulp )
!>
!> Packed storage means that, for example, if UPLO='U', then the columns
!> of the upper triangle of A are stored one after another, so that
!> A(1,j+1) immediately follows A(j,j) in the array AP.  Similarly, if
!> UPLO='L', then the columns of the lower triangle of A are stored one
!> after another in AP, so that A(j+1,j+1) immediately follows A(n,j)
!> in the array AP.  This means that A(i,j) is stored in:
!>
!>    AP( i + j*(j-1)/2 )                 if UPLO='U'
!>
!>    AP( i + (2*n-j)*(j-1)/2 )           if UPLO='L'
!>
!> The array VP bears the same relation to the matrix V that A does to
!> AP.
!>
!> For ITYPE > 1, the transformation U is expressed as a product
!> of Householder transformations:
!>
!>    If UPLO='U', then  V = H(n-1)...H(1),  where
!>
!>        H(j) = I  -  tau(j) v(j) v(j)**T
!>
!>    and the first j-1 elements of v(j) are stored in V(1:j-1,j+1),
!>    (i.e., VP( j*(j+1)/2 + 1 : j*(j+1)/2 + j-1 ) ),
!>    the j-th element is 1, and the last n-j elements are 0.
!>
!>    If UPLO='L', then  V = H(1)...H(n-1),  where
!>
!>        H(j) = I  -  tau(j) v(j) v(j)**T
!>
!>    and the first j elements of v(j) are 0, the (j+1)-st is 1, and the
!>    (j+2)-nd through n-th elements are stored in V(j+2:n,j) (i.e.,
!>    in VP( (2*n-j)*(j-1)/2 + j+2 : (2*n-j)*(j-1)/2 + n ) .)
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
!>             RESULT(1) = | A - U S U**T | / ( |A| n ulp ) and
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
!>          If UPLO='U', AP and VP are considered to contain the upper
!>          triangle of A and V.
!>          If UPLO='L', AP and VP are considered to contain the lower
!>          triangle of A and V.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The size of the matrix.  If it is zero, SSPT21 does nothing.
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
!> \param[in] AP
!> \verbatim
!>          AP is REAL array, dimension (N*(N+1)/2)
!>          The original (unfactored) matrix.  It is assumed to be
!>          symmetric, and contains the columns of just the upper
!>          triangle (UPLO='U') or only the lower triangle (UPLO='L'),
!>          packed one after another.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The diagonal of the (symmetric tri-) diagonal matrix.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          The off-diagonal of the (symmetric tri-) diagonal matrix.
!>          E(1) is the (1,2) and (2,1) element, E(2) is the (2,3) and
!>          (3,2) element, etc.
!>          Not referenced if KBAND=0.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is REAL array, dimension (LDU, N)
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
!> \param[in] VP
!> \verbatim
!>          VP is REAL array, dimension (N*(N+1)/2)
!>          If ITYPE=2 or 3, the columns of this array contain the
!>          Householder vectors used to describe the orthogonal matrix
!>          in the decomposition, as described in purpose.
!>          *NOTE* If ITYPE=2 or 3, V is modified and restored.  The
!>          subdiagonal (if UPLO='L') or the superdiagonal (if UPLO='U')
!>          is set to one, and later reset to its original value, during
!>          the course of the calculation.
!>          If ITYPE=1, then it is neither referenced nor modified.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is REAL array, dimension (N)
!>          If ITYPE >= 2, then TAU(j) is the scalar factor of
!>          v(j) v(j)**T in the Householder transformation H(j) of
!>          the product  U = H(1)...H(n-2)
!>          If ITYPE < 2, then TAU is not referenced.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (N**2+N)
!>          Workspace.
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (2)
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
!> \ingroup single_eig
!
!  =====================================================================
   SUBROUTINE SSPT21( ITYPE, UPLO, N, KBAND, AP, D, E, U, LDU, VP, &
                      TAU, WORK, RESULT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            ITYPE, KBAND, LDU, N
!     ..
!     .. Array Arguments ..
   REAL               AP( * ), D( * ), E( * ), RESULT( 2 ), TAU( * ), &
                      U( LDU, * ), VP( * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            LOWER
   CHARACTER          CUPLO
   INTEGER            IINFO, J, JP, JP1, JR, LAP
   REAL               ANORM, TEMP, ULP, UNFL, VSAVE, WNORM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               SDOT, SLAMCH, SLANGE, SLANSP
   EXTERNAL           LSAME, SDOT, SLAMCH, SLANGE, SLANSP
!     ..
!     .. External Subroutines ..
   EXTERNAL           SAXPY, SCOPY, SGEMM, SLACPY, SLASET, SOPMTR, &
                      SSPMV, SSPR, SSPR2
!     ..
!     .. Executable Statements ..
!
!     1)      Constants
!
   RESULT( 1 ) = 0.0E+0
   IF( ITYPE == 1 ) RESULT( 2 ) = 0.0E+0
   IF( N <= 0 ) RETURN
!
   LAP = ( N*( N+1 ) ) / 2
!
   IF( LSAME( UPLO, 'U' ) ) THEN
      LOWER = .FALSE.
      CUPLO = 'U'
   ELSE
      LOWER = .TRUE.
      CUPLO = 'L'
   END IF
!
   UNFL = SLAMCH( 'Safe minimum' )
   ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
!
!     Some Error Checks
!
   IF( ITYPE < 1 .OR. ITYPE > 3 ) THEN
      RESULT( 1 ) = 10.0E+0 / ULP
      RETURN
   END IF
!
!     Do Test 1
!
!     Norm of A:
!
   IF( ITYPE == 3 ) THEN
      ANORM = 1.0E+0
   ELSE
      ANORM = MAX( SLANSP( '1', CUPLO, N, AP, WORK ), UNFL )
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
      CALL SLASET( 'Full', N, N, 0.0E+0, 0.0E+0, WORK, N )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL SCOPY( LAP, AP, 1, WORK, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
      DO J = 1, N
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SSPR( CUPLO, N, -D( J ), U( 1, J ), 1, WORK )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SSPR : ',&
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
            CALL SSPR2( CUPLO, N, -E( J ), U( 1, J ), 1, U( 1, J+1 ), &
                        1, WORK )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SSPR2 : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ENDDO
      END IF
      WNORM = SLANSP( '1', CUPLO, N, WORK, WORK( N**2+1 ) )
!
   ELSE IF( ITYPE == 2 ) THEN
!
!        ITYPE=2: error = V S V**T - A
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL SLASET( 'Full', N, N, 0.0E+0, 0.0E+0, WORK, N )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
      IF( LOWER ) THEN
         WORK( LAP ) = D( N )
         DO J = N - 1, 1, -1
            JP = ( ( 2*N-J )*( J-1 ) ) / 2
            JP1 = JP + N - J
            IF( KBAND == 1 ) THEN
               WORK( JP+J+1 ) = ( 1.0E+0-TAU( J ) )*E( J )
               DO JR = J + 2, N
                  WORK( JP+JR ) = -TAU( J )*E( J )*VP( JP+JR )
               ENDDO
            END IF
!
            IF( TAU( J ) /= 0.0E+0 ) THEN
               VSAVE = VP( JP+J+1 )
               VP( JP+J+1 ) = 1.0E+0
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SSPMV( 'L', N-J, 1.0E+0, WORK( JP1+J+1 ), &
                           VP( JP+J+1 ), 1, 0.0E+0, WORK( LAP+1 ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SSPMV : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               TEMP = -0.5E+0*TAU( J )*SDOT( N-J, WORK( LAP+1 ), 1, &
                      VP( JP+J+1 ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SAXPY( N-J, TEMP, VP( JP+J+1 ), 1, WORK( LAP+1 ), &
                           1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SAXPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SSPR2( 'L', N-J, -TAU( J ), VP( JP+J+1 ), 1, &
                           WORK( LAP+1 ), 1, WORK( JP1+J+1 ) )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SSPR2 : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               VP( JP+J+1 ) = VSAVE
            END IF
            WORK( JP+J ) = D( J )
         ENDDO
      ELSE
         WORK( 1 ) = D( 1 )
         DO J = 1, N - 1
            JP = ( J*( J-1 ) ) / 2
            JP1 = JP + J
            IF( KBAND == 1 ) THEN
               WORK( JP1+J ) = ( 1.0E+0-TAU( J ) )*E( J )
               DO JR = 1, J - 1
                  WORK( JP1+JR ) = -TAU( J )*E( J )*VP( JP1+JR )
               ENDDO
            END IF
!
            IF( TAU( J ) /= 0.0E+0 ) THEN
               VSAVE = VP( JP1+J )
               VP( JP1+J ) = 1.0E+0
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SSPMV( 'U', J, 1.0E+0, WORK, VP( JP1+1 ), 1, 0.0E+0, &
                           WORK( LAP+1 ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SSPMV : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               TEMP = -0.5E+0*TAU( J )*SDOT( J, WORK( LAP+1 ), 1, &
                      VP( JP1+1 ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SAXPY( J, TEMP, VP( JP1+1 ), 1, WORK( LAP+1 ), &
                           1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SAXPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SSPR2( 'U', J, -TAU( J ), VP( JP1+1 ), 1, &
                           WORK( LAP+1 ), 1, WORK )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SSPR2 : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               VP( JP1+J ) = VSAVE
            END IF
            WORK( JP1+J+1 ) = D( J+1 )
         ENDDO
      END IF
!
      DO J = 1, LAP
         WORK( J ) = WORK( J ) - AP( J )
      ENDDO
      WNORM = SLANSP( '1', CUPLO, N, WORK, WORK( LAP+1 ) )
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
      CALL SLACPY( ' ', N, N, U, LDU, WORK, N )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL SOPMTR( 'R', CUPLO, 'T', N, N, VP, TAU, WORK, N, &
                   WORK( N**2+1 ), IINFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : SOPMTR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( IINFO /= 0 ) THEN
         RESULT( 1 ) = 10.0E+0 / ULP
         RETURN
      END IF
!
      DO J = 1, N
         WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - 1.0E+0
      ENDDO
!
      WNORM = SLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) )
   END IF
!
   IF( ANORM > WNORM ) THEN
      RESULT( 1 ) = ( WNORM / ANORM ) / ( N*ULP )
   ELSE
      IF( ANORM < 1.0E+0 ) THEN
         RESULT( 1 ) = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
      ELSE
         RESULT( 1 ) = MIN( WNORM / ANORM, REAL( N ) ) / ( N*ULP )
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
      CALL SGEMM( 'N', 'C', N, N, N, 1.0E+0, U, LDU, U, LDU, 0.0E+0, WORK, &
                  N )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
      DO J = 1, N
         WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - 1.0E+0
      ENDDO
!
      RESULT( 2 ) = MIN( SLANGE( '1', N, N, WORK, N, &
                    WORK( N**2+1 ) ), REAL( N ) ) / ( N*ULP )
   END IF
!
   RETURN
!
!     End of SSPT21
!
END




