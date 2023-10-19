!> \brief \b SGET51
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGET51( ITYPE, N, A, LDA, B, LDB, U, LDU, V, LDV, WORK,
!                          RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            ITYPE, LDA, LDB, LDU, LDV, N
!       REAL               RESULT
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), B( LDB, * ), U( LDU, * ),
!      $                   V( LDV, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      SGET51  generally checks a decomposition of the form
!>
!>              A = U B V'
!>
!>      where ' means transpose and U and V are orthogonal.
!>
!>      Specifically, if ITYPE=1
!>
!>              RESULT = | A - U B V' | / ( |A| n ulp )
!>
!>      If ITYPE=2, then:
!>
!>              RESULT = | A - B | / ( |A| n ulp )
!>
!>      If ITYPE=3, then:
!>
!>              RESULT = | I - UU' | / ( n ulp )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          Specifies the type of tests to be performed.
!>          =1: RESULT = | A - U B V' | / ( |A| n ulp )
!>          =2: RESULT = | A - B | / ( |A| n ulp )
!>          =3: RESULT = | I - UU' | / ( n ulp )
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The size of the matrix.  If it is zero, SGET51 does nothing.
!>          It must be at least zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA, N)
!>          The original (unfactored) matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  It must be at least 1
!>          and at least N.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (LDB, N)
!>          The factored matrix.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of B.  It must be at least 1
!>          and at least N.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is REAL array, dimension (LDU, N)
!>          The orthogonal matrix on the left-hand side in the
!>          decomposition.
!>          Not referenced if ITYPE=2
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
!>          V is REAL array, dimension (LDV, N)
!>          The orthogonal matrix on the left-hand side in the
!>          decomposition.
!>          Not referenced if ITYPE=2
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of V.  LDV must be at least N and
!>          at least 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (2*N**2)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL
!>          The values computed by the test specified by ITYPE.  The
!>          value is currently limited to 1/ulp, to avoid overflow.
!>          Errors are flagged by RESULT=10/ulp.
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
   SUBROUTINE SGET51( ITYPE, N, A, LDA, B, LDB, U, LDU, V, LDV, WORK, &
                      RESULT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            ITYPE, LDA, LDB, LDU, LDV, N
   REAL               RESULT
!     ..
!     .. Array Arguments ..
   REAL               A( LDA, * ), B( LDB, * ), U( LDU, * ), &
                      V( LDV, * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            JCOL, JDIAG, JROW
   REAL               ANORM, ULP, UNFL, WNORM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   REAL               SLAMCH, SLANGE
   EXTERNAL           SLAMCH, SLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           SGEMM, SLACPY
!     ..
!     .. Executable Statements ..
!
   RESULT = 0.0E+0
   IF( N <= 0 ) RETURN
!
!     Constants
!
   UNFL = SLAMCH( 'Safe minimum' )
   ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
!
!     Some Error Checks
!
   IF( ITYPE < 1 .OR. ITYPE > 3 ) THEN
      RESULT = 10.0E+0 / ULP
      RETURN
   END IF
!
   IF( ITYPE <= 2 ) THEN
!
!        Tests scaled by the norm(A)
!
      ANORM = MAX( SLANGE( '1', N, N, A, LDA, WORK ), UNFL )
!
      IF( ITYPE == 1 ) THEN
!
!           ITYPE=1: Compute W = A - UBV'
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SLACPY( ' ', N, N, A, LDA, WORK, N )
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
         CALL SGEMM( 'N', 'N', N, N, N, 1.0E+0, U, LDU, B, LDB, 0.0E+0, &
                     WORK( N**2+1 ), N )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SGEMM( 'N', 'C', N, N, N, -1.0E+0, WORK( N**2+1 ), N, V, &
                     LDV, 1.0E+0, WORK, N )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
      ELSE
!
!           ITYPE=2: Compute W = A - B
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SLACPY( ' ', N, N, B, LDB, WORK, N )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
         DO JCOL = 1, N
            DO JROW = 1, N
               WORK( JROW+N*( JCOL-1 ) ) = WORK( JROW+N*( JCOL-1 ) ) &
                   - A( JROW, JCOL )
            ENDDO
         ENDDO
      END IF
!
!        Compute norm(W)/ ( ulp*norm(A) )
!
      WNORM = SLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) )
!
      IF( ANORM > WNORM ) THEN
         RESULT = ( WNORM / ANORM ) / ( N*ULP )
      ELSE
         IF( ANORM < 1.0E+0 ) THEN
            RESULT = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
         ELSE
            RESULT = MIN( WNORM / ANORM, REAL( N ) ) / ( N*ULP )
         END IF
      END IF
!
   ELSE
!
!        Tests not scaled by norm(A)
!
!        ITYPE=3: Compute  UU' - I
!
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
      DO JDIAG = 1, N
         WORK( ( N+1 )*( JDIAG-1 )+1 ) = WORK( ( N+1 )*( JDIAG-1 )+ &
            1 ) - 1.0E+0
      ENDDO
!
      RESULT = MIN( SLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) ), &
               REAL( N ) ) / ( N*ULP )
   END IF
!
   RETURN
!
!     End of SGET51
!
END




