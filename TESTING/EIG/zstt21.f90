!> \brief \b ZSTT21
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZSTT21( N, KBAND, AD, AE, SD, SE, U, LDU, WORK, RWORK,
!                          RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            KBAND, LDU, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AD( * ), AE( * ), RESULT( 2 ), RWORK( * ),
!      $                   SD( * ), SE( * )
!       COMPLEX*16         U( LDU, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZSTT21  checks a decomposition of the form
!>
!>    A = U S U**H
!>
!> where **H means conjugate transpose, A is real symmetric tridiagonal,
!> U is unitary, and S is real and diagonal (if KBAND=0) or symmetric
!> tridiagonal (if KBAND=1).  Two tests are performed:
!>
!>    RESULT(1) = | A - U S U**H | / ( |A| n ulp )
!>
!>    RESULT(2) = | I - U U**H | / ( n ulp )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The size of the matrix.  If it is zero, ZSTT21 does nothing.
!>          It must be at least zero.
!> \endverbatim
!>
!> \param[in] KBAND
!> \verbatim
!>          KBAND is INTEGER
!>          The bandwidth of the matrix S.  It may only be zero or one.
!>          If zero, then S is diagonal, and SE is not referenced.  If
!>          one, then S is symmetric tri-diagonal.
!> \endverbatim
!>
!> \param[in] AD
!> \verbatim
!>          AD is DOUBLE PRECISION array, dimension (N)
!>          The diagonal of the original (unfactored) matrix A.  A is
!>          assumed to be real symmetric tridiagonal.
!> \endverbatim
!>
!> \param[in] AE
!> \verbatim
!>          AE is DOUBLE PRECISION array, dimension (N-1)
!>          The off-diagonal of the original (unfactored) matrix A.  A
!>          is assumed to be symmetric tridiagonal.  AE(1) is the (1,2)
!>          and (2,1) element, AE(2) is the (2,3) and (3,2) element, etc.
!> \endverbatim
!>
!> \param[in] SD
!> \verbatim
!>          SD is DOUBLE PRECISION array, dimension (N)
!>          The diagonal of the real (symmetric tri-) diagonal matrix S.
!> \endverbatim
!>
!> \param[in] SE
!> \verbatim
!>          SE is DOUBLE PRECISION array, dimension (N-1)
!>          The off-diagonal of the (symmetric tri-) diagonal matrix S.
!>          Not referenced if KBSND=0.  If KBAND=1, then AE(1) is the
!>          (1,2) and (2,1) element, SE(2) is the (2,3) and (3,2)
!>          element, etc.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is COMPLEX*16 array, dimension (LDU, N)
!>          The unitary matrix in the decomposition.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U.  LDU must be at least N.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (N**2)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (2)
!>          The values computed by the two tests described above.  The
!>          values are currently limited to 1/ulp, to avoid overflow.
!>          RESULT(1) is always modified.
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
   SUBROUTINE ZSTT21( N, KBAND, AD, AE, SD, SE, U, LDU, WORK, RWORK, &
                      RESULT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            KBAND, LDU, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   AD( * ), AE( * ), RESULT( 2 ), RWORK( * ), &
                      SD( * ), SE( * )
   COMPLEX*16         U( LDU, * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            J
   DOUBLE PRECISION   ANORM, TEMP1, TEMP2, ULP, UNFL, WNORM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH, ZLANGE, ZLANHE
   EXTERNAL           DLAMCH, ZLANGE, ZLANHE
!     ..
!     .. External Subroutines ..
   EXTERNAL           ZGEMM, ZHER, ZHER2, ZLASET
!     ..
!     .. Executable Statements ..
!
!     1)      Constants
!
   RESULT( 1:2 ) = 0.0D0
   IF( N <= 0 ) RETURN
!
   UNFL = DLAMCH( 'Safe minimum' )
   ULP = DLAMCH( 'Precision' )
!
!     Do Test 1
!
!     Copy A & Compute its 1-Norm:
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', N, N, (0.0D+0,0.0D+0), (0.0D+0,0.0D+0), WORK, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   ANORM = 0.0D0
   TEMP1 = 0.0D0
!
   DO J = 1, N - 1
      WORK( ( N+1 )*( J-1 )+1 ) = AD( J )
      WORK( ( N+1 )*( J-1 )+2 ) = AE( J )
      TEMP2 = ABS( AE( J ) )
      ANORM = MAX( ANORM, ABS( AD( J ) )+TEMP1+TEMP2 )
      TEMP1 = TEMP2
   ENDDO
!
   WORK( N**2 ) = AD( N )
   ANORM = MAX( ANORM, ABS( AD( N ) )+TEMP1, UNFL )
!
!     Norm of A - USU*
!
   DO J = 1, N
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZHER( 'L', N, -SD( J ), U( 1, J ), 1, WORK, N )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZHER : ',&
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
         CALL ZHER2( 'L', N, -DCMPLX( SE( J ) ), U( 1, J ), 1, &
                     U( 1, J+1 ), 1, WORK, N )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZHER2 : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
   END IF
!
   WNORM = ZLANHE( '1', 'L', N, WORK, N, RWORK )
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
!     Compute  U U**H - I
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'N', 'C', N, N, N, (1.0D0,0.0D0), U, LDU, U, LDU, (0.0D+0,0.0D+0), WORK, &
               N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   DO J = 1, N
      WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - (1.0D0,0.0D0)
   ENDDO
!
   RESULT( 2 ) = MIN( DBLE( N ), ZLANGE( '1', N, N, WORK, N, &
                 RWORK ) ) / ( N*ULP )
!
   RETURN
!
!     End of ZSTT21
!
END




