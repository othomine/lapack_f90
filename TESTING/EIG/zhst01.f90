!> \brief \b ZHST01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHST01( N, ILO, IHI, A, LDA, H, LDH, Q, LDQ, WORK,
!                          LWORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, LDA, LDH, LDQ, LWORK, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RESULT( 2 ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), H( LDH, * ), Q( LDQ, * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHST01 tests the reduction of a general matrix A to upper Hessenberg
!> form:  A = Q*H*Q'.  Two test ratios are computed;
!>
!> RESULT(1) = norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )
!> RESULT(2) = norm( I - Q'*Q ) / ( N * EPS )
!>
!> The matrix Q is assumed to be given explicitly as it would be
!> following ZGEHRD + ZUNGHR.
!>
!> In this version, ILO and IHI are not used, but they could be used
!> to save some work if this is desired.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>
!>          A is assumed to be upper triangular in rows and columns
!>          1:ILO-1 and IHI+1:N, so Q differs from the identity only in
!>          rows and columns ILO+1:IHI.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The original n by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] H
!> \verbatim
!>          H is COMPLEX*16 array, dimension (LDH,N)
!>          The upper Hessenberg matrix H from the reduction A = Q*H*Q'
!>          as computed by ZGEHRD.  H is assumed to be zero below the
!>          first subdiagonal.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          The leading dimension of the array H.  LDH >= max(1,N).
!> \endverbatim
!>
!> \param[in] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDQ,N)
!>          The orthogonal matrix Q from the reduction A = Q*H*Q' as
!>          computed by ZGEHRD + ZUNGHR.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= max(1,N).
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
!>          The length of the array WORK.  LWORK >= 2*N*N.
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
!>          RESULT(1) = norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )
!>          RESULT(2) = norm( I - Q'*Q ) / ( N * EPS )
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
   SUBROUTINE ZHST01( N, ILO, IHI, A, LDA, H, LDH, Q, LDQ, WORK, &
                      LWORK, RWORK, RESULT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            IHI, ILO, LDA, LDH, LDQ, LWORK, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   RESULT( 2 ), RWORK( * )
   COMPLEX*16         A( LDA, * ), H( LDH, * ), Q( LDQ, * ), &
                      WORK( LWORK )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            LDWORK
   DOUBLE PRECISION   ANORM, EPS, OVFL, SMLNUM, UNFL, WNORM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH, ZLANGE
   EXTERNAL           DLAMCH, ZLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           ZGEMM, ZLACPY, ZUNT01
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
   IF( N <= 0 ) THEN
      RESULT( 1:2 ) = 0.0D0
      RETURN
   END IF
!
   UNFL = DLAMCH( 'Safe minimum' )
   EPS = DLAMCH( 'Precision' )
   OVFL = 1.0D0 / UNFL
   SMLNUM = UNFL*N / EPS
!
!     Test 1:  Compute norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )
!
!     Copy A to WORK
!
   LDWORK = MAX( 1, N )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( ' ', N, N, A, LDA, WORK, LDWORK )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute Q*H
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'No transpose', 'No transpose', N, N, N, &
               DCMPLX( 1.0D0 ), Q, LDQ, H, LDH, DCMPLX( 0.0D0 ), &
               WORK( LDWORK*N+1 ), LDWORK )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compute A - Q*H*Q'
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEMM( 'No transpose', 'Conjugate transpose', N, N, N, &
               DCMPLX( -1.0D0 ), WORK( LDWORK*N+1 ), LDWORK, Q, LDQ, &
               DCMPLX( 1.0D0 ), WORK, LDWORK )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   ANORM = MAX( ZLANGE( '1', N, N, A, LDA, RWORK ), UNFL )
   WNORM = ZLANGE( '1', N, N, WORK, LDWORK, RWORK )
!
!     Note that RESULT(1) cannot overflow and is bounded by 1/(N*EPS)
!
   RESULT( 1 ) = MIN( WNORM, ANORM ) / MAX( SMLNUM, ANORM*EPS ) / N
!
!     Test 2:  Compute norm( I - Q'*Q ) / ( N * EPS )
!
   CALL ZUNT01( 'Columns', N, N, Q, LDQ, WORK, LWORK, RWORK, &
                RESULT( 2 ) )
!
   RETURN
!
!     End of ZHST01
!
END




