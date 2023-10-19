!> \brief \b DGET10
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGET10( M, N, A, LDA, B, LDB, WORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, M, N
!       DOUBLE PRECISION   RESULT
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGET10 compares two matrices A and B and computes the ratio
!> RESULT = norm( A - B ) / ( norm(A) * M * EPS )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrices A and B.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrices A and B.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>          The m by n matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION
!>          RESULT = norm( A - B ) / ( norm(A) * M * EPS )
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
   SUBROUTINE DGET10( M, N, A, LDA, B, LDB, WORK, RESULT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, LDB, M, N
   DOUBLE PRECISION   RESULT
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            J
   DOUBLE PRECISION   ANORM, EPS, UNFL, WNORM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DASUM, DLAMCH, DLANGE
   EXTERNAL           DASUM, DLAMCH, DLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           DAXPY, DCOPY
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
   IF( M <= 0 .OR. N <= 0 ) THEN
      RESULT = 0.0D+0
      RETURN
   END IF
!
   UNFL = DLAMCH( 'Safe minimum' )
   EPS = DLAMCH( 'Precision' )
!
   WNORM = 0.0D+0
   DO J = 1, N
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( M, A( 1, J ), 1, WORK, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DAXPY( M, -1.0D+0, B( 1, J ), 1, WORK, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DAXPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      WNORM = MAX( WNORM, DASUM( N, WORK, 1 ) )
   ENDDO
!
   ANORM = MAX( DLANGE( '1', M, N, A, LDA, WORK ), UNFL )
!
   IF( ANORM > WNORM ) THEN
      RESULT = ( WNORM / ANORM ) / ( M*EPS )
   ELSE
      IF( ANORM < 1.0D+0 ) THEN
         RESULT = ( MIN( WNORM, M*ANORM ) / ANORM ) / ( M*EPS )
      ELSE
         RESULT = MIN( WNORM / ANORM, DBLE( M ) ) / ( M*EPS )
      END IF
   END IF
!
   RETURN
!
!     End of DGET10
!
END




