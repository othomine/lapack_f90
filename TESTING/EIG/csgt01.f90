!> \brief \b CSGT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSGT01( ITYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D,
!                          WORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            ITYPE, LDA, LDB, LDZ, M, N
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), RESULT( * ), RWORK( * )
!       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSGT01 checks a decomposition of the form
!>
!>    A Z   =  B Z D or
!>    A B Z =  Z D or
!>    B A Z =  Z D
!>
!> where A is a Hermitian matrix, B is Hermitian positive definite,
!> Z is unitary, and D is diagonal.
!>
!> One of the following test ratios is computed:
!>
!> ITYPE = 1:  RESULT(1) = | A Z - B Z D | / ( |A| |Z| n ulp )
!>
!> ITYPE = 2:  RESULT(1) = | A B Z - Z D | / ( |A| |Z| n ulp )
!>
!> ITYPE = 3:  RESULT(1) = | B A Z - Z D | / ( |A| |Z| n ulp )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          The form of the Hermitian generalized eigenproblem.
!>          = 1:  A*z = (lambda)*B*z
!>          = 2:  A*B*z = (lambda)*z
!>          = 3:  B*A*z = (lambda)*z
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          Hermitian matrices A and B is stored.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of eigenvalues found.  M >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA, N)
!>          The original Hermitian matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB, N)
!>          The original Hermitian positive definite matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ, M)
!>          The computed eigenvectors of the generalized eigenproblem.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (M)
!>          The computed eigenvalues of the generalized eigenproblem.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (N*N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (1)
!>          The test ratio as described above.
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
   SUBROUTINE CSGT01( ITYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, &
                      WORK, RWORK, RESULT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            ITYPE, LDA, LDB, LDZ, M, N
!     ..
!     .. Array Arguments ..
   REAL               D( * ), RESULT( * ), RWORK( * )
   COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ), &
                      Z( LDZ, * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I
   REAL               ANORM, ULP
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   REAL               CLANGE, CLANHE, SLAMCH
   EXTERNAL           CLANGE, CLANHE, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           CHEMM, CSSCAL
!     ..
!     .. Executable Statements ..
!
   RESULT( 1 ) = 0.0E+0
   IF( N <= 0 ) RETURN
!
   ULP = SLAMCH( 'Epsilon' )
!
!     Compute product of 1-norms of A and Z.
!
   ANORM = CLANHE( '1', UPLO, N, A, LDA, RWORK )* &
           CLANGE( '1', N, M, Z, LDZ, RWORK )
   IF( ANORM == 0.0E+0 ) &
      ANORM = 1.0E+0
!
   IF( ITYPE == 1 ) THEN
!
!        Norm of AZ - BZD
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEMM( 'Left', UPLO, N, M, (1.0E+0,0.0E+0), A, LDA, Z, LDZ, (0.0E+0,0.0E+0), WORK, N )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEMM : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      DO I = 1, M
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL CSSCAL( N, D( I ), Z( 1, I ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CSSCAL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEMM( 'Left', UPLO, N, M, (1.0E+0,0.0E+0), B, LDB, Z, LDZ, -(1.0E+0,0.0E+0), WORK, N )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEMM : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
      RESULT( 1 ) = ( CLANGE( '1', N, M, WORK, N, RWORK ) / ANORM ) / ( N*ULP )
!
   ELSE IF( ITYPE == 2 ) THEN
!
!        Norm of ABZ - ZD
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEMM( 'Left', UPLO, N, M, (1.0E+0,0.0E+0), B, LDB, Z, LDZ, (0.0E+0,0.0E+0), WORK, N )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEMM : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      DO I = 1, M
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL CSSCAL( N, D( I ), Z( 1, I ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CSSCAL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEMM( 'Left', UPLO, N, M, (1.0E+0,0.0E+0), A, LDA, WORK, N, -(1.0E+0,0.0E+0), Z, LDZ )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEMM : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
      RESULT( 1 ) = ( CLANGE( '1', N, M, Z, LDZ, RWORK ) / ANORM ) / ( N*ULP )
!
   ELSE IF( ITYPE == 3 ) THEN
!
!        Norm of BAZ - ZD
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEMM( 'Left', UPLO, N, M, (1.0E+0,0.0E+0), A, LDA, Z, LDZ, (0.0E+0,0.0E+0), WORK, N )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEMM : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      DO I = 1, M
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL CSSCAL( N, D( I ), Z( 1, I ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : CSSCAL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CHEMM( 'Left', UPLO, N, M, (1.0E+0,0.0E+0), B, LDB, WORK, N, -(1.0E+0,0.0E+0), Z, LDZ )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CHEMM : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
      RESULT( 1 ) = ( CLANGE( '1', N, M, Z, LDZ, RWORK ) / ANORM ) / ( N*ULP )
   END IF
!
   RETURN
!
!     End of CSGT01
!
END




