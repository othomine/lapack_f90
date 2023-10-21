!> \brief \b ZLATM6
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLATM6( TYPE, N, A, LDA, B, X, LDX, Y, LDY, ALPHA,
!                          BETA, WX, WY, S, DIF )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDX, LDY, N, TYPE
!       COMPLEX*16         ALPHA, BETA, WX, WY
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   DIF( * ), S( * )
!       COMPLEX*16         A( LDA, * ), B( LDA, * ), X( LDX, * ),
!      $                   Y( LDY, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLATM6 generates test matrices for the generalized eigenvalue
!> problem, their corresponding right and left eigenvector matrices,
!> and also reciprocal condition numbers for all eigenvalues and
!> the reciprocal condition numbers of eigenvectors corresponding to
!> the 1th and 5th eigenvalues.
!>
!> Test Matrices
!> =============
!>
!> Two kinds of test matrix pairs
!>          (A, B) = inverse(YH) * (Da, Db) * inverse(X)
!> are used in the tests:
!>
!> Type 1:
!>    Da = 1+a   0    0    0    0    Db = 1   0   0   0   0
!>          0   2+a   0    0    0         0   1   0   0   0
!>          0    0   3+a   0    0         0   0   1   0   0
!>          0    0    0   4+a   0         0   0   0   1   0
!>          0    0    0    0   5+a ,      0   0   0   0   1
!> and Type 2:
!>    Da = 1+i   0    0       0       0    Db = 1   0   0   0   0
!>          0   1-i   0       0       0         0   1   0   0   0
!>          0    0    1       0       0         0   0   1   0   0
!>          0    0    0 (1+a)+(1+b)i  0         0   0   0   1   0
!>          0    0    0       0 (1+a)-(1+b)i,   0   0   0   0   1 .
!>
!> In both cases the same inverse(YH) and inverse(X) are used to compute
!> (A, B), giving the exact eigenvectors to (A,B) as (YH, X):
!>
!> YH:  =  1    0   -y    y   -y    X =  1   0  -x  -x   x
!>         0    1   -y    y   -y         0   1   x  -x  -x
!>         0    0    1    0    0         0   0   1   0   0
!>         0    0    0    1    0         0   0   0   1   0
!>         0    0    0    0    1,        0   0   0   0   1 , where
!>
!> a, b, x and y will have all values independently of each other.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TYPE
!> \verbatim
!>          TYPE is INTEGER
!>          Specifies the problem type (see further details).
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          Size of the matrices A and B.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA, N).
!>          On exit A N-by-N is initialized according to TYPE.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A and of B.
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDA, N).
!>          On exit B N-by-N is initialized according to TYPE.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (LDX, N).
!>          On exit X is the N-by-N matrix of right eigenvectors.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of X.
!> \endverbatim
!>
!> \param[out] Y
!> \verbatim
!>          Y is COMPLEX*16 array, dimension (LDY, N).
!>          On exit Y is the N-by-N matrix of left eigenvectors.
!> \endverbatim
!>
!> \param[in] LDY
!> \verbatim
!>          LDY is INTEGER
!>          The leading dimension of Y.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX*16
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is COMPLEX*16
!> \verbatim
!>          Weighting constants for matrix A.
!> \endverbatim
!>
!> \param[in] WX
!> \verbatim
!>          WX is COMPLEX*16
!>          Constant for right eigenvector matrix.
!> \endverbatim
!>
!> \param[in] WY
!> \verbatim
!>          WY is COMPLEX*16
!>          Constant for left eigenvector matrix.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (N)
!>          S(i) is the reciprocal condition number for eigenvalue i.
!> \endverbatim
!>
!> \param[out] DIF
!> \verbatim
!>          DIF is DOUBLE PRECISION array, dimension (N)
!>          DIF(i) is the reciprocal condition number for eigenvector i.
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
!> \ingroup complex16_matgen
!
!  =====================================================================
   SUBROUTINE ZLATM6( TYPE, N, A, LDA, B, X, LDX, Y, LDY, ALPHA, &
                      BETA, WX, WY, S, DIF )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, LDX, LDY, N, TYPE
   COMPLEX*16         ALPHA, BETA, WX, WY
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   DIF( * ), S( * )
   COMPLEX*16         A( LDA, * ), B( LDA, * ), X( LDX, * ), &
                      Y( LDY, * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, INFO, J
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   RWORK( 50 )
   COMPLEX*16         WORK( 26 ), Z( 8, 8 )
!     ..
!     .. External Subroutines ..
   EXTERNAL           ZGESVD, ZLACPY, ZLAKF2
!     ..
!     .. Executable Statements ..
!
!     Generate test problem ...
!     (Da, Db) ...
!
   A(1:N,1:N) = (0.0D+0,0.0D+0)
   B(1:N,1:N) = (0.0D+0,0.0D+0)
   FORALL (I=1:N) A( I, I ) = DCMPLX( I ) + ALPHA
   FORALL (I=1:N) B( I, I ) = (1.0E+0,0.0E+0)

   IF( TYPE == 2 ) THEN
      A( 1, 1 ) = DCMPLX( 1.0D+0, 1.0D+0 )
      A( 2, 2 ) = DCONJG( A( 1, 1 ) )
      A( 3, 3 ) = (1.0D+0,0.0D+0)
      A( 4, 4 ) = DCMPLX( DBLE( (1.0D+0,0.0D+0)+ALPHA ), DBLE( (1.0D+0,0.0D+0)+BETA ) )
      A( 5, 5 ) = DCONJG( A( 4, 4 ) )
   END IF
!
!     Form X and Y
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'F', N, N, B, LDA, Y, LDY )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   Y( 3, 1 ) = -DCONJG( WY )
   Y( 4, 1 ) = DCONJG( WY )
   Y( 5, 1 ) = -DCONJG( WY )
   Y( 3, 2 ) = -DCONJG( WY )
   Y( 4, 2 ) = DCONJG( WY )
   Y( 5, 2 ) = -DCONJG( WY )
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'F', N, N, B, LDA, X, LDX )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   X( 1, 3 ) = -WX
   X( 1, 4 ) = -WX
   X( 1, 5 ) = WX
   X( 2, 3 ) = WX
   X( 2, 4 ) = -WX
   X( 2, 5 ) = -WX
!
!     Form (A, B)
!
   B( 1, 3 ) = WX + WY
   B( 2, 3 ) = -WX + WY
   B( 1, 4 ) = WX - WY
   B( 2, 4 ) = WX - WY
   B( 1, 5 ) = -WX + WY
   B( 2, 5 ) = WX + WY
   A( 1, 3 ) = WX*A( 1, 1 ) + WY*A( 3, 3 )
   A( 2, 3 ) = -WX*A( 2, 2 ) + WY*A( 3, 3 )
   A( 1, 4 ) = WX*A( 1, 1 ) - WY*A( 4, 4 )
   A( 2, 4 ) = WX*A( 2, 2 ) - WY*A( 4, 4 )
   A( 1, 5 ) = -WX*A( 1, 1 ) + WY*A( 5, 5 )
   A( 2, 5 ) = WX*A( 2, 2 ) + WY*A( 5, 5 )
!
!     Compute condition numbers
!
   S( 1 ) = 1.0D+0 / SQRT( ( 1.0D+0+3.0D+0*CDABS( WY )*CDABS( WY ) ) / &
            ( 1.0D+0+CDABS( A( 1, 1 ) )*CDABS( A( 1, 1 ) ) ) )
   S( 2 ) = 1.0D+0 / SQRT( ( 1.0D+0+3.0D+0*CDABS( WY )*CDABS( WY ) ) / &
            ( 1.0D+0+CDABS( A( 2, 2 ) )*CDABS( A( 2, 2 ) ) ) )
   S( 3 ) = 1.0D+0 / SQRT( ( 1.0D+0+2.0D+0*CDABS( WX )*CDABS( WX ) ) / &
            ( 1.0D+0+CDABS( A( 3, 3 ) )*CDABS( A( 3, 3 ) ) ) )
   S( 4 ) = 1.0D+0 / SQRT( ( 1.0D+0+2.0D+0*CDABS( WX )*CDABS( WX ) ) / &
            ( 1.0D+0+CDABS( A( 4, 4 ) )*CDABS( A( 4, 4 ) ) ) )
   S( 5 ) = 1.0D+0 / SQRT( ( 1.0D+0+2.0D+0*CDABS( WX )*CDABS( WX ) ) / &
            ( 1.0D+0+CDABS( A( 5, 5 ) )*CDABS( A( 5, 5 ) ) ) )
!
   CALL ZLAKF2( 1, 4, A, LDA, A( 2, 2 ), B, B( 2, 2 ), Z, 8 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGESVD( 'N', 'N', 8, 8, Z, 8, RWORK, WORK, 1, WORK( 2 ), 1, &
                WORK( 3 ), 24, RWORK( 9 ), INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGESVD : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   DIF( 1 ) = RWORK( 8 )
!
   CALL ZLAKF2( 4, 1, A, LDA, A( 5, 5 ), B, B( 5, 5 ), Z, 8 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGESVD( 'N', 'N', 8, 8, Z, 8, RWORK, WORK, 1, WORK( 2 ), 1, &
                WORK( 3 ), 24, RWORK( 9 ), INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGESVD : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   DIF( 5 ) = RWORK( 8 )
!
   RETURN
!
!     End of ZLATM6
!
END

