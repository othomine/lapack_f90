!> \brief \b SBDT04
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SBDT04( UPLO, N, D, E, S, NS, U, LDU, VT, LDVT,
!                          WORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDU, LDVT, N, NS
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * ), S( * ), U( LDU, * ),
!      $                   VT( LDVT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SBDT04 reconstructs a bidiagonal matrix B from its (partial) SVD:
!>    S = U' * B * V
!> where U and V are orthogonal matrices and S is diagonal.
!>
!> The test ratio to test the singular value decomposition is
!>    RESID = norm( S - U' * B * V ) / ( n * norm(B) * EPS )
!> where VT = V' and EPS is the machine precision.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix B is upper or lower bidiagonal.
!>          = 'U':  Upper bidiagonal
!>          = 'L':  Lower bidiagonal
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix B.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The n diagonal elements of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          The (n-1) superdiagonal elements of the bidiagonal matrix B
!>          if UPLO = 'U', or the (n-1) subdiagonal elements of B if
!>          UPLO = 'L'.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is REAL array, dimension (NS)
!>          The singular values from the (partial) SVD of B, sorted in
!>          decreasing order.
!> \endverbatim
!>
!> \param[in] NS
!> \verbatim
!>          NS is INTEGER
!>          The number of singular values/vectors from the (partial)
!>          SVD of B.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is REAL array, dimension (LDU,NS)
!>          The n by ns orthogonal matrix U in S = U' * B * V.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= max(1,N)
!> \endverbatim
!>
!> \param[in] VT
!> \verbatim
!>          VT is REAL array, dimension (LDVT,N)
!>          The n by ns orthogonal matrix V in S = U' * B * V.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>          The leading dimension of the array VT.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
!>          The test ratio:  norm(S - U' * B * V) / ( n * norm(B) * EPS )
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
   SUBROUTINE SBDT04( UPLO, N, D, E, S, NS, U, LDU, VT, LDVT, WORK, &
                      RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            LDU, LDVT, N, NS
   REAL               RESID
!     ..
!     .. Array Arguments ..
   REAL               D( * ), E( * ), S( * ), U( LDU, * ), &
                      VT( LDVT, * ), WORK( * )
!     ..
!
! ======================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, J, K
   REAL               BNORM, EPS
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ISAMAX
   REAL               SASUM, SLAMCH
   EXTERNAL           LSAME, ISAMAX, SASUM, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           SGEMM
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible.
!
   RESID = 0.0E+0
   IF( N <= 0 .OR. NS <= 0 ) RETURN
!
   EPS = SLAMCH( 'Precision' )
!
!     Compute S - U' * B * V.
!
   BNORM = 0.0E+0
!
   IF( LSAME( UPLO, 'U' ) ) THEN
!
!        B is upper bidiagonal.
!
      K = 0
      DO I = 1, NS
         DO J = 1, N-1
            K = K + 1
            WORK( K ) = D( J )*VT( I, J ) + E( J )*VT( I, J+1 )
         ENDDO
         K = K + 1
         WORK( K ) = D( N )*VT( I, N )
      ENDDO
      BNORM = ABS( D( 1 ) )
      DO I = 2, N
         BNORM = MAX( BNORM, ABS( D( I ) )+ABS( E( I-1 ) ) )
      ENDDO
   ELSE
!
!        B is lower bidiagonal.
!
      K = 0
      DO I = 1, NS
         K = K + 1
         WORK( K ) = D( 1 )*VT( I, 1 )
         DO J = 1, N-1
            K = K + 1
            WORK( K ) = E( J )*VT( I, J ) + D( J+1 )*VT( I, J+1 )
         ENDDO
      ENDDO
      BNORM = ABS( D( N ) )
      DO I = 1, N-1
         BNORM = MAX( BNORM, ABS( D( I ) )+ABS( E( I ) ) )
      ENDDO
   END IF
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SGEMM( 'T', 'N', NS, NS, N, -1.0E+0, U, LDU, WORK( 1 ), &
               N, 0.0E+0, WORK( 1+N*NS ), NS )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     norm(S - U' * B * V)
!
   K = N*NS
   DO I = 1, NS
      WORK( K+I ) =  WORK( K+I ) + S( I )
      RESID = MAX( RESID, SASUM( NS, WORK( K+1 ), 1 ) )
      K = K + NS
   ENDDO
!
   IF( BNORM <= 0.0E+0 ) THEN
      IF( RESID /= 0.0E+0 ) &
         RESID = 1.0E+0 / EPS
   ELSE
      IF( BNORM >= RESID ) THEN
         RESID = ( RESID / BNORM ) / ( REAL( N )*EPS )
      ELSE
         IF( BNORM < 1.0E+0 ) THEN
            RESID = ( MIN( RESID, REAL( N )*BNORM ) / BNORM ) / &
                    ( REAL( N )*EPS )
         ELSE
            RESID = MIN( RESID / BNORM, REAL( N ) ) / &
                    ( REAL( N )*EPS )
         END IF
      END IF
   END IF
!
   RETURN
!
!     End of SBDT04
!
END




