!> \brief \b SLAGTS solves the system of equations (T-λI)x = y
!> or (T-λI)^Tx = y, where T is a general tridiagonal matrix
!> and λ a scalar, using the LU factorization computed by slagtf.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAGTS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slagts.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slagts.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slagts.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAGTS( JOB, N, A, B, C, D, IN, Y, TOL, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, JOB, N
!       REAL               TOL
!       ..
!       .. Array Arguments ..
!       INTEGER            IN( * )
!       REAL               A( * ), B( * ), C( * ), D( * ), Y( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAGTS may be used to solve one of the systems of equations
!>
!>    (T - lambda*I)*x = y   or   (T - lambda*I)**T*x = y,
!>
!> where T is an n by n tridiagonal matrix, for x, following the
!> factorization of (T - lambda*I) as
!>
!>    (T - lambda*I) = P*L*U ,
!>
!> by routine SLAGTF. The choice of equation to be solved is
!> controlled by the argument JOB, and in each case there is an option
!> to perturb zero or very small diagonal elements of U, this option
!> being intended for use in applications such as inverse iteration.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is INTEGER
!>          Specifies the job to be performed by SLAGTS as follows:
!>          =  1: The equations  (T - lambda*I)x = y  are to be solved,
!>                but diagonal elements of U are not to be perturbed.
!>          = -1: The equations  (T - lambda*I)x = y  are to be solved
!>                and, if overflow would otherwise occur, the diagonal
!>                elements of U are to be perturbed. See argument TOL
!>                below.
!>          =  2: The equations  (T - lambda*I)**Tx = y  are to be solved,
!>                but diagonal elements of U are not to be perturbed.
!>          = -2: The equations  (T - lambda*I)**Tx = y  are to be solved
!>                and, if overflow would otherwise occur, the diagonal
!>                elements of U are to be perturbed. See argument TOL
!>                below.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix T.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (N)
!>          On entry, A must contain the diagonal elements of U as
!>          returned from SLAGTF.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (N-1)
!>          On entry, B must contain the first super-diagonal elements of
!>          U as returned from SLAGTF.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is REAL array, dimension (N-1)
!>          On entry, C must contain the sub-diagonal elements of L as
!>          returned from SLAGTF.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N-2)
!>          On entry, D must contain the second super-diagonal elements
!>          of U as returned from SLAGTF.
!> \endverbatim
!>
!> \param[in] IN
!> \verbatim
!>          IN is INTEGER array, dimension (N)
!>          On entry, IN must contain details of the matrix P as returned
!>          from SLAGTF.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is REAL array, dimension (N)
!>          On entry, the right hand side vector y.
!>          On exit, Y is overwritten by the solution vector x.
!> \endverbatim
!>
!> \param[in,out] TOL
!> \verbatim
!>          TOL is REAL
!>          On entry, with  JOB < 0, TOL should be the minimum
!>          perturbation to be made to very small diagonal elements of U.
!>          TOL should normally be chosen as about eps*norm(U), where eps
!>          is the relative machine precision, but if TOL is supplied as
!>          non-positive, then it is reset to eps*max( abs( u(i,j) ) ).
!>          If  JOB > 0  then TOL is not referenced.
!>
!>          On exit, TOL is changed as described above, only if TOL is
!>          non-positive on entry. Otherwise TOL is unchanged.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: overflow would occur when computing the INFO(th)
!>               element of the solution vector x. This can only occur
!>               when JOB is supplied as positive and either means
!>               that a diagonal element of U is very small, or that
!>               the elements of the right-hand side vector y are very
!>               large.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup lagts
!
!  =====================================================================
   SUBROUTINE SLAGTS( JOB, N, A, B, C, D, IN, Y, TOL, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, JOB, N
   REAL               TOL
!     ..
!     .. Array Arguments ..
   INTEGER            IN( * )
   REAL               A( * ), B( * ), C( * ), D( * ), Y( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ONE, ZERO
   PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            K
   REAL               ABSAK, AK, BIGNUM, EPS, PERT, SFMIN, TEMP
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, MAX, SIGN
!     ..
!     .. External Functions ..
   REAL               SLAMCH
   EXTERNAL           SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
   INFO = 0
   IF( ( ABS( JOB ) > 2 ) .OR. ( JOB == 0 ) ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'SLAGTS', -INFO )
      RETURN
   END IF
!
   IF( N == 0 ) &
      RETURN
!
   EPS = SLAMCH( 'Epsilon' )
   SFMIN = SLAMCH( 'Safe minimum' )
   BIGNUM = ONE / SFMIN
!
   IF( JOB < 0 ) THEN
      IF( TOL <= ZERO ) THEN
         TOL = ABS( A( 1 ) )
         IF( N > 1 ) &
            TOL = MAX( TOL, ABS( A( 2 ) ), ABS( B( 1 ) ) )
         DO K = 3, N
            TOL = MAX( TOL, ABS( A( K ) ), ABS( B( K-1 ) ), &
                  ABS( D( K-2 ) ) )
         ENDDO
         TOL = TOL*EPS
         IF( TOL == ZERO ) &
            TOL = EPS
      END IF
   END IF
!
   IF( ABS( JOB ) == 1 ) THEN
      DO K = 2, N
         IF( IN( K-1 ) == 0 ) THEN
            Y( K ) = Y( K ) - C( K-1 )*Y( K-1 )
         ELSE
            TEMP = Y( K-1 )
            Y( K-1 ) = Y( K )
            Y( K ) = TEMP - C( K-1 )*Y( K )
         END IF
      ENDDO
      IF( JOB == 1 ) THEN
         DO K = N, 1, -1
            IF( K <= N-2 ) THEN
               TEMP = Y( K ) - B( K )*Y( K+1 ) - D( K )*Y( K+2 )
            ELSE IF( K == N-1 ) THEN
               TEMP = Y( K ) - B( K )*Y( K+1 )
            ELSE
               TEMP = Y( K )
            END IF
            AK = A( K )
            ABSAK = ABS( AK )
            IF( ABSAK < ONE ) THEN
               IF( ABSAK < SFMIN ) THEN
                  IF( ABSAK == ZERO .OR. ABS( TEMP )*SFMIN > ABSAK ) &
                       THEN
                     INFO = K
                     RETURN
                  ELSE
                     TEMP = TEMP*BIGNUM
                     AK = AK*BIGNUM
                  END IF
               ELSE IF( ABS( TEMP ) > ABSAK*BIGNUM ) THEN
                  INFO = K
                  RETURN
               END IF
            END IF
            Y( K ) = TEMP / AK
         ENDDO
      ELSE
         DO K = N, 1, -1
            IF( K <= N-2 ) THEN
               TEMP = Y( K ) - B( K )*Y( K+1 ) - D( K )*Y( K+2 )
            ELSE IF( K == N-1 ) THEN
               TEMP = Y( K ) - B( K )*Y( K+1 )
            ELSE
               TEMP = Y( K )
            END IF
            AK = A( K )
            PERT = SIGN( TOL, AK )
40          CONTINUE
            ABSAK = ABS( AK )
            IF( ABSAK < ONE ) THEN
               IF( ABSAK < SFMIN ) THEN
                  IF( ABSAK == ZERO .OR. ABS( TEMP )*SFMIN > ABSAK ) &
                       THEN
                     AK = AK + PERT
                     PERT = 2*PERT
                     GO TO 40
                  ELSE
                     TEMP = TEMP*BIGNUM
                     AK = AK*BIGNUM
                  END IF
               ELSE IF( ABS( TEMP ) > ABSAK*BIGNUM ) THEN
                  AK = AK + PERT
                  PERT = 2*PERT
                  GO TO 40
               END IF
            END IF
            Y( K ) = TEMP / AK
         ENDDO
      END IF
   ELSE
!
!        Come to here if  JOB = 2 or -2
!
      IF( JOB == 2 ) THEN
         DO K = 1, N
            IF( K >= 3 ) THEN
               TEMP = Y( K ) - B( K-1 )*Y( K-1 ) - D( K-2 )*Y( K-2 )
            ELSE IF( K == 2 ) THEN
               TEMP = Y( K ) - B( K-1 )*Y( K-1 )
            ELSE
               TEMP = Y( K )
            END IF
            AK = A( K )
            ABSAK = ABS( AK )
            IF( ABSAK < ONE ) THEN
               IF( ABSAK < SFMIN ) THEN
                  IF( ABSAK == ZERO .OR. ABS( TEMP )*SFMIN > ABSAK ) &
                       THEN
                     INFO = K
                     RETURN
                  ELSE
                     TEMP = TEMP*BIGNUM
                     AK = AK*BIGNUM
                  END IF
               ELSE IF( ABS( TEMP ) > ABSAK*BIGNUM ) THEN
                  INFO = K
                  RETURN
               END IF
            END IF
            Y( K ) = TEMP / AK
         ENDDO
      ELSE
         DO K = 1, N
            IF( K >= 3 ) THEN
               TEMP = Y( K ) - B( K-1 )*Y( K-1 ) - D( K-2 )*Y( K-2 )
            ELSE IF( K == 2 ) THEN
               TEMP = Y( K ) - B( K-1 )*Y( K-1 )
            ELSE
               TEMP = Y( K )
            END IF
            AK = A( K )
            PERT = SIGN( TOL, AK )
70          CONTINUE
            ABSAK = ABS( AK )
            IF( ABSAK < ONE ) THEN
               IF( ABSAK < SFMIN ) THEN
                  IF( ABSAK == ZERO .OR. ABS( TEMP )*SFMIN > ABSAK ) &
                       THEN
                     AK = AK + PERT
                     PERT = 2*PERT
                     GO TO 70
                  ELSE
                     TEMP = TEMP*BIGNUM
                     AK = AK*BIGNUM
                  END IF
               ELSE IF( ABS( TEMP ) > ABSAK*BIGNUM ) THEN
                  AK = AK + PERT
                  PERT = 2*PERT
                  GO TO 70
               END IF
            END IF
            Y( K ) = TEMP / AK
         ENDDO
      END IF
!
      DO K = N, 2, -1
         IF( IN( K-1 ) == 0 ) THEN
            Y( K-1 ) = Y( K-1 ) - C( K-1 )*Y( K )
         ELSE
            TEMP = Y( K-1 )
            Y( K-1 ) = Y( K )
            Y( K ) = TEMP - C( K-1 )*Y( K )
         END IF
      ENDDO
   END IF
!
!     End of SLAGTS
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
