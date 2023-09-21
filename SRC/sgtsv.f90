!> \brief <b> SGTSV computes the solution to system of linear equations A * X = B for GT matrices </b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGTSV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgtsv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgtsv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgtsv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       REAL               B( LDB, * ), D( * ), DL( * ), DU( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGTSV  solves the equation
!>
!>    A*X = B,
!>
!> where A is an n by n tridiagonal matrix, by Gaussian elimination with
!> partial pivoting.
!>
!> Note that the equation  A**T*X = B  may be solved by interchanging the
!> order of the arguments DU and DL.
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
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in,out] DL
!> \verbatim
!>          DL is REAL array, dimension (N-1)
!>          On entry, DL must contain the (n-1) sub-diagonal elements of
!>          A.
!>
!>          On exit, DL is overwritten by the (n-2) elements of the
!>          second super-diagonal of the upper triangular matrix U from
!>          the LU factorization of A, in DL(1), ..., DL(n-2).
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          On entry, D must contain the diagonal elements of A.
!>
!>          On exit, D is overwritten by the n diagonal elements of U.
!> \endverbatim
!>
!> \param[in,out] DU
!> \verbatim
!>          DU is REAL array, dimension (N-1)
!>          On entry, DU must contain the (n-1) super-diagonal elements
!>          of A.
!>
!>          On exit, DU is overwritten by the (n-1) elements of the first
!>          super-diagonal of U.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,NRHS)
!>          On entry, the N by NRHS matrix of right hand side matrix B.
!>          On exit, if INFO = 0, the N by NRHS solution matrix X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, U(i,i) is exactly zero, and the solution
!>               has not been computed.  The factorization has not been
!>               completed unless i = N.
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
!> \ingroup gtsv
!
!  =====================================================================
   SUBROUTINE SGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
!
!  -- LAPACK driver routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
   REAL               B( LDB, * ), D( * ), DL( * ), DU( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO
   PARAMETER          ( ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, J
   REAL               FACT, TEMP
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, MAX
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
   INFO = 0
   IF( N < 0 ) THEN
      INFO = -1
   ELSE IF( NRHS < 0 ) THEN
      INFO = -2
   ELSE IF( LDB < MAX( 1, N ) ) THEN
      INFO = -7
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'SGTSV ', -INFO )
      RETURN
   END IF
!
   IF( N == 0 ) &
      RETURN
!
   IF( NRHS == 1 ) THEN
      DO I = 1, N - 2
         IF( ABS( D( I ) ) >= ABS( DL( I ) ) ) THEN
!
!              No row interchange required
!
            IF( D( I ) /= ZERO ) THEN
               FACT = DL( I ) / D( I )
               D( I+1 ) = D( I+1 ) - FACT*DU( I )
               B( I+1, 1 ) = B( I+1, 1 ) - FACT*B( I, 1 )
            ELSE
               INFO = I
               RETURN
            END IF
            DL( I ) = ZERO
         ELSE
!
!              Interchange rows I and I+1
!
            FACT = D( I ) / DL( I )
            D( I ) = DL( I )
            TEMP = D( I+1 )
            D( I+1 ) = DU( I ) - FACT*TEMP
            DL( I ) = DU( I+1 )
            DU( I+1 ) = -FACT*DL( I )
            DU( I ) = TEMP
            TEMP = B( I, 1 )
            B( I, 1 ) = B( I+1, 1 )
            B( I+1, 1 ) = TEMP - FACT*B( I+1, 1 )
         END IF
      ENDDO
      IF( N > 1 ) THEN
         I = N - 1
         IF( ABS( D( I ) ) >= ABS( DL( I ) ) ) THEN
            IF( D( I ) /= ZERO ) THEN
               FACT = DL( I ) / D( I )
               D( I+1 ) = D( I+1 ) - FACT*DU( I )
               B( I+1, 1 ) = B( I+1, 1 ) - FACT*B( I, 1 )
            ELSE
               INFO = I
               RETURN
            END IF
         ELSE
            FACT = D( I ) / DL( I )
            D( I ) = DL( I )
            TEMP = D( I+1 )
            D( I+1 ) = DU( I ) - FACT*TEMP
            DU( I ) = TEMP
            TEMP = B( I, 1 )
            B( I, 1 ) = B( I+1, 1 )
            B( I+1, 1 ) = TEMP - FACT*B( I+1, 1 )
         END IF
      END IF
      IF( D( N ) == ZERO ) THEN
         INFO = N
         RETURN
      END IF
   ELSE
      DO I = 1, N - 2
         IF( ABS( D( I ) ) >= ABS( DL( I ) ) ) THEN
!
!              No row interchange required
!
            IF( D( I ) /= ZERO ) THEN
               FACT = DL( I ) / D( I )
               D( I+1 ) = D( I+1 ) - FACT*DU( I )
               DO J = 1, NRHS
                  B( I+1, J ) = B( I+1, J ) - FACT*B( I, J )
               ENDDO
            ELSE
               INFO = I
               RETURN
            END IF
            DL( I ) = ZERO
         ELSE
!
!              Interchange rows I and I+1
!
            FACT = D( I ) / DL( I )
            D( I ) = DL( I )
            TEMP = D( I+1 )
            D( I+1 ) = DU( I ) - FACT*TEMP
            DL( I ) = DU( I+1 )
            DU( I+1 ) = -FACT*DL( I )
            DU( I ) = TEMP
            DO J = 1, NRHS
               TEMP = B( I, J )
               B( I, J ) = B( I+1, J )
               B( I+1, J ) = TEMP - FACT*B( I+1, J )
            ENDDO
         END IF
      ENDDO
      IF( N > 1 ) THEN
         I = N - 1
         IF( ABS( D( I ) ) >= ABS( DL( I ) ) ) THEN
            IF( D( I ) /= ZERO ) THEN
               FACT = DL( I ) / D( I )
               D( I+1 ) = D( I+1 ) - FACT*DU( I )
               DO J = 1, NRHS
                  B( I+1, J ) = B( I+1, J ) - FACT*B( I, J )
               ENDDO
            ELSE
               INFO = I
               RETURN
            END IF
         ELSE
            FACT = D( I ) / DL( I )
            D( I ) = DL( I )
            TEMP = D( I+1 )
            D( I+1 ) = DU( I ) - FACT*TEMP
            DU( I ) = TEMP
            DO J = 1, NRHS
               TEMP = B( I, J )
               B( I, J ) = B( I+1, J )
               B( I+1, J ) = TEMP - FACT*B( I+1, J )
            ENDDO
         END IF
      END IF
      IF( D( N ) == ZERO ) THEN
         INFO = N
         RETURN
      END IF
   END IF
!
!     Back solve with the matrix U from the factorization.
!
   IF( NRHS <= 2 ) THEN
      J = 1
70    CONTINUE
      B( N, J ) = B( N, J ) / D( N )
      IF( N > 1 ) &
         B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 )
      DO I = N - 2, 1, -1
         B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )* &
                     B( I+2, J ) ) / D( I )
      ENDDO
      IF( J < NRHS ) THEN
         J = J + 1
         GO TO 70
      END IF
   ELSE
      DO J = 1, NRHS
         B( N, J ) = B( N, J ) / D( N )
         IF( N > 1 ) &
            B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / &
                          D( N-1 )
         DO I = N - 2, 1, -1
            B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )* &
                        B( I+2, J ) ) / D( I )
         ENDDO
         ENDDO
   END IF
!
   RETURN
!
!     End of SGTSV
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        