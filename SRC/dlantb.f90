!> \brief \b DLANTB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a triangular band matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLANTB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlantb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlantb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlantb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DLANTB( NORM, UPLO, DIAG, N, K, AB,
!                        LDAB, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, NORM, UPLO
!       INTEGER            K, LDAB, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AB( LDAB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLANTB  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the element of  largest absolute value  of an
!> n by n triangular band matrix A,  with ( k + 1 ) diagonals.
!> \endverbatim
!>
!> \return DLANTB
!> \verbatim
!>
!>    DLANTB = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!>             (
!>             ( norm1(A),         NORM = '1', 'O' or 'o'
!>             (
!>             ( normI(A),         NORM = 'I' or 'i'
!>             (
!>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!>
!> where  norm1  denotes the  one norm of a matrix (maximum column sum),
!> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!> normF  denotes the  Frobenius norm of a matrix (square root of sum of
!> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NORM
!> \verbatim
!>          NORM is CHARACTER*1
!>          Specifies the value to be returned in DLANTB as described
!>          above.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A is upper or lower triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A is unit triangular.
!>          = 'N':  Non-unit triangular
!>          = 'U':  Unit triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, DLANTB is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of super-diagonals of the matrix A if UPLO = 'U',
!>          or the number of sub-diagonals of the matrix A if UPLO = 'L'.
!>          K >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
!>          The upper or lower triangular band matrix A, stored in the
!>          first k+1 rows of AB.  The j-th column of A is stored
!>          in the j-th column of the array AB as follows:
!>          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k).
!>          Note that when DIAG = 'U', the elements of the array AB
!>          corresponding to the diagonal elements of the matrix A are
!>          not referenced, but are assumed to be one.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= K+1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
!>          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
!>          referenced.
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
!> \ingroup lantb
!
!  =====================================================================
   DOUBLE PRECISION FUNCTION DLANTB( NORM, UPLO, DIAG, N, K, AB, &
                    LDAB, WORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          DIAG, NORM, UPLO
   INTEGER            K, LDAB, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   AB( LDAB, * ), WORK( * )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ONE, ZERO
   PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            UDIAG
   INTEGER            I, J, L
   DOUBLE PRECISION   SCALE, SUM, VALUE
!     ..
!     .. External Subroutines ..
   EXTERNAL           DLASSQ
!     ..
!     .. External Functions ..
   LOGICAL            LSAME, DISNAN
   EXTERNAL           LSAME, DISNAN
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
   IF( N == 0 ) THEN
      VALUE = ZERO
   ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
      IF( LSAME( DIAG, 'U' ) ) THEN
         VALUE = ONE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO J = 1, N
               DO I = MAX( K+2-J, 1 ), K
                  SUM = ABS( AB( I, J ) )
                  IF( VALUE  <  SUM .OR. DISNAN( SUM ) ) VALUE = SUM
               ENDDO
            ENDDO
         ELSE
            DO J = 1, N
               DO I = 2, MIN( N+1-J, K+1 )
                  SUM = ABS( AB( I, J ) )
                  IF( VALUE  <  SUM .OR. DISNAN( SUM ) ) VALUE = SUM
               ENDDO
            ENDDO
         END IF
      ELSE
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO J = 1, N
               DO I = MAX( K+2-J, 1 ), K + 1
                  SUM = ABS( AB( I, J ) )
                  IF( VALUE  <  SUM .OR. DISNAN( SUM ) ) VALUE = SUM
               ENDDO
            ENDDO
         ELSE
            DO J = 1, N
               DO I = 1, MIN( N+1-J, K+1 )
                  SUM = ABS( AB( I, J ) )
                  IF( VALUE  <  SUM .OR. DISNAN( SUM ) ) VALUE = SUM
               ENDDO
            ENDDO
         END IF
      END IF
   ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM == '1' ) ) THEN
!
!        Find norm1(A).
!
      VALUE = ZERO
      UDIAG = LSAME( DIAG, 'U' )
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO J = 1, N
            IF( UDIAG ) THEN
               SUM = ONE
               DO I = MAX( K+2-J, 1 ), K
                  SUM = SUM + ABS( AB( I, J ) )
               ENDDO
            ELSE
               SUM = ZERO
               DO I = MAX( K+2-J, 1 ), K + 1
                  SUM = SUM + ABS( AB( I, J ) )
                  ENDDO
            END IF
            IF( VALUE  <  SUM .OR. DISNAN( SUM ) ) VALUE = SUM
            ENDDO
      ELSE
         DO J = 1, N
            IF( UDIAG ) THEN
               SUM = ONE
               DO I = 2, MIN( N+1-J, K+1 )
                  SUM = SUM + ABS( AB( I, J ) )
                  ENDDO
            ELSE
               SUM = ZERO
               DO I = 1, MIN( N+1-J, K+1 )
                  SUM = SUM + ABS( AB( I, J ) )
                  ENDDO
            END IF
            IF( VALUE  <  SUM .OR. DISNAN( SUM ) ) VALUE = SUM
            ENDDO
      END IF
   ELSE IF( LSAME( NORM, 'I' ) ) THEN
!
!        Find normI(A).
!
      VALUE = ZERO
      IF( LSAME( UPLO, 'U' ) ) THEN
         IF( LSAME( DIAG, 'U' ) ) THEN
            DO I = 1, N
               WORK( I ) = ONE
               ENDDO
            DO J = 1, N
               L = K + 1 - J
               DO I = MAX( 1, J-K ), J - 1
                  WORK( I ) = WORK( I ) + ABS( AB( L+I, J ) )
                  ENDDO
               ENDDO
         ELSE
            DO I = 1, N
               WORK( I ) = ZERO
               ENDDO
            DO J = 1, N
               L = K + 1 - J
               DO I = MAX( 1, J-K ), J
                  WORK( I ) = WORK( I ) + ABS( AB( L+I, J ) )
                  ENDDO
               ENDDO
         END IF
      ELSE
         IF( LSAME( DIAG, 'U' ) ) THEN
            DO I = 1, N
               WORK( I ) = ONE
               ENDDO
            DO J = 1, N
               L = 1 - J
               DO I = J + 1, MIN( N, J+K )
                  WORK( I ) = WORK( I ) + ABS( AB( L+I, J ) )
                  ENDDO
               ENDDO
         ELSE
            DO I = 1, N
               WORK( I ) = ZERO
               ENDDO
            DO J = 1, N
               L = 1 - J
               DO I = J, MIN( N, J+K )
                  WORK( I ) = WORK( I ) + ABS( AB( L+I, J ) )
                  ENDDO
               ENDDO
         END IF
      END IF
      DO I = 1, N
         SUM = WORK( I )
         IF( VALUE  <  SUM .OR. DISNAN( SUM ) ) VALUE = SUM
         ENDDO
   ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
      IF( LSAME( UPLO, 'U' ) ) THEN
         IF( LSAME( DIAG, 'U' ) ) THEN
            SCALE = ONE
            SUM = N
            IF( K > 0 ) THEN
               DO J = 2, N
                  CALL DLASSQ( MIN( J-1, K ), &
                               AB( MAX( K+2-J, 1 ), J ), 1, SCALE, &
                               SUM )
                  ENDDO
            END IF
         ELSE
            SCALE = ZERO
            SUM = ONE
            DO J = 1, N
               CALL DLASSQ( MIN( J, K+1 ), AB( MAX( K+2-J, 1 ), J ), &
                            1, SCALE, SUM )
               ENDDO
         END IF
      ELSE
         IF( LSAME( DIAG, 'U' ) ) THEN
            SCALE = ONE
            SUM = N
            IF( K > 0 ) THEN
               DO J = 1, N - 1
                  CALL DLASSQ( MIN( N-J, K ), AB( 2, J ), 1, SCALE, &
                               SUM )
                  ENDDO
            END IF
         ELSE
            SCALE = ZERO
            SUM = ONE
            DO J = 1, N
               CALL DLASSQ( MIN( N-J+1, K+1 ), AB( 1, J ), 1, SCALE, &
                            SUM )
               ENDDO
         END IF
      END IF
      VALUE = SCALE*SQRT( SUM )
   END IF
!
   DLANTB = VALUE
   RETURN
!
!     End of DLANTB
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

