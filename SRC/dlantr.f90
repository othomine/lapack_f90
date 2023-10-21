!> \brief \b DLANTR returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a trapezoidal or triangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLANTR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlantr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlantr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlantr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DLANTR( NORM, UPLO, DIAG, M, N, A, LDA,
!                        WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, NORM, UPLO
!       INTEGER            LDA, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLANTR  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> trapezoidal or triangular matrix A.
!> \endverbatim
!>
!> \return DLANTR
!> \verbatim
!>
!>    DLANTR = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in DLANTR as described
!>          above.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A is upper or lower trapezoidal.
!>          = 'U':  Upper trapezoidal
!>          = 'L':  Lower trapezoidal
!>          Note that A is triangular instead of trapezoidal if M = N.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A has unit diagonal.
!>          = 'N':  Non-unit diagonal
!>          = 'U':  Unit diagonal
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0, and if
!>          UPLO = 'U', M <= N.  When M = 0, DLANTR is set to zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0, and if
!>          UPLO = 'L', N <= M.  When N = 0, DLANTR is set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The trapezoidal matrix A (A is triangular if M = N).
!>          If UPLO = 'U', the leading m by n upper trapezoidal part of
!>          the array A contains the upper trapezoidal matrix, and the
!>          strictly lower triangular part of A is not referenced.
!>          If UPLO = 'L', the leading m by n lower trapezoidal part of
!>          the array A contains the lower trapezoidal matrix, and the
!>          strictly upper triangular part of A is not referenced.  Note
!>          that when DIAG = 'U', the diagonal elements of A are not
!>          referenced and are assumed to be one.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(M,1).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
!>          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
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
!> \ingroup lantr
!
!  =====================================================================
   DOUBLE PRECISION FUNCTION DLANTR( NORM, UPLO, DIAG, M, N, A, LDA, &
                    WORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          DIAG, NORM, UPLO
   INTEGER            LDA, M, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   A( LDA, * ), WORK( * )
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
   INTEGER            I, J
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
   INTRINSIC          ABS, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
   IF( MIN( M, N ) == 0 ) THEN
      VALUE = ZERO
   ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
      IF( LSAME( DIAG, 'U' ) ) THEN
         VALUE = ONE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO J = 1, N
               DO I = 1, MIN( M, J-1 )
                  SUM = ABS( A( I, J ) )
                  IF( VALUE  <  SUM .OR. DISNAN( SUM ) ) VALUE = SUM
               ENDDO
            ENDDO
         ELSE
            DO J = 1, N
               DO I = J + 1, M
                  SUM = ABS( A( I, J ) )
                  IF( VALUE  <  SUM .OR. DISNAN( SUM ) ) VALUE = SUM
               ENDDO
            ENDDO
         END IF
      ELSE
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO J = 1, N
               DO I = 1, MIN( M, J )
                  SUM = ABS( A( I, J ) )
                  IF( VALUE  <  SUM .OR. DISNAN( SUM ) ) VALUE = SUM
               ENDDO
            ENDDO
         ELSE
            DO J = 1, N
               DO I = J, M
                  SUM = ABS( A( I, J ) )
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
            IF( ( UDIAG ) .AND. ( J <= M ) ) THEN
               SUM = ONE
               DO I = 1, J - 1
                  SUM = SUM + ABS( A( I, J ) )
               ENDDO
            ELSE
               SUM = ZERO
               DO I = 1, MIN( M, J )
                  SUM = SUM + ABS( A( I, J ) )
                  ENDDO
            END IF
            IF( VALUE  <  SUM .OR. DISNAN( SUM ) ) VALUE = SUM
            ENDDO
      ELSE
         DO J = 1, N
            IF( UDIAG ) THEN
               SUM = ONE
               DO I = J + 1, M
                  SUM = SUM + ABS( A( I, J ) )
                  ENDDO
            ELSE
               SUM = ZERO
               DO I = J, M
                  SUM = SUM + ABS( A( I, J ) )
                  ENDDO
            END IF
            IF( VALUE  <  SUM .OR. DISNAN( SUM ) ) VALUE = SUM
            ENDDO
      END IF
   ELSE IF( LSAME( NORM, 'I' ) ) THEN
!
!        Find normI(A).
!
      IF( LSAME( UPLO, 'U' ) ) THEN
         IF( LSAME( DIAG, 'U' ) ) THEN
            DO I = 1, M
               WORK( I ) = ONE
               ENDDO
            DO J = 1, N
               DO I = 1, MIN( M, J-1 )
                  WORK( I ) = WORK( I ) + ABS( A( I, J ) )
                  ENDDO
               ENDDO
         ELSE
            DO I = 1, M
               WORK( I ) = ZERO
               ENDDO
            DO J = 1, N
               DO I = 1, MIN( M, J )
                  WORK( I ) = WORK( I ) + ABS( A( I, J ) )
                  ENDDO
               ENDDO
         END IF
      ELSE
         IF( LSAME( DIAG, 'U' ) ) THEN
            DO I = 1, MIN( M, N )
               WORK( I ) = ONE
               ENDDO
            DO I = N + 1, M
               WORK( I ) = ZERO
               ENDDO
            DO J = 1, N
               DO I = J + 1, M
                  WORK( I ) = WORK( I ) + ABS( A( I, J ) )
                  ENDDO
               ENDDO
         ELSE
            DO I = 1, M
               WORK( I ) = ZERO
               ENDDO
            DO J = 1, N
               DO I = J, M
                  WORK( I ) = WORK( I ) + ABS( A( I, J ) )
                  ENDDO
               ENDDO
         END IF
      END IF
      VALUE = ZERO
      DO I = 1, M
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
            SUM = MIN( M, N )
            DO J = 2, N
               CALL DLASSQ( MIN( M, J-1 ), A( 1, J ), 1, SCALE, SUM )
               ENDDO
         ELSE
            SCALE = ZERO
            SUM = ONE
            DO J = 1, N
               CALL DLASSQ( MIN( M, J ), A( 1, J ), 1, SCALE, SUM )
               ENDDO
         END IF
      ELSE
         IF( LSAME( DIAG, 'U' ) ) THEN
            SCALE = ONE
            SUM = MIN( M, N )
            DO J = 1, N
               CALL DLASSQ( M-J, A( MIN( M, J+1 ), J ), 1, SCALE, &
                            SUM )
               ENDDO
         ELSE
            SCALE = ZERO
            SUM = ONE
            DO J = 1, N
               CALL DLASSQ( M-J+1, A( J, J ), 1, SCALE, SUM )
               ENDDO
         END IF
      END IF
      VALUE = SCALE*SQRT( SUM )
   END IF
!
   DLANTR = VALUE
   RETURN
!
!     End of DLANTR
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

