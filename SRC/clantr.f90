!> \brief \b CLANTR returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a trapezoidal or triangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLANTR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clantr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clantr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clantr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION CLANTR( NORM, UPLO, DIAG, M, N, A, LDA,
!                        WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, NORM, UPLO
!       INTEGER            LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL               WORK( * )
!       COMPLEX            A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLANTR  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> trapezoidal or triangular matrix A.
!> \endverbatim
!>
!> \return CLANTR
!> \verbatim
!>
!>    CLANTR = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in CLANTR as described
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
!>          UPLO = 'U', M <= N.  When M = 0, CLANTR is set to zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0, and if
!>          UPLO = 'L', N <= M.  When N = 0, CLANTR is set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
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
!>          WORK is REAL array, dimension (MAX(1,LWORK)),
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
   REAL             FUNCTION CLANTR( NORM, UPLO, DIAG, M, N, A, LDA, WORK )
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
   REAL               WORK( * )
   COMPLEX            A( LDA, * )
!     ..
!
! =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            UDIAG
   INTEGER            I, J
   REAL               SCALE, SOMME, VALUE
!     ..
!     .. External Functions ..
   LOGICAL            LSAME, SISNAN
   EXTERNAL           LSAME, SISNAN
!     ..
!     .. External Subroutines ..
   EXTERNAL           CLASSQ
!     ..
!     .. Executable Statements ..
!
   IF( MIN( M, N ) == 0 ) THEN
      VALUE = 0.0E+0
   ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
      IF( LSAME( DIAG, 'U' ) ) THEN
         VALUE = 1.0E+0
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO J = 1, N
               DO I = 1, MIN( M, J-1 )
                  SOMME = ABS( A( I, J ) )
                  IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
               ENDDO
            ENDDO
         ELSE
            DO J = 1, N
               DO I = J + 1, M
                  SOMME = ABS( A( I, J ) )
                  IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
               ENDDO
            ENDDO
         END IF
      ELSE
         VALUE = 0.0E+0
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO J = 1, N
               DO I = 1, MIN( M, J )
                  SOMME = ABS( A( I, J ) )
                  IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
               ENDDO
            ENDDO
         ELSE
            DO J = 1, N
               DO I = J, M
                  SOMME = ABS( A( I, J ) )
                  IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
               ENDDO
            ENDDO
         END IF
      END IF
   ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM == '1' ) ) THEN
!
!        Find norm1(A).
!
      VALUE = 0.0E+0
      UDIAG = LSAME( DIAG, 'U' )
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO J = 1, N
            IF( ( UDIAG ) .AND. ( J <= M ) ) THEN
               SOMME = 1.0E+0 + SUM(ABS(A(1:J-1,J)))
            ELSE
               SOMME = SUM(ABS(A(1:MIN(M,J),J)))
            END IF
            IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
         ENDDO
      ELSE
         DO J = 1, N
            IF( UDIAG ) THEN
               SOMME = 1.0E+0 + SUM(ABS(A(J+1:M,J)))
            ELSE
               SOMME = SUM(ABS(A(J:M,J)))
            END IF
            IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
         ENDDO
      END IF
   ELSE IF( LSAME( NORM, 'I' ) ) THEN
!
!        Find normI(A).
!
      IF( LSAME( UPLO, 'U' ) ) THEN
         IF( LSAME( DIAG, 'U' ) ) THEN
            WORK(1:M) = 1.0E+0
            DO J = 1, N
               WORK(1:MIN(M,J-1)) = WORK(1:MIN(M,J-1)) + ABS(A(1:MIN(M,J-1),J))
            ENDDO
         ELSE
            WORK(1:M) = 0.0E+0
            DO J = 1, N
               WORK(1:MIN(M,J)) = WORK(1:MIN(M,J)) + ABS(A(1:MIN(M,J),J))
            ENDDO
         END IF
      ELSE
         IF( LSAME( DIAG, 'U' ) ) THEN
            DO I = 1, MIN( M, N )
               WORK( I ) = 1.0E+0
            ENDDO
            WORK(N+1:M) = 0.0E+0
            DO J = 1, N
               WORK(J+1:M) = WORK(J+1:M) + ABS( A(J+1:M, J ) )
            ENDDO
         ELSE
            WORK(1:M) = 0.0E+0
            DO J = 1, N
               WORK(J:M) = WORK(J:M) + ABS( A(J:M, J ) )
            ENDDO
         END IF
      END IF
      VALUE = 0.0E+0
      DO I = 1, M
         SOMME = WORK( I )
         IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
      ENDDO
   ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
      IF( LSAME( UPLO, 'U' ) ) THEN
         IF( LSAME( DIAG, 'U' ) ) THEN
            SCALE = 1.0E+0
            SOMME = MIN( M, N )
            DO J = 2, N
               CALL CLASSQ( MIN( M, J-1 ), A( 1, J ), 1, SCALE, SOMME )
            ENDDO
         ELSE
            SCALE = 0.0E+0
            SOMME = 1.0E+0
            DO J = 1, N
               CALL CLASSQ( MIN( M, J ), A( 1, J ), 1, SCALE, SOMME )
            ENDDO
         END IF
      ELSE
         IF( LSAME( DIAG, 'U' ) ) THEN
            SCALE = 1.0E+0
            SOMME = MIN( M, N )
            DO J = 1, N
               CALL CLASSQ( M-J, A( MIN( M, J+1 ), J ), 1, SCALE, SOMME )
            ENDDO
         ELSE
            SCALE = 0.0E+0
            SOMME = 1.0E+0
            DO J = 1, N
               CALL CLASSQ( M-J+1, A( J, J ), 1, SCALE, SOMME )
            ENDDO
         END IF
      END IF
      VALUE = SCALE*SQRT( SOMME )
   END IF
!
   CLANTR = VALUE
   RETURN
!
!     End of CLANTR
!
END
