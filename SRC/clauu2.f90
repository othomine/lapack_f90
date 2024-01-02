!> \brief \b CLAUU2 computes the product UUH or LHL, where U and L are upper or lower triangular matrices (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAUU2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clauu2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clauu2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clauu2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAUU2( UPLO, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAUU2 computes the product U * U**H or L**H * L, where the triangular
!> factor U or L is stored in the upper or lower triangular part of
!> the array A.
!>
!> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
!> overwriting the factor U in A.
!> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
!> overwriting the factor L in A.
!>
!> This is the unblocked form of the algorithm, calling Level 2 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the triangular factor stored in the array A
!>          is upper or lower triangular:
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the triangular factor U or L.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the triangular factor U or L.
!>          On exit, if UPLO = 'U', the upper triangle of A is
!>          overwritten with the upper triangle of the product U * U**H;
!>          if UPLO = 'L', the lower triangle of A is overwritten with
!>          the lower triangle of the product L**H * L.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -k, the k-th argument had an illegal value
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
!> \ingroup lauu2
!
!  =====================================================================
   SUBROUTINE CLAUU2( UPLO, N, A, LDA, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
   COMPLEX            A( LDA, * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            UPPER
   INTEGER            I
   REAL               AII
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMV, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   UPPER = LSAME( UPLO, 'U' )
   IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -4
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CLAUU2', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) RETURN
!
   IF( UPPER ) THEN
!
!        Compute the product U * U**H.
!
      DO I = 1, N
         AII = REAL( A( I, I ) )
         IF( I < N ) THEN
            A( I, I ) = AII*AII + REAL(SUM(CONJG(A(I,I+1:N))*A(I,I+1:N)))
            A(I,I+1:N) = CONJG(A(I,I+1:N))
            CALL CGEMV( 'No transpose', I-1, N-I, (1.0E+0,0.0E+0), A( 1, I+1 ), &
                        LDA, A( I, I+1 ), LDA, CMPLX( AII ), &
                        A( 1, I ), 1 )
            A(I,I+1:N) = CONJG(A(I,I+1:N))
         ELSE
            A(1:I,I) = AII*A(1:I,I)
         END IF
      ENDDO
!
   ELSE
!
!        Compute the product L**H * L.
!
      DO I = 1, N
         AII = REAL( A( I, I ) )
         IF( I < N ) THEN
            A( I, I ) = AII*AII + REAL(SUM(CONJG(A(I+1:N,I))*A(I+1:N,I)))
            A(I,1:I-1) = CONJG(A(I,1:I-1))
            CALL CGEMV( 'Conjugate transpose', N-I, I-1, (1.0E+0,0.0E+0), &
                        A( I+1, 1 ), LDA, A( I+1, I ), 1, &
                        CMPLX( AII ), A( I, 1 ), LDA )
            A(I,1:I-1) = CONJG(A(I,1:I-1))
         ELSE
            
            A(I,1:I) = AII*A(I,1:I)
         END IF
      ENDDO
   END IF
!
   RETURN
!
!     End of CLAUU2
!
END
