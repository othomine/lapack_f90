!> \brief \b CHETRI_ROOK computes the inverse of HE matrix using the factorization obtained with the bounded Bunch-Kaufman ("rook") diagonal pivoting method.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHETRI_ROOK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetri_rook.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetri_rook.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetri_rook.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHETRI_ROOK( UPLO, N, A, LDA, IPIV, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHETRI_ROOK computes the inverse of a complex Hermitian indefinite matrix
!> A using the factorization A = U*D*U**H or A = L*D*L**H computed by
!> CHETRF_ROOK.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the details of the factorization are stored
!>          as an upper or lower triangular matrix.
!>          = 'U':  Upper triangular, form is A = U*D*U**H;
!>          = 'L':  Lower triangular, form is A = L*D*L**H.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the block diagonal matrix D and the multipliers
!>          used to obtain the factor U or L as computed by CHETRF_ROOK.
!>
!>          On exit, if INFO = 0, the (Hermitian) inverse of the original
!>          matrix.  If UPLO = 'U', the upper triangular part of the
!>          inverse is formed and the part of A below the diagonal is not
!>          referenced; if UPLO = 'L' the lower triangular part of the
!>          inverse is formed and the part of A above the diagonal is
!>          not referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D
!>          as determined by CHETRF_ROOK.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
!>               inverse could not be computed.
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
!> \ingroup hetri_rook
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>  November 2013,  Igor Kozachenko,
!>                  Computer Science Division,
!>                  University of California, Berkeley
!>
!>  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas,
!>                  School of Mathematics,
!>                  University of Manchester
!> \endverbatim
!
!  =====================================================================
   SUBROUTINE CHETRI_ROOK( UPLO, N, A, LDA, IPIV, WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * )
   COMPLEX            A( LDA, * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            UPPER
   INTEGER            J, K, KP, KSTEP
   REAL               AK, AKP1, D, T
   COMPLEX            AKKP1, TEMP, A_TMP( LDA )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   COMPLEX            CDOTC
   EXTERNAL           LSAME, CDOTC
!     ..
!     .. External Subroutines ..
   EXTERNAL           CHEMV, XERBLA
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
      CALL XERBLA( 'CHETRI_ROOK', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) RETURN
!
!     Check that the diagonal matrix D is nonsingular.
!
   IF( UPPER ) THEN
!
!        Upper triangular storage: examine D from bottom to top
!
      DO INFO = N, 1, -1
         IF( IPIV( INFO ) > 0 .AND. A( INFO, INFO ) == (0.0E+0,0.0E+0) ) RETURN
      ENDDO
   ELSE
!
!        Lower triangular storage: examine D from top to bottom.
!
      DO INFO = 1, N
         IF( IPIV( INFO ) > 0 .AND. A( INFO, INFO ) == (0.0E+0,0.0E+0) ) RETURN
      ENDDO
   END IF
   INFO = 0
!
   IF( UPPER ) THEN
!
!        Compute inv(A) from the factorization A = U*D*U**H.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
      K = 1
30    CONTINUE
!
!        If K > N, exit from loop.
!
      IF( K > N ) GO TO 70
!
      IF( IPIV( K ) > 0 ) THEN
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
         A( K, K ) = 1.0E+0 / REAL( A( K, K ) )
!
!           Compute column K of the inverse.
!
         IF( K > 1 ) THEN
            WORK(1:K-1) = A(1:K-1,K)
            CALL CHEMV( UPLO, K-1, -(1.0E+0,0.0E+0), A, LDA, WORK, 1, (0.0E+0,0.0E+0), &
                        A( 1, K ), 1 )
            A( K, K ) = A( K, K ) - REAL( CDOTC( K-1, WORK, 1, A( 1, K ), 1 ) )
         END IF
         KSTEP = 1
      ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
         T = ABS( A( K, K+1 ) )
         AK = REAL( A( K, K ) ) / T
         AKP1 = REAL( A( K+1, K+1 ) ) / T
         AKKP1 = A( K, K+1 ) / T
         D = T*( AK*AKP1-1.0E+0 )
         A( K, K ) = AKP1 / D
         A( K+1, K+1 ) = AK / D
         A( K, K+1 ) = -AKKP1 / D
!
!           Compute columns K and K+1 of the inverse.
!
         IF( K > 1 ) THEN
            WORK(1:K-1) = A(1:K-1,K)
            CALL CHEMV( UPLO, K-1, -(1.0E+0,0.0E+0), A, LDA, WORK, 1, (0.0E+0,0.0E+0), &
                        A( 1, K ), 1 )
            A( K, K ) = A( K, K ) - REAL( CDOTC( K-1, WORK, 1, A( 1, &
                        K ), 1 ) )
            A( K, K+1 ) = A( K, K+1 ) - &
                          CDOTC( K-1, A( 1, K ), 1, A( 1, K+1 ), 1 )
            WORK(1:K-1) = A(1:K-1,K+1)
            CALL CHEMV( UPLO, K-1, -(1.0E+0,0.0E+0), A, LDA, WORK, 1, (0.0E+0,0.0E+0), &
                        A( 1, K+1 ), 1 )
            A( K+1, K+1 ) = A( K+1, K+1 ) - &
                            REAL( CDOTC( K-1, WORK, 1, A( 1, K+1 ), 1 ) )
         END IF
         KSTEP = 2
      END IF
!
      IF( KSTEP == 1 ) THEN
!
!           Interchange rows and columns K and IPIV(K) in the leading
!           submatrix A(1:k,1:k)
!
         KP = IPIV( K )
         IF( KP /= K ) THEN
!
            IF( KP > 1 ) THEN
               A_TMP(1:KP-1) = A(1:KP-1,K)
               A(1:KP-1,K) = A(1:KP-1,KP)
               A(1:KP-1,KP) = A_TMP(1:KP-1)
            ENDIF
!
            A_TMP(1:K-KP-1) = CONJG(A(KP+1:K-1,K))
            A(KP+1:K-1,K) = CONJG(A(KP,KP+1:K-1))
            A(KP,KP+1:K-1) = A_TMP(1:K-KP-1)
!
            A( KP, K ) = CONJG( A( KP, K ) )
!
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
         END IF
      ELSE
!
!           Interchange rows and columns K and K+1 with -IPIV(K) and
!           -IPIV(K+1) in the leading submatrix A(k+1:n,k+1:n)
!
!           (1) Interchange rows and columns K and -IPIV(K)
!
         KP = -IPIV( K )
         IF( KP /= K ) THEN
!
            IF( KP > 1 ) THEN
               A_TMP(1:KP-1) = A(1:KP-1,K)
               A(1:KP-1,K) = A(1:KP-1,KP)
               A(1:KP-1,KP) = A_TMP(1:KP-1)
            ENDIF
!
            A_TMP(1:K-KP-1) = CONJG(A(KP+1:K-1,K))
            A(KP+1:K-1,K) = CONJG(A(KP,KP+1:K-1))
            A(KP,KP+1:K-1) = A_TMP(1:K-KP-1)
!
            A( KP, K ) = CONJG( A( KP, K ) )
!
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
!
            TEMP = A( K, K+1 )
            A( K, K+1 ) = A( KP, K+1 )
            A( KP, K+1 ) = TEMP
         END IF
!
!           (2) Interchange rows and columns K+1 and -IPIV(K+1)
!
         K = K + 1
         KP = -IPIV( K )
         IF( KP /= K ) THEN
!
            IF( KP > 1 ) THEN
               A_TMP(1:KP-1) = A(1:KP-1,K)
               A(1:KP-1,K) = A(1:KP-1,KP)
               A(1:KP-1,KP) = A_TMP(1:KP-1)
            ENDIF
!
            A_TMP(1:K-KP-1) = CONJG(A(KP+1:K-1,K))
            A(KP+1:K-1,K) = CONJG(A(KP,KP+1:K-1))
            A(KP,KP+1:K-1) = A_TMP(1:K-KP-1)
!
            A( KP, K ) = CONJG( A( KP, K ) )
!
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
         END IF
      END IF
!
      K = K + 1
      GO TO 30
70    CONTINUE
!
   ELSE
!
!        Compute inv(A) from the factorization A = L*D*L**H.
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
      K = N
80    CONTINUE
!
!        If K < 1, exit from loop.
!
      IF( K < 1 ) GO TO 120
!
      IF( IPIV( K ) > 0 ) THEN
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
         A( K, K ) = 1.0E+0 / REAL( A( K, K ) )
!
!           Compute column K of the inverse.
!
         IF( K < N ) THEN
            WORK(1:N-K) = A(K+1:N,K)
            CALL CHEMV( UPLO, N-K, -(1.0E+0,0.0E+0), A( K+1, K+1 ), LDA, WORK, &
                        1, (0.0E+0,0.0E+0), A( K+1, K ), 1 )
            A( K, K ) = A( K, K ) - REAL( CDOTC( N-K, WORK, 1, A( K+1, K ), 1 ) )
         END IF
         KSTEP = 1
      ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
         T = ABS( A( K, K-1 ) )
         AK = REAL( A( K-1, K-1 ) ) / T
         AKP1 = REAL( A( K, K ) ) / T
         AKKP1 = A( K, K-1 ) / T
         D = T*( AK*AKP1-1.0E+0 )
         A( K-1, K-1 ) = AKP1 / D
         A( K, K ) = AK / D
         A( K, K-1 ) = -AKKP1 / D
!
!           Compute columns K-1 and K of the inverse.
!
         IF( K < N ) THEN
            WORK(1:N-K) = A(K+1:N,K)
            CALL CHEMV( UPLO, N-K, -(1.0E+0,0.0E+0), A( K+1, K+1 ), LDA, WORK, &
                        1, (0.0E+0,0.0E+0), A( K+1, K ), 1 )
            A( K, K ) = A( K, K ) - REAL( CDOTC( N-K, WORK, 1, &
                        A( K+1, K ), 1 ) )
            A( K, K-1 ) = A( K, K-1 ) - &
                          CDOTC( N-K, A( K+1, K ), 1, A( K+1, K-1 ), 1 )
            WORK(1:N-K) = A(K+1:N,K-1)
            CALL CHEMV( UPLO, N-K, -(1.0E+0,0.0E+0), A( K+1, K+1 ), LDA, WORK, &
                        1, (0.0E+0,0.0E+0), A( K+1, K-1 ), 1 )
            A( K-1, K-1 ) = A( K-1, K-1 ) - &
                            REAL( CDOTC( N-K, WORK, 1, A( K+1, K-1 ), 1 ) )
         END IF
         KSTEP = 2
      END IF
!
      IF( KSTEP == 1 ) THEN
!
!           Interchange rows and columns K and IPIV(K) in the trailing
!           submatrix A(k:n,k:n)
!
         KP = IPIV( K )
         IF( KP /= K ) THEN
!
            IF( KP < N ) THEN
               A_TMP(1:N-KP) = A(KP+1:N,K)
               A(KP+1:N,K) = A(KP+1:N,KP)
               A(KP+1:N,KP) = A_TMP(1:N-KP)
            ENDIF
!
            A_TMP(1:KP-K-1) = CONJG(A(K+1:KP-1,K))
            A(K+1:KP-1,K) = CONJG(A(KP,K+1:KP-1))
            A(KP,K+1:KP-1) = A_TMP(1:KP-K-1)
!
            A( KP, K ) = CONJG( A( KP, K ) )
!
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
         END IF
      ELSE
!
!           Interchange rows and columns K and K-1 with -IPIV(K) and
!           -IPIV(K-1) in the trailing submatrix A(k-1:n,k-1:n)
!
!           (1) Interchange rows and columns K and -IPIV(K)
!
         KP = -IPIV( K )
         IF( KP /= K ) THEN
!
            IF( KP < N ) THEN
               A_TMP(1:N-KP) = A(KP+1:N,K)
               A(KP+1:N,K) = A(KP+1:N,KP)
               A(KP+1:N,KP) = A_TMP(1:N-KP)
            ENDIF
!
            A_TMP(1:KP-K-1) = CONJG(A(K+1:KP-1,K))
            A(K+1:KP-1,K) = CONJG(A(KP,K+1:KP-1))
            A(KP,K+1:KP-1) = A_TMP(1:KP-K-1)
!
            A( KP, K ) = CONJG( A( KP, K ) )
!
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
!
            TEMP = A( K, K-1 )
            A( K, K-1 ) = A( KP, K-1 )
            A( KP, K-1 ) = TEMP
         END IF
!
!           (2) Interchange rows and columns K-1 and -IPIV(K-1)
!
         K = K - 1
         KP = -IPIV( K )
         IF( KP /= K ) THEN
!
            IF( KP < N ) THEN
               A_TMP(1:N-KP) = A(KP+1:N,K)
               A(KP+1:N,K) = A(KP+1:N,KP)
               A(KP+1:N,KP) = A_TMP(1:N-KP)
            ENDIF
!
            A_TMP(1:KP-K-1) = CONJG(A(K+1:KP-1,K))
            A(K+1:KP-1,K) = CONJG(A(KP,K+1:KP-1))
            A(KP,K+1:KP-1) = A_TMP(1:KP-K-1)
!
            A( KP, K ) = CONJG( A( KP, K ) )
!
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
         END IF
      END IF
!
      K = K - 1
      GO TO 80
  120    CONTINUE
   END IF
!
   RETURN
!
!     End of CHETRI_ROOK
!
END

