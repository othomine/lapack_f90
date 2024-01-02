!> \brief \b CSYTRI
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CSYTRI + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytri.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytri.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytri.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
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
!> CSYTRI computes the inverse of a complex symmetric indefinite matrix
!> A using the factorization A = U*D*U**T or A = L*D*L**T computed by
!> CSYTRF.
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
!>          = 'U':  Upper triangular, form is A = U*D*U**T;
!>          = 'L':  Lower triangular, form is A = L*D*L**T.
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
!>          used to obtain the factor U or L as computed by CSYTRF.
!>
!>          On exit, if INFO = 0, the (symmetric) inverse of the original
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
!>          as determined by CSYTRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (2*N)
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
!> \ingroup hetri
!
!  =====================================================================
   SUBROUTINE CSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
   IMPLICIT NONE
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
!     .. Local Array ..
   COMPLEX            A_tmp(LDA)
!     .. Local Scalars ..
   LOGICAL            UPPER
   INTEGER            K, KP, KSTEP
   COMPLEX            AK, AKKP1, AKP1, uoD, uoT, TEMP
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           CSWAP, CSYMV, XERBLA
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
      CALL XERBLA( 'CSYTRI', -INFO )
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
!        Compute inv(A) from the factorization A = U*D*U**T.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
      K = 1
30    CONTINUE
!
!        If K > N, exit from loop.
!
      IF( K > N ) GO TO 40
!
      IF( IPIV( K ) > 0 ) THEN
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
         A( K, K ) = (1.0E+0,0.0E+0) / A( K, K )
!
!           Compute column K of the inverse.
!
         IF( K > 1 ) THEN
            WORK(1:K-1) = A(1:K-1,K)
            CALL CSYMV( UPLO, K-1, -(1.0E+0,0.0E+0), A, LDA, WORK, 1, (0.0E+0,0.0E+0), A( 1, K ), 1 )
            A( K, K ) = A( K, K ) - SUM(WORK(1:K-1)*A(1:K-1,K))
         END IF
         KSTEP = 1
      ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
         uoT = (1.0E+0,0.0E+0) / A( K, K+1 )
         AK = A( K, K ) * uoT
         AKP1 = A( K+1, K+1 ) * uoT
         AKKP1 = A( K, K+1 ) * uoT
         uoD = uoT/( AK*AKP1-(1.0E+0,0.0E+0) )
         A( K, K ) = AKP1 * uoD
         A( K+1, K+1 ) = AK * uoD
         A( K, K+1 ) = -AKKP1 * uoD
!
!           Compute columns K and K+1 of the inverse.
!
         IF( K > 1 ) THEN
            WORK(1:K-1) = A(1:K-1,K)
            CALL CSYMV( UPLO, K-1, -(1.0E+0,0.0E+0), A, LDA, WORK, 1, (0.0E+0,0.0E+0), &
                        A( 1, K ), 1 )
            A( K, K ) = A( K, K ) - SUM(WORK(1:K-1)*A(1:K-1,K))
            A( K, K+1 ) = A( K, K+1 ) - SUM(A(1:K-1,K)*A(1:K-1,K+1))
            WORK(1:K-1) = A(1:K-1,K+1)
            CALL CSYMV( UPLO, K-1, -(1.0E+0,0.0E+0), A, LDA, WORK, 1, (0.0E+0,0.0E+0), &
                        A( 1, K+1 ), 1 )
            A( K+1, K+1 ) = A( K+1, K+1 ) - SUM(WORK(1:K-1)*A(1:K-1,K+1))
         END IF
         KSTEP = 2
      END IF
!
      KP = ABS( IPIV( K ) )
      IF( KP /= K ) THEN
!
!           Interchange rows and columns K and KP in the leading
!           submatrix A(1:k+1,1:k+1)
!
         A_tmp(1:KP-1) = A(1:KP-1,K)
         A(1:KP-1,K) = A(1:KP-1,KP)
         A(1:KP-1,KP) = A_tmp(1:KP-1)
         CALL CSWAP( K-KP-1, A( KP+1, K ), 1, A( KP, KP+1 ), LDA )
         TEMP = A( K, K )
         A( K, K ) = A( KP, KP )
         A( KP, KP ) = TEMP
         IF( KSTEP == 2 ) THEN
            TEMP = A( K, K+1 )
            A( K, K+1 ) = A( KP, K+1 )
            A( KP, K+1 ) = TEMP
         END IF
      END IF
!
      K = K + KSTEP
      GO TO 30
40    CONTINUE
!
   ELSE
!
!        Compute inv(A) from the factorization A = L*D*L**T.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
      K = N
50    CONTINUE
!
!        If K < 1, exit from loop.
!
      IF( K < 1 ) GO TO 60
!
      IF( IPIV( K ) > 0 ) THEN
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
         A( K, K ) = (1.0E+0,0.0E+0) / A( K, K )
!
!           Compute column K of the inverse.
!
         IF( K < N ) THEN
            WORK(1:N-K) = A(K+1:N,K)
            CALL CSYMV( UPLO, N-K,-(1.0E+0,0.0E+0), A( K+1, K+1 ), LDA, WORK, 1, &
                        (0.0E+0,0.0E+0), A( K+1, K ), 1 )
            A( K, K ) = A( K, K ) - SUM(WORK(1:N-K)*A(K+1:N,K))
         END IF
         KSTEP = 1
      ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
         uoT = (1.0E+0,0.0E+0) / A( K, K-1 )
         AK = A( K-1, K-1 ) * uoT
         AKP1 = A( K, K ) * uoT
         AKKP1 = A( K, K-1 ) * uoT
         uoD = uoT/( AK*AKP1-(1.0E+0,0.0E+0) )
         A( K-1, K-1 ) = AKP1 * uoD
         A( K, K ) = AK * uoD
         A( K, K-1 ) = -AKKP1 * uoD
!
!           Compute columns K-1 and K of the inverse.
!
         IF( K < N ) THEN
            WORK(1:N-K) = A(K+1:N,K)
            CALL CSYMV( UPLO, N-K,-(1.0E+0,0.0E+0), A( K+1, K+1 ), LDA, WORK, 1, &
                        (0.0E+0,0.0E+0), A( K+1, K ), 1 )
            A( K, K ) = A( K, K ) - SUM(WORK(1:N-K)*A(K+1:N,K))
            A( K, K-1 ) = A( K, K-1 ) - SUM(A(K+1:N,K)*A(K+1:N,K-1))
            WORK(1:N-K) = A(K+1:N,K-1)
            CALL CSYMV( UPLO, N-K,-(1.0E+0,0.0E+0), A( K+1, K+1 ), LDA, WORK, 1, &
                        (0.0E+0,0.0E+0), A( K+1, K-1 ), 1 )
            A( K-1, K-1 ) = A( K-1, K-1 ) - SUM(WORK(1:N-K)*A(K+1:N,K-1))
         END IF
         KSTEP = 2
      END IF
!
      KP = ABS( IPIV( K ) )
      IF( KP /= K ) THEN
!
!           Interchange rows and columns K and KP in the trailing
!           submatrix A(k-1:n,k-1:n)
!
         IF( KP < N ) THEN
            A_tmp(1:N-KP) = A(KP+1:N,K)
            A(KP+1:N,K) = A(KP+1:N,KP)
            A(KP+1:N,KP) = A_tmp(1:N-KP)
         ENDIF
         CALL CSWAP( KP-K-1, A( K+1, K ), 1, A( KP, K+1 ), LDA )
         TEMP = A( K, K )
         A( K, K ) = A( KP, KP )
         A( KP, KP ) = TEMP
         IF( KSTEP == 2 ) THEN
            TEMP = A( K, K-1 )
            A( K, K-1 ) = A( KP, K-1 )
            A( KP, K-1 ) = TEMP
         END IF
      END IF
!
      K = K - KSTEP
      GO TO 50
60    CONTINUE
   END IF
!
   RETURN
!
!     End of CSYTRI
!
END
