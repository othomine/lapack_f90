!> \brief \b CSPTRI
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CSPTRI + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csptri.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csptri.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csptri.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            AP( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSPTRI computes the inverse of a complex symmetric indefinite matrix
!> A in packed storage using the factorization A = U*D*U**T or
!> A = L*D*L**T computed by CSPTRF.
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
!> \param[in,out] AP
!> \verbatim
!>          AP is COMPLEX array, dimension (N*(N+1)/2)
!>          On entry, the block diagonal matrix D and the multipliers
!>          used to obtain the factor U or L as computed by CSPTRF,
!>          stored as a packed triangular matrix.
!>
!>          On exit, if INFO = 0, the (symmetric) inverse of the original
!>          matrix, stored as a packed triangular matrix. The j-th column
!>          of inv(A) is stored in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = inv(A)(i,j) for 1<=i<=j;
!>          if UPLO = 'L',
!>             AP(i + (j-1)*(2n-j)/2) = inv(A)(i,j) for j<=i<=n.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D
!>          as determined by CSPTRF.
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
!> \ingroup hptri
!
!  =====================================================================
   SUBROUTINE CSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
   IMPLICIT NONE
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * )
   COMPLEX            AP( * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..

!     .. Local Array ..
   COMPLEX            AP_tmp(N*(N+1)/2)
!     .. Local Scalars ..
   LOGICAL            UPPER
   INTEGER            J, K, KC, KCNEXT, KP, KPC, KSTEP, KX, NPP
   COMPLEX            AK, AKKP1, AKP1, uoD, uoT, TEMP
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           CSPMV, XERBLA
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
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CSPTRI', -INFO )
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
      KP = N*( N+1 ) / 2
      DO INFO = N, 1, -1
         IF( IPIV( INFO ) > 0 .AND. AP( KP ) == (0.0E+0,0.0E+0 ) ) RETURN
         KP = KP - INFO
      ENDDO
   ELSE
!
!        Lower triangular storage: examine D from top to bottom.
!
      KP = 1
      DO INFO = 1, N
         IF( IPIV( INFO ) > 0 .AND. AP( KP ) == (0.0E+0,0.0E+0 ) ) RETURN
         KP = KP + N - INFO + 1
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
      KC = 1
30    CONTINUE
!
!        If K > N, exit from loop.
!
      IF( K > N ) GO TO 50
!
      KCNEXT = KC + K
      IF( IPIV( K ) > 0 ) THEN
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
         AP( KC+K-1 ) = (1.0E+0,0.0E+0 ) / AP( KC+K-1 )
!
!           Compute column K of the inverse.
!
         IF( K > 1 ) THEN
            WORK(1:K-1) = AP(KC:KC+K-2)
            CALL CSPMV( UPLO, K-1, -(1.0E+0,0.0E+0 ), AP, WORK, 1, (0.0E+0,0.0E+0 ), AP( KC ), 1 )
            AP( KC+K-1 ) = AP( KC+K-1 ) - sum(AP(KC:KC+K-2)*WORK(1:K-1))
         END IF
         KSTEP = 1
      ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
         uoT = (1.0E+0,0.0E+0)/AP( KCNEXT+K-1 )
         AK = AP( KC+K-1 ) * uoT
         AKP1 = AP( KCNEXT+K ) * uoT
         AKKP1 = AP( KCNEXT+K-1 ) * uoT
         uoD = uoT/(AK*AKP1-(1.0E+0,0.0E+0))
         AP( KC+K-1 ) = AKP1 * uoD
         AP( KCNEXT+K ) = AK * uoD
         AP( KCNEXT+K-1 ) = -AKKP1 * uoD
!
!           Compute columns K and K+1 of the inverse.
!
         IF( K > 1 ) THEN
            WORK(1:K-1) = AP(KC:KC+K-2)
            CALL CSPMV( UPLO, K-1, -(1.0E+0,0.0E+0 ), AP, WORK, 1, (0.0E+0,0.0E+0 ), AP( KC ), &
                        1 )
            AP( KC+K-1 ) = AP( KC+K-1 ) - sum(WORK(1:K-1)*AP(KC:KC+K-2))
            AP( KCNEXT+K-1 ) = AP( KCNEXT+K-1 ) - sum(AP(KC:KC+K-2)*AP(KCNEXT:KCNEXT+K-2))
            WORK(1:K-1) = AP(KCNEXT:KCNEXT+K-2)
            CALL CSPMV( UPLO, K-1, -(1.0E+0,0.0E+0 ), AP, WORK, 1, (0.0E+0,0.0E+0 ), AP( KCNEXT ), 1 )
            AP( KCNEXT+K ) = AP( KCNEXT+K ) - sum(WORK(1:K-1)*AP(KCNEXT:KCNEXT+K-2))
         END IF
         KSTEP = 2
         KCNEXT = KCNEXT + K + 1
      END IF
!
      KP = ABS( IPIV( K ) )
      IF( KP /= K ) THEN
!
!           Interchange rows and columns K and KP in the leading
!           submatrix A(1:k+1,1:k+1)
!
         KPC = ( KP-1 )*KP / 2 + 1
         AP_tmp(1:KP-1) = AP(KC:KC+KP-2)
         AP(KC:KC+KP-2) = AP(KPC:KPC+KP-2)
         AP(KPC:KPC+KP-2) = AP_tmp(1:KP-1)
         KX = KPC + KP - 1
         DO J = KP + 1, K - 1
            KX = KX + J - 1
            TEMP = AP( KC+J-1 )
            AP( KC+J-1 ) = AP( KX )
            AP( KX ) = TEMP
         ENDDO
         TEMP = AP( KC+K-1 )
         AP( KC+K-1 ) = AP( KPC+KP-1 )
         AP( KPC+KP-1 ) = TEMP
         IF( KSTEP == 2 ) THEN
            TEMP = AP( KC+K+K-1 )
            AP( KC+K+K-1 ) = AP( KC+K+KP-1 )
            AP( KC+K+KP-1 ) = TEMP
         END IF
      END IF
!
      K = K + KSTEP
      KC = KCNEXT
      GO TO 30
50    CONTINUE
!
   ELSE
!
!        Compute inv(A) from the factorization A = L*D*L**T.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
      NPP = N*( N+1 ) / 2
      K = N
      KC = NPP
60    CONTINUE
!
!        If K < 1, exit from loop.
!
      IF( K < 1 ) GO TO 80
!
      KCNEXT = KC - ( N-K+2 )
      IF( IPIV( K ) > 0 ) THEN
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
         AP( KC ) = (1.0E+0,0.0E+0 ) / AP( KC )
!
!           Compute column K of the inverse.
!
         IF( K < N ) THEN
            WORK(1:N-K) = AP(KC+1:KC+N-K)
            CALL CSPMV( UPLO, N-K, -(1.0E+0,0.0E+0 ), AP(KC+N-K+1), WORK, 1, &
                        (0.0E+0,0.0E+0 ), AP( KC+1 ), 1 )
            AP( KC ) = AP( KC ) - sum(AP(KC+1:KC+N-K)*WORK(1:N-K))
         END IF
         KSTEP = 1
      ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
         uoT = (1.0E+0,0.0E+0)/AP( KCNEXT+1 )
         AK = AP( KCNEXT ) * uoT
         AKP1 = AP( KC ) * uoT
         AKKP1 = AP( KCNEXT+1 ) * uoT
         uoD = uoT/(AK*AKP1-(1.0E+0,0.0E+0 ))
         AP( KCNEXT ) = AKP1 * uoD
         AP( KC ) = AK * uoD
         AP( KCNEXT+1 ) = -AKKP1 * uoD
!
!           Compute columns K-1 and K of the inverse.
!
         IF( K < N ) THEN
            WORK(1:N-K) = AP(KC+1:KC+N-K)
            CALL CSPMV( UPLO, N-K, -(1.0E+0,0.0E+0 ), AP( KC+( N-K+1 ) ), WORK, 1, &
                        (0.0E+0,0.0E+0 ), AP( KC+1 ), 1 )
            AP( KC ) = AP( KC ) - sum(WORK(1:N-K)*AP(KC+1:KC+N-K))
            AP( KCNEXT+1 ) = AP( KCNEXT+1 ) - sum(AP(KC+1:KC+N-K)*AP(KCNEXT+2:KCNEXT+1+N-K))
            WORK(1:N-K) = AP(KCNEXT+2:KCNEXT+1+N-K)
            CALL CSPMV( UPLO, N-K, -(1.0E+0,0.0E+0 ), AP( KC+( N-K+1 ) ), WORK, 1, &
                        (0.0E+0,0.0E+0 ), AP( KCNEXT+2 ), 1 )
            AP( KCNEXT ) = AP( KCNEXT ) - sum(WORK(1:N-K)*AP(KCNEXT+2:KCNEXT+1+N-K))
         END IF
         KSTEP = 2
         KCNEXT = KCNEXT - ( N-K+3 )
      END IF
!
      KP = ABS( IPIV( K ) )
      IF( KP /= K ) THEN
!
!           Interchange rows and columns K and KP in the trailing
!           submatrix A(k-1:n,k-1:n)
!
         KPC = NPP - ( N-KP+1 )*( N-KP+2 ) / 2 + 1
         IF( KP < N ) THEN
           AP_TMP(1:N-KP) = AP(KC+KP-K+1:KC-K+N)
           AP(KC+KP-K+1:KC-K+N) = AP(KPC+1:KPC+N-KP)
           AP(KPC+1:KPC+N-KP) = AP_TMP(1:N-KP)
         ENDIF
         KX = KC + KP - K
         DO J = K + 1, KP - 1
            KX = KX + N - J + 1
            TEMP = AP( KC+J-K )
            AP( KC+J-K ) = AP( KX )
            AP( KX ) = TEMP
         ENDDO
         TEMP = AP( KC )
         AP( KC ) = AP( KPC )
         AP( KPC ) = TEMP
         IF( KSTEP == 2 ) THEN
            TEMP = AP( KC-N+K-1 )
            AP( KC-N+K-1 ) = AP( KC-N+KP-1 )
            AP( KC-N+KP-1 ) = TEMP
         END IF
      END IF
!
      K = K - KSTEP
      KC = KCNEXT
      GO TO 60
80    CONTINUE
   END IF
!
   RETURN
!
!     End of CSPTRI
!
END
