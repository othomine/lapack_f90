!> \brief \b CLARZB applies a block reflector or its conjugate-transpose to a general matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLARZB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarzb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarzb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarzb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLARZB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, V,
!                          LDV, T, LDT, C, LDC, WORK, LDWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, SIDE, STOREV, TRANS
!       INTEGER            K, L, LDC, LDT, LDV, LDWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            C( LDC, * ), T( LDT, * ), V( LDV, * ),
!      $                   WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLARZB applies a complex block reflector H or its transpose H**H
!> to a complex distributed M-by-N  C from the left or the right.
!>
!> Currently, only STOREV = 'R' and DIRECT = 'B' are supported.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply H or H**H from the Left
!>          = 'R': apply H or H**H from the Right
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': apply H (No transpose)
!>          = 'C': apply H**H (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] DIRECT
!> \verbatim
!>          DIRECT is CHARACTER*1
!>          Indicates how H is formed from a product of elementary
!>          reflectors
!>          = 'F': H = H(1) H(2) . . . H(k) (Forward, not supported yet)
!>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!> \endverbatim
!>
!> \param[in] STOREV
!> \verbatim
!>          STOREV is CHARACTER*1
!>          Indicates how the vectors which define the elementary
!>          reflectors are stored:
!>          = 'C': Columnwise                        (not supported yet)
!>          = 'R': Rowwise
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The order of the matrix T (= the number of elementary
!>          reflectors whose product defines the block reflector).
!> \endverbatim
!>
!> \param[in] L
!> \verbatim
!>          L is INTEGER
!>          The number of columns of the matrix V containing the
!>          meaningful part of the Householder reflectors.
!>          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is COMPLEX array, dimension (LDV,NV).
!>          If STOREV = 'C', NV = K; if STOREV = 'R', NV = L.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V.
!>          If STOREV = 'C', LDV >= L; if STOREV = 'R', LDV >= K.
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is COMPLEX array, dimension (LDT,K)
!>          The triangular K-by-K matrix T in the representation of the
!>          block reflector.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= K.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by H*C or H**H*C or C*H or C*H**H.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (LDWORK,K)
!> \endverbatim
!>
!> \param[in] LDWORK
!> \verbatim
!>          LDWORK is INTEGER
!>          The leading dimension of the array WORK.
!>          If SIDE = 'L', LDWORK >= max(1,N);
!>          if SIDE = 'R', LDWORK >= max(1,M).
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
!> \ingroup larzb
!
!> \par Contributors:
!  ==================
!>
!>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE CLARZB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, V, &
                      LDV, T, LDT, C, LDC, WORK, LDWORK )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          DIRECT, SIDE, STOREV, TRANS
   INTEGER            K, L, LDC, LDT, LDV, LDWORK, M, N
!     ..
!     .. Array Arguments ..
   COMPLEX            C( LDC, * ), T( LDT, * ), V( LDV, * ), &
                      WORK( LDWORK, * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   CHARACTER          TRANST
   INTEGER            I, INFO, J
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMM, CTRMM, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
   IF( M <= 0 .OR. N <= 0 ) RETURN
!
!     Check for currently supported options
!
   INFO = 0
   IF( .NOT.LSAME( DIRECT, 'B' ) ) THEN
      INFO = -3
   ELSE IF( .NOT.LSAME( STOREV, 'R' ) ) THEN
      INFO = -4
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CLARZB', -INFO )
      RETURN
   END IF
!
   IF( LSAME( TRANS, 'N' ) ) THEN
      TRANST = 'C'
   ELSE
      TRANST = 'N'
   END IF
!
   IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  H * C  or  H**H * C
!
!        W( 1:n, 1:k ) = C( 1:k, 1:n )**H
!
      DO J = 1, K
         WORK(1:N,J) = C(J,1:N)
      ENDDO
!
!        W( 1:n, 1:k ) = W( 1:n, 1:k ) + ...
!                        C( m-l+1:m, 1:n )**H * V( 1:k, 1:l )**T
!
      IF( L > 0 ) &
         CALL CGEMM( 'Transpose', 'Conjugate transpose', N, K, L, &
                     (1.0E+0,0.0E+0), C( M-L+1, 1 ), LDC, V, LDV, (1.0E+0,0.0E+0), WORK, &
                     LDWORK )
!
!        W( 1:n, 1:k ) = W( 1:n, 1:k ) * T**T  or  W( 1:m, 1:k ) * T
!
      CALL CTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, (1.0E+0,0.0E+0), T, &
                  LDT, WORK, LDWORK )
!
!        C( 1:k, 1:n ) = C( 1:k, 1:n ) - W( 1:n, 1:k )**H
!
      DO J = 1, N
         C(1:K,J) = C(1:K,J) - WORK( J,1:K)
      ENDDO
!
!        C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
!                            V( 1:k, 1:l )**H * W( 1:n, 1:k )**H
!
      IF( L > 0 ) &
         CALL CGEMM( 'Transpose', 'Transpose', L, N, K, -(1.0E+0,0.0E+0), V, LDV, &
                     WORK, LDWORK, (1.0E+0,0.0E+0), C( M-L+1, 1 ), LDC )
!
   ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!        Form  C * H  or  C * H**H
!
!        W( 1:m, 1:k ) = C( 1:m, 1:k )
!
      DO J = 1, K
         WORK(1:M,J) = C(1:M,J)
      ENDDO
!
!        W( 1:m, 1:k ) = W( 1:m, 1:k ) + ...
!                        C( 1:m, n-l+1:n ) * V( 1:k, 1:l )**H
!
      IF( L > 0 ) &
         CALL CGEMM( 'No transpose', 'Transpose', M, K, L, (1.0E+0,0.0E+0), &
                     C( 1, N-L+1 ), LDC, V, LDV, (1.0E+0,0.0E+0), WORK, LDWORK )
!
!        W( 1:m, 1:k ) = W( 1:m, 1:k ) * conjg( T )  or
!                        W( 1:m, 1:k ) * T**H
!
      DO J = 1, K
         T(J:K,J) = CONJG(T(J:K,J))
      ENDDO
      CALL CTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, (1.0E+0,0.0E+0), T, &
                  LDT, WORK, LDWORK )
      DO J = 1, K
         T(J:K,J) = CONJG(T(J:K,J))
      ENDDO
!
!        C( 1:m, 1:k ) = C( 1:m, 1:k ) - W( 1:m, 1:k )
!
      C(1:M,1:K) = C(1:M,1:K) - WORK(1:M,1:K)
!
!        C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
!                            W( 1:m, 1:k ) * conjg( V( 1:k, 1:l ) )
!
      DO J = 1, L
         V(1:K,J) = CONJG(V(1:K,J))
      ENDDO
      IF( L > 0 ) &
         CALL CGEMM( 'No transpose', 'No transpose', M, L, K, -(1.0E+0,0.0E+0), &
                     WORK, LDWORK, V, LDV, (1.0E+0,0.0E+0), C( 1, N-L+1 ), LDC )
      DO J = 1, L
         V(1:K,J) = CONJG(V(1:K,J))
      ENDDO
!
   END IF
!
   RETURN
!
!     End of CLARZB
!
END
