!> \brief \b CLARFB applies a block reflector or its conjugate-transpose to a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLARFB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarfb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarfb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarfb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
!                          T, LDT, C, LDC, WORK, LDWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, SIDE, STOREV, TRANS
!       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
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
!> CLARFB applies a complex block reflector H or its transpose H**H to a
!> complex M-by-N matrix C, from either the left or the right.
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
!>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!> \endverbatim
!>
!> \param[in] STOREV
!> \verbatim
!>          STOREV is CHARACTER*1
!>          Indicates how the vectors which define the elementary
!>          reflectors are stored:
!>          = 'C': Columnwise
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
!>          If SIDE = 'L', M >= K >= 0;
!>          if SIDE = 'R', N >= K >= 0.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is COMPLEX array, dimension
!>                                (LDV,K) if STOREV = 'C'
!>                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!>                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!>          The matrix V. See Further Details.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V.
!>          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
!>          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
!>          if STOREV = 'R', LDV >= K.
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
!> \ingroup larfb
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The shape of the matrix V and the storage of the vectors which define
!>  the H(i) is best illustrated by the following example with n = 5 and
!>  k = 3. The elements equal to 1 are not stored; the corresponding
!>  array elements are modified but restored on exit. The rest of the
!>  array is not used.
!>
!>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!>
!>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!>                   ( v1  1    )                     (     1 v2 v2 v2 )
!>                   ( v1 v2  1 )                     (        1 v3 v3 )
!>                   ( v1 v2 v3 )
!>                   ( v1 v2 v3 )
!>
!>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!>
!>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!>                   (     1 v3 )
!>                   (        1 )
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE CLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, &
                      T, LDT, C, LDC, WORK, LDWORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          DIRECT, SIDE, STOREV, TRANS
   INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
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
   INTEGER            I, J
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMM, CTRMM
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
   IF( M <= 0 .OR. N <= 0 ) RETURN
!
   IF( LSAME( TRANS, 'N' ) ) THEN
      TRANST = 'C'
   ELSE
      TRANST = 'N'
   END IF
!
   IF( LSAME( STOREV, 'C' ) ) THEN
!
      IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1 )    (first K rows)
!                     ( V2 )
!           where  V1  is unit lower triangular.
!
         IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H**H * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
!
!              W := C1**H
!
            DO J = 1, K
               WORK(1:N,J) = CONJG(C(J,1:N))
            ENDDO
!
!              W := W * V1
!
            CALL CTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                        K, (1.0E+0,0.0E+0), V, LDV, WORK, LDWORK )
            IF( M > K ) THEN
!
!                 W := W + C2**H *V2
!
               CALL CGEMM( 'Conjugate transpose', 'No transpose', N, &
                           K, M-K, (1.0E+0,0.0E+0), C( K+1, 1 ), LDC, &
                           V( K+1, 1 ), LDV, (1.0E+0,0.0E+0), WORK, LDWORK )
            END IF
!
!              W := W * T**H  or  W * T
!
            CALL CTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                        (1.0E+0,0.0E+0), T, LDT, WORK, LDWORK )
!
!              C := C - V * W**H
!
            IF( M > K ) THEN
!
!                 C2 := C2 - V2 * W**H
!
               CALL CGEMM( 'No transpose', 'Conjugate transpose', &
                           M-K, N, K, -(1.0E+0,0.0E+0), V( K+1, 1 ), LDV, WORK, &
                           LDWORK, (1.0E+0,0.0E+0), C( K+1, 1 ), LDC )
            END IF
!
!              W := W * V1**H
!
            CALL CTRMM( 'Right', 'Lower', 'Conjugate transpose', &
                        'Unit', N, K, (1.0E+0,0.0E+0), V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W**H
!
            DO J = 1, K
               C(J,1:N) = C(J,1:N) - CONJG(WORK(1:N,J))
            ENDDO
!
         ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C1
!
            DO J = 1, K
               WORK(1:M,J) = C(1:M,J)
            ENDDO
!
!              W := W * V1
!
            CALL CTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                        K, (1.0E+0,0.0E+0), V, LDV, WORK, LDWORK )
            IF( N > K ) THEN
!
!                 W := W + C2 * V2
!
               CALL CGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                           (1.0E+0,0.0E+0), C( 1, K+1 ), LDC, V( K+1, 1 ), LDV, &
                           (1.0E+0,0.0E+0), WORK, LDWORK )
            END IF
!
!              W := W * T  or  W * T**H
!
            CALL CTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                        (1.0E+0,0.0E+0), T, LDT, WORK, LDWORK )
!
!              C := C - W * V**H
!
            IF( N > K ) THEN
!
!                 C2 := C2 - W * V2**H
!
               CALL CGEMM( 'No transpose', 'Conjugate transpose', M, &
                           N-K, K, -(1.0E+0,0.0E+0), WORK, LDWORK, V( K+1, 1 ), &
                           LDV, (1.0E+0,0.0E+0), C( 1, K+1 ), LDC )
            END IF
!
!              W := W * V1**H
!
            CALL CTRMM( 'Right', 'Lower', 'Conjugate transpose', &
                        'Unit', M, K, (1.0E+0,0.0E+0), V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
            C(1:M,1:K) = C(1:M,1:K) - WORK(1:M,1:K)
         END IF
!
      ELSE
!
!           Let  V =  ( V1 )
!                     ( V2 )    (last K rows)
!           where  V2  is unit upper triangular.
!
         IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H**H * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
!
!              W := C2**H
!
            DO J = 1, K
               WORK(1:N,J) = CONJG(C(M-K+J,1:N))
            ENDDO
!
!              W := W * V2
!
            CALL CTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                        K, (1.0E+0,0.0E+0), V( M-K+1, 1 ), LDV, WORK, LDWORK )
            IF( M > K ) THEN
!
!                 W := W + C1**H * V1
!
               CALL CGEMM( 'Conjugate transpose', 'No transpose', N, &
                           K, M-K, (1.0E+0,0.0E+0), C, LDC, V, LDV, (1.0E+0,0.0E+0), WORK, &
                           LDWORK )
            END IF
!
!              W := W * T**H  or  W * T
!
            CALL CTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                        (1.0E+0,0.0E+0), T, LDT, WORK, LDWORK )
!
!              C := C - V * W**H
!
            IF( M > K ) THEN
!
!                 C1 := C1 - V1 * W**H
!
               CALL CGEMM( 'No transpose', 'Conjugate transpose', &
                           M-K, N, K, -(1.0E+0,0.0E+0), V, LDV, WORK, LDWORK, &
                           (1.0E+0,0.0E+0), C, LDC )
            END IF
!
!              W := W * V2**H
!
            CALL CTRMM( 'Right', 'Upper', 'Conjugate transpose', &
                        'Unit', N, K, (1.0E+0,0.0E+0), V( M-K+1, 1 ), LDV, WORK, &
                        LDWORK )
!
!              C2 := C2 - W**H
!
            DO J = 1, K
               C(M-K+J,1:N) = C(M-K+J,1:N) - CONJG(WORK(1:N,J))
            ENDDO
!
         ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C2
!
            DO J = 1, K
               WORK(1:M,J) = C(1:M,N-K+J)
            ENDDO
!
!              W := W * V2
!
            CALL CTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                        K, (1.0E+0,0.0E+0), V( N-K+1, 1 ), LDV, WORK, LDWORK )
            IF( N > K ) THEN
!
!                 W := W + C1 * V1
!
               CALL CGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                           (1.0E+0,0.0E+0), C, LDC, V, LDV, (1.0E+0,0.0E+0), WORK, LDWORK )
            END IF
!
!              W := W * T  or  W * T**H
!
            CALL CTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                        (1.0E+0,0.0E+0), T, LDT, WORK, LDWORK )
!
!              C := C - W * V**H
!
            IF( N > K ) THEN
!
!                 C1 := C1 - W * V1**H
!
               CALL CGEMM( 'No transpose', 'Conjugate transpose', M, &
                           N-K, K, -(1.0E+0,0.0E+0), WORK, LDWORK, V, LDV, (1.0E+0,0.0E+0), &
                           C, LDC )
            END IF
!
!              W := W * V2**H
!
            CALL CTRMM( 'Right', 'Upper', 'Conjugate transpose', &
                        'Unit', M, K, (1.0E+0,0.0E+0), V( N-K+1, 1 ), LDV, WORK, &
                        LDWORK )
!
!              C2 := C2 - W
!
            C(1:M,N-K+1:N ) = C(1:M,N-K+1:N) - WORK(1:M,1:K)
         END IF
      END IF
!
   ELSE IF( LSAME( STOREV, 'R' ) ) THEN
!
      IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1  V2 )    (V1: first K columns)
!           where  V1  is unit upper triangular.
!
         IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H**H * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
!
!              W := C1**H
!
            DO J = 1, K
               WORK(1:N,J) = CONJG(C(J,1:N))
            ENDDO
!
!              W := W * V1**H
!
            CALL CTRMM( 'Right', 'Upper', 'Conjugate transpose', &
                        'Unit', N, K, (1.0E+0,0.0E+0), V, LDV, WORK, LDWORK )
            IF( M > K ) THEN
!
!                 W := W + C2**H * V2**H
!
               CALL CGEMM( 'Conjugate transpose', &
                           'Conjugate transpose', N, K, M-K, (1.0E+0,0.0E+0), &
                           C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, (1.0E+0,0.0E+0), &
                           WORK, LDWORK )
            END IF
!
!              W := W * T**H  or  W * T
!
            CALL CTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                        (1.0E+0,0.0E+0), T, LDT, WORK, LDWORK )
!
!              C := C - V**H * W**H
!
            IF( M > K ) THEN
!
!                 C2 := C2 - V2**H * W**H
!
               CALL CGEMM( 'Conjugate transpose', &
                           'Conjugate transpose', M-K, N, K, -(1.0E+0,0.0E+0), &
                           V( 1, K+1 ), LDV, WORK, LDWORK, (1.0E+0,0.0E+0), &
                           C( K+1, 1 ), LDC )
            END IF
!
!              W := W * V1
!
            CALL CTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                        K, (1.0E+0,0.0E+0), V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W**H
!
            DO J = 1, K
               C(J,1:N) = C(J,1:N) - CONJG(WORK(1:N,J))
            ENDDO
!
         ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
!              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
!
!              W := C1
!
            DO J = 1, K
               WORK(1:M,J) = C(1:M,J)
            ENDDO
!
!              W := W * V1**H
!
            CALL CTRMM( 'Right', 'Upper', 'Conjugate transpose', &
                        'Unit', M, K, (1.0E+0,0.0E+0), V, LDV, WORK, LDWORK )
            IF( N > K ) THEN
!
!                 W := W + C2 * V2**H
!
               CALL CGEMM( 'No transpose', 'Conjugate transpose', M, &
                           K, N-K, (1.0E+0,0.0E+0), C( 1, K+1 ), LDC, &
                           V( 1, K+1 ), LDV, (1.0E+0,0.0E+0), WORK, LDWORK )
            END IF
!
!              W := W * T  or  W * T**H
!
            CALL CTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                        (1.0E+0,0.0E+0), T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
            IF( N > K ) THEN
!
!                 C2 := C2 - W * V2
!
               CALL CGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                           -(1.0E+0,0.0E+0), WORK, LDWORK, V( 1, K+1 ), LDV, (1.0E+0,0.0E+0), &
                           C( 1, K+1 ), LDC )
            END IF
!
!              W := W * V1
!
            CALL CTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                        K, (1.0E+0,0.0E+0), V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
            C(1:M,1:K) = C(1:M,1:K) - WORK(1:M,1:K)
!
         END IF
!
      ELSE
!
!           Let  V =  ( V1  V2 )    (V2: last K columns)
!           where  V2  is unit lower triangular.
!
         IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H**H * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
!
!              W := C2**H
!
            DO J = 1, K
               WORK(1:N,J) = CONJG(C(M-K+J,1:N))
            ENDDO
!
!              W := W * V2**H
!
            CALL CTRMM( 'Right', 'Lower', 'Conjugate transpose', &
                        'Unit', N, K, (1.0E+0,0.0E+0), V( 1, M-K+1 ), LDV, WORK, &
                        LDWORK )
            IF( M > K ) THEN
!
!                 W := W + C1**H * V1**H
!
               CALL CGEMM( 'Conjugate transpose', &
                           'Conjugate transpose', N, K, M-K, (1.0E+0,0.0E+0), C, &
                           LDC, V, LDV, (1.0E+0,0.0E+0), WORK, LDWORK )
            END IF
!
!              W := W * T**H  or  W * T
!
            CALL CTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                        (1.0E+0,0.0E+0), T, LDT, WORK, LDWORK )
!
!              C := C - V**H * W**H
!
            IF( M > K ) THEN
!
!                 C1 := C1 - V1**H * W**H
!
               CALL CGEMM( 'Conjugate transpose', &
                           'Conjugate transpose', M-K, N, K, -(1.0E+0,0.0E+0), V, &
                           LDV, WORK, LDWORK, (1.0E+0,0.0E+0), C, LDC )
            END IF
!
!              W := W * V2
!
            CALL CTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                        K, (1.0E+0,0.0E+0), V( 1, M-K+1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W**H
!
            DO J = 1, K
               C(M-K+J,1:N) = C(M-K+J,1:N) - CONJG(WORK(1:N,J))
            ENDDO
!
         ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
!              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
!
!              W := C2
!
            DO J = 1, K
               WORK(1:M,J) = C(1:M,N-K+J)
            ENDDO
!
!              W := W * V2**H
!
            CALL CTRMM( 'Right', 'Lower', 'Conjugate transpose', &
                        'Unit', M, K, (1.0E+0,0.0E+0), V( 1, N-K+1 ), LDV, WORK, &
                        LDWORK )
            IF( N > K ) THEN
!
!                 W := W + C1 * V1**H
!
               CALL CGEMM( 'No transpose', 'Conjugate transpose', M, &
                           K, N-K, (1.0E+0,0.0E+0), C, LDC, V, LDV, (1.0E+0,0.0E+0), WORK, &
                           LDWORK )
            END IF
!
!              W := W * T  or  W * T**H
!
            CALL CTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                        (1.0E+0,0.0E+0), T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
            IF( N > K ) THEN
!
!                 C1 := C1 - W * V1
!
               CALL CGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                           -(1.0E+0,0.0E+0), WORK, LDWORK, V, LDV, (1.0E+0,0.0E+0), C, LDC )
            END IF
!
!              W := W * V2
!
            CALL CTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                        K, (1.0E+0,0.0E+0), V( 1, N-K+1 ), LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
            C(1:M,N-K+1:N) = C(1:M,N-K+1:N) - WORK(1:M,1:K)
!
         END IF
!
      END IF
   END IF
!
   RETURN
!
!     End of CLARFB
!
END
