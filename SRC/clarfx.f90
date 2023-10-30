!> \brief \b CLARFX applies an elementary reflector to a general rectangular matrix, with loop unrolling when the reflector has order ≤ 10.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLARFX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarfx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarfx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarfx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE
!       INTEGER            LDC, M, N
!       COMPLEX            TAU
!       ..
!       .. Array Arguments ..
!       COMPLEX            C( LDC, * ), V( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLARFX applies a complex elementary reflector H to a complex m by n
!> matrix C, from either the left or the right. H is represented in the
!> form
!>
!>       H = I - tau * v * v**H
!>
!> where tau is a complex scalar and v is a complex vector.
!>
!> If tau = 0, then H is taken to be the unit matrix
!>
!> This version uses inline code if H has order < 11.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': form  H * C
!>          = 'R': form  C * H
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
!> \param[in] V
!> \verbatim
!>          V is COMPLEX array, dimension (M) if SIDE = 'L'
!>                                        or (N) if SIDE = 'R'
!>          The vector v in the representation of H.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX
!>          The value tau in the representation of H.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC,N)
!>          On entry, the m by n matrix C.
!>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!>          or C * H if SIDE = 'R'.
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
!>          WORK is COMPLEX array, dimension (N) if SIDE = 'L'
!>                                            or (M) if SIDE = 'R'
!>          WORK is not referenced if H has order < 11.
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
!> \ingroup larfx
!
!  =====================================================================
   SUBROUTINE CLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          SIDE
   INTEGER            LDC, M, N
   COMPLEX            TAU
!     ..
!     .. Array Arguments ..
   COMPLEX            C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            J
   COMPLEX            SOMME, TT( 10 ), VV( 10 )
   COMPLEX            T1, T10, T2, T3, T4, T5, T6, T7, T8, T9, &
                      V1, V10, V2, V3, V4, V5, V6, V7, V8, V9
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           CLARF
!     ..
!     .. Executable Statements ..
!
   IF( TAU == (0.0E+0,0.0E+0) ) RETURN
   IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  H * C, where H has order m.
!
      IF (M == 1) THEN
!
!        Special code for 1 x 1 Householder
!
         C(1,1:N) = ((1.0E+0,0.0E+0)-TAU*V(1)*CONJG(V(1)))*C(1,1:N)
      ELSEIF (M <= 10) THEN
!
!        Special code for M x M Householder (1 < M <= 10)
!
         VV(1:M) = CONJG(V(1:M))
         TT(1:M) = TAU*CONJG(VV(1:M))
         DO J = 1, N
            SOMME = SUM(VV(1:M)*C(1:M,J))
            C(1:M,J) = C(1:M,J) - SOMME*TT(1:M)
         ENDDO
      ELSE
!
!        Code for general M, serialized
!
         CALL CLARF( SIDE, M, N, V, 1, TAU, C, LDC, WORK )
      ENDIF
   ELSE
!
!        Form  C * H, where H has order n.
!
      IF (N == 1) THEN
!
!        Special code for 1 x 1 Householder
!
         C(1:M,1) = ((1.0E+0,0.0E+0)-TAU*V(1)*CONJG(V(1)))*C(1:M,1)
      ELSEIF (N <= 10) THEN
!
!        Special code for M x M Householder (1 < M <= 10)
!
         TT(1:N) = TAU*CONJG(V(1:N))
         DO J = 1, M
            SOMME = SUM(V(1:N)*C(J,1:N))
            C(J,1:N) = C(J,1:N) - SOMME*TT(1:N)
         ENDDO
      ELSE
!
!        Code for general M, serialized
!
         CALL CLARF( SIDE, M, N, V, 1, TAU, C, LDC, WORK )
      ENDIF
   END IF
   RETURN
!
!     End of CLARFX
!
END
