!> \brief \b CLASCL multiplies a general rectangular matrix by a real scalar defined as cto/cfrom.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLASCL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clascl.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clascl.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clascl.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TYPE
!       INTEGER            INFO, KL, KU, LDA, M, N
!       REAL               CFROM, CTO
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
!> CLASCL multiplies the M by N complex matrix A by the real scalar
!> CTO/CFROM.  This is done without over/underflow as long as the final
!> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
!> A may be full, upper triangular, lower triangular, upper Hessenberg,
!> or banded.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TYPE
!> \verbatim
!>          TYPE is CHARACTER*1
!>          TYPE indices the storage type of the input matrix.
!>          = 'G':  A is a full matrix.
!>          = 'L':  A is a lower triangular matrix.
!>          = 'U':  A is an upper triangular matrix.
!>          = 'H':  A is an upper Hessenberg matrix.
!>          = 'B':  A is a symmetric band matrix with lower bandwidth KL
!>                  and upper bandwidth KU and with the only the lower
!>                  half stored.
!>          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
!>                  and upper bandwidth KU and with the only the upper
!>                  half stored.
!>          = 'Z':  A is a band matrix with lower bandwidth KL and upper
!>                  bandwidth KU. See CGBTRF for storage details.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The lower bandwidth of A.  Referenced only if TYPE = 'B',
!>          'Q' or 'Z'.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The upper bandwidth of A.  Referenced only if TYPE = 'B',
!>          'Q' or 'Z'.
!> \endverbatim
!>
!> \param[in] CFROM
!> \verbatim
!>          CFROM is REAL
!> \endverbatim
!>
!> \param[in] CTO
!> \verbatim
!>          CTO is REAL
!>
!>          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
!>          without over/underflow if the final result CTO*A(I,J)/CFROM
!>          can be represented without over/underflow.  CFROM must be
!>          nonzero.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
!>          storage type.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          If TYPE = 'G', 'L', 'U', 'H', LDA >= max(1,M);
!>             TYPE = 'B', LDA >= KL+1;
!>             TYPE = 'Q', LDA >= KU+1;
!>             TYPE = 'Z', LDA >= 2*KL+KU+1.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          0  - successful exit
!>          <0 - if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup lascl
!
!  =====================================================================
   SUBROUTINE CLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          TYPE
   INTEGER            INFO, KL, KU, LDA, M, N
   REAL               CFROM, CTO
!     ..
!     .. Array Arguments ..
   COMPLEX            A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO, ONE
   PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            DONE
   INTEGER            I, ITYPE, J, K1, K2, K3, K4
   REAL               BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
!     ..
!     .. External Functions ..
   LOGICAL            LSAME, SISNAN
   REAL               SLAMCH
   EXTERNAL           LSAME, SLAMCH, SISNAN
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, MAX, MIN
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
   INFO = 0
!
   IF( LSAME( TYPE, 'G' ) ) THEN
      ITYPE = 0
   ELSE IF( LSAME( TYPE, 'L' ) ) THEN
      ITYPE = 1
   ELSE IF( LSAME( TYPE, 'U' ) ) THEN
      ITYPE = 2
   ELSE IF( LSAME( TYPE, 'H' ) ) THEN
      ITYPE = 3
   ELSE IF( LSAME( TYPE, 'B' ) ) THEN
      ITYPE = 4
   ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
      ITYPE = 5
   ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
      ITYPE = 6
   ELSE
      ITYPE = -1
   END IF
!
   IF( ITYPE == -1 ) THEN
      INFO = -1
   ELSE IF( CFROM == ZERO .OR. SISNAN(CFROM) ) THEN
      INFO = -4
   ELSE IF( SISNAN(CTO) ) THEN
      INFO = -5
   ELSE IF( M < 0 ) THEN
      INFO = -6
   ELSE IF( N < 0 .OR. ( ITYPE == 4 .AND. N /= M ) .OR. &
            ( ITYPE == 5 .AND. N /= M ) ) THEN
      INFO = -7
   ELSE IF( ITYPE <= 3 .AND. LDA < MAX( 1, M ) ) THEN
      INFO = -9
   ELSE IF( ITYPE >= 4 ) THEN
      IF( KL < 0 .OR. KL > MAX( M-1, 0 ) ) THEN
         INFO = -2
      ELSE IF( KU < 0 .OR. KU > MAX( N-1, 0 ) .OR. &
               ( ( ITYPE == 4 .OR. ITYPE == 5 ) .AND. KL /= KU ) ) &
                THEN
         INFO = -3
      ELSE IF( ( ITYPE == 4 .AND. LDA < KL+1 ) .OR. &
               ( ITYPE == 5 .AND. LDA < KU+1 ) .OR. &
               ( ITYPE == 6 .AND. LDA < 2*KL+KU+1 ) ) THEN
         INFO = -9
      END IF
   END IF
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CLASCL', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 .OR. M == 0 ) &
      RETURN
!
!     Get machine parameters
!
   SMLNUM = SLAMCH( 'S' )
   BIGNUM = ONE / SMLNUM
!
   CFROMC = CFROM
   CTOC = CTO
!
10 CONTINUE
   CFROM1 = CFROMC*SMLNUM
   IF( CFROM1 == CFROMC ) THEN
!        CFROMC is an inf.  Multiply by a correctly signed zero for
!        finite CTOC, or a NaN if CTOC is infinite.
      MUL = CTOC / CFROMC
      DONE = .TRUE.
      CTO1 = CTOC
   ELSE
      CTO1 = CTOC / BIGNUM
      IF( CTO1 == CTOC ) THEN
!           CTOC is either 0 or an inf.  In both cases, CTOC itself
!           serves as the correct multiplication factor.
         MUL = CTOC
         DONE = .TRUE.
         CFROMC = ONE
      ELSE IF( ABS( CFROM1 ) > ABS( CTOC ) .AND. CTOC /= ZERO ) THEN
         MUL = SMLNUM
         DONE = .FALSE.
         CFROMC = CFROM1
      ELSE IF( ABS( CTO1 ) > ABS( CFROMC ) ) THEN
         MUL = BIGNUM
         DONE = .FALSE.
         CTOC = CTO1
      ELSE
         MUL = CTOC / CFROMC
         DONE = .TRUE.
         IF (MUL  ==  ONE) &
            RETURN
      END IF
   END IF
!
   IF( ITYPE == 0 ) THEN
!
!        Full matrix
!
      DO J = 1, N
         DO I = 1, M
            A( I, J ) = A( I, J )*MUL
         ENDDO
      ENDDO
!
   ELSE IF( ITYPE == 1 ) THEN
!
!        Lower triangular matrix
!
      DO J = 1, N
         DO I = J, M
            A( I, J ) = A( I, J )*MUL
         ENDDO
      ENDDO
!
   ELSE IF( ITYPE == 2 ) THEN
!
!        Upper triangular matrix
!
      DO J = 1, N
         DO I = 1, MIN( J, M )
            A( I, J ) = A( I, J )*MUL
         ENDDO
      ENDDO
!
   ELSE IF( ITYPE == 3 ) THEN
!
!        Upper Hessenberg matrix
!
      DO J = 1, N
         DO I = 1, MIN( J+1, M )
            A( I, J ) = A( I, J )*MUL
         ENDDO
      ENDDO
!
   ELSE IF( ITYPE == 4 ) THEN
!
!        Lower half of a symmetric band matrix
!
      K3 = KL + 1
      K4 = N + 1
      DO J = 1, N
         DO I = 1, MIN( K3, K4-J )
            A( I, J ) = A( I, J )*MUL
            ENDDO
         ENDDO
!
   ELSE IF( ITYPE == 5 ) THEN
!
!        Upper half of a symmetric band matrix
!
      K1 = KU + 2
      K3 = KU + 1
      DO J = 1, N
         DO I = MAX( K1-J, 1 ), K3
            A( I, J ) = A( I, J )*MUL
            ENDDO
         ENDDO
!
   ELSE IF( ITYPE == 6 ) THEN
!
!        Band matrix
!
      K1 = KL + KU + 2
      K2 = KL + 1
      K3 = 2*KL + KU + 1
      K4 = KL + KU + 1 + M
      DO J = 1, N
         DO I = MAX( K1-J, K2 ), MIN( K3, K4-J )
            A( I, J ) = A( I, J )*MUL
            ENDDO
         ENDDO
!
   END IF
!
   IF( .NOT.DONE ) &
      GO TO 10
!
   RETURN
!
!     End of CLASCL
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

