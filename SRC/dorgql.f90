!> \brief \b DORGQL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DORGQL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgql.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgql.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgql.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DORGQL generates an M-by-N real matrix Q with orthonormal columns,
!> which is defined as the last N columns of a product of K elementary
!> reflectors of order M
!>
!>       Q  =  H(k) . . . H(2) H(1)
!>
!> as returned by DGEQLF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix Q. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix Q. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. N >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the (n-k+i)-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by DGEQLF in the last k columns of its array
!>          argument A.
!>          On exit, the M-by-N matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The first dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DGEQLF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= max(1,N).
!>          For optimum performance LWORK >= N*NB, where NB is the
!>          optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument has an illegal value
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
!> \ingroup ungql
!
!  =====================================================================
   SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO
   PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            LQUERY
   INTEGER            I, IB, IINFO, IWS, J, KK, L, LDWORK, LWKOPT, &
                      NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
   EXTERNAL           DLARFB, DLARFT, DORG2L, XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
   INTEGER            ILAENV
   EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
   INFO = 0
   LQUERY = ( LWORK == -1 )
   IF( M < 0 ) THEN
      INFO = -1
   ELSE IF( N < 0 .OR. N > M ) THEN
      INFO = -2
   ELSE IF( K < 0 .OR. K > N ) THEN
      INFO = -3
   ELSE IF( LDA < MAX( 1, M ) ) THEN
      INFO = -5
   END IF
!
   IF( INFO == 0 ) THEN
      IF( N == 0 ) THEN
         LWKOPT = 1
      ELSE
         NB = ILAENV( 1, 'DORGQL', ' ', M, N, K, -1 )
         LWKOPT = N*NB
      END IF
      WORK( 1 ) = LWKOPT
!
      IF( LWORK < MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
   END IF
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'DORGQL', -INFO )
      RETURN
   ELSE IF( LQUERY ) THEN
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N <= 0 ) THEN
      RETURN
   END IF
!
   NBMIN = 2
   NX = 0
   IWS = N
   IF( NB > 1 .AND. NB < K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
      NX = MAX( 0, ILAENV( 3, 'DORGQL', ' ', M, N, K, -1 ) )
      IF( NX < K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
         LDWORK = N
         IWS = LDWORK*NB
         IF( LWORK < IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORGQL', ' ', M, N, K, -1 ) )
         END IF
      END IF
   END IF
!
   IF( NB >= NBMIN .AND. NB < K .AND. NX < K ) THEN
!
!        Use blocked code after the first block.
!        The last kk columns are handled by the block method.
!
      KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )
!
!        Set A(m-kk+1:m,1:n-kk) to zero.
!
      DO J = 1, N - KK
         DO I = M - KK + 1, M
            A( I, J ) = ZERO
         ENDDO
      ENDDO
   ELSE
      KK = 0
   END IF
!
!     Use unblocked code for the first or only block.
!
   CALL DORG2L( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
!
   IF( KK > 0 ) THEN
!
!        Use blocked code
!
      DO I = K - KK + 1, K, NB
         IB = MIN( NB, K-I+1 )
         IF( N-K+I > 1 ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i+ib-1) . . . H(i+1) H(i)
!
            CALL DLARFT( 'Backward', 'Columnwise', M-K+I+IB-1, IB, &
                         A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
!
            CALL DLARFB( 'Left', 'No transpose', 'Backward', &
                         'Columnwise', M-K+I+IB-1, N-K+I-1, IB, &
                         A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA, &
                         WORK( IB+1 ), LDWORK )
         END IF
!
!           Apply H to rows 1:m-k+i+ib-1 of current block
!
         CALL DORG2L( M-K+I+IB-1, IB, IB, A( 1, N-K+I ), LDA, &
                      TAU( I ), WORK, IINFO )
!
!           Set rows m-k+i+ib:m of current block to zero
!
         DO J = N - K + I, N - K + I + IB - 1
            DO L = M - K + I + IB, M
               A( L, J ) = ZERO
            ENDDO
         ENDDO
      ENDDO
   END IF
!
   WORK( 1 ) = IWS
   RETURN
!
!     End of DORGQL
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

