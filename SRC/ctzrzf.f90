!> \brief \b CTZRZF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CTZRZF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctzrzf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctzrzf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctzrzf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTZRZF reduces the M-by-N ( M<=N ) complex upper trapezoidal matrix A
!> to upper triangular form by means of unitary transformations.
!>
!> The upper trapezoidal matrix A is factored as
!>
!>    A = ( R  0 ) * Z,
!>
!> where Z is an N-by-N unitary matrix and R is an M-by-M upper
!> triangular matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= M.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the leading M-by-N upper trapezoidal part of the
!>          array A must contain the matrix to be factorized.
!>          On exit, the leading M-by-M upper triangular part of A
!>          contains the upper triangular matrix R, and elements M+1 to
!>          N of the first M rows of A, with the array TAU, represent the
!>          unitary matrix Z as a product of M elementary reflectors.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (M)
!>          The scalar factors of the elementary reflectors.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= max(1,M).
!>          For optimum performance LWORK >= M*NB, where NB is
!>          the optimal blocksize.
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
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup tzrzf
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
!>
!>  The N-by-N matrix Z can be computed by
!>
!>     Z =  Z(1)*Z(2)* ... *Z(M)
!>
!>  where each N-by-N Z(k) is given by
!>
!>     Z(k) = I - tau(k)*v(k)*v(k)**H
!>
!>  with v(k) is the kth row vector of the M-by-N matrix
!>
!>     V = ( I   A(:,M+1:N) )
!>
!>  I is the M-by-M identity matrix, A(:,M+1:N)
!>  is the output stored in A on exit from CTZRZF,
!>  and tau(k) is the kth element of the array TAU.
!>
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE CTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
   COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   COMPLEX            ZERO
   PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
   LOGICAL            LQUERY
   INTEGER            I, IB, IWS, KI, KK, LDWORK, LWKMIN, LWKOPT, &
                      M1, MU, NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA, CLARZB, CLARZT, CLATRZ
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
   ELSE IF( N < M ) THEN
      INFO = -2
   ELSE IF( LDA < MAX( 1, M ) ) THEN
      INFO = -4
   END IF
!
   IF( INFO == 0 ) THEN
      IF( M == 0 .OR. M == N ) THEN
         LWKOPT = 1
         LWKMIN = 1
      ELSE
!
!           Determine the block size.
!
         NB = ILAENV( 1, 'CGERQF', ' ', M, N, -1, -1 )
         LWKOPT = M*NB
         LWKMIN = MAX( 1, M )
      END IF
      WORK( 1 ) = LWKOPT
!
      IF( LWORK < LWKMIN .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
   END IF
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CTZRZF', -INFO )
      RETURN
   ELSE IF( LQUERY ) THEN
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( M == 0 ) THEN
      RETURN
   ELSE IF( M == N ) THEN
      DO I = 1, N
         TAU( I ) = ZERO
      ENDDO
      RETURN
   END IF
!
   NBMIN = 2
   NX = 1
   IWS = M
   IF( NB > 1 .AND. NB < M ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
      NX = MAX( 0, ILAENV( 3, 'CGERQF', ' ', M, N, -1, -1 ) )
      IF( NX < M ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
         LDWORK = M
         IWS = LDWORK*NB
         IF( LWORK < IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'CGERQF', ' ', M, N, -1, &
                    -1 ) )
         END IF
      END IF
   END IF
!
   IF( NB >= NBMIN .AND. NB < M .AND. NX < M ) THEN
!
!        Use blocked code initially.
!        The last kk rows are handled by the block method.
!
      M1 = MIN( M+1, N )
      KI = ( ( M-NX-1 ) / NB )*NB
      KK = MIN( M, KI+NB )
!
      DO I = M - KK + KI + 1, M - KK + 1, -NB
         IB = MIN( M-I+1, NB )
!
!           Compute the TZ factorization of the current block
!           A(i:i+ib-1,i:n)
!
         CALL CLATRZ( IB, N-I+1, N-M, A( I, I ), LDA, TAU( I ), &
                      WORK )
         IF( I > 1 ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i+ib-1) . . . H(i+1) H(i)
!
            CALL CLARZT( 'Backward', 'Rowwise', N-M, IB, A( I, M1 ), &
                         LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(1:i-1,i:n) from the right
!
            CALL CLARZB( 'Right', 'No transpose', 'Backward', &
                         'Rowwise', I-1, N-I+1, IB, N-M, A( I, M1 ), &
                         LDA, WORK, LDWORK, A( 1, I ), LDA, &
                         WORK( IB+1 ), LDWORK )
         END IF
      ENDDO
      MU = I + NB - 1
   ELSE
      MU = M
   END IF
!
!     Use unblocked code to factor the last or only block
!
   IF( MU > 0 ) &
      CALL CLATRZ( MU, N, N-M, A, LDA, TAU, WORK )
!
   WORK( 1 ) = LWKOPT
!
   RETURN
!
!     End of CTZRZF
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        