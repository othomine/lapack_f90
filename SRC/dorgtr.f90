!> \brief \b DORGTR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DORGTR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgtr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgtr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgtr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LWORK, N
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
!> DORGTR generates a real orthogonal matrix Q which is defined as the
!> product of n-1 elementary reflectors of order N, as returned by
!> DSYTRD:
!>
!> if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
!>
!> if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U': Upper triangle of A contains elementary reflectors
!>                 from DSYTRD;
!>          = 'L': Lower triangle of A contains elementary reflectors
!>                 from DSYTRD.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix Q. N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the vectors which define the elementary reflectors,
!>          as returned by DSYTRD.
!>          On exit, the N-by-N orthogonal matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (N-1)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DSYTRD.
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
!>          The dimension of the array WORK. LWORK >= max(1,N-1).
!>          For optimum performance LWORK >= (N-1)*NB, where NB is
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
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
!
!> \ingroup ungtr
!
!  =====================================================================
   SUBROUTINE DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            LQUERY, UPPER
   INTEGER            I, IINFO, J, LWKOPT, NB
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ILAENV
   EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
   EXTERNAL           DORGQL, DORGQR, XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
   INFO = 0
   LQUERY = ( LWORK == -1 )
   UPPER = LSAME( UPLO, 'U' )
   IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -4
   ELSE IF( LWORK < MAX( 1, N-1 ) .AND. .NOT.LQUERY ) THEN
      INFO = -7
   END IF
!
   IF( INFO == 0 ) THEN
      IF( UPPER ) THEN
         NB = ILAENV( 1, 'DORGQL', ' ', N-1, N-1, N-1, -1 )
      ELSE
         NB = ILAENV( 1, 'DORGQR', ' ', N-1, N-1, N-1, -1 )
      END IF
      LWKOPT = MAX( 1, N-1 )*NB
      WORK( 1 ) = LWKOPT
   END IF
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'DORGTR', -INFO )
      RETURN
   ELSE IF( LQUERY ) THEN
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) THEN
      WORK( 1 ) = 1
      RETURN
   END IF
!
   IF( UPPER ) THEN
!
!        Q was determined by a call to DSYTRD with UPLO = 'U'
!
!        Shift the vectors which define the elementary reflectors one
!        column to the left, and set the last row and column of Q to
!        those of the unit matrix
!
      DO J = 1, N - 1
         DO I = 1, J - 1
            A( I, J ) = A( I, J+1 )
         ENDDO
         A( N, J ) = ZERO
      ENDDO
      DO I = 1, N - 1
         A( I, N ) = ZERO
      ENDDO
      A( N, N ) = ONE
!
!        Generate Q(1:n-1,1:n-1)
!
      CALL DORGQL( N-1, N-1, N-1, A, LDA, TAU, WORK, LWORK, IINFO )
!
   ELSE
!
!        Q was determined by a call to DSYTRD with UPLO = 'L'.
!
!        Shift the vectors which define the elementary reflectors one
!        column to the right, and set the first row and column of Q to
!        those of the unit matrix
!
      DO J = N, 2, -1
         A( 1, J ) = ZERO
         DO I = J + 1, N
            A( I, J ) = A( I, J-1 )
         ENDDO
      ENDDO
      A( 1, 1 ) = ONE
      DO I = 2, N
         A( I, 1 ) = ZERO
      ENDDO
      IF( N > 1 ) THEN
!
!           Generate Q(2:n,2:n)
!
         CALL DORGQR( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK, &
                      LWORK, IINFO )
      END IF
   END IF
   WORK( 1 ) = LWKOPT
   RETURN
!
!     End of DORGTR
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

