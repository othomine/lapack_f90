!> \brief \b DORGTSQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DORGTSQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgtsqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgtsqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgtsqr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORGTSQR( M, N, MB, NB, A, LDA, T, LDT, WORK, LWORK,
!      $                     INFO )
!
!       .. Scalar Arguments ..
!       INTEGER           INFO, LDA, LDT, LWORK, M, N, MB, NB
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION  A( LDA, * ), T( LDT, * ), WORK( * )
!       ..
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DORGTSQR generates an M-by-N real matrix Q_out with orthonormal columns,
!> which are the first N columns of a product of real orthogonal
!> matrices of order M which are returned by DLATSQR
!>
!>      Q_out = first_N_columns_of( Q(1)_in * Q(2)_in * ... * Q(k)_in ).
!>
!> See the documentation for DLATSQR.
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
!>          The number of columns of the matrix A. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] MB
!> \verbatim
!>          MB is INTEGER
!>          The row block size used by DLATSQR to return
!>          arrays A and T. MB > N.
!>          (Note that if MB > M, then M is used instead of MB
!>          as the row block size).
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The column block size used by DLATSQR to return
!>          arrays A and T. NB >= 1.
!>          (Note that if NB > N, then N is used instead of NB
!>          as the column block size).
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>
!>          On entry:
!>
!>             The elements on and above the diagonal are not accessed.
!>             The elements below the diagonal represent the unit
!>             lower-trapezoidal blocked matrix V computed by DLATSQR
!>             that defines the input matrices Q_in(k) (ones on the
!>             diagonal are not stored) (same format as the output A
!>             below the diagonal in DLATSQR).
!>
!>          On exit:
!>
!>             The array A contains an M-by-N orthonormal matrix Q_out,
!>             i.e the columns of A are orthogonal unit vectors.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is DOUBLE PRECISION array,
!>          dimension (LDT, N * NIRB)
!>          where NIRB = Number_of_input_row_blocks
!>                     = MAX( 1, CEIL((M-N)/(MB-N)) )
!>          Let NICB = Number_of_input_col_blocks
!>                   = CEIL(N/NB)
!>
!>          The upper-triangular block reflectors used to define the
!>          input matrices Q_in(k), k=(1:NIRB*NICB). The block
!>          reflectors are stored in compact form in NIRB block
!>          reflector sequences. Each of NIRB block reflector sequences
!>          is stored in a larger NB-by-N column block of T and consists
!>          of NICB smaller NB-by-NB upper-triangular column blocks.
!>          (same format as the output T in DLATSQR).
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.
!>          LDT >= max(1,min(NB1,N)).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          (workspace) DOUBLE PRECISION array, dimension (MAX(2,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= (M+NB)*N.
!>          If LWORK = -1, then a workspace query is assumed.
!>          The routine only calculates the optimal size of the WORK
!>          array, returns this value as the first entry of the WORK
!>          array, and no error message related to LWORK is issued
!>          by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!>
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup ungtsqr
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!> November 2019, Igor Kozachenko,
!>                Computer Science Division,
!>                University of California, Berkeley
!>
!> \endverbatim
!
!  =====================================================================
   SUBROUTINE DORGTSQR( M, N, MB, NB, A, LDA, T, LDT, WORK, LWORK, &
                        INFO )
   IMPLICIT NONE
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER           INFO, LDA, LDT, LWORK, M, N, MB, NB
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION  A( LDA, * ), T( LDT, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ONE, ZERO
   PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            LQUERY
   INTEGER            IINFO, LDC, LWORKOPT, LC, LW, NBLOCAL, J
!     ..
!     .. External Subroutines ..
   EXTERNAL           DCOPY, DLAMTSQR, DLASET, XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          DBLE, MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
   LQUERY  = LWORK == -1
   INFO = 0
   IF( M < 0 ) THEN
      INFO = -1
   ELSE IF( N < 0 .OR. M < N ) THEN
      INFO = -2
   ELSE IF( MB <= N ) THEN
      INFO = -3
   ELSE IF( NB < 1 ) THEN
      INFO = -4
   ELSE IF( LDA < MAX( 1, M ) ) THEN
      INFO = -6
   ELSE IF( LDT < MAX( 1, MIN( NB, N ) ) ) THEN
      INFO = -8
   ELSE
!
!        Test the input LWORK for the dimension of the array WORK.
!        This workspace is used to store array C(LDC, N) and WORK(LWORK)
!        in the call to DLAMTSQR. See the documentation for DLAMTSQR.
!
      IF( LWORK < 2 .AND. (.NOT.LQUERY) ) THEN
         INFO = -10
      ELSE
!
!           Set block size for column blocks
!
         NBLOCAL = MIN( NB, N )
!
!           LWORK = -1, then set the size for the array C(LDC,N)
!           in DLAMTSQR call and set the optimal size of the work array
!           WORK(LWORK) in DLAMTSQR call.
!
         LDC = M
         LC = LDC*N
         LW = N * NBLOCAL
!
         LWORKOPT = LC+LW
!
         IF( ( LWORK < MAX( 1, LWORKOPT ) ).AND.(.NOT.LQUERY) ) THEN
            INFO = -10
         END IF
      END IF
!
   END IF
!
!     Handle error in the input parameters and return workspace query.
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'DORGTSQR', -INFO )
      RETURN
   ELSE IF ( LQUERY ) THEN
      WORK( 1 ) = DBLE( LWORKOPT )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( MIN( M, N ) == 0 ) THEN
      WORK( 1 ) = DBLE( LWORKOPT )
      RETURN
   END IF
!
!     (1) Form explicitly the tall-skinny M-by-N left submatrix Q1_in
!     of M-by-M orthogonal matrix Q_in, which is implicitly stored in
!     the subdiagonal part of input array A and in the input array T.
!     Perform by the following operation using the routine DLAMTSQR.
!
!         Q1_in = Q_in * ( I ), where I is a N-by-N identity matrix,
!                        ( 0 )        0 is a (M-N)-by-N zero matrix.
!
!     (1a) Form M-by-N matrix in the array WORK(1:LDC*N) with ones
!     on the diagonal and zeros elsewhere.
!
   CALL DLASET( 'F', M, N, ZERO, ONE, WORK, LDC )
!
!     (1b)  On input, WORK(1:LDC*N) stores ( I );
!                                          ( 0 )
!
!           On output, WORK(1:LDC*N) stores Q1_in.
!
   CALL DLAMTSQR( 'L', 'N', M, N, N, MB, NBLOCAL, A, LDA, T, LDT, &
                  WORK, LDC, WORK( LC+1 ), LW, IINFO )
!
!     (2) Copy the result from the part of the work array (1:M,1:N)
!     with the leading dimension LDC that starts at WORK(1) into
!     the output array A(1:M,1:N) column-by-column.
!
   DO J = 1, N
      CALL DCOPY( M, WORK( (J-1)*LDC + 1 ), 1, A( 1, J ), 1 )
   END DO
!
   WORK( 1 ) = DBLE( LWORKOPT )
   RETURN
!
!     End of DORGTSQR
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        