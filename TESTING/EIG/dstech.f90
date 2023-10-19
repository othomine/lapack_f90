!> \brief \b DSTECH
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSTECH( N, A, B, EIG, TOL, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, N
!       DOUBLE PRECISION   TOL
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( * ), B( * ), EIG( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    Let T be the tridiagonal matrix with diagonal entries A(1) ,...,
!>    A(N) and offdiagonal entries B(1) ,..., B(N-1)).  DSTECH checks to
!>    see if EIG(1) ,..., EIG(N) are indeed accurate eigenvalues of T.
!>    It does this by expanding each EIG(I) into an interval
!>    [SVD(I) - EPS, SVD(I) + EPS], merging overlapping intervals if
!>    any, and using Sturm sequences to count and verify whether each
!>    resulting interval has the correct number of eigenvalues (using
!>    DSTECT).  Here EPS = TOL*MAZHEPS*MAXEIG, where MACHEPS is the
!>    machine precision and MAXEIG is the absolute value of the largest
!>    eigenvalue. If each interval contains the correct number of
!>    eigenvalues, INFO = 0 is returned, otherwise INFO is the index of
!>    the first eigenvalue in the first bad interval.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The dimension of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (N)
!>          The diagonal entries of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (N-1)
!>          The offdiagonal entries of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] EIG
!> \verbatim
!>          EIG is DOUBLE PRECISION array, dimension (N)
!>          The purported eigenvalues to be checked.
!> \endverbatim
!>
!> \param[in] TOL
!> \verbatim
!>          TOL is DOUBLE PRECISION
!>          Error tolerance for checking, a multiple of the
!>          machine precision.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          0  if the eigenvalues are all correct (to within
!>             1 +- TOL*MAZHEPS*MAXEIG)
!>          >0 if the interval containing the INFO-th eigenvalue
!>             contains the incorrect number of eigenvalues.
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
!> \ingroup double_eig
!
!  =====================================================================
   SUBROUTINE DSTECH( N, A, B, EIG, TOL, WORK, INFO )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, N
   DOUBLE PRECISION   TOL
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   A( * ), B( * ), EIG( * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            BPNT, COUNT, I, ISUB, J, NUML, NUMU, TPNT
   DOUBLE PRECISION   EMIN, EPS, LOWER, MX, TUPPR, UNFLEP, UPPER
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH
   EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           DSTECT
!     ..
!     .. Executable Statements ..
!
!     Check input parameters
!
   INFO = 0
   IF( N == 0 ) RETURN
   IF( N < 0 ) THEN
      INFO = -1
      RETURN
   END IF
   IF( TOL < 0.0D+0 ) THEN
      INFO = -5
      RETURN
   END IF
!
!     Get machine constants
!
   EPS = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
   UNFLEP = DLAMCH( 'Safe minimum' ) / EPS
   EPS = TOL*EPS
!
!     Compute maximum absolute eigenvalue, error tolerance
!
   MX = MAXVAL(ABS(EIG(1:N)))
   EPS = MAX( EPS*MX, UNFLEP )
!
!     Sort eigenvalues from EIG into WORK
!
   WORK(1:N) = EIG(1:N)
   DO I = 1, N - 1
      ISUB = 1
      EMIN = WORK( 1 )
      DO J = 2, N + 1 - I
         IF( WORK( J ) < EMIN ) THEN
            ISUB = J
            EMIN = WORK( J )
         END IF
      ENDDO
      IF( ISUB /= N+1-I ) THEN
         WORK( ISUB ) = WORK( N+1-I )
         WORK( N+1-I ) = EMIN
      END IF
   ENDDO
!
!     TPNT points to singular value at right endpoint of interval
!     BPNT points to singular value at left  endpoint of interval
!
   TPNT = 1
   BPNT = 1
!
!     Begin loop over all intervals
!
   DO
   UPPER = WORK( TPNT ) + EPS
   LOWER = WORK( BPNT ) - EPS
!
!     Begin loop merging overlapping intervals
!
   DO
     IF( BPNT == N ) EXIT
     TUPPR = WORK( BPNT+1 ) + EPS
     IF( TUPPR < LOWER ) EXIT
!
!     Merge
!
     BPNT = BPNT + 1
     LOWER = WORK( BPNT ) - EPS
   ENDDO
!
!     Count singular values in interval [ LOWER, UPPER ]
!
   CALL DSTECT( N, A, B, LOWER, NUML )
   CALL DSTECT( N, A, B, UPPER, NUMU )
   COUNT = NUMU - NUML
   IF( COUNT /= BPNT-TPNT+1 ) THEN
!
!        Wrong number of singular values in interval
!
      INFO = TPNT
      RETURN
   END IF
   TPNT = BPNT + 1
   BPNT = TPNT
   IF( TPNT > N ) RETURN
   ENDDO
!
!     End of DSTECH
!
END



