!> \brief \b SSVDCH
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSVDCH( N, S, E, SVD, TOL, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, N
!       REAL               TOL
!       ..
!       .. Array Arguments ..
!       REAL               E( * ), S( * ), SVD( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSVDCH checks to see if SVD(1) ,..., SVD(N) are accurate singular
!> values of the bidiagonal matrix B with diagonal entries
!> S(1) ,..., S(N) and superdiagonal entries E(1) ,..., E(N-1)).
!> It does this by expanding each SVD(I) into an interval
!> [SVD(I) * (1-EPS) , SVD(I) * (1+EPS)], merging overlapping intervals
!> if any, and using Sturm sequences to count and verify whether each
!> resulting interval has the correct number of singular values (using
!> SSVDCT). Here EPS=TOL*MAX(N/10,1)*MACHEP, where MACHEP is the
!> machine precision. The routine assumes the singular values are sorted
!> with SVD(1) the largest and SVD(N) smallest.  If each interval
!> contains the correct number of singular values, INFO = 0 is returned,
!> otherwise INFO is the index of the first singular value in the first
!> bad interval.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The dimension of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is REAL array, dimension (N)
!>          The diagonal entries of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          The superdiagonal entries of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] SVD
!> \verbatim
!>          SVD is REAL array, dimension (N)
!>          The computed singular values to be checked.
!> \endverbatim
!>
!> \param[in] TOL
!> \verbatim
!>          TOL is REAL
!>          Error tolerance for checking, a multiplier of the
!>          machine precision.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          =0 if the singular values are all correct (to within
!>             1 +- TOL*MACHEPS)
!>          >0 if the interval containing the INFO-th singular value
!>             contains the incorrect number of singular values.
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
!> \ingroup single_eig
!
!  =====================================================================
   SUBROUTINE SSVDCH( N, S, E, SVD, TOL, INFO )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, N
   REAL               TOL
!     ..
!     .. Array Arguments ..
   REAL               E( * ), S( * ), SVD( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            BPNT, COUNT, NUML, NUMU, TPNT
   REAL               EPS, LOWER, OVFL, TUPPR, UNFL, UNFLEP, UPPER
!     ..
!     .. External Functions ..
   REAL               SLAMCH
   EXTERNAL           SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           SSVDCT
!     ..
!     .. Executable Statements ..
!
!     Get machine constants
!
   INFO = 0
   IF( N <= 0 ) RETURN
   UNFL = SLAMCH( 'Safe minimum' )
   OVFL = SLAMCH( 'Overflow' )
   EPS = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
!
!     UNFLEP is chosen so that when an eigenvalue is multiplied by the
!     scale factor sqrt(OVFL)*sqrt(sqrt(UNFL))/MX in SSVDCT, it exceeds
!     sqrt(UNFL), which is the lower limit for SSVDCT.
!
   UNFLEP = ( SQRT( SQRT( UNFL ) ) / SQRT( OVFL ) )*SVD( 1 ) + UNFL / EPS
!
!     The value of EPS works best when TOL  >=  10.
!
   EPS = TOL*MAX( N / 10, 1 )*EPS
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
   UPPER = ( 1.0E+0+EPS )*SVD( TPNT ) + UNFLEP
   LOWER = ( 1.0E+0-EPS )*SVD( BPNT ) - UNFLEP
   IF( LOWER <= UNFLEP ) LOWER = -UPPER
!
!     Begin loop merging overlapping intervals
!
   DO
     IF( BPNT == N ) EXIT
     TUPPR = ( 1.0E+0+EPS )*SVD( BPNT+1 ) + UNFLEP
     IF( TUPPR < LOWER ) EXIT
!
!     Merge
!
     BPNT = BPNT + 1
     LOWER = ( 1.0E+0-EPS )*SVD( BPNT ) - UNFLEP
     IF( LOWER <= UNFLEP ) LOWER = -UPPER
   ENDDO

!
!     Count singular values in interval [ LOWER, UPPER ]
!
   CALL SSVDCT( N, S, E, LOWER, NUML )
   CALL SSVDCT( N, S, E, UPPER, NUMU )
   COUNT = NUMU - NUML
   IF( LOWER < 0.0E+0 ) COUNT = COUNT / 2
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
!     End of SSVDCH
!
END



