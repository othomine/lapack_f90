!> \brief \b CSRSCL multiplies a vector by the reciprocal of a real scalar.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CSRSCL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csrscl.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csrscl.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csrscl.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSRSCL( N, SA, SX, INCX )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       REAL               SA
!       ..
!       .. Array Arguments ..
!       COMPLEX            SX( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSRSCL multiplies an n-element complex vector x by the real scalar
!> 1/a.  This is done without overflow or underflow as long as
!> the final result x/a does not overflow or underflow.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of components of the vector x.
!> \endverbatim
!>
!> \param[in] SA
!> \verbatim
!>          SA is REAL
!>          The scalar a which is used to divide each component of x.
!>          SA must be >= 0, or the subroutine will divide by zero.
!> \endverbatim
!>
!> \param[in,out] SX
!> \verbatim
!>          SX is COMPLEX array, dimension
!>                         (1+(N-1)*abs(INCX))
!>          The n-element vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between successive values of the vector SX.
!>          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n
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
!> \ingroup rscl
!
!  =====================================================================
   SUBROUTINE CSRSCL( N, SA, SX, INCX )
   IMPLICIT NONE
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INCX, N
   REAL               SA
!     ..
!     .. Array Arguments ..
   COMPLEX            SX( * )
!     ..
!
! =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            DONE
   REAL               BIGNUM, CDEN, CDEN1, CNUM, CNUM1, MUL, SMLNUM
!     ..
!     .. External Functions ..
   REAL               SLAMCH
   EXTERNAL           SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           CSSCAL
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
   IF( N <= 0 ) RETURN
!
!     Get machine parameters
!
   SMLNUM = SLAMCH( 'S' )
   BIGNUM = 1.0E+0 / SMLNUM
!
!     Initialize the denominator to SA and the numerator to 1.
!
   CDEN = SA
   CNUM = 1.0E+0
!
10 CONTINUE
   CDEN1 = CDEN*SMLNUM
   CNUM1 = CNUM / BIGNUM
   IF( ABS( CDEN1 ) > ABS( CNUM ) .AND. CNUM /= 0.0E+0 ) THEN
!
!        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM.
!
      MUL = SMLNUM
      DONE = .FALSE.
      CDEN = CDEN1
   ELSE IF( ABS( CNUM1 ) > ABS( CDEN ) ) THEN
!
!        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM.
!
      MUL = BIGNUM
      DONE = .FALSE.
      CNUM = CNUM1
   ELSE
!
!        Multiply X by CNUM / CDEN and return.
!
      MUL = CNUM / CDEN
      DONE = .TRUE.
   END IF
!
!     Scale the vector X by MUL
!
   SX(1:N*INCX:INCX) = CMPLX(MUL*REAL(SX(1:N*INCX:INCX)),MUL*AIMAG(SX(1:N*INCX:INCX)))
!
   IF( .NOT.DONE ) GO TO 10
!
   RETURN
!
!     End of CSRSCL
!
END
