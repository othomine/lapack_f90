!> \brief \b CLAESY computes the eigenvalues and eigenvectors of a 2-by-2 complex symmetric matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAESY + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claesy.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claesy.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claesy.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAESY( A, B, C, RT1, RT2, EVSCAL, CS1, SN1 )
!
!       .. Scalar Arguments ..
!       COMPLEX            A, B, C, CS1, EVSCAL, RT1, RT2, SN1
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAESY computes the eigendecomposition of a 2-by-2 symmetric matrix
!>    ( ( A, B );( B, C ) )
!> provided the norm of the matrix of eigenvectors is larger than
!> some threshold value.
!>
!> RT1 is the eigenvalue of larger absolute value, and RT2 of
!> smaller absolute value.  If the eigenvectors are computed, then
!> on return ( CS1, SN1 ) is the unit eigenvector for RT1, hence
!>
!> [  CS1     SN1   ] . [ A  B ] . [ CS1    -SN1   ] = [ RT1  0  ]
!> [ -SN1     CS1   ]   [ B  C ]   [ SN1     CS1   ]   [  0  RT2 ]
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] A
!> \verbatim
!>          A is COMPLEX
!>          The ( 1, 1 ) element of input matrix.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX
!>          The ( 1, 2 ) element of input matrix.  The ( 2, 1 ) element
!>          is also given by B, since the 2-by-2 matrix is symmetric.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is COMPLEX
!>          The ( 2, 2 ) element of input matrix.
!> \endverbatim
!>
!> \param[out] RT1
!> \verbatim
!>          RT1 is COMPLEX
!>          The eigenvalue of larger modulus.
!> \endverbatim
!>
!> \param[out] RT2
!> \verbatim
!>          RT2 is COMPLEX
!>          The eigenvalue of smaller modulus.
!> \endverbatim
!>
!> \param[out] EVSCAL
!> \verbatim
!>          EVSCAL is COMPLEX
!>          The complex value by which the eigenvector matrix was scaled
!>          to make it orthonormal.  If EVSCAL is zero, the eigenvectors
!>          were not computed.  This means one of two things:  the 2-by-2
!>          matrix could not be diagonalized, or the norm of the matrix
!>          of eigenvectors before scaling was larger than the threshold
!>          value 0.1E+0 (set below).
!> \endverbatim
!>
!> \param[out] CS1
!> \verbatim
!>          CS1 is COMPLEX
!> \endverbatim
!>
!> \param[out] SN1
!> \verbatim
!>          SN1 is COMPLEX
!>          If EVSCAL  /=  0,  ( CS1, SN1 ) is the unit right eigenvector
!>          for RT1.
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
!> \ingroup laesy
!
!  =====================================================================
   SUBROUTINE CLAESY( A, B, C, RT1, RT2, EVSCAL, CS1, SN1 )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   COMPLEX            A, B, C, CS1, EVSCAL, RT1, RT2, SN1
!     ..
!
! =====================================================================
!     ..
!     .. Local Scalars ..
   REAL               BABS, EVNORM, TABS, Z
   COMPLEX            S, T, TMP
!     ..
!     .. Executable Statements ..
!
!
!     Special case:  The matrix is actually diagonal.
!     To avoid divide by zero later, we treat this case separately.
!
   IF( ABS( B ) == 0.0E+0 ) THEN
      RT1 = A
      RT2 = C
      IF( ABS( RT1 ) < ABS( RT2 ) ) THEN
         TMP = RT1
         RT1 = RT2
         RT2 = TMP
         CS1 = 0.0E+0
         SN1 = 1.0E+0
      ELSE
         CS1 = 1.0E+0
         SN1 = 0.0E+0
      END IF
   ELSE
!
!        Compute the eigenvalues and eigenvectors.
!        The characteristic equation is
!           lambda **2 - (A+C) lambda + (A*C - B*B)
!        and we solve it using the quadratic formula.
!
      S = ( A+C )*(0.5E+0,0.0E+0)
      T = ( A-C )*(0.5E+0,0.0E+0)
!
!        Take the square root carefully to avoid over/under flow.
!
      BABS = ABS( B )
      TABS = ABS( T )
      Z = MAX( BABS, TABS )
      IF( Z > 0.0E+0 ) T = Z*SQRT( ( T / Z )**2+( B / Z )**2 )
!
!        Compute the two eigenvalues.  RT1 and RT2 are exchanged
!        if necessary so that RT1 will have the greater magnitude.
!
      RT1 = S + T
      RT2 = S - T
      IF( ABS( RT1 ) < ABS( RT2 ) ) THEN
         TMP = RT1
         RT1 = RT2
         RT2 = TMP
      END IF
!
!        Choose CS1 = 1 and SN1 to satisfy the first equation, then
!        scale the components of this eigenvector so that the matrix
!        of eigenvectors X satisfies  X * X**T = I .  (No scaling is
!        done if the norm of the eigenvalue matrix is less than 0.1E+0.)
!
      SN1 = ( RT1-A ) / B
      TABS = ABS( SN1 )
      IF( TABS > 1.0E+0 ) THEN
         T = TABS*SQRT( ( 1.0E+0 / TABS )**2+( SN1 / TABS )**2 )
      ELSE
         T = SQRT( (1.0E+0,0.0E+0)+SN1*SN1 )
      END IF
      EVNORM = ABS( T )
      IF( EVNORM >= 0.1E+0 ) THEN
         EVSCAL = (1.0E+0,0.0E+0) / T
         CS1 = EVSCAL
         SN1 = SN1*EVSCAL
      ELSE
         EVSCAL = 0.0E+0
      END IF
   END IF
   RETURN
!
!     End of CLAESY
!
END
