!> \brief \b SSXT1
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION SSXT1( IJOB, D1, N1, D2, N2, ABSTOL,
!                        ULP, UNFL )
!
!       .. Scalar Arguments ..
!       INTEGER            IJOB, N1, N2
!       REAL               ABSTOL, ULP, UNFL
!       ..
!       .. Array Arguments ..
!       REAL               D1( * ), D2( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSXT1  computes the difference between a set of eigenvalues.
!> The result is returned as the function value.
!>
!> IJOB = 1:   Computes   max { min | D1(i)-D2(j) | }
!>                         i     j
!>
!> IJOB = 2:   Computes   max { min | D1(i)-D2(j) | /
!>                         i     j
!>                              ( ABSTOL + |D1(i)|*ULP ) }
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] IJOB
!> \verbatim
!>          IJOB is INTEGER
!>          Specifies the type of tests to be performed.  (See above.)
!> \endverbatim
!>
!> \param[in] D1
!> \verbatim
!>          D1 is REAL array, dimension (N1)
!>          The first array.  D1 should be in increasing order, i.e.,
!>          D1(j) <= D1(j+1).
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!>          The length of D1.
!> \endverbatim
!>
!> \param[in] D2
!> \verbatim
!>          D2 is REAL array, dimension (N2)
!>          The second array.  D2 should be in increasing order, i.e.,
!>          D2(j) <= D2(j+1).
!> \endverbatim
!>
!> \param[in] N2
!> \verbatim
!>          N2 is INTEGER
!>          The length of D2.
!> \endverbatim
!>
!> \param[in] ABSTOL
!> \verbatim
!>          ABSTOL is REAL
!>          The absolute tolerance, used as a measure of the error.
!> \endverbatim
!>
!> \param[in] ULP
!> \verbatim
!>          ULP is REAL
!>          Machine precision.
!> \endverbatim
!>
!> \param[in] UNFL
!> \verbatim
!>          UNFL is REAL
!>          The smallest positive number whose reciprocal does not
!>          overflow.
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
   REAL             FUNCTION SSXT1( IJOB, D1, N1, D2, N2, ABSTOL, &
                    ULP, UNFL )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            IJOB, N1, N2
   REAL               ABSTOL, ULP, UNFL
!     ..
!     .. Array Arguments ..
   REAL               D1( * ), D2( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, J
   REAL               TEMP1, TEMP2
!     ..
!     .. Executable Statements ..
!
   TEMP1 = 0.0E+0
!
   J = 1
   DO I = 1, N1
      DO WHILE ( D2( J ) < D1( I ) .AND. J < N2 )
         J = J + 1
      ENDDO
      IF( J == 1 ) THEN
         TEMP2 = ABS( D2( J )-D1( I ) )
         IF( IJOB == 2 ) TEMP2 = TEMP2 / MAX( UNFL, ABSTOL+ULP*ABS( D1( I ) ) )
      ELSE
         TEMP2 = MIN( ABS( D2( J )-D1( I ) ), ABS( D1( I )-D2( J-1 ) ) )
         IF( IJOB == 2 ) TEMP2 = TEMP2 / MAX( UNFL, ABSTOL+ULP*ABS( D1( I ) ) )
      END IF
      TEMP1 = MAX( TEMP1, TEMP2 )
   ENDDO
!
   SSXT1 = TEMP1
   RETURN
!
!     End of SSXT1
!
END



