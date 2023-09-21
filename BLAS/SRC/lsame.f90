!> \brief \b LSAME
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION LSAME(CA,CB)
!
!       .. Scalar Arguments ..
!       CHARACTER CA,CB
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> LSAME returns .TRUE. if CA is the same letter as CB regardless of
!> case.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] CA
!> \verbatim
!>          CA is CHARACTER*1
!> \endverbatim
!>
!> \param[in] CB
!> \verbatim
!>          CB is CHARACTER*1
!>          CA and CB specify the single characters to be compared.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!> \author Olivier Thomine
!
!>
!>     converted to F90 and optimized 2023, Olivier Thomine
!> \ingroup lsame
!
!  =====================================================================
   LOGICAL PURE FUNCTION LSAME(CA,CB)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER, INTENT(IN) :: CA,CB
!     ..
!
! =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER INTA,INTB,ZCODE
!     ..
!
!     Test if the characters are equal
!
   LSAME = CA  ==  CB
   IF (LSAME) RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
   ZCODE = ICHAR('Z')
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
   INTA = ICHAR(CA)
   INTB = ICHAR(CB)
!
   IF (ZCODE == 90 .OR. ZCODE == 122) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
       IF (INTA >= 97 .AND. INTA <= 122) INTA = INTA - 32
       IF (INTB >= 97 .AND. INTB <= 122) INTB = INTB - 32
!
   ELSE IF (ZCODE == 233 .OR. ZCODE == 169) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
       IF (INTA >= 129 .AND. INTA <= 137 .OR. &
           INTA >= 145 .AND. INTA <= 153 .OR. &
           INTA >= 162 .AND. INTA <= 169) INTA = INTA + 64
       IF (INTB >= 129 .AND. INTB <= 137 .OR. &
           INTB >= 145 .AND. INTB <= 153 .OR. &
           INTB >= 162 .AND. INTB <= 169) INTB = INTB + 64
!
   ELSE IF (ZCODE == 218 .OR. ZCODE == 250) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
       IF (INTA >= 225 .AND. INTA <= 250) INTA = INTA - 32
       IF (INTB >= 225 .AND. INTB <= 250) INTB = INTB - 32
   END IF
   LSAME = INTA  ==  INTB
!
!     RETURN
!
!     End of LSAME
!
END
