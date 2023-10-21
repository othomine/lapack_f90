!> \brief \b CLAPMR rearranges rows of a matrix as specified by a permutation vector.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAPMR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clapmr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clapmr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clapmr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAPMR( FORWRD, M, N, X, LDX, K )
!
!       .. Scalar Arguments ..
!       LOGICAL            FORWRD
!       INTEGER            LDX, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            K( * )
!       COMPLEX            X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAPMR rearranges the rows of the M by N matrix X as specified
!> by the permutation K(1),K(2),...,K(M) of the integers 1,...,M.
!> If FORWRD = .TRUE.,  forward permutation:
!>
!>      X(K(I),*) is moved X(I,*) for I = 1,2,...,M.
!>
!> If FORWRD = .FALSE., backward permutation:
!>
!>      X(I,*) is moved to X(K(I),*) for I = 1,2,...,M.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] FORWRD
!> \verbatim
!>          FORWRD is LOGICAL
!>          = .TRUE., forward permutation
!>          = .FALSE., backward permutation
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix X. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix X. N >= 0.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX array, dimension (LDX,N)
!>          On entry, the M by N matrix X.
!>          On exit, X contains the permuted matrix X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X, LDX >= MAX(1,M).
!> \endverbatim
!>
!> \param[in,out] K
!> \verbatim
!>          K is INTEGER array, dimension (M)
!>          On entry, K contains the permutation vector. K is used as
!>          internal workspace, but reset to its original value on
!>          output.
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
!> \ingroup lapmr
!
!  =====================================================================
   SUBROUTINE CLAPMR( FORWRD, M, N, X, LDX, K )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            FORWRD
   INTEGER            LDX, M, N
!     ..
!     .. Array Arguments ..
   INTEGER            K( * )
   COMPLEX            X( LDX, * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
   INTEGER            I, IN, J, JJ
   COMPLEX            TEMP
!     ..
!     .. Executable Statements ..
!
   IF( M <= 1 ) &
      RETURN
!
   DO I = 1, M
      K( I ) = -K( I )
   ENDDO
!
   IF( FORWRD ) THEN
!
!        Forward permutation
!
      DO I = 1, M
!
         IF( K( I ) > 0 ) &
            GO TO 40
!
         J = I
         K( J ) = -K( J )
         IN = K( J )
!
20       CONTINUE
         IF( K( IN ) > 0 ) &
            GO TO 40
!
         DO JJ = 1, N
            TEMP = X( J, JJ )
            X( J, JJ ) = X( IN, JJ )
            X( IN, JJ ) = TEMP
         ENDDO
!
         K( IN ) = -K( IN )
         J = IN
         IN = K( IN )
         GO TO 20
!
40       CONTINUE
!
      ENDDO
!
   ELSE
!
!        Backward permutation
!
      DO I = 1, M
!
         IF( K( I ) > 0 ) &
            GO TO 80
!
         K( I ) = -K( I )
         J = K( I )
60       CONTINUE
         IF( J == I ) &
            GO TO 80
!
         DO JJ = 1, N
            TEMP = X( I, JJ )
            X( I, JJ ) = X( J, JJ )
            X( J, JJ ) = TEMP
         ENDDO
!
         K( J ) = -K( J )
         J = K( J )
         GO TO 60
!
80       CONTINUE
!
      ENDDO
!
   END IF
!
   RETURN
!
!     End of CLAPMR
!
   END

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

