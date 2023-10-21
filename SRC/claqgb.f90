!> \brief \b CLAQGB scales a general band matrix, using row and column scaling factors computed by sgbequ.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAQGB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqgb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqgb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqgb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAQGB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND,
!                          AMAX, EQUED )
!
!       .. Scalar Arguments ..
!       CHARACTER          EQUED
!       INTEGER            KL, KU, LDAB, M, N
!       REAL               AMAX, COLCND, ROWCND
!       ..
!       .. Array Arguments ..
!       REAL               C( * ), R( * )
!       COMPLEX            AB( LDAB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAQGB equilibrates a general M by N band matrix A with KL
!> subdiagonals and KU superdiagonals using the row and scaling factors
!> in the vectors R and C.
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
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is COMPLEX array, dimension (LDAB,N)
!>          On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
!>          The j-th column of A is stored in the j-th column of the
!>          array AB as follows:
!>          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!>
!>          On exit, the equilibrated matrix, in the same storage format
!>          as A.  See EQUED for the form of the equilibrated matrix.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDA >= KL+KU+1.
!> \endverbatim
!>
!> \param[in] R
!> \verbatim
!>          R is REAL array, dimension (M)
!>          The row scale factors for A.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is REAL array, dimension (N)
!>          The column scale factors for A.
!> \endverbatim
!>
!> \param[in] ROWCND
!> \verbatim
!>          ROWCND is REAL
!>          Ratio of the smallest R(i) to the largest R(i).
!> \endverbatim
!>
!> \param[in] COLCND
!> \verbatim
!>          COLCND is REAL
!>          Ratio of the smallest C(i) to the largest C(i).
!> \endverbatim
!>
!> \param[in] AMAX
!> \verbatim
!>          AMAX is REAL
!>          Absolute value of largest matrix entry.
!> \endverbatim
!>
!> \param[out] EQUED
!> \verbatim
!>          EQUED is CHARACTER*1
!>          Specifies the form of equilibration that was done.
!>          = 'N':  No equilibration
!>          = 'R':  Row equilibration, i.e., A has been premultiplied by
!>                  diag(R).
!>          = 'C':  Column equilibration, i.e., A has been postmultiplied
!>                  by diag(C).
!>          = 'B':  Both row and column equilibration, i.e., A has been
!>                  replaced by diag(R) * A * diag(C).
!> \endverbatim
!
!> \par Internal Parameters:
!  =========================
!>
!> \verbatim
!>  THRESH is a threshold value used to decide if row or column scaling
!>  should be done based on the ratio of the row or column scaling
!>  factors.  If ROWCND < THRESH, row scaling is done, and if
!>  COLCND < THRESH, column scaling is done.
!>
!>  LARGE and SMALL are threshold values used to decide if row scaling
!>  should be done based on the absolute size of the largest matrix
!>  element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.
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
!> \ingroup laqgb
!
!  =====================================================================
   SUBROUTINE CLAQGB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, &
                      AMAX, EQUED )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          EQUED
   INTEGER            KL, KU, LDAB, M, N
   REAL               AMAX, COLCND, ROWCND
!     ..
!     .. Array Arguments ..
   REAL               C( * ), R( * )
   COMPLEX            AB( LDAB, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ONE, THRESH
   PARAMETER          ( ONE = 1.0E+0, THRESH = 0.1E+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, J
   REAL               CJ, LARGE, SMALL
!     ..
!     .. External Functions ..
   REAL               SLAMCH
   EXTERNAL           SLAMCH
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
   IF( M <= 0 .OR. N <= 0 ) THEN
      EQUED = 'N'
      RETURN
   END IF
!
!     Initialize LARGE and SMALL.
!
   SMALL = SLAMCH( 'Safe minimum' ) / SLAMCH( 'Precision' )
   LARGE = ONE / SMALL
!
   IF( ROWCND >= THRESH .AND. AMAX >= SMALL .AND. AMAX <= LARGE ) &
        THEN
!
!        No row scaling
!
      IF( COLCND >= THRESH ) THEN
!
!           No column scaling
!
         EQUED = 'N'
      ELSE
!
!           Column scaling
!
         DO J = 1, N
            CJ = C( J )
            DO I = MAX( 1, J-KU ), MIN( M, J+KL )
               AB( KU+1+I-J, J ) = CJ*AB( KU+1+I-J, J )
            ENDDO
         ENDDO
         EQUED = 'C'
      END IF
   ELSE IF( COLCND >= THRESH ) THEN
!
!        Row scaling, no column scaling
!
      DO J = 1, N
         DO I = MAX( 1, J-KU ), MIN( M, J+KL )
            AB( KU+1+I-J, J ) = R( I )*AB( KU+1+I-J, J )
         ENDDO
      ENDDO
      EQUED = 'R'
   ELSE
!
!        Row and column scaling
!
      DO J = 1, N
         CJ = C( J )
         DO I = MAX( 1, J-KU ), MIN( M, J+KL )
            AB( KU+1+I-J, J ) = CJ*R( I )*AB( KU+1+I-J, J )
         ENDDO
      ENDDO
      EQUED = 'B'
   END IF
!
   RETURN
!
!     End of CLAQGB
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

