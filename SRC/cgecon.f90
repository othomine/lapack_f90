!> \brief \b CGECON
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGECON + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgecon.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgecon.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgecon.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, RWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            INFO, LDA, N
!       REAL               ANORM, RCOND
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGECON estimates the reciprocal of the condition number of a general
!> complex matrix A, in either the 1-norm or the infinity-norm, using
!> the LU factorization computed by CGETRF.
!>
!> An estimate is obtained for norm(inv(A)), and the reciprocal of the
!> condition number is computed as
!>    RCOND = 1 / ( norm(A) * norm(inv(A)) ).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NORM
!> \verbatim
!>          NORM is CHARACTER*1
!>          Specifies whether the 1-norm condition number or the
!>          infinity-norm condition number is required:
!>          = '1' or 'O':  1-norm;
!>          = 'I':         Infinity-norm.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The factors L and U from the factorization A = P*L*U
!>          as computed by CGETRF.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] ANORM
!> \verbatim
!>          ANORM is REAL
!>          If NORM = '1' or 'O', the 1-norm of the original matrix A.
!>          If NORM = 'I', the infinity-norm of the original matrix A.
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is REAL
!>          The reciprocal of the condition number of the matrix A,
!>          computed as RCOND = 1/(norm(A) * norm(inv(A))).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          =-5:  if ANORM is NAN or negative.
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
!> \ingroup gecon
!
!  =====================================================================
   SUBROUTINE CGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, RWORK, &
                      INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          NORM
   INTEGER            INFO, LDA, N
   REAL               ANORM, RCOND
!     ..
!     .. Array Arguments ..
   REAL               RWORK( * )
   COMPLEX            A( LDA, * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            ONENRM
   CHARACTER          NORMIN
   INTEGER            IX, KASE, KASE1
   REAL               AINVNM, SCALE, SL, SMLNUM, SU
   COMPLEX            ZDUM
!     ..
!     .. Local Arrays ..
   INTEGER            ISAVE( 3 )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME, SISNAN
   INTEGER            ICAMAX
   REAL               SLAMCH
   EXTERNAL           LSAME, ICAMAX, SLAMCH, SISNAN
!     ..
!     .. External Subroutines ..
   EXTERNAL           CLACN2, CLATRS, CSRSCL, XERBLA
!     ..
!     .. Statement Functions ..
   REAL               CABS1
!     ..
!     .. Statement Function definitions ..
   CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   ONENRM = NORM == '1' .OR. LSAME( NORM, 'O' )
   IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -4
   ELSE IF( ANORM < 0.0E+0 .OR. SISNAN( ANORM ) ) THEN
      INFO = -5
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CGECON', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   RCOND = 0.0E+0
   IF( N == 0 ) THEN
      RCOND = 1.0E+0
      RETURN
   ELSE IF( ANORM == 0.0E+0 ) THEN
      RETURN
   END IF
!
   SMLNUM = SLAMCH( 'Safe minimum' )
!
!     Estimate the norm of inv(A).
!
   AINVNM = 0.0E+0
   NORMIN = 'N'
   IF( ONENRM ) THEN
      KASE1 = 1
   ELSE
      KASE1 = 2
   END IF
   KASE = 0
10 CONTINUE
   CALL CLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
   IF( KASE /= 0 ) THEN
      IF( KASE == KASE1 ) THEN
!
!           Multiply by inv(L).
!
         CALL CLATRS( 'Lower', 'No transpose', 'Unit', NORMIN, N, A, &
                      LDA, WORK, SL, RWORK, INFO )
!
!           Multiply by inv(U).
!
         CALL CLATRS( 'Upper', 'No transpose', 'Non-unit', NORMIN, N, &
                      A, LDA, WORK, SU, RWORK( N+1 ), INFO )
      ELSE
!
!           Multiply by inv(U**H).
!
         CALL CLATRS( 'Upper', 'Conjugate transpose', 'Non-unit', &
                      NORMIN, N, A, LDA, WORK, SU, RWORK( N+1 ), &
                      INFO )
!
!           Multiply by inv(L**H).
!
         CALL CLATRS( 'Lower', 'Conjugate transpose', 'Unit', NORMIN, &
                      N, A, LDA, WORK, SL, RWORK, INFO )
      END IF
!
!        Divide X by 1/(SL*SU) if doing so will not cause overflow.
!
      SCALE = SL*SU
      NORMIN = 'Y'
      IF( SCALE /= 1.0E+0 ) THEN
         IX = ICAMAX( N, WORK, 1 )
         IF( SCALE < CABS1( WORK( IX ) )*SMLNUM .OR. SCALE == 0.0E+0 ) RETURN
         CALL CSRSCL( N, SCALE, WORK, 1 )
      END IF
      GO TO 10
   END IF
!
!     Compute the estimate of the reciprocal condition number.
!
   IF( AINVNM /= 0.0E+0 ) RCOND = ( 1.0E+0 / AINVNM ) / ANORM
!
   RETURN
!
!     End of CGECON
!
END

