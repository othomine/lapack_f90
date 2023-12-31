!> \brief \b CTBCON
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CTBCON + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctbcon.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctbcon.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctbcon.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK,
!                          RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, NORM, UPLO
!       INTEGER            INFO, KD, LDAB, N
!       REAL               RCOND
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            AB( LDAB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTBCON estimates the reciprocal of the condition number of a
!> triangular band matrix A, in either the 1-norm or the infinity-norm.
!>
!> The norm of A is computed and an estimate is obtained for
!> norm(inv(A)), then the reciprocal of the condition number is
!> computed as
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
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  A is upper triangular;
!>          = 'L':  A is lower triangular.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          = 'N':  A is non-unit triangular;
!>          = 'U':  A is unit triangular.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of superdiagonals or subdiagonals of the
!>          triangular band matrix A.  KD >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is COMPLEX array, dimension (LDAB,N)
!>          The upper or lower triangular band matrix A, stored in the
!>          first kd+1 rows of the array. The j-th column of A is stored
!>          in the j-th column of the array AB as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>          If DIAG = 'U', the diagonal elements of A are not referenced
!>          and are assumed to be 1.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
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
!>          RWORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup tbcon
!
!  =====================================================================
   SUBROUTINE CTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, &
                      RWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          DIAG, NORM, UPLO
   INTEGER            INFO, KD, LDAB, N
   REAL               RCOND
!     ..
!     .. Array Arguments ..
   REAL               RWORK( * )
   COMPLEX            AB( LDAB, * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            NOUNIT, ONENRM, UPPER
   CHARACTER          NORMIN
   INTEGER            IX, KASE, KASE1
   REAL               AINVNM, ANORM, SCALE, SMLNUM, XNORM
   COMPLEX            ZDUM
!     ..
!     .. Local Arrays ..
   INTEGER            ISAVE( 3 )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ICAMAX
   REAL               CLANTB, SLAMCH, CABS1
   EXTERNAL           LSAME, ICAMAX, CLANTB, SLAMCH, CABS1
!     ..
!     .. External Subroutines ..
   EXTERNAL           CLACN2, CLATBS, CSRSCL, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   UPPER = LSAME( UPLO, 'U' )
   ONENRM = NORM == '1' .OR. LSAME( NORM, 'O' )
   NOUNIT = LSAME( DIAG, 'N' )
!
   IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
      INFO = -1
   ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -2
   ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
      INFO = -3
   ELSE IF( N < 0 ) THEN
      INFO = -4
   ELSE IF( KD < 0 ) THEN
      INFO = -5
   ELSE IF( LDAB < KD+1 ) THEN
      INFO = -7
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CTBCON', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) THEN
      RCOND = 1.0E+0
      RETURN
   END IF
!
   RCOND = 0.0E+0
   SMLNUM = SLAMCH( 'Safe minimum' )*REAL( MAX( N, 1 ) )
!
!     Compute the 1-norm of the triangular matrix A or A**H.
!
   ANORM = CLANTB( NORM, UPLO, DIAG, N, KD, AB, LDAB, RWORK )
!
!     Continue only if ANORM > 0.
!
   IF( ANORM > 0.0E+0 ) THEN
!
!        Estimate the 1-norm of the inverse of A.
!
      AINVNM = 0.0E+0
      NORMIN = 'N'
      IF( ONENRM ) THEN
         KASE1 = 1
      ELSE
         KASE1 = 2
      END IF
      KASE = 0
10    CONTINUE
      CALL CLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
      IF( KASE /= 0 ) THEN
         IF( KASE == KASE1 ) THEN
!
!              Multiply by inv(A).
!
            CALL CLATBS( UPLO, 'No transpose', DIAG, NORMIN, N, KD, &
                         AB, LDAB, WORK, SCALE, RWORK, INFO )
         ELSE
!
!              Multiply by inv(A**H).
!
            CALL CLATBS( UPLO, 'Conjugate transpose', DIAG, NORMIN, &
                         N, KD, AB, LDAB, WORK, SCALE, RWORK, INFO )
         END IF
         NORMIN = 'Y'
!
!           Multiply by 1/SCALE if doing so will not cause overflow.
!
         IF( SCALE /= 1.0E+0 ) THEN
            IX = ICAMAX( N, WORK, 1 )
            XNORM = CABS1( WORK( IX ) )
            IF( SCALE < XNORM*SMLNUM .OR. SCALE == 0.0E+0 ) GO TO 20
            CALL CSRSCL( N, SCALE, WORK, 1 )
         END IF
         GO TO 10
      END IF
!
!        Compute the estimate of the reciprocal condition number.
!
      IF( AINVNM /= 0.0E+0 ) RCOND = ( 1.0E+0 / ANORM ) / AINVNM
   END IF
!
20 CONTINUE
   RETURN
!
!     End of CTBCON
!
END
