!> \brief \b CPBTRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CPBTRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbtrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbtrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbtrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CPBTRF( UPLO, N, KD, AB, LDAB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, KD, LDAB, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            AB( LDAB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CPBTRF computes the Cholesky factorization of a complex Hermitian
!> positive definite band matrix A.
!>
!> The factorization has the form
!>    A = U**H * U,  if UPLO = 'U', or
!>    A = L  * L**H,  if UPLO = 'L',
!> where U is an upper triangular matrix and L is lower triangular.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
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
!>          The number of superdiagonals of the matrix A if UPLO = 'U',
!>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is COMPLEX array, dimension (LDAB,N)
!>          On entry, the upper or lower triangle of the Hermitian band
!>          matrix A, stored in the first KD+1 rows of the array.  The
!>          j-th column of A is stored in the j-th column of the array AB
!>          as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>
!>          On exit, if INFO = 0, the triangular factor U or L from the
!>          Cholesky factorization A = U**H*U or A = L*L**H of the band
!>          matrix A, in the same storage format as A.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, the leading principal minor of order i
!>                is not positive, and the factorization could not be
!>                completed.
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
!> \ingroup pbtrf
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The band storage scheme is illustrated by the following example, when
!>  N = 6, KD = 2, and UPLO = 'U':
!>
!>  On entry:                       On exit:
!>
!>      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
!>      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!>     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!>
!>  Similarly, if UPLO = 'L' the format of A is as follows:
!>
!>  On entry:                       On exit:
!>
!>     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
!>     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
!>     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
!>
!>  Array elements marked * are not used by the routine.
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!>  Peter Mayes and Giuseppe Radicati, IBM ECSEC, Rome, March 23, 1989
!
!  =====================================================================
   SUBROUTINE CPBTRF( UPLO, N, KD, AB, LDAB, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INFO, KD, LDAB, N
!     ..
!     .. Array Arguments ..
   COMPLEX            AB( LDAB, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            NBMAX, LDWORK
   PARAMETER          ( NBMAX = 32, LDWORK = NBMAX+1 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, I2, I3, IB, II, J, JJ, NB
!     ..
!     .. Local Arrays ..
   COMPLEX            WORK( LDWORK, NBMAX )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ILAENV
   EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMM, CHERK, CPBTF2, CPOTF2, CTRSM, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   IF( ( .NOT.LSAME( UPLO, 'U' ) ) .AND. &
       ( .NOT.LSAME( UPLO, 'L' ) ) ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( KD < 0 ) THEN
      INFO = -3
   ELSE IF( LDAB < KD+1 ) THEN
      INFO = -5
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CPBTRF', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) RETURN
!
!     Determine the block size for this environment
!
   NB = ILAENV( 1, 'CPBTRF', UPLO, N, KD, -1, -1 )
!
!     The block size must not exceed the semi-bandwidth KD, and must not
!     exceed the limit set by the size of the local array WORK.
!
   NB = MIN( NB, NBMAX )
!
   IF( NB <= 1 .OR. NB > KD ) THEN
!
!        Use unblocked code
!
      CALL CPBTF2( UPLO, N, KD, AB, LDAB, INFO )
   ELSE
!
!        Use blocked code
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!           Compute the Cholesky factorization of a Hermitian band
!           matrix, given the upper triangle of the matrix in band
!           storage.
!
!           Zero the upper triangle of the work array.
!
         DO J = 1, NB
            WORK(1:J-1,J) = 0.0E+0
         ENDDO
!
!           Process the band matrix one diagonal block at a time.
!
         DO I = 1, N, NB
            IB = MIN( NB, N-I+1 )
!
!              Factorize the diagonal block
!
            CALL CPOTF2( UPLO, IB, AB( KD+1, I ), LDAB-1, II )
            IF( II /= 0 ) THEN
               INFO = I + II - 1
               GO TO 150
            END IF
            IF( I+IB <= N ) THEN
!
!                 Update the relevant part of the trailing submatrix.
!                 If A11 denotes the diagonal block which has just been
!                 factorized, then we need to update the remaining
!                 blocks in the diagram:
!
!                    A11   A12   A13
!                          A22   A23
!                                A33
!
!                 The numbers of rows and columns in the partitioning
!                 are IB, I2, I3 respectively. The blocks A12, A22 and
!                 A23 are empty if IB = KD. The upper triangle of A13
!                 lies outside the band.
!
               I2 = MIN( KD-IB, N-I-IB+1 )
               I3 = MIN( IB, N-I-KD+1 )
!
               IF( I2 > 0 ) THEN
!
!                    Update A12
!
                  CALL CTRSM( 'Left', 'Upper', 'Conjugate transpose', &
                              'Non-unit', IB, I2, (1.0E+0,0.0E+0), &
                              AB( KD+1, I ), LDAB-1, &
                              AB( KD+1-IB, I+IB ), LDAB-1 )
!
!                    Update A22
!
                  CALL CHERK( 'Upper', 'Conjugate transpose', I2, IB, &
                              -1.0E+0, AB( KD+1-IB, I+IB ), LDAB-1, 1.0E+0, &
                              AB( KD+1, I+IB ), LDAB-1 )
               END IF
!
               IF( I3 > 0 ) THEN
!
!                    Copy the lower triangle of A13 into the work array.
!
                  DO JJ = 1, I3
                     WORK(JJ:IB,JJ) = AB(1:1-JJ+IB,JJ+I+KD-1)
                  ENDDO
!
!                    Update A13 (in the work array).
!
                  CALL CTRSM( 'Left', 'Upper', 'Conjugate transpose', &
                              'Non-unit', IB, I3, (1.0E+0,0.0E+0), &
                              AB( KD+1, I ), LDAB-1, WORK, LDWORK )
!
!                    Update A23
!
                  IF( I2 > 0 ) &
                     CALL CGEMM( 'Conjugate transpose', &
                                 'No transpose', I2, I3, IB, -(1.0E+0,0.0E+0), &
                                 AB( KD+1-IB, I+IB ), LDAB-1, WORK, &
                                 LDWORK, (1.0E+0,0.0E+0), AB( 1+IB, I+KD ), &
                                 LDAB-1 )
!
!                    Update A33
!
                  CALL CHERK( 'Upper', 'Conjugate transpose', I3, IB, &
                              -1.0E+0, WORK, LDWORK, 1.0E+0, &
                              AB( KD+1, I+KD ), LDAB-1 )
!
!                    Copy the lower triangle of A13 back into place.
!
                  DO JJ = 1, I3
                     AB(1:1-JJ+IB, JJ+I+KD-1) = WORK(JJ:IB,JJ)
                  ENDDO
               END IF
            END IF
         ENDDO
      ELSE
!
!           Compute the Cholesky factorization of a Hermitian band
!           matrix, given the lower triangle of the matrix in band
!           storage.
!
!           Zero the lower triangle of the work array.
!
         DO J = 1, NB
            WORK(J+1:NB,J) = 0.0E+0
         ENDDO
!
!           Process the band matrix one diagonal block at a time.
!
         DO I = 1, N, NB
            IB = MIN( NB, N-I+1 )
!
!              Factorize the diagonal block
!
            CALL CPOTF2( UPLO, IB, AB( 1, I ), LDAB-1, II )
            IF( II /= 0 ) THEN
               INFO = I + II - 1
               GO TO 150
            END IF
            IF( I+IB <= N ) THEN
!
!                 Update the relevant part of the trailing submatrix.
!                 If A11 denotes the diagonal block which has just been
!                 factorized, then we need to update the remaining
!                 blocks in the diagram:
!
!                    A11
!                    A21   A22
!                    A31   A32   A33
!
!                 The numbers of rows and columns in the partitioning
!                 are IB, I2, I3 respectively. The blocks A21, A22 and
!                 A32 are empty if IB = KD. The lower triangle of A31
!                 lies outside the band.
!
               I2 = MIN( KD-IB, N-I-IB+1 )
               I3 = MIN( IB, N-I-KD+1 )
!
               IF( I2 > 0 ) THEN
!
!                    Update A21
!
                  CALL CTRSM( 'Right', 'Lower', &
                              'Conjugate transpose', 'Non-unit', I2, &
                              IB, (1.0E+0,0.0E+0), AB( 1, I ), LDAB-1, &
                              AB( 1+IB, I ), LDAB-1 )
!
!                    Update A22
!
                  CALL CHERK( 'Lower', 'No transpose', I2, IB, -1.0E+0, &
                              AB( 1+IB, I ), LDAB-1, 1.0E+0, &
                              AB( 1, I+IB ), LDAB-1 )
               END IF
!
               IF( I3 > 0 ) THEN
!
!                    Copy the upper triangle of A31 into the work array.
!
                  DO JJ = 1, IB
                     DO II = 1, MIN( JJ, I3 )
                        WORK( II, JJ ) = AB( KD+1-JJ+II, JJ+I-1 )
                     ENDDO
                  ENDDO
!
!                    Update A31 (in the work array).
!
                  CALL CTRSM( 'Right', 'Lower', &
                              'Conjugate transpose', 'Non-unit', I3, &
                              IB, (1.0E+0,0.0E+0), AB( 1, I ), LDAB-1, WORK, &
                              LDWORK )
!
!                    Update A32
!
                  IF( I2 > 0 ) &
                     CALL CGEMM( 'No transpose', &
                                 'Conjugate transpose', I3, I2, IB, &
                                 -(1.0E+0,0.0E+0), WORK, LDWORK, AB( 1+IB, I ), &
                                 LDAB-1, (1.0E+0,0.0E+0), AB( 1+KD-IB, I+IB ), &
                                 LDAB-1 )
!
!                    Update A33
!
                  CALL CHERK( 'Lower', 'No transpose', I3, IB, -1.0E+0, &
                              WORK, LDWORK, 1.0E+0, AB( 1, I+KD ), &
                              LDAB-1 )
!
!                    Copy the upper triangle of A31 back into place.
!
                  DO JJ = 1, IB
                     DO II = 1, MIN( JJ, I3 )
                        AB( KD+1-JJ+II, JJ+I-1 ) = WORK( II, JJ )
                     ENDDO
                  ENDDO
               END IF
            END IF
         ENDDO
      END IF
   END IF
   RETURN
!
  150 CONTINUE
   RETURN
!
!     End of CPBTRF
!
END
