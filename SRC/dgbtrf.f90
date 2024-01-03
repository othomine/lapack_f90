!> \brief \b DGBTRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGBTRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbtrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbtrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbtrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, KL, KU, LDAB, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   AB( LDAB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGBTRF computes an LU factorization of a real m-by-n band matrix A
!> using partial pivoting with row interchanges.
!>
!> This is the blocked version of the algorithm, calling Level 3 BLAS.
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
!>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
!>          On entry, the matrix A in band storage, in rows KL+1 to
!>          2*KL+KU+1; rows 1 to KL of the array need not be set.
!>          The j-th column of A is stored in the j-th column of the
!>          array AB as follows:
!>          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!>
!>          On exit, details of the factorization: U is stored as an
!>          upper triangular band matrix with KL+KU superdiagonals in
!>          rows 1 to KL+KU+1, and the multipliers used during the
!>          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!>          See below for further details.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (min(M,N))
!>          The pivot indices; for 1 <= i <= min(M,N), row i of the
!>          matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!>               has been completed, but the factor U is exactly
!>               singular, and division by zero will occur if it is used
!>               to solve a system of equations.
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
!> \ingroup gbtrf
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The band storage scheme is illustrated by the following example, when
!>  M = N = 6, KL = 2, KU = 1:
!>
!>  On entry:                       On exit:
!>
!>      *    *    *    +    +    +       *    *    *   u14  u25  u36
!>      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!>      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!>     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!>     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!>     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!>
!>  Array elements marked * are not used by the routine; elements marked
!>  + need not be set on entry, but are required by the routine to store
!>  elements of U because of fill-in resulting from the row interchanges.
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * )
   DOUBLE PRECISION   AB( LDAB, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            NBMAX, LDWORK
   PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP, &
                      JU, K2, KM, KV, NB, NW
   DOUBLE PRECISION   TEMP
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   WORK13( LDWORK, NBMAX ), WORK31( LDWORK, NBMAX )
!     ..
!     .. External Functions ..
   INTEGER            IDAMAX, ILAENV
   EXTERNAL           IDAMAX, ILAENV
!     ..
!     .. External Subroutines ..
   EXTERNAL           DCOPY, DGBTF2, DGEMM, DGER, DLASWP, DSCAL, &
                      DSWAP, DTRSM, XERBLA
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in
!
   KV = KU + KL
!
!     Test the input parameters.
!
   INFO = 0
   IF( M < 0 ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( KL < 0 ) THEN
      INFO = -3
   ELSE IF( KU < 0 ) THEN
      INFO = -4
   ELSE IF( LDAB < KL+KV+1 ) THEN
      INFO = -6
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'DGBTRF', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( M == 0 .OR. N == 0 ) RETURN
!
!     Determine the block size for this environment
!
   NB = ILAENV( 1, 'DGBTRF', ' ', M, N, KL, KU )
!
!     The block size must not exceed the limit set by the size of the
!     local arrays WORK13 and WORK31.
!
   NB = MIN( NB, NBMAX )
!
   IF( NB <= 1 .OR. NB > KL ) THEN
!
!        Use unblocked code
!
      CALL DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
   ELSE
!
!        Use blocked code
!
!        Zero the superdiagonal elements of the work array WORK13
!
      DO J = 1, NB
         WORK13(1:J-1,J) = 0.0D0
      ENDDO
!
!        Zero the subdiagonal elements of the work array WORK31
!
      DO J = 1, NB
         WORK31(J+1:NB,J) = 0.0D0
      ENDDO
!
!        Gaussian elimination with partial pivoting
!
!        Set fill-in elements in columns KU+2 to KV to zero
!
      DO J = KU + 2, MIN( KV, N )
         AB(KV-J+2:KL,J) = 0.0D0
      ENDDO
!
!        JU is the index of the last column affected by the current
!        stage of the factorization
!
      JU = 1
!
      DO J = 1, MIN( M, N ), NB
         JB = MIN( NB, MIN( M, N )-J+1 )
!
!           The active part of the matrix is partitioned
!
!              A11   A12   A13
!              A21   A22   A23
!              A31   A32   A33
!
!           Here A11, A21 and A31 denote the current block of JB columns
!           which is about to be factorized. The number of rows in the
!           partitioning are JB, I2, I3 respectively, and the numbers
!           of columns are JB, J2, J3. The superdiagonal elements of A13
!           and the subdiagonal elements of A31 lie outside the band.
!
         I2 = MIN( KL-JB, M-J-JB+1 )
         I3 = MIN( JB, M-J-KL+1 )
!
!           J2 and J3 are computed after JU has been updated.
!
!           Factorize the current block of JB columns
!
         DO JJ = J, J + JB - 1
!
!              Set fill-in elements in column JJ+KV to zero
!
            IF( JJ+KV <= N ) AB(1:KL,JJ+KV) = 0.0D0
!
!              Find pivot and test for singularity. KM is the number of
!              subdiagonal elements in the current column.
!
            KM = MIN( KL, M-JJ )
            JP = IDAMAX( KM+1, AB( KV+1, JJ ), 1 )
            IPIV( JJ ) = JP + JJ - J
            IF( AB( KV+JP, JJ ) /= 0.0D0 ) THEN
               JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
               IF( JP /= 1 ) THEN
!
!                    Apply interchange to columns J to J+JB-1
!
                  IF( JP+JJ-1 < J+KL ) THEN
!
                     CALL DSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1, &
                                 AB( KV+JP+JJ-J, J ), LDAB-1 )
                  ELSE
!
!                       The interchange affects columns J to JJ-1 of A31
!                       which are stored in the work array WORK31
!
                     CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                                 WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                     CALL DSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1, &
                                 AB( KV+JP, JJ ), LDAB-1 )
                  END IF
               END IF
!
!                 Compute multipliers
!
               CALL DSCAL( KM, 1.0D0 / AB( KV+1, JJ ), AB( KV+2, JJ ), 1 )
!
!                 Update trailing submatrix within the band and within
!                 the current block. JM is the index of the last column
!                 which needs to be updated.
!
               JM = MIN( JU, J+JB-1 )
               IF( JM > JJ ) &
                  CALL DGER( KM, JM-JJ, -1.0D0, AB( KV+2, JJ ), 1, &
                             AB( KV, JJ+1 ), LDAB-1, &
                             AB( KV+1, JJ+1 ), LDAB-1 )
            ELSE
!
!                 If pivot is zero, set INFO to the index of the pivot
!                 unless a zero pivot has already been found.
!
               IF( INFO == 0 ) INFO = JJ
            END IF
!
!              Copy current column of A31 into the work array WORK31
!
            NW = MIN( JJ-J+1, I3 )
            IF( NW > 0 ) &
               CALL DCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1, WORK31( 1, JJ-J+1 ), 1 )
         ENDDO
         IF( J+JB <= N ) THEN
!
!              Apply the row interchanges to the other blocks.
!
            J2 = MIN( JU-J+1, KV ) - JB
            J3 = MAX( 0, JU-J-KV+1 )
!
!              Use DLASWP to apply the row interchanges to A12, A22, and
!              A32.
!
            CALL DLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB, IPIV( J ), 1 )
!
!              Adjust the pivot indices.
!
            IPIV(J:J+JB-1) = IPIV(J:J+JB-1) + J - 1
!
!              Apply the row interchanges to A13, A23, and A33
!              columnwise.
!
            K2 = J - 1 + JB + J2
            DO I = 1, J3
               JJ = K2 + I
               DO II = J + I - 1, J + JB - 1
                  IP = IPIV( II )
                  IF( IP /= II ) THEN
                     TEMP = AB( KV+1+II-JJ, JJ )
                     AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                     AB( KV+1+IP-JJ, JJ ) = TEMP
                  END IF
               ENDDO
            ENDDO
!
!              Update the relevant part of the trailing submatrix
!
            IF( J2 > 0 ) THEN
!
!                 Update A12
!
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', &
                           JB, J2, 1.0D0, AB( KV+1, J ), LDAB-1, &
                           AB( KV+1-JB, J+JB ), LDAB-1 )
!
               IF( I2 > 0 ) THEN
!
!                    Update A22
!
                  CALL DGEMM( 'No transpose', 'No transpose', I2, J2, &
                              JB, -1.0D0, AB( KV+1+JB, J ), LDAB-1, &
                              AB( KV+1-JB, J+JB ), LDAB-1, 1.0D0, &
                              AB( KV+1, J+JB ), LDAB-1 )
               END IF
!
               IF( I3 > 0 ) THEN
!
!                    Update A32
!
                  CALL DGEMM( 'No transpose', 'No transpose', I3, J2, &
                              JB, -1.0D0, WORK31, LDWORK, &
                              AB( KV+1-JB, J+JB ), LDAB-1, 1.0D0, &
                              AB( KV+KL+1-JB, J+JB ), LDAB-1 )
               END IF
            END IF
!
            IF( J3 > 0 ) THEN
!
!                 Copy the lower triangle of A13 into the work array
!                 WORK13
!
               DO JJ = 1, J3
                  WORK13(JJ:JB, JJ ) = AB(1:JB-JJ+1, JJ+J+KV-1 )
               ENDDO
!
!                 Update A13 in the work array
!
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', &
                           JB, J3, 1.0D0, AB( KV+1, J ), LDAB-1, WORK13, LDWORK )
!
               IF( I2 > 0 ) THEN
!
!                    Update A23
!
                  CALL DGEMM( 'No transpose', 'No transpose', I2, J3, &
                              JB, -1.0D0, AB( KV+1+JB, J ), LDAB-1, &
                              WORK13, LDWORK, 1.0D0, AB( 1+JB, J+KV ), &
                              LDAB-1 )
               END IF
!
               IF( I3 > 0 ) THEN
!
!                    Update A33
!
                  CALL DGEMM( 'No transpose', 'No transpose', I3, J3, &
                              JB, -1.0D0, WORK31, LDWORK, WORK13, &
                              LDWORK, 1.0D0, AB( 1+KL, J+KV ), LDAB-1 )
               END IF
!
!                 Copy the lower triangle of A13 back into place
!
               DO JJ = 1, J3
                  AB(JJ-JJ+1:JB-JJ+1, JJ+J+KV-1 ) = WORK13(JJ:JB,JJ)
               ENDDO
            END IF
         ELSE
!
!              Adjust the pivot indices.
!
            IPIV(J:J+JB-1) = IPIV(J:J+JB-1) + J - 1
         END IF
!
!           Partially undo the interchanges in the current block to
!           restore the upper triangular form of A31 and copy the upper
!           triangle of A31 back into place
!
         DO JJ = J + JB - 1, J, -1
            JP = IPIV( JJ ) - JJ + 1
            IF( JP /= 1 ) THEN
!
!                 Apply interchange to columns J to JJ-1
!
               IF( JP+JJ-1 < J+KL ) THEN
!
!                    The interchange does not affect A31
!
                  CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                              AB( KV+JP+JJ-J, J ), LDAB-1 )
               ELSE
!
!                    The interchange does affect A31
!
                  CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                              WORK31( JP+JJ-J-KL, 1 ), LDWORK )
               END IF
            END IF
!
!              Copy the current column of A31 back into place
!
            NW = MIN( I3, JJ-J+1 )
            IF( NW > 0 ) AB(KV+KL+1-JJ+J:KV+KL-JJ+J+NW,JJ) = WORK31(1:NW,JJ-J+1)
         ENDDO
      ENDDO
   END IF
!
   RETURN
!
!     End of DGBTRF
!
END
