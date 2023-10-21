!> \brief \b ZHBTRD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHBTRD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhbtrd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhbtrd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhbtrd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ,
!                          WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO, VECT
!       INTEGER            INFO, KD, LDAB, LDQ, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), E( * )
!       COMPLEX*16         AB( LDAB, * ), Q( LDQ, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHBTRD reduces a complex Hermitian band matrix A to real symmetric
!> tridiagonal form T by a unitary similarity transformation:
!> Q**H * A * Q = T.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] VECT
!> \verbatim
!>          VECT is CHARACTER*1
!>          = 'N':  do not form Q;
!>          = 'V':  form Q;
!>          = 'U':  update a matrix X, by forming X*Q.
!> \endverbatim
!>
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
!>          AB is COMPLEX*16 array, dimension (LDAB,N)
!>          On entry, the upper or lower triangle of the Hermitian band
!>          matrix A, stored in the first KD+1 rows of the array.  The
!>          j-th column of A is stored in the j-th column of the array AB
!>          as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>          On exit, the diagonal elements of AB are overwritten by the
!>          diagonal elements of the tridiagonal matrix T; if KD > 0, the
!>          elements on the first superdiagonal (if UPLO = 'U') or the
!>          first subdiagonal (if UPLO = 'L') are overwritten by the
!>          off-diagonal elements of T; the rest of AB is overwritten by
!>          values generated during the reduction.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The diagonal elements of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          The off-diagonal elements of the tridiagonal matrix T:
!>          E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDQ,N)
!>          On entry, if VECT = 'U', then Q must contain an N-by-N
!>          matrix X; if VECT = 'N' or 'V', then Q need not be set.
!>
!>          On exit:
!>          if VECT = 'V', Q contains the N-by-N unitary matrix Q;
!>          if VECT = 'U', Q contains the product X*Q;
!>          if VECT = 'N', the array Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.
!>          LDQ >= 1, and LDQ >= N if VECT = 'V' or 'U'.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (N)
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
!> \ingroup hbtrd
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Modified by Linda Kaufman, Bell Labs.
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE ZHBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ, &
                      WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO, VECT
   INTEGER            INFO, KD, LDAB, LDQ, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   D( * ), E( * )
   COMPLEX*16         AB( LDAB, * ), Q( LDQ, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO
   PARAMETER          ( ZERO = 0.0D+0 )
   COMPLEX*16         CZERO, CONE
   PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), &
                      CONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
   LOGICAL            INITQ, UPPER, WANTQ
   INTEGER            I, I2, IBL, INCA, INCX, IQAEND, IQB, IQEND, J, &
                      J1, J1END, J1INC, J2, JEND, JIN, JINC, K, KD1, &
                      KDM1, KDN, L, LAST, LEND, NQ, NR, NRT
   DOUBLE PRECISION   ABST
   COMPLEX*16         T, TEMP
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA, ZLACGV, ZLAR2V, ZLARGV, ZLARTG, ZLARTV, &
                      ZLASET, ZROT, ZSCAL
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, DBLE, DCONJG, MAX, MIN
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
   INITQ = LSAME( VECT, 'V' )
   WANTQ = INITQ .OR. LSAME( VECT, 'U' )
   UPPER = LSAME( UPLO, 'U' )
   KD1 = KD + 1
   KDM1 = KD - 1
   INCX = LDAB - 1
   IQEND = 1
!
   INFO = 0
   IF( .NOT.WANTQ .AND. .NOT.LSAME( VECT, 'N' ) ) THEN
      INFO = -1
   ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -2
   ELSE IF( N < 0 ) THEN
      INFO = -3
   ELSE IF( KD < 0 ) THEN
      INFO = -4
   ELSE IF( LDAB < KD1 ) THEN
      INFO = -6
   ELSE IF( LDQ < MAX( 1, N ) .AND. WANTQ ) THEN
      INFO = -10
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'ZHBTRD', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) &
      RETURN
!
!     Initialize Q to the unit matrix, if needed
!
   IF( INITQ ) &
      CALL ZLASET( 'Full', N, N, CZERO, CONE, Q, LDQ )
!
!     Wherever possible, plane rotations are generated and applied in
!     vector operations of length NR over the index set J1:J2:KD1.
!
!     The real cosines and complex sines of the plane rotations are
!     stored in the arrays D and WORK.
!
   INCA = KD1*LDAB
   KDN = MIN( N-1, KD )
   IF( UPPER ) THEN
!
      IF( KD > 1 ) THEN
!
!           Reduce to complex Hermitian tridiagonal form, working with
!           the upper triangle
!
         NR = 0
         J1 = KDN + 2
         J2 = 1
!
         AB( KD1, 1 ) = DBLE( AB( KD1, 1 ) )
         DO I = 1, N - 2
!
!              Reduce i-th row of matrix to tridiagonal form
!
            DO K = KDN + 1, 2, -1
               J1 = J1 + KDN
               J2 = J2 + KDN
!
               IF( NR > 0 ) THEN
!
!                    generate plane rotations to annihilate nonzero
!                    elements which have been created outside the band
!
                  CALL ZLARGV( NR, AB( 1, J1-1 ), INCA, WORK( J1 ), &
                               KD1, D( J1 ), KD1 )
!
!                    apply rotations from the right
!
!
!                    Dependent on the the number of diagonals either
!                    ZLARTV or ZROT is used
!
                  IF( NR >= 2*KD-1 ) THEN
                     DO L = 1, KD - 1
                        CALL ZLARTV( NR, AB( L+1, J1-1 ), INCA, &
                                     AB( L, J1 ), INCA, D( J1 ), &
                                     WORK( J1 ), KD1 )
                     ENDDO
!
                  ELSE
                     JEND = J1 + ( NR-1 )*KD1
                     DO JINC = J1, JEND, KD1
                        CALL ZROT( KDM1, AB( 2, JINC-1 ), 1, &
                                   AB( 1, JINC ), 1, D( JINC ), &
                                   WORK( JINC ) )
                     ENDDO
                  END IF
               END IF
!
!
               IF( K > 2 ) THEN
                  IF( K <= N-I+1 ) THEN
!
!                       generate plane rotation to annihilate a(i,i+k-1)
!                       within the band
!
                     CALL ZLARTG( AB( KD-K+3, I+K-2 ), &
                                  AB( KD-K+2, I+K-1 ), D( I+K-1 ), &
                                  WORK( I+K-1 ), TEMP )
                     AB( KD-K+3, I+K-2 ) = TEMP
!
!                       apply rotation from the right
!
                     CALL ZROT( K-3, AB( KD-K+4, I+K-2 ), 1, &
                                AB( KD-K+3, I+K-1 ), 1, D( I+K-1 ), &
                                WORK( I+K-1 ) )
                  END IF
                  NR = NR + 1
                  J1 = J1 - KDN - 1
               END IF
!
!                 apply plane rotations from both sides to diagonal
!                 blocks
!
               IF( NR > 0 ) &
                  CALL ZLAR2V( NR, AB( KD1, J1-1 ), AB( KD1, J1 ), &
                               AB( KD, J1 ), INCA, D( J1 ), &
                               WORK( J1 ), KD1 )
!
!                 apply plane rotations from the left
!
               IF( NR > 0 ) THEN
                  CALL ZLACGV( NR, WORK( J1 ), KD1 )
                  IF( 2*KD-1 < NR ) THEN
!
!                    Dependent on the the number of diagonals either
!                    ZLARTV or ZROT is used
!
                     DO L = 1, KD - 1
                        IF( J2+L > N ) THEN
                           NRT = NR - 1
                        ELSE
                           NRT = NR
                        END IF
                        IF( NRT > 0 ) &
                           CALL ZLARTV( NRT, AB( KD-L, J1+L ), INCA, &
                                        AB( KD-L+1, J1+L ), INCA, &
                                        D( J1 ), WORK( J1 ), KD1 )
                     ENDDO
                  ELSE
                     J1END = J1 + KD1*( NR-2 )
                     IF( J1END >= J1 ) THEN
                        DO JIN = J1, J1END, KD1
                           CALL ZROT( KD-1, AB( KD-1, JIN+1 ), INCX, &
                                      AB( KD, JIN+1 ), INCX, &
                                      D( JIN ), WORK( JIN ) )
                        ENDDO
                     END IF
                     LEND = MIN( KDM1, N-J2 )
                     LAST = J1END + KD1
                     IF( LEND > 0 ) &
                        CALL ZROT( LEND, AB( KD-1, LAST+1 ), INCX, &
                                   AB( KD, LAST+1 ), INCX, D( LAST ), &
                                   WORK( LAST ) )
                  END IF
               END IF
!
               IF( WANTQ ) THEN
!
!                    accumulate product of plane rotations in Q
!
                  IF( INITQ ) THEN
!
!                 take advantage of the fact that Q was
!                 initially the Identity matrix
!
                     IQEND = MAX( IQEND, J2 )
                     I2 = MAX( 0, K-3 )
                     IQAEND = 1 + I*KD
                     IF( K == 2 ) &
                        IQAEND = IQAEND + KD
                     IQAEND = MIN( IQAEND, IQEND )
                     DO J = J1, J2, KD1
                        IBL = I - I2 / KDM1
                        I2 = I2 + 1
                        IQB = MAX( 1, J-IBL )
                        NQ = 1 + IQAEND - IQB
                        IQAEND = MIN( IQAEND+KD, IQEND )
                        CALL ZROT( NQ, Q( IQB, J-1 ), 1, Q( IQB, J ), &
                                   1, D( J ), DCONJG( WORK( J ) ) )
                     ENDDO
                  ELSE
!
                     DO J = J1, J2, KD1
                        CALL ZROT( N, Q( 1, J-1 ), 1, Q( 1, J ), 1, &
                                   D( J ), DCONJG( WORK( J ) ) )
                     ENDDO
                  END IF
!
               END IF
!
               IF( J2+KDN > N ) THEN
!
!                    adjust J2 to keep within the bounds of the matrix
!
                  NR = NR - 1
                  J2 = J2 - KDN - 1
               END IF
!
               DO J = J1, J2, KD1
!
!                    create nonzero element a(j-1,j+kd) outside the band
!                    and store it in WORK
!
                  WORK( J+KD ) = WORK( J )*AB( 1, J+KD )
                  AB( 1, J+KD ) = D( J )*AB( 1, J+KD )
               ENDDO
            ENDDO
         ENDDO
      END IF
!
      IF( KD > 0 ) THEN
!
!           make off-diagonal elements real and copy them to E
!
         DO I = 1, N - 1
            T = AB( KD, I+1 )
            ABST = ABS( T )
            AB( KD, I+1 ) = ABST
            E( I ) = ABST
            IF( ABST /= ZERO ) THEN
               T = T / ABST
            ELSE
               T = CONE
            END IF
            IF( I < N-1 ) &
               AB( KD, I+2 ) = AB( KD, I+2 )*T
            IF( WANTQ ) THEN
               CALL ZSCAL( N, DCONJG( T ), Q( 1, I+1 ), 1 )
            END IF
            ENDDO
      ELSE
!
!           set E to zero if original matrix was diagonal
!
         DO I = 1, N - 1
            E( I ) = ZERO
            ENDDO
      END IF
!
!        copy diagonal elements to D
!
      DO I = 1, N
         D( I ) = DBLE( AB( KD1, I ) )
         ENDDO
!
   ELSE
!
      IF( KD > 1 ) THEN
!
!           Reduce to complex Hermitian tridiagonal form, working with
!           the lower triangle
!
         NR = 0
         J1 = KDN + 2
         J2 = 1
!
         AB( 1, 1 ) = DBLE( AB( 1, 1 ) )
         DO I = 1, N - 2
!
!              Reduce i-th column of matrix to tridiagonal form
!
            DO K = KDN + 1, 2, -1
               J1 = J1 + KDN
               J2 = J2 + KDN
!
               IF( NR > 0 ) THEN
!
!                    generate plane rotations to annihilate nonzero
!                    elements which have been created outside the band
!
                  CALL ZLARGV( NR, AB( KD1, J1-KD1 ), INCA, &
                               WORK( J1 ), KD1, D( J1 ), KD1 )
!
!                    apply plane rotations from one side
!
!
!                    Dependent on the the number of diagonals either
!                    ZLARTV or ZROT is used
!
                  IF( NR > 2*KD-1 ) THEN
                     DO L = 1, KD - 1
                        CALL ZLARTV( NR, AB( KD1-L, J1-KD1+L ), INCA, &
                                     AB( KD1-L+1, J1-KD1+L ), INCA, &
                                     D( J1 ), WORK( J1 ), KD1 )
                        ENDDO
                  ELSE
                     JEND = J1 + KD1*( NR-1 )
                     DO JINC = J1, JEND, KD1
                        CALL ZROT( KDM1, AB( KD, JINC-KD ), INCX, &
                                   AB( KD1, JINC-KD ), INCX, &
                                   D( JINC ), WORK( JINC ) )
                        ENDDO
                  END IF
!
               END IF
!
               IF( K > 2 ) THEN
                  IF( K <= N-I+1 ) THEN
!
!                       generate plane rotation to annihilate a(i+k-1,i)
!                       within the band
!
                     CALL ZLARTG( AB( K-1, I ), AB( K, I ), &
                                  D( I+K-1 ), WORK( I+K-1 ), TEMP )
                     AB( K-1, I ) = TEMP
!
!                       apply rotation from the left
!
                     CALL ZROT( K-3, AB( K-2, I+1 ), LDAB-1, &
                                AB( K-1, I+1 ), LDAB-1, D( I+K-1 ), &
                                WORK( I+K-1 ) )
                  END IF
                  NR = NR + 1
                  J1 = J1 - KDN - 1
               END IF
!
!                 apply plane rotations from both sides to diagonal
!                 blocks
!
               IF( NR > 0 ) &
                  CALL ZLAR2V( NR, AB( 1, J1-1 ), AB( 1, J1 ), &
                               AB( 2, J1-1 ), INCA, D( J1 ), &
                               WORK( J1 ), KD1 )
!
!                 apply plane rotations from the right
!
!
!                    Dependent on the the number of diagonals either
!                    ZLARTV or ZROT is used
!
               IF( NR > 0 ) THEN
                  CALL ZLACGV( NR, WORK( J1 ), KD1 )
                  IF( NR > 2*KD-1 ) THEN
                     DO L = 1, KD - 1
                        IF( J2+L > N ) THEN
                           NRT = NR - 1
                        ELSE
                           NRT = NR
                        END IF
                        IF( NRT > 0 ) &
                           CALL ZLARTV( NRT, AB( L+2, J1-1 ), INCA, &
                                        AB( L+1, J1 ), INCA, D( J1 ), &
                                        WORK( J1 ), KD1 )
                        ENDDO
                  ELSE
                     J1END = J1 + KD1*( NR-2 )
                     IF( J1END >= J1 ) THEN
                        DO J1INC = J1, J1END, KD1
                           CALL ZROT( KDM1, AB( 3, J1INC-1 ), 1, &
                                      AB( 2, J1INC ), 1, D( J1INC ), &
                                      WORK( J1INC ) )
                           ENDDO
                     END IF
                     LEND = MIN( KDM1, N-J2 )
                     LAST = J1END + KD1
                     IF( LEND > 0 ) &
                        CALL ZROT( LEND, AB( 3, LAST-1 ), 1, &
                                   AB( 2, LAST ), 1, D( LAST ), &
                                   WORK( LAST ) )
                  END IF
               END IF
!
!
!
               IF( WANTQ ) THEN
!
!                    accumulate product of plane rotations in Q
!
                  IF( INITQ ) THEN
!
!                 take advantage of the fact that Q was
!                 initially the Identity matrix
!
                     IQEND = MAX( IQEND, J2 )
                     I2 = MAX( 0, K-3 )
                     IQAEND = 1 + I*KD
                     IF( K == 2 ) &
                        IQAEND = IQAEND + KD
                     IQAEND = MIN( IQAEND, IQEND )
                     DO J = J1, J2, KD1
                        IBL = I - I2 / KDM1
                        I2 = I2 + 1
                        IQB = MAX( 1, J-IBL )
                        NQ = 1 + IQAEND - IQB
                        IQAEND = MIN( IQAEND+KD, IQEND )
                        CALL ZROT( NQ, Q( IQB, J-1 ), 1, Q( IQB, J ), &
                                   1, D( J ), WORK( J ) )
                        ENDDO
                  ELSE
!
                     DO J = J1, J2, KD1
                        CALL ZROT( N, Q( 1, J-1 ), 1, Q( 1, J ), 1, &
                                   D( J ), WORK( J ) )
                        ENDDO
                  END IF
               END IF
!
               IF( J2+KDN > N ) THEN
!
!                    adjust J2 to keep within the bounds of the matrix
!
                  NR = NR - 1
                  J2 = J2 - KDN - 1
               END IF
!
               DO J = J1, J2, KD1
!
!                    create nonzero element a(j+kd,j-1) outside the
!                    band and store it in WORK
!
                  WORK( J+KD ) = WORK( J )*AB( KD1, J )
                  AB( KD1, J ) = D( J )*AB( KD1, J )
                  ENDDO
               ENDDO
            ENDDO
      END IF
!
      IF( KD > 0 ) THEN
!
!           make off-diagonal elements real and copy them to E
!
         DO I = 1, N - 1
            T = AB( 2, I )
            ABST = ABS( T )
            AB( 2, I ) = ABST
            E( I ) = ABST
            IF( ABST /= ZERO ) THEN
               T = T / ABST
            ELSE
               T = CONE
            END IF
            IF( I < N-1 ) &
               AB( 2, I+1 ) = AB( 2, I+1 )*T
            IF( WANTQ ) THEN
               CALL ZSCAL( N, T, Q( 1, I+1 ), 1 )
            END IF
            ENDDO
      ELSE
!
!           set E to zero if original matrix was diagonal
!
         DO I = 1, N - 1
            E( I ) = ZERO
            ENDDO
      END IF
!
!        copy diagonal elements to D
!
      DO I = 1, N
         D( I ) = DBLE( AB( 1, I ) )
         ENDDO
   END IF
!
   RETURN
!
!     End of ZHBTRD
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

