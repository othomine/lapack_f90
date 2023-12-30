!> \brief \b CLAHEF_AA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAHEF_AA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahef_aa.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahef_aa.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahef_aa.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAHEF_AA( UPLO, J1, M, NB, A, LDA, IPIV,
!                             H, LDH, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER    UPLO
!       INTEGER      J1, M, NB, LDA, LDH
!       ..
!       .. Array Arguments ..
!       INTEGER      IPIV( * )
!       COMPLEX      A( LDA, * ), H( LDH, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAHEF_AA factorizes a panel of a complex hermitian matrix A using
!> the Aasen's algorithm. The panel consists of a set of NB rows of A
!> when UPLO is U, or a set of NB columns when UPLO is L.
!>
!> In order to factorize the panel, the Aasen's algorithm requires the
!> last row, or column, of the previous panel. The first row, or column,
!> of A is set to be the first row, or column, of an identity matrix,
!> which is used to factorize the first panel.
!>
!> The resulting J-th row of U, or J-th column of L, is stored in the
!> (J-1)-th row, or column, of A (without the unit diagonals), while
!> the diagonal and subdiagonal of A are overwritten by those of T.
!>
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
!> \param[in] J1
!> \verbatim
!>          J1 is INTEGER
!>          The location of the first row, or column, of the panel
!>          within the submatrix of A, passed to this routine, e.g.,
!>          when called by CHETRF_AA, for the first panel, J1 is 1,
!>          while for the remaining panels, J1 is 2.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The dimension of the submatrix. M >= 0.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The dimension of the panel to be facotorized.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,M) for
!>          the first panel, while dimension (LDA,M+1) for the
!>          remaining panels.
!>
!>          On entry, A contains the last row, or column, of
!>          the previous panel, and the trailing submatrix of A
!>          to be factorized, except for the first panel, only
!>          the panel is passed.
!>
!>          On exit, the leading panel is factorized.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the row and column interchanges,
!>          the row and column k were interchanged with the row and
!>          column IPIV(k).
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX workspace, dimension (LDH,NB).
!>
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          The leading dimension of the workspace H. LDH >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX workspace, dimension (M).
!> \endverbatim
!>
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
!> \ingroup lahef_aa
!
!  =====================================================================
   SUBROUTINE CLAHEF_AA( UPLO, J1, M, NB, A, LDA, IPIV, &
                         H, LDH, WORK )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
   IMPLICIT NONE
!
!     .. Scalar Arguments ..
   CHARACTER    UPLO
   INTEGER      M, NB, J1, LDA, LDH
!     ..
!     .. Array Arguments ..
   INTEGER      IPIV( * )
   COMPLEX      A( LDA, * ), H( LDH, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
   INTEGER      J, K, K1, I1, I2, MJ
   COMPLEX      PIV, ALPHA, A_TMP( LDA ), H_TMP( NB )
!     ..
!     .. External Functions ..
   LOGICAL      LSAME
   INTEGER      ICAMAX, ILAENV
   EXTERNAL     LSAME, ILAENV, ICAMAX
!     ..
!     .. External Subroutines ..
   EXTERNAL     CLACGV, CGEMV, XERBLA
!     ..
!     .. Executable Statements ..
!
   J = 1
!
!     K1 is the first column of the panel to be factorized
!     i.e.,  K1 is 2 for the first block column, and 1 for the rest of the blocks
!
   K1 = (2-J1)+1
!
   IF( LSAME( UPLO, 'U' ) ) THEN
!
!        .....................................................
!        Factorize A as U**T*D*U using the upper triangle of A
!        .....................................................
!
 10      CONTINUE
      IF ( J > MIN(M, NB) ) GO TO 20
!
!        K is the column to be factorized
!         when being called from CHETRF_AA,
!         > for the first block column, J1 is 1, hence J1+J-1 is J,
!         > for the rest of the columns, J1 is 2, and J1+J-1 is J+1,
!
      K = J1+J-1
      IF( J == M ) THEN
!
!            Only need to compute T(J, J)
!
          MJ = 1
      ELSE
          MJ = M-J+1
      END IF
!
!        H(J:N, J) := A(J, J:N) - H(J:N, 1:(J-1)) * L(J1:(J-1), J),
!         where H(J:N, J) has been initialized to be A(J, J:N)
!
      IF( K > 2 ) THEN
!
!        K is the column to be factorized
!         > for the first block column, K is J, skipping the first two
!           columns
!         > for the rest of the columns, K is J+1, skipping only the
!           first column
!
         CALL CLACGV( J-K1, A( 1, J ), 1 )
         CALL CGEMV( 'No transpose', MJ, J-K1, &
                    -(1.0E+0,0.0E+0), H( J, K1 ), LDH, &
                          A( 1, J ), 1, &
                     (1.0E+0,0.0E+0), H( J, J ), 1 )
         CALL CLACGV( J-K1, A( 1, J ), 1 )
      END IF
!
!        Copy H(i:n, i) into WORK
!
      WORK(1:MJ) = H(J:J+MJ-1,J)
!
      IF( J > K1 ) THEN
!
!           Compute WORK := WORK - L(J-1, J:N) * T(J-1,J),
!            where A(J-1, J) stores T(J-1, J) and A(J-2, J:N) stores U(J-1, J:N)
!
         ALPHA = -CONJG( A( K-1, J ) )
         WORK(1:MJ) = WORK(1:MJ) + ALPHA*A(K-2,J:J+MJ-1)
      END IF
!
!        Set A(J, J) = T(J, J)
!
      A( K, J ) = REAL( WORK( 1 ) )
!
      IF( J < M ) THEN
!
!           Compute WORK(2:N) = T(J, J) L(J, (J+1):N)
!            where A(J, J) stores T(J, J) and A(J-1, (J+1):N) stores U(J, (J+1):N)
!
         IF( K > 1 ) THEN
            ALPHA = -A( K, J )
            WORK(2:1+M-J) = WORK(2:1+M-J) + ALPHA*A(K-1,J+1:M)
         ENDIF
!
!           Find max(|WORK(2:n)|)
!
         I2 = ICAMAX( M-J, WORK( 2 ), 1 ) + 1
         PIV = WORK( I2 )
!
!           Apply hermitian pivot
!
         IF( (I2 /= 2) .AND. (PIV /= 0) ) THEN
!
!              Swap WORK(I1) and WORK(I2)
!
            I1 = 2
            WORK( I2 ) = WORK( I1 )
            WORK( I1 ) = PIV
!
!              Swap A(I1, I1+1:N) with A(I1+1:N, I2)
!
            I1 = I1+J-1
            I2 = I2+J-1
            A_TMP(1:I2-I1-1) = A(J1+I1-1,I1+1:I2-1)
            A(J1+I1-1,I1+1:I2-1) = A(J1+I1:J1+I2-2,I2)
            A(J1+I1:J1+I2-2,I2) = A_TMP(1:I2-I1-1)
            CALL CLACGV( I2-I1, A( J1+I1-1, I1+1 ), LDA )
            CALL CLACGV( I2-I1-1, A( J1+I1, I2 ), 1 )
!
!              Swap A(I1, I2+1:N) with A(I2, I2+1:N)
!
            IF( I2 < M ) THEN
               A_TMP(1:M-I2) = A(J1+I1-1,I2+1:M)
               A(J1+I1-1,I2+1:M) = A(J1+I2-1,I2+1:M)
               A(J1+I2-1,I2+1:M) = A_TMP(1:M-I2)
            ENDIF
!
!              Swap A(I1, I1) with A(I2,I2)
!
            PIV = A( I1+J1-1, I1 )
            A( J1+I1-1, I1 ) = A( J1+I2-1, I2 )
            A( J1+I2-1, I2 ) = PIV
!
!              Swap H(I1, 1:J1) with H(I2, 1:J1)
!
            H_TMP(1:I1-1) = H(I1,1:I1-1)
            H(I1,1:I1-1) = H(I2,1:I1-1)
            H(I2,1:I1-1) = H_TMP(1:I1-1)
            IPIV( I1 ) = I2
!
            IF( I1 > (K1-1) ) THEN
!
!                 Swap L(1:I1-1, I1) with L(1:I1-1, I2),
!                  skipping the first column
!
               A_TMP(1:I1-K1+1) = A(1:1+I1-K1,I1)
               A(1:1+I1-K1,I1) = A(1:1+I1-K1,I2)
               A(1:1+I1-K1,I2) = A_TMP(1:I1-K1+1)
            END IF
         ELSE
            IPIV( J+1 ) = J+1
         ENDIF
!
!           Set A(J, J+1) = T(J, J+1)
!
         A( K, J+1 ) = WORK( 2 )
!
!
!              Copy A(J+1:N, J+1) into H(J:N, J),
!
         IF( J < NB ) H(J+1:M,J+1) = A(K+1,J+1:M)
!
!           Compute L(J+2, J+1) = WORK( 3:N ) / T(J, J+1),
!            where A(J, J+1) = T(J, J+1) and A(J+2:N, J) = L(J+2:N, J+1)
!
         IF( J < (M-1) ) THEN
            IF( A( K, J+1 ) /= (0.0E+0,0.0E+0) ) THEN
               ALPHA = (1.0E+0,0.0E+0) / A( K, J+1 )
               A(K,J+2:M) = WORK(3:1+M-J)
               A(K,J+2:M) = ALPHA*A(K,J+2:M)
            ELSE
               A(K,J+2:M) = (0.0E+0,0.0E+0)
            END IF
         END IF
      END IF
      J = J + 1
      GO TO 10
 20      CONTINUE
!
   ELSE
!
!        .....................................................
!        Factorize A as L*D*L**T using the lower triangle of A
!        .....................................................
!
 30      CONTINUE
      IF( J > MIN( M, NB ) ) GO TO 40
!
!        K is the column to be factorized
!         when being called from CHETRF_AA,
!         > for the first block column, J1 is 1, hence J1+J-1 is J,
!         > for the rest of the columns, J1 is 2, and J1+J-1 is J+1,
!
      K = J1+J-1
      IF( J == M ) THEN
!
!            Only need to compute T(J, J)
!
          MJ = 1
      ELSE
          MJ = M-J+1
      END IF
!
!        H(J:N, J) := A(J:N, J) - H(J:N, 1:(J-1)) * L(J, J1:(J-1))^T,
!         where H(J:N, J) has been initialized to be A(J:N, J)
!
      IF( K > 2 ) THEN
!
!        K is the column to be factorized
!         > for the first block column, K is J, skipping the first two
!           columns
!         > for the rest of the columns, K is J+1, skipping only the
!           first column
!
         CALL CLACGV( J-K1, A( J, 1 ), LDA )
         CALL CGEMV( 'No transpose', MJ, J-K1, &
                    -(1.0E+0,0.0E+0), H( J, K1 ), LDH, &
                          A( J, 1 ), LDA, &
                     (1.0E+0,0.0E+0), H( J, J ), 1 )
         CALL CLACGV( J-K1, A( J, 1 ), LDA )
      END IF
!
!        Copy H(J:N, J) into WORK
!
      WORK(1:MJ) = H(J:J+MJ-1,J)
!
      IF( J > K1 ) THEN
!
!           Compute WORK := WORK - L(J:N, J-1) * T(J-1,J),
!            where A(J-1, J) = T(J-1, J) and A(J, J-2) = L(J, J-1)
!
         ALPHA = -CONJG( A( J, K-1 ) )
         WORK(1:MJ) = WORK(1:MJ) + ALPHA*A(J:J+MJ-1,K-2)
      END IF
!
!        Set A(J, J) = T(J, J)
!
      A( J, K ) = REAL( WORK( 1 ) )
!
      IF( J < M ) THEN
!
!           Compute WORK(2:N) = T(J, J) L((J+1):N, J)
!            where A(J, J) = T(J, J) and A((J+1):N, J-1) = L((J+1):N, J)
!
         IF( K > 1 ) WORK(2:1+M-J) = WORK(2:1+M-J) -A(J,K)*A(J+1:M,K-1)
!
!           Find max(|WORK(2:n)|)
!
         I2 = ICAMAX( M-J, WORK( 2 ), 1 ) + 1
         PIV = WORK( I2 )
!
!           Apply hermitian pivot
!
         IF( (I2 /= 2) .AND. (PIV /= 0) ) THEN
!
!              Swap WORK(I1) and WORK(I2)
!
            I1 = 2
            WORK( I2 ) = WORK( I1 )
            WORK( I1 ) = PIV
!
!              Swap A(I1+1:N, I1) with A(I2, I1+1:N)
!
            I1 = I1+J-1
            I2 = I2+J-1
            A_TMP(1:I2-I1-1) = A(I1+1:I2-1,J1+I1-1)
            A(I1+1:I2-1,J1+I1-1) = A(I2,J1+I1:J1+I2-2)
            A(I2,J1+I1:J1+I2-2) = A_TMP(1:I2-I1-1)
            CALL CLACGV( I2-I1, A( I1+1, J1+I1-1 ), 1 )
            CALL CLACGV( I2-I1-1, A( I2, J1+I1 ), LDA )
!
!              Swap A(I2+1:N, I1) with A(I2+1:N, I2)
!
            IF( I2 < M ) THEN
               A_TMP(1:M-I2) = A(I2+1:M,J1+I1-1)
               A(I2+1:M,J1+I1-1) = A(I2+1:M,J1+I2-1)
               A(I2+1:M,J1+I2-1) = A_TMP(1:M-I2)
            ENDIF
!
!              Swap A(I1, I1) with A(I2, I2)
!
            PIV = A( I1, J1+I1-1 )
            A( I1, J1+I1-1 ) = A( I2, J1+I2-1 )
            A( I2, J1+I2-1 ) = PIV
!
!              Swap H(I1, I1:J1) with H(I2, I2:J1)
!
            H_TMP(1:I1-1) = H(I1,1:I1-1)
            H(I1,1:I1-1) = H(I2,1:I1-1)
            H(I2,1:I1-1) = H_TMP(1:I1-1)
            IPIV( I1 ) = I2
!
            IF( I1 > (K1-1) ) THEN
!
!                 Swap L(1:I1-1, I1) with L(1:I1-1, I2),
!                  skipping the first column
!
               A_TMP(1:I1-K1+1) = A(I1,1:1+I1-K1)
               A(I1,1:I1-K1+1) = A(I2,1:I1-K1+1)
               A(I2,1:I1-K1+1) = A_TMP(1:I1-K1+1)
            END IF
         ELSE
            IPIV( J+1 ) = J+1
         ENDIF
!
!           Set A(J+1, J) = T(J+1, J)
!
         A( J+1, K ) = WORK( 2 )
!
!
!              Copy A(J+1:N, J+1) into H(J+1:N, J),
!
         IF( J < NB ) H(J+1:M,J+1) = A(J+1:M,K+1)
!
!           Compute L(J+2, J+1) = WORK( 3:N ) / T(J, J+1),
!            where A(J, J+1) = T(J, J+1) and A(J+2:N, J) = L(J+2:N, J+1)
!
         IF( J < (M-1) ) THEN
            IF( A( J+1, K ) /= (0.0E+0,0.0E+0) ) THEN
               ALPHA = (1.0E+0,0.0E+0) / A( J+1, K )
               A(J+2:M,K) = WORK(3:1+M-J)
               A(J+2:M,K) = ALPHA*A(J+2:M,K)
            ELSE
               A(J+2:M,K) = (0.0E+0,0.0E+0)
            END IF
         END IF
      END IF
      J = J + 1
      GO TO 30
 40      CONTINUE
   END IF
   RETURN
!
!     End of CLAHEF_AA
!
END
