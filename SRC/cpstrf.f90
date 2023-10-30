!> \brief \b CPSTRF computes the Cholesky factorization with complete pivoting of complex Hermitian positive semidefinite matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CPSTRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpstrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpstrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpstrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CPSTRF( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )
!
!       .. Scalar Arguments ..
!       REAL               TOL
!       INTEGER            INFO, LDA, N, RANK
!       CHARACTER          UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * )
!       REAL               WORK( 2*N )
!       INTEGER            PIV( N )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CPSTRF computes the Cholesky factorization with complete
!> pivoting of a complex Hermitian positive semidefinite matrix A.
!>
!> The factorization has the form
!>    P**T * A * P = U**H * U ,  if UPLO = 'U',
!>    P**T * A * P = L  * L**H,  if UPLO = 'L',
!> where U is an upper triangular matrix and L is lower triangular, and
!> P is stored as vector PIV.
!>
!> This algorithm does not attempt to check that A is positive
!> semidefinite. This version of the algorithm calls level 3 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix A is stored.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!>          n by n upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading n by n lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, if INFO = 0, the factor U or L from the Cholesky
!>          factorization as above.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] PIV
!> \verbatim
!>          PIV is INTEGER array, dimension (N)
!>          PIV is such that the nonzero entries are P( PIV(K), K ) = 1.
!> \endverbatim
!>
!> \param[out] RANK
!> \verbatim
!>          RANK is INTEGER
!>          The rank of A given by the number of steps the algorithm
!>          completed.
!> \endverbatim
!>
!> \param[in] TOL
!> \verbatim
!>          TOL is REAL
!>          User defined tolerance. If TOL < 0, then N*U*MAX( A(K,K) )
!>          will be used. The algorithm terminates at the (K-1)st step
!>          if the pivot <= TOL.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (2*N)
!>          Work space.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          < 0: If INFO = -K, the K-th argument had an illegal value,
!>          = 0: algorithm completed successfully, and
!>          > 0: the matrix A is either rank deficient with computed rank
!>               as returned in RANK, or is not positive semidefinite. See
!>               Section 7 of LAPACK Working Note #161 for further
!>               information.
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
!> \ingroup pstrf
!
!  =====================================================================
   SUBROUTINE CPSTRF( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   REAL               TOL
   INTEGER            INFO, LDA, N, RANK
   CHARACTER          UPLO
!     ..
!     .. Array Arguments ..
   COMPLEX            A( LDA, * )
   REAL               WORK( 2*N )
   INTEGER            PIV( N )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   COMPLEX            CTEMP, A_TMP( LDA ), A_TMP2( N )
   REAL               AJJ, SSTOP, STEMP
   INTEGER            I, ITEMP, J, JB, K, NB, PVT
   LOGICAL            UPPER
!     ..
!     .. External Functions ..
   REAL               SLAMCH
   INTEGER            ILAENV
   LOGICAL            LSAME, SISNAN
   EXTERNAL           SLAMCH, ILAENV, LSAME, SISNAN
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMV, CHERK, CLACGV, CPSTF2, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   UPPER = LSAME( UPLO, 'U' )
   IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -4
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CPSTRF', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) RETURN
!
!     Get block size
!
   NB = ILAENV( 1, 'CPOTRF', UPLO, N, -1, -1, -1 )
   IF( NB <= 1 .OR. NB >= N ) THEN
!
!        Use unblocked code
!
      CALL CPSTF2( UPLO, N, A( 1, 1 ), LDA, PIV, RANK, TOL, WORK, INFO )
      GO TO 230
!
   ELSE
!
!     Initialize PIV
!
      DO I = 1, N
         PIV( I ) = I
      ENDDO
!
!     Compute stopping value
!
      DO I = 1, N
         WORK( I ) = REAL( A( I, I ) )
      ENDDO
      PVT = MAXLOC( WORK( 1:N ), 1 )
      AJJ = REAL( A( PVT, PVT ) )
      IF( AJJ <= 0.0E+0.OR.SISNAN( AJJ ) ) THEN
         RANK = 0
         INFO = 1
         GO TO 230
      END IF
!
!     Compute stopping value if not supplied
!
      IF( TOL < 0.0E+0 ) THEN
         SSTOP = N * SLAMCH( 'Epsilon' ) * AJJ
      ELSE
         SSTOP = TOL
      END IF
!
!
      IF( UPPER ) THEN
!
!           Compute the Cholesky factorization P**T * A * P = U**H * U
!
         DO K = 1, N, NB
!
!              Account for last block not being NB wide
!
            JB = MIN( NB, N-K+1 )
!
!              Set relevant part of first half of WORK to zero,
!              holds dot products
!
            WORK(K:N) = 0
!
            DO J = K, K + JB - 1
!
!              Find pivot, test for exit, else swap rows and columns
!              Update dot products, compute possible pivots which are
!              stored in the second half of WORK
!
               DO I = J, N
!
                  IF( J > K ) THEN
                     WORK( I ) = WORK( I ) + &
                                 REAL( CONJG( A( J-1, I ) )* &
                                       A( J-1, I ) )
                  END IF
                  WORK( N+I ) = REAL( A( I, I ) ) - WORK( I )
!
               ENDDO
!
               IF( J > 1 ) THEN
                  ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 )
                  PVT = ITEMP + J - 1
                  AJJ = WORK( N+PVT )
                  IF( AJJ <= SSTOP.OR.SISNAN( AJJ ) ) THEN
                     A( J, J ) = AJJ
                     GO TO 220
                  END IF
               END IF
!
               IF( J /= PVT ) THEN
!
!                    Pivot OK, so can now swap pivot rows and columns
!
                  A( PVT, PVT ) = A( J, J )
                  A_TMP(1:J-1) = A(1:J-1,J)
                  A(1:J-1,J) = A(1:J-1,PVT)
                  A(1:J-1,PVT) = A_TMP(1:J-1)
                  IF( PVT < N ) THEN
                     A_TMP2(1:N-PVT) = A(J,PVT+1:N)
                     A(J,PVT+1:N) = A(PVT,PVT+1:N)
                     A(PVT,PVT+1:N) = A_TMP2(1:N-PVT)
                  ENDIF
                  A_TMP2(J+1:PVT-1) = CONJG(A(J,J+1:PVT-1))
                  A(J,J+1:PVT-1) = CONJG(A(J+1:PVT-1,PVT))
                  A(J+1:PVT-1,PVT) = A_TMP2(J+1:PVT-1)
                  A( J, PVT ) = CONJG( A( J, PVT ) )
!
!                    Swap dot products and PIV
!
                  STEMP = WORK( J )
                  WORK( J ) = WORK( PVT )
                  WORK( PVT ) = STEMP
                  ITEMP = PIV( PVT )
                  PIV( PVT ) = PIV( J )
                  PIV( J ) = ITEMP
               END IF
!
               AJJ = SQRT( AJJ )
               A( J, J ) = AJJ
!
!                 Compute elements J+1:N of row J.
!
               IF( J < N ) THEN
                  A(1:J-1,J) = CONJG(A(1:J-1,J))
                  CALL CGEMV( 'Trans', J-K, N-J, -(1.0E+0,0.0E+0), A( K, J+1 ), &
                              LDA, A( K, J ), 1, (1.0E+0,0.0E+0), A( J, J+1 ), &
                              LDA )
                  A(1:J-1,J) = CONJG(A(1:J-1,J))
                  A(J,J+1:N) = A(J,J+1:N) / AJJ
               END IF
!
            ENDDO
!
!              Update trailing matrix, J already incremented
!
            IF( K+JB <= N ) THEN
               CALL CHERK( 'Upper', 'Conj Trans', N-J+1, JB, -1.0E+0, &
                           A( K, J ), LDA, 1.0E+0, A( J, J ), LDA )
            END IF
!
         ENDDO
!
      ELSE
!
!        Compute the Cholesky factorization P**T * A * P = L * L**H
!
         DO K = 1, N, NB
!
!              Account for last block not being NB wide
!
            JB = MIN( NB, N-K+1 )
!
!              Set relevant part of first half of WORK to zero,
!              holds dot products
!
            WORK(K:N) = 0
!
            DO J = K, K + JB - 1
!
!              Find pivot, test for exit, else swap rows and columns
!              Update dot products, compute possible pivots which are
!              stored in the second half of WORK
!
               DO I = J, N
!
                  IF( J > K ) THEN
                     WORK( I ) = WORK( I ) + &
                                 REAL( CONJG( A( I, J-1 ) )* &
                                       A( I, J-1 ) )
                  END IF
                  WORK( N+I ) = REAL( A( I, I ) ) - WORK( I )
!
               ENDDO
!
               IF( J > 1 ) THEN
                  ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 )
                  PVT = ITEMP + J - 1
                  AJJ = WORK( N+PVT )
                  IF( AJJ <= SSTOP.OR.SISNAN( AJJ ) ) THEN
                     A( J, J ) = AJJ
                     GO TO 220
                  END IF
               END IF
!
               IF( J /= PVT ) THEN
!
!                    Pivot OK, so can now swap pivot rows and columns
!
                  A( PVT, PVT ) = A( J, J )
                  A_TMP(1:J-1) = A(J,1:J-1)
                  A(J,1:J-1) = A(PVT,1:J-1)
                  A(PVT,1:J-1) = A_TMP(1:J-1)
                  IF( PVT < N ) THEN
                     A_TMP(1:N-PVT) = A(PVT+1:N,J)
                     A(PVT+1:N,J) = A(PVT+1:N,PVT)
                     A(PVT+1:N,PVT) = A_TMP(1:N-PVT)
                  ENDIF
                  A_TMP(J+1:PVT-1) = CONJG(A(J+1:PVT-1,J))
                  A(J+1:PVT-1,J) = CONJG(A(PVT,J+1:PVT-1))
                  A(PVT,J+1:PVT-1) = A_TMP(J+1:PVT-1) 
                  A( PVT, J ) = CONJG( A( PVT, J ) )
!
!                    Swap dot products and PIV
!
                  STEMP = WORK( J )
                  WORK( J ) = WORK( PVT )
                  WORK( PVT ) = STEMP
                  ITEMP = PIV( PVT )
                  PIV( PVT ) = PIV( J )
                  PIV( J ) = ITEMP
               END IF
!
               AJJ = SQRT( AJJ )
               A( J, J ) = AJJ
!
!                 Compute elements J+1:N of column J.
!
               IF( J < N ) THEN
                  A(J,1:J-1) = CONJG(A(J,1:J-1))
                  CALL CGEMV( 'No Trans', N-J, J-K, -(1.0E+0,0.0E+0), &
                              A( J+1, K ), LDA, A( J, K ), LDA, (1.0E+0,0.0E+0), &
                              A( J+1, J ), 1 )
                  A(J,1:J-1) = CONJG(A(J,1:J-1))
                  A(J+1:N,J) = A(J+1:N,J) / AJJ
               END IF
!
            ENDDO
!
!              Update trailing matrix, J already incremented
!
            IF( K+JB <= N ) THEN
               CALL CHERK( 'Lower', 'No Trans', N-J+1, JB, -1.0E+0, &
                           A( J, K ), LDA, 1.0E+0, A( J, J ), LDA )
            END IF
!
         ENDDO
!
      END IF
   END IF
!
!     Ran to completion, A has full rank
!
   RANK = N
!
   GO TO 230
  220 CONTINUE
!
!     Rank is the number of steps completed.  Set INFO = 1 to signal
!     that the factorization cannot be used to solve a system.
!
   RANK = J - 1
   INFO = 1
!
  230 CONTINUE
   RETURN
!
!     End of CPSTRF
!
END
