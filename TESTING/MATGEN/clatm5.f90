!> \brief \b CLATM5
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLATM5( PRTYPE, M, N, A, LDA, B, LDB, C, LDC, D, LDD,
!                          E, LDE, F, LDF, R, LDR, L, LDL, ALPHA, QBLCKA,
!                          QBLCKB )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LDC, LDD, LDE, LDF, LDL, LDR, M, N,
!      $                   PRTYPE, QBLCKA, QBLCKB
!       REAL               ALPHA
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * ),
!      $                   D( LDD, * ), E( LDE, * ), F( LDF, * ),
!      $                   L( LDL, * ), R( LDR, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLATM5 generates matrices involved in the Generalized Sylvester
!> equation:
!>
!>     A * R - L * B = C
!>     D * R - L * E = F
!>
!> They also satisfy (the diagonalization condition)
!>
!>  [ I -L ] ( [ A  -C ], [ D -F ] ) [ I  R ] = ( [ A    ], [ D    ] )
!>  [    I ] ( [     B ]  [    E ] ) [    I ]   ( [    B ]  [    E ] )
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] PRTYPE
!> \verbatim
!>          PRTYPE is INTEGER
!>          "Points" to a certain type of the matrices to generate
!>          (see further details).
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          Specifies the order of A and D and the number of rows in
!>          C, F,  R and L.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          Specifies the order of B and E and the number of columns in
!>          C, F, R and L.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA, M).
!>          On exit A M-by-M is initialized according to PRTYPE.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB, N).
!>          On exit B N-by-N is initialized according to PRTYPE.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of B.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC, N).
!>          On exit C M-by-N is initialized according to PRTYPE.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of C.
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is COMPLEX array, dimension (LDD, M).
!>          On exit D M-by-M is initialized according to PRTYPE.
!> \endverbatim
!>
!> \param[in] LDD
!> \verbatim
!>          LDD is INTEGER
!>          The leading dimension of D.
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is COMPLEX array, dimension (LDE, N).
!>          On exit E N-by-N is initialized according to PRTYPE.
!> \endverbatim
!>
!> \param[in] LDE
!> \verbatim
!>          LDE is INTEGER
!>          The leading dimension of E.
!> \endverbatim
!>
!> \param[out] F
!> \verbatim
!>          F is COMPLEX array, dimension (LDF, N).
!>          On exit F M-by-N is initialized according to PRTYPE.
!> \endverbatim
!>
!> \param[in] LDF
!> \verbatim
!>          LDF is INTEGER
!>          The leading dimension of F.
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is COMPLEX array, dimension (LDR, N).
!>          On exit R M-by-N is initialized according to PRTYPE.
!> \endverbatim
!>
!> \param[in] LDR
!> \verbatim
!>          LDR is INTEGER
!>          The leading dimension of R.
!> \endverbatim
!>
!> \param[out] L
!> \verbatim
!>          L is COMPLEX array, dimension (LDL, N).
!>          On exit L M-by-N is initialized according to PRTYPE.
!> \endverbatim
!>
!> \param[in] LDL
!> \verbatim
!>          LDL is INTEGER
!>          The leading dimension of L.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is REAL
!>          Parameter used in generating PRTYPE = 1 and 5 matrices.
!> \endverbatim
!>
!> \param[in] QBLCKA
!> \verbatim
!>          QBLCKA is INTEGER
!>          When PRTYPE = 3, specifies the distance between 2-by-2
!>          blocks on the diagonal in A. Otherwise, QBLCKA is not
!>          referenced. QBLCKA > 1.
!> \endverbatim
!>
!> \param[in] QBLCKB
!> \verbatim
!>          QBLCKB is INTEGER
!>          When PRTYPE = 3, specifies the distance between 2-by-2
!>          blocks on the diagonal in B. Otherwise, QBLCKB is not
!>          referenced. QBLCKB > 1.
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
!> \ingroup complex_matgen
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  PRTYPE = 1: A and B are Jordan blocks, D and E are identity matrices
!>
!>             A : if (i == j) then A(i, j) = 1.0
!>                 if (j == i + 1) then A(i, j) = -1.0
!>                 else A(i, j) = 0.0,            i, j = 1...M
!>
!>             B : if (i == j) then B(i, j) = 1.0 - ALPHA
!>                 if (j == i + 1) then B(i, j) = 1.0
!>                 else B(i, j) = 0.0,            i, j = 1...N
!>
!>             D : if (i == j) then D(i, j) = 1.0
!>                 else D(i, j) = 0.0,            i, j = 1...M
!>
!>             E : if (i == j) then E(i, j) = 1.0
!>                 else E(i, j) = 0.0,            i, j = 1...N
!>
!>             L =  R are chosen from [-10...10],
!>                  which specifies the right hand sides (C, F).
!>
!>  PRTYPE = 2 or 3: Triangular and/or quasi- triangular.
!>
!>             A : if (i <= j) then A(i, j) = [-1...1]
!>                 else A(i, j) = 0.0,             i, j = 1...M
!>
!>                 if (PRTYPE = 3) then
!>                    A(k + 1, k + 1) = A(k, k)
!>                    A(k + 1, k) = [-1...1]
!>                    sign(A(k, k + 1) = -(sin(A(k + 1, k))
!>                        k = 1, M - 1, QBLCKA
!>
!>             B : if (i <= j) then B(i, j) = [-1...1]
!>                 else B(i, j) = 0.0,            i, j = 1...N
!>
!>                 if (PRTYPE = 3) then
!>                    B(k + 1, k + 1) = B(k, k)
!>                    B(k + 1, k) = [-1...1]
!>                    sign(B(k, k + 1) = -(sign(B(k + 1, k))
!>                        k = 1, N - 1, QBLCKB
!>
!>             D : if (i <= j) then D(i, j) = [-1...1].
!>                 else D(i, j) = 0.0,            i, j = 1...M
!>
!>
!>             E : if (i <= j) then D(i, j) = [-1...1]
!>                 else E(i, j) = 0.0,            i, j = 1...N
!>
!>                 L, R are chosen from [-10...10],
!>                 which specifies the right hand sides (C, F).
!>
!>  PRTYPE = 4 Full
!>             A(i, j) = [-10...10]
!>             D(i, j) = [-1...1]    i,j = 1...M
!>             B(i, j) = [-10...10]
!>             E(i, j) = [-1...1]    i,j = 1...N
!>             R(i, j) = [-10...10]
!>             L(i, j) = [-1...1]    i = 1..M ,j = 1...N
!>
!>             L, R specifies the right hand sides (C, F).
!>
!>  PRTYPE = 5 special case common and/or close eigs.
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE CLATM5( PRTYPE, M, N, A, LDA, B, LDB, C, LDC, D, LDD, &
                      E, LDE, F, LDF, R, LDR, L, LDL, ALPHA, QBLCKA, &
                      QBLCKB )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, LDB, LDC, LDD, LDE, LDF, LDL, LDR, M, N, &
                      PRTYPE, QBLCKA, QBLCKB
   REAL               ALPHA
!     ..
!     .. Array Arguments ..
   COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * ), &
                      D( LDD, * ), E( LDE, * ), F( LDF, * ), &
                      L( LDL, * ), R( LDR, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   COMPLEX            ONE, TWO, ZERO, HALF, TWENTY
   PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ), &
                      TWO = ( 2.0E+0, 0.0E+0 ), &
                      ZERO = ( 0.0E+0, 0.0E+0 ), &
                      HALF = ( 0.5E+0, 0.0E+0 ), &
                      TWENTY = ( 2.0E+1, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
   INTEGER            I, J, K
   COMPLEX            IMEPS, REEPS
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          CMPLX, MOD, SIN
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMM
!     ..
!     .. Executable Statements ..
!
   IF( PRTYPE == 1 ) THEN
      DO I = 1, M
         DO J = 1, M
            IF( I == J ) THEN
               A( I, J ) = ONE
               D( I, J ) = ONE
            ELSE IF( I == J-1 ) THEN
               A( I, J ) = -ONE
               D( I, J ) = ZERO
            ELSE
               A( I, J ) = ZERO
               D( I, J ) = ZERO
            END IF
         ENDDO
      ENDDO
!
      DO I = 1, N
         DO J = 1, N
            IF( I == J ) THEN
               B( I, J ) = ONE - ALPHA
               E( I, J ) = ONE
            ELSE IF( I == J-1 ) THEN
               B( I, J ) = ONE
               E( I, J ) = ZERO
            ELSE
               B( I, J ) = ZERO
               E( I, J ) = ZERO
            END IF
         ENDDO
      ENDDO
!
      DO I = 1, M
         DO J = 1, N
            R( I, J ) = ( HALF-SIN( CMPLX( I / J ) ) )*TWENTY
            L( I, J ) = R( I, J )
         ENDDO
      ENDDO
!
   ELSE IF( PRTYPE == 2 .OR. PRTYPE == 3 ) THEN
      DO I = 1, M
         DO J = 1, M
            IF( I <= J ) THEN
               A( I, J ) = ( HALF-SIN( CMPLX( I ) ) )*TWO
               D( I, J ) = ( HALF-SIN( CMPLX( I*J ) ) )*TWO
            ELSE
               A( I, J ) = ZERO
               D( I, J ) = ZERO
            END IF
         ENDDO
      ENDDO
!
      DO I = 1, N
         DO J = 1, N
            IF( I <= J ) THEN
               B( I, J ) = ( HALF-SIN( CMPLX( I+J ) ) )*TWO
               E( I, J ) = ( HALF-SIN( CMPLX( J ) ) )*TWO
            ELSE
               B( I, J ) = ZERO
               E( I, J ) = ZERO
            END IF
         ENDDO
         ENDDO
!
      DO I = 1, M
         DO J = 1, N
            R( I, J ) = ( HALF-SIN( CMPLX( I*J ) ) )*TWENTY
            L( I, J ) = ( HALF-SIN( CMPLX( I+J ) ) )*TWENTY
            ENDDO
         ENDDO
!
      IF( PRTYPE == 3 ) THEN
         IF( QBLCKA <= 1 ) &
            QBLCKA = 2
         DO K = 1, M - 1, QBLCKA
            A( K+1, K+1 ) = A( K, K )
            A( K+1, K ) = -SIN( A( K, K+1 ) )
            ENDDO
!
         IF( QBLCKB <= 1 ) &
            QBLCKB = 2
         DO K = 1, N - 1, QBLCKB
            B( K+1, K+1 ) = B( K, K )
            B( K+1, K ) = -SIN( B( K, K+1 ) )
            ENDDO
      END IF
!
   ELSE IF( PRTYPE == 4 ) THEN
      DO I = 1, M
         DO J = 1, M
            A( I, J ) = ( HALF-SIN( CMPLX( I*J ) ) )*TWENTY
            D( I, J ) = ( HALF-SIN( CMPLX( I+J ) ) )*TWO
            ENDDO
         ENDDO
!
      DO I = 1, N
         DO J = 1, N
            B( I, J ) = ( HALF-SIN( CMPLX( I+J ) ) )*TWENTY
            E( I, J ) = ( HALF-SIN( CMPLX( I*J ) ) )*TWO
            ENDDO
         ENDDO
!
      DO I = 1, M
         DO J = 1, N
            R( I, J ) = ( HALF-SIN( CMPLX( J / I ) ) )*TWENTY
            L( I, J ) = ( HALF-SIN( CMPLX( I*J ) ) )*TWO
            ENDDO
         ENDDO
!
   ELSE IF( PRTYPE >= 5 ) THEN
      REEPS = HALF*TWO*TWENTY / ALPHA
      IMEPS = ( HALF-TWO ) / ALPHA
      DO I = 1, M
         DO J = 1, N
            R( I, J ) = ( HALF-SIN( CMPLX( I*J ) ) )*ALPHA / TWENTY
            L( I, J ) = ( HALF-SIN( CMPLX( I+J ) ) )*ALPHA / TWENTY
            ENDDO
         ENDDO
!
      DO I = 1, M
         D( I, I ) = ONE
         ENDDO
!
      DO I = 1, M
         IF( I <= 4 ) THEN
            A( I, I ) = ONE
            IF( I > 2 ) &
               A( I, I ) = ONE + REEPS
            IF( MOD( I, 2 ) /= 0 .AND. I < M ) THEN
               A( I, I+1 ) = IMEPS
            ELSE IF( I > 1 ) THEN
               A( I, I-1 ) = -IMEPS
            END IF
         ELSE IF( I <= 8 ) THEN
            IF( I <= 6 ) THEN
               A( I, I ) = REEPS
            ELSE
               A( I, I ) = -REEPS
            END IF
            IF( MOD( I, 2 ) /= 0 .AND. I < M ) THEN
               A( I, I+1 ) = ONE
            ELSE IF( I > 1 ) THEN
               A( I, I-1 ) = -ONE
            END IF
         ELSE
            A( I, I ) = ONE
            IF( MOD( I, 2 ) /= 0 .AND. I < M ) THEN
               A( I, I+1 ) = IMEPS*2
            ELSE IF( I > 1 ) THEN
               A( I, I-1 ) = -IMEPS*2
            END IF
         END IF
         ENDDO
!
      DO I = 1, N
         E( I, I ) = ONE
         IF( I <= 4 ) THEN
            B( I, I ) = -ONE
            IF( I > 2 ) &
               B( I, I ) = ONE - REEPS
            IF( MOD( I, 2 ) /= 0 .AND. I < N ) THEN
               B( I, I+1 ) = IMEPS
            ELSE IF( I > 1 ) THEN
               B( I, I-1 ) = -IMEPS
            END IF
         ELSE IF( I <= 8 ) THEN
            IF( I <= 6 ) THEN
               B( I, I ) = REEPS
            ELSE
               B( I, I ) = -REEPS
            END IF
            IF( MOD( I, 2 ) /= 0 .AND. I < N ) THEN
               B( I, I+1 ) = ONE + IMEPS
            ELSE IF( I > 1 ) THEN
               B( I, I-1 ) = -ONE - IMEPS
            END IF
         ELSE
            B( I, I ) = ONE - REEPS
            IF( MOD( I, 2 ) /= 0 .AND. I < N ) THEN
               B( I, I+1 ) = IMEPS*2
            ELSE IF( I > 1 ) THEN
               B( I, I-1 ) = -IMEPS*2
            END IF
         END IF
         ENDDO
   END IF
!
!     Compute rhs (C, F)
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CGEMM( 'N', 'N', M, N, M, ONE, A, LDA, R, LDR, ZERO, C, LDC )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CGEMM( 'N', 'N', M, N, N, -ONE, L, LDL, B, LDB, ONE, C, LDC )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CGEMM( 'N', 'N', M, N, M, ONE, D, LDD, R, LDR, ZERO, F, LDF )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CGEMM( 'N', 'N', M, N, N, -ONE, L, LDL, E, LDE, ONE, F, LDF )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     End of CLATM5
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        


