!> \brief \b DLATTP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLATTP( IMAT, UPLO, TRANS, DIAG, ISEED, N, A, B, WORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, TRANS, UPLO
!       INTEGER            IMAT, INFO, N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       DOUBLE PRECISION   A( * ), B( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLATTP generates a triangular test matrix in packed storage.
!> IMAT and UPLO uniquely specify the properties of the test
!> matrix, which is returned in the array AP.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] IMAT
!> \verbatim
!>          IMAT is INTEGER
!>          An integer key describing which matrix to generate for this
!>          path.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A will be upper or lower
!>          triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies whether the matrix or its transpose will be used.
!>          = 'N':  No transpose
!>          = 'T':  Transpose
!>          = 'C':  Conjugate transpose (= Transpose)
!> \endverbatim
!>
!> \param[out] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A is unit triangular.
!>          = 'N':  Non-unit triangular
!>          = 'U':  Unit triangular
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          The seed vector for the random number generator (used in
!>          DLATMS).  Modified on exit.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix to be generated.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          The upper or lower triangular matrix A, packed columnwise in
!>          a linear array.  The j-th column of A is stored in the array
!>          AP as follows:
!>          if UPLO = 'U', AP((j-1)*j/2 + i) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L',
!>             AP((j-1)*(n-j) + j*(j+1)/2 + i-j) = A(i,j) for j<=i<=n.
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (N)
!>          The right hand side vector, if IMAT > 10.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (3*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0: if INFO = -k, the k-th argument had an illegal value
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
!> \ingroup double_lin
!
!  =====================================================================
   SUBROUTINE DLATTP( IMAT, UPLO, TRANS, DIAG, ISEED, N, A, B, WORK, &
                      INFO )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          DIAG, TRANS, UPLO
   INTEGER            IMAT, INFO, N
!     ..
!     .. Array Arguments ..
   INTEGER            ISEED( 4 )
   DOUBLE PRECISION   A( * ), B( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ONE, TWO, ZERO
   PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            UPPER
   CHARACTER          DIST, PACKIT, TYPE
   CHARACTER*3        PATH
   INTEGER            I, IY, J, JC, JCNEXT, JCOUNT, JJ, JL, JR, JX, &
                      KL, KU, MODE
   DOUBLE PRECISION   ANORM, BIGNUM, BNORM, BSCAL, C, CNDNUM, PLUS1, &
                      PLUS2, RA, RB, REXP, S, SFAC, SMLNUM, STAR1, &
                      STEMP, T, TEXP, TLEFT, TSCAL, ULP, UNFL, X, Y, &
                      Z
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            IDAMAX
   DOUBLE PRECISION   DLAMCH, DLARND
   EXTERNAL           LSAME, IDAMAX, DLAMCH, DLARND
!     ..
!     .. External Subroutines ..
   EXTERNAL           DLARNV, DLATB4, DLATMS, DROT, DROTG, DSCAL
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, DBLE, MAX, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
   PATH( 1: 1 ) = 'Double precision'
   PATH( 2: 3 ) = 'TP'
   UNFL = DLAMCH( 'Safe minimum' )
   ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
   SMLNUM = UNFL
   BIGNUM = ( ONE-ULP ) / SMLNUM
   IF( ( IMAT >= 7 .AND. IMAT <= 10 ) .OR. IMAT == 18 ) THEN
      DIAG = 'U'
   ELSE
      DIAG = 'N'
   END IF
   INFO = 0
!
!     Quick return if N <= 0.
!
   IF( N <= 0 ) &
      RETURN
!
!     Call DLATB4 to set parameters for DLATMS.
!
   UPPER = LSAME( UPLO, 'U' )
   IF( UPPER ) THEN
      CALL DLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, &
                   CNDNUM, DIST )
      PACKIT = 'C'
   ELSE
      CALL DLATB4( PATH, -IMAT, N, N, TYPE, KL, KU, ANORM, MODE, &
                   CNDNUM, DIST )
      PACKIT = 'R'
   END IF
!
!     IMAT <= 6:  Non-unit triangular matrix
!
   IF( IMAT <= 6 ) THEN
      CALL DLATMS( N, N, DIST, ISEED, TYPE, B, MODE, CNDNUM, ANORM, &
                   KL, KU, PACKIT, A, N, WORK, INFO )
!
!     IMAT > 6:  Unit triangular matrix
!     The diagonal is deliberately set to something other than 1.
!
!     IMAT = 7:  Matrix is the identity
!
   ELSE IF( IMAT == 7 ) THEN
      IF( UPPER ) THEN
         JC = 1
         DO J = 1, N
            DO I = 1, J - 1
               A( JC+I-1 ) = ZERO
            ENDDO
            A( JC+J-1 ) = J
            JC = JC + J
         ENDDO
      ELSE
         JC = 1
         DO J = 1, N
            A( JC ) = J
            DO I = J + 1, N
               A( JC+I-J ) = ZERO
            ENDDO
            JC = JC + N - J + 1
         ENDDO
      END IF
!
!     IMAT > 7:  Non-trivial unit triangular matrix
!
!     Generate a unit triangular matrix T with condition CNDNUM by
!     forming a triangular matrix with known singular values and
!     filling in the zero entries with Givens rotations.
!
   ELSE IF( IMAT <= 10 ) THEN
      IF( UPPER ) THEN
         JC = 0
         DO J = 1, N
            DO I = 1, J - 1
               A( JC+I ) = ZERO
            ENDDO
            A( JC+J ) = J
            JC = JC + J
         ENDDO
      ELSE
         JC = 1
         DO J = 1, N
            A( JC ) = J
            DO I = J + 1, N
               A( JC+I-J ) = ZERO
            ENDDO
            JC = JC + N - J + 1
         ENDDO
      END IF
!
!        Since the trace of a unit triangular matrix is 1, the product
!        of its singular values must be 1.  Let s = sqrt(CNDNUM),
!        x = sqrt(s) - 1/sqrt(s), y = sqrt(2/(n-2))*x, and z = x**2.
!        The following triangular matrix has singular values s, 1, 1,
!        ..., 1, 1/s:
!
!        1  y  y  y  ...  y  y  z
!           1  0  0  ...  0  0  y
!              1  0  ...  0  0  y
!                 .  ...  .  .  .
!                     .   .  .  .
!                         1  0  y
!                            1  y
!                               1
!
!        To fill in the zeros, we first multiply by a matrix with small
!        condition number of the form
!
!        1  0  0  0  0  ...
!           1  +  *  0  0  ...
!              1  +  0  0  0
!                 1  +  *  0  0
!                    1  +  0  0
!                       ...
!                          1  +  0
!                             1  0
!                                1
!
!        Each element marked with a '*' is formed by taking the product
!        of the adjacent elements marked with '+'.  The '*'s can be
!        chosen freely, and the '+'s are chosen so that the inverse of
!        T will have elements of the same magnitude as T.  If the *'s in
!        both T and inv(T) have small magnitude, T is well conditioned.
!        The two offdiagonals of T are stored in WORK.
!
!        The product of these two matrices has the form
!
!        1  y  y  y  y  y  .  y  y  z
!           1  +  *  0  0  .  0  0  y
!              1  +  0  0  .  0  0  y
!                 1  +  *  .  .  .  .
!                    1  +  .  .  .  .
!                       .  .  .  .  .
!                          .  .  .  .
!                             1  +  y
!                                1  y
!                                   1
!
!        Now we multiply by Givens rotations, using the fact that
!
!              [  c   s ] [  1   w ] [ -c  -s ] =  [  1  -w ]
!              [ -s   c ] [  0   1 ] [  s  -c ]    [  0   1 ]
!        and
!              [ -c  -s ] [  1   0 ] [  c   s ] =  [  1   0 ]
!              [  s  -c ] [  w   1 ] [ -s   c ]    [ -w   1 ]
!
!        where c = w / sqrt(w**2+4) and s = 2 / sqrt(w**2+4).
!
      STAR1 = 0.25D0
      SFAC = 0.5D0
      PLUS1 = SFAC
      DO J = 1, N, 2
         PLUS2 = STAR1 / PLUS1
         WORK( J ) = PLUS1
         WORK( N+J ) = STAR1
         IF( J+1 <= N ) THEN
            WORK( J+1 ) = PLUS2
            WORK( N+J+1 ) = ZERO
            PLUS1 = STAR1 / PLUS2
            REXP = DLARND( 2, ISEED )
            STAR1 = STAR1*( SFAC**REXP )
            IF( REXP < ZERO ) THEN
               STAR1 = -SFAC**( ONE-REXP )
            ELSE
               STAR1 = SFAC**( ONE+REXP )
            END IF
         END IF
      ENDDO
!
      X = SQRT( CNDNUM ) - ONE / SQRT( CNDNUM )
      IF( N > 2 ) THEN
         Y = SQRT( TWO / DBLE( N-2 ) )*X
      ELSE
         Y = ZERO
      END IF
      Z = X*X
!
      IF( UPPER ) THEN
!
!           Set the upper triangle of A with a unit triangular matrix
!           of known condition number.
!
         JC = 1
         DO J = 2, N
            A( JC+1 ) = Y
            IF( J > 2 ) &
               A( JC+J-1 ) = WORK( J-2 )
            IF( J > 3 ) &
               A( JC+J-2 ) = WORK( N+J-3 )
            JC = JC + J
            ENDDO
         JC = JC - N
         A( JC+1 ) = Z
         DO J = 2, N - 1
            A( JC+J ) = Y
            ENDDO
      ELSE
!
!           Set the lower triangle of A with a unit triangular matrix
!           of known condition number.
!
         DO I = 2, N - 1
            A( I ) = Y
            ENDDO
         A( N ) = Z
         JC = N + 1
         DO J = 2, N - 1
            A( JC+1 ) = WORK( J-1 )
            IF( J < N-1 ) &
               A( JC+2 ) = WORK( N+J-1 )
            A( JC+N-J ) = Y
            JC = JC + N - J + 1
            ENDDO
      END IF
!
!        Fill in the zeros using Givens rotations
!
      IF( UPPER ) THEN
         JC = 1
         DO J = 1, N - 1
            JCNEXT = JC + J
            RA = A( JCNEXT+J-1 )
            RB = TWO
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DROTG( RA, RB, C, S )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DROTG : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
!              Multiply by [ c  s; -s  c] on the left.
!
            IF( N > J+1 ) THEN
               JX = JCNEXT + J
               DO I = J + 2, N
                  STEMP = C*A( JX+J ) + S*A( JX+J+1 )
                  A( JX+J+1 ) = -S*A( JX+J ) + C*A( JX+J+1 )
                  A( JX+J ) = STEMP
                  JX = JX + I
                  ENDDO
            END IF
!
!              Multiply by [-c -s;  s -c] on the right.
!
            IF( J > 1 )  THEN
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DROT( J-1, A( JCNEXT ), 1, A( JC ), 1, -C, -S )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DROT : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
            ENDIF
!
!              Negate A(J,J+1).
!
            A( JCNEXT+J-1 ) = -A( JCNEXT+J-1 )
            JC = JCNEXT
            ENDDO
      ELSE
         JC = 1
         DO J = 1, N - 1
            JCNEXT = JC + N - J + 1
            RA = A( JC+1 )
            RB = TWO
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DROTG( RA, RB, C, S )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DROTG : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
!              Multiply by [ c -s;  s  c] on the right.
!
            IF( N > J+1 )  THEN
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DROT( N-J-1, A( JCNEXT+1 ), 1, A( JC+2 ), 1, C, &
                          -S )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DROT : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
            ENDIF
!
!              Multiply by [-c  s; -s -c] on the left.
!
            IF( J > 1 ) THEN
               JX = 1
               DO I = 1, J - 1
                  STEMP = -C*A( JX+J-I ) + S*A( JX+J-I+1 )
                  A( JX+J-I+1 ) = -S*A( JX+J-I ) - C*A( JX+J-I+1 )
                  A( JX+J-I ) = STEMP
                  JX = JX + N - I + 1
                  ENDDO
            END IF
!
!              Negate A(J+1,J).
!
            A( JC+1 ) = -A( JC+1 )
            JC = JCNEXT
            ENDDO
      END IF
!
!     IMAT > 10:  Pathological test cases.  These triangular matrices
!     are badly scaled or badly conditioned, so when used in solving a
!     triangular system they may cause overflow in the solution vector.
!
   ELSE IF( IMAT == 11 ) THEN
!
!        Type 11:  Generate a triangular matrix with elements between
!        -1 and 1. Give the diagonal norm 2 to make it well-conditioned.
!        Make the right hand side large so that it requires scaling.
!
      IF( UPPER ) THEN
         JC = 1
         DO J = 1, N
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLARNV( 2, ISEED, J, A( JC ) )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            A( JC+J-1 ) = SIGN( TWO, A( JC+J-1 ) )
            JC = JC + J
            ENDDO
      ELSE
         JC = 1
         DO J = 1, N
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLARNV( 2, ISEED, N-J+1, A( JC ) )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            A( JC ) = SIGN( TWO, A( JC ) )
            JC = JC + N - J + 1
            ENDDO
      END IF
!
!        Set the right hand side so that the largest value is BIGNUM.
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLARNV( 2, ISEED, N, B )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IY = IDAMAX( N, B, 1 )
      BNORM = ABS( B( IY ) )
      BSCAL = BIGNUM / MAX( ONE, BNORM )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DSCAL( N, BSCAL, B, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
   ELSE IF( IMAT == 12 ) THEN
!
!        Type 12:  Make the first diagonal element in the solve small to
!        cause immediate overflow when dividing by T(j,j).
!        In type 12, the offdiagonal elements are small (CNORM(j) < 1).
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLARNV( 2, ISEED, N, B )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      TSCAL = ONE / MAX( ONE, DBLE( N-1 ) )
      IF( UPPER ) THEN
         JC = 1
         DO J = 1, N
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLARNV( 2, ISEED, J-1, A( JC ) )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DSCAL( J-1, TSCAL, A( JC ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            A( JC+J-1 ) = SIGN( ONE, DLARND( 2, ISEED ) )
            JC = JC + J
            ENDDO
         A( N*( N+1 ) / 2 ) = SMLNUM
      ELSE
         JC = 1
         DO J = 1, N
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLARNV( 2, ISEED, N-J, A( JC+1 ) )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DSCAL( N-J, TSCAL, A( JC+1 ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            A( JC ) = SIGN( ONE, DLARND( 2, ISEED ) )
            JC = JC + N - J + 1
            ENDDO
         A( 1 ) = SMLNUM
      END IF
!
   ELSE IF( IMAT == 13 ) THEN
!
!        Type 13:  Make the first diagonal element in the solve small to
!        cause immediate overflow when dividing by T(j,j).
!        In type 13, the offdiagonal elements are O(1) (CNORM(j) > 1).
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLARNV( 2, ISEED, N, B )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( UPPER ) THEN
         JC = 1
         DO J = 1, N
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLARNV( 2, ISEED, J-1, A( JC ) )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            A( JC+J-1 ) = SIGN( ONE, DLARND( 2, ISEED ) )
            JC = JC + J
            ENDDO
         A( N*( N+1 ) / 2 ) = SMLNUM
      ELSE
         JC = 1
         DO J = 1, N
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLARNV( 2, ISEED, N-J, A( JC+1 ) )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            A( JC ) = SIGN( ONE, DLARND( 2, ISEED ) )
            JC = JC + N - J + 1
            ENDDO
         A( 1 ) = SMLNUM
      END IF
!
   ELSE IF( IMAT == 14 ) THEN
!
!        Type 14:  T is diagonal with small numbers on the diagonal to
!        make the growth factor underflow, but a small right hand side
!        chosen so that the solution does not overflow.
!
      IF( UPPER ) THEN
         JCOUNT = 1
         JC = ( N-1 )*N / 2 + 1
         DO J = N, 1, -1
            DO I = 1, J - 1
               A( JC+I-1 ) = ZERO
               ENDDO
            IF( JCOUNT <= 2 ) THEN
               A( JC+J-1 ) = SMLNUM
            ELSE
               A( JC+J-1 ) = ONE
            END IF
            JCOUNT = JCOUNT + 1
            IF( JCOUNT > 4 ) &
               JCOUNT = 1
            JC = JC - J + 1
            ENDDO
      ELSE
         JCOUNT = 1
         JC = 1
         DO J = 1, N
            DO I = J + 1, N
               A( JC+I-J ) = ZERO
               ENDDO
            IF( JCOUNT <= 2 ) THEN
               A( JC ) = SMLNUM
            ELSE
               A( JC ) = ONE
            END IF
            JCOUNT = JCOUNT + 1
            IF( JCOUNT > 4 ) &
               JCOUNT = 1
            JC = JC + N - J + 1
            ENDDO
      END IF
!
!        Set the right hand side alternately zero and small.
!
      IF( UPPER ) THEN
         B( 1 ) = ZERO
         DO I = N, 2, -2
            B( I ) = ZERO
            B( I-1 ) = SMLNUM
            ENDDO
      ELSE
         B( N ) = ZERO
         DO I = 1, N - 1, 2
            B( I ) = ZERO
            B( I+1 ) = SMLNUM
            ENDDO
      END IF
!
   ELSE IF( IMAT == 15 ) THEN
!
!        Type 15:  Make the diagonal elements small to cause gradual
!        overflow when dividing by T(j,j).  To control the amount of
!        scaling needed, the matrix is bidiagonal.
!
      TEXP = ONE / MAX( ONE, DBLE( N-1 ) )
      TSCAL = SMLNUM**TEXP
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLARNV( 2, ISEED, N, B )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( UPPER ) THEN
         JC = 1
         DO J = 1, N
            DO I = 1, J - 2
               A( JC+I-1 ) = ZERO
               ENDDO
            IF( J > 1 ) &
               A( JC+J-2 ) = -ONE
            A( JC+J-1 ) = TSCAL
            JC = JC + J
            ENDDO
         B( N ) = ONE
      ELSE
         JC = 1
         DO J = 1, N
            DO I = J + 2, N
               A( JC+I-J ) = ZERO
               ENDDO
            IF( J < N ) &
               A( JC+1 ) = -ONE
            A( JC ) = TSCAL
            JC = JC + N - J + 1
            ENDDO
         B( 1 ) = ONE
      END IF
!
   ELSE IF( IMAT == 16 ) THEN
!
!        Type 16:  One zero diagonal element.
!
      IY = N / 2 + 1
      IF( UPPER ) THEN
         JC = 1
         DO J = 1, N
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLARNV( 2, ISEED, J, A( JC ) )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            IF( J /= IY ) THEN
               A( JC+J-1 ) = SIGN( TWO, A( JC+J-1 ) )
            ELSE
               A( JC+J-1 ) = ZERO
            END IF
            JC = JC + J
            ENDDO
      ELSE
         JC = 1
         DO J = 1, N
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLARNV( 2, ISEED, N-J+1, A( JC ) )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            IF( J /= IY ) THEN
               A( JC ) = SIGN( TWO, A( JC ) )
            ELSE
               A( JC ) = ZERO
            END IF
            JC = JC + N - J + 1
            ENDDO
      END IF
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLARNV( 2, ISEED, N, B )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DSCAL( N, TWO, B, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
   ELSE IF( IMAT == 17 ) THEN
!
!        Type 17:  Make the offdiagonal elements large to cause overflow
!        when adding a column of T.  In the non-transposed case, the
!        matrix is constructed to cause overflow when adding a column in
!        every other step.
!
      TSCAL = UNFL / ULP
      TSCAL = ( ONE-ULP ) / TSCAL
      DO J = 1, N*( N+1 ) / 2
         A( J ) = ZERO
         ENDDO
      TEXP = ONE
      IF( UPPER ) THEN
         JC = ( N-1 )*N / 2 + 1
         DO J = N, 2, -2
            A( JC ) = -TSCAL / DBLE( N+1 )
            A( JC+J-1 ) = ONE
            B( J ) = TEXP*( ONE-ULP )
            JC = JC - J + 1
            A( JC ) = -( TSCAL / DBLE( N+1 ) ) / DBLE( N+2 )
            A( JC+J-2 ) = ONE
            B( J-1 ) = TEXP*DBLE( N*N+N-1 )
            TEXP = TEXP*TWO
            JC = JC - J + 2
            ENDDO
         B( 1 ) = ( DBLE( N+1 ) / DBLE( N+2 ) )*TSCAL
      ELSE
         JC = 1
         DO J = 1, N - 1, 2
            A( JC+N-J ) = -TSCAL / DBLE( N+1 )
            A( JC ) = ONE
            B( J ) = TEXP*( ONE-ULP )
            JC = JC + N - J + 1
            A( JC+N-J-1 ) = -( TSCAL / DBLE( N+1 ) ) / DBLE( N+2 )
            A( JC ) = ONE
            B( J+1 ) = TEXP*DBLE( N*N+N-1 )
            TEXP = TEXP*TWO
            JC = JC + N - J
            ENDDO
         B( N ) = ( DBLE( N+1 ) / DBLE( N+2 ) )*TSCAL
      END IF
!
   ELSE IF( IMAT == 18 ) THEN
!
!        Type 18:  Generate a unit triangular matrix with elements
!        between -1 and 1, and make the right hand side large so that it
!        requires scaling.
!
      IF( UPPER ) THEN
         JC = 1
         DO J = 1, N
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLARNV( 2, ISEED, J-1, A( JC ) )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            A( JC+J-1 ) = ZERO
            JC = JC + J
            ENDDO
      ELSE
         JC = 1
         DO J = 1, N
            IF( J < N )  THEN
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DLARNV( 2, ISEED, N-J, A( JC+1 ) )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
            ENDIF
            A( JC ) = ZERO
            JC = JC + N - J + 1
            ENDDO
      END IF
!
!        Set the right hand side so that the largest value is BIGNUM.
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLARNV( 2, ISEED, N, B )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IY = IDAMAX( N, B, 1 )
      BNORM = ABS( B( IY ) )
      BSCAL = BIGNUM / MAX( ONE, BNORM )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DSCAL( N, BSCAL, B, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
   ELSE IF( IMAT == 19 ) THEN
!
!        Type 19:  Generate a triangular matrix with elements between
!        BIGNUM/(n-1) and BIGNUM so that at least one of the column
!        norms will exceed BIGNUM.
!
      TLEFT = BIGNUM / MAX( ONE, DBLE( N-1 ) )
      TSCAL = BIGNUM*( DBLE( N-1 ) / MAX( ONE, DBLE( N ) ) )
      IF( UPPER ) THEN
         JC = 1
         DO J = 1, N
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLARNV( 2, ISEED, J, A( JC ) )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            DO I = 1, J
               A( JC+I-1 ) = SIGN( TLEFT, A( JC+I-1 ) ) + &
                             TSCAL*A( JC+I-1 )
               ENDDO
            JC = JC + J
            ENDDO
      ELSE
         JC = 1
         DO J = 1, N
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLARNV( 2, ISEED, N-J+1, A( JC ) )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            DO I = J, N
               A( JC+I-J ) = SIGN( TLEFT, A( JC+I-J ) ) + &
                             TSCAL*A( JC+I-J )
               ENDDO
            JC = JC + N - J + 1
            ENDDO
      END IF
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLARNV( 2, ISEED, N, B )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DSCAL( N, TWO, B, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   END IF
!
!     Flip the matrix across its counter-diagonal if the transpose will
!     be used.
!
   IF( .NOT.LSAME( TRANS, 'N' ) ) THEN
      IF( UPPER ) THEN
         JJ = 1
         JR = N*( N+1 ) / 2
         DO J = 1, N / 2
            JL = JJ
            DO I = J, N - J
               T = A( JR-I+J )
               A( JR-I+J ) = A( JL )
               A( JL ) = T
               JL = JL + I
               ENDDO
            JJ = JJ + J + 1
            JR = JR - ( N-J+1 )
            ENDDO
      ELSE
         JL = 1
         JJ = N*( N+1 ) / 2
         DO J = 1, N / 2
            JR = JJ
            DO I = J, N - J
               T = A( JL+I-J )
               A( JL+I-J ) = A( JR )
               A( JR ) = T
               JR = JR - I
               ENDDO
            JL = JL + N - J + 1
            JJ = JJ - J - 1
            ENDDO
      END IF
   END IF
!
   RETURN
!
!     End of DLATTP
!
END
                                                                                                                                                                                                                                                                                                            




