!> \brief \b CGEBAL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGEBAL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgebal.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgebal.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgebal.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOB
!       INTEGER            IHI, ILO, INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       REAL               SCALE( * )
!       COMPLEX            A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGEBAL balances a general complex matrix A.  This involves, first,
!> permuting A by a similarity transformation to isolate eigenvalues
!> in the first 1 to ILO-1 and last IHI+1 to N elements on the
!> diagonal; and second, applying a diagonal similarity transformation
!> to rows and columns ILO to IHI to make the rows and columns as
!> close in norm as possible.  Both steps are optional.
!>
!> Balancing may reduce the 1-norm of the matrix, and improve the
!> accuracy of the computed eigenvalues and/or eigenvectors.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>          Specifies the operations to be performed on A:
!>          = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0
!>                  for i = 1,...,N;
!>          = 'P':  permute only;
!>          = 'S':  scale only;
!>          = 'B':  both permute and scale.
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
!>          On entry, the input matrix A.
!>          On exit,  A is overwritten by the balanced matrix.
!>          If JOB = 'N', A is not referenced.
!>          See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[out] IHI
!> \verbatim
!>          IHI is INTEGER
!>          ILO and IHI are set to integers such that on exit
!>          A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N.
!>          If JOB = 'N' or 'S', ILO = 1 and IHI = N.
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is REAL array, dimension (N)
!>          Details of the permutations and scaling factors applied to
!>          A.  If P(j) is the index of the row and column interchanged
!>          with row and column j and D(j) is the scaling factor
!>          applied to row and column j, then
!>          SCALE(j) = P(j)    for j = 1,...,ILO-1
!>                   = D(j)    for j = ILO,...,IHI
!>                   = P(j)    for j = IHI+1,...,N.
!>          The order in which the interchanges are made is N to IHI+1,
!>          then 1 to ILO-1.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup gebal
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The permutations consist of row and column interchanges which put
!>  the matrix in the form
!>
!>             ( T1   X   Y  )
!>     P A P = (  0   B   Z  )
!>             (  0   0   T2 )
!>
!>  where T1 and T2 are upper triangular matrices whose eigenvalues lie
!>  along the diagonal.  The column indices ILO and IHI mark the starting
!>  and ending columns of the submatrix B. Balancing consists of applying
!>  a diagonal similarity transformation inv(D) * B * D to make the
!>  1-norms of each row of B and its corresponding column nearly equal.
!>  The output matrix is
!>
!>     ( T1     X*D          Y    )
!>     (  0  inv(D)*B*D  inv(D)*Z ).
!>     (  0      0           T2   )
!>
!>  Information about the permutations P and the diagonal matrix D is
!>  returned in the vector SCALE.
!>
!>  This subroutine is based on the EISPACK routine CBAL.
!>
!>  Modified by Tzu-Yi Chen, Computer Science Division, University of
!>    California at Berkeley, USA
!>
!>  Refactored by Evert Provoost, Department of Computer Science,
!>    KU Leuven, Belgium
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE CGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          JOB
   INTEGER            IHI, ILO, INFO, LDA, N
!     ..
!     .. Array Arguments ..
   REAL               SCALE( * )
   COMPLEX            A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO, ONE
   PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
   REAL               SCLFAC
   PARAMETER          ( SCLFAC = 2.0E+0 )
   REAL               FACTOR
   PARAMETER          ( FACTOR = 0.95E+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            NOCONV, CANSWAP
   INTEGER            I, ICA, IRA, J, K, L
   REAL               C, CA, F, G, R, RA, S, SFMAX1, SFMAX2, SFMIN1, &
                      SFMIN2
!     ..
!     .. External Functions ..
   LOGICAL            SISNAN, LSAME
   INTEGER            ICAMAX
   REAL               SLAMCH, SCNRM2
   EXTERNAL           SISNAN, LSAME, ICAMAX, SLAMCH, SCNRM2
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA, CSSCAL, CSWAP
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, REAL, AIMAG, MAX, MIN
!
!     Test the input parameters
!
   INFO = 0
   IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND. &
       .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -4
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CGEBAL', -INFO )
      RETURN
   END IF
!
!     Quick returns.
!
   IF( N == 0 ) THEN
      ILO = 1
      IHI = 0
      RETURN
   END IF
!
   IF( LSAME( JOB, 'N' ) ) THEN
      DO I = 1, N
         SCALE( I ) = ONE
      END DO
      ILO = 1
      IHI = N
      RETURN
   END IF
!
!     Permutation to isolate eigenvalues if possible.
!
   K = 1
   L = N
!
   IF( .NOT.LSAME( JOB, 'S' ) ) THEN
!
!        Row and column exchange.
!
      NOCONV = .TRUE.
      DO WHILE( NOCONV )
!
!           Search for rows isolating an eigenvalue and push them down.
!
         NOCONV = .FALSE.
         DO I = L, 1, -1
            CANSWAP = .TRUE.
            DO J = 1, L
               IF( I /= J .AND. ( REAL( A( I, J ) ) /= ZERO .OR. &
                   AIMAG( A( I, J ) ) /= ZERO ) ) THEN
                  CANSWAP = .FALSE.
                  EXIT
               END IF
            END DO
!
            IF( CANSWAP ) THEN
               SCALE( L ) = I
               IF( I /= L ) THEN
                  CALL CSWAP( L, A( 1, I ), 1, A( 1, L ), 1 )
                  CALL CSWAP( N-K+1, A( I, K ), LDA, A( L, K ), LDA )
               END IF
               NOCONV = .TRUE.
!
               IF( L == 1 ) THEN
                  ILO = 1
                  IHI = 1
                  RETURN
               END IF
!
               L = L - 1
            END IF
         END DO
!
      END DO

      NOCONV = .TRUE.
      DO WHILE( NOCONV )
!
!           Search for columns isolating an eigenvalue and push them left.
!
         NOCONV = .FALSE.
         DO J = K, L
            CANSWAP = .TRUE.
            DO I = K, L
               IF( I /= J .AND. ( REAL( A( I, J ) ) /= ZERO .OR. &
                   AIMAG( A( I, J ) ) /= ZERO ) ) THEN
                  CANSWAP = .FALSE.
                  EXIT
               END IF
            END DO
!
            IF( CANSWAP ) THEN
               SCALE( K ) = J
               IF( J /= K ) THEN
                  CALL CSWAP( L, A( 1, J ), 1, A( 1, K ), 1 )
                  CALL CSWAP( N-K+1, A( J, K ), LDA, A( K, K ), LDA )
               END IF
               NOCONV = .TRUE.
!
               K = K + 1
            END IF
         END DO
!
      END DO
!
   END IF
!
!     Initialize SCALE for non-permuted submatrix.
!
   DO I = K, L
      SCALE( I ) = ONE
   END DO
!
!     If we only had to permute, we are done.
!
   IF( LSAME( JOB, 'P' ) ) THEN
      ILO = K
      IHI = L
      RETURN
   END IF
!
!     Balance the submatrix in rows K to L.
!
!     Iterative loop for norm reduction.
!
   SFMIN1 = SLAMCH( 'S' ) / SLAMCH( 'P' )
   SFMAX1 = ONE / SFMIN1
   SFMIN2 = SFMIN1*SCLFAC
   SFMAX2 = ONE / SFMIN2
!
   NOCONV = .TRUE.
   DO WHILE( NOCONV )
      NOCONV = .FALSE.
!
      DO I = K, L
!
         C = SCNRM2( L-K+1, A( K, I ), 1 )
         R = SCNRM2( L-K+1, A( I, K ), LDA )
         ICA = ICAMAX( L, A( 1, I ), 1 )
         CA = ABS( A( ICA, I ) )
         IRA = ICAMAX( N-K+1, A( I, K ), LDA )
         RA = ABS( A( I, IRA+K-1 ) )
!
!           Guard against zero C or R due to underflow.
!
         IF( C == ZERO .OR. R == ZERO ) CYCLE
!
!           Exit if NaN to avoid infinite loop
!
         IF( SISNAN( C+CA+R+RA ) ) THEN
            INFO = -3
            CALL XERBLA( 'CGEBAL', -INFO )
            RETURN
         END IF
!
         G = R / SCLFAC
         F = ONE
         S = C + R
!
         DO WHILE( C < G .AND. MAX( F, C, CA ) < SFMAX2 .AND. &
                   MIN( R, G, RA ) > SFMIN2 )
            F = F*SCLFAC
            C = C*SCLFAC
            CA = CA*SCLFAC
            R = R / SCLFAC
            G = G / SCLFAC
            RA = RA / SCLFAC
         END DO
!
         G = C / SCLFAC
!
         DO WHILE( G >= R .AND. MAX( R, RA ) < SFMAX2 .AND. &
                   MIN( F, C, G, CA ) > SFMIN2 )
            F = F / SCLFAC
            C = C / SCLFAC
            G = G / SCLFAC
            CA = CA / SCLFAC
            R = R*SCLFAC
            RA = RA*SCLFAC
         END DO
!
!           Now balance.
!
         IF( ( C+R ) >= FACTOR*S ) CYCLE
         IF( F < ONE .AND. SCALE( I ) < ONE ) THEN
            IF( F*SCALE( I ) <= SFMIN1 ) CYCLE
         END IF
         IF( F > ONE .AND. SCALE( I ) > ONE ) THEN
            IF( SCALE( I ) >= SFMAX1 / F ) CYCLE
         END IF
         G = ONE / F
         SCALE( I ) = SCALE( I )*F
         NOCONV = .TRUE.
!
         CALL CSSCAL( N-K+1, G, A( I, K ), LDA )
         CALL CSSCAL( L, F, A( 1, I ), 1 )
!
      END DO
!
   END DO
!
   ILO = K
   IHI = L
!
   RETURN
!
!     End of CGEBAL
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
