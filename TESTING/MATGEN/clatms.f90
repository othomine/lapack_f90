!> \brief \b CLATMS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLATMS( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX,
!                          KL, KU, PACK, A, LDA, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIST, PACK, SYM
!       INTEGER            INFO, KL, KU, LDA, M, MODE, N
!       REAL               COND, DMAX
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       REAL               D( * )
!       COMPLEX            A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CLATMS generates random matrices with specified singular values
!>    (or hermitian with specified eigenvalues)
!>    for testing LAPACK programs.
!>
!>    CLATMS operates by applying the following sequence of
!>    operations:
!>
!>      Set the diagonal to D, where D may be input or
!>         computed according to MODE, COND, DMAX, and SYM
!>         as described below.
!>
!>      Generate a matrix with the appropriate band structure, by one
!>         of two methods:
!>
!>      Method A:
!>          Generate a dense M x N matrix by multiplying D on the left
!>              and the right by random unitary matrices, then:
!>
!>          Reduce the bandwidth according to KL and KU, using
!>              Householder transformations.
!>
!>      Method B:
!>          Convert the bandwidth-0 (i.e., diagonal) matrix to a
!>              bandwidth-1 matrix using Givens rotations, "chasing"
!>              out-of-band elements back, much as in QR; then convert
!>              the bandwidth-1 to a bandwidth-2 matrix, etc.  Note
!>              that for reasonably small bandwidths (relative to M and
!>              N) this requires less storage, as a dense matrix is not
!>              generated.  Also, for hermitian or symmetric matrices,
!>              only one triangle is generated.
!>
!>      Method A is chosen if the bandwidth is a large fraction of the
!>          order of the matrix, and LDA is at least M (so a dense
!>          matrix can be stored.)  Method B is chosen if the bandwidth
!>          is small (< 1/2 N for hermitian or symmetric, < .3 N+M for
!>          non-symmetric), or LDA is less than M and not less than the
!>          bandwidth.
!>
!>      Pack the matrix if desired. Options specified by PACK are:
!>         no packing
!>         zero out upper half (if hermitian)
!>         zero out lower half (if hermitian)
!>         store the upper half columnwise (if hermitian or upper
!>               triangular)
!>         store the lower half columnwise (if hermitian or lower
!>               triangular)
!>         store the lower triangle in banded format (if hermitian or
!>               lower triangular)
!>         store the upper triangle in banded format (if hermitian or
!>               upper triangular)
!>         store the entire matrix in banded format
!>      If Method B is chosen, and band format is specified, then the
!>         matrix will be generated in the band format, so no repacking
!>         will be necessary.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           The number of rows of A. Not modified.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           The number of columns of A. N must equal M if the matrix
!>           is symmetric or hermitian (i.e., if SYM is not 'N')
!>           Not modified.
!> \endverbatim
!>
!> \param[in] DIST
!> \verbatim
!>          DIST is CHARACTER*1
!>           On entry, DIST specifies the type of distribution to be used
!>           to generate the random eigen-/singular values.
!>           'U' => UNIFORM( 0, 1 )  ( 'U' for uniform )
!>           'S' => UNIFORM( -1, 1 ) ( 'S' for symmetric )
!>           'N' => NORMAL( 0, 1 )   ( 'N' for normal )
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension ( 4 )
!>           On entry ISEED specifies the seed of the random number
!>           generator. They should lie between 0 and 4095 inclusive,
!>           and ISEED(4) should be odd. The random number generator
!>           uses a linear congruential sequence limited to small
!>           integers, and so should produce machine independent
!>           random numbers. The values of ISEED are changed on
!>           exit, and can be used in the next call to CLATMS
!>           to continue the same random number sequence.
!>           Changed on exit.
!> \endverbatim
!>
!> \param[in] SYM
!> \verbatim
!>          SYM is CHARACTER*1
!>           If SYM='H', the generated matrix is hermitian, with
!>             eigenvalues specified by D, COND, MODE, and DMAX; they
!>             may be positive, negative, or zero.
!>           If SYM='P', the generated matrix is hermitian, with
!>             eigenvalues (= singular values) specified by D, COND,
!>             MODE, and DMAX; they will not be negative.
!>           If SYM='N', the generated matrix is nonsymmetric, with
!>             singular values specified by D, COND, MODE, and DMAX;
!>             they will not be negative.
!>           If SYM='S', the generated matrix is (complex) symmetric,
!>             with singular values specified by D, COND, MODE, and
!>             DMAX; they will not be negative.
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension ( MIN( M, N ) )
!>           This array is used to specify the singular values or
!>           eigenvalues of A (see SYM, above.)  If MODE=0, then D is
!>           assumed to contain the singular/eigenvalues, otherwise
!>           they will be computed according to MODE, COND, and DMAX,
!>           and placed in D.
!>           Modified if MODE is nonzero.
!> \endverbatim
!>
!> \param[in] MODE
!> \verbatim
!>          MODE is INTEGER
!>           On entry this describes how the singular/eigenvalues are to
!>           be specified:
!>           MODE = 0 means use D as input
!>           MODE = 1 sets D(1)=1 and D(2:N)=1.0/COND
!>           MODE = 2 sets D(1:N-1)=1 and D(N)=1.0/COND
!>           MODE = 3 sets D(I)=COND**(-(I-1)/(N-1))
!>           MODE = 4 sets D(i)=1 - (i-1)/(N-1)*(1 - 1/COND)
!>           MODE = 5 sets D to random numbers in the range
!>                    ( 1/COND , 1 ) such that their logarithms
!>                    are uniformly distributed.
!>           MODE = 6 set D to random numbers from same distribution
!>                    as the rest of the matrix.
!>           MODE < 0 has the same meaning as ABS(MODE), except that
!>              the order of the elements of D is reversed.
!>           Thus if MODE is positive, D has entries ranging from
!>              1 to 1/COND, if negative, from 1/COND to 1,
!>           If SYM='H', and MODE is neither 0, 6, nor -6, then
!>              the elements of D will also be multiplied by a random
!>              sign (i.e., +1 or -1.)
!>           Not modified.
!> \endverbatim
!>
!> \param[in] COND
!> \verbatim
!>          COND is REAL
!>           On entry, this is used as described under MODE above.
!>           If used, it must be >= 1. Not modified.
!> \endverbatim
!>
!> \param[in] DMAX
!> \verbatim
!>          DMAX is REAL
!>           If MODE is neither -6, 0 nor 6, the contents of D, as
!>           computed according to MODE and COND, will be scaled by
!>           DMAX / max(abs(D(i))); thus, the maximum absolute eigen- or
!>           singular value (which is to say the norm) will be abs(DMAX).
!>           Note that DMAX need not be positive: if DMAX is negative
!>           (or zero), D will be scaled by a negative number (or zero).
!>           Not modified.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>           This specifies the lower bandwidth of the  matrix. For
!>           example, KL=0 implies upper triangular, KL=1 implies upper
!>           Hessenberg, and KL being at least M-1 means that the matrix
!>           has full lower bandwidth.  KL must equal KU if the matrix
!>           is symmetric or hermitian.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>           This specifies the upper bandwidth of the  matrix. For
!>           example, KU=0 implies lower triangular, KU=1 implies lower
!>           Hessenberg, and KU being at least N-1 means that the matrix
!>           has full upper bandwidth.  KL must equal KU if the matrix
!>           is symmetric or hermitian.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] PACK
!> \verbatim
!>          PACK is CHARACTER*1
!>           This specifies packing of matrix as follows:
!>           'N' => no packing
!>           'U' => zero out all subdiagonal entries (if symmetric
!>                  or hermitian)
!>           'L' => zero out all superdiagonal entries (if symmetric
!>                  or hermitian)
!>           'C' => store the upper triangle columnwise (only if the
!>                  matrix is symmetric, hermitian, or upper triangular)
!>           'R' => store the lower triangle columnwise (only if the
!>                  matrix is symmetric, hermitian, or lower triangular)
!>           'B' => store the lower triangle in band storage scheme
!>                  (only if the matrix is symmetric, hermitian, or
!>                  lower triangular)
!>           'Q' => store the upper triangle in band storage scheme
!>                  (only if the matrix is symmetric, hermitian, or
!>                  upper triangular)
!>           'Z' => store the entire matrix in band storage scheme
!>                      (pivoting can be provided for by using this
!>                      option to store A in the trailing rows of
!>                      the allocated storage)
!>
!>           Using these options, the various LAPACK packed and banded
!>           storage schemes can be obtained:
!>           GB                    - use 'Z'
!>           PB, SB, HB, or TB     - use 'B' or 'Q'
!>           PP, SP, HB, or TP     - use 'C' or 'R'
!>
!>           If two calls to CLATMS differ only in the PACK parameter,
!>           they will generate mathematically equivalent matrices.
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension ( LDA, N )
!>           On exit A is the desired test matrix.  A is first generated
!>           in full (unpacked) form, and then packed, if so specified
!>           by PACK.  Thus, the first M elements of the first N
!>           columns will always be modified.  If PACK specifies a
!>           packed or banded storage scheme, all LDA elements of the
!>           first N columns will be modified; the elements of the
!>           array which do not correspond to elements of the generated
!>           matrix are set to zero.
!>           Modified.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           LDA specifies the first dimension of A as declared in the
!>           calling program.  If PACK='N', 'U', 'L', 'C', or 'R', then
!>           LDA must be at least M.  If PACK='B' or 'Q', then LDA must
!>           be at least MIN( KL, M-1) (which is equal to MIN(KU,N-1)).
!>           If PACK='Z', LDA must be large enough to hold the packed
!>           array: MIN( KU, N-1) + MIN( KL, M-1) + 1.
!>           Not modified.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension ( 3*MAX( N, M ) )
!>           Workspace.
!>           Modified.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           Error code.  On exit, INFO will be set to one of the
!>           following values:
!>             0 => normal return
!>            -1 => M negative or unequal to N and SYM='S', 'H', or 'P'
!>            -2 => N negative
!>            -3 => DIST illegal string
!>            -5 => SYM illegal string
!>            -7 => MODE not in range -6 to 6
!>            -8 => COND less than 1.0, and MODE neither -6, 0 nor 6
!>           -10 => KL negative
!>           -11 => KU negative, or SYM is not 'N' and KU is not equal to
!>                  KL
!>           -12 => PACK illegal string, or PACK='U' or 'L', and SYM='N';
!>                  or PACK='C' or 'Q' and SYM='N' and KL is not zero;
!>                  or PACK='R' or 'B' and SYM='N' and KU is not zero;
!>                  or PACK='U', 'L', 'C', 'R', 'B', or 'Q', and M is not
!>                  N.
!>           -14 => LDA is less than M, or PACK='Z' and LDA is less than
!>                  MIN(KU,N-1) + MIN(KL,M-1) + 1.
!>            1  => Error return from SLATM1
!>            2  => Cannot scale to DMAX (max. sing. value is 0)
!>            3  => Error return from CLAGGE, CLAGHE or CLAGSY
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
!  =====================================================================
   SUBROUTINE CLATMS( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX, &
                      KL, KU, PACK, A, LDA, WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          DIST, PACK, SYM
   INTEGER            INFO, KL, KU, LDA, M, MODE, N
   REAL               COND, DMAX
!     ..
!     .. Array Arguments ..
   INTEGER            ISEED( 4 )
   REAL               D( * )
   COMPLEX            A( LDA, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               TWOPI
   PARAMETER  ( TWOPI = 6.28318530717958647692528676655900576839E+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            CSYM, GIVENS, ILEXTR, ILTEMP, TOPDWN
   INTEGER            I, IC, ICOL, IDIST, IENDCH, IINFO, IL, ILDA, &
                      IOFFG, IOFFST, IPACK, IPACKG, IR, IR1, IR2, &
                      IROW, IRSIGN, ISKEW, ISYM, ISYMPK, J, JC, JCH, &
                      JKL, JKU, JR, K, LLB, MINLDA, MNMIN, MR, NC, &
                      UUB
   REAL               ALPHA, ANGLE, REALC, TEMP
   COMPLEX            C, CT, CTEMP, DUMMY, EXTRA, S, ST
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               SLARND
   COMPLEX            CLARND
   EXTERNAL           LSAME, SLARND, CLARND
!     ..
!     .. External Subroutines ..
   EXTERNAL           CLAGGE, CLAGHE, CLAGSY, CLAROT, CLARTG, CLASET, &
                      SLATM1, SSCAL, XERBLA
!     ..
!     .. Executable Statements ..
!
!     1)      Decode and Test the input parameters.
!             Initialize flags & seed.
!
   INFO = 0
!
!     Quick return if possible
!
   IF( M == 0 .OR. N == 0 ) RETURN
!
!     Decode DIST
!
   IF( LSAME( DIST, 'U' ) ) THEN
      IDIST = 1
   ELSE IF( LSAME( DIST, 'S' ) ) THEN
      IDIST = 2
   ELSE IF( LSAME( DIST, 'N' ) ) THEN
      IDIST = 3
   ELSE
      IDIST = -1
   END IF
!
!     Decode SYM
!
   IF( LSAME( SYM, 'N' ) ) THEN
      ISYM = 1
      IRSIGN = 0
      CSYM = .FALSE.
   ELSE IF( LSAME( SYM, 'P' ) ) THEN
      ISYM = 2
      IRSIGN = 0
      CSYM = .FALSE.
   ELSE IF( LSAME( SYM, 'S' ) ) THEN
      ISYM = 2
      IRSIGN = 0
      CSYM = .TRUE.
   ELSE IF( LSAME( SYM, 'H' ) ) THEN
      ISYM = 2
      IRSIGN = 1
      CSYM = .FALSE.
   ELSE
      ISYM = -1
   END IF
!
!     Decode PACK
!
   ISYMPK = 0
   IF( LSAME( PACK, 'N' ) ) THEN
      IPACK = 0
   ELSE IF( LSAME( PACK, 'U' ) ) THEN
      IPACK = 1
      ISYMPK = 1
   ELSE IF( LSAME( PACK, 'L' ) ) THEN
      IPACK = 2
      ISYMPK = 1
   ELSE IF( LSAME( PACK, 'C' ) ) THEN
      IPACK = 3
      ISYMPK = 2
   ELSE IF( LSAME( PACK, 'R' ) ) THEN
      IPACK = 4
      ISYMPK = 3
   ELSE IF( LSAME( PACK, 'B' ) ) THEN
      IPACK = 5
      ISYMPK = 3
   ELSE IF( LSAME( PACK, 'Q' ) ) THEN
      IPACK = 6
      ISYMPK = 2
   ELSE IF( LSAME( PACK, 'Z' ) ) THEN
      IPACK = 7
   ELSE
      IPACK = -1
   END IF
!
!     Set certain internal parameters
!
   MNMIN = MIN( M, N )
   LLB = MIN( KL, M-1 )
   UUB = MIN( KU, N-1 )
   MR = MIN( M, N+LLB )
   NC = MIN( N, M+UUB )
!
   IF( IPACK == 5 .OR. IPACK == 6 ) THEN
      MINLDA = UUB + 1
   ELSE IF( IPACK == 7 ) THEN
      MINLDA = LLB + UUB + 1
   ELSE
      MINLDA = M
   END IF
!
!     Use Givens rotation method if bandwidth small enough,
!     or if LDA is too small to store the matrix unpacked.
!
   IF( ISYM == 1 ) THEN
      GIVENS = ( REAL( LLB+UUB ) < 0.3*REAL( MAX( 1, MR+NC ) ) )
   ELSE
      GIVENS = ( 2*LLB < M )
   END IF
   IF( LDA < M .AND. LDA >= MINLDA ) GIVENS = .TRUE.
!
!     Set INFO if an error
!
   IF( M < 0 ) THEN
      INFO = -1
   ELSE IF( M /= N .AND. ISYM /= 1 ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( IDIST == -1 ) THEN
      INFO = -3
   ELSE IF( ISYM == -1 ) THEN
      INFO = -5
   ELSE IF( ABS( MODE ) > 6 ) THEN
      INFO = -7
   ELSE IF( ( MODE /= 0 .AND. ABS( MODE ) /= 6 ) .AND. COND < 1.0E+0 ) &
             THEN
      INFO = -8
   ELSE IF( KL < 0 ) THEN
      INFO = -10
   ELSE IF( KU < 0 .OR. ( ISYM /= 1 .AND. KL /= KU ) ) THEN
      INFO = -11
   ELSE IF( IPACK == -1 .OR. ( ISYMPK == 1 .AND. ISYM == 1 ) .OR. &
            ( ISYMPK == 2 .AND. ISYM == 1 .AND. KL > 0 ) .OR. &
            ( ISYMPK == 3 .AND. ISYM == 1 .AND. KU > 0 ) .OR. &
            ( ISYMPK /= 0 .AND. M /= N ) ) THEN
      INFO = -12
   ELSE IF( LDA < MAX( 1, MINLDA ) ) THEN
      INFO = -14
   END IF
!
   IF( INFO /= 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'CLATMS', -INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : XERBLA : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      RETURN
   END IF
!
!     Initialize random number generator
!
   ISEED(1:4) = MOD( ABS( ISEED(1:4) ), 4096 )
!
   IF( MOD( ISEED( 4 ), 2 ) /= 1 ) ISEED( 4 ) = ISEED( 4 ) + 1
!
!     2)      Set up D  if indicated.
!
!             Compute D according to COND and MODE
!
   CALL SLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, IINFO )
   IF( IINFO /= 0 ) THEN
      INFO = 1
      RETURN
   END IF
!
!     Choose Top-Down if D is (apparently) increasing,
!     Bottom-Up if D is (apparently) decreasing.
!
   TOPDWN = ( ABS( D( 1 ) ) <= ABS( D( MNMIN ) ) )
!
   IF( MODE /= 0 .AND. ABS( MODE ) /= 6 ) THEN
!
!        Scale by DMAX
!
      TEMP = MAXVAL(ABS( D(1:MNMIN) ) )
!
      IF( TEMP > 0.0E+0 ) THEN
         ALPHA = DMAX / TEMP
      ELSE
         INFO = 2
         RETURN
      END IF
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL SSCAL( MNMIN, ALPHA, D, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : SSCAL : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
   END IF
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CLASET( 'Full', LDA, N, (0.0E+0,0.0E+0), (0.0E+0,0.0E+0), A, LDA )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     3)      Generate Banded Matrix using Givens rotations.
!             Also the special case of UUB=LLB=0
!
!               Compute Addressing constants to cover all
!               storage formats.  Whether GE, HE, SY, GB, HB, or SB,
!               upper or lower triangle or both,
!               the (i,j)-th element is in
!               A( i - ISKEW*j + IOFFST, j )
!
   IF( IPACK > 4 ) THEN
      ILDA = LDA - 1
      ISKEW = 1
      IF( IPACK > 5 ) THEN
         IOFFST = UUB + 1
      ELSE
         IOFFST = 1
      END IF
   ELSE
      ILDA = LDA
      ISKEW = 0
      IOFFST = 0
   END IF
!
!     IPACKG is the format that the matrix is generated in. If this is
!     different from IPACK, then the matrix must be repacked at the
!     end.  It also signals how to compute the norm, for scaling.
!
   IPACKG = 0
!
!     Diagonal Matrix -- We are done, unless it
!     is to be stored HP/SP/PP/TP (PACK='R' or 'C')
!
   IF( LLB == 0 .AND. UUB == 0 ) THEN
      DO J = 1, MNMIN
         A( ( 1-ISKEW )*J+IOFFST, J ) = CMPLX( D( J ) )
      ENDDO
!
      IF( IPACK <= 2 .OR. IPACK >= 5 ) IPACKG = IPACK
!
   ELSE IF( GIVENS ) THEN
!
!        Check whether to use Givens rotations,
!        Householder transformations, or nothing.
!
      IF( ISYM == 1 ) THEN
!
!           Non-symmetric -- A = U D V
!
         IF( IPACK > 4 ) THEN
            IPACKG = IPACK
         ELSE
            IPACKG = 0
         END IF
!
         DO J = 1, MNMIN
            A( ( 1-ISKEW )*J+IOFFST, J ) = CMPLX( D( J ) )
         ENDDO
!
         IF( TOPDWN ) THEN
            JKL = 0
            DO JKU = 1, UUB
!
!                 Transform from bandwidth JKL, JKU-1 to JKL, JKU
!
!                 Last row actually rotated is M
!                 Last column actually rotated is MIN( M+JKU, N )
!
               DO JR = 1, MIN( M+JKU, N ) + JKL - 1
                  EXTRA = (0.0E+0,0.0E+0)
                  ANGLE = TWOPI*SLARND( 1, ISEED )
                  C = COS( ANGLE )*CLARND( 5, ISEED )
                  S = SIN( ANGLE )*CLARND( 5, ISEED )
                  ICOL = MAX( 1, JR-JKL )
                  IF( JR < M ) THEN
                     IL = MIN( N, JR+JKU ) + 1 - ICOL
                     CALL CLAROT( .TRUE., JR > JKL, .FALSE., IL, C, &
                                  S, A( JR-ISKEW*ICOL+IOFFST, ICOL ), &
                                  ILDA, EXTRA, DUMMY )
                  END IF
!
!                    Chase "EXTRA" back up
!
                  IR = JR
                  IC = ICOL
                  DO JCH = JR - JKL, 1, -JKL - JKU
                     IF( IR < M ) THEN
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CLARTG( A( IR+1-ISKEW*( IC+1 )+IOFFST, &
                                     IC+1 ), EXTRA, REALC, S, DUMMY )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CLARTG : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                        DUMMY = CLARND( 5, ISEED )
                        C = CONJG( REALC*DUMMY )
                        S = CONJG( -S*DUMMY )
                     END IF
                     IROW = MAX( 1, JCH-JKU )
                     IL = IR + 2 - IROW
                     CTEMP = (0.0E+0,0.0E+0)
                     ILTEMP = JCH > JKU
                     CALL CLAROT( .FALSE., ILTEMP, .TRUE., IL, C, S, &
                                  A( IROW-ISKEW*IC+IOFFST, IC ), &
                                  ILDA, CTEMP, EXTRA )
                     IF( ILTEMP ) THEN
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CLARTG( A( IROW+1-ISKEW*( IC+1 )+IOFFST, &
                                     IC+1 ), CTEMP, REALC, S, DUMMY )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CLARTG : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                        DUMMY = CLARND( 5, ISEED )
                        C = CONJG( REALC*DUMMY )
                        S = CONJG( -S*DUMMY )
!
                        ICOL = MAX( 1, JCH-JKU-JKL )
                        IL = IC + 2 - ICOL
                        EXTRA = (0.0E+0,0.0E+0)
                        CALL CLAROT( .TRUE., JCH > JKU+JKL, .TRUE., &
                                     IL, C, S, A( IROW-ISKEW*ICOL+ &
                                     IOFFST, ICOL ), ILDA, EXTRA, &
                                     CTEMP )
                        IC = ICOL
                        IR = IROW
                     END IF
                  ENDDO
               ENDDO
            ENDDO
!
            JKU = UUB
            DO JKL = 1, LLB
!
!                 Transform from bandwidth JKL-1, JKU to JKL, JKU
!
               DO JC = 1, MIN( N+JKL, M ) + JKU - 1
                  EXTRA = (0.0E+0,0.0E+0)
                  ANGLE = TWOPI*SLARND( 1, ISEED )
                  C = COS( ANGLE )*CLARND( 5, ISEED )
                  S = SIN( ANGLE )*CLARND( 5, ISEED )
                  IROW = MAX( 1, JC-JKU )
                  IF( JC < N ) THEN
                     IL = MIN( M, JC+JKL ) + 1 - IROW
                     CALL CLAROT( .FALSE., JC > JKU, .FALSE., IL, C, &
                                  S, A( IROW-ISKEW*JC+IOFFST, JC ), &
                                  ILDA, EXTRA, DUMMY )
                  END IF
!
!                    Chase "EXTRA" back up
!
                  IC = JC
                  IR = IROW
                  DO JCH = JC - JKU, 1, -JKL - JKU
                     IF( IC < N ) THEN
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CLARTG( A( IR+1-ISKEW*( IC+1 )+IOFFST, &
                                     IC+1 ), EXTRA, REALC, S, DUMMY )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CLARTG : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                        DUMMY = CLARND( 5, ISEED )
                        C = CONJG( REALC*DUMMY )
                        S = CONJG( -S*DUMMY )
                     END IF
                     ICOL = MAX( 1, JCH-JKL )
                     IL = IC + 2 - ICOL
                     CTEMP = (0.0E+0,0.0E+0)
                     ILTEMP = JCH > JKL
                     CALL CLAROT( .TRUE., ILTEMP, .TRUE., IL, C, S, &
                                  A( IR-ISKEW*ICOL+IOFFST, ICOL ), &
                                  ILDA, CTEMP, EXTRA )
                     IF( ILTEMP ) THEN
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CLARTG( A( IR+1-ISKEW*( ICOL+1 )+IOFFST, &
                                     ICOL+1 ), CTEMP, REALC, S, &
                                     DUMMY )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CLARTG : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                        DUMMY = CLARND( 5, ISEED )
                        C = CONJG( REALC*DUMMY )
                        S = CONJG( -S*DUMMY )
                        IROW = MAX( 1, JCH-JKL-JKU )
                        IL = IR + 2 - IROW
                        EXTRA = (0.0E+0,0.0E+0)
                        CALL CLAROT( .FALSE., JCH > JKL+JKU, .TRUE., &
                                     IL, C, S, A( IROW-ISKEW*ICOL+ &
                                     IOFFST, ICOL ), ILDA, EXTRA, &
                                     CTEMP )
                        IC = ICOL
                        IR = IROW
                     END IF
                  ENDDO
               ENDDO
            ENDDO
!
         ELSE
!
!              Bottom-Up -- Start at the bottom right.
!
            JKL = 0
            DO JKU = 1, UUB
!
!                 Transform from bandwidth JKL, JKU-1 to JKL, JKU
!
!                 First row actually rotated is M
!                 First column actually rotated is MIN( M+JKU, N )
!
               IENDCH = MIN( M, N+JKL ) - 1
               DO JC = MIN( M+JKU, N ) - 1, 1 - JKL, -1
                  EXTRA = (0.0E+0,0.0E+0)
                  ANGLE = TWOPI*SLARND( 1, ISEED )
                  C = COS( ANGLE )*CLARND( 5, ISEED )
                  S = SIN( ANGLE )*CLARND( 5, ISEED )
                  IROW = MAX( 1, JC-JKU+1 )
                  IF( JC > 0 ) THEN
                     IL = MIN( M, JC+JKL+1 ) + 1 - IROW
                     CALL CLAROT( .FALSE., .FALSE., JC+JKL < M, IL, &
                                  C, S, A( IROW-ISKEW*JC+IOFFST, &
                                  JC ), ILDA, DUMMY, EXTRA )
                  END IF
!
!                    Chase "EXTRA" back down
!
                  IC = JC
                  DO JCH = JC + JKL, IENDCH, JKL + JKU
                     ILEXTR = IC > 0
                     IF( ILEXTR ) THEN
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CLARTG( A( JCH-ISKEW*IC+IOFFST, IC ), &
                                     EXTRA, REALC, S, DUMMY )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CLARTG : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                        DUMMY = CLARND( 5, ISEED )
                        C = REALC*DUMMY
                        S = S*DUMMY
                     END IF
                     IC = MAX( 1, IC )
                     ICOL = MIN( N-1, JCH+JKU )
                     ILTEMP = JCH + JKU < N
                     CTEMP = (0.0E+0,0.0E+0)
                     CALL CLAROT( .TRUE., ILEXTR, ILTEMP, ICOL+2-IC, &
                                  C, S, A( JCH-ISKEW*IC+IOFFST, IC ), &
                                  ILDA, EXTRA, CTEMP )
                     IF( ILTEMP ) THEN
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CLARTG( A( JCH-ISKEW*ICOL+IOFFST, &
                                     ICOL ), CTEMP, REALC, S, DUMMY )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CLARTG : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                        DUMMY = CLARND( 5, ISEED )
                        C = REALC*DUMMY
                        S = S*DUMMY
                        IL = MIN( IENDCH, JCH+JKL+JKU ) + 2 - JCH
                        EXTRA = (0.0E+0,0.0E+0)
                        CALL CLAROT( .FALSE., .TRUE., &
                                     JCH+JKL+JKU <= IENDCH, IL, C, S, &
                                     A( JCH-ISKEW*ICOL+IOFFST, &
                                     ICOL ), ILDA, CTEMP, EXTRA )
                        IC = ICOL
                     END IF
                  ENDDO
               ENDDO
            ENDDO
!
            JKU = UUB
            DO JKL = 1, LLB
!
!                 Transform from bandwidth JKL-1, JKU to JKL, JKU
!
!                 First row actually rotated is MIN( N+JKL, M )
!                 First column actually rotated is N
!
               IENDCH = MIN( N, M+JKU ) - 1
               DO JR = MIN( N+JKL, M ) - 1, 1 - JKU, -1
                  EXTRA = (0.0E+0,0.0E+0)
                  ANGLE = TWOPI*SLARND( 1, ISEED )
                  C = COS( ANGLE )*CLARND( 5, ISEED )
                  S = SIN( ANGLE )*CLARND( 5, ISEED )
                  ICOL = MAX( 1, JR-JKL+1 )
                  IF( JR > 0 ) THEN
                     IL = MIN( N, JR+JKU+1 ) + 1 - ICOL
                     CALL CLAROT( .TRUE., .FALSE., JR+JKU < N, IL, &
                                  C, S, A( JR-ISKEW*ICOL+IOFFST, &
                                  ICOL ), ILDA, DUMMY, EXTRA )
                  END IF
!
!                    Chase "EXTRA" back down
!
                  IR = JR
                  DO JCH = JR + JKU, IENDCH, JKL + JKU
                     ILEXTR = IR > 0
                     IF( ILEXTR ) THEN
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CLARTG( A( IR-ISKEW*JCH+IOFFST, JCH ), &
                                     EXTRA, REALC, S, DUMMY )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CLARTG : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                        DUMMY = CLARND( 5, ISEED )
                        C = REALC*DUMMY
                        S = S*DUMMY
                     END IF
                     IR = MAX( 1, IR )
                     IROW = MIN( M-1, JCH+JKL )
                     ILTEMP = JCH + JKL < M
                     CTEMP = (0.0E+0,0.0E+0)
                     CALL CLAROT( .FALSE., ILEXTR, ILTEMP, IROW+2-IR, &
                                  C, S, A( IR-ISKEW*JCH+IOFFST, &
                                  JCH ), ILDA, EXTRA, CTEMP )
                     IF( ILTEMP ) THEN
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL CLARTG( A( IROW-ISKEW*JCH+IOFFST, JCH ), &
                                     CTEMP, REALC, S, DUMMY )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : CLARTG : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                        DUMMY = CLARND( 5, ISEED )
                        C = REALC*DUMMY
                        S = S*DUMMY
                        IL = MIN( IENDCH, JCH+JKL+JKU ) + 2 - JCH
                        EXTRA = (0.0E+0,0.0E+0)
                        CALL CLAROT( .TRUE., .TRUE., &
                                     JCH+JKL+JKU <= IENDCH, IL, C, S, &
                                     A( IROW-ISKEW*JCH+IOFFST, JCH ), &
                                     ILDA, CTEMP, EXTRA )
                        IR = IROW
                     END IF
                  ENDDO
               ENDDO
            ENDDO
!
         END IF
!
      ELSE
!
!           Symmetric -- A = U D U'
!           Hermitian -- A = U D U*
!
         IPACKG = IPACK
         IOFFG = IOFFST
!
         IF( TOPDWN ) THEN
!
!              Top-Down -- Generate Upper triangle only
!
            IF( IPACK >= 5 ) THEN
               IPACKG = 6
               IOFFG = UUB + 1
            ELSE
               IPACKG = 1
            END IF
!
            DO J = 1, MNMIN
               A( ( 1-ISKEW )*J+IOFFG, J ) = CMPLX( D( J ) )
            ENDDO
!
            DO K = 1, UUB
               DO JC = 1, N - 1
                  IROW = MAX( 1, JC-K )
                  IL = MIN( JC+1, K+2 )
                  EXTRA = (0.0E+0,0.0E+0)
                  CTEMP = A( JC-ISKEW*( JC+1 )+IOFFG, JC+1 )
                  ANGLE = TWOPI*SLARND( 1, ISEED )
                  C = COS( ANGLE )*CLARND( 5, ISEED )
                  S = SIN( ANGLE )*CLARND( 5, ISEED )
                  IF( CSYM ) THEN
                     CT = C
                     ST = S
                  ELSE
                     CTEMP = CONJG( CTEMP )
                     CT = CONJG( C )
                     ST = CONJG( S )
                  END IF
                  CALL CLAROT( .FALSE., JC > K, .TRUE., IL, C, S, &
                               A( IROW-ISKEW*JC+IOFFG, JC ), ILDA, &
                               EXTRA, CTEMP )
                  CALL CLAROT( .TRUE., .TRUE., .FALSE., &
                               MIN( K, N-JC )+1, CT, ST, &
                               A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA, &
                               CTEMP, DUMMY )
!
!                    Chase EXTRA back up the matrix
!
                  ICOL = JC
                  DO JCH = JC - K, 1, -K
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CLARTG( A( JCH+1-ISKEW*( ICOL+1 )+IOFFG, &
                                  ICOL+1 ), EXTRA, REALC, S, DUMMY )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CLARTG : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     DUMMY = CLARND( 5, ISEED )
                     C = CONJG( REALC*DUMMY )
                     S = CONJG( -S*DUMMY )
                     CTEMP = A( JCH-ISKEW*( JCH+1 )+IOFFG, JCH+1 )
                     IF( CSYM ) THEN
                        CT = C
                        ST = S
                     ELSE
                        CTEMP = CONJG( CTEMP )
                        CT = CONJG( C )
                        ST = CONJG( S )
                     END IF
                     CALL CLAROT( .TRUE., .TRUE., .TRUE., K+2, C, S, &
                                  A( ( 1-ISKEW )*JCH+IOFFG, JCH ), &
                                  ILDA, CTEMP, EXTRA )
                     IROW = MAX( 1, JCH-K )
                     IL = MIN( JCH+1, K+2 )
                     EXTRA = (0.0E+0,0.0E+0)
                     CALL CLAROT( .FALSE., JCH > K, .TRUE., IL, CT, &
                                  ST, A( IROW-ISKEW*JCH+IOFFG, JCH ), &
                                  ILDA, EXTRA, CTEMP )
                     ICOL = JCH
                  ENDDO
               ENDDO
            ENDDO
!
!              If we need lower triangle, copy from upper. Note that
!              the order of copying is chosen to work for 'q' -> 'b'
!
            IF( IPACK /= IPACKG .AND. IPACK /= 3 ) THEN
               DO JC = 1, N
                  IROW = IOFFST - ISKEW*JC
                  IF( CSYM ) THEN
                     DO JR = JC, MIN( N, JC+UUB )
                        A( JR+IROW, JC ) = A( JC-ISKEW*JR+IOFFG, JR )
                     ENDDO
                  ELSE
                     DO JR = JC, MIN( N, JC+UUB )
                        A( JR+IROW, JC ) = CONJG( A( JC-ISKEW*JR+ &
                                           IOFFG, JR ) )
                     ENDDO
                  END IF
               ENDDO
               IF( IPACK == 5 ) THEN
                  DO JC = N - UUB + 1, N
                     A(N+2-JC:UUB+1, JC ) = (0.0E+0,0.0E+0)
                  ENDDO
               END IF
               IF( IPACKG == 6 ) THEN
                  IPACKG = IPACK
               ELSE
                  IPACKG = 0
               END IF
            END IF
         ELSE
!
!              Bottom-Up -- Generate Lower triangle only
!
            IF( IPACK >= 5 ) THEN
               IPACKG = 5
               IF( IPACK == 6 ) IOFFG = 1
            ELSE
               IPACKG = 2
            END IF
!
            DO J = 1, MNMIN
               A( ( 1-ISKEW )*J+IOFFG, J ) = CMPLX( D( J ) )
            ENDDO
!
            DO K = 1, UUB
               DO JC = N - 1, 1, -1
                  IL = MIN( N+1-JC, K+2 )
                  EXTRA = (0.0E+0,0.0E+0)
                  CTEMP = A( 1+( 1-ISKEW )*JC+IOFFG, JC )
                  ANGLE = TWOPI*SLARND( 1, ISEED )
                  C = COS( ANGLE )*CLARND( 5, ISEED )
                  S = SIN( ANGLE )*CLARND( 5, ISEED )
                  IF( CSYM ) THEN
                     CT = C
                     ST = S
                  ELSE
                     CTEMP = CONJG( CTEMP )
                     CT = CONJG( C )
                     ST = CONJG( S )
                  END IF
                  CALL CLAROT( .FALSE., .TRUE., N-JC > K, IL, C, S, &
                               A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA, &
                               CTEMP, EXTRA )
                  ICOL = MAX( 1, JC-K+1 )
                  CALL CLAROT( .TRUE., .FALSE., .TRUE., JC+2-ICOL, &
                               CT, ST, A( JC-ISKEW*ICOL+IOFFG, &
                               ICOL ), ILDA, DUMMY, CTEMP )
!
!                    Chase EXTRA back down the matrix
!
                  ICOL = JC
                  DO JCH = JC + K, N - 1, K
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CLARTG( A( JCH-ISKEW*ICOL+IOFFG, ICOL ), &
                                  EXTRA, REALC, S, DUMMY )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CLARTG : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     DUMMY = CLARND( 5, ISEED )
                     C = REALC*DUMMY
                     S = S*DUMMY
                     CTEMP = A( 1+( 1-ISKEW )*JCH+IOFFG, JCH )
                     IF( CSYM ) THEN
                        CT = C
                        ST = S
                     ELSE
                        CTEMP = CONJG( CTEMP )
                        CT = CONJG( C )
                        ST = CONJG( S )
                     END IF
                     CALL CLAROT( .TRUE., .TRUE., .TRUE., K+2, C, S, &
                                  A( JCH-ISKEW*ICOL+IOFFG, ICOL ), &
                                  ILDA, EXTRA, CTEMP )
                     IL = MIN( N+1-JCH, K+2 )
                     EXTRA = (0.0E+0,0.0E+0)
                     CALL CLAROT( .FALSE., .TRUE., N-JCH > K, IL, &
                                  CT, ST, A( ( 1-ISKEW )*JCH+IOFFG, &
                                  JCH ), ILDA, CTEMP, EXTRA )
                     ICOL = JCH
                  ENDDO
               ENDDO
            ENDDO
!
!              If we need upper triangle, copy from lower. Note that
!              the order of copying is chosen to work for 'b' -> 'q'
!
            IF( IPACK /= IPACKG .AND. IPACK /= 4 ) THEN
               DO JC = N, 1, -1
                  IROW = IOFFST - ISKEW*JC
                  IF( CSYM ) THEN
                     DO JR = JC, MAX( 1, JC-UUB ), -1
                        A( JR+IROW, JC ) = A( JC-ISKEW*JR+IOFFG, JR )
                     ENDDO
                  ELSE
                     DO JR = JC, MAX( 1, JC-UUB ), -1
                        A( JR+IROW, JC ) = CONJG( A( JC-ISKEW*JR+ &
                                           IOFFG, JR ) )
                     ENDDO
                  END IF
               ENDDO
               IF( IPACK == 6 ) THEN
                  DO JC = 1, UUB
                     A(1:UUB+1-JC, JC ) = (0.0E+0,0.0E+0)
                  ENDDO
               END IF
               IF( IPACKG == 5 ) THEN
                  IPACKG = IPACK
               ELSE
                  IPACKG = 0
               END IF
            END IF
         END IF
!
!           Ensure that the diagonal is real if Hermitian
!
         IF( .NOT.CSYM ) THEN
            DO JC = 1, N
               IROW = IOFFST + ( 1-ISKEW )*JC
               A( IROW, JC ) = CMPLX( REAL( A( IROW, JC ) ) )
            ENDDO
         END IF
!
      END IF
!
   ELSE
!
!        4)      Generate Banded Matrix by first
!                Rotating by random Unitary matrices,
!                then reducing the bandwidth using Householder
!                transformations.
!
!                Note: we should get here only if LDA .ge. N
!
      IF( ISYM == 1 ) THEN
!
!           Non-symmetric -- A = U D V
!
         CALL CLAGGE( MR, NC, LLB, UUB, D, A, LDA, ISEED, WORK, &
                      IINFO )
      ELSE
!
!           Symmetric -- A = U D U' or
!           Hermitian -- A = U D U*
!
         IF( CSYM ) THEN
            CALL CLAGSY( M, LLB, D, A, LDA, ISEED, WORK, IINFO )
         ELSE
            CALL CLAGHE( M, LLB, D, A, LDA, ISEED, WORK, IINFO )
         END IF
      END IF
!
      IF( IINFO /= 0 ) THEN
         INFO = 3
         RETURN
      END IF
   END IF
!
!     5)      Pack the matrix
!
   IF( IPACK /= IPACKG ) THEN
      IF( IPACK == 1 ) THEN
!
!           'U' -- Upper triangular, not packed
!
         DO J = 1, M
            A(J+1:M, J ) = (0.0E+0,0.0E+0)
         ENDDO
!
      ELSE IF( IPACK == 2 ) THEN
!
!           'L' -- Lower triangular, not packed
!
         DO J = 2, M
            A(1:J-1, J ) = (0.0E+0,0.0E+0)
         ENDDO
!
      ELSE IF( IPACK == 3 ) THEN
!
!           'C' -- Upper triangle packed Columnwise.
!
         ICOL = 1
         IROW = 0
         DO J = 1, M
            DO I = 1, J
               IROW = IROW + 1
               IF( IROW > LDA ) THEN
                  IROW = 1
                  ICOL = ICOL + 1
               END IF
               A( IROW, ICOL ) = A( I, J )
            ENDDO
         ENDDO
!
      ELSE IF( IPACK == 4 ) THEN
!
!           'R' -- Lower triangle packed Columnwise.
!
         ICOL = 1
         IROW = 0
         DO J = 1, M
            DO I = J, M
               IROW = IROW + 1
               IF( IROW > LDA ) THEN
                  IROW = 1
                  ICOL = ICOL + 1
               END IF
               A( IROW, ICOL ) = A( I, J )
            ENDDO
         ENDDO
!
      ELSE IF( IPACK >= 5 ) THEN
!
!           'B' -- The lower triangle is packed as a band matrix.
!           'Q' -- The upper triangle is packed as a band matrix.
!           'Z' -- The whole matrix is packed as a band matrix.
!
         IF( IPACK == 5 ) UUB = 0
         IF( IPACK == 6 ) LLB = 0
!
         DO J = 1, UUB
            DO I = MIN( J+LLB, M ), 1, -1
               A( I-J+UUB+1, J ) = A( I, J )
            ENDDO
         ENDDO
!
         DO J = UUB + 2, N
            DO I = J - UUB, MIN( J+LLB, M )
               A( I-J+UUB+1, J ) = A( I, J )
            ENDDO
         ENDDO
      END IF
!
!        If packed, zero out extraneous elements.
!
!        Symmetric/Triangular Packed --
!        zero out everything after A(IROW,ICOL)
!
      IF( IPACK == 3 .OR. IPACK == 4 ) THEN
         DO JC = ICOL, M
            DO JR = IROW + 1, LDA
               A( JR, JC ) = (0.0E+0,0.0E+0)
            ENDDO
            IROW = 0
         ENDDO
!
      ELSE IF( IPACK >= 5 ) THEN
!
!           Packed Band --
!              1st row is now in A( UUB+2-j, j), zero above it
!              m-th row is now in A( M+UUB-j,j), zero below it
!              last non-zero diagonal is now in A( UUB+LLB+1,j ),
!                 zero below it, too.
!
         IR1 = UUB + LLB + 2
         IR2 = UUB + M + 2
         DO JC = 1, N
            DO JR = 1, UUB + 1 - JC
               A( JR, JC ) = (0.0E+0,0.0E+0)
            ENDDO
            DO JR = MAX( 1, MIN( IR1, IR2-JC ) ), LDA
               A( JR, JC ) = (0.0E+0,0.0E+0)
            ENDDO
         ENDDO
      END IF
   END IF
!
   RETURN
!
!     End of CLATMS
!
END

