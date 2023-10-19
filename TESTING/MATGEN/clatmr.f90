!> \brief \b CLATMR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLATMR( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX,
!                          RSIGN, GRADE, DL, MODEL, CONDL, DR, MODER,
!                          CONDR, PIVTNG, IPIVOT, KL, KU, SPARSE, ANORM,
!                          PACK, A, LDA, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIST, GRADE, PACK, PIVTNG, RSIGN, SYM
!       INTEGER            INFO, KL, KU, LDA, M, MODE, MODEL, MODER, N
!       REAL               ANORM, COND, CONDL, CONDR, SPARSE
!       COMPLEX            DMAX
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIVOT( * ), ISEED( 4 ), IWORK( * )
!       COMPLEX            A( LDA, * ), D( * ), DL( * ), DR( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CLATMR generates random matrices of various types for testing
!>    LAPACK programs.
!>
!>    CLATMR operates by applying the following sequence of
!>    operations:
!>
!>      Generate a matrix A with random entries of distribution DIST
!>         which is symmetric if SYM='S', Hermitian if SYM='H', and
!>         nonsymmetric if SYM='N'.
!>
!>      Set the diagonal to D, where D may be input or
!>         computed according to MODE, COND, DMAX and RSIGN
!>         as described below.
!>
!>      Grade the matrix, if desired, from the left and/or right
!>         as specified by GRADE. The inputs DL, MODEL, CONDL, DR,
!>         MODER and CONDR also determine the grading as described
!>         below.
!>
!>      Permute, if desired, the rows and/or columns as specified by
!>         PIVTNG and IPIVOT.
!>
!>      Set random entries to zero, if desired, to get a random sparse
!>         matrix as specified by SPARSE.
!>
!>      Make A a band matrix, if desired, by zeroing out the matrix
!>         outside a band of lower bandwidth KL and upper bandwidth KU.
!>
!>      Scale A, if desired, to have maximum entry ANORM.
!>
!>      Pack the matrix if desired. Options specified by PACK are:
!>         no packing
!>         zero out upper half (if symmetric or Hermitian)
!>         zero out lower half (if symmetric or Hermitian)
!>         store the upper half columnwise (if symmetric or Hermitian
!>             or square upper triangular)
!>         store the lower half columnwise (if symmetric or Hermitian
!>             or square lower triangular)
!>             same as upper half rowwise if symmetric
!>             same as conjugate upper half rowwise if Hermitian
!>         store the lower triangle in banded format
!>             (if symmetric or Hermitian)
!>         store the upper triangle in banded format
!>             (if symmetric or Hermitian)
!>         store the entire matrix in banded format
!>
!>    Note: If two calls to CLATMR differ only in the PACK parameter,
!>          they will generate mathematically equivalent matrices.
!>
!>          If two calls to CLATMR both have full bandwidth (KL = M-1
!>          and KU = N-1), and differ only in the PIVTNG and PACK
!>          parameters, then the matrices generated will differ only
!>          in the order of the rows and/or columns, and otherwise
!>          contain the same data. This consistency cannot be and
!>          is not maintained with less than full bandwidth.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           Number of rows of A. Not modified.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           Number of columns of A. Not modified.
!> \endverbatim
!>
!> \param[in] DIST
!> \verbatim
!>          DIST is CHARACTER*1
!>           On entry, DIST specifies the type of distribution to be used
!>           to generate a random matrix .
!>           'U' => real and imaginary parts are independent
!>                  UNIFORM( 0, 1 )  ( 'U' for uniform )
!>           'S' => real and imaginary parts are independent
!>                  UNIFORM( -1, 1 ) ( 'S' for symmetric )
!>           'N' => real and imaginary parts are independent
!>                  NORMAL( 0, 1 )   ( 'N' for normal )
!>           'D' => uniform on interior of unit disk ( 'D' for disk )
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>           On entry ISEED specifies the seed of the random number
!>           generator. They should lie between 0 and 4095 inclusive,
!>           and ISEED(4) should be odd. The random number generator
!>           uses a linear congruential sequence limited to small
!>           integers, and so should produce machine independent
!>           random numbers. The values of ISEED are changed on
!>           exit, and can be used in the next call to CLATMR
!>           to continue the same random number sequence.
!>           Changed on exit.
!> \endverbatim
!>
!> \param[in] SYM
!> \verbatim
!>          SYM is CHARACTER*1
!>           If SYM='S', generated matrix is symmetric.
!>           If SYM='H', generated matrix is Hermitian.
!>           If SYM='N', generated matrix is nonsymmetric.
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is COMPLEX array, dimension (min(M,N))
!>           On entry this array specifies the diagonal entries
!>           of the diagonal of A.  D may either be specified
!>           on entry, or set according to MODE and COND as described
!>           below. If the matrix is Hermitian, the real part of D
!>           will be taken. May be changed on exit if MODE is nonzero.
!> \endverbatim
!>
!> \param[in] MODE
!> \verbatim
!>          MODE is INTEGER
!>           On entry describes how D is to be used:
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
!>           Not modified.
!> \endverbatim
!>
!> \param[in] COND
!> \verbatim
!>          COND is REAL
!>           On entry, used as described under MODE above.
!>           If used, it must be >= 1. Not modified.
!> \endverbatim
!>
!> \param[in] DMAX
!> \verbatim
!>          DMAX is COMPLEX
!>           If MODE neither -6, 0 nor 6, the diagonal is scaled by
!>           DMAX / max(abs(D(i))), so that maximum absolute entry
!>           of diagonal is abs(DMAX). If DMAX is complex (or zero),
!>           diagonal will be scaled by a complex number (or zero).
!> \endverbatim
!>
!> \param[in] RSIGN
!> \verbatim
!>          RSIGN is CHARACTER*1
!>           If MODE neither -6, 0 nor 6, specifies sign of diagonal
!>           as follows:
!>           'T' => diagonal entries are multiplied by a random complex
!>                  number uniformly distributed with absolute value 1
!>           'F' => diagonal unchanged
!>           Not modified.
!> \endverbatim
!>
!> \param[in] GRADE
!> \verbatim
!>          GRADE is CHARACTER*1
!>           Specifies grading of matrix as follows:
!>           'N'  => no grading
!>           'L'  => matrix premultiplied by diag( DL )
!>                   (only if matrix nonsymmetric)
!>           'R'  => matrix postmultiplied by diag( DR )
!>                   (only if matrix nonsymmetric)
!>           'B'  => matrix premultiplied by diag( DL ) and
!>                         postmultiplied by diag( DR )
!>                   (only if matrix nonsymmetric)
!>           'H'  => matrix premultiplied by diag( DL ) and
!>                         postmultiplied by diag( CONJG(DL) )
!>                   (only if matrix Hermitian or nonsymmetric)
!>           'S'  => matrix premultiplied by diag( DL ) and
!>                         postmultiplied by diag( DL )
!>                   (only if matrix symmetric or nonsymmetric)
!>           'E'  => matrix premultiplied by diag( DL ) and
!>                         postmultiplied by inv( diag( DL ) )
!>                         ( 'S' for similarity )
!>                   (only if matrix nonsymmetric)
!>                   Note: if GRADE='S', then M must equal N.
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] DL
!> \verbatim
!>          DL is COMPLEX array, dimension (M)
!>           If MODEL=0, then on entry this array specifies the diagonal
!>           entries of a diagonal matrix used as described under GRADE
!>           above. If MODEL is not zero, then DL will be set according
!>           to MODEL and CONDL, analogous to the way D is set according
!>           to MODE and COND (except there is no DMAX parameter for DL).
!>           If GRADE='E', then DL cannot have zero entries.
!>           Not referenced if GRADE = 'N' or 'R'. Changed on exit.
!> \endverbatim
!>
!> \param[in] MODEL
!> \verbatim
!>          MODEL is INTEGER
!>           This specifies how the diagonal array DL is to be computed,
!>           just as MODE specifies how D is to be computed.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] CONDL
!> \verbatim
!>          CONDL is REAL
!>           When MODEL is not zero, this specifies the condition number
!>           of the computed DL.  Not modified.
!> \endverbatim
!>
!> \param[in,out] DR
!> \verbatim
!>          DR is COMPLEX array, dimension (N)
!>           If MODER=0, then on entry this array specifies the diagonal
!>           entries of a diagonal matrix used as described under GRADE
!>           above. If MODER is not zero, then DR will be set according
!>           to MODER and CONDR, analogous to the way D is set according
!>           to MODE and COND (except there is no DMAX parameter for DR).
!>           Not referenced if GRADE = 'N', 'L', 'H' or 'S'.
!>           Changed on exit.
!> \endverbatim
!>
!> \param[in] MODER
!> \verbatim
!>          MODER is INTEGER
!>           This specifies how the diagonal array DR is to be computed,
!>           just as MODE specifies how D is to be computed.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] CONDR
!> \verbatim
!>          CONDR is REAL
!>           When MODER is not zero, this specifies the condition number
!>           of the computed DR.  Not modified.
!> \endverbatim
!>
!> \param[in] PIVTNG
!> \verbatim
!>          PIVTNG is CHARACTER*1
!>           On entry specifies pivoting permutations as follows:
!>           'N' or ' ' => none.
!>           'L' => left or row pivoting (matrix must be nonsymmetric).
!>           'R' => right or column pivoting (matrix must be
!>                  nonsymmetric).
!>           'B' or 'F' => both or full pivoting, i.e., on both sides.
!>                         In this case, M must equal N
!>
!>           If two calls to CLATMR both have full bandwidth (KL = M-1
!>           and KU = N-1), and differ only in the PIVTNG and PACK
!>           parameters, then the matrices generated will differ only
!>           in the order of the rows and/or columns, and otherwise
!>           contain the same data. This consistency cannot be
!>           maintained with less than full bandwidth.
!> \endverbatim
!>
!> \param[in] IPIVOT
!> \verbatim
!>          IPIVOT is INTEGER array, dimension (N or M)
!>           This array specifies the permutation used.  After the
!>           basic matrix is generated, the rows, columns, or both
!>           are permuted.   If, say, row pivoting is selected, CLATMR
!>           starts with the *last* row and interchanges the M-th and
!>           IPIVOT(M)-th rows, then moves to the next-to-last row,
!>           interchanging the (M-1)-th and the IPIVOT(M-1)-th rows,
!>           and so on.  In terms of "2-cycles", the permutation is
!>           (1 IPIVOT(1)) (2 IPIVOT(2)) ... (M IPIVOT(M))
!>           where the rightmost cycle is applied first.  This is the
!>           *inverse* of the effect of pivoting in LINPACK.  The idea
!>           is that factoring (with pivoting) an identity matrix
!>           which has been inverse-pivoted in this way should
!>           result in a pivot vector identical to IPIVOT.
!>           Not referenced if PIVTNG = 'N'. Not modified.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>           On entry specifies the lower bandwidth of the  matrix. For
!>           example, KL=0 implies upper triangular, KL=1 implies upper
!>           Hessenberg, and KL at least M-1 implies the matrix is not
!>           banded. Must equal KU if matrix is symmetric or Hermitian.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>           On entry specifies the upper bandwidth of the  matrix. For
!>           example, KU=0 implies lower triangular, KU=1 implies lower
!>           Hessenberg, and KU at least N-1 implies the matrix is not
!>           banded. Must equal KL if matrix is symmetric or Hermitian.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] SPARSE
!> \verbatim
!>          SPARSE is REAL
!>           On entry specifies the sparsity of the matrix if a sparse
!>           matrix is to be generated. SPARSE should lie between
!>           0 and 1. To generate a sparse matrix, for each matrix entry
!>           a uniform ( 0, 1 ) random number x is generated and
!>           compared to SPARSE; if x is larger the matrix entry
!>           is unchanged and if x is smaller the entry is set
!>           to zero. Thus on the average a fraction SPARSE of the
!>           entries will be set to zero.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] ANORM
!> \verbatim
!>          ANORM is REAL
!>           On entry specifies maximum entry of output matrix
!>           (output matrix will by multiplied by a constant so that
!>           its largest absolute entry equal ANORM)
!>           if ANORM is nonnegative. If ANORM is negative no scaling
!>           is done. Not modified.
!> \endverbatim
!>
!> \param[in] PACK
!> \verbatim
!>          PACK is CHARACTER*1
!>           On entry specifies packing of matrix as follows:
!>           'N' => no packing
!>           'U' => zero out all subdiagonal entries
!>                  (if symmetric or Hermitian)
!>           'L' => zero out all superdiagonal entries
!>                  (if symmetric or Hermitian)
!>           'C' => store the upper triangle columnwise
!>                  (only if matrix symmetric or Hermitian or
!>                   square upper triangular)
!>           'R' => store the lower triangle columnwise
!>                  (only if matrix symmetric or Hermitian or
!>                   square lower triangular)
!>                  (same as upper half rowwise if symmetric)
!>                  (same as conjugate upper half rowwise if Hermitian)
!>           'B' => store the lower triangle in band storage scheme
!>                  (only if matrix symmetric or Hermitian)
!>           'Q' => store the upper triangle in band storage scheme
!>                  (only if matrix symmetric or Hermitian)
!>           'Z' => store the entire matrix in band storage scheme
!>                      (pivoting can be provided for by using this
!>                      option to store A in the trailing rows of
!>                      the allocated storage)
!>
!>           Using these options, the various LAPACK packed and banded
!>           storage schemes can be obtained:
!>           GB               - use 'Z'
!>           PB, HB or TB     - use 'B' or 'Q'
!>           PP, HP or TP     - use 'C' or 'R'
!>
!>           If two calls to CLATMR differ only in the PACK parameter,
!>           they will generate mathematically equivalent matrices.
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>           On exit A is the desired test matrix. Only those
!>           entries of A which are significant on output
!>           will be referenced (even if A is in packed or band
!>           storage format). The 'unoccupied corners' of A in
!>           band format will be zeroed out.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           on entry LDA specifies the first dimension of A as
!>           declared in the calling program.
!>           If PACK='N', 'U' or 'L', LDA must be at least max ( 1, M ).
!>           If PACK='C' or 'R', LDA must be at least 1.
!>           If PACK='B', or 'Q', LDA must be MIN ( KU+1, N )
!>           If PACK='Z', LDA must be at least KUU+KLL+1, where
!>           KUU = MIN ( KU, N-1 ) and KLL = MIN ( KL, M-1 )
!>           Not modified.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N or M)
!>           Workspace. Not referenced if PIVTNG = 'N'. Changed on exit.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           Error parameter on exit:
!>             0 => normal return
!>            -1 => M negative or unequal to N and SYM='S' or 'H'
!>            -2 => N negative
!>            -3 => DIST illegal string
!>            -5 => SYM illegal string
!>            -7 => MODE not in range -6 to 6
!>            -8 => COND less than 1.0, and MODE neither -6, 0 nor 6
!>           -10 => MODE neither -6, 0 nor 6 and RSIGN illegal string
!>           -11 => GRADE illegal string, or GRADE='E' and
!>                  M not equal to N, or GRADE='L', 'R', 'B', 'S' or 'E'
!>                  and SYM = 'H', or GRADE='L', 'R', 'B', 'H' or 'E'
!>                  and SYM = 'S'
!>           -12 => GRADE = 'E' and DL contains zero
!>           -13 => MODEL not in range -6 to 6 and GRADE= 'L', 'B', 'H',
!>                  'S' or 'E'
!>           -14 => CONDL less than 1.0, GRADE='L', 'B', 'H', 'S' or 'E',
!>                  and MODEL neither -6, 0 nor 6
!>           -16 => MODER not in range -6 to 6 and GRADE= 'R' or 'B'
!>           -17 => CONDR less than 1.0, GRADE='R' or 'B', and
!>                  MODER neither -6, 0 nor 6
!>           -18 => PIVTNG illegal string, or PIVTNG='B' or 'F' and
!>                  M not equal to N, or PIVTNG='L' or 'R' and SYM='S'
!>                  or 'H'
!>           -19 => IPIVOT contains out of range number and
!>                  PIVTNG not equal to 'N'
!>           -20 => KL negative
!>           -21 => KU negative, or SYM='S' or 'H' and KU not equal to KL
!>           -22 => SPARSE not in range 0. to 1.
!>           -24 => PACK illegal string, or PACK='U', 'L', 'B' or 'Q'
!>                  and SYM='N', or PACK='C' and SYM='N' and either KL
!>                  not equal to 0 or N not equal to M, or PACK='R' and
!>                  SYM='N', and either KU not equal to 0 or N not equal
!>                  to M
!>           -26 => LDA too small
!>             1 => Error return from CLATM1 (computing D)
!>             2 => Cannot scale diagonal to DMAX (max. entry is 0)
!>             3 => Error return from CLATM1 (computing DL)
!>             4 => Error return from CLATM1 (computing DR)
!>             5 => ANORM is positive, but matrix constructed prior to
!>                  attempting to scale it to have norm ANORM, is zero
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
   SUBROUTINE CLATMR( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX, &
                      RSIGN, GRADE, DL, MODEL, CONDL, DR, MODER, &
                      CONDR, PIVTNG, IPIVOT, KL, KU, SPARSE, ANORM, &
                      PACK, A, LDA, IWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          DIST, GRADE, PACK, PIVTNG, RSIGN, SYM
   INTEGER            INFO, KL, KU, LDA, M, MODE, MODEL, MODER, N
   REAL               ANORM, COND, CONDL, CONDR, SPARSE
   COMPLEX            DMAX
!     ..
!     .. Array Arguments ..
   INTEGER            IPIVOT( * ), ISEED( 4 ), IWORK( * )
   COMPLEX            A( LDA, * ), D( * ), DL( * ), DR( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO
   PARAMETER          ( ZERO = 0.0E0 )
   REAL               ONE
   PARAMETER          ( ONE = 1.0E0 )
   COMPLEX            CONE
   PARAMETER          ( CONE = ( 1.0E0, 0.0E0 ) )
   COMPLEX            CZERO
   PARAMETER          ( CZERO = ( 0.0E0, 0.0E0 ) )
!     ..
!     .. Local Scalars ..
   LOGICAL            BADPVT, DZERO, FULBND
   INTEGER            I, IDIST, IGRADE, IISUB, IPACK, IPVTNG, IRSIGN, &
                      ISUB, ISYM, J, JJSUB, JSUB, K, KLL, KUU, MNMIN, &
                      MNSUB, MXSUB, NPVTS
   REAL               ONORM, TEMP
   COMPLEX            CALPHA, CTEMP
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   REAL               TEMPA( 1 )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               CLANGB, CLANGE, CLANSB, CLANSP, CLANSY
   COMPLEX            CLATM2, CLATM3
   EXTERNAL           LSAME, CLANGB, CLANGE, CLANSB, CLANSP, CLANSY, &
                      CLATM2, CLATM3
!     ..
!     .. External Subroutines ..
   EXTERNAL           CLATM1, CSSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, CONJG, MAX, MIN, MOD, REAL
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
   IF( M == 0 .OR. N == 0 ) &
      RETURN
!
!     Decode DIST
!
   IF( LSAME( DIST, 'U' ) ) THEN
      IDIST = 1
   ELSE IF( LSAME( DIST, 'S' ) ) THEN
      IDIST = 2
   ELSE IF( LSAME( DIST, 'N' ) ) THEN
      IDIST = 3
   ELSE IF( LSAME( DIST, 'D' ) ) THEN
      IDIST = 4
   ELSE
      IDIST = -1
   END IF
!
!     Decode SYM
!
   IF( LSAME( SYM, 'H' ) ) THEN
      ISYM = 0
   ELSE IF( LSAME( SYM, 'N' ) ) THEN
      ISYM = 1
   ELSE IF( LSAME( SYM, 'S' ) ) THEN
      ISYM = 2
   ELSE
      ISYM = -1
   END IF
!
!     Decode RSIGN
!
   IF( LSAME( RSIGN, 'F' ) ) THEN
      IRSIGN = 0
   ELSE IF( LSAME( RSIGN, 'T' ) ) THEN
      IRSIGN = 1
   ELSE
      IRSIGN = -1
   END IF
!
!     Decode PIVTNG
!
   IF( LSAME( PIVTNG, 'N' ) ) THEN
      IPVTNG = 0
   ELSE IF( LSAME( PIVTNG, ' ' ) ) THEN
      IPVTNG = 0
   ELSE IF( LSAME( PIVTNG, 'L' ) ) THEN
      IPVTNG = 1
      NPVTS = M
   ELSE IF( LSAME( PIVTNG, 'R' ) ) THEN
      IPVTNG = 2
      NPVTS = N
   ELSE IF( LSAME( PIVTNG, 'B' ) ) THEN
      IPVTNG = 3
      NPVTS = MIN( N, M )
   ELSE IF( LSAME( PIVTNG, 'F' ) ) THEN
      IPVTNG = 3
      NPVTS = MIN( N, M )
   ELSE
      IPVTNG = -1
   END IF
!
!     Decode GRADE
!
   IF( LSAME( GRADE, 'N' ) ) THEN
      IGRADE = 0
   ELSE IF( LSAME( GRADE, 'L' ) ) THEN
      IGRADE = 1
   ELSE IF( LSAME( GRADE, 'R' ) ) THEN
      IGRADE = 2
   ELSE IF( LSAME( GRADE, 'B' ) ) THEN
      IGRADE = 3
   ELSE IF( LSAME( GRADE, 'E' ) ) THEN
      IGRADE = 4
   ELSE IF( LSAME( GRADE, 'H' ) ) THEN
      IGRADE = 5
   ELSE IF( LSAME( GRADE, 'S' ) ) THEN
      IGRADE = 6
   ELSE
      IGRADE = -1
   END IF
!
!     Decode PACK
!
   IF( LSAME( PACK, 'N' ) ) THEN
      IPACK = 0
   ELSE IF( LSAME( PACK, 'U' ) ) THEN
      IPACK = 1
   ELSE IF( LSAME( PACK, 'L' ) ) THEN
      IPACK = 2
   ELSE IF( LSAME( PACK, 'C' ) ) THEN
      IPACK = 3
   ELSE IF( LSAME( PACK, 'R' ) ) THEN
      IPACK = 4
   ELSE IF( LSAME( PACK, 'B' ) ) THEN
      IPACK = 5
   ELSE IF( LSAME( PACK, 'Q' ) ) THEN
      IPACK = 6
   ELSE IF( LSAME( PACK, 'Z' ) ) THEN
      IPACK = 7
   ELSE
      IPACK = -1
   END IF
!
!     Set certain internal parameters
!
   MNMIN = MIN( M, N )
   KLL = MIN( KL, M-1 )
   KUU = MIN( KU, N-1 )
!
!     If inv(DL) is used, check to see if DL has a zero entry.
!
   DZERO = .FALSE.
   IF( IGRADE == 4 .AND. MODEL == 0 ) THEN
      DO I = 1, M
         IF( DL( I ) == CZERO ) &
            DZERO = .TRUE.
      ENDDO
   END IF
!
!     Check values in IPIVOT
!
   BADPVT = .FALSE.
   IF( IPVTNG > 0 ) THEN
      DO J = 1, NPVTS
         IF( IPIVOT( J ) <= 0 .OR. IPIVOT( J ) > NPVTS ) &
            BADPVT = .TRUE.
      ENDDO
   END IF
!
!     Set INFO if an error
!
   IF( M < 0 ) THEN
      INFO = -1
   ELSE IF( M /= N .AND. ( ISYM == 0 .OR. ISYM == 2 ) ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( IDIST == -1 ) THEN
      INFO = -3
   ELSE IF( ISYM == -1 ) THEN
      INFO = -5
   ELSE IF( MODE < -6 .OR. MODE > 6 ) THEN
      INFO = -7
   ELSE IF( ( MODE /= -6 .AND. MODE /= 0 .AND. MODE /= 6 ) .AND. &
            COND < ONE ) THEN
      INFO = -8
   ELSE IF( ( MODE /= -6 .AND. MODE /= 0 .AND. MODE /= 6 ) .AND. &
            IRSIGN == -1 ) THEN
      INFO = -10
   ELSE IF( IGRADE == -1 .OR. ( IGRADE == 4 .AND. M /= N ) .OR. &
            ( ( IGRADE == 1 .OR. IGRADE == 2 .OR. IGRADE == 3 .OR. &
            IGRADE == 4 .OR. IGRADE == 6 ) .AND. ISYM == 0 ) .OR. &
            ( ( IGRADE == 1 .OR. IGRADE == 2 .OR. IGRADE == 3 .OR. &
            IGRADE == 4 .OR. IGRADE == 5 ) .AND. ISYM == 2 ) ) THEN
      INFO = -11
   ELSE IF( IGRADE == 4 .AND. DZERO ) THEN
      INFO = -12
   ELSE IF( ( IGRADE == 1 .OR. IGRADE == 3 .OR. IGRADE == 4 .OR. &
            IGRADE == 5 .OR. IGRADE == 6 ) .AND. &
            ( MODEL < -6 .OR. MODEL > 6 ) ) THEN
      INFO = -13
   ELSE IF( ( IGRADE == 1 .OR. IGRADE == 3 .OR. IGRADE == 4 .OR. &
            IGRADE == 5 .OR. IGRADE == 6 ) .AND. &
            ( MODEL /= -6 .AND. MODEL /= 0 .AND. MODEL /= 6 ) .AND. &
            CONDL < ONE ) THEN
      INFO = -14
   ELSE IF( ( IGRADE == 2 .OR. IGRADE == 3 ) .AND. &
            ( MODER < -6 .OR. MODER > 6 ) ) THEN
      INFO = -16
   ELSE IF( ( IGRADE == 2 .OR. IGRADE == 3 ) .AND. &
            ( MODER /= -6 .AND. MODER /= 0 .AND. MODER /= 6 ) .AND. &
            CONDR < ONE ) THEN
      INFO = -17
   ELSE IF( IPVTNG == -1 .OR. ( IPVTNG == 3 .AND. M /= N ) .OR. &
            ( ( IPVTNG == 1 .OR. IPVTNG == 2 ) .AND. ( ISYM == 0 .OR. &
            ISYM == 2 ) ) ) THEN
      INFO = -18
   ELSE IF( IPVTNG /= 0 .AND. BADPVT ) THEN
      INFO = -19
   ELSE IF( KL < 0 ) THEN
      INFO = -20
   ELSE IF( KU < 0 .OR. ( ( ISYM == 0 .OR. ISYM == 2 ) .AND. KL /= &
            KU ) ) THEN
      INFO = -21
   ELSE IF( SPARSE < ZERO .OR. SPARSE > ONE ) THEN
      INFO = -22
   ELSE IF( IPACK == -1 .OR. ( ( IPACK == 1 .OR. IPACK == 2 .OR. &
            IPACK == 5 .OR. IPACK == 6 ) .AND. ISYM == 1 ) .OR. &
            ( IPACK == 3 .AND. ISYM == 1 .AND. ( KL /= 0 .OR. M /= &
            N ) ) .OR. ( IPACK == 4 .AND. ISYM == 1 .AND. ( KU /= &
            0 .OR. M /= N ) ) ) THEN
      INFO = -24
   ELSE IF( ( ( IPACK == 0 .OR. IPACK == 1 .OR. IPACK == 2 ) .AND. &
            LDA < MAX( 1, M ) ) .OR. ( ( IPACK == 3 .OR. IPACK == &
            4 ) .AND. LDA < 1 ) .OR. ( ( IPACK == 5 .OR. IPACK == &
            6 ) .AND. LDA < KUU+1 ) .OR. &
            ( IPACK == 7 .AND. LDA < KLL+KUU+1 ) ) THEN
      INFO = -26
   END IF
!
   IF( INFO /= 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'CLATMR', -INFO )
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
!     Decide if we can pivot consistently
!
   FULBND = .FALSE.
   IF( KUU == N-1 .AND. KLL == M-1 ) &
      FULBND = .TRUE.
!
!     Initialize random number generator
!
   DO I = 1, 4
      ISEED( I ) = MOD( ABS( ISEED( I ) ), 4096 )
   ENDDO
!
   ISEED( 4 ) = 2*( ISEED( 4 ) / 2 ) + 1
!
!     2)      Set up D, DL, and DR, if indicated.
!
!             Compute D according to COND and MODE
!
   CALL CLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, INFO )
   IF( INFO /= 0 ) THEN
      INFO = 1
      RETURN
   END IF
   IF( MODE /= 0 .AND. MODE /= -6 .AND. MODE /= 6 ) THEN
!
!        Scale by DMAX
!
      TEMP = ABS( D( 1 ) )
      DO I = 2, MNMIN
         TEMP = MAX( TEMP, ABS( D( I ) ) )
      ENDDO
      IF( TEMP == ZERO .AND. DMAX /= CZERO ) THEN
         INFO = 2
         RETURN
      END IF
      IF( TEMP /= ZERO ) THEN
         CALPHA = DMAX / TEMP
      ELSE
         CALPHA = CONE
      END IF
      DO I = 1, MNMIN
         D( I ) = CALPHA*D( I )
      ENDDO
!
   END IF
!
!     If matrix Hermitian, make D real
!
   IF( ISYM == 0 ) THEN
      DO I = 1, MNMIN
         D( I ) = REAL( D( I ) )
      ENDDO
   END IF
!
!     Compute DL if grading set
!
   IF( IGRADE == 1 .OR. IGRADE == 3 .OR. IGRADE == 4 .OR. IGRADE == &
       5 .OR. IGRADE == 6 ) THEN
      CALL CLATM1( MODEL, CONDL, 0, IDIST, ISEED, DL, M, INFO )
      IF( INFO /= 0 ) THEN
         INFO = 3
         RETURN
      END IF
   END IF
!
!     Compute DR if grading set
!
   IF( IGRADE == 2 .OR. IGRADE == 3 ) THEN
      CALL CLATM1( MODER, CONDR, 0, IDIST, ISEED, DR, N, INFO )
      IF( INFO /= 0 ) THEN
         INFO = 4
         RETURN
      END IF
   END IF
!
!     3)     Generate IWORK if pivoting
!
   IF( IPVTNG > 0 ) THEN
      DO I = 1, NPVTS
         IWORK( I ) = I
      ENDDO
      IF( FULBND ) THEN
         DO I = 1, NPVTS
            K = IPIVOT( I )
            J = IWORK( I )
            IWORK( I ) = IWORK( K )
            IWORK( K ) = J
         ENDDO
      ELSE
         DO I = NPVTS, 1, -1
            K = IPIVOT( I )
            J = IWORK( I )
            IWORK( I ) = IWORK( K )
            IWORK( K ) = J
         ENDDO
      END IF
   END IF
!
!     4)      Generate matrices for each kind of PACKing
!             Always sweep matrix columnwise (if symmetric, upper
!             half only) so that matrix generated does not depend
!             on PACK
!
   IF( FULBND ) THEN
!
!        Use CLATM3 so matrices generated with differing PIVOTing only
!        differ only in the order of their rows and/or columns.
!
      IF( IPACK == 0 ) THEN
         IF( ISYM == 0 ) THEN
            DO J = 1, N
               DO I = 1, J
                  CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, &
                          IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, &
                          IWORK, SPARSE )
                  A( ISUB, JSUB ) = CTEMP
                  A( JSUB, ISUB ) = CONJG( CTEMP )
                  ENDDO
               ENDDO
         ELSE IF( ISYM == 1 ) THEN
            DO J = 1, N
               DO I = 1, M
                  CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, &
                          IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, &
                          IWORK, SPARSE )
                  A( ISUB, JSUB ) = CTEMP
                  ENDDO
               ENDDO
         ELSE IF( ISYM == 2 ) THEN
            DO J = 1, N
               DO I = 1, J
                  CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, &
                          IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, &
                          IWORK, SPARSE )
                  A( ISUB, JSUB ) = CTEMP
                  A( JSUB, ISUB ) = CTEMP
                  ENDDO
               ENDDO
         END IF
!
      ELSE IF( IPACK == 1 ) THEN
!
         DO J = 1, N
            DO I = 1, J
               CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, &
                       ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, &
                       SPARSE )
               MNSUB = MIN( ISUB, JSUB )
               MXSUB = MAX( ISUB, JSUB )
               IF( MXSUB == ISUB .AND. ISYM == 0 ) THEN
                  A( MNSUB, MXSUB ) = CONJG( CTEMP )
               ELSE
                  A( MNSUB, MXSUB ) = CTEMP
               END IF
               IF( MNSUB /= MXSUB ) &
                  A( MXSUB, MNSUB ) = CZERO
               ENDDO
            ENDDO
!
      ELSE IF( IPACK == 2 ) THEN
!
         DO J = 1, N
            DO I = 1, J
               CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, &
                       ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, &
                       SPARSE )
               MNSUB = MIN( ISUB, JSUB )
               MXSUB = MAX( ISUB, JSUB )
               IF( MXSUB == JSUB .AND. ISYM == 0 ) THEN
                  A( MXSUB, MNSUB ) = CONJG( CTEMP )
               ELSE
                  A( MXSUB, MNSUB ) = CTEMP
               END IF
               IF( MNSUB /= MXSUB ) &
                  A( MNSUB, MXSUB ) = CZERO
               ENDDO
            ENDDO
!
      ELSE IF( IPACK == 3 ) THEN
!
         DO J = 1, N
            DO I = 1, J
               CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, &
                       ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, &
                       SPARSE )
!
!                 Compute K = location of (ISUB,JSUB) entry in packed
!                 array
!
               MNSUB = MIN( ISUB, JSUB )
               MXSUB = MAX( ISUB, JSUB )
               K = MXSUB*( MXSUB-1 ) / 2 + MNSUB
!
!                 Convert K to (IISUB,JJSUB) location
!
               JJSUB = ( K-1 ) / LDA + 1
               IISUB = K - LDA*( JJSUB-1 )
!
               IF( MXSUB == ISUB .AND. ISYM == 0 ) THEN
                  A( IISUB, JJSUB ) = CONJG( CTEMP )
               ELSE
                  A( IISUB, JJSUB ) = CTEMP
               END IF
               ENDDO
            ENDDO
!
      ELSE IF( IPACK == 4 ) THEN
!
         DO J = 1, N
            DO I = 1, J
               CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, &
                       ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, &
                       SPARSE )
!
!                 Compute K = location of (I,J) entry in packed array
!
               MNSUB = MIN( ISUB, JSUB )
               MXSUB = MAX( ISUB, JSUB )
               IF( MNSUB == 1 ) THEN
                  K = MXSUB
               ELSE
                  K = N*( N+1 ) / 2 - ( N-MNSUB+1 )*( N-MNSUB+2 ) / &
                      2 + MXSUB - MNSUB + 1
               END IF
!
!                 Convert K to (IISUB,JJSUB) location
!
               JJSUB = ( K-1 ) / LDA + 1
               IISUB = K - LDA*( JJSUB-1 )
!
               IF( MXSUB == JSUB .AND. ISYM == 0 ) THEN
                  A( IISUB, JJSUB ) = CONJG( CTEMP )
               ELSE
                  A( IISUB, JJSUB ) = CTEMP
               END IF
               ENDDO
            ENDDO
!
      ELSE IF( IPACK == 5 ) THEN
!
         DO J = 1, N
            DO I = J - KUU, J
               IF( I < 1 ) THEN
                  A( J-I+1, I+N ) = CZERO
               ELSE
                  CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, &
                          IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, &
                          IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  IF( MXSUB == JSUB .AND. ISYM == 0 ) THEN
                     A( MXSUB-MNSUB+1, MNSUB ) = CONJG( CTEMP )
                  ELSE
                     A( MXSUB-MNSUB+1, MNSUB ) = CTEMP
                  END IF
               END IF
               ENDDO
            ENDDO
!
      ELSE IF( IPACK == 6 ) THEN
!
         DO J = 1, N
            DO I = J - KUU, J
               CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, &
                       ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, &
                       SPARSE )
               MNSUB = MIN( ISUB, JSUB )
               MXSUB = MAX( ISUB, JSUB )
               IF( MXSUB == ISUB .AND. ISYM == 0 ) THEN
                  A( MNSUB-MXSUB+KUU+1, MXSUB ) = CONJG( CTEMP )
               ELSE
                  A( MNSUB-MXSUB+KUU+1, MXSUB ) = CTEMP
               END IF
               ENDDO
            ENDDO
!
      ELSE IF( IPACK == 7 ) THEN
!
         IF( ISYM /= 1 ) THEN
            DO J = 1, N
               DO I = J - KUU, J
                  CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, &
                          IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, &
                          IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  IF( I < 1 ) &
                     A( J-I+1+KUU, I+N ) = CZERO
                  IF( MXSUB == ISUB .AND. ISYM == 0 ) THEN
                     A( MNSUB-MXSUB+KUU+1, MXSUB ) = CONJG( CTEMP )
                  ELSE
                     A( MNSUB-MXSUB+KUU+1, MXSUB ) = CTEMP
                  END IF
                  IF( I >= 1 .AND. MNSUB /= MXSUB ) THEN
                     IF( MNSUB == ISUB .AND. ISYM == 0 ) THEN
                        A( MXSUB-MNSUB+1+KUU, &
                           MNSUB ) = CONJG( CTEMP )
                     ELSE
                        A( MXSUB-MNSUB+1+KUU, MNSUB ) = CTEMP
                     END IF
                  END IF
                  ENDDO
               ENDDO
         ELSE IF( ISYM == 1 ) THEN
            DO J = 1, N
               DO I = J - KUU, J + KLL
                  CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, &
                          IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, &
                          IWORK, SPARSE )
                  A( ISUB-JSUB+KUU+1, JSUB ) = CTEMP
                  ENDDO
               ENDDO
         END IF
!
      END IF
!
   ELSE
!
!        Use CLATM2
!
      IF( IPACK == 0 ) THEN
         IF( ISYM == 0 ) THEN
            DO J = 1, N
               DO I = 1, J
                  A( I, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, &
                              ISEED, D, IGRADE, DL, DR, IPVTNG, &
                              IWORK, SPARSE )
                  A( J, I ) = CONJG( A( I, J ) )
                  ENDDO
               ENDDO
         ELSE IF( ISYM == 1 ) THEN
            DO J = 1, N
               DO I = 1, M
                  A( I, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, &
                              ISEED, D, IGRADE, DL, DR, IPVTNG, &
                              IWORK, SPARSE )
                  ENDDO
               ENDDO
         ELSE IF( ISYM == 2 ) THEN
            DO J = 1, N
               DO I = 1, J
                  A( I, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, &
                              ISEED, D, IGRADE, DL, DR, IPVTNG, &
                              IWORK, SPARSE )
                  A( J, I ) = A( I, J )
                  ENDDO
               ENDDO
         END IF
!
      ELSE IF( IPACK == 1 ) THEN
!
         DO J = 1, N
            DO I = 1, J
               A( I, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, &
                           D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
               IF( I /= J ) &
                  A( J, I ) = CZERO
               ENDDO
            ENDDO
!
      ELSE IF( IPACK == 2 ) THEN
!
         DO J = 1, N
            DO I = 1, J
               IF( ISYM == 0 ) THEN
                  A( J, I ) = CONJG( CLATM2( M, N, I, J, KL, KU, &
                              IDIST, ISEED, D, IGRADE, DL, DR, &
                              IPVTNG, IWORK, SPARSE ) )
               ELSE
                  A( J, I ) = CLATM2( M, N, I, J, KL, KU, IDIST, &
                              ISEED, D, IGRADE, DL, DR, IPVTNG, &
                              IWORK, SPARSE )
               END IF
               IF( I /= J ) &
                  A( I, J ) = CZERO
               ENDDO
            ENDDO
!
      ELSE IF( IPACK == 3 ) THEN
!
         ISUB = 0
         JSUB = 1
         DO J = 1, N
            DO I = 1, J
               ISUB = ISUB + 1
               IF( ISUB > LDA ) THEN
                  ISUB = 1
                  JSUB = JSUB + 1
               END IF
               A( ISUB, JSUB ) = CLATM2( M, N, I, J, KL, KU, IDIST, &
                                 ISEED, D, IGRADE, DL, DR, IPVTNG, &
                                 IWORK, SPARSE )
               ENDDO
            ENDDO
!
      ELSE IF( IPACK == 4 ) THEN
!
         IF( ISYM == 0 .OR. ISYM == 2 ) THEN
            DO J = 1, N
               DO I = 1, J
!
!                    Compute K = location of (I,J) entry in packed array
!
                  IF( I == 1 ) THEN
                     K = J
                  ELSE
                     K = N*( N+1 ) / 2 - ( N-I+1 )*( N-I+2 ) / 2 + &
                         J - I + 1
                  END IF
!
!                    Convert K to (ISUB,JSUB) location
!
                  JSUB = ( K-1 ) / LDA + 1
                  ISUB = K - LDA*( JSUB-1 )
!
                  A( ISUB, JSUB ) = CLATM2( M, N, I, J, KL, KU, &
                                    IDIST, ISEED, D, IGRADE, DL, DR, &
                                    IPVTNG, IWORK, SPARSE )
                  IF( ISYM == 0 ) &
                     A( ISUB, JSUB ) = CONJG( A( ISUB, JSUB ) )
                  ENDDO
               ENDDO
         ELSE
            ISUB = 0
            JSUB = 1
            DO J = 1, N
               DO I = J, M
                  ISUB = ISUB + 1
                  IF( ISUB > LDA ) THEN
                     ISUB = 1
                     JSUB = JSUB + 1
                  END IF
                  A( ISUB, JSUB ) = CLATM2( M, N, I, J, KL, KU, &
                                    IDIST, ISEED, D, IGRADE, DL, DR, &
                                    IPVTNG, IWORK, SPARSE )
                  ENDDO
               ENDDO
         END IF
!
      ELSE IF( IPACK == 5 ) THEN
!
         DO J = 1, N
            DO I = J - KUU, J
               IF( I < 1 ) THEN
                  A( J-I+1, I+N ) = CZERO
               ELSE
                  IF( ISYM == 0 ) THEN
                     A( J-I+1, I ) = CONJG( CLATM2( M, N, I, J, KL, &
                                     KU, IDIST, ISEED, D, IGRADE, DL, &
                                     DR, IPVTNG, IWORK, SPARSE ) )
                  ELSE
                     A( J-I+1, I ) = CLATM2( M, N, I, J, KL, KU, &
                                     IDIST, ISEED, D, IGRADE, DL, DR, &
                                     IPVTNG, IWORK, SPARSE )
                  END IF
               END IF
               ENDDO
            ENDDO
!
      ELSE IF( IPACK == 6 ) THEN
!
         DO J = 1, N
            DO I = J - KUU, J
               A( I-J+KUU+1, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, &
                                   ISEED, D, IGRADE, DL, DR, IPVTNG, &
                                   IWORK, SPARSE )
               ENDDO
            ENDDO
!
      ELSE IF( IPACK == 7 ) THEN
!
         IF( ISYM /= 1 ) THEN
            DO J = 1, N
               DO I = J - KUU, J
                  A( I-J+KUU+1, J ) = CLATM2( M, N, I, J, KL, KU, &
                                      IDIST, ISEED, D, IGRADE, DL, &
                                      DR, IPVTNG, IWORK, SPARSE )
                  IF( I < 1 ) &
                     A( J-I+1+KUU, I+N ) = CZERO
                  IF( I >= 1 .AND. I /= J ) THEN
                     IF( ISYM == 0 ) THEN
                        A( J-I+1+KUU, I ) = CONJG( A( I-J+KUU+1, &
                                            J ) )
                     ELSE
                        A( J-I+1+KUU, I ) = A( I-J+KUU+1, J )
                     END IF
                  END IF
                  ENDDO
               ENDDO
         ELSE IF( ISYM == 1 ) THEN
            DO J = 1, N
               DO I = J - KUU, J + KLL
                  A( I-J+KUU+1, J ) = CLATM2( M, N, I, J, KL, KU, &
                                      IDIST, ISEED, D, IGRADE, DL, &
                                      DR, IPVTNG, IWORK, SPARSE )
                  ENDDO
               ENDDO
         END IF
!
      END IF
!
   END IF
!
!     5)      Scaling the norm
!
   IF( IPACK == 0 ) THEN
      ONORM = CLANGE( 'M', M, N, A, LDA, TEMPA )
   ELSE IF( IPACK == 1 ) THEN
      ONORM = CLANSY( 'M', 'U', N, A, LDA, TEMPA )
   ELSE IF( IPACK == 2 ) THEN
      ONORM = CLANSY( 'M', 'L', N, A, LDA, TEMPA )
   ELSE IF( IPACK == 3 ) THEN
      ONORM = CLANSP( 'M', 'U', N, A, TEMPA )
   ELSE IF( IPACK == 4 ) THEN
      ONORM = CLANSP( 'M', 'L', N, A, TEMPA )
   ELSE IF( IPACK == 5 ) THEN
      ONORM = CLANSB( 'M', 'L', N, KLL, A, LDA, TEMPA )
   ELSE IF( IPACK == 6 ) THEN
      ONORM = CLANSB( 'M', 'U', N, KUU, A, LDA, TEMPA )
   ELSE IF( IPACK == 7 ) THEN
      ONORM = CLANGB( 'M', N, KLL, KUU, A, LDA, TEMPA )
   END IF
!
   IF( ANORM >= ZERO ) THEN
!
      IF( ANORM > ZERO .AND. ONORM == ZERO ) THEN
!
!           Desired scaling impossible
!
         INFO = 5
         RETURN
!
      ELSE IF( ( ANORM > ONE .AND. ONORM < ONE ) .OR. &
               ( ANORM < ONE .AND. ONORM > ONE ) ) THEN
!
!           Scale carefully to avoid over / underflow
!
         IF( IPACK <= 2 ) THEN
            DO J = 1, N
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL CSSCAL( M, ONE / ONORM, A( 1, J ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : CSSCAL : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL CSSCAL( M, ANORM, A( 1, J ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : CSSCAL : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               ENDDO
!
         ELSE IF( IPACK == 3 .OR. IPACK == 4 ) THEN
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL CSSCAL( N*( N+1 ) / 2, ONE / ONORM, A, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : CSSCAL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL CSSCAL( N*( N+1 ) / 2, ANORM, A, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : CSSCAL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
         ELSE IF( IPACK >= 5 ) THEN
!
            DO J = 1, N
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL CSSCAL( KLL+KUU+1, ONE / ONORM, A( 1, J ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : CSSCAL : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL CSSCAL( KLL+KUU+1, ANORM, A( 1, J ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : CSSCAL : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               ENDDO
!
         END IF
!
      ELSE
!
!           Scale straightforwardly
!
         IF( IPACK <= 2 ) THEN
            DO J = 1, N
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL CSSCAL( M, ANORM / ONORM, A( 1, J ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : CSSCAL : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               ENDDO
!
         ELSE IF( IPACK == 3 .OR. IPACK == 4 ) THEN
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL CSSCAL( N*( N+1 ) / 2, ANORM / ONORM, A, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : CSSCAL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
         ELSE IF( IPACK >= 5 ) THEN
!
            DO J = 1, N
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL CSSCAL( KLL+KUU+1, ANORM / ONORM, A( 1, J ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : CSSCAL : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               ENDDO
         END IF
!
      END IF
!
   END IF
!
!     End of CLATMR
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        


