!> \brief \b SCHKBD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SCHKBD( NSIZES, MVAL, NVAL, NTYPES, DOTYPE, NRHS,
!                          ISEED, THRESH, A, LDA, BD, BE, S1, S2, X, LDX,
!                          Y, Z, Q, LDQ, PT, LDPT, U, VT, WORK, LWORK,
!                          IWORK, NOUT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDPT, LDQ, LDX, LWORK, NOUT, NRHS,
!      $                   NSIZES, NTYPES
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), IWORK( * ), MVAL( * ), NVAL( * )
!       REAL               A( LDA, * ), BD( * ), BE( * ), PT( LDPT, * ),
!      $                   Q( LDQ, * ), S1( * ), S2( * ), U( LDPT, * ),
!      $                   VT( LDPT, * ), WORK( * ), X( LDX, * ),
!      $                   Y( LDX, * ), Z( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SCHKBD checks the singular value decomposition (SVD) routines.
!>
!> SGEBRD reduces a real general m by n matrix A to upper or lower
!> bidiagonal form B by an orthogonal transformation:  Q' * A * P = B
!> (or A = Q * B * P').  The matrix B is upper bidiagonal if m >= n
!> and lower bidiagonal if m < n.
!>
!> SORGBR generates the orthogonal matrices Q and P' from SGEBRD.
!> Note that Q and P are not necessarily square.
!>
!> SBDSQR computes the singular value decomposition of the bidiagonal
!> matrix B as B = U S V'.  It is called three times to compute
!>    1)  B = U S1 V', where S1 is the diagonal matrix of singular
!>        values and the columns of the matrices U and V are the left
!>        and right singular vectors, respectively, of B.
!>    2)  Same as 1), but the singular values are stored in S2 and the
!>        singular vectors are not computed.
!>    3)  A = (UQ) S (P'V'), the SVD of the original matrix A.
!> In addition, SBDSQR has an option to apply the left orthogonal matrix
!> U to a matrix X, useful in least squares applications.
!>
!> SBDSDC computes the singular value decomposition of the bidiagonal
!> matrix B as B = U S V' using divide-and-conquer. It is called twice
!> to compute
!>    1) B = U S1 V', where S1 is the diagonal matrix of singular
!>        values and the columns of the matrices U and V are the left
!>        and right singular vectors, respectively, of B.
!>    2) Same as 1), but the singular values are stored in S2 and the
!>        singular vectors are not computed.
!>
!>  SBDSVDX computes the singular value decomposition of the bidiagonal
!>  matrix B as B = U S V' using bisection and inverse iteration. It is
!>  called six times to compute
!>     1) B = U S1 V', RANGE='A', where S1 is the diagonal matrix of singular
!>         values and the columns of the matrices U and V are the left
!>         and right singular vectors, respectively, of B.
!>     2) Same as 1), but the singular values are stored in S2 and the
!>         singular vectors are not computed.
!>     3) B = U S1 V', RANGE='I', with where S1 is the diagonal matrix of singular
!>         values and the columns of the matrices U and V are the left
!>         and right singular vectors, respectively, of B
!>     4) Same as 3), but the singular values are stored in S2 and the
!>         singular vectors are not computed.
!>     5) B = U S1 V', RANGE='V', with where S1 is the diagonal matrix of singular
!>         values and the columns of the matrices U and V are the left
!>         and right singular vectors, respectively, of B
!>     6) Same as 5), but the singular values are stored in S2 and the
!>         singular vectors are not computed.
!>
!> For each pair of matrix dimensions (M,N) and each selected matrix
!> type, an M by N matrix A and an M by NRHS matrix X are generated.
!> The problem dimensions are as follows
!>    A:          M x N
!>    Q:          M x min(M,N) (but M x M if NRHS > 0)
!>    P:          min(M,N) x N
!>    B:          min(M,N) x min(M,N)
!>    U, V:       min(M,N) x min(M,N)
!>    S1, S2      diagonal, order min(M,N)
!>    X:          M x NRHS
!>
!> For each generated matrix, 14 tests are performed:
!>
!> Test SGEBRD and SORGBR
!>
!> (1)   | A - Q B PT | / ( |A| max(M,N) ulp ), PT = P'
!>
!> (2)   | I - Q' Q | / ( M ulp )
!>
!> (3)   | I - PT PT' | / ( N ulp )
!>
!> Test SBDSQR on bidiagonal matrix B
!>
!> (4)   | B - U S1 VT | / ( |B| min(M,N) ulp ), VT = V'
!>
!> (5)   | Y - U Z | / ( |Y| max(min(M,N),k) ulp ), where Y = Q' X
!>                                                  and   Z = U' Y.
!> (6)   | I - U' U | / ( min(M,N) ulp )
!>
!> (7)   | I - VT VT' | / ( min(M,N) ulp )
!>
!> (8)   S1 contains min(M,N) nonnegative values in decreasing order.
!>       (Return 0 if true, 1/ULP if false.)
!>
!> (9)   | S1 - S2 | / ( |S1| ulp ), where S2 is computed without
!>                                   computing U and V.
!>
!> (10)  0 if the true singular values of B are within THRESH of
!>       those in S1.  2*THRESH if they are not.  (Tested using
!>       SSVDCH)
!>
!> Test SBDSQR on matrix A
!>
!> (11)  | A - (QU) S (VT PT) | / ( |A| max(M,N) ulp )
!>
!> (12)  | X - (QU) Z | / ( |X| max(M,k) ulp )
!>
!> (13)  | I - (QU)'(QU) | / ( M ulp )
!>
!> (14)  | I - (VT PT) (PT'VT') | / ( N ulp )
!>
!> Test SBDSDC on bidiagonal matrix B
!>
!> (15)  | B - U S1 VT | / ( |B| min(M,N) ulp ), VT = V'
!>
!> (16)  | I - U' U | / ( min(M,N) ulp )
!>
!> (17)  | I - VT VT' | / ( min(M,N) ulp )
!>
!> (18)  S1 contains min(M,N) nonnegative values in decreasing order.
!>       (Return 0 if true, 1/ULP if false.)
!>
!> (19)  | S1 - S2 | / ( |S1| ulp ), where S2 is computed without
!>                                   computing U and V.
!>  Test SBDSVDX on bidiagonal matrix B
!>
!>  (20)  | B - U S1 VT | / ( |B| min(M,N) ulp ), VT = V'
!>
!>  (21)  | I - U' U | / ( min(M,N) ulp )
!>
!>  (22)  | I - VT VT' | / ( min(M,N) ulp )
!>
!>  (23)  S1 contains min(M,N) nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!>  (24)  | S1 - S2 | / ( |S1| ulp ), where S2 is computed without
!>                                    computing U and V.
!>
!>  (25)  | S1 - U' B VT' | / ( |S| n ulp )    SBDSVDX('V', 'I')
!>
!>  (26)  | I - U' U | / ( min(M,N) ulp )
!>
!>  (27)  | I - VT VT' | / ( min(M,N) ulp )
!>
!>  (28)  S1 contains min(M,N) nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!>  (29)  | S1 - S2 | / ( |S1| ulp ), where S2 is computed without
!>                                    computing U and V.
!>
!>  (30)  | S1 - U' B VT' | / ( |S1| n ulp )   SBDSVDX('V', 'V')
!>
!>  (31)  | I - U' U | / ( min(M,N) ulp )
!>
!>  (32)  | I - VT VT' | / ( min(M,N) ulp )
!>
!>  (33)  S1 contains min(M,N) nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!>  (34)  | S1 - S2 | / ( |S1| ulp ), where S2 is computed without
!>                                    computing U and V.
!>
!> The possible matrix types are
!>
!> (1)  The zero matrix.
!> (2)  The identity matrix.
!>
!> (3)  A diagonal matrix with evenly spaced entries
!>      1, ..., ULP  and random signs.
!>      (ULP = (first number larger than 1) - 1 )
!> (4)  A diagonal matrix with geometrically spaced entries
!>      1, ..., ULP  and random signs.
!> (5)  A diagonal matrix with "clustered" entries 1, ULP, ..., ULP
!>      and random signs.
!>
!> (6)  Same as (3), but multiplied by SQRT( overflow threshold )
!> (7)  Same as (3), but multiplied by SQRT( underflow threshold )
!>
!> (8)  A matrix of the form  U D V, where U and V are orthogonal and
!>      D has evenly spaced entries 1, ..., ULP with random signs
!>      on the diagonal.
!>
!> (9)  A matrix of the form  U D V, where U and V are orthogonal and
!>      D has geometrically spaced entries 1, ..., ULP with random
!>      signs on the diagonal.
!>
!> (10) A matrix of the form  U D V, where U and V are orthogonal and
!>      D has "clustered" entries 1, ULP,..., ULP with random
!>      signs on the diagonal.
!>
!> (11) Same as (8), but multiplied by SQRT( overflow threshold )
!> (12) Same as (8), but multiplied by SQRT( underflow threshold )
!>
!> (13) Rectangular matrix with random entries chosen from (-1,1).
!> (14) Same as (13), but multiplied by SQRT( overflow threshold )
!> (15) Same as (13), but multiplied by SQRT( underflow threshold )
!>
!> Special case:
!> (16) A bidiagonal matrix with random entries chosen from a
!>      logarithmic distribution on [ulp^2,ulp^(-2)]  (I.e., each
!>      entry is  e^x, where x is chosen uniformly on
!>      [ 2 log(ulp), -2 log(ulp) ] .)  For *this* type:
!>      (a) SGEBRD is not called to reduce it to bidiagonal form.
!>      (b) the bidiagonal is  min(M,N) x min(M,N); if M<N, the
!>          matrix will be lower bidiagonal, otherwise upper.
!>      (c) only tests 5--8 and 14 are performed.
!>
!> A subset of the full set of matrix types may be selected through
!> the logical array DOTYPE.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NSIZES
!> \verbatim
!>          NSIZES is INTEGER
!>          The number of values of M and N contained in the vectors
!>          MVAL and NVAL.  The matrix sizes are used in pairs (M,N).
!> \endverbatim
!>
!> \param[in] MVAL
!> \verbatim
!>          MVAL is INTEGER array, dimension (NM)
!>          The values of the matrix row dimension M.
!> \endverbatim
!>
!> \param[in] NVAL
!> \verbatim
!>          NVAL is INTEGER array, dimension (NM)
!>          The values of the matrix column dimension N.
!> \endverbatim
!>
!> \param[in] NTYPES
!> \verbatim
!>          NTYPES is INTEGER
!>          The number of elements in DOTYPE.   If it is zero, SCHKBD
!>          does nothing.  It must be at least zero.  If it is MAXTYP+1
!>          and NSIZES is 1, then an additional type, MAXTYP+1 is
!>          defined, which is to use whatever matrices are in A and B.
!>          This is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
!>          DOTYPE(MAXTYP+1) is .TRUE. .
!> \endverbatim
!>
!> \param[in] DOTYPE
!> \verbatim
!>          DOTYPE is LOGICAL array, dimension (NTYPES)
!>          If DOTYPE(j) is .TRUE., then for each size (m,n), a matrix
!>          of type j will be generated.  If NTYPES is smaller than the
!>          maximum number of types defined (PARAMETER MAXTYP), then
!>          types NTYPES+1 through MAXTYP will not be generated.  If
!>          NTYPES is larger than MAXTYP, DOTYPE(MAXTYP+1) through
!>          DOTYPE(NTYPES) will be ignored.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns in the "right-hand side" matrices X, Y,
!>          and Z, used in testing SBDSQR.  If NRHS = 0, then the
!>          operations on the right-hand side will not be tested.
!>          NRHS must be at least 0.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry ISEED specifies the seed of the random number
!>          generator. The array elements should be between 0 and 4095;
!>          if not they will be reduced mod 4096.  Also, ISEED(4) must
!>          be odd.  The values of ISEED are changed on exit, and can be
!>          used in the next call to SCHKBD to continue the same random
!>          number sequence.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is REAL
!>          The threshold value for the test ratios.  A result is
!>          included in the output file if RESULT >= THRESH.  To have
!>          every test ratio printed, use THRESH = 0.  Note that the
!>          expected value of the test ratios is O(1), so THRESH should
!>          be a reasonably small multiple of 1, e.g., 10 or 100.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,NMAX)
!>          where NMAX is the maximum value of N in NVAL.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,MMAX),
!>          where MMAX is the maximum value of M in MVAL.
!> \endverbatim
!>
!> \param[out] BD
!> \verbatim
!>          BD is REAL array, dimension
!>                      (max(min(MVAL(j),NVAL(j))))
!> \endverbatim
!>
!> \param[out] BE
!> \verbatim
!>          BE is REAL array, dimension
!>                      (max(min(MVAL(j),NVAL(j))))
!> \endverbatim
!>
!> \param[out] S1
!> \verbatim
!>          S1 is REAL array, dimension
!>                      (max(min(MVAL(j),NVAL(j))))
!> \endverbatim
!>
!> \param[out] S2
!> \verbatim
!>          S2 is REAL array, dimension
!>                      (max(min(MVAL(j),NVAL(j))))
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is REAL array, dimension (LDX,NRHS)
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the arrays X, Y, and Z.
!>          LDX >= max(1,MMAX)
!> \endverbatim
!>
!> \param[out] Y
!> \verbatim
!>          Y is REAL array, dimension (LDX,NRHS)
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is REAL array, dimension (LDX,NRHS)
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDQ,MMAX)
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= max(1,MMAX).
!> \endverbatim
!>
!> \param[out] PT
!> \verbatim
!>          PT is REAL array, dimension (LDPT,NMAX)
!> \endverbatim
!>
!> \param[in] LDPT
!> \verbatim
!>          LDPT is INTEGER
!>          The leading dimension of the arrays PT, U, and V.
!>          LDPT >= max(1, max(min(MVAL(j),NVAL(j)))).
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is REAL array, dimension
!>                      (LDPT,max(min(MVAL(j),NVAL(j))))
!> \endverbatim
!>
!> \param[out] VT
!> \verbatim
!>          VT is REAL array, dimension
!>                      (LDPT,max(min(MVAL(j),NVAL(j))))
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The number of entries in WORK.  This must be at least
!>          3(M+N) and  M(M + max(M,N,k) + 1) + N*min(M,N)  for all
!>          pairs  (M,N)=(MM(j),NN(j))
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension at least 8*min(M,N)
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The FORTRAN unit number for printing out error messages
!>          (e.g., if a routine returns IINFO not equal to 0.)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          If 0, then everything ran OK.
!>           -1: NSIZES < 0
!>           -2: Some MM(j) < 0
!>           -3: Some NN(j) < 0
!>           -4: NTYPES < 0
!>           -6: NRHS  < 0
!>           -8: THRESH < 0
!>          -11: LDA < 1 or LDA < MMAX, where MMAX is max( MM(j) ).
!>          -17: LDB < 1 or LDB < MMAX.
!>          -21: LDQ < 1 or LDQ < MMAX.
!>          -23: LDPT< 1 or LDPT< MNMAX.
!>          -27: LWORK too small.
!>          If  SLATMR, SLATMS, SGEBRD, SORGBR, or SBDSQR,
!>              returns an error code, the
!>              absolute value of it is returned.
!>
!>-----------------------------------------------------------------------
!>
!>     Some Local Variables and Parameters:
!>     ---- ----- --------- --- ----------
!>
!>     0.0E+0, 1.0E+0       Real 0 and 1.
!>     MAXTYP          The number of types defined.
!>     NTEST           The number of tests performed, or which can
!>                     be performed so far, for the current matrix.
!>     MMAX            Largest value in NN.
!>     NMAX            Largest value in NN.
!>     MNMIN           min(MM(j), NN(j)) (the dimension of the bidiagonal
!>                     matrix.)
!>     MNMAX           The maximum value of MNMIN for j=1,...,NSIZES.
!>     NFAIL           The number of tests which have exceeded THRESH
!>     COND, IMODE     Values to be passed to the matrix generators.
!>     ANORM           Norm of A; passed to matrix generators.
!>
!>     OVFL, UNFL      Overflow and underflow thresholds.
!>     RTOVFL, RTUNFL  Square roots of the previous 2 values.
!>     ULP, ULPINV     Finest relative precision and its inverse.
!>
!>             The following four arrays decode JTYPE:
!>     KTYPE(j)        The general type (1-10) for type "j".
!>     KMODE(j)        The MODE value to be passed to the matrix
!>                     generator for type "j".
!>     KMAGN(j)        The order of magnitude ( O(1),
!>                     O(overflow^(1/2) ), O(underflow^(1/2) )
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
!> \ingroup single_eig
!
!  =====================================================================
   SUBROUTINE SCHKBD( NSIZES, MVAL, NVAL, NTYPES, DOTYPE, NRHS, &
                      ISEED, THRESH, A, LDA, BD, BE, S1, S2, X, LDX, &
                      Y, Z, Q, LDQ, PT, LDPT, U, VT, WORK, LWORK, &
                      IWORK, NOUT, INFO )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, LDA, LDPT, LDQ, LDX, LWORK, NOUT, NRHS, &
                      NSIZES, NTYPES
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            ISEED( 4 ), IWORK( * ), MVAL( * ), NVAL( * )
   REAL               A( LDA, * ), BD( * ), BE( * ), PT( LDPT, * ), &
                      Q( LDQ, * ), S1( * ), S2( * ), U( LDPT, * ), &
                      VT( LDPT, * ), WORK( * ), X( LDX, * ), &
                      Y( LDX, * ), Z( LDX, * )
!     ..
!
! ======================================================================
!
!     .. Parameters ..
   INTEGER            MAXTYP
   PARAMETER          ( MAXTYP = 16 )
!     ..
!     .. Local Scalars ..
   LOGICAL            BADMM, BADNN, BIDIAG
   CHARACTER          UPLO
   CHARACTER*3        PATH
   INTEGER            I, IINFO, IL, IMODE, ITEMP, ITYPE, IU, IWBD, &
                      IWBE, IWBS, IWBZ, IWWORK, J, JCOL, JSIZE, &
                      JTYPE, LOG2UI, M, MINWRK, MMAX, MNMAX, MNMIN, &
                      MNMIN2, MQ, MTYPES, N, NFAIL, NMAX, &
                      NS1, NS2, NTEST
   REAL               ABSTOL, AMNINV, ANORM, COND, OVFL, RTOVFL, &
                      RTUNFL, TEMP1, TEMP2, ULP, ULPINV, UNFL, &
                      VL, VU
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   INTEGER            IDUM( 1 ), IOLDSD( 4 ), ISEED2( 4 ), &
                      KMAGN( MAXTYP ), KMODE( MAXTYP ), &
                      KTYPE( MAXTYP )
   REAL               DUM( 1 ), DUMMA( 1 ), RESULT( 40 )
!     ..
!     .. External Functions ..
   REAL               SLAMCH, SLARND, SSXT1
   EXTERNAL           SLAMCH, SLARND, SSXT1
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALASUM, SBDSDC, SBDSQR, SBDSVDX, SBDT01, &
                      SBDT02, SBDT03, SBDT04, SCOPY, SGEBRD, &
                      SGEMM, SLACPY, SLAHD2, SLASET, SLATMR, &
                      SLATMS, SORGBR, SORT01, XERBLA
!     ..
!     .. Scalars in Common ..
   LOGICAL            LERR, OK
   CHARACTER*32       SRNAMT
   INTEGER            INFOT, NUNIT
!     ..
!     .. Common blocks ..
   COMMON             / INFOC / INFOT, NUNIT, OK, LERR
   COMMON             / SRNAMC / SRNAMT
!     ..
!     .. Data statements ..
   DATA            KTYPE / 1, 2, 5*4, 5*6, 3*9, 10 /
   DATA            KMAGN / 2*1, 3*1, 2, 3, 3*1, 2, 3, 1, 2, 3, 0 /
   DATA            KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 0 /
!     ..
!     .. Executable Statements ..
!
!     Check for errors
!
   INFO = 0
!
   BADMM = ANY(MVAL(1:NSIZES) < 0)
   BADNN = ANY(NVAL(1:NSIZES) < 0)
   MMAX = MAX(1, MAXVAL(MVAL(1:NSIZES)))
   NMAX = MAX(1, MAXVAL(NVAL(1:NSIZES)))
   MNMAX = MAX(1, MINVAL(MVAL(1:NSIZES)), MINVAL(NVAL(1:NSIZES)))
   MINWRK = 1
   DO J = 1, NSIZES
      MINWRK = MAX( MINWRK, 3*( MVAL( J )+NVAL( J ) ), &
               MVAL( J )*( MVAL( J )+MAX( MVAL( J ), NVAL( J ), &
               NRHS )+1 )+NVAL( J )*MIN( NVAL( J ), MVAL( J ) ) )
   ENDDO
!
!     Check for errors
!
   IF( NSIZES < 0 ) THEN
      INFO = -1
   ELSE IF( BADMM ) THEN
      INFO = -2
   ELSE IF( BADNN ) THEN
      INFO = -3
   ELSE IF( NTYPES < 0 ) THEN
      INFO = -4
   ELSE IF( NRHS < 0 ) THEN
      INFO = -6
   ELSE IF( LDA < MMAX ) THEN
      INFO = -11
   ELSE IF( LDX < MMAX ) THEN
      INFO = -17
   ELSE IF( LDQ < MMAX ) THEN
      INFO = -21
   ELSE IF( LDPT < MNMAX ) THEN
      INFO = -23
   ELSE IF( MINWRK > LWORK ) THEN
      INFO = -27
   END IF
!
   IF( INFO /= 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'SCHKBD', -INFO )
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
!     Initialize constants
!
   PATH( 1: 1 ) = 'Single precision'
   PATH( 2: 3 ) = 'BD'
   NFAIL = 0
   NTEST = 0
   UNFL = SLAMCH( 'Safe minimum' )
   OVFL = SLAMCH( 'Overflow' )
   ULP = SLAMCH( 'Precision' )
   ULPINV = 1.0E+0 / ULP
   LOG2UI = INT( LOG( ULPINV ) / LOG( 2.0E+0 ) )
   RTUNFL = SQRT( UNFL )
   RTOVFL = SQRT( OVFL )
   INFOT = 0
   ABSTOL = 2*UNFL
!
!     Loop over sizes, types
!
   DO JSIZE = 1, NSIZES
      M = MVAL( JSIZE )
      N = NVAL( JSIZE )
      MNMIN = MIN( M, N )
      AMNINV = 1.0E+0 / MAX( M, N, 1 )
!
      IF( NSIZES /= 1 ) THEN
         MTYPES = MIN( MAXTYP, NTYPES )
      ELSE
         MTYPES = MIN( MAXTYP+1, NTYPES )
      END IF
!
      DO JTYPE = 1, MTYPES
         IF(DOTYPE( JTYPE ) ) THEN
!
         IOLDSD(1:4) = ISEED(1:4)
!
         RESULT(1:34) = -1.0E0
!
         UPLO = ' '
!
!           Compute "A"
!
!           Control parameters:
!
!           KMAGN  KMODE        KTYPE
!       =1  O(1)   clustered 1  zero
!       =2  large  clustered 2  identity
!       =3  small  exponential  (none)
!       =4         arithmetic   diagonal, (w/ eigenvalues)
!       =5         random       symmetric, w/ eigenvalues
!       =6                      nonsymmetric, w/ singular values
!       =7                      random diagonal
!       =8                      random symmetric
!       =9                      random nonsymmetric
!       =10                     random bidiagonal (log. distrib.)
!
         IF( MTYPES > MAXTYP ) GO TO 100
!
         ITYPE = KTYPE( JTYPE )
         IMODE = KMODE( JTYPE )
!
!           Compute norm
!
         SELECT CASE (KMAGN(JTYPE))
          CASE (1)
           ANORM = 1.0E0
          CASE (2)
           ANORM = ( RTOVFL*ULP )*AMNINV
          CASE (3)
           ANORM = RTUNFL*MAX( M, N )*ULPINV
         END SELECT
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SLASET( 'Full', LDA, N, 0.0E+0, 0.0E+0, A, LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IINFO = 0
         COND = ULPINV
!
         BIDIAG = .FALSE.
         IF( ITYPE == 1 ) THEN
!
!              Zero matrix
!
            IINFO = 0
!
         ELSE IF( ITYPE == 2 ) THEN
!
!              Identity
!
            DO JCOL = 1, MNMIN
               A( JCOL, JCOL ) = ANORM
            ENDDO
!
         ELSE IF( ITYPE == 4 ) THEN
!
!              Diagonal Matrix, [Eigen]values Specified
!
            CALL SLATMS( MNMIN, MNMIN, 'S', ISEED, 'N', WORK, IMODE, &
                         COND, ANORM, 0, 0, 'N', A, LDA, &
                         WORK( MNMIN+1 ), IINFO )
!
         ELSE IF( ITYPE == 5 ) THEN
!
!              Symmetric, eigenvalues specified
!
            CALL SLATMS( MNMIN, MNMIN, 'S', ISEED, 'S', WORK, IMODE, &
                         COND, ANORM, M, N, 'N', A, LDA, &
                         WORK( MNMIN+1 ), IINFO )
!
         ELSE IF( ITYPE == 6 ) THEN
!
!              Nonsymmetric, singular values specified
!
            CALL SLATMS( M, N, 'S', ISEED, 'N', WORK, IMODE, COND, &
                         ANORM, M, N, 'N', A, LDA, WORK( MNMIN+1 ), &
                         IINFO )
!
         ELSE IF( ITYPE == 7 ) THEN
!
!              Diagonal, random entries
!
            CALL SLATMR( MNMIN, MNMIN, 'S', ISEED, 'N', WORK, 6, 1.0E+0, &
                         1.0E+0, 'T', 'N', WORK( MNMIN+1 ), 1, 1.0E+0, &
                         WORK( 2*MNMIN+1 ), 1, 1.0E+0, 'N', IWORK, 0, 0, &
                         0.0E+0, ANORM, 'NO', A, LDA, IWORK, IINFO )
!
         ELSE IF( ITYPE == 8 ) THEN
!
!              Symmetric, random entries
!
            CALL SLATMR( MNMIN, MNMIN, 'S', ISEED, 'S', WORK, 6, 1.0E+0, &
                         1.0E+0, 'T', 'N', WORK( MNMIN+1 ), 1, 1.0E+0, &
                         WORK( M+MNMIN+1 ), 1, 1.0E+0, 'N', IWORK, M, N, &
                         0.0E+0, ANORM, 'NO', A, LDA, IWORK, IINFO )
!
         ELSE IF( ITYPE == 9 ) THEN
!
!              Nonsymmetric, random entries
!
            CALL SLATMR( M, N, 'S', ISEED, 'N', WORK, 6, 1.0E+0, 1.0E+0, &
                         'T', 'N', WORK( MNMIN+1 ), 1, 1.0E+0, &
                         WORK( M+MNMIN+1 ), 1, 1.0E+0, 'N', IWORK, M, N, &
                         0.0E+0, ANORM, 'NO', A, LDA, IWORK, IINFO )
!
         ELSE IF( ITYPE == 10 ) THEN
!
!              Bidiagonal, random entries
!
            TEMP1 = -2.0E+0*LOG( ULP )
            DO J = 1, MNMIN
               BD( J ) = EXP( TEMP1*SLARND( 2, ISEED ) )
               IF( J < MNMIN ) BE( J ) = EXP( TEMP1*SLARND( 2, ISEED ) )
            ENDDO
!
            IINFO = 0
            BIDIAG = .TRUE.
            IF( M >= N ) THEN
               UPLO = 'U'
            ELSE
               UPLO = 'L'
            END IF
         ELSE
            IINFO = 1
         END IF
!
         IF( IINFO == 0 ) THEN
!
!              Generate Right-Hand Side
!
            IF( BIDIAG ) THEN
               CALL SLATMR( MNMIN, NRHS, 'S', ISEED, 'N', WORK, 6, &
                            1.0E+0, 1.0E+0, 'T', 'N', WORK( MNMIN+1 ), 1, &
                            1.0E+0, WORK( 2*MNMIN+1 ), 1, 1.0E+0, 'N', &
                            IWORK, MNMIN, NRHS, 0.0E+0, 1.0E+0, 'NO', Y, &
                            LDX, IWORK, IINFO )
            ELSE
               CALL SLATMR( M, NRHS, 'S', ISEED, 'N', WORK, 6, 1.0E+0, &
                            1.0E+0, 'T', 'N', WORK( M+1 ), 1, 1.0E+0, &
                            WORK( 2*M+1 ), 1, 1.0E+0, 'N', IWORK, M, &
                            NRHS, 0.0E+0, 1.0E+0, 'NO', X, LDX, IWORK, &
                            IINFO )
            END IF
         END IF
!
!           Error Exit
!
         IF( IINFO /= 0 ) THEN
            WRITE( NOUT, FMT = 9998 )'Generator', IINFO, M, N, &
               JTYPE, IOLDSD
            INFO = ABS( IINFO )
            RETURN
         END IF
!
  100       CONTINUE
!
!           Call SGEBRD and SORGBR to compute B, Q, and P, do tests.
!
         IF( .NOT.BIDIAG ) THEN
!
!              Compute transformations to reduce A to bidiagonal form:
!              B := Q' * A * P.
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SLACPY( ' ', M, N, A, LDA, Q, LDQ )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGEBRD( M, N, Q, LDQ, BD, BE, WORK, WORK( MNMIN+1 ), &
                         WORK( 2*MNMIN+1 ), LWORK-2*MNMIN, IINFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGEBRD : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
!              Check error code from SGEBRD.
!
            IF( IINFO /= 0 ) THEN
               WRITE( NOUT, FMT = 9998 )'SGEBRD', IINFO, M, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SLACPY( ' ', M, N, Q, LDQ, PT, LDPT )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            IF( M >= N ) THEN
               UPLO = 'U'
            ELSE
               UPLO = 'L'
            END IF
!
!              Generate Q
!
            MQ = M
            IF( NRHS <= 0 ) MQ = MNMIN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SORGBR( 'Q', M, MQ, N, Q, LDQ, WORK, &
                         WORK( 2*MNMIN+1 ), LWORK-2*MNMIN, IINFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SORGBR : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
!              Check error code from SORGBR.
!
            IF( IINFO /= 0 ) THEN
               WRITE( NOUT, FMT = 9998 )'SORGBR(Q)', IINFO, M, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
!
!              Generate P'
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SORGBR( 'P', MNMIN, N, M, PT, LDPT, WORK( MNMIN+1 ), &
                         WORK( 2*MNMIN+1 ), LWORK-2*MNMIN, IINFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SORGBR : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
!              Check error code from SORGBR.
!
            IF( IINFO /= 0 ) THEN
               WRITE( NOUT, FMT = 9998 )'SORGBR(P)', IINFO, M, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
!
!              Apply Q' to an M by NRHS matrix X:  Y := Q' * X.
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGEMM( 'Transpose', 'No transpose', M, NRHS, M, 1.0E+0, &
                        Q, LDQ, X, LDX, 0.0E+0, Y, LDX )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
!              Test 1:  Check the decomposition A := Q * B * PT
!                   2:  Check the orthogonality of Q
!                   3:  Check the orthogonality of PT
!
            CALL SBDT01( M, N, 1, A, LDA, Q, LDQ, BD, BE, PT, LDPT, WORK, RESULT( 1 ) )
            CALL SORT01( 'Columns', M, MQ, Q, LDQ, WORK, LWORK, RESULT( 2 ) )
            CALL SORT01( 'Rows', MNMIN, N, PT, LDPT, WORK, LWORK, RESULT( 3 ) )
         END IF
!
!           Use SBDSQR to form the SVD of the bidiagonal matrix B:
!           B := U * S1 * VT, and compute Z = U' * Y.
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SCOPY( MNMIN, BD, 1, S1, 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( MNMIN > 0 ) CALL SCOPY( MNMIN-1, BE, 1, WORK, 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SLACPY( ' ', M, NRHS, Y, LDX, Z, LDX )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SLASET( 'Full', MNMIN, MNMIN, 0.0E+0, 1.0E+0, U, LDPT )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SLASET( 'Full', MNMIN, MNMIN, 0.0E+0, 1.0E+0, VT, LDPT )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SBDSQR( UPLO, MNMIN, MNMIN, MNMIN, NRHS, S1, WORK, VT, &
                      LDPT, U, LDPT, Z, LDX, WORK( MNMIN+1 ), IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SBDSQR : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
!           Check error code from SBDSQR.
!
         IF( IINFO /= 0 ) THEN
            WRITE( NOUT, FMT = 9998 )'SBDSQR(vects)', IINFO, M, N, JTYPE, IOLDSD
            INFO = ABS( IINFO )
            IF( IINFO < 0 ) THEN
               RETURN
            ELSE
               RESULT( 4 ) = ULPINV
               GO TO 270
            END IF
         END IF
!
!           Use SBDSQR to compute only the singular values of the
!           bidiagonal matrix B;  U, VT, and Z should not be modified.
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SCOPY( MNMIN, BD, 1, S2, 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( MNMIN > 0 ) CALL SCOPY( MNMIN-1, BE, 1, WORK, 1 )
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SBDSQR( UPLO, MNMIN, 0, 0, 0, S2, WORK, VT, LDPT, U, &
                      LDPT, Z, LDX, WORK( MNMIN+1 ), IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SBDSQR : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
!           Check error code from SBDSQR.
!
         IF( IINFO /= 0 ) THEN
            WRITE( NOUT, FMT = 9998 )'SBDSQR(values)', IINFO, M, N, JTYPE, IOLDSD
            INFO = ABS( IINFO )
            IF( IINFO < 0 ) THEN
               RETURN
            ELSE
               RESULT( 9 ) = ULPINV
               GO TO 270
            END IF
         END IF
!
!           Test 4:  Check the decomposition B := U * S1 * VT
!                5:  Check the computation Z := U' * Y
!                6:  Check the orthogonality of U
!                7:  Check the orthogonality of VT
!
         CALL SBDT03( UPLO, MNMIN, 1, BD, BE, U, LDPT, S1, VT, LDPT, WORK, RESULT( 4 ) )
         CALL SBDT02( MNMIN, NRHS, Y, LDX, Z, LDX, U, LDPT, WORK, RESULT( 5 ) )
         CALL SORT01( 'Columns', MNMIN, MNMIN, U, LDPT, WORK, LWORK, RESULT( 6 ) )
         CALL SORT01( 'Rows', MNMIN, MNMIN, VT, LDPT, WORK, LWORK, RESULT( 7 ) )
!
!           Test 8:  Check that the singular values are sorted in
!                    non-increasing order and are non-negative
!
         RESULT( 8 ) = 0.0E+0
         DO I = 1, MNMIN - 1
            IF( S1( I ) < S1( I+1 ) ) &
               RESULT( 8 ) = ULPINV
            IF( S1( I ) < 0.0E+0 ) &
               RESULT( 8 ) = ULPINV
            ENDDO
         IF( MNMIN >= 1 ) THEN
            IF(S1( MNMIN ) < 0.0E+0 ) RESULT( 8 ) = ULPINV
         END IF
!
!           Test 9:  Compare SBDSQR with and without singular vectors
!
         TEMP2 = 0.0E+0
!
         DO J = 1, MNMIN
            TEMP1 = ABS( S1( J )-S2( J ) ) / MAX( SQRT( UNFL )*MAX( S1( 1 ), 1.0E+0 ), &
                    ULP*MAX( ABS( S1( J ) ), ABS( S2( J ) ) ) )
            TEMP2 = MAX( TEMP1, TEMP2 )
            ENDDO
!
         RESULT( 9 ) = TEMP2
!
!           Test 10:  Sturm sequence test of singular values
!                     Go up by factors of two until it succeeds
!
         TEMP1 = THRESH*( 0.5E+0-ULP )
!
         DO J = 0, LOG2UI
!               CALL SSVDCH( MNMIN, BD, BE, S1, TEMP1, IINFO )
            IF( IINFO == 0 ) &
               GO TO 140
            TEMP1 = TEMP1*2.0E+0
            ENDDO
!
  140       CONTINUE
         RESULT( 10 ) = TEMP1
!
!           Use SBDSQR to form the decomposition A := (QU) S (VT PT)
!           from the bidiagonal form A := Q B PT.
!
         IF( .NOT.BIDIAG ) THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SCOPY( MNMIN, BD, 1, S2, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            IF( MNMIN > 0 )  THEN
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SCOPY( MNMIN-1, BE, 1, WORK, 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
            ENDIF
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SBDSQR( UPLO, MNMIN, N, M, NRHS, S2, WORK, PT, LDPT, &
                         Q, LDQ, Y, LDX, WORK( MNMIN+1 ), IINFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SBDSQR : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
!              Test 11:  Check the decomposition A := Q*U * S2 * VT*PT
!                   12:  Check the computation Z := U' * Q' * X
!                   13:  Check the orthogonality of Q*U
!                   14:  Check the orthogonality of VT*PT
!
            CALL SBDT01( M, N, 0, A, LDA, Q, LDQ, S2, DUMMA, PT, &
                         LDPT, WORK, RESULT( 11 ) )
            CALL SBDT02( M, NRHS, X, LDX, Y, LDX, Q, LDQ, WORK, &
                         RESULT( 12 ) )
            CALL SORT01( 'Columns', M, MQ, Q, LDQ, WORK, LWORK, &
                         RESULT( 13 ) )
            CALL SORT01( 'Rows', MNMIN, N, PT, LDPT, WORK, LWORK, &
                         RESULT( 14 ) )
         END IF
!
!           Use SBDSDC to form the SVD of the bidiagonal matrix B:
!           B := U * S1 * VT
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SCOPY( MNMIN, BD, 1, S1, 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( MNMIN > 0 )  THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SCOPY( MNMIN-1, BE, 1, WORK, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ENDIF
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SLASET( 'Full', MNMIN, MNMIN, 0.0E+0, 1.0E+0, U, LDPT )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SLASET( 'Full', MNMIN, MNMIN, 0.0E+0, 1.0E+0, VT, LDPT )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SBDSDC( UPLO, 'I', MNMIN, S1, WORK, U, LDPT, VT, LDPT, &
                      DUM, IDUM, WORK( MNMIN+1 ), IWORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SBDSDC : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
!           Check error code from SBDSDC.
!
         IF( IINFO /= 0 ) THEN
            WRITE( NOUT, FMT = 9998 )'SBDSDC(vects)', IINFO, M, N, &
               JTYPE, IOLDSD
            INFO = ABS( IINFO )
            IF( IINFO < 0 ) THEN
               RETURN
            ELSE
               RESULT( 15 ) = ULPINV
               GO TO 270
            END IF
         END IF
!
!           Use SBDSDC to compute only the singular values of the
!           bidiagonal matrix B;  U and VT should not be modified.
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SCOPY( MNMIN, BD, 1, S2, 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( MNMIN > 0 )  THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SCOPY( MNMIN-1, BE, 1, WORK, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ENDIF
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SBDSDC( UPLO, 'N', MNMIN, S2, WORK, DUM, 1, DUM, 1, &
                      DUM, IDUM, WORK( MNMIN+1 ), IWORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SBDSDC : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
!           Check error code from SBDSDC.
!
         IF( IINFO /= 0 ) THEN
            WRITE( NOUT, FMT = 9998 )'SBDSDC(values)', IINFO, M, N, &
               JTYPE, IOLDSD
            INFO = ABS( IINFO )
            IF( IINFO < 0 ) THEN
               RETURN
            ELSE
               RESULT( 18 ) = ULPINV
               GO TO 270
            END IF
         END IF
!
!           Test 15:  Check the decomposition B := U * S1 * VT
!                16:  Check the orthogonality of U
!                17:  Check the orthogonality of VT
!
         CALL SBDT03( UPLO, MNMIN, 1, BD, BE, U, LDPT, S1, VT, LDPT, &
                      WORK, RESULT( 15 ) )
         CALL SORT01( 'Columns', MNMIN, MNMIN, U, LDPT, WORK, LWORK, &
                      RESULT( 16 ) )
         CALL SORT01( 'Rows', MNMIN, MNMIN, VT, LDPT, WORK, LWORK, &
                      RESULT( 17 ) )
!
!           Test 18:  Check that the singular values are sorted in
!                     non-increasing order and are non-negative
!
         RESULT( 18 ) = 0.0E+0
         DO I = 1, MNMIN - 1
            IF( S1( I ) < S1( I+1 ) ) &
               RESULT( 18 ) = ULPINV
            IF( S1( I ) < 0.0E+0 ) &
               RESULT( 18 ) = ULPINV
            ENDDO
         IF( MNMIN >= 1 ) THEN
            IF( S1( MNMIN ) < 0.0E+0 ) RESULT( 18 ) = ULPINV
         END IF
!
!           Test 19:  Compare SBDSQR with and without singular vectors
!
         TEMP2 = 0.0E+0
!
         DO J = 1, MNMIN
            TEMP1 = ABS( S1( J )-S2( J ) ) / MAX( SQRT( UNFL )*MAX( S1( 1 ), 1.0E+0 ), &
                    ULP*MAX( ABS( S1( 1 ) ), ABS( S2( 1 ) ) ) )
            TEMP2 = MAX( TEMP1, TEMP2 )
         ENDDO
!
         RESULT( 19 ) = TEMP2
!
!
!           Use SBDSVDX to compute the SVD of the bidiagonal matrix B:
!           B := U * S1 * VT
!
         IF( JTYPE == 10 .OR. JTYPE == 16 ) THEN
!              =================================
!              Matrix types temporarily disabled
!              =================================
            RESULT( 20:34 ) = 0.0E+0
            GO TO 270
         END IF
!
         IWBS = 1
         IWBD = IWBS + MNMIN
         IWBE = IWBD + MNMIN
         IWBZ = IWBE + MNMIN
         IWWORK = IWBZ + 2*MNMIN*(MNMIN+1)
         MNMIN2 = MAX( 1,MNMIN*2 )
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SCOPY( MNMIN, BD, 1, WORK( IWBD ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( MNMIN > 0 ) CALL SCOPY( MNMIN-1, BE, 1, WORK( IWBE ), 1 )
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SBDSVDX( UPLO, 'V', 'A', MNMIN, WORK( IWBD ), &
                       WORK( IWBE ), 0.0E+0, 0.0E+0, 0, 0, NS1, S1, &
                       WORK( IWBZ ), MNMIN2, WORK( IWWORK ), &
                       IWORK, IINFO)
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SBDSVDX : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
!           Check error code from SBDSVDX.
!
         IF( IINFO /= 0 ) THEN
            WRITE( NOUT, FMT = 9998 )'SBDSVDX(vects,A)', IINFO, M, N, &
               JTYPE, IOLDSD
            INFO = ABS( IINFO )
            IF( IINFO < 0 ) THEN
               RETURN
            ELSE
               RESULT( 20 ) = ULPINV
               GO TO 270
            END IF
         END IF
!
         J = IWBZ
         DO I = 1, NS1
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SCOPY( MNMIN, WORK( J ), 1, U( 1,I ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            J = J + MNMIN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SCOPY( MNMIN, WORK( J ), 1, VT( I,1 ), LDPT )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            J = J + MNMIN
         ENDDO
!
!           Use SBDSVDX to compute only the singular values of the
!           bidiagonal matrix B;  U and VT should not be modified.
!
         IF( JTYPE == 9 ) THEN
!              =================================
!              Matrix types temporarily disabled
!              =================================
            RESULT( 24 ) = 0.0E+0
            GO TO 270
         END IF
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SCOPY( MNMIN, BD, 1, WORK( IWBD ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( MNMIN > 0 ) CALL SCOPY( MNMIN-1, BE, 1, WORK( IWBE ), 1 )
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SBDSVDX( UPLO, 'N', 'A', MNMIN, WORK( IWBD ), &
                       WORK( IWBE ), 0.0E+0, 0.0E+0, 0, 0, NS2, S2, &
                       WORK( IWBZ ), MNMIN2, WORK( IWWORK ), &
                       IWORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SBDSVDX : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
!           Check error code from SBDSVDX.
!
         IF( IINFO /= 0 ) THEN
            WRITE( NOUT, FMT = 9998 )'SBDSVDX(values,A)', IINFO, &
               M, N, JTYPE, IOLDSD
            INFO = ABS( IINFO )
            IF( IINFO < 0 ) THEN
               RETURN
            ELSE
               RESULT( 24 ) = ULPINV
               GO TO 270
            END IF
         END IF
!
!           Save S1 for tests 30-34.
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SCOPY( MNMIN, S1, 1, WORK( IWBS ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
!           Test 20:  Check the decomposition B := U * S1 * VT
!                21:  Check the orthogonality of U
!                22:  Check the orthogonality of VT
!                23:  Check that the singular values are sorted in
!                     non-increasing order and are non-negative
!                24:  Compare SBDSVDX with and without singular vectors
!
         CALL SBDT03( UPLO, MNMIN, 1, BD, BE, U, LDPT, S1, VT, &
                      LDPT, WORK( IWBS+MNMIN ), RESULT( 20 ) )
         CALL SORT01( 'Columns', MNMIN, MNMIN, U, LDPT, &
                      WORK( IWBS+MNMIN ), LWORK-MNMIN, &
                      RESULT( 21 ) )
         CALL SORT01( 'Rows', MNMIN, MNMIN, VT, LDPT, &
                      WORK( IWBS+MNMIN ), LWORK-MNMIN, &
                      RESULT( 22) )
!
         RESULT( 23 ) = 0.0E+0
         DO I = 1, MNMIN - 1
            IF( S1( I ) < S1( I+1 ) ) RESULT( 23 ) = ULPINV
            IF( S1( I ) < 0.0E+0 ) RESULT( 23 ) = ULPINV
         ENDDO
         IF( MNMIN >= 1 ) THEN
            IF( S1( MNMIN ) < 0.0E+0 ) RESULT( 23 ) = ULPINV
         END IF
!
         TEMP2 = 0.0E+0
         DO J = 1, MNMIN
            TEMP1 = ABS( S1( J )-S2( J ) ) / MAX( SQRT( UNFL )*MAX( S1( 1 ), 1.0E+0 ), &
                    ULP*MAX( ABS( S1( 1 ) ), ABS( S2( 1 ) ) ) )
            TEMP2 = MAX( TEMP1, TEMP2 )
         ENDDO
         RESULT( 24 ) = TEMP2
         ANORM = S1( 1 )
!
!           Use SBDSVDX with RANGE='I': choose random values for IL and
!           IU, and ask for the IL-th through IU-th singular values
!           and corresponding vectors.
!
         ISEED2(1:4) = ISEED(1:4)
         IF( MNMIN <= 1 ) THEN
            IL = 1
            IU = MNMIN
         ELSE
            IL = 1 + INT( ( MNMIN-1 )*SLARND( 1, ISEED2 ) )
            IU = 1 + INT( ( MNMIN-1 )*SLARND( 1, ISEED2 ) )
            IF( IU < IL ) THEN
               ITEMP = IU
               IU = IL
               IL = ITEMP
            END IF
         END IF
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SCOPY( MNMIN, BD, 1, WORK( IWBD ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( MNMIN > 0 )  THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SCOPY( MNMIN-1, BE, 1, WORK( IWBE ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ENDIF
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SBDSVDX( UPLO, 'V', 'I', MNMIN, WORK( IWBD ), &
                       WORK( IWBE ), 0.0E+0, 0.0E+0, IL, IU, NS1, S1, &
                       WORK( IWBZ ), MNMIN2, WORK( IWWORK ), &
                       IWORK, IINFO)
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SBDSVDX : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
!           Check error code from SBDSVDX.
!
         IF( IINFO /= 0 ) THEN
            WRITE( NOUT, FMT = 9998 )'SBDSVDX(vects,I)', IINFO, &
               M, N, JTYPE, IOLDSD
            INFO = ABS( IINFO )
            IF( IINFO < 0 ) THEN
               RETURN
            ELSE
               RESULT( 25 ) = ULPINV
               GO TO 270
            END IF
         END IF
!
         J = IWBZ
         DO I = 1, NS1
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SCOPY( MNMIN, WORK( J ), 1, U( 1,I ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            J = J + MNMIN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SCOPY( MNMIN, WORK( J ), 1, VT( I,1 ), LDPT )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            J = J + MNMIN
         ENDDO
!
!           Use SBDSVDX to compute only the singular values of the
!           bidiagonal matrix B;  U and VT should not be modified.
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SCOPY( MNMIN, BD, 1, WORK( IWBD ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( MNMIN > 0 ) CALL SCOPY( MNMIN-1, BE, 1, WORK( IWBE ), 1 )
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SBDSVDX( UPLO, 'N', 'I', MNMIN, WORK( IWBD ), &
                       WORK( IWBE ), 0.0E+0, 0.0E+0, IL, IU, NS2, S2, &
                       WORK( IWBZ ), MNMIN2, WORK( IWWORK ), &
                       IWORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SBDSVDX : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
!           Check error code from SBDSVDX.
!
         IF( IINFO /= 0 ) THEN
            WRITE( NOUT, FMT = 9998 )'SBDSVDX(values,I)', IINFO, &
               M, N, JTYPE, IOLDSD
            INFO = ABS( IINFO )
            IF( IINFO < 0 ) THEN
               RETURN
            ELSE
               RESULT( 29 ) = ULPINV
               GO TO 270
            END IF
         END IF
!
!           Test 25:  Check S1 - U' * B * VT'
!                26:  Check the orthogonality of U
!                27:  Check the orthogonality of VT
!                28:  Check that the singular values are sorted in
!                     non-increasing order and are non-negative
!                29:  Compare SBDSVDX with and without singular vectors
!
         CALL SBDT04( UPLO, MNMIN, BD, BE, S1, NS1, U, LDPT, VT, LDPT, WORK( IWBS+MNMIN ), &
                      RESULT( 25 ) )
         CALL SORT01( 'Columns', MNMIN, NS1, U, LDPT, WORK( IWBS+MNMIN ), LWORK-MNMIN, &
                      RESULT( 26 ) )
         CALL SORT01( 'Rows', NS1, MNMIN, VT, LDPT,  WORK( IWBS+MNMIN ), LWORK-MNMIN, &
                      RESULT( 27 ) )
!
         RESULT( 28 ) = 0.0E+0
         DO I = 1, NS1 - 1
            IF( S1( I ) < S1( I+1 ) ) &
               RESULT( 28 ) = ULPINV
            IF( S1( I ) < 0.0E+0 ) &
               RESULT( 28 ) = ULPINV
            ENDDO
         IF( NS1 >= 1 ) THEN
            IF( S1( NS1 ) < 0.0E+0 ) RESULT( 28 ) = ULPINV
         END IF
!
         TEMP2 = 0.0E+0
         DO J = 1, NS1
            TEMP1 = ABS( S1( J )-S2( J ) ) / &
                    MAX( SQRT( UNFL )*MAX( S1( 1 ), 1.0E+0 ), &
                    ULP*MAX( ABS( S1( 1 ) ), ABS( S2( 1 ) ) ) )
            TEMP2 = MAX( TEMP1, TEMP2 )
            ENDDO
         RESULT( 29 ) = TEMP2
!
!           Use SBDSVDX with RANGE='V': determine the values VL and VU
!           of the IL-th and IU-th singular values and ask for all
!           singular values in this range.
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SCOPY( MNMIN, WORK( IWBS ), 1, S1, 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
         IF( MNMIN > 0 ) THEN
            IF( IL /= 1 ) THEN
               VU = S1( IL ) + MAX( 0.5E+0*ABS( S1( IL )-S1( IL-1 ) ), &
                    ULP*ANORM, 2.0E+0*RTUNFL )
            ELSE
               VU = S1( 1 ) + MAX( 0.5E+0*ABS( S1( MNMIN )-S1( 1 ) ), &
                    ULP*ANORM, 2.0E+0*RTUNFL )
            END IF
            IF( IU /= NS1 ) THEN
               VL = S1( IU ) - MAX( ULP*ANORM, 2.0E+0*RTUNFL, &
                    0.5E+0*ABS( S1( IU+1 )-S1( IU ) ) )
            ELSE
               VL = S1( NS1 ) - MAX( ULP*ANORM, 2.0E+0*RTUNFL, &
                    0.5E+0*ABS( S1( MNMIN )-S1( 1 ) ) )
            END IF
            VL = MAX( VL,0.0E+0 )
            VU = MAX( VU,0.0E+0 )
            IF( VL >= VU ) VU = MAX( VU*2, VU+VL+0.5E+0 )
         ELSE
            VL = 0.0E+0
            VU = 1.0E+0
         END IF
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SCOPY( MNMIN, BD, 1, WORK( IWBD ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( MNMIN > 0 )  THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SCOPY( MNMIN-1, BE, 1, WORK( IWBE ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ENDIF
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SBDSVDX( UPLO, 'V', 'V', MNMIN, WORK( IWBD ), &
                       WORK( IWBE ), VL, VU, 0, 0, NS1, S1, &
                       WORK( IWBZ ), MNMIN2, WORK( IWWORK ), &
                       IWORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SBDSVDX : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
!           Check error code from SBDSVDX.
!
         IF( IINFO /= 0 ) THEN
            WRITE( NOUT, FMT = 9998 )'SBDSVDX(vects,V)', IINFO, M, N, JTYPE, IOLDSD
            INFO = ABS( IINFO )
            IF( IINFO < 0 ) THEN
               RETURN
            ELSE
               RESULT( 30 ) = ULPINV
               GO TO 270
            END IF
         END IF
!
         J = IWBZ
         DO I = 1, NS1
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SCOPY( MNMIN, WORK( J ), 1, U( 1,I ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            J = J + MNMIN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SCOPY( MNMIN, WORK( J ), 1, VT( I,1 ), LDPT )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            J = J + MNMIN
         ENDDO
!
!           Use SBDSVDX to compute only the singular values of the
!           bidiagonal matrix B;  U and VT should not be modified.
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SCOPY( MNMIN, BD, 1, WORK( IWBD ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( MNMIN > 0 )  THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SCOPY( MNMIN-1, BE, 1, WORK( IWBE ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ENDIF
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SBDSVDX( UPLO, 'N', 'V', MNMIN, WORK( IWBD ), &
                       WORK( IWBE ), VL, VU, 0, 0, NS2, S2, &
                       WORK( IWBZ ), MNMIN2, WORK( IWWORK ), &
                       IWORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SBDSVDX : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
!           Check error code from SBDSVDX.
!
         IF( IINFO /= 0 ) THEN
            WRITE( NOUT, FMT = 9998 )'SBDSVDX(values,V)', IINFO, &
               M, N, JTYPE, IOLDSD
            INFO = ABS( IINFO )
            IF( IINFO < 0 ) THEN
               RETURN
            ELSE
               RESULT( 34 ) = ULPINV
               GO TO 270
            END IF
         END IF
!
!           Test 30:  Check S1 - U' * B * VT'
!                31:  Check the orthogonality of U
!                32:  Check the orthogonality of VT
!                33:  Check that the singular values are sorted in
!                     non-increasing order and are non-negative
!                34:  Compare SBDSVDX with and without singular vectors
!
         CALL SBDT04( UPLO, MNMIN, BD, BE, S1, NS1, U, &
                      LDPT, VT, LDPT, WORK( IWBS+MNMIN ), &
                      RESULT( 30 ) )
         CALL SORT01( 'Columns', MNMIN, NS1, U, LDPT, &
                      WORK( IWBS+MNMIN ), LWORK-MNMIN, &
                      RESULT( 31 ) )
         CALL SORT01( 'Rows', NS1, MNMIN, VT, LDPT, &
                      WORK( IWBS+MNMIN ), LWORK-MNMIN, &
                      RESULT( 32 ) )
!
         RESULT( 33 ) = 0.0E+0
         DO I = 1, NS1 - 1
            IF( S1( I ) < S1( I+1 ) ) RESULT( 28 ) = ULPINV
            IF( S1( I ) < 0.0E+0 ) RESULT( 28 ) = ULPINV
            ENDDO
         IF( NS1 >= 1 ) THEN
            IF( S1( NS1 ) < 0.0E+0 ) RESULT( 28 ) = ULPINV
         END IF
!
         TEMP2 = 0.0E+0
         DO J = 1, NS1
            TEMP1 = ABS( S1( J )-S2( J ) ) / &
                    MAX( SQRT( UNFL )*MAX( S1( 1 ), 1.0E+0 ), &
                    ULP*MAX( ABS( S1( 1 ) ), ABS( S2( 1 ) ) ) )
            TEMP2 = MAX( TEMP1, TEMP2 )
            ENDDO
         RESULT( 34 ) = TEMP2
!
!           End of Loop -- Check for RESULT(j) > THRESH
!
  270       CONTINUE
!
         DO J = 1, 34
            IF( RESULT( J ) >= THRESH ) THEN
               IF( NFAIL == 0 ) CALL SLAHD2( NOUT, PATH )
               WRITE( NOUT, FMT = 9999 )M, N, JTYPE, IOLDSD, J, RESULT( J )
               NFAIL = NFAIL + 1
            END IF
         ENDDO
         IF( .NOT.BIDIAG ) THEN
            NTEST = NTEST + 34
         ELSE
            NTEST = NTEST + 30
         END IF
!
         ENDIF
         ENDDO
      ENDDO
!
!     Summary
!
   CALL ALASUM( PATH, NOUT, NFAIL, NTEST, 0 )
!
   RETURN
!
!     End of SCHKBD
!
 9999 FORMAT( ' M=', I5, ', N=', I5, ', type ', I2, ', seed=', &
         4( I4, ',' ), ' test(', I2, ')=', G11.4 )
 9998 FORMAT( ' SCHKBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=', &
         I6, ', N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), &
         I5, ')' )
!
END



