!> \brief \b SDRVST2STG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SDRVST2STG( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                          NOUNIT, A, LDA, D1, D2, D3, D4, EVEIGS, WA1,
!                          WA2, WA3, U, LDU, V, TAU, Z, WORK, LWORK,
!                          IWORK, LIWORK, RESULT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDU, LIWORK, LWORK, NOUNIT, NSIZES,
!      $                   NTYPES
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
!       REAL               A( LDA, * ), D1( * ), D2( * ), D3( * ),
!      $                   D4( * ), EVEIGS( * ), RESULT( * ), TAU( * ),
!      $                   U( LDU, * ), V( LDU, * ), WA1( * ), WA2( * ),
!      $                   WA3( * ), WORK( * ), Z( LDU, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      SDRVST2STG  checks the symmetric eigenvalue problem drivers.
!>
!>              SSTEV computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric tridiagonal matrix.
!>
!>              SSTEVX computes selected eigenvalues and, optionally,
!>              eigenvectors of a real symmetric tridiagonal matrix.
!>
!>              SSTEVR computes selected eigenvalues and, optionally,
!>              eigenvectors of a real symmetric tridiagonal matrix
!>              using the Relatively Robust Representation where it can.
!>
!>              SSYEV computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric matrix.
!>
!>              SSYEVX computes selected eigenvalues and, optionally,
!>              eigenvectors of a real symmetric matrix.
!>
!>              SSYEVR computes selected eigenvalues and, optionally,
!>              eigenvectors of a real symmetric matrix
!>              using the Relatively Robust Representation where it can.
!>
!>              SSPEV computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric matrix in packed
!>              storage.
!>
!>              SSPEVX computes selected eigenvalues and, optionally,
!>              eigenvectors of a real symmetric matrix in packed
!>              storage.
!>
!>              SSBEV computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric band matrix.
!>
!>              SSBEVX computes selected eigenvalues and, optionally,
!>              eigenvectors of a real symmetric band matrix.
!>
!>              SSYEVD computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric matrix using
!>              a divide and conquer algorithm.
!>
!>              SSPEVD computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric matrix in packed
!>              storage, using a divide and conquer algorithm.
!>
!>              SSBEVD computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric band matrix,
!>              using a divide and conquer algorithm.
!>
!>      When SDRVST2STG is called, a number of matrix "sizes" ("n's") and a
!>      number of matrix "types" are specified.  For each size ("n")
!>      and each type of matrix, one matrix will be generated and used
!>      to test the appropriate drivers.  For each matrix and each
!>      driver routine called, the following tests will be performed:
!>
!>      (1)     | A - Z D Z' | / ( |A| n ulp )
!>
!>      (2)     | I - Z Z' | / ( n ulp )
!>
!>      (3)     | D1 - D2 | / ( |D1| ulp )
!>
!>      where Z is the matrix of eigenvectors returned when the
!>      eigenvector option is given and D1 and D2 are the eigenvalues
!>      returned with and without the eigenvector option.
!>
!>      The "sizes" are specified by an array NN(1:NSIZES); the value of
!>      each element NN(j) specifies one size.
!>      The "types" are specified by a logical array DOTYPE( 1:NTYPES );
!>      if DOTYPE(j) is .TRUE., then matrix type "j" will be generated.
!>      Currently, the list of possible types is:
!>
!>      (1)  The zero matrix.
!>      (2)  The identity matrix.
!>
!>      (3)  A diagonal matrix with evenly spaced eigenvalues
!>           1, ..., ULP  and random signs.
!>           (ULP = (first number larger than 1) - 1 )
!>      (4)  A diagonal matrix with geometrically spaced eigenvalues
!>           1, ..., ULP  and random signs.
!>      (5)  A diagonal matrix with "clustered" eigenvalues
!>           1, ULP, ..., ULP and random signs.
!>
!>      (6)  Same as (4), but multiplied by SQRT( overflow threshold )
!>      (7)  Same as (4), but multiplied by SQRT( underflow threshold )
!>
!>      (8)  A matrix of the form  U' D U, where U is orthogonal and
!>           D has evenly spaced entries 1, ..., ULP with random signs
!>           on the diagonal.
!>
!>      (9)  A matrix of the form  U' D U, where U is orthogonal and
!>           D has geometrically spaced entries 1, ..., ULP with random
!>           signs on the diagonal.
!>
!>      (10) A matrix of the form  U' D U, where U is orthogonal and
!>           D has "clustered" entries 1, ULP,..., ULP with random
!>           signs on the diagonal.
!>
!>      (11) Same as (8), but multiplied by SQRT( overflow threshold )
!>      (12) Same as (8), but multiplied by SQRT( underflow threshold )
!>
!>      (13) Symmetric matrix with random entries chosen from (-1,1).
!>      (14) Same as (13), but multiplied by SQRT( overflow threshold )
!>      (15) Same as (13), but multiplied by SQRT( underflow threshold )
!>      (16) A band matrix with half bandwidth randomly chosen between
!>           0 and N-1, with evenly spaced eigenvalues 1, ..., ULP
!>           with random signs.
!>      (17) Same as (16), but multiplied by SQRT( overflow threshold )
!>      (18) Same as (16), but multiplied by SQRT( underflow threshold )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \verbatim
!>  NSIZES  INTEGER
!>          The number of sizes of matrices to use.  If it is zero,
!>          SDRVST2STG does nothing.  It must be at least zero.
!>          Not modified.
!>
!>  NN      INTEGER array, dimension (NSIZES)
!>          An array containing the sizes to be used for the matrices.
!>          Zero values will be skipped.  The values must be at least
!>          zero.
!>          Not modified.
!>
!>  NTYPES  INTEGER
!>          The number of elements in DOTYPE.   If it is zero, SDRVST2STG
!>          does nothing.  It must be at least zero.  If it is MAXTYP+1
!>          and NSIZES is 1, then an additional type, MAXTYP+1 is
!>          defined, which is to use whatever matrix is in A.  This
!>          is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
!>          DOTYPE(MAXTYP+1) is .TRUE. .
!>          Not modified.
!>
!>  DOTYPE  LOGICAL array, dimension (NTYPES)
!>          If DOTYPE(j) is .TRUE., then for each size in NN a
!>          matrix of that size and of type j will be generated.
!>          If NTYPES is smaller than the maximum number of types
!>          defined (PARAMETER MAXTYP), then types NTYPES+1 through
!>          MAXTYP will not be generated.  If NTYPES is larger
!>          than MAXTYP, DOTYPE(MAXTYP+1) through DOTYPE(NTYPES)
!>          will be ignored.
!>          Not modified.
!>
!>  ISEED   INTEGER array, dimension (4)
!>          On entry ISEED specifies the seed of the random number
!>          generator. The array elements should be between 0 and 4095;
!>          if not they will be reduced mod 4096.  Also, ISEED(4) must
!>          be odd.  The random number generator uses a linear
!>          congruential sequence limited to small integers, and so
!>          should produce machine independent random numbers. The
!>          values of ISEED are changed on exit, and can be used in the
!>          next call to SDRVST2STG to continue the same random number
!>          sequence.
!>          Modified.
!>
!>  THRESH  REAL
!>          A test will count as "failed" if the "error", computed as
!>          described above, exceeds THRESH.  Note that the error
!>          is scaled to be O(1), so THRESH should be a reasonably
!>          small multiple of 1, e.g., 10 or 100.  In particular,
!>          it should not depend on the precision (single vs. double)
!>          or the size of the matrix.  It must be at least zero.
!>          Not modified.
!>
!>  NOUNIT  INTEGER
!>          The FORTRAN unit number for printing out error messages
!>          (e.g., if a routine returns IINFO not equal to 0.)
!>          Not modified.
!>
!>  A       REAL             array, dimension (LDA , max(NN))
!>          Used to hold the matrix whose eigenvalues are to be
!>          computed.  On exit, A contains the last matrix actually
!>          used.
!>          Modified.
!>
!>  LDA     INTEGER
!>          The leading dimension of A.  It must be at
!>          least 1 and at least max( NN ).
!>          Not modified.
!>
!>  D1      REAL             array, dimension (max(NN))
!>          The eigenvalues of A, as computed by SSTEQR simultaneously
!>          with Z.  On exit, the eigenvalues in D1 correspond with the
!>          matrix in A.
!>          Modified.
!>
!>  D2      REAL             array, dimension (max(NN))
!>          The eigenvalues of A, as computed by SSTEQR if Z is not
!>          computed.  On exit, the eigenvalues in D2 correspond with
!>          the matrix in A.
!>          Modified.
!>
!>  D3      REAL             array, dimension (max(NN))
!>          The eigenvalues of A, as computed by SSTERF.  On exit, the
!>          eigenvalues in D3 correspond with the matrix in A.
!>          Modified.
!>
!>  D4      REAL             array, dimension
!>
!>  EVEIGS  REAL array, dimension (max(NN))
!>          The eigenvalues as computed by SSTEV('N', ... )
!>          (I reserve the right to change this to the output of
!>          whichever algorithm computes the most accurate eigenvalues).
!>
!>  WA1     REAL array, dimension
!>
!>  WA2     REAL array, dimension
!>
!>  WA3     REAL array, dimension
!>
!>  U       REAL             array, dimension (LDU, max(NN))
!>          The orthogonal matrix computed by SSYTRD + SORGTR.
!>          Modified.
!>
!>  LDU     INTEGER
!>          The leading dimension of U, Z, and V.  It must be at
!>          least 1 and at least max( NN ).
!>          Not modified.
!>
!>  V       REAL             array, dimension (LDU, max(NN))
!>          The Housholder vectors computed by SSYTRD in reducing A to
!>          tridiagonal form.
!>          Modified.
!>
!>  TAU     REAL array, dimension (max(NN))
!>          The Householder factors computed by SSYTRD in reducing A
!>          to tridiagonal form.
!>          Modified.
!>
!>  Z       REAL             array, dimension (LDU, max(NN))
!>          The orthogonal matrix of eigenvectors computed by SSTEQR,
!>          SPTEQR, and SSTEIN.
!>          Modified.
!>
!>  WORK    REAL array, dimension (LWORK)
!>          Workspace.
!>          Modified.
!>
!>  LWORK   INTEGER
!>          The number of entries in WORK.  This must be at least
!>          1 + 4 * Nmax + 2 * Nmax * lg Nmax + 4 * Nmax**2
!>          where Nmax = max( NN(j), 2 ) and lg = log base 2.
!>          Not modified.
!>
!>  IWORK   INTEGER array,
!>             dimension (6 + 6*Nmax + 5 * Nmax * lg Nmax )
!>          where Nmax = max( NN(j), 2 ) and lg = log base 2.
!>          Workspace.
!>          Modified.
!>
!>  RESULT  REAL array, dimension (105)
!>          The values computed by the tests described above.
!>          The values are currently limited to 1/ulp, to avoid
!>          overflow.
!>          Modified.
!>
!>  INFO    INTEGER
!>          If 0, then everything ran OK.
!>           -1: NSIZES < 0
!>           -2: Some NN(j) < 0
!>           -3: NTYPES < 0
!>           -5: THRESH < 0
!>           -9: LDA < 1 or LDA < NMAX, where NMAX is max( NN(j) ).
!>          -16: LDU < 1 or LDU < NMAX.
!>          -21: LWORK too small.
!>          If  SLATMR, SLATMS, SSYTRD, SORGTR, SSTEQR, SSTERF,
!>              or SORMTR returns an error code, the
!>              absolute value of it is returned.
!>          Modified.
!>
!>-----------------------------------------------------------------------
!>
!>       Some Local Variables and Parameters:
!>       ---- ----- --------- --- ----------
!>       0.0E+0, 1.0E+0       Real 0 and 1.
!>       MAXTYP          The number of types defined.
!>       NTEST           The number of tests performed, or which can
!>                       be performed so far, for the current matrix.
!>       NTESTT          The total number of tests performed so far.
!>       NMAX            Largest value in NN.
!>       NMATS           The number of matrices generated so far.
!>       NERRS           The number of tests which have exceeded THRESH
!>                       so far (computed by SLAFTS).
!>       COND, IMODE     Values to be passed to the matrix generators.
!>       ANORM           Norm of A; passed to matrix generators.
!>
!>       OVFL, UNFL      Overflow and underflow thresholds.
!>       ULP, ULPINV     Finest relative precision and its inverse.
!>       RTOVFL, RTUNFL  Square roots of the previous 2 values.
!>               The following four arrays decode JTYPE:
!>       KTYPE(j)        The general type (1-10) for type "j".
!>       KMODE(j)        The MODE value to be passed to the matrix
!>                       generator for type "j".
!>       KMAGN(j)        The order of magnitude ( O(1),
!>                       O(overflow^(1/2) ), O(underflow^(1/2) )
!>
!>     The tests performed are:                 Routine tested
!>    1= | A - U S U' | / ( |A| n ulp )         SSTEV('V', ... )
!>    2= | I - U U' | / ( n ulp )               SSTEV('V', ... )
!>    3= |D(with Z) - D(w/o Z)| / (|D| ulp)     SSTEV('N', ... )
!>    4= | A - U S U' | / ( |A| n ulp )         SSTEVX('V','A', ... )
!>    5= | I - U U' | / ( n ulp )               SSTEVX('V','A', ... )
!>    6= |D(with Z) - EVEIGS| / (|D| ulp)       SSTEVX('N','A', ... )
!>    7= | A - U S U' | / ( |A| n ulp )         SSTEVR('V','A', ... )
!>    8= | I - U U' | / ( n ulp )               SSTEVR('V','A', ... )
!>    9= |D(with Z) - EVEIGS| / (|D| ulp)       SSTEVR('N','A', ... )
!>    10= | A - U S U' | / ( |A| n ulp )        SSTEVX('V','I', ... )
!>    11= | I - U U' | / ( n ulp )              SSTEVX('V','I', ... )
!>    12= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSTEVX('N','I', ... )
!>    13= | A - U S U' | / ( |A| n ulp )        SSTEVX('V','V', ... )
!>    14= | I - U U' | / ( n ulp )              SSTEVX('V','V', ... )
!>    15= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSTEVX('N','V', ... )
!>    16= | A - U S U' | / ( |A| n ulp )        SSTEVD('V', ... )
!>    17= | I - U U' | / ( n ulp )              SSTEVD('V', ... )
!>    18= |D(with Z) - EVEIGS| / (|D| ulp)      SSTEVD('N', ... )
!>    19= | A - U S U' | / ( |A| n ulp )        SSTEVR('V','I', ... )
!>    20= | I - U U' | / ( n ulp )              SSTEVR('V','I', ... )
!>    21= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSTEVR('N','I', ... )
!>    22= | A - U S U' | / ( |A| n ulp )        SSTEVR('V','V', ... )
!>    23= | I - U U' | / ( n ulp )              SSTEVR('V','V', ... )
!>    24= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSTEVR('N','V', ... )
!>
!>    25= | A - U S U' | / ( |A| n ulp )        SSYEV('L','V', ... )
!>    26= | I - U U' | / ( n ulp )              SSYEV('L','V', ... )
!>    27= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSYEV_2STAGE('L','N', ... )
!>    28= | A - U S U' | / ( |A| n ulp )        SSYEVX('L','V','A', ... )
!>    29= | I - U U' | / ( n ulp )              SSYEVX('L','V','A', ... )
!>    30= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSYEVX_2STAGE('L','N','A', ... )
!>    31= | A - U S U' | / ( |A| n ulp )        SSYEVX('L','V','I', ... )
!>    32= | I - U U' | / ( n ulp )              SSYEVX('L','V','I', ... )
!>    33= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSYEVX_2STAGE('L','N','I', ... )
!>    34= | A - U S U' | / ( |A| n ulp )        SSYEVX('L','V','V', ... )
!>    35= | I - U U' | / ( n ulp )              SSYEVX('L','V','V', ... )
!>    36= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSYEVX_2STAGE('L','N','V', ... )
!>    37= | A - U S U' | / ( |A| n ulp )        SSPEV('L','V', ... )
!>    38= | I - U U' | / ( n ulp )              SSPEV('L','V', ... )
!>    39= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSPEV('L','N', ... )
!>    40= | A - U S U' | / ( |A| n ulp )        SSPEVX('L','V','A', ... )
!>    41= | I - U U' | / ( n ulp )              SSPEVX('L','V','A', ... )
!>    42= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSPEVX('L','N','A', ... )
!>    43= | A - U S U' | / ( |A| n ulp )        SSPEVX('L','V','I', ... )
!>    44= | I - U U' | / ( n ulp )              SSPEVX('L','V','I', ... )
!>    45= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSPEVX('L','N','I', ... )
!>    46= | A - U S U' | / ( |A| n ulp )        SSPEVX('L','V','V', ... )
!>    47= | I - U U' | / ( n ulp )              SSPEVX('L','V','V', ... )
!>    48= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSPEVX('L','N','V', ... )
!>    49= | A - U S U' | / ( |A| n ulp )        SSBEV('L','V', ... )
!>    50= | I - U U' | / ( n ulp )              SSBEV('L','V', ... )
!>    51= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSBEV_2STAGE('L','N', ... )
!>    52= | A - U S U' | / ( |A| n ulp )        SSBEVX('L','V','A', ... )
!>    53= | I - U U' | / ( n ulp )              SSBEVX('L','V','A', ... )
!>    54= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSBEVX_2STAGE('L','N','A', ... )
!>    55= | A - U S U' | / ( |A| n ulp )        SSBEVX('L','V','I', ... )
!>    56= | I - U U' | / ( n ulp )              SSBEVX('L','V','I', ... )
!>    57= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSBEVX_2STAGE('L','N','I', ... )
!>    58= | A - U S U' | / ( |A| n ulp )        SSBEVX('L','V','V', ... )
!>    59= | I - U U' | / ( n ulp )              SSBEVX('L','V','V', ... )
!>    60= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSBEVX_2STAGE('L','N','V', ... )
!>    61= | A - U S U' | / ( |A| n ulp )        SSYEVD('L','V', ... )
!>    62= | I - U U' | / ( n ulp )              SSYEVD('L','V', ... )
!>    63= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSYEVD_2STAGE('L','N', ... )
!>    64= | A - U S U' | / ( |A| n ulp )        SSPEVD('L','V', ... )
!>    65= | I - U U' | / ( n ulp )              SSPEVD('L','V', ... )
!>    66= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSPEVD('L','N', ... )
!>    67= | A - U S U' | / ( |A| n ulp )        SSBEVD('L','V', ... )
!>    68= | I - U U' | / ( n ulp )              SSBEVD('L','V', ... )
!>    69= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSBEVD_2STAGE('L','N', ... )
!>    70= | A - U S U' | / ( |A| n ulp )        SSYEVR('L','V','A', ... )
!>    71= | I - U U' | / ( n ulp )              SSYEVR('L','V','A', ... )
!>    72= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSYEVR_2STAGE('L','N','A', ... )
!>    73= | A - U S U' | / ( |A| n ulp )        SSYEVR('L','V','I', ... )
!>    74= | I - U U' | / ( n ulp )              SSYEVR('L','V','I', ... )
!>    75= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSYEVR_2STAGE('L','N','I', ... )
!>    76= | A - U S U' | / ( |A| n ulp )        SSYEVR('L','V','V', ... )
!>    77= | I - U U' | / ( n ulp )              SSYEVR('L','V','V', ... )
!>    78= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSYEVR_2STAGE('L','N','V', ... )
!>
!>    Tests 25 through 78 are repeated (as tests 79 through 132)
!>    with UPLO='U'
!>
!>    To be added in 1999
!>
!>    79= | A - U S U' | / ( |A| n ulp )        SSPEVR('L','V','A', ... )
!>    80= | I - U U' | / ( n ulp )              SSPEVR('L','V','A', ... )
!>    81= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSPEVR('L','N','A', ... )
!>    82= | A - U S U' | / ( |A| n ulp )        SSPEVR('L','V','I', ... )
!>    83= | I - U U' | / ( n ulp )              SSPEVR('L','V','I', ... )
!>    84= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSPEVR('L','N','I', ... )
!>    85= | A - U S U' | / ( |A| n ulp )        SSPEVR('L','V','V', ... )
!>    86= | I - U U' | / ( n ulp )              SSPEVR('L','V','V', ... )
!>    87= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSPEVR('L','N','V', ... )
!>    88= | A - U S U' | / ( |A| n ulp )        SSBEVR('L','V','A', ... )
!>    89= | I - U U' | / ( n ulp )              SSBEVR('L','V','A', ... )
!>    90= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSBEVR('L','N','A', ... )
!>    91= | A - U S U' | / ( |A| n ulp )        SSBEVR('L','V','I', ... )
!>    92= | I - U U' | / ( n ulp )              SSBEVR('L','V','I', ... )
!>    93= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSBEVR('L','N','I', ... )
!>    94= | A - U S U' | / ( |A| n ulp )        SSBEVR('L','V','V', ... )
!>    95= | I - U U' | / ( n ulp )              SSBEVR('L','V','V', ... )
!>    96= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSBEVR('L','N','V', ... )
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
!> \ingroup single_eig
!
!  =====================================================================
   SUBROUTINE SDRVST2STG( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, &
                      NOUNIT, A, LDA, D1, D2, D3, D4, EVEIGS, WA1, &
                      WA2, WA3, U, LDU, V, TAU, Z, WORK, LWORK, &
                      IWORK, LIWORK, RESULT, INFO )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, LDA, LDU, LIWORK, LWORK, NOUNIT, NSIZES, &
                      NTYPES
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
   REAL               A( LDA, * ), D1( * ), D2( * ), D3( * ), &
                      D4( * ), EVEIGS( * ), RESULT( * ), TAU( * ), &
                      U( LDU, * ), V( LDU, * ), WA1( * ), WA2( * ), &
                      WA3( * ), WORK( * ), Z( LDU, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            MAXTYP
   PARAMETER          ( MAXTYP = 18 )
!     ..
!     .. Local Scalars ..
   LOGICAL            BADNN
   CHARACTER          UPLO
   INTEGER            I, IDIAG, IHBW, IINFO, IL, IMODE, INDX, IROW, &
                      ITEMP, ITYPE, IU, IUPLO, J, J1, J2, JCOL, &
                      JSIZE, JTYPE, KD, LGN, LIWEDC, LWEDC, M, M2, &
                      M3, MTYPES, N, NERRS, NMATS, NMAX, NTEST, &
                      NTESTT
   REAL               ABSTOL, ANINV, ANORM, COND, OVFL, RTOVFL, &
                      RTUNFL, TEMP1, TEMP2, TEMP3, ULP, ULPINV, UNFL, &
                      VL, VU
!     ..
!     .. Local Arrays ..
   INTEGER            IDUMMA( 1 ), IOLDSD( 4 ), ISEED2( 4 ), &
                      ISEED3( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), &
                      KTYPE( MAXTYP )
!     ..
!     .. External Functions ..
   REAL               SLAMCH, SLARND, SSXT1
   EXTERNAL           SLAMCH, SLARND, SSXT1
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALASVM, SLACPY, SLAFTS, SLASET, SLATMR, &
                      SLATMS, SSBEV, SSBEVD, SSBEVX, SSPEV, SSPEVD, &
                      SSPEVX, SSTEV, SSTEVD, SSTEVR, SSTEVX, SSTT21, &
                      SSTT22, SSYEV, SSYEVD, SSYEVR, SSYEVX, SSYT21, &
                      SSYEVD_2STAGE, SSYEVR_2STAGE, SSYEVX_2STAGE, &
                      SSYEV_2STAGE, SSBEV_2STAGE, SSBEVD_2STAGE, &
                      SSBEVX_2STAGE, SSYTRD_2STAGE, SSYTRD_SY2SB, &
                      SSYTRD_SB2ST, SSYT22, XERBLA
!     ..
!     .. Scalars in Common ..
   CHARACTER*32       SRNAMT
!     ..
!     .. Common blocks ..
   COMMON             / SRNAMC / SRNAMT
!     ..
!     .. Data statements ..
   DATA               KTYPE / 1, 2, 5*4, 5*5, 3*8, 3*9 /
   DATA               KMAGN / 2*1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 2, 3 /
   DATA               KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 4, 4 /
!     ..
!     .. Executable Statements ..
!
!     Keep ftrnchek happy
!
   VL = 0.0E+0
   VU = 0.0E+0
!
!     1)      Check for errors
!
   NTESTT = 0
   INFO = 0
!
   BADNN = ANY(NN(1:NSIZES) < 0)
   NMAX = MAXVAL(NN(1:NSIZES))
!
!     Check for errors
!
   IF( NSIZES < 0 ) THEN
      INFO = -1
   ELSE IF( BADNN ) THEN
      INFO = -2
   ELSE IF( NTYPES < 0 ) THEN
      INFO = -3
   ELSE IF( LDA < NMAX ) THEN
      INFO = -9
   ELSE IF( LDU < NMAX ) THEN
      INFO = -16
   ELSE IF( 2*MAX( 2, NMAX )**2 > LWORK ) THEN
      INFO = -21
   END IF
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'SDRVST2STG', -INFO )
      RETURN
   END IF
!
!     Quick return if nothing to do
!
   IF( NSIZES == 0 .OR. NTYPES == 0 ) RETURN
!
!     More Important constants
!
   UNFL = SLAMCH( 'Safe minimum' )
   OVFL = SLAMCH( 'Overflow' )
   ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
   ULPINV = 1.0E+0 / ULP
   RTUNFL = SQRT( UNFL )
   RTOVFL = SQRT( OVFL )
!
!     Loop over sizes, types
!
   ISEED2(1:4) = ISEED(1:4)
   ISEED3(1:4) = ISEED(1:4)
!
   NERRS = 0
   NMATS = 0
!
   DO JSIZE = 1, NSIZES
      N = NN( JSIZE )
      IF( N > 0 ) THEN
         LGN = INT( LOG( REAL( N ) ) / LOG( 2.0E+0 ) )
         IF( 2**LGN < N ) LGN = LGN + 1
         IF( 2**LGN < N ) LGN = LGN + 1
         LWEDC = 1 + 4*N + 2*N*LGN + 4*N**2
!           LIWEDC = 6 + 6*N + 5*N*LGN
         LIWEDC = 3 + 5*N
      ELSE
         LWEDC = 9
!           LIWEDC = 12
         LIWEDC = 8
      END IF
      ANINV = 1.0E+0 / REAL( MAX( 1, N ) )
!
      IF( NSIZES /= 1 ) THEN
         MTYPES = MIN( MAXTYP, NTYPES )
      ELSE
         MTYPES = MIN( MAXTYP+1, NTYPES )
      END IF
!
      DO JTYPE = 1, MTYPES
!
         IF( .NOT.DOTYPE( JTYPE ) ) GO TO 1730
         NMATS = NMATS + 1
         NTEST = 0
!
         IOLDSD(1:4) = ISEED(1:4)
!
!           2)      Compute "A"
!
!                   Control parameters:
!
!               KMAGN  KMODE        KTYPE
!           =1  O(1)   clustered 1  zero
!           =2  large  clustered 2  identity
!           =3  small  exponential  (none)
!           =4         arithmetic   diagonal, (w/ eigenvalues)
!           =5         random log   symmetric, w/ eigenvalues
!           =6         random       (none)
!           =7                      random diagonal
!           =8                      random symmetric
!           =9                      band symmetric, w/ eigenvalues
!
         IF( MTYPES > MAXTYP ) &
            GO TO 110
!
         ITYPE = KTYPE( JTYPE )
         IMODE = KMODE( JTYPE )
!
!           Compute norm
!
         SELECT CASE (KMAGN( JTYPE ))
          CASE(1)
           ANORM = 1.0E+0
          CASE(2)
           ANORM = ( RTOVFL*ULP )*ANINV
          CASE(3)
           ANORM = RTUNFL*N*ULPINV
         END SELECT
!
         CALL SLASET( 'Full', LDA, N, 0.0E+0, 0.0E+0, A, LDA )
         IINFO = 0
         COND = ULPINV
!
!           Special Matrices -- Identity & Jordan block
!
!                   Zero
!
         IF( ITYPE == 1 ) THEN
            IINFO = 0
!
         ELSE IF( ITYPE == 2 ) THEN
!
!              Identity
!
            FORALL (JCOL = 1:N) A( JCOL, JCOL ) = ANORM
!
         ELSE IF( ITYPE == 4 ) THEN
!
!              Diagonal Matrix, [Eigen]values Specified
!
            CALL SLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, &
                         ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ), &
                         IINFO )
!
         ELSE IF( ITYPE == 5 ) THEN
!
!              Symmetric, eigenvalues specified
!
            CALL SLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, &
                         ANORM, N, N, 'N', A, LDA, WORK( N+1 ), &
                         IINFO )
!
         ELSE IF( ITYPE == 7 ) THEN
!
!              Diagonal, random eigenvalues
!
            IDUMMA( 1 ) = 1
            CALL SLATMR( N, N, 'S', ISEED, 'S', WORK, 6, 1.0E+0, 1.0E+0, &
                         'T', 'N', WORK( N+1 ), 1, 1.0E+0, &
                         WORK( 2*N+1 ), 1, 1.0E+0, 'N', IDUMMA, 0, 0, &
                         0.0E+0, ANORM, 'NO', A, LDA, IWORK, IINFO )
!
         ELSE IF( ITYPE == 8 ) THEN
!
!              Symmetric, random eigenvalues
!
            IDUMMA( 1 ) = 1
            CALL SLATMR( N, N, 'S', ISEED, 'S', WORK, 6, 1.0E+0, 1.0E+0, &
                         'T', 'N', WORK( N+1 ), 1, 1.0E+0, &
                         WORK( 2*N+1 ), 1, 1.0E+0, 'N', IDUMMA, N, N, &
                         0.0E+0, ANORM, 'NO', A, LDA, IWORK, IINFO )
!
         ELSE IF( ITYPE == 9 ) THEN
!
!              Symmetric banded, eigenvalues specified
!
            IHBW = INT( ( N-1 )*SLARND( 1, ISEED3 ) )
            CALL SLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, &
                         ANORM, IHBW, IHBW, 'Z', U, LDU, WORK( N+1 ), &
                         IINFO )
!
!              Store as dense matrix for most routines.
!
            CALL SLASET( 'Full', LDA, N, 0.0E+0, 0.0E+0, A, LDA )
            DO IDIAG = -IHBW, IHBW
               IROW = IHBW - IDIAG + 1
               J1 = MAX( 1, IDIAG+1 )
               J2 = MIN( N, N+IDIAG )
               DO J = J1, J2
                  I = J - IDIAG
                  A( I, J ) = U( IROW, J )
               ENDDO
            ENDDO
         ELSE
            IINFO = 1
         END IF
!
         IF( IINFO /= 0 ) THEN
            WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, &
               IOLDSD
            INFO = ABS( IINFO )
            RETURN
         END IF
!
  110       CONTINUE
!
         ABSTOL = UNFL + UNFL
         IF( N <= 1 ) THEN
            IL = 1
            IU = N
         ELSE
            IL = 1 + INT( ( N-1 )*SLARND( 1, ISEED2 ) )
            IU = 1 + INT( ( N-1 )*SLARND( 1, ISEED2 ) )
            IF( IL > IU ) THEN
               ITEMP = IL
               IL = IU
               IU = ITEMP
            END IF
         END IF
!
!           3)      If matrix is tridiagonal, call SSTEV and SSTEVX.
!
         IF( JTYPE <= 7 ) THEN
            NTEST = 1
            DO I = 1, N
               D1( I ) = REAL( A( I, I ) )
               ENDDO
            DO I = 1, N - 1
               D2( I ) = REAL( A( I+1, I ) )
               ENDDO
            SRNAMT = 'SSTEV'
            CALL SSTEV( 'V', N, D1, D2, Z, LDU, WORK, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSTEV(V)', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 1 ) = ULPINV
                  RESULT( 2 ) = ULPINV
                  RESULT( 3 ) = ULPINV
                  GO TO 180
               END IF
            END IF
!
!              Do tests 1 and 2.
!
            DO I = 1, N
               D3( I ) = REAL( A( I, I ) )
            ENDDO
            DO I = 1, N - 1
               D4( I ) = REAL( A( I+1, I ) )
            ENDDO
            CALL SSTT21( N, 0, D3, D4, D1, D2, Z, LDU, WORK, &
                         RESULT( 1 ) )
!
            NTEST = 3
            DO I = 1, N - 1
               D4( I ) = REAL( A( I+1, I ) )
               ENDDO
            SRNAMT = 'SSTEV'
            CALL SSTEV( 'N', N, D3, D4, Z, LDU, WORK, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSTEV(N)', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 3 ) = ULPINV
                  GO TO 180
               END IF
            END IF
!
!              Do test 3.
!
            TEMP1 = 0.0E+0
            TEMP2 = 0.0E+0
            DO J = 1, N
               TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
               TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
               ENDDO
            RESULT( 3 ) = TEMP2 / MAX( UNFL, &
                          ULP*MAX( TEMP1, TEMP2 ) )
!
  180          CONTINUE
!
            NTEST = 4
            DO I = 1, N
               EVEIGS( I ) = D3( I )
               D1( I ) = REAL( A( I, I ) )
            ENDDO
            DO I = 1, N - 1
               D2( I ) = REAL( A( I+1, I ) )
            ENDDO
            SRNAMT = 'SSTEVX'
            CALL SSTEVX( 'V', 'A', N, D1, D2, VL, VU, IL, IU, ABSTOL, &
                         M, WA1, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), &
                         IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSTEVX(V,A)', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 4 ) = ULPINV
                  RESULT( 5 ) = ULPINV
                  RESULT( 6 ) = ULPINV
                  GO TO 250
               END IF
            END IF
            IF( N > 0 ) THEN
               TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
            ELSE
               TEMP3 = 0.0E+0
            END IF
!
!              Do tests 4 and 5.
!
            DO I = 1, N
               D3( I ) = REAL( A( I, I ) )
            ENDDO
            DO I = 1, N - 1
               D4( I ) = REAL( A( I+1, I ) )
            ENDDO
            CALL SSTT21( N, 0, D3, D4, WA1, D2, Z, LDU, WORK, &
                         RESULT( 4 ) )
!
            NTEST = 6
            DO I = 1, N - 1
               D4( I ) = REAL( A( I+1, I ) )
            ENDDO
            SRNAMT = 'SSTEVX'
            CALL SSTEVX( 'N', 'A', N, D3, D4, VL, VU, IL, IU, ABSTOL, &
                         M2, WA2, Z, LDU, WORK, IWORK, &
                         IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSTEVX(N,A)', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 6 ) = ULPINV
                  GO TO 250
               END IF
            END IF
!
!              Do test 6.
!
            TEMP1 = 0.0E+0
            TEMP2 = 0.0E+0
            DO J = 1, N
               TEMP1 = MAX( TEMP1, ABS( WA2( J ) ), &
                       ABS( EVEIGS( J ) ) )
               TEMP2 = MAX( TEMP2, ABS( WA2( J )-EVEIGS( J ) ) )
            ENDDO
            RESULT( 6 ) = TEMP2 / MAX( UNFL, &
                          ULP*MAX( TEMP1, TEMP2 ) )
!
  250          CONTINUE
!
            NTEST = 7
            DO I = 1, N
               D1( I ) = REAL( A( I, I ) )
            ENDDO
            DO I = 1, N - 1
               D2( I ) = REAL( A( I+1, I ) )
            ENDDO
            SRNAMT = 'SSTEVR'
            CALL SSTEVR( 'V', 'A', N, D1, D2, VL, VU, IL, IU, ABSTOL, &
                         M, WA1, Z, LDU, IWORK, WORK, LWORK, &
                         IWORK(2*N+1), LIWORK-2*N, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSTEVR(V,A)', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 7 ) = ULPINV
                  RESULT( 8 ) = ULPINV
                  GO TO 320
               END IF
            END IF
            IF( N > 0 ) THEN
               TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
            ELSE
               TEMP3 = 0.0E+0
            END IF
!
!              Do tests 7 and 8.
!
            DO I = 1, N
               D3( I ) = REAL( A( I, I ) )
            ENDDO
            DO I = 1, N - 1
               D4( I ) = REAL( A( I+1, I ) )
            ENDDO
            CALL SSTT21( N, 0, D3, D4, WA1, D2, Z, LDU, WORK, &
                         RESULT( 7 ) )
!
            NTEST = 9
            DO I = 1, N - 1
               D4( I ) = REAL( A( I+1, I ) )
            ENDDO
            SRNAMT = 'SSTEVR'
            CALL SSTEVR( 'N', 'A', N, D3, D4, VL, VU, IL, IU, ABSTOL, &
                         M2, WA2, Z, LDU, IWORK, WORK, LWORK, &
                         IWORK(2*N+1), LIWORK-2*N, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSTEVR(N,A)', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 9 ) = ULPINV
                  GO TO 320
               END IF
            END IF
!
!              Do test 9.
!
            TEMP1 = 0.0E+0
            TEMP2 = 0.0E+0
            DO J = 1, N
               TEMP1 = MAX( TEMP1, ABS( WA2( J ) ), &
                       ABS( EVEIGS( J ) ) )
               TEMP2 = MAX( TEMP2, ABS( WA2( J )-EVEIGS( J ) ) )
            ENDDO
            RESULT( 9 ) = TEMP2 / MAX( UNFL, &
                          ULP*MAX( TEMP1, TEMP2 ) )
!
  320          CONTINUE
!
!
            NTEST = 10
            DO I = 1, N
               D1( I ) = REAL( A( I, I ) )
            ENDDO
            DO I = 1, N - 1
               D2( I ) = REAL( A( I+1, I ) )
            ENDDO
            SRNAMT = 'SSTEVX'
            CALL SSTEVX( 'V', 'I', N, D1, D2, VL, VU, IL, IU, ABSTOL, &
                         M2, WA2, Z, LDU, WORK, IWORK, &
                         IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSTEVX(V,I)', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 10 ) = ULPINV
                  RESULT( 11 ) = ULPINV
                  RESULT( 12 ) = ULPINV
                  GO TO 380
               END IF
            END IF
!
!              Do tests 10 and 11.
!
            DO I = 1, N
               D3( I ) = REAL( A( I, I ) )
            ENDDO
            DO I = 1, N - 1
               D4( I ) = REAL( A( I+1, I ) )
            ENDDO
            CALL SSTT22( N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK, &
                         MAX( 1, M2 ), RESULT( 10 ) )
!
!
            NTEST = 12
            DO I = 1, N - 1
               D4( I ) = REAL( A( I+1, I ) )
               ENDDO
            SRNAMT = 'SSTEVX'
            CALL SSTEVX( 'N', 'I', N, D3, D4, VL, VU, IL, IU, ABSTOL, &
                         M3, WA3, Z, LDU, WORK, IWORK, &
                         IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSTEVX(N,I)', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 12 ) = ULPINV
                  GO TO 380
               END IF
            END IF
!
!              Do test 12.
!
            TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
            TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
            RESULT( 12 ) = ( TEMP1+TEMP2 ) / MAX( UNFL, ULP*TEMP3 )
!
  380          CONTINUE
!
            NTEST = 12
            IF( N > 0 ) THEN
               IF( IL /= 1 ) THEN
                  VL = WA1( IL ) - MAX( 0.5E+0* &
                       ( WA1( IL )-WA1( IL-1 ) ), 10.0E+0*ULP*TEMP3, &
                       10.0E+0*RTUNFL )
               ELSE
                  VL = WA1( 1 ) - MAX( 0.5E+0*( WA1( N )-WA1( 1 ) ), &
                       10.0E+0*ULP*TEMP3, 10.0E+0*RTUNFL )
               END IF
               IF( IU /= N ) THEN
                  VU = WA1( IU ) + MAX( 0.5E+0* &
                       ( WA1( IU+1 )-WA1( IU ) ), 10.0E+0*ULP*TEMP3, &
                       10.0E+0*RTUNFL )
               ELSE
                  VU = WA1( N ) + MAX( 0.5E+0*( WA1( N )-WA1( 1 ) ), &
                       10.0E+0*ULP*TEMP3, 10.0E+0*RTUNFL )
               END IF
            ELSE
               VL = 0.0E+0
               VU = 1.0E+0
            END IF
!
            DO I = 1, N
               D1( I ) = REAL( A( I, I ) )
               ENDDO
            DO I = 1, N - 1
               D2( I ) = REAL( A( I+1, I ) )
               ENDDO
            SRNAMT = 'SSTEVX'
            CALL SSTEVX( 'V', 'V', N, D1, D2, VL, VU, IL, IU, ABSTOL, &
                         M2, WA2, Z, LDU, WORK, IWORK, &
                         IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSTEVX(V,V)', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 13 ) = ULPINV
                  RESULT( 14 ) = ULPINV
                  RESULT( 15 ) = ULPINV
                  GO TO 440
               END IF
            END IF
!
            IF( M2 == 0 .AND. N > 0 ) THEN
               RESULT( 13 ) = ULPINV
               RESULT( 14 ) = ULPINV
               RESULT( 15 ) = ULPINV
               GO TO 440
            END IF
!
!              Do tests 13 and 14.
!
            DO I = 1, N
               D3( I ) = REAL( A( I, I ) )
               ENDDO
            DO I = 1, N - 1
               D4( I ) = REAL( A( I+1, I ) )
               ENDDO
            CALL SSTT22( N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK, &
                         MAX( 1, M2 ), RESULT( 13 ) )
!
            NTEST = 15
            DO I = 1, N - 1
               D4( I ) = REAL( A( I+1, I ) )
               ENDDO
            SRNAMT = 'SSTEVX'
            CALL SSTEVX( 'N', 'V', N, D3, D4, VL, VU, IL, IU, ABSTOL, &
                         M3, WA3, Z, LDU, WORK, IWORK, &
                         IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSTEVX(N,V)', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 15 ) = ULPINV
                  GO TO 440
               END IF
            END IF
!
!              Do test 15.
!
            TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
            TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
            RESULT( 15 ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )
!
  440          CONTINUE
!
            NTEST = 16
            DO I = 1, N
               D1( I ) = REAL( A( I, I ) )
            ENDDO
            DO I = 1, N - 1
               D2( I ) = REAL( A( I+1, I ) )
            ENDDO
            SRNAMT = 'SSTEVD'
            CALL SSTEVD( 'V', N, D1, D2, Z, LDU, WORK, LWEDC, IWORK, &
                         LIWEDC, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSTEVD(V)', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 16 ) = ULPINV
                  RESULT( 17 ) = ULPINV
                  RESULT( 18 ) = ULPINV
                  GO TO 510
               END IF
            END IF
!
!              Do tests 16 and 17.
!
            DO I = 1, N
               D3( I ) = REAL( A( I, I ) )
            ENDDO
            DO I = 1, N - 1
               D4( I ) = REAL( A( I+1, I ) )
            ENDDO
            CALL SSTT21( N, 0, D3, D4, D1, D2, Z, LDU, WORK, &
                         RESULT( 16 ) )
!
            NTEST = 18
            DO I = 1, N - 1
               D4( I ) = REAL( A( I+1, I ) )
               ENDDO
            SRNAMT = 'SSTEVD'
            CALL SSTEVD( 'N', N, D3, D4, Z, LDU, WORK, LWEDC, IWORK, &
                         LIWEDC, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSTEVD(N)', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 18 ) = ULPINV
                  GO TO 510
               END IF
            END IF
!
!              Do test 18.
!
            TEMP1 = 0.0E+0
            TEMP2 = 0.0E+0
            DO J = 1, N
               TEMP1 = MAX( TEMP1, ABS( EVEIGS( J ) ), &
                       ABS( D3( J ) ) )
               TEMP2 = MAX( TEMP2, ABS( EVEIGS( J )-D3( J ) ) )
               ENDDO
            RESULT( 18 ) = TEMP2 / MAX( UNFL, &
                           ULP*MAX( TEMP1, TEMP2 ) )
!
  510          CONTINUE
!
            NTEST = 19
            DO I = 1, N
               D1( I ) = REAL( A( I, I ) )
               ENDDO
            DO I = 1, N - 1
               D2( I ) = REAL( A( I+1, I ) )
               ENDDO
            SRNAMT = 'SSTEVR'
            CALL SSTEVR( 'V', 'I', N, D1, D2, VL, VU, IL, IU, ABSTOL, &
                         M2, WA2, Z, LDU, IWORK, WORK, LWORK, &
                         IWORK(2*N+1), LIWORK-2*N, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSTEVR(V,I)', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 19 ) = ULPINV
                  RESULT( 20 ) = ULPINV
                  RESULT( 21 ) = ULPINV
                  GO TO 570
               END IF
            END IF
!
!              DO tests 19 and 20.
!
            DO I = 1, N
               D3( I ) = REAL( A( I, I ) )
            ENDDO
            DO I = 1, N - 1
               D4( I ) = REAL( A( I+1, I ) )
            ENDDO
            CALL SSTT22( N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK, &
                         MAX( 1, M2 ), RESULT( 19 ) )
!
!
            NTEST = 21
            DO I = 1, N - 1
               D4( I ) = REAL( A( I+1, I ) )
               ENDDO
            SRNAMT = 'SSTEVR'
            CALL SSTEVR( 'N', 'I', N, D3, D4, VL, VU, IL, IU, ABSTOL, &
                         M3, WA3, Z, LDU, IWORK, WORK, LWORK, &
                         IWORK(2*N+1), LIWORK-2*N, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSTEVR(N,I)', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 21 ) = ULPINV
                  GO TO 570
               END IF
            END IF
!
!              Do test 21.
!
            TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
            TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
            RESULT( 21 ) = ( TEMP1+TEMP2 ) / MAX( UNFL, ULP*TEMP3 )
!
  570          CONTINUE
!
            NTEST = 21
            IF( N > 0 ) THEN
               IF( IL /= 1 ) THEN
                  VL = WA1( IL ) - MAX( 0.5E+0* &
                       ( WA1( IL )-WA1( IL-1 ) ), 10.0E+0*ULP*TEMP3, &
                       10.0E+0*RTUNFL )
               ELSE
                  VL = WA1( 1 ) - MAX( 0.5E+0*( WA1( N )-WA1( 1 ) ), &
                       10.0E+0*ULP*TEMP3, 10.0E+0*RTUNFL )
               END IF
               IF( IU /= N ) THEN
                  VU = WA1( IU ) + MAX( 0.5E+0* &
                       ( WA1( IU+1 )-WA1( IU ) ), 10.0E+0*ULP*TEMP3, &
                       10.0E+0*RTUNFL )
               ELSE
                  VU = WA1( N ) + MAX( 0.5E+0*( WA1( N )-WA1( 1 ) ), &
                       10.0E+0*ULP*TEMP3, 10.0E+0*RTUNFL )
               END IF
            ELSE
               VL = 0.0E+0
               VU = 1.0E+0
            END IF
!
            DO I = 1, N
               D1( I ) = REAL( A( I, I ) )
               ENDDO
            DO I = 1, N - 1
               D2( I ) = REAL( A( I+1, I ) )
               ENDDO
            SRNAMT = 'SSTEVR'
            CALL SSTEVR( 'V', 'V', N, D1, D2, VL, VU, IL, IU, ABSTOL, &
                         M2, WA2, Z, LDU, IWORK, WORK, LWORK, &
                         IWORK(2*N+1), LIWORK-2*N, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSTEVR(V,V)', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 22 ) = ULPINV
                  RESULT( 23 ) = ULPINV
                  RESULT( 24 ) = ULPINV
                  GO TO 630
               END IF
            END IF
!
            IF( M2 == 0 .AND. N > 0 ) THEN
               RESULT( 22 ) = ULPINV
               RESULT( 23 ) = ULPINV
               RESULT( 24 ) = ULPINV
               GO TO 630
            END IF
!
!              Do tests 22 and 23.
!
            DO I = 1, N
               D3( I ) = REAL( A( I, I ) )
            ENDDO
            DO I = 1, N - 1
               D4( I ) = REAL( A( I+1, I ) )
            ENDDO
            CALL SSTT22( N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK, &
                         MAX( 1, M2 ), RESULT( 22 ) )
!
            NTEST = 24
            DO I = 1, N - 1
               D4( I ) = REAL( A( I+1, I ) )
               ENDDO
            SRNAMT = 'SSTEVR'
            CALL SSTEVR( 'N', 'V', N, D3, D4, VL, VU, IL, IU, ABSTOL, &
                         M3, WA3, Z, LDU, IWORK, WORK, LWORK, &
                         IWORK(2*N+1), LIWORK-2*N, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSTEVR(N,V)', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 24 ) = ULPINV
                  GO TO 630
               END IF
            END IF
!
!              Do test 24.
!
            TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
            TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
            RESULT( 24 ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )
!
  630          CONTINUE
!
!
!
         ELSE
!
            RESULT(1:24) = 0.0E+0
            NTEST = 24
         END IF
!
!           Perform remaining tests storing upper or lower triangular
!           part of matrix.
!
         DO IUPLO = 0, 1
            IF( IUPLO == 0 ) THEN
               UPLO = 'L'
            ELSE
               UPLO = 'U'
            END IF
!
!              4)      Call SSYEV and SSYEVX.
!
            CALL SLACPY( ' ', N, N, A, LDA, V, LDU )
!
            NTEST = NTEST + 1
            SRNAMT = 'SSYEV'
            CALL SSYEV( 'V', UPLO, N, A, LDU, D1, WORK, LWORK, &
                        IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSYEV(V,' // UPLO // ')', &
                  IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  RESULT( NTEST+1 ) = ULPINV
                  RESULT( NTEST+2 ) = ULPINV
                  GO TO 660
               END IF
            END IF
!
!              Do tests 25 and 26 (or +54)
!
            CALL SSYT21( 1, UPLO, N, 0, V, LDU, D1, D2, A, LDU, Z, &
                         LDU, TAU, WORK, RESULT( NTEST ) )
!
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
!
            NTEST = NTEST + 2
            SRNAMT = 'SSYEV_2STAGE'
            CALL SSYEV_2STAGE( 'N', UPLO, N, A, LDU, D3, WORK, LWORK, &
                        IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 ) &
                  'SSYEV_2STAGE(N,' // UPLO // ')', &
                  IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  GO TO 660
               END IF
            END IF
!
!              Do test 27 (or +54)
!
            TEMP1 = MAX(MAXVAL(ABS(D1(1:N))), MAXVAL(ABS(D3(1:N))))
            TEMP2 = MAXVAL(ABS(D1(1:N)-D3(1:N)))
            RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
!
  660          CONTINUE
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
!
            NTEST = NTEST + 1
!
            IF( N > 0 ) THEN
               TEMP3 = MAX( ABS( D1( 1 ) ), ABS( D1( N ) ) )
               IF( IL /= 1 ) THEN
                  VL = D1( IL ) - MAX( 0.5E+0*( D1( IL )-D1( IL-1 ) ), &
                       10.0E+0*ULP*TEMP3, 10.0E+0*RTUNFL )
               ELSE IF( N > 0 ) THEN
                  VL = D1( 1 ) - MAX( 0.5E+0*( D1( N )-D1( 1 ) ), &
                       10.0E+0*ULP*TEMP3, 10.0E+0*RTUNFL )
               END IF
               IF( IU /= N ) THEN
                  VU = D1( IU ) + MAX( 0.5E+0*( D1( IU+1 )-D1( IU ) ), &
                       10.0E+0*ULP*TEMP3, 10.0E+0*RTUNFL )
               ELSE IF( N > 0 ) THEN
                  VU = D1( N ) + MAX( 0.5E+0*( D1( N )-D1( 1 ) ), &
                       10.0E+0*ULP*TEMP3, 10.0E+0*RTUNFL )
               END IF
            ELSE
               TEMP3 = 0.0E+0
               VL = 0.0E+0
               VU = 1.0E+0
            END IF
!
            SRNAMT = 'SSYEVX'
            CALL SSYEVX( 'V', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, &
                         ABSTOL, M, WA1, Z, LDU, WORK, LWORK, IWORK, &
                         IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSYEVX(V,A,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  RESULT( NTEST+1 ) = ULPINV
                  RESULT( NTEST+2 ) = ULPINV
                  GO TO 680
               END IF
            END IF
!
!              Do tests 28 and 29 (or +54)
!
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
!
            CALL SSYT21( 1, UPLO, N, 0, A, LDU, D1, D2, Z, LDU, V, &
                         LDU, TAU, WORK, RESULT( NTEST ) )
!
            NTEST = NTEST + 2
            SRNAMT = 'SSYEVX_2STAGE'
            CALL SSYEVX_2STAGE( 'N', 'A', UPLO, N, A, LDU, VL, VU, &
                         IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, &
                         LWORK, IWORK, IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 ) &
                  'SSYEVX_2STAGE(N,A,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  GO TO 680
               END IF
            END IF
!
!              Do test 30 (or +54)
!
            TEMP1 = MAX(MAXVAL(ABS(WA1(1:N))), MAXVAL(ABS(WA2(1:N))))
            TEMP2 = MAXVAL(ABS(WA1(1:N)-WA2(1:N)))
            RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
!
  680          CONTINUE
!
            NTEST = NTEST + 1
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
            SRNAMT = 'SSYEVX'
            CALL SSYEVX( 'V', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, &
                         ABSTOL, M2, WA2, Z, LDU, WORK, LWORK, IWORK, &
                         IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSYEVX(V,I,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST:NTEST+2 ) = ULPINV
                  GO TO 690
               END IF
            END IF
!
!              Do tests 31 and 32 (or +54)
!
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
!
            CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, &
                         V, LDU, TAU, WORK, RESULT( NTEST ) )
!
            NTEST = NTEST + 2
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
            SRNAMT = 'SSYEVX_2STAGE'
            CALL SSYEVX_2STAGE( 'N', 'I', UPLO, N, A, LDU, VL, VU, &
                         IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, &
                         LWORK, IWORK, IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 ) &
                  'SSYEVX_2STAGE(N,I,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  GO TO 690
               END IF
            END IF
!
!              Do test 33 (or +54)
!
            TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
            TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
            RESULT( NTEST ) = ( TEMP1+TEMP2 ) / &
                              MAX( UNFL, ULP*TEMP3 )
  690          CONTINUE
!
            NTEST = NTEST + 1
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
            SRNAMT = 'SSYEVX'
            CALL SSYEVX( 'V', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, &
                         ABSTOL, M2, WA2, Z, LDU, WORK, LWORK, IWORK, &
                         IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSYEVX(V,V,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST:NTEST+2 ) = ULPINV
                  GO TO 700
               END IF
            END IF
!
!              Do tests 34 and 35 (or +54)
!
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
!
            CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, &
                         V, LDU, TAU, WORK, RESULT( NTEST ) )
!
            NTEST = NTEST + 2
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
            SRNAMT = 'SSYEVX_2STAGE'
            CALL SSYEVX_2STAGE( 'N', 'V', UPLO, N, A, LDU, VL, VU, &
                         IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, &
                         LWORK, IWORK, IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 ) &
                  'SSYEVX_2STAGE(N,V,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  GO TO 700
               END IF
            END IF
!
            IF( M3 == 0 .AND. N > 0 ) THEN
               RESULT( NTEST ) = ULPINV
               GO TO 700
            END IF
!
!              Do test 36 (or +54)
!
            TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
            TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
            IF( N > 0 ) THEN
               TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
            ELSE
               TEMP3 = 0.0E+0
            END IF
            RESULT( NTEST ) = ( TEMP1+TEMP2 ) / &
                              MAX( UNFL, TEMP3*ULP )
!
  700          CONTINUE
!
!              5)      Call SSPEV and SSPEVX.
!
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
!
!              Load array WORK with the upper or lower triangular
!              part of the matrix in packed form.
!
            IF( IUPLO == 1 ) THEN
               INDX = 1
               DO J = 1, N
                  DO I = 1, J
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            ELSE
               INDX = 1
               DO J = 1, N
                  DO I = J, N
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            END IF
!
            NTEST = NTEST + 1
            SRNAMT = 'SSPEV'
            CALL SSPEV( 'V', UPLO, N, WORK, D1, Z, LDU, V, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSPEV(V,' // UPLO // ')', &
                  IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  RESULT( NTEST+1 ) = ULPINV
                  RESULT( NTEST+2 ) = ULPINV
                  GO TO 800
               END IF
            END IF
!
!              Do tests 37 and 38 (or +54)
!
            CALL SSYT21( 1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, &
                         LDU, TAU, WORK, RESULT( NTEST ) )
!
            IF( IUPLO == 1 ) THEN
               INDX = 1
               DO J = 1, N
                  DO I = 1, J
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            ELSE
               INDX = 1
               DO J = 1, N
                  DO I = J, N
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            END IF
!
            NTEST = NTEST + 2
            SRNAMT = 'SSPEV'
            CALL SSPEV( 'N', UPLO, N, WORK, D3, Z, LDU, V, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSPEV(N,' // UPLO // ')', &
                  IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  GO TO 800
               END IF
            END IF
!
!              Do test 39 (or +54)
!
            TEMP1 = MAX(MAXVAL(ABS(D1(1:N))), MAXVAL(ABS(D3(1:N))))
            TEMP2 = MAXVAL(ABS(D1(1:N)-D3(1:N)))
            RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
!
!              Load array WORK with the upper or lower triangular part
!              of the matrix in packed form.
!
  800          CONTINUE
            IF( IUPLO == 1 ) THEN
               INDX = 1
               DO J = 1, N
                  DO I = 1, J
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            ELSE
               INDX = 1
               DO J = 1, N
                  DO I = J, N
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            END IF
!
            NTEST = NTEST + 1
!
            IF( N > 0 ) THEN
               TEMP3 = MAX( ABS( D1( 1 ) ), ABS( D1( N ) ) )
               IF( IL /= 1 ) THEN
                  VL = D1( IL ) - MAX( 0.5E+0*( D1( IL )-D1( IL-1 ) ), &
                       10.0E+0*ULP*TEMP3, 10.0E+0*RTUNFL )
               ELSE IF( N > 0 ) THEN
                  VL = D1( 1 ) - MAX( 0.5E+0*( D1( N )-D1( 1 ) ), &
                       10.0E+0*ULP*TEMP3, 10.0E+0*RTUNFL )
               END IF
               IF( IU /= N ) THEN
                  VU = D1( IU ) + MAX( 0.5E+0*( D1( IU+1 )-D1( IU ) ), &
                       10.0E+0*ULP*TEMP3, 10.0E+0*RTUNFL )
               ELSE IF( N > 0 ) THEN
                  VU = D1( N ) + MAX( 0.5E+0*( D1( N )-D1( 1 ) ), &
                       10.0E+0*ULP*TEMP3, 10.0E+0*RTUNFL )
               END IF
            ELSE
               TEMP3 = 0.0E+0
               VL = 0.0E+0
               VU = 1.0E+0
            END IF
!
            SRNAMT = 'SSPEVX'
            CALL SSPEVX( 'V', 'A', UPLO, N, WORK, VL, VU, IL, IU, &
                         ABSTOL, M, WA1, Z, LDU, V, IWORK, &
                         IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSPEVX(V,A,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  RESULT( NTEST+1 ) = ULPINV
                  RESULT( NTEST+2 ) = ULPINV
                  GO TO 900
               END IF
            END IF
!
!              Do tests 40 and 41 (or +54)
!
            CALL SSYT21( 1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, &
                         LDU, TAU, WORK, RESULT( NTEST ) )
!
            NTEST = NTEST + 2
!
            IF( IUPLO == 1 ) THEN
               INDX = 1
               DO J = 1, N
                  DO I = 1, J
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            ELSE
               INDX = 1
               DO J = 1, N
                  DO I = J, N
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            END IF
!
            SRNAMT = 'SSPEVX'
            CALL SSPEVX( 'N', 'A', UPLO, N, WORK, VL, VU, IL, IU, &
                         ABSTOL, M2, WA2, Z, LDU, V, IWORK, &
                         IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSPEVX(N,A,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  GO TO 900
               END IF
            END IF
!
!              Do test 42 (or +54)
!
            TEMP1 = MAX(MAXVAL(ABS(WA1(1:N))), MAXVAL(ABS(WA2(1:N))))
            TEMP2 = MAXVAL(ABS(WA1(1:N)-WA2(1:N)))
            RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
!
  900          CONTINUE
            IF( IUPLO == 1 ) THEN
               INDX = 1
               DO J = 1, N
                  DO I = 1, J
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            ELSE
               INDX = 1
               DO J = 1, N
                  DO I = J, N
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            END IF
!
            NTEST = NTEST + 1
!
            SRNAMT = 'SSPEVX'
            CALL SSPEVX( 'V', 'I', UPLO, N, WORK, VL, VU, IL, IU, &
                         ABSTOL, M2, WA2, Z, LDU, V, IWORK, &
                         IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSPEVX(V,I,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  RESULT( NTEST+1 ) = ULPINV
                  RESULT( NTEST+2 ) = ULPINV
                  GO TO 990
               END IF
            END IF
!
!              Do tests 43 and 44 (or +54)
!
            CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, &
                         V, LDU, TAU, WORK, RESULT( NTEST ) )
!
            NTEST = NTEST + 2
!
            IF( IUPLO == 1 ) THEN
               INDX = 1
               DO J = 1, N
                  DO I = 1, J
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            ELSE
               INDX = 1
               DO J = 1, N
                  DO I = J, N
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            END IF
!
            SRNAMT = 'SSPEVX'
            CALL SSPEVX( 'N', 'I', UPLO, N, WORK, VL, VU, IL, IU, &
                         ABSTOL, M3, WA3, Z, LDU, V, IWORK, &
                         IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSPEVX(N,I,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  GO TO 990
               END IF
            END IF
!
            IF( M3 == 0 .AND. N > 0 ) THEN
               RESULT( NTEST ) = ULPINV
               GO TO 990
            END IF
!
!              Do test 45 (or +54)
!
            TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
            TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
            IF( N > 0 ) THEN
               TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
            ELSE
               TEMP3 = 0.0E+0
            END IF
            RESULT( NTEST ) = ( TEMP1+TEMP2 ) / &
                              MAX( UNFL, TEMP3*ULP )
!
  990          CONTINUE
            IF( IUPLO == 1 ) THEN
               INDX = 1
               DO J = 1, N
                  DO I = 1, J
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            ELSE
               INDX = 1
               DO J = 1, N
                  DO I = J, N
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            END IF
!
            NTEST = NTEST + 1
!
            SRNAMT = 'SSPEVX'
            CALL SSPEVX( 'V', 'V', UPLO, N, WORK, VL, VU, IL, IU, &
                         ABSTOL, M2, WA2, Z, LDU, V, IWORK, &
                         IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSPEVX(V,V,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  RESULT( NTEST+1 ) = ULPINV
                  RESULT( NTEST+2 ) = ULPINV
                  GO TO 1080
               END IF
            END IF
!
!              Do tests 46 and 47 (or +54)
!
            CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, &
                         V, LDU, TAU, WORK, RESULT( NTEST ) )
!
            NTEST = NTEST + 2
!
            IF( IUPLO == 1 ) THEN
               INDX = 1
               DO J = 1, N
                  DO I = 1, J
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            ELSE
               INDX = 1
               DO J = 1, N
                  DO I = J, N
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            END IF
!
            SRNAMT = 'SSPEVX'
            CALL SSPEVX( 'N', 'V', UPLO, N, WORK, VL, VU, IL, IU, &
                         ABSTOL, M3, WA3, Z, LDU, V, IWORK, &
                         IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSPEVX(N,V,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  GO TO 1080
               END IF
            END IF
!
            IF( M3 == 0 .AND. N > 0 ) THEN
               RESULT( NTEST ) = ULPINV
               GO TO 1080
            END IF
!
!              Do test 48 (or +54)
!
            TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
            TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
            IF( N > 0 ) THEN
               TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
            ELSE
               TEMP3 = 0.0E+0
            END IF
            RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )
!
 1080          CONTINUE
!
!              6)      Call SSBEV and SSBEVX.
!
            IF( JTYPE <= 7 ) THEN
               KD = 1
            ELSE IF( JTYPE >= 8 .AND. JTYPE <= 15 ) THEN
               KD = MAX( N-1, 0 )
            ELSE
               KD = IHBW
            END IF
!
!              Load array V with the upper or lower triangular part
!              of the matrix in band form.
!
            IF( IUPLO == 1 ) THEN
               DO J = 1, N
                  DO I = MAX( 1, J-KD ), J
                     V( KD+1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            ELSE
               DO J = 1, N
                  DO I = J, MIN( N, J+KD )
                     V( 1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            END IF
!
            NTEST = NTEST + 1
            SRNAMT = 'SSBEV'
            CALL SSBEV( 'V', UPLO, N, KD, V, LDU, D1, Z, LDU, WORK, &
                        IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSBEV(V,' // UPLO // ')', &
                  IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  RESULT( NTEST+1 ) = ULPINV
                  RESULT( NTEST+2 ) = ULPINV
                  GO TO 1180
               END IF
            END IF
!
!              Do tests 49 and 50 (or ... )
!
            CALL SSYT21( 1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, &
                         LDU, TAU, WORK, RESULT( NTEST ) )
!
            IF( IUPLO == 1 ) THEN
               DO J = 1, N
                  DO I = MAX( 1, J-KD ), J
                     V( KD+1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            ELSE
               DO J = 1, N
                  DO I = J, MIN( N, J+KD )
                     V( 1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            END IF
!
            NTEST = NTEST + 2
            SRNAMT = 'SSBEV_2STAGE'
            CALL SSBEV_2STAGE( 'N', UPLO, N, KD, V, LDU, D3, Z, LDU, &
                        WORK, LWORK, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 ) 'SSBEV_2STAGE(N,' // UPLO // ')', &
                  IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  GO TO 1180
               END IF
            END IF
!
!              Do test 51 (or +54)
!
            TEMP1 = MAX(MAXVAL(ABS(D1(1:N))), MAXVAL(ABS(D3(1:N))))
            TEMP2 = MAXVAL(ABS(D1(1:N)-D3(1:N)))
            RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
!
!              Load array V with the upper or lower triangular part
!              of the matrix in band form.
!
 1180          CONTINUE
            IF( IUPLO == 1 ) THEN
               DO J = 1, N
                  DO I = MAX( 1, J-KD ), J
                     V( KD+1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            ELSE
               DO J = 1, N
                  DO I = J, MIN( N, J+KD )
                     V( 1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            END IF
!
            NTEST = NTEST + 1
            SRNAMT = 'SSBEVX'
            CALL SSBEVX( 'V', 'A', UPLO, N, KD, V, LDU, U, LDU, VL, &
                         VU, IL, IU, ABSTOL, M, WA2, Z, LDU, WORK, &
                         IWORK, IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSBEVX(V,A,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  RESULT( NTEST+1 ) = ULPINV
                  RESULT( NTEST+2 ) = ULPINV
                  GO TO 1280
               END IF
            END IF
!
!              Do tests 52 and 53 (or +54)
!
            CALL SSYT21( 1, UPLO, N, 0, A, LDU, WA2, D2, Z, LDU, V, &
                         LDU, TAU, WORK, RESULT( NTEST ) )
!
            NTEST = NTEST + 2
!
            IF( IUPLO == 1 ) THEN
               DO J = 1, N
                  DO I = MAX( 1, J-KD ), J
                     V( KD+1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            ELSE
               DO J = 1, N
                  DO I = J, MIN( N, J+KD )
                     V( 1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            END IF
!
            SRNAMT = 'SSBEVX_2STAGE'
            CALL SSBEVX_2STAGE( 'N', 'A', UPLO, N, KD, V, LDU, &
                         U, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, &
                         Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ), &
                         IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 ) &
                  'SSBEVX_2STAGE(N,A,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  GO TO 1280
               END IF
            END IF
!
!              Do test 54 (or +54)
!
            TEMP1 = MAX(MAXVAL(ABS(WA2(1:N))), MAXVAL(ABS(WA3(1:N))))
            TEMP2 = MAXVAL(ABS(WA2(1:N)-WA3(1:N)))
            RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
!
 1280          CONTINUE
            NTEST = NTEST + 1
            IF( IUPLO == 1 ) THEN
               DO J = 1, N
                  DO I = MAX( 1, J-KD ), J
                     V( KD+1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            ELSE
               DO J = 1, N
                  DO I = J, MIN( N, J+KD )
                     V( 1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            END IF
!
            SRNAMT = 'SSBEVX'
            CALL SSBEVX( 'V', 'I', UPLO, N, KD, V, LDU, U, LDU, VL, &
                         VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, &
                         IWORK, IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSBEVX(V,I,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  RESULT( NTEST+1 ) = ULPINV
                  RESULT( NTEST+2 ) = ULPINV
                  GO TO 1370
               END IF
            END IF
!
!              Do tests 55 and 56 (or +54)
!
            CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, &
                         V, LDU, TAU, WORK, RESULT( NTEST ) )
!
            NTEST = NTEST + 2
!
            IF( IUPLO == 1 ) THEN
               DO J = 1, N
                  DO I = MAX( 1, J-KD ), J
                     V( KD+1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            ELSE
               DO J = 1, N
                  DO I = J, MIN( N, J+KD )
                     V( 1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            END IF
!
            SRNAMT = 'SSBEVX_2STAGE'
            CALL SSBEVX_2STAGE( 'N', 'I', UPLO, N, KD, V, LDU, &
                         U, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, &
                         Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ), &
                         IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 ) &
                  'SSBEVX_2STAGE(N,I,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  GO TO 1370
               END IF
            END IF
!
!              Do test 57 (or +54)
!
            TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
            TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
            IF( N > 0 ) THEN
               TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
            ELSE
               TEMP3 = 0.0E+0
            END IF
            RESULT( NTEST ) = ( TEMP1+TEMP2 ) / &
                              MAX( UNFL, TEMP3*ULP )
!
 1370          CONTINUE
            NTEST = NTEST + 1
            IF( IUPLO == 1 ) THEN
               DO J = 1, N
                  DO I = MAX( 1, J-KD ), J
                     V( KD+1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            ELSE
               DO J = 1, N
                  DO I = J, MIN( N, J+KD )
                     V( 1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            END IF
!
            SRNAMT = 'SSBEVX'
            CALL SSBEVX( 'V', 'V', UPLO, N, KD, V, LDU, U, LDU, VL, &
                         VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, &
                         IWORK, IWORK( 5*N+1 ), IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSBEVX(V,V,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  RESULT( NTEST+1 ) = ULPINV
                  RESULT( NTEST+2 ) = ULPINV
                  GO TO 1460
               END IF
            END IF
!
!              Do tests 58 and 59 (or +54)
!
            CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, &
                         V, LDU, TAU, WORK, RESULT( NTEST ) )
!
            NTEST = NTEST + 2
!
            IF( IUPLO == 1 ) THEN
               DO J = 1, N
                  DO I = MAX( 1, J-KD ), J
                     V( KD+1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            ELSE
               DO J = 1, N
                  DO I = J, MIN( N, J+KD )
                     V( 1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            END IF
!
            SRNAMT = 'SSBEVX_2STAGE'
            CALL SSBEVX_2STAGE( 'N', 'V', UPLO, N, KD, V, LDU, &
                         U, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, &
                         Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ), &
                         IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 ) &
                  'SSBEVX_2STAGE(N,V,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  GO TO 1460
               END IF
            END IF
!
            IF( M3 == 0 .AND. N > 0 ) THEN
               RESULT( NTEST ) = ULPINV
               GO TO 1460
            END IF
!
!              Do test 60 (or +54)
!
            TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
            TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
            IF( N > 0 ) THEN
               TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
            ELSE
               TEMP3 = 0.0E+0
            END IF
            RESULT( NTEST ) = ( TEMP1+TEMP2 ) / &
                              MAX( UNFL, TEMP3*ULP )
!
 1460          CONTINUE
!
!              7)      Call SSYEVD
!
            CALL SLACPY( ' ', N, N, A, LDA, V, LDU )
!
            NTEST = NTEST + 1
            SRNAMT = 'SSYEVD'
            CALL SSYEVD( 'V', UPLO, N, A, LDU, D1, WORK, LWEDC, &
                         IWORK, LIWEDC, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSYEVD(V,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  RESULT( NTEST+1 ) = ULPINV
                  RESULT( NTEST+2 ) = ULPINV
                  GO TO 1480
               END IF
            END IF
!
!              Do tests 61 and 62 (or +54)
!
            CALL SSYT21( 1, UPLO, N, 0, V, LDU, D1, D2, A, LDU, Z, &
                         LDU, TAU, WORK, RESULT( NTEST ) )
!
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
!
            NTEST = NTEST + 2
            SRNAMT = 'SSYEVD_2STAGE'
            CALL SSYEVD_2STAGE( 'N', UPLO, N, A, LDU, D3, WORK, &
                                 LWORK, IWORK, LIWEDC, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 ) &
                  'SSYEVD_2STAGE(N,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  GO TO 1480
               END IF
            END IF
!
!              Do test 63 (or +54)
!
            TEMP1 = MAX(MAXVAL(ABS(D1(1:N))), MAXVAL(ABS(D3(1:N))))
            TEMP2 = MAXVAL(ABS(D1(1:N)-D3(1:N)))
            RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
!
 1480          CONTINUE
!
!              8)      Call SSPEVD.
!
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
!
!              Load array WORK with the upper or lower triangular
!              part of the matrix in packed form.
!
            IF( IUPLO == 1 ) THEN
               INDX = 1
               DO J = 1, N
                  DO I = 1, J
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            ELSE
               INDX = 1
               DO J = 1, N
                  DO I = J, N
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            END IF
!
            NTEST = NTEST + 1
            SRNAMT = 'SSPEVD'
            CALL SSPEVD( 'V', UPLO, N, WORK, D1, Z, LDU, &
                         WORK( INDX ), LWEDC-INDX+1, IWORK, LIWEDC, &
                         IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSPEVD(V,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  RESULT( NTEST+1 ) = ULPINV
                  RESULT( NTEST+2 ) = ULPINV
                  GO TO 1580
               END IF
            END IF
!
!              Do tests 64 and 65 (or +54)
!
            CALL SSYT21( 1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, &
                         LDU, TAU, WORK, RESULT( NTEST ) )
!
            IF( IUPLO == 1 ) THEN
               INDX = 1
               DO J = 1, N
                  DO I = 1, J
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            ELSE
               INDX = 1
               DO J = 1, N
                  DO I = J, N
                     WORK( INDX ) = A( I, J )
                     INDX = INDX + 1
                  ENDDO
               ENDDO
            END IF
!
            NTEST = NTEST + 2
            SRNAMT = 'SSPEVD'
            CALL SSPEVD( 'N', UPLO, N, WORK, D3, Z, LDU, &
                         WORK( INDX ), LWEDC-INDX+1, IWORK, LIWEDC, &
                         IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSPEVD(N,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  GO TO 1580
               END IF
            END IF
!
!              Do test 66 (or +54)
!
            TEMP1 = MAX(MAXVAL(ABS(D1(1:N))), MAXVAL(ABS(D3(1:N))))
            TEMP2 = MAXVAL(ABS(D1(1:N)-D3(1:N)))
            RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
 1580          CONTINUE
!
!              9)      Call SSBEVD.
!
            IF( JTYPE <= 7 ) THEN
               KD = 1
            ELSE IF( JTYPE >= 8 .AND. JTYPE <= 15 ) THEN
               KD = MAX( N-1, 0 )
            ELSE
               KD = IHBW
            END IF
!
!              Load array V with the upper or lower triangular part
!              of the matrix in band form.
!
            IF( IUPLO == 1 ) THEN
               DO J = 1, N
                  DO I = MAX( 1, J-KD ), J
                     V( KD+1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            ELSE
               DO J = 1, N
                  DO I = J, MIN( N, J+KD )
                     V( 1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            END IF
!
            NTEST = NTEST + 1
            SRNAMT = 'SSBEVD'
            CALL SSBEVD( 'V', UPLO, N, KD, V, LDU, D1, Z, LDU, WORK, &
                         LWEDC, IWORK, LIWEDC, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSBEVD(V,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  RESULT( NTEST+1 ) = ULPINV
                  RESULT( NTEST+2 ) = ULPINV
                  GO TO 1680
               END IF
            END IF
!
!              Do tests 67 and 68 (or +54)
!
            CALL SSYT21( 1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, &
                         LDU, TAU, WORK, RESULT( NTEST ) )
!
            IF( IUPLO == 1 ) THEN
               DO J = 1, N
                  DO I = MAX( 1, J-KD ), J
                     V( KD+1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            ELSE
               DO J = 1, N
                  DO I = J, MIN( N, J+KD )
                     V( 1+I-J, J ) = A( I, J )
                  ENDDO
               ENDDO
            END IF
!
            NTEST = NTEST + 2
            SRNAMT = 'SSBEVD_2STAGE'
            CALL SSBEVD_2STAGE( 'N', UPLO, N, KD, V, LDU, D3, Z, LDU, &
                                WORK, LWORK, IWORK, LIWEDC, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 ) &
                  'SSBEVD_2STAGE(N,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  GO TO 1680
               END IF
            END IF
!
!              Do test 69 (or +54)
!
            TEMP1 = MAX(MAXVAL(ABS(D1(1:N))), MAXVAL(ABS(D3(1:N))))
            TEMP2 = MAXVAL(ABS(D1(1:N)-D3(1:N)))
            RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
!
 1680          CONTINUE
!
!
            CALL SLACPY( ' ', N, N, A, LDA, V, LDU )
            NTEST = NTEST + 1
            SRNAMT = 'SSYEVR'
            CALL SSYEVR( 'V', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, &
                         ABSTOL, M, WA1, Z, LDU, IWORK, WORK, LWORK, &
                         IWORK(2*N+1), LIWORK-2*N, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSYEVR(V,A,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  RESULT( NTEST+1 ) = ULPINV
                  RESULT( NTEST+2 ) = ULPINV
                  GO TO 1700
               END IF
            END IF
!
!              Do tests 70 and 71 (or ... )
!
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
!
            CALL SSYT21( 1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, &
                         LDU, TAU, WORK, RESULT( NTEST ) )
!
            NTEST = NTEST + 2
            SRNAMT = 'SSYEVR_2STAGE'
            CALL SSYEVR_2STAGE( 'N', 'A', UPLO, N, A, LDU, VL, VU, &
                         IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, &
                         WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, &
                         IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 ) &
                  'SSYEVR_2STAGE(N,A,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  GO TO 1700
               END IF
            END IF
!
!              Do test 72 (or ... )
!
            TEMP1 = MAX(MAXVAL(ABS(WA1(1:N))), MAXVAL(ABS(WA2(1:N))))
            TEMP2 = MAXVAL(ABS(WA1(1:N)-WA2(1:N)))
            RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
!
 1700          CONTINUE
!
            NTEST = NTEST + 1
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
            SRNAMT = 'SSYEVR'
            CALL SSYEVR( 'V', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, &
                         ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, &
                         IWORK(2*N+1), LIWORK-2*N, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSYEVR(V,I,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  RESULT( NTEST+1 ) = ULPINV
                  RESULT( NTEST+2 ) = ULPINV
                  GO TO 1710
               END IF
            END IF
!
!              Do tests 73 and 74 (or +54)
!
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
!
            CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, &
                         V, LDU, TAU, WORK, RESULT( NTEST ) )
!
            NTEST = NTEST + 2
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
            SRNAMT = 'SSYEVR_2STAGE'
            CALL SSYEVR_2STAGE( 'N', 'I', UPLO, N, A, LDU, VL, VU, &
                         IL, IU, ABSTOL, M3, WA3, Z, LDU, IWORK, &
                         WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, &
                         IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 ) &
                  'SSYEVR_2STAGE(N,I,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  GO TO 1710
               END IF
            END IF
!
!              Do test 75 (or +54)
!
            TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
            TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
            RESULT( NTEST ) = ( TEMP1+TEMP2 ) / &
                              MAX( UNFL, ULP*TEMP3 )
 1710          CONTINUE
!
            NTEST = NTEST + 1
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
            SRNAMT = 'SSYEVR'
            CALL SSYEVR( 'V', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, &
                         ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, &
                         IWORK(2*N+1), LIWORK-2*N, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SSYEVR(V,V,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST:NTEST+2  ) = ULPINV
                  GO TO 700
               END IF
            END IF
!
!              Do tests 76 and 77 (or +54)
!
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
!
            CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, &
                         V, LDU, TAU, WORK, RESULT( NTEST ) )
!
            NTEST = NTEST + 2
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
            SRNAMT = 'SSYEVR_2STAGE'
            CALL SSYEVR_2STAGE( 'N', 'V', UPLO, N, A, LDU, VL, VU, &
                         IL, IU, ABSTOL, M3, WA3, Z, LDU, IWORK, &
                         WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, &
                         IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 ) &
                  'SSYEVR_2STAGE(N,V,' // UPLO // &
                  ')', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( NTEST ) = ULPINV
                  GO TO 700
               END IF
            END IF
!
            IF( M3 == 0 .AND. N > 0 ) THEN
               RESULT( NTEST ) = ULPINV
               GO TO 700
            END IF
!
!              Do test 78 (or +54)
!
            TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
            TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
            IF( N > 0 ) THEN
               TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
            ELSE
               TEMP3 = 0.0E+0
            END IF
            RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )
!
            CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
!
            ENDDO
!
!           End of Loop -- Check for RESULT(j) > THRESH
!
         NTESTT = NTESTT + NTEST
         CALL SLAFTS( 'SST', N, N, JTYPE, NTEST, RESULT, IOLDSD, &
                      THRESH, NOUNIT, NERRS )
!
 1730    CONTINUE
         ENDDO
      ENDDO
!
!     Summary
!
   CALL ALASVM( 'SST', NOUNIT, NERRS, NTESTT, 0 )
!
 9999 FORMAT( ' SDRVST2STG: ', A, ' returned INFO=', I6, '.', / 9X, &
       'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
!
   RETURN
!
!     End of SDRVST2STG
!
END
