!> \brief \b DCHKHS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCHKHS( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                          NOUNIT, A, LDA, H, T1, T2, U, LDU, Z, UZ, WR1,
!                          WI1, WR2, WI2, WR3, WI3, EVECTL, EVECTR, EVECTY,
!                          EVECTX, UU, TAU, WORK, NWORK, IWORK, SELECT,
!                          RESULT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDU, NOUNIT, NSIZES, NTYPES, NWORK
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * ), SELECT( * )
!       INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
!       DOUBLE PRECISION   A( LDA, * ), EVECTL( LDU, * ),
!      $                   EVECTR( LDU, * ), EVECTX( LDU, * ),
!      $                   EVECTY( LDU, * ), H( LDA, * ), RESULT( 16 ),
!      $                   T1( LDA, * ), T2( LDA, * ), TAU( * ),
!      $                   U( LDU, * ), UU( LDU, * ), UZ( LDU, * ),
!      $                   WI1( * ), WI2( * ), WI3( * ), WORK( * ),
!      $                   WR1( * ), WR2( * ), WR3( * ), Z( LDU, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DCHKHS  checks the nonsymmetric eigenvalue problem routines.
!>
!>            DGEHRD factors A as  U H U' , where ' means transpose,
!>            H is hessenberg, and U is an orthogonal matrix.
!>
!>            DORGHR generates the orthogonal matrix U.
!>
!>            DORMHR multiplies a matrix by the orthogonal matrix U.
!>
!>            DHSEQR factors H as  Z T Z' , where Z is orthogonal and
!>            T is "quasi-triangular", and the eigenvalue vector W.
!>
!>            DTREVC computes the left and right eigenvector matrices
!>            L and R for T. L is lower quasi-triangular, and R is
!>            upper quasi-triangular.
!>
!>            DHSEIN computes the left and right eigenvector matrices
!>            Y and X for H, using inverse iteration.
!>
!>            DTREVC3 computes left and right eigenvector matrices
!>            from a Schur matrix T and backtransforms them with Z
!>            to eigenvector matrices L and R for A. L and R are
!>            GE matrices.
!>
!>    When DCHKHS is called, a number of matrix "sizes" ("n's") and a
!>    number of matrix "types" are specified.  For each size ("n")
!>    and each type of matrix, one matrix will be generated and used
!>    to test the nonsymmetric eigenroutines.  For each matrix, 16
!>    tests will be performed:
!>
!>    (1)     | A - U H U**T | / ( |A| n ulp )
!>
!>    (2)     | I - UU**T | / ( n ulp )
!>
!>    (3)     | H - Z T Z**T | / ( |H| n ulp )
!>
!>    (4)     | I - ZZ**T | / ( n ulp )
!>
!>    (5)     | A - UZ H (UZ)**T | / ( |A| n ulp )
!>
!>    (6)     | I - UZ (UZ)**T | / ( n ulp )
!>
!>    (7)     | T(Z computed) - T(Z not computed) | / ( |T| ulp )
!>
!>    (8)     | W(Z computed) - W(Z not computed) | / ( |W| ulp )
!>
!>    (9)     | TR - RW | / ( |T| |R| ulp )
!>
!>    (10)    | L**H T - W**H L | / ( |T| |L| ulp )
!>
!>    (11)    | HX - XW | / ( |H| |X| ulp )
!>
!>    (12)    | Y**H H - W**H Y | / ( |H| |Y| ulp )
!>
!>    (13)    | AX - XW | / ( |A| |X| ulp )
!>
!>    (14)    | Y**H A - W**H Y | / ( |A| |Y| ulp )
!>
!>    (15)    | AR - RW | / ( |A| |R| ulp )
!>
!>    (16)    | LA - WL | / ( |A| |L| ulp )
!>
!>    The "sizes" are specified by an array NN(1:NSIZES); the value of
!>    each element NN(j) specifies one size.
!>    The "types" are specified by a logical array DOTYPE( 1:NTYPES );
!>    if DOTYPE(j) is .TRUE., then matrix type "j" will be generated.
!>    Currently, the list of possible types is:
!>
!>    (1)  The zero matrix.
!>    (2)  The identity matrix.
!>    (3)  A (transposed) Jordan block, with 1's on the diagonal.
!>
!>    (4)  A diagonal matrix with evenly spaced entries
!>         1, ..., ULP  and random signs.
!>         (ULP = (first number larger than 1) - 1 )
!>    (5)  A diagonal matrix with geometrically spaced entries
!>         1, ..., ULP  and random signs.
!>    (6)  A diagonal matrix with "clustered" entries 1, ULP, ..., ULP
!>         and random signs.
!>
!>    (7)  Same as (4), but multiplied by SQRT( overflow threshold )
!>    (8)  Same as (4), but multiplied by SQRT( underflow threshold )
!>
!>    (9)  A matrix of the form  U' T U, where U is orthogonal and
!>         T has evenly spaced entries 1, ..., ULP with random signs
!>         on the diagonal and random O(1) entries in the upper
!>         triangle.
!>
!>    (10) A matrix of the form  U' T U, where U is orthogonal and
!>         T has geometrically spaced entries 1, ..., ULP with random
!>         signs on the diagonal and random O(1) entries in the upper
!>         triangle.
!>
!>    (11) A matrix of the form  U' T U, where U is orthogonal and
!>         T has "clustered" entries 1, ULP,..., ULP with random
!>         signs on the diagonal and random O(1) entries in the upper
!>         triangle.
!>
!>    (12) A matrix of the form  U' T U, where U is orthogonal and
!>         T has real or complex conjugate paired eigenvalues randomly
!>         chosen from ( ULP, 1 ) and random O(1) entries in the upper
!>         triangle.
!>
!>    (13) A matrix of the form  X' T X, where X has condition
!>         SQRT( ULP ) and T has evenly spaced entries 1, ..., ULP
!>         with random signs on the diagonal and random O(1) entries
!>         in the upper triangle.
!>
!>    (14) A matrix of the form  X' T X, where X has condition
!>         SQRT( ULP ) and T has geometrically spaced entries
!>         1, ..., ULP with random signs on the diagonal and random
!>         O(1) entries in the upper triangle.
!>
!>    (15) A matrix of the form  X' T X, where X has condition
!>         SQRT( ULP ) and T has "clustered" entries 1, ULP,..., ULP
!>         with random signs on the diagonal and random O(1) entries
!>         in the upper triangle.
!>
!>    (16) A matrix of the form  X' T X, where X has condition
!>         SQRT( ULP ) and T has real or complex conjugate paired
!>         eigenvalues randomly chosen from ( ULP, 1 ) and random
!>         O(1) entries in the upper triangle.
!>
!>    (17) Same as (16), but multiplied by SQRT( overflow threshold )
!>    (18) Same as (16), but multiplied by SQRT( underflow threshold )
!>
!>    (19) Nonsymmetric matrix with random entries chosen from (-1,1).
!>    (20) Same as (19), but multiplied by SQRT( overflow threshold )
!>    (21) Same as (19), but multiplied by SQRT( underflow threshold )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \verbatim
!>  NSIZES - INTEGER
!>           The number of sizes of matrices to use.  If it is zero,
!>           DCHKHS does nothing.  It must be at least zero.
!>           Not modified.
!>
!>  NN     - INTEGER array, dimension (NSIZES)
!>           An array containing the sizes to be used for the matrices.
!>           Zero values will be skipped.  The values must be at least
!>           zero.
!>           Not modified.
!>
!>  NTYPES - INTEGER
!>           The number of elements in DOTYPE.   If it is zero, DCHKHS
!>           does nothing.  It must be at least zero.  If it is MAXTYP+1
!>           and NSIZES is 1, then an additional type, MAXTYP+1 is
!>           defined, which is to use whatever matrix is in A.  This
!>           is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
!>           DOTYPE(MAXTYP+1) is .TRUE. .
!>           Not modified.
!>
!>  DOTYPE - LOGICAL array, dimension (NTYPES)
!>           If DOTYPE(j) is .TRUE., then for each size in NN a
!>           matrix of that size and of type j will be generated.
!>           If NTYPES is smaller than the maximum number of types
!>           defined (PARAMETER MAXTYP), then types NTYPES+1 through
!>           MAXTYP will not be generated.  If NTYPES is larger
!>           than MAXTYP, DOTYPE(MAXTYP+1) through DOTYPE(NTYPES)
!>           will be ignored.
!>           Not modified.
!>
!>  ISEED  - INTEGER array, dimension (4)
!>           On entry ISEED specifies the seed of the random number
!>           generator. The array elements should be between 0 and 4095;
!>           if not they will be reduced mod 4096.  Also, ISEED(4) must
!>           be odd.  The random number generator uses a linear
!>           congruential sequence limited to small integers, and so
!>           should produce machine independent random numbers. The
!>           values of ISEED are changed on exit, and can be used in the
!>           next call to DCHKHS to continue the same random number
!>           sequence.
!>           Modified.
!>
!>  THRESH - DOUBLE PRECISION
!>           A test will count as "failed" if the "error", computed as
!>           described above, exceeds THRESH.  Note that the error
!>           is scaled to be O(1), so THRESH should be a reasonably
!>           small multiple of 1, e.g., 10 or 100.  In particular,
!>           it should not depend on the precision (single vs. double)
!>           or the size of the matrix.  It must be at least zero.
!>           Not modified.
!>
!>  NOUNIT - INTEGER
!>           The FORTRAN unit number for printing out error messages
!>           (e.g., if a routine returns IINFO not equal to 0.)
!>           Not modified.
!>
!>  A      - DOUBLE PRECISION array, dimension (LDA,max(NN))
!>           Used to hold the matrix whose eigenvalues are to be
!>           computed.  On exit, A contains the last matrix actually
!>           used.
!>           Modified.
!>
!>  LDA    - INTEGER
!>           The leading dimension of A, H, T1 and T2.  It must be at
!>           least 1 and at least max( NN ).
!>           Not modified.
!>
!>  H      - DOUBLE PRECISION array, dimension (LDA,max(NN))
!>           The upper hessenberg matrix computed by DGEHRD.  On exit,
!>           H contains the Hessenberg form of the matrix in A.
!>           Modified.
!>
!>  T1     - DOUBLE PRECISION array, dimension (LDA,max(NN))
!>           The Schur (="quasi-triangular") matrix computed by DHSEQR
!>           if Z is computed.  On exit, T1 contains the Schur form of
!>           the matrix in A.
!>           Modified.
!>
!>  T2     - DOUBLE PRECISION array, dimension (LDA,max(NN))
!>           The Schur matrix computed by DHSEQR when Z is not computed.
!>           This should be identical to T1.
!>           Modified.
!>
!>  LDU    - INTEGER
!>           The leading dimension of U, Z, UZ and UU.  It must be at
!>           least 1 and at least max( NN ).
!>           Not modified.
!>
!>  U      - DOUBLE PRECISION array, dimension (LDU,max(NN))
!>           The orthogonal matrix computed by DGEHRD.
!>           Modified.
!>
!>  Z      - DOUBLE PRECISION array, dimension (LDU,max(NN))
!>           The orthogonal matrix computed by DHSEQR.
!>           Modified.
!>
!>  UZ     - DOUBLE PRECISION array, dimension (LDU,max(NN))
!>           The product of U times Z.
!>           Modified.
!>
!>  WR1    - DOUBLE PRECISION array, dimension (max(NN))
!>  WI1    - DOUBLE PRECISION array, dimension (max(NN))
!>           The real and imaginary parts of the eigenvalues of A,
!>           as computed when Z is computed.
!>           On exit, WR1 + WI1*i are the eigenvalues of the matrix in A.
!>           Modified.
!>
!>  WR2    - DOUBLE PRECISION array, dimension (max(NN))
!>  WI2    - DOUBLE PRECISION array, dimension (max(NN))
!>           The real and imaginary parts of the eigenvalues of A,
!>           as computed when T is computed but not Z.
!>           On exit, WR2 + WI2*i are the eigenvalues of the matrix in A.
!>           Modified.
!>
!>  WR3    - DOUBLE PRECISION array, dimension (max(NN))
!>  WI3    - DOUBLE PRECISION array, dimension (max(NN))
!>           Like WR1, WI1, these arrays contain the eigenvalues of A,
!>           but those computed when DHSEQR only computes the
!>           eigenvalues, i.e., not the Schur vectors and no more of the
!>           Schur form than is necessary for computing the
!>           eigenvalues.
!>           Modified.
!>
!>  EVECTL - DOUBLE PRECISION array, dimension (LDU,max(NN))
!>           The (upper triangular) left eigenvector matrix for the
!>           matrix in T1.  For complex conjugate pairs, the real part
!>           is stored in one row and the imaginary part in the next.
!>           Modified.
!>
!>  EVEZTR - DOUBLE PRECISION array, dimension (LDU,max(NN))
!>           The (upper triangular) right eigenvector matrix for the
!>           matrix in T1.  For complex conjugate pairs, the real part
!>           is stored in one column and the imaginary part in the next.
!>           Modified.
!>
!>  EVECTY - DOUBLE PRECISION array, dimension (LDU,max(NN))
!>           The left eigenvector matrix for the
!>           matrix in H.  For complex conjugate pairs, the real part
!>           is stored in one row and the imaginary part in the next.
!>           Modified.
!>
!>  EVECTX - DOUBLE PRECISION array, dimension (LDU,max(NN))
!>           The right eigenvector matrix for the
!>           matrix in H.  For complex conjugate pairs, the real part
!>           is stored in one column and the imaginary part in the next.
!>           Modified.
!>
!>  UU     - DOUBLE PRECISION array, dimension (LDU,max(NN))
!>           Details of the orthogonal matrix computed by DGEHRD.
!>           Modified.
!>
!>  TAU    - DOUBLE PRECISION array, dimension(max(NN))
!>           Further details of the orthogonal matrix computed by DGEHRD.
!>           Modified.
!>
!>  WORK   - DOUBLE PRECISION array, dimension (NWORK)
!>           Workspace.
!>           Modified.
!>
!>  NWORK  - INTEGER
!>           The number of entries in WORK.  NWORK >= 4*NN(j)*NN(j) + 2.
!>
!>  IWORK  - INTEGER array, dimension (max(NN))
!>           Workspace.
!>           Modified.
!>
!>  SELECT - LOGICAL array, dimension (max(NN))
!>           Workspace.
!>           Modified.
!>
!>  RESULT - DOUBLE PRECISION array, dimension (16)
!>           The values computed by the fourteen tests described above.
!>           The values are currently limited to 1/ulp, to avoid
!>           overflow.
!>           Modified.
!>
!>  INFO   - INTEGER
!>           If 0, then everything ran OK.
!>            -1: NSIZES < 0
!>            -2: Some NN(j) < 0
!>            -3: NTYPES < 0
!>            -6: THRESH < 0
!>            -9: LDA < 1 or LDA < NMAX, where NMAX is max( NN(j) ).
!>           -14: LDU < 1 or LDU < NMAX.
!>           -28: NWORK too small.
!>           If  DLATMR, SLATMS, or SLATME returns an error code, the
!>               absolute value of it is returned.
!>           If 1, then DHSEQR could not find all the shifts.
!>           If 2, then the EISPACK code (for small blocks) failed.
!>           If >2, then 30*N iterations were not enough to find an
!>               eigenvalue or to decompose the problem.
!>           Modified.
!>
!>-----------------------------------------------------------------------
!>
!>     Some Local Variables and Parameters:
!>     ---- ----- --------- --- ----------
!>
!>     0.0D0, 1.0D0       Real 0 and 1.
!>     MAXTYP          The number of types defined.
!>     MTEST           The number of tests defined: care must be taken
!>                     that (1) the size of RESULT, (2) the number of
!>                     tests actually performed, and (3) MTEST agree.
!>     NTEST           The number of tests performed on this matrix
!>                     so far.  This should be less than MTEST, and
!>                     equal to it by the last test.  It will be less
!>                     if any of the routines being tested indicates
!>                     that it could not compute the matrices that
!>                     would be tested.
!>     NMAX            Largest value in NN.
!>     NMATS           The number of matrices generated so far.
!>     NERRS           The number of tests which have exceeded THRESH
!>                     so far (computed by DLAFTS).
!>     COND, CONDS,
!>     IMODE           Values to be passed to the matrix generators.
!>     ANORM           Norm of A; passed to matrix generators.
!>
!>     OVFL, UNFL      Overflow and underflow thresholds.
!>     ULP, ULPINV     Finest relative precision and its inverse.
!>     RTOVFL, RTUNFL,
!>     RTULP, RTULPI   Square roots of the previous 4 values.
!>
!>             The following four arrays decode JTYPE:
!>     KTYPE(j)        The general type (1-10) for type "j".
!>     KMODE(j)        The MODE value to be passed to the matrix
!>                     generator for type "j".
!>     KMAGN(j)        The order of magnitude ( O(1),
!>                     O(overflow^(1/2) ), O(underflow^(1/2) )
!>     KCONDS(j)       Selects whether CONDS is to be 1 or
!>                     1/sqrt(ulp).  (0 means irrelevant.)
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
!> \ingroup double_eig
!
!  =====================================================================
   SUBROUTINE DCHKHS( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, &
                      NOUNIT, A, LDA, H, T1, T2, U, LDU, Z, UZ, WR1, &
                      WI1, WR2, WI2, WR3, WI3, EVECTL, EVECTR, &
                      EVECTY, EVECTX, UU, TAU, WORK, NWORK, IWORK, &
                      SELECT, RESULT, INFO )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, LDA, LDU, NOUNIT, NSIZES, NTYPES, NWORK
   DOUBLE PRECISION   THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * ), SELECT( * )
   INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
   DOUBLE PRECISION   A( LDA, * ), EVECTL( LDU, * ), &
                      EVECTR( LDU, * ), EVECTX( LDU, * ), &
                      EVECTY( LDU, * ), H( LDA, * ), RESULT( 16 ), &
                      T1( LDA, * ), T2( LDA, * ), TAU( * ), &
                      U( LDU, * ), UU( LDU, * ), UZ( LDU, * ), &
                      WI1( * ), WI2( * ), WI3( * ), WORK( * ), &
                      WR1( * ), WR2( * ), WR3( * ), Z( LDU, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            MAXTYP
   PARAMETER          ( MAXTYP = 21 )
!     ..
!     .. Local Scalars ..
   LOGICAL            BADNN, MATCH
   INTEGER            I, IHI, IINFO, ILO, IMODE, IN, ITYPE, J, JCOL, &
                      JJ, JSIZE, JTYPE, K, MTYPES, N, N1, NERRS, &
                      NMATS, NMAX, NSELC, NSELR, NTEST, NTESTT
   DOUBLE PRECISION   ANINV, ANORM, COND, CONDS, OVFL, RTOVFL, RTULP, &
                      RTULPI, RTUNFL, TEMP1, TEMP2, ULP, ULPINV, UNFL
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          ADUMMA( 1 )
   INTEGER            IDUMMA( 1 ), IOLDSD( 4 ), KCONDS( MAXTYP ), &
                      KMAGN( MAXTYP ), KMODE( MAXTYP ), &
                      KTYPE( MAXTYP )
   DOUBLE PRECISION   DUMMA( 6 )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH
   EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           DCOPY, DGEHRD, DGEMM, DGET10, DGET22, DHSEIN, &
                      DHSEQR, DHST01, DLACPY, DLAFTS, DLASET, DLASUM, &
                      DLATME, DLATMR, DLATMS, DORGHR, DORMHR, DTREVC, &
                      DTREVC3, XERBLA
!     ..
!     .. Data statements ..
   DATA               KTYPE / 1, 2, 3, 5*4, 4*6, 6*6, 3*9 /
   DATA               KMAGN / 3*1, 1, 1, 1, 2, 3, 4*1, 1, 1, 1, 1, 2, 3, 1, 2, 3 /
   DATA               KMODE / 3*0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1 /
   DATA               KCONDS / 3*0, 5*0, 4*1, 6*2, 3*0 /
!     ..
!     .. Executable Statements ..
!
!     Check for errors
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
   ELSE IF( THRESH < 0.0D0 ) THEN
      INFO = -6
   ELSE IF( LDA <= 1 .OR. LDA < NMAX ) THEN
      INFO = -9
   ELSE IF( LDU <= 1 .OR. LDU < NMAX ) THEN
      INFO = -14
   ELSE IF( 4*NMAX*NMAX+2 > NWORK ) THEN
      INFO = -28
   END IF
!
   IF( INFO /= 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'DCHKHS', -INFO )
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
!     Quick return if possible
!
   IF( NSIZES == 0 .OR. NTYPES == 0 ) RETURN
!
!     More important constants
!
   UNFL = DLAMCH( 'Safe minimum' )
   OVFL = DLAMCH( 'Overflow' )
   ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
   ULPINV = 1.0D0 / ULP
   RTUNFL = SQRT( UNFL )
   RTOVFL = SQRT( OVFL )
   RTULP = SQRT( ULP )
   RTULPI = 1.0D0 / RTULP
!
!     Loop over sizes, types
!
   NERRS = 0
   NMATS = 0
!
   DO JSIZE = 1, NSIZES
      N = NN( JSIZE )
      IF( N == 0 ) GO TO 270
      N1 = MAX( 1, N )
      ANINV = 1.0D0 / DBLE( N1 )
!
      IF( NSIZES /= 1 ) THEN
         MTYPES = MIN( MAXTYP, NTYPES )
      ELSE
         MTYPES = MIN( MAXTYP+1, NTYPES )
      END IF
!
      DO JTYPE = 1, MTYPES
         IF( .NOT.DOTYPE( JTYPE ) ) GO TO 260
         NMATS = NMATS + 1
         NTEST = 0
!
!           Save ISEED in case of an error.
!
         IOLDSD(1:4) = ISEED(1:4)
!
!           Initialize RESULT
!
         RESULT(1:14) = 0.0D+0
!
!           Compute "A"
!
!           Control parameters:
!
!           KMAGN  KCONDS  KMODE        KTYPE
!       =1  O(1)   1       clustered 1  zero
!       =2  large  large   clustered 2  identity
!       =3  small          exponential  Jordan
!       =4                 arithmetic   diagonal, (w/ eigenvalues)
!       =5                 random log   symmetric, w/ eigenvalues
!       =6                 random       general, w/ eigenvalues
!       =7                              random diagonal
!       =8                              random symmetric
!       =9                              random general
!       =10                             random triangular
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
           ANORM = 1.0D+0
          CASE (2)
           ANORM = ( RTOVFL*ULP )*ANINV
          CASE (3)
           ANORM = RTUNFL*N*ULPINV
         END SELECT
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DLASET( 'Full', LDA, N, 0.0D0, 0.0D0, A, LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IINFO = 0
         COND = ULPINV
!
!           Special Matrices
!
         IF( ITYPE == 1 ) THEN
!
!              Zero
!
            IINFO = 0
         ELSE IF( ITYPE == 2 ) THEN
!
!              Identity
!
            DO JCOL = 1, N
               A( JCOL, JCOL ) = ANORM
            ENDDO
!
         ELSE IF( ITYPE == 3 ) THEN
!
!              Jordan Block
!
            DO JCOL = 1, N
               A( JCOL, JCOL ) = ANORM
               IF( JCOL > 1 ) A( JCOL, JCOL-1 ) = 1.0D0
            ENDDO
!
         ELSE IF( ITYPE == 4 ) THEN
!
!              Diagonal Matrix, [Eigen]values Specified
!
            CALL DLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, &
                         ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ), &
                         IINFO )
!
         ELSE IF( ITYPE == 5 ) THEN
!
!              Symmetric, eigenvalues specified
!
            CALL DLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, &
                         ANORM, N, N, 'N', A, LDA, WORK( N+1 ), &
                         IINFO )
!
         ELSE IF( ITYPE == 6 ) THEN
!
!              General, eigenvalues specified
!
            IF( KCONDS( JTYPE ) == 1 ) THEN
               CONDS = 1.0D0
            ELSE IF( KCONDS( JTYPE ) == 2 ) THEN
               CONDS = RTULPI
            ELSE
               CONDS = 0.0D0
            END IF
!
            ADUMMA( 1 ) = ' '
            CALL DLATME( N, 'S', ISEED, WORK, IMODE, COND, 1.0D0, &
                         ADUMMA, 'T', 'T', 'T', WORK( N+1 ), 4, &
                         CONDS, N, N, ANORM, A, LDA, WORK( 2*N+1 ), &
                         IINFO )
!
         ELSE IF( ITYPE == 7 ) THEN
!
!              Diagonal, random eigenvalues
!
            CALL DLATMR( N, N, 'S', ISEED, 'S', WORK, 6, 1.0D0, 1.0D0, &
                         'T', 'N', WORK( N+1 ), 1, 1.0D0, &
                         WORK( 2*N+1 ), 1, 1.0D0, 'N', IDUMMA, 0, 0, &
                         0.0D0, ANORM, 'NO', A, LDA, IWORK, IINFO )
!
         ELSE IF( ITYPE == 8 ) THEN
!
!              Symmetric, random eigenvalues
!
            CALL DLATMR( N, N, 'S', ISEED, 'S', WORK, 6, 1.0D0, 1.0D0, &
                         'T', 'N', WORK( N+1 ), 1, 1.0D0, &
                         WORK( 2*N+1 ), 1, 1.0D0, 'N', IDUMMA, N, N, &
                         0.0D0, ANORM, 'NO', A, LDA, IWORK, IINFO )
!
         ELSE IF( ITYPE == 9 ) THEN
!
!              General, random eigenvalues
!
            CALL DLATMR( N, N, 'S', ISEED, 'N', WORK, 6, 1.0D0, 1.0D0, &
                         'T', 'N', WORK( N+1 ), 1, 1.0D0, &
                         WORK( 2*N+1 ), 1, 1.0D0, 'N', IDUMMA, N, N, &
                         0.0D0, ANORM, 'NO', A, LDA, IWORK, IINFO )
!
         ELSE IF( ITYPE == 10 ) THEN
!
!              Triangular, random eigenvalues
!
            CALL DLATMR( N, N, 'S', ISEED, 'N', WORK, 6, 1.0D0, 1.0D0, &
                         'T', 'N', WORK( N+1 ), 1, 1.0D0, &
                         WORK( 2*N+1 ), 1, 1.0D0, 'N', IDUMMA, N, 0, &
                         0.0D0, ANORM, 'NO', A, LDA, IWORK, IINFO )
!
         ELSE
!
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
  100       CONTINUE
!
!           Call DGEHRD to compute H and U, do tests.
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DLACPY( ' ', N, N, A, LDA, H, LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
         NTEST = 1
!
         ILO = 1
         IHI = N
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DGEHRD( N, ILO, IHI, H, LDA, WORK, WORK( N+1 ), &
                      NWORK-N, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DGEHRD : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
         IF( IINFO /= 0 ) THEN
            RESULT( 1 ) = ULPINV
            WRITE( NOUNIT, FMT = 9999 )'DGEHRD', IINFO, N, JTYPE, &
               IOLDSD
            INFO = ABS( IINFO )
            GO TO 250
         END IF
!
         DO J = 1, N - 1
            UU(J+1, J ) = 0.0D+0
            UU(J+2:N,J) = H(J+2:N,J)
            U(J+2:N,J) = H(J+2:N,J)
            H(J+2:N,J) = 0.0D+0
            ENDDO
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DCOPY( N-1, WORK, 1, TAU, 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DORGHR( N, ILO, IHI, U, LDU, WORK, WORK( N+1 ), NWORK-N, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DORGHR : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         NTEST = 2
!
         CALL DHST01( N, ILO, IHI, A, LDA, H, LDA, U, LDU, WORK, &
                      NWORK, RESULT( 1 ) )
!
!           Call DHSEQR to compute T1, T2 and Z, do tests.
!
!           Eigenvalues only (WR3,WI3)
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DLACPY( ' ', N, N, H, LDA, T2, LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         NTEST = 3
         RESULT( 3 ) = ULPINV
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DHSEQR( 'E', 'N', N, ILO, IHI, T2, LDA, WR3, WI3, UZ, &
                      LDU, WORK, NWORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DHSEQR : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( IINFO /= 0 ) THEN
            WRITE( NOUNIT, FMT = 9999 )'DHSEQR(E)', IINFO, N, JTYPE, &
               IOLDSD
            IF( IINFO <= N+2 ) THEN
               INFO = ABS( IINFO )
               GO TO 250
            END IF
         END IF
!
!           Eigenvalues (WR2,WI2) and Full Schur Form (T2)
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DLACPY( ' ', N, N, H, LDA, T2, LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DHSEQR( 'S', 'N', N, ILO, IHI, T2, LDA, WR2, WI2, UZ, &
                      LDU, WORK, NWORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DHSEQR : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( IINFO /= 0 .AND. IINFO <= N+2 ) THEN
            WRITE( NOUNIT, FMT = 9999 )'DHSEQR(S)', IINFO, N, JTYPE, &
               IOLDSD
            INFO = ABS( IINFO )
            GO TO 250
         END IF
!
!           Eigenvalues (WR1,WI1), Schur Form (T1), and Schur vectors
!           (UZ)
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DLACPY( ' ', N, N, H, LDA, T1, LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DLACPY( ' ', N, N, U, LDU, UZ, LDU )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DHSEQR( 'S', 'V', N, ILO, IHI, T1, LDA, WR1, WI1, UZ, &
                      LDU, WORK, NWORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DHSEQR : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( IINFO /= 0 .AND. IINFO <= N+2 ) THEN
            WRITE( NOUNIT, FMT = 9999 )'DHSEQR(V)', IINFO, N, JTYPE, &
               IOLDSD
            INFO = ABS( IINFO )
            GO TO 250
         END IF
!
!           Compute Z = U' UZ
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DGEMM( 'T', 'N', N, N, N, 1.0D0, U, LDU, UZ, LDU, 0.0D0, &
                     Z, LDU )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DGEMM : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         NTEST = 8
!
!           Do Tests 3: | H - Z T Z' | / ( |H| n ulp )
!                and 4: | I - Z Z' | / ( n ulp )
!
         CALL DHST01( N, ILO, IHI, H, LDA, T1, LDA, Z, LDU, WORK, &
                      NWORK, RESULT( 3 ) )
!
!           Do Tests 5: | A - UZ T (UZ)' | / ( |A| n ulp )
!                and 6: | I - UZ (UZ)' | / ( n ulp )
!
         CALL DHST01( N, ILO, IHI, A, LDA, T1, LDA, UZ, LDU, WORK, &
                      NWORK, RESULT( 5 ) )
!
!           Do Test 7: | T2 - T1 | / ( |T| n ulp )
!
         CALL DGET10( N, N, T2, LDA, T1, LDA, WORK, RESULT( 7 ) )
!
!           Do Test 8: | W2 - W1 | / ( max(|W1|,|W2|) ulp )
!
         TEMP1 = MAX(MAXVAL(ABS(WR1(1:N))+ABS(WI1(1:N))),MAXVAL(ABS(WR2(1:N))+ABS(WI2(1:N))))
         TEMP2 = MAXVAL(ABS(WR1(1:N)-WR2(1:N))+ABS(WI1(1:N)-WI2(1:N)))
!
         RESULT( 8 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
!
!           Compute the Left and Right Eigenvectors of T
!
!           Compute the Right eigenvector Matrix:
!
         NTEST = 9
         RESULT( 9 ) = ULPINV
!
!           Select last max(N/4,1) real, max(N/4,1) complex eigenvectors
!
         NSELC = 0
         NSELR = 0
         J = N
  140       CONTINUE
         IF( WI1( J ) == 0.0D0 ) THEN
            IF( NSELR < MAX( N / 4, 1 ) ) THEN
               NSELR = NSELR + 1
               SELECT( J ) = .TRUE.
            ELSE
               SELECT( J ) = .FALSE.
            END IF
            J = J - 1
         ELSE
            IF( NSELC < MAX( N / 4, 1 ) ) THEN
               NSELC = NSELC + 1
               SELECT( J ) = .TRUE.
               SELECT( J-1 ) = .FALSE.
            ELSE
               SELECT( J ) = .FALSE.
               SELECT( J-1 ) = .FALSE.
            END IF
            J = J - 2
         END IF
         IF( J > 0 ) &
            GO TO 140
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DTREVC( 'Right', 'All', SELECT, N, T1, LDA, DUMMA, LDU, &
                      EVECTR, LDU, N, IN, WORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DTREVC : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( IINFO /= 0 ) THEN
            WRITE( NOUNIT, FMT = 9999 )'DTREVC(R,A)', IINFO, N, &
               JTYPE, IOLDSD
            INFO = ABS( IINFO )
            GO TO 250
         END IF
!
!           Test 9:  | TR - RW | / ( |T| |R| ulp )
!
         CALL DGET22( 'N', 'N', 'N', N, T1, LDA, EVECTR, LDU, WR1, &
                      WI1, WORK, DUMMA( 1 ) )
         RESULT( 9 ) = DUMMA( 1 )
         IF( DUMMA( 2 ) > THRESH ) THEN
            WRITE( NOUNIT, FMT = 9998 )'Right', 'DTREVC', &
               DUMMA( 2 ), N, JTYPE, IOLDSD
         END IF
!
!           Compute selected right eigenvectors and confirm that
!           they agree with previous right eigenvectors
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DTREVC( 'Right', 'Some', SELECT, N, T1, LDA, DUMMA, &
                      LDU, EVECTL, LDU, N, IN, WORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DTREVC : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( IINFO /= 0 ) THEN
            WRITE( NOUNIT, FMT = 9999 )'DTREVC(R,S)', IINFO, N, &
               JTYPE, IOLDSD
            INFO = ABS( IINFO )
            GO TO 250
         END IF
!
         K = 1
         MATCH = .TRUE.
         DO J = 1, N
            IF( SELECT( J ) .AND. WI1( J ) == 0.0D0 ) THEN
               IF (ANY(EVECTR(1:N,J) /= EVECTL(1:N,K))) THEN
                  MATCH = .FALSE.
                  GO TO 180
               END IF
               K = K + 1
            ELSE IF( SELECT( J ) .AND. WI1( J ) /= 0.0D0 ) THEN
               IF (ANY(EVECTR(1:N,J) /= EVECTL(1:N,K)) .or. ANY(EVECTR(1:N, J+1 ) /= EVECTL(1:N,K+1))) THEN
                  MATCH = .FALSE.
                  GO TO 180
               END IF
               K = K + 2
            END IF
            ENDDO
  180       CONTINUE
         IF( .NOT.MATCH ) WRITE( NOUNIT, FMT = 9997 )'Right', 'DTREVC', N, JTYPE, IOLDSD
!
!           Compute the Left eigenvector Matrix:
!
         NTEST = 10
         RESULT( 10 ) = ULPINV
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DTREVC( 'Left', 'All', SELECT, N, T1, LDA, EVECTL, LDU, &
                      DUMMA, LDU, N, IN, WORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DTREVC : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( IINFO /= 0 ) THEN
            WRITE( NOUNIT, FMT = 9999 )'DTREVC(L,A)', IINFO, N, JTYPE, IOLDSD
            INFO = ABS( IINFO )
            GO TO 250
         END IF
!
!           Test 10:  | LT - WL | / ( |T| |L| ulp )
!
         CALL DGET22( 'Trans', 'N', 'Conj', N, T1, LDA, EVECTL, LDU, &
                      WR1, WI1, WORK, DUMMA( 3 ) )
         RESULT( 10 ) = DUMMA( 3 )
         IF( DUMMA( 4 ) > THRESH ) THEN
            WRITE( NOUNIT, FMT = 9998 )'Left', 'DTREVC', DUMMA( 4 ), &
               N, JTYPE, IOLDSD
         END IF
!
!           Compute selected left eigenvectors and confirm that
!           they agree with previous left eigenvectors
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DTREVC( 'Left', 'Some', SELECT, N, T1, LDA, EVECTR, &
                      LDU, DUMMA, LDU, N, IN, WORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DTREVC : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( IINFO /= 0 ) THEN
            WRITE( NOUNIT, FMT = 9999 )'DTREVC(L,S)', IINFO, N, &
               JTYPE, IOLDSD
            INFO = ABS( IINFO )
            GO TO 250
         END IF
!
         K = 1
         MATCH = .TRUE.
         DO J = 1, N
            IF( SELECT( J ) .AND. WI1( J ) == 0.0D0 ) THEN
               IF (ANY(EVECTL(1:N,J) /= EVECTR(1:N,K))) THEN
                  MATCH = .FALSE.
                  GO TO 220
               END IF

               K = K + 1
            ELSE IF( SELECT( J ) .AND. WI1( J ) /= 0.0D0 ) THEN
               IF (ANY(EVECTL(1:N,J) /= EVECTR(1:N,K)) .or. ANY(EVECTL(1:N,J+1) /= EVECTR(1:N,K+1))) THEN
                  MATCH = .FALSE.
                  GO TO 220
               END IF
               K = K + 2
            END IF
            ENDDO
  220       CONTINUE
         IF( .NOT.MATCH ) WRITE( NOUNIT, FMT = 9997 )'Left', 'DTREVC', N, JTYPE, &
            IOLDSD
!
!           Call DHSEIN for Right eigenvectors of H, do test 11
!
         NTEST = 11
         RESULT( 11 ) = ULPINV
         SELECT(1:N) = .TRUE.
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DHSEIN( 'Right', 'Qr', 'Ninitv', SELECT, N, H, LDA, &
                      WR3, WI3, DUMMA, LDU, EVECTX, LDU, N1, IN, &
                      WORK, IWORK, IWORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DHSEIN : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( IINFO /= 0 ) THEN
            WRITE( NOUNIT, FMT = 9999 )'DHSEIN(R)', IINFO, N, JTYPE, &
               IOLDSD
            INFO = ABS( IINFO )
            IF( IINFO < 0 ) GO TO 250
         ELSE
!
!              Test 11:  | HX - XW | / ( |H| |X| ulp )
!
!                        (from inverse iteration)
!
            CALL DGET22( 'N', 'N', 'N', N, H, LDA, EVECTX, LDU, WR3, &
                         WI3, WORK, DUMMA( 1 ) )
            IF( DUMMA( 1 ) < ULPINV ) &
               RESULT( 11 ) = DUMMA( 1 )*ANINV
            IF( DUMMA( 2 ) > THRESH ) THEN
               WRITE( NOUNIT, FMT = 9998 )'Right', 'DHSEIN', &
                  DUMMA( 2 ), N, JTYPE, IOLDSD
            END IF
         END IF
!
!           Call DHSEIN for Left eigenvectors of H, do test 12
!
         NTEST = 12
         RESULT( 12 ) = ULPINV
         SELECT(1:N) = .TRUE.
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DHSEIN( 'Left', 'Qr', 'Ninitv', SELECT, N, H, LDA, WR3, &
                      WI3, EVECTY, LDU, DUMMA, LDU, N1, IN, WORK, &
                      IWORK, IWORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DHSEIN : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( IINFO /= 0 ) THEN
            WRITE( NOUNIT, FMT = 9999 )'DHSEIN(L)', IINFO, N, JTYPE, &
               IOLDSD
            INFO = ABS( IINFO )
            IF( IINFO < 0 ) GO TO 250
         ELSE
!
!              Test 12:  | YH - WY | / ( |H| |Y| ulp )
!
!                        (from inverse iteration)
!
            CALL DGET22( 'C', 'N', 'C', N, H, LDA, EVECTY, LDU, WR3, &
                         WI3, WORK, DUMMA( 3 ) )
            IF( DUMMA( 3 ) < ULPINV ) &
               RESULT( 12 ) = DUMMA( 3 )*ANINV
            IF( DUMMA( 4 ) > THRESH ) THEN
               WRITE( NOUNIT, FMT = 9998 )'Left', 'DHSEIN', &
                  DUMMA( 4 ), N, JTYPE, IOLDSD
            END IF
         END IF
!
!           Call DORMHR for Right eigenvectors of A, do test 13
!
         NTEST = 13
         RESULT( 13 ) = ULPINV
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DORMHR( 'Left', 'No transpose', N, N, ILO, IHI, UU, &
                      LDU, TAU, EVECTX, LDU, WORK, NWORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DORMHR : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( IINFO /= 0 ) THEN
            WRITE( NOUNIT, FMT = 9999 )'DORMHR(R)', IINFO, N, JTYPE, &
               IOLDSD
            INFO = ABS( IINFO )
            IF( IINFO < 0 ) GO TO 250
         ELSE
!
!              Test 13:  | AX - XW | / ( |A| |X| ulp )
!
!                        (from inverse iteration)
!
            CALL DGET22( 'N', 'N', 'N', N, A, LDA, EVECTX, LDU, WR3, &
                         WI3, WORK, DUMMA( 1 ) )
            IF( DUMMA( 1 ) < ULPINV ) RESULT( 13 ) = DUMMA( 1 )*ANINV
         END IF
!
!           Call DORMHR for Left eigenvectors of A, do test 14
!
         NTEST = 14
         RESULT( 14 ) = ULPINV
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DORMHR( 'Left', 'No transpose', N, N, ILO, IHI, UU, &
                      LDU, TAU, EVECTY, LDU, WORK, NWORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DORMHR : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( IINFO /= 0 ) THEN
            WRITE( NOUNIT, FMT = 9999 )'DORMHR(L)', IINFO, N, JTYPE, &
               IOLDSD
            INFO = ABS( IINFO )
            IF( IINFO < 0 ) GO TO 250
         ELSE
!
!              Test 14:  | YA - WY | / ( |A| |Y| ulp )
!
!                        (from inverse iteration)
!
            CALL DGET22( 'C', 'N', 'C', N, A, LDA, EVECTY, LDU, WR3, &
                         WI3, WORK, DUMMA( 3 ) )
            IF( DUMMA( 3 ) < ULPINV ) &
               RESULT( 14 ) = DUMMA( 3 )*ANINV
         END IF
!
!           Compute Left and Right Eigenvectors of A
!
!           Compute a Right eigenvector matrix:
!
         NTEST = 15
         RESULT( 15 ) = ULPINV
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DLACPY( ' ', N, N, UZ, LDU, EVECTR, LDU )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DTREVC3( 'Right', 'Back', SELECT, N, T1, LDA, DUMMA, &
                       LDU, EVECTR, LDU, N, IN, WORK, NWORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DTREVC3 : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( IINFO /= 0 ) THEN
            WRITE( NOUNIT, FMT = 9999 )'DTREVC3(R,B)', IINFO, N, &
               JTYPE, IOLDSD
            INFO = ABS( IINFO )
            GO TO 250
         END IF
!
!           Test 15:  | AR - RW | / ( |A| |R| ulp )
!
!                     (from Schur decomposition)
!
         CALL DGET22( 'N', 'N', 'N', N, A, LDA, EVECTR, LDU, WR1, &
                      WI1, WORK, DUMMA( 1 ) )
         RESULT( 15 ) = DUMMA( 1 )
         IF( DUMMA( 2 ) > THRESH ) THEN
            WRITE( NOUNIT, FMT = 9998 )'Right', 'DTREVC3', &
               DUMMA( 2 ), N, JTYPE, IOLDSD
         END IF
!
!           Compute a Left eigenvector matrix:
!
         NTEST = 16
         RESULT( 16 ) = ULPINV
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DLACPY( ' ', N, N, UZ, LDU, EVECTL, LDU )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DTREVC3( 'Left', 'Back', SELECT, N, T1, LDA, EVECTL, &
                       LDU, DUMMA, LDU, N, IN, WORK, NWORK, IINFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DTREVC3 : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( IINFO /= 0 ) THEN
            WRITE( NOUNIT, FMT = 9999 )'DTREVC3(L,B)', IINFO, N, &
               JTYPE, IOLDSD
            INFO = ABS( IINFO )
            GO TO 250
         END IF
!
!           Test 16:  | LA - WL | / ( |A| |L| ulp )
!
!                     (from Schur decomposition)
!
         CALL DGET22( 'Trans', 'N', 'Conj', N, A, LDA, EVECTL, LDU, &
                      WR1, WI1, WORK, DUMMA( 3 ) )
         RESULT( 16 ) = DUMMA( 3 )
         IF( DUMMA( 4 ) > THRESH ) THEN
            WRITE( NOUNIT, FMT = 9998 )'Left', 'DTREVC3', DUMMA( 4 ), &
               N, JTYPE, IOLDSD
         END IF
!
!           End of Loop -- Check for RESULT(j) > THRESH
!
  250       CONTINUE
!
         NTESTT = NTESTT + NTEST
         CALL DLAFTS( 'DHS', N, N, JTYPE, NTEST, RESULT, IOLDSD, &
                      THRESH, NOUNIT, NERRS )
!
  260    CONTINUE
         ENDDO
  270 CONTINUE
      ENDDO
!
!     Summary
!
   CALL DLASUM( 'DHS', NOUNIT, NERRS, NTESTT )
!
   RETURN
!
 9999 FORMAT( ' DCHKHS: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', &
         I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
 9998 FORMAT( ' DCHKHS: ', A, ' Eigenvectors from ', A, ' incorrectly ', &
         'normalized.', / ' Bits of error=', 0P, G10.3, ',', 9X, &
         'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, &
         ')' )
 9997 FORMAT( ' DCHKHS: Selected ', A, ' Eigenvectors from ', A, &
         ' do not match other eigenvectors ', 9X, 'N=', I6, &
         ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
!
!     End of DCHKHS
!
END



