!> \brief \b ZDRVSG2STG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZDRVSG2STG( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                              NOUNIT, A, LDA, B, LDB, D, D2, Z, LDZ, AB,
!                              BB, AP, BP, WORK, NWORK, RWORK, LRWORK,
!                              IWORK, LIWORK, RESULT, INFO )
!
!       IMPLICIT N1.0D0
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB, LDZ, LIWORK, LRWORK, NOUNIT,
!      $                   NSIZES, NTYPES, NWORK
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
!       DOUBLE PRECISION   D( * ), RESULT( * ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), AB( LDA, * ), AP( * ),
!      $                   B( LDB, * ), BB( LDB, * ), BP( * ), WORK( * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      ZDRVSG2STG checks the complex Hermitian generalized eigenproblem
!>      drivers.
!>
!>              ZHEGV computes all eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian-definite generalized
!>              eigenproblem.
!>
!>              ZHEGVD computes all eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian-definite generalized
!>              eigenproblem using a divide and conquer algorithm.
!>
!>              ZHEGVX computes selected eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian-definite generalized
!>              eigenproblem.
!>
!>              ZHPGV computes all eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian-definite generalized
!>              eigenproblem in packed storage.
!>
!>              ZHPGVD computes all eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian-definite generalized
!>              eigenproblem in packed storage using a divide and
!>              conquer algorithm.
!>
!>              ZHPGVX computes selected eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian-definite generalized
!>              eigenproblem in packed storage.
!>
!>              ZHBGV computes all eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian-definite banded
!>              generalized eigenproblem.
!>
!>              ZHBGVD computes all eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian-definite banded
!>              generalized eigenproblem using a divide and conquer
!>              algorithm.
!>
!>              ZHBGVX computes selected eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian-definite banded
!>              generalized eigenproblem.
!>
!>      When ZDRVSG2STG is called, a number of matrix "sizes" ("n's") and a
!>      number of matrix "types" are specified.  For each size ("n")
!>      and each type of matrix, one matrix A of the given type will be
!>      generated; a random well-conditioned matrix B is also generated
!>      and the pair (A,B) is used to test the drivers.
!>
!>      For each pair (A,B), the following tests are performed:
!>
!>      (1) ZHEGV with ITYPE = 1 and UPLO ='U':
!>
!>              | A Z - B Z D | / ( |A| |Z| n ulp )
!>              | D - D2 | / ( |D| ulp )   where D is computed by
!>                                         ZHEGV and  D2 is computed by
!>                                         ZHEGV_2STAGE. This test is
!>                                         only performed for DSYGV
!>
!>      (2) as (1) but calling ZHPGV
!>      (3) as (1) but calling ZHBGV
!>      (4) as (1) but with UPLO = 'L'
!>      (5) as (4) but calling ZHPGV
!>      (6) as (4) but calling ZHBGV
!>
!>      (7) ZHEGV with ITYPE = 2 and UPLO ='U':
!>
!>              | A B Z - Z D | / ( |A| |Z| n ulp )
!>
!>      (8) as (7) but calling ZHPGV
!>      (9) as (7) but with UPLO = 'L'
!>      (10) as (9) but calling ZHPGV
!>
!>      (11) ZHEGV with ITYPE = 3 and UPLO ='U':
!>
!>              | B A Z - Z D | / ( |A| |Z| n ulp )
!>
!>      (12) as (11) but calling ZHPGV
!>      (13) as (11) but with UPLO = 'L'
!>      (14) as (13) but calling ZHPGV
!>
!>      ZHEGVD, ZHPGVD and ZHBGVD performed the same 14 tests.
!>
!>      ZHEGVX, ZHPGVX and ZHBGVX performed the above 14 tests with
!>      the parameter RANGE = 'A', 'N' and 'I', respectively.
!>
!>      The "sizes" are specified by an array NN(1:NSIZES); the value of
!>      each element NN(j) specifies one size.
!>      The "types" are specified by a logical array DOTYPE( 1:NTYPES );
!>      if DOTYPE(j) is .TRUE., then matrix type "j" will be generated.
!>      This type is used for the matrix A which has half-bandwidth KA.
!>      B is generated as a well-conditioned positive definite matrix
!>      with half-bandwidth KB (<= KA).
!>      Currently, the list of possible types for A is:
!>
!>      (1)  The zero matrix.
!>      (2)  The identity matrix.
!>
!>      (3)  A diagonal matrix with evenly spaced entries
!>           1, ..., ULP  and random signs.
!>           (ULP = (first number larger than 1) - 1 )
!>      (4)  A diagonal matrix with geometrically spaced entries
!>           1, ..., ULP  and random signs.
!>      (5)  A diagonal matrix with "clustered" entries 1, ULP, ..., ULP
!>           and random signs.
!>
!>      (6)  Same as (4), but multiplied by SQRT( overflow threshold )
!>      (7)  Same as (4), but multiplied by SQRT( underflow threshold )
!>
!>      (8)  A matrix of the form  U* D U, where U is unitary and
!>           D has evenly spaced entries 1, ..., ULP with random signs
!>           on the diagonal.
!>
!>      (9)  A matrix of the form  U* D U, where U is unitary and
!>           D has geometrically spaced entries 1, ..., ULP with random
!>           signs on the diagonal.
!>
!>      (10) A matrix of the form  U* D U, where U is unitary and
!>           D has "clustered" entries 1, ULP,..., ULP with random
!>           signs on the diagonal.
!>
!>      (11) Same as (8), but multiplied by SQRT( overflow threshold )
!>      (12) Same as (8), but multiplied by SQRT( underflow threshold )
!>
!>      (13) Hermitian matrix with random entries chosen from (-1,1).
!>      (14) Same as (13), but multiplied by SQRT( overflow threshold )
!>      (15) Same as (13), but multiplied by SQRT( underflow threshold )
!>
!>      (16) Same as (8), but with KA = 1 and KB = 1
!>      (17) Same as (8), but with KA = 2 and KB = 1
!>      (18) Same as (8), but with KA = 2 and KB = 2
!>      (19) Same as (8), but with KA = 3 and KB = 1
!>      (20) Same as (8), but with KA = 3 and KB = 2
!>      (21) Same as (8), but with KA = 3 and KB = 3
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \verbatim
!>  NSIZES  INTEGER
!>          The number of sizes of matrices to use.  If it is zero,
!>          ZDRVSG2STG does nothing.  It must be at least zero.
!>          Not modified.
!>
!>  NN      INTEGER array, dimension (NSIZES)
!>          An array containing the sizes to be used for the matrices.
!>          Zero values will be skipped.  The values must be at least
!>          zero.
!>          Not modified.
!>
!>  NTYPES  INTEGER
!>          The number of elements in DOTYPE.   If it is zero, ZDRVSG2STG
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
!>          next call to ZDRVSG2STG to continue the same random number
!>          sequence.
!>          Modified.
!>
!>  THRESH  DOUBLE PRECISION
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
!>  A       COMPLEX*16 array, dimension (LDA , max(NN))
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
!>  B       COMPLEX*16 array, dimension (LDB , max(NN))
!>          Used to hold the Hermitian positive definite matrix for
!>          the generalized problem.
!>          On exit, B contains the last matrix actually
!>          used.
!>          Modified.
!>
!>  LDB     INTEGER
!>          The leading dimension of B.  It must be at
!>          least 1 and at least max( NN ).
!>          Not modified.
!>
!>  D       DOUBLE PRECISION array, dimension (max(NN))
!>          The eigenvalues of A. On exit, the eigenvalues in D
!>          correspond with the matrix in A.
!>          Modified.
!>
!>  Z       COMPLEX*16 array, dimension (LDZ, max(NN))
!>          The matrix of eigenvectors.
!>          Modified.
!>
!>  LDZ     INTEGER
!>          The leading dimension of ZZ.  It must be at least 1 and
!>          at least max( NN ).
!>          Not modified.
!>
!>  AB      COMPLEX*16 array, dimension (LDA, max(NN))
!>          Workspace.
!>          Modified.
!>
!>  BB      COMPLEX*16 array, dimension (LDB, max(NN))
!>          Workspace.
!>          Modified.
!>
!>  AP      COMPLEX*16 array, dimension (max(NN)**2)
!>          Workspace.
!>          Modified.
!>
!>  BP      COMPLEX*16 array, dimension (max(NN)**2)
!>          Workspace.
!>          Modified.
!>
!>  WORK    COMPLEX*16 array, dimension (NWORK)
!>          Workspace.
!>          Modified.
!>
!>  NWORK   INTEGER
!>          The number of entries in WORK.  This must be at least
!>          2*N + N**2  where  N = max( NN(j), 2 ).
!>          Not modified.
!>
!>  RWORK   DOUBLE PRECISION array, dimension (LRWORK)
!>          Workspace.
!>          Modified.
!>
!>  LRWORK  INTEGER
!>          The number of entries in RWORK.  This must be at least
!>          max( 7*N, 1 + 4*N + 2*N*lg(N) + 3*N**2 ) where
!>          N = max( NN(j) ) and lg( N ) = smallest integer k such
!>          that 2**k >= N .
!>          Not modified.
!>
!>  IWORK   INTEGER array, dimension (LIWORK))
!>          Workspace.
!>          Modified.
!>
!>  LIWORK  INTEGER
!>          The number of entries in IWORK.  This must be at least
!>          2 + 5*max( NN(j) ).
!>          Not modified.
!>
!>  RESULT  DOUBLE PRECISION array, dimension (70)
!>          The values computed by the 70 tests described above.
!>          Modified.
!>
!>  INFO    INTEGER
!>          If 0, then everything ran OK.
!>           -1: NSIZES < 0
!>           -2: Some NN(j) < 0
!>           -3: NTYPES < 0
!>           -5: THRESH < 0
!>           -9: LDA < 1 or LDA < NMAX, where NMAX is max( NN(j) ).
!>          -16: LDZ < 1 or LDZ < NMAX.
!>          -21: NWORK too small.
!>          -23: LRWORK too small.
!>          -25: LIWORK too small.
!>          If  ZLATMR, CLATMS, ZHEGV, ZHPGV, ZHBGV, CHEGVD, CHPGVD,
!>              ZHPGVD, ZHEGVX, CHPGVX, ZHBGVX returns an error code,
!>              the absolute value of it is returned.
!>          Modified.
!>
!>-----------------------------------------------------------------------
!>
!>       Some Local Variables and Parameters:
!>       ---- ----- --------- --- ----------
!>       0.0D0, 1.0D0       Real 0 and 1.
!>       MAXTYP          The number of types defined.
!>       NTEST           The number of tests that have been run
!>                       on this matrix.
!>       NTESTT          The total number of tests for this call.
!>       NMAX            Largest value in NN.
!>       NMATS           The number of matrices generated so far.
!>       NERRS           The number of tests which have exceeded THRESH
!>                       so far (computed by DLAFTS).
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
!> \ingroup complex16_eig
!
!  =====================================================================
   SUBROUTINE ZDRVSG2STG( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, &
                          NOUNIT, A, LDA, B, LDB, D, D2, Z, LDZ, AB, &
                          BB, AP, BP, WORK, NWORK, RWORK, LRWORK, &
                          IWORK, LIWORK, RESULT, INFO )
!
   IMPLICIT NONE
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, LDA, LDB, LDZ, LIWORK, LRWORK, NOUNIT, &
                      NSIZES, NTYPES, NWORK
   DOUBLE PRECISION   THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
   DOUBLE PRECISION   D( * ), D2( * ), RESULT( * ), RWORK( * )
   COMPLEX*16         A( LDA, * ), AB( LDA, * ), AP( * ), &
                      B( LDB, * ), BB( LDB, * ), BP( * ), WORK( * ), &
                      Z( LDZ, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            MAXTYP
   PARAMETER          ( MAXTYP = 21 )
!     ..
!     .. Local Scalars ..
   LOGICAL            BADNN
   CHARACTER          UPLO
   INTEGER            I, IBTYPE, IBUPLO, IINFO, IJ, IL, IMODE, ITEMP, &
                      ITYPE, IU, J, JCOL, JSIZE, JTYPE, KA, KA9, KB, &
                      KB9, M, MTYPES, N, NERRS, NMATS, NMAX, NTEST, &
                      NTESTT
   DOUBLE PRECISION   ABSTOL, ANINV, ANORM, COND, OVFL, RTOVFL, &
                      RTUNFL, ULP, ULPINV, UNFL, VL, VU, TEMP1, TEMP2
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   INTEGER            IDUMMA( 1 ), IOLDSD( 4 ), ISEED2( 4 ), &
                      KMAGN( MAXTYP ), KMODE( MAXTYP ), &
                      KTYPE( MAXTYP )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   DOUBLE PRECISION   DLAMCH, DLARND
   EXTERNAL           LSAME, DLAMCH, DLARND
!     ..
!     .. External Subroutines ..
   EXTERNAL           DLAFTS, DLASUM, XERBLA, ZHBGV, ZHBGVD, &
                      ZHBGVX, ZHEGV, ZHEGVD, ZHEGVX, ZHPGV, ZHPGVD, &
                      ZHPGVX, ZLACPY, ZLASET, ZLATMR, ZLATMS, ZSGT01, &
                      ZHEGV_2STAGE

!     ..
!     .. Data statements ..
   DATA               KTYPE / 1, 2, 5*4, 5*5, 3*8, 6*9 /
   DATA               KMAGN / 2*1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 6*1 /
   DATA               KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 6*4 /
!     ..
!     .. Executable Statements ..
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
   ELSE IF( LDA <= 1 .OR. LDA < NMAX ) THEN
      INFO = -9
   ELSE IF( LDZ <= 1 .OR. LDZ < NMAX ) THEN
      INFO = -16
   ELSE IF( 2*MAX( NMAX, 2 )**2 > NWORK ) THEN
      INFO = -21
   ELSE IF( 2*MAX( NMAX, 2 )**2 > LRWORK ) THEN
      INFO = -23
   ELSE IF( 2*MAX( NMAX, 2 )**2 > LIWORK ) THEN
      INFO = -25
   END IF
!
   IF( INFO /= 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'ZDRVSG2STG', -INFO )
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
!     More Important constants
!
   UNFL = DLAMCH( 'Safe minimum' )
   OVFL = DLAMCH( 'Overflow' )
   ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
   ULPINV = 1.0D0 / ULP
   RTUNFL = SQRT( UNFL )
   RTOVFL = SQRT( OVFL )
!
   ISEED2(1:4) = ISEED(1:4)
!
!     Loop over sizes, types
!
   NERRS = 0
   NMATS = 0
!
   DO JSIZE = 1, NSIZES
      N = NN( JSIZE )
      ANINV = 1.0D0 / DBLE( MAX( 1, N ) )
!
      IF( NSIZES /= 1 ) THEN
         MTYPES = MIN( MAXTYP, NTYPES )
      ELSE
         MTYPES = MIN( MAXTYP+1, NTYPES )
      END IF
!
      KA9 = 0
      KB9 = 0
      DO JTYPE = 1, MTYPES
         IF (DOTYPE( JTYPE ) ) THEN
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
!           =4         arithmetic   diagonal, w/ eigenvalues
!           =5         random log   hermitian, w/ eigenvalues
!           =6         random       (none)
!           =7                      random diagonal
!           =8                      random hermitian
!           =9                      banded, w/ eigenvalues
!
         IF( MTYPES > MAXTYP ) GO TO 90
!
         ITYPE = KTYPE( JTYPE )
         IMODE = KMODE( JTYPE )
!
!           Compute norm
!
         SELECT CASE (KMAGN(JTYPE))
          CASE (1)
           ANORM = (1.0D+0,0.0D+0)
          CASE (2)
           ANORM = ( RTOVFL*ULP )*ANINV
          CASE (3)
           ANORM = RTUNFL*N*ULPINV
         END SELECT
!
         IINFO = 0
         COND = ULPINV
!
!           Special Matrices -- Identity & Jordan block
!
         IF( ITYPE == 1 ) THEN
!
!              Zero
!
            KA = 0
            KB = 0
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZLASET( 'Full', LDA, N, (0.0D+0,0.0D+0), (0.0D+0,0.0D+0), A, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
         ELSE IF( ITYPE == 2 ) THEN
!
!              Identity
!
            KA = 0
            KB = 0
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZLASET( 'Full', LDA, N, (0.0D+0,0.0D+0), (0.0D+0,0.0D+0), A, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            DO JCOL = 1, N
               A( JCOL, JCOL ) = ANORM
            ENDDO
!
         ELSE IF( ITYPE == 4 ) THEN
!
!              Diagonal Matrix, [Eigen]values Specified
!
            KA = 0
            KB = 0
            CALL ZLATMS( N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, &
                         ANORM, 0, 0, 'N', A, LDA, WORK, IINFO )
!
         ELSE IF( ITYPE == 5 ) THEN
!
!              Hermitian, eigenvalues specified
!
            KA = MAX( 0, N-1 )
            KB = KA
            CALL ZLATMS( N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, &
                         ANORM, N, N, 'N', A, LDA, WORK, IINFO )
!
         ELSE IF( ITYPE == 7 ) THEN
!
!              Diagonal, random eigenvalues
!
            KA = 0
            KB = 0
            CALL ZLATMR( N, N, 'S', ISEED, 'H', WORK, 6, 1.0D0, (1.0D+0,0.0D+0), &
                         'T', 'N', WORK( N+1 ), 1, 1.0D0, &
                         WORK( 2*N+1 ), 1, 1.0D0, 'N', IDUMMA, 0, 0, &
                         0.0D0, ANORM, 'NO', A, LDA, IWORK, IINFO )
!
         ELSE IF( ITYPE == 8 ) THEN
!
!              Hermitian, random eigenvalues
!
            KA = MAX( 0, N-1 )
            KB = KA
            CALL ZLATMR( N, N, 'S', ISEED, 'H', WORK, 6, 1.0D0, (1.0D+0,0.0D+0), &
                         'T', 'N', WORK( N+1 ), 1, 1.0D0, &
                         WORK( 2*N+1 ), 1, 1.0D0, 'N', IDUMMA, N, N, &
                         0.0D0, ANORM, 'NO', A, LDA, IWORK, IINFO )
!
         ELSE IF( ITYPE == 9 ) THEN
!
!              Hermitian banded, eigenvalues specified
!
!              The following values are used for the half-bandwidths:
!
!                ka = 1   kb = 1
!                ka = 2   kb = 1
!                ka = 2   kb = 2
!                ka = 3   kb = 1
!                ka = 3   kb = 2
!                ka = 3   kb = 3
!
            KB9 = KB9 + 1
            IF( KB9 > KA9 ) THEN
               KA9 = KA9 + 1
               KB9 = 1
            END IF
            KA = MAX( 0, MIN( N-1, KA9 ) )
            KB = MAX( 0, MIN( N-1, KB9 ) )
            CALL ZLATMS( N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, &
                         ANORM, KA, KA, 'N', A, LDA, WORK, IINFO )
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
90       CONTINUE
!
         ABSTOL = UNFL + UNFL
         IF( N <= 1 ) THEN
            IL = 1
            IU = N
         ELSE
            IL = 1 + INT( ( N-1 )*DLARND( 1, ISEED2 ) )
            IU = 1 + INT( ( N-1 )*DLARND( 1, ISEED2 ) )
            IF( IL > IU ) THEN
               ITEMP = IL
               IL = IU
               IU = ITEMP
            END IF
         END IF
!
!           3) Call ZHEGV, ZHPGV, ZHBGV, CHEGVD, CHPGVD, CHBGVD,
!              ZHEGVX, ZHPGVX and ZHBGVX, do tests.
!
!           loop over the three generalized problems
!                 IBTYPE = 1: A*x = (lambda)*B*x
!                 IBTYPE = 2: A*B*x = (lambda)*x
!                 IBTYPE = 3: B*A*x = (lambda)*x
!
         DO IBTYPE = 1, 3
!
!              loop over the setting UPLO
!
            DO IBUPLO = 1, 2
               IF( IBUPLO == 1 ) UPLO = 'U'
               IF( IBUPLO == 2 ) UPLO = 'L'
!
!                 Generate random well-conditioned positive definite
!                 matrix B, of bandwidth not greater than that of A.
!
               CALL ZLATMS( N, N, 'U', ISEED, 'P', RWORK, 5, 10.0D0, &
                            1.0D0, KB, KB, UPLO, B, LDB, WORK( N+1 ), &
                            IINFO )
!
!                 Test ZHEGV
!
               NTEST = NTEST + 1
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLACPY( ' ', N, N, A, LDA, Z, LDZ )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLACPY( UPLO, N, N, B, LDB, BB, LDB )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZHEGV( IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, &
                           WORK, NWORK, RWORK, IINFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZHEGV : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               IF( IINFO /= 0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'ZHEGV(V,' // UPLO // &
                     ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO < 0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 100
                  END IF
               END IF
!
!                 Do Test
!
               CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, &
                            LDZ, D, WORK, RWORK, RESULT( NTEST ) )
!
!                 Test ZHEGV_2STAGE
!
               NTEST = NTEST + 1
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLACPY( ' ', N, N, A, LDA, Z, LDZ )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLACPY( UPLO, N, N, B, LDB, BB, LDB )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZHEGV_2STAGE( IBTYPE, 'N', UPLO, N, Z, LDZ, &
                                  BB, LDB, D2, WORK, NWORK, RWORK, &
                                  IINFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZHEGV_2STAGE : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               IF( IINFO /= 0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 ) &
                     'ZHEGV_2STAGE(V,' // UPLO // &
                     ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO < 0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 100
                  END IF
               END IF
!
!                 Do Test
!
!                  CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z,
!     $                         LDZ, D, WORK, RWORK, RESULT( NTEST ) )
!
!                 Do Tests | D1 - D2 | / ( |D1| ulp )
!                 D1 computed using the standard 1-stage reduction as reference
!                 D2 computed using the 2-stage reduction
!
               TEMP1 = MAX(MAXVAL(ABS(D(1:N))), MAXVAL(ABS(D2(1:N))))
               TEMP2 = MAXVAL(ABS(D(1:N)-D2(1:N)))
               RESULT( NTEST ) = TEMP2 / &
                                 MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
!
!                 Test ZHEGVD
!
               NTEST = NTEST + 1
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLACPY( ' ', N, N, A, LDA, Z, LDZ )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLACPY( UPLO, N, N, B, LDB, BB, LDB )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZHEGVD( IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, &
                            WORK, NWORK, RWORK, LRWORK, IWORK, &
                            LIWORK, IINFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZHEGVD : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               IF( IINFO /= 0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'ZHEGVD(V,' // UPLO // &
                     ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO < 0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 100
                  END IF
               END IF
!
!                 Do Test
!
               CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, &
                            LDZ, D, WORK, RWORK, RESULT( NTEST ) )
!
!                 Test ZHEGVX
!
               NTEST = NTEST + 1
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLACPY( ' ', N, N, A, LDA, AB, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLACPY( UPLO, N, N, B, LDB, BB, LDB )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZHEGVX( IBTYPE, 'V', 'A', UPLO, N, AB, LDA, BB, &
                            LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, &
                            LDZ, WORK, NWORK, RWORK, IWORK( N+1 ), &
                            IWORK, IINFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZHEGVX : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               IF( IINFO /= 0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'ZHEGVX(V,A' // UPLO // &
                     ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO < 0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 100
                  END IF
               END IF
!
!                 Do Test
!
               CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, &
                            LDZ, D, WORK, RWORK, RESULT( NTEST ) )
!
               NTEST = NTEST + 1
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLACPY( ' ', N, N, A, LDA, AB, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLACPY( UPLO, N, N, B, LDB, BB, LDB )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 since we do not know the exact eigenvalues of this
!                 eigenpair, we just set VL and VU as constants.
!                 It is quite possible that there are no eigenvalues
!                 in this interval.
!
               VL = 0.0D0
               VU = ANORM
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZHEGVX( IBTYPE, 'V', 'V', UPLO, N, AB, LDA, BB, &
                            LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, &
                            LDZ, WORK, NWORK, RWORK, IWORK( N+1 ), &
                            IWORK, IINFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZHEGVX : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               IF( IINFO /= 0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'ZHEGVX(V,V,' // &
                     UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO < 0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 100
                  END IF
               END IF
!
!                 Do Test
!
               CALL ZSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, &
                            LDZ, D, WORK, RWORK, RESULT( NTEST ) )
!
               NTEST = NTEST + 1
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLACPY( ' ', N, N, A, LDA, AB, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLACPY( UPLO, N, N, B, LDB, BB, LDB )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZHEGVX( IBTYPE, 'V', 'I', UPLO, N, AB, LDA, BB, &
                            LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, &
                            LDZ, WORK, NWORK, RWORK, IWORK( N+1 ), &
                            IWORK, IINFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZHEGVX : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               IF( IINFO /= 0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'ZHEGVX(V,I,' // &
                     UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO < 0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 100
                  END IF
               END IF
!
!                 Do Test
!
               CALL ZSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, &
                            LDZ, D, WORK, RWORK, RESULT( NTEST ) )
!
  100             CONTINUE
!
!                 Test ZHPGV
!
               NTEST = NTEST + 1
!
!                 Copy the matrices into packed storage.
!
               IF( LSAME( UPLO, 'U' ) ) THEN
                  IJ = 1
                  DO J = 1, N
                     DO I = 1, J
                        AP( IJ ) = A( I, J )
                        BP( IJ ) = B( I, J )
                        IJ = IJ + 1
                     ENDDO
                  ENDDO
               ELSE
                  IJ = 1
                  DO J = 1, N
                     DO I = J, N
                        AP( IJ ) = A( I, J )
                        BP( IJ ) = B( I, J )
                        IJ = IJ + 1
                     ENDDO
                  ENDDO
               END IF
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZHPGV( IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, &
                           WORK, RWORK, IINFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZHPGV : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               IF( IINFO /= 0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'ZHPGV(V,' // UPLO // &
                     ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO < 0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 310
                  END IF
               END IF
!
!                 Do Test
!
               CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, &
                            LDZ, D, WORK, RWORK, RESULT( NTEST ) )
!
!                 Test ZHPGVD
!
               NTEST = NTEST + 1
!
!                 Copy the matrices into packed storage.
!
               IF( LSAME( UPLO, 'U' ) ) THEN
                  IJ = 1
                  DO J = 1, N
                     DO I = 1, J
                        AP( IJ ) = A( I, J )
                        BP( IJ ) = B( I, J )
                        IJ = IJ + 1
                     ENDDO
                  ENDDO
               ELSE
                  IJ = 1
                  DO J = 1, N
                     DO I = J, N
                        AP( IJ ) = A( I, J )
                        BP( IJ ) = B( I, J )
                        IJ = IJ + 1
                     ENDDO
                  ENDDO
               END IF
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZHPGVD( IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, &
                            WORK, NWORK, RWORK, LRWORK, IWORK, &
                            LIWORK, IINFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZHPGVD : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               IF( IINFO /= 0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'ZHPGVD(V,' // UPLO // &
                     ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO < 0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 310
                  END IF
               END IF
!
!                 Do Test
!
               CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, &
                            LDZ, D, WORK, RWORK, RESULT( NTEST ) )
!
!                 Test ZHPGVX
!
               NTEST = NTEST + 1
!
!                 Copy the matrices into packed storage.
!
               IF( LSAME( UPLO, 'U' ) ) THEN
                  IJ = 1
                  DO J = 1, N
                     DO I = 1, J
                        AP( IJ ) = A( I, J )
                        BP( IJ ) = B( I, J )
                        IJ = IJ + 1
                     ENDDO
                  ENDDO
               ELSE
                  IJ = 1
                  DO J = 1, N
                     DO I = J, N
                        AP( IJ ) = A( I, J )
                        BP( IJ ) = B( I, J )
                        IJ = IJ + 1
                     ENDDO
                  ENDDO
               END IF
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZHPGVX( IBTYPE, 'V', 'A', UPLO, N, AP, BP, VL, &
                            VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, &
                            RWORK, IWORK( N+1 ), IWORK, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZHPGVX : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               IF( IINFO /= 0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'ZHPGVX(V,A' // UPLO // &
                     ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO < 0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 310
                  END IF
               END IF
!
!                 Do Test
!
               CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, &
                            LDZ, D, WORK, RWORK, RESULT( NTEST ) )
!
               NTEST = NTEST + 1
!
!                 Copy the matrices into packed storage.
!
               IF( LSAME( UPLO, 'U' ) ) THEN
                  IJ = 1
                  DO J = 1, N
                     DO I = 1, J
                        AP( IJ ) = A( I, J )
                        BP( IJ ) = B( I, J )
                        IJ = IJ + 1
                     ENDDO
                  ENDDO
               ELSE
                  IJ = 1
                  DO J = 1, N
                     DO I = J, N
                        AP( IJ ) = A( I, J )
                        BP( IJ ) = B( I, J )
                        IJ = IJ + 1
                     ENDDO
                  ENDDO
               END IF
!
               VL = 0.0D0
               VU = ANORM
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZHPGVX( IBTYPE, 'V', 'V', UPLO, N, AP, BP, VL, &
                            VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, &
                            RWORK, IWORK( N+1 ), IWORK, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZHPGVX : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               IF( IINFO /= 0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'ZHPGVX(V,V' // UPLO // &
                     ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO < 0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 310
                  END IF
               END IF
!
!                 Do Test
!
               CALL ZSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, &
                            LDZ, D, WORK, RWORK, RESULT( NTEST ) )
!
               NTEST = NTEST + 1
!
!                 Copy the matrices into packed storage.
!
               IF( LSAME( UPLO, 'U' ) ) THEN
                  IJ = 1
                  DO J = 1, N
                     DO I = 1, J
                        AP( IJ ) = A( I, J )
                        BP( IJ ) = B( I, J )
                        IJ = IJ + 1
                     ENDDO
                  ENDDO
               ELSE
                  IJ = 1
                  DO J = 1, N
                     DO I = J, N
                        AP( IJ ) = A( I, J )
                        BP( IJ ) = B( I, J )
                        IJ = IJ + 1
                     ENDDO
                  ENDDO
               END IF
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZHPGVX( IBTYPE, 'V', 'I', UPLO, N, AP, BP, VL, &
                            VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, &
                            RWORK, IWORK( N+1 ), IWORK, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZHPGVX : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               IF( IINFO /= 0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'ZHPGVX(V,I' // UPLO // &
                     ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO < 0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 310
                  END IF
               END IF
!
!                 Do Test
!
               CALL ZSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, &
                            LDZ, D, WORK, RWORK, RESULT( NTEST ) )
!
  310             CONTINUE
!
               IF( IBTYPE == 1 ) THEN
!
!                    TEST ZHBGV
!
                  NTEST = NTEST + 1
!
!                    Copy the matrices into band storage.
!
                  IF( LSAME( UPLO, 'U' ) ) THEN
                     DO J = 1, N
                        DO I = MAX( 1, J-KA ), J
                           AB( KA+1+I-J, J ) = A( I, J )
                        ENDDO
                        DO I = MAX( 1, J-KB ), J
                           BB( KB+1+I-J, J ) = B( I, J )
                        ENDDO
                     ENDDO
                  ELSE
                     DO J = 1, N
                        DO I = J, MIN( N, J+KA )
                           AB( 1+I-J, J ) = A( I, J )
                        ENDDO
                        DO I = J, MIN( N, J+KB )
                           BB( 1+I-J, J ) = B( I, J )
                        ENDDO
                     ENDDO
                  END IF
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZHBGV( 'V', UPLO, N, KA, KB, AB, LDA, BB, LDB, &
                              D, Z, LDZ, WORK, RWORK, IINFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZHBGV : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  IF( IINFO /= 0 ) THEN
                     WRITE( NOUNIT, FMT = 9999 )'ZHBGV(V,' // &
                        UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     IF( IINFO < 0 ) THEN
                        RETURN
                     ELSE
                        RESULT( NTEST ) = ULPINV
                        GO TO 620
                     END IF
                  END IF
!
!                    Do Test
!
                  CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, &
                               LDZ, D, WORK, RWORK, RESULT( NTEST ) )
!
!                    TEST ZHBGVD
!
                  NTEST = NTEST + 1
!
!                    Copy the matrices into band storage.
!
                  IF( LSAME( UPLO, 'U' ) ) THEN
                     DO J = 1, N
                        DO I = MAX( 1, J-KA ), J
                           AB( KA+1+I-J, J ) = A( I, J )
                        ENDDO
                        DO I = MAX( 1, J-KB ), J
                           BB( KB+1+I-J, J ) = B( I, J )
                        ENDDO
                     ENDDO
                  ELSE
                     DO J = 1, N
                        DO I = J, MIN( N, J+KA )
                           AB( 1+I-J, J ) = A( I, J )
                        ENDDO
                        DO I = J, MIN( N, J+KB )
                           BB( 1+I-J, J ) = B( I, J )
                        ENDDO
                     ENDDO
                  END IF
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZHBGVD( 'V', UPLO, N, KA, KB, AB, LDA, BB, &
                               LDB, D, Z, LDZ, WORK, NWORK, RWORK, &
                               LRWORK, IWORK, LIWORK, IINFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZHBGVD : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  IF( IINFO /= 0 ) THEN
                     WRITE( NOUNIT, FMT = 9999 )'ZHBGVD(V,' // &
                        UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     IF( IINFO < 0 ) THEN
                        RETURN
                     ELSE
                        RESULT( NTEST ) = ULPINV
                        GO TO 620
                     END IF
                  END IF
!
!                    Do Test
!
                  CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, &
                               LDZ, D, WORK, RWORK, RESULT( NTEST ) )
!
!                    Test ZHBGVX
!
                  NTEST = NTEST + 1
!
!                    Copy the matrices into band storage.
!
                  IF( LSAME( UPLO, 'U' ) ) THEN
                     DO J = 1, N
                        DO I = MAX( 1, J-KA ), J
                           AB( KA+1+I-J, J ) = A( I, J )
                        ENDDO
                        DO I = MAX( 1, J-KB ), J
                           BB( KB+1+I-J, J ) = B( I, J )
                        ENDDO
                     ENDDO
                  ELSE
                     DO J = 1, N
                        DO I = J, MIN( N, J+KA )
                           AB( 1+I-J, J ) = A( I, J )
                        ENDDO
                        DO I = J, MIN( N, J+KB )
                           BB( 1+I-J, J ) = B( I, J )
                        ENDDO
                     ENDDO
                  END IF
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZHBGVX( 'V', 'A', UPLO, N, KA, KB, AB, LDA, &
                               BB, LDB, BP, MAX( 1, N ), VL, VU, IL, &
                               IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, &
                               IWORK( N+1 ), IWORK, IINFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZHBGVX : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  IF( IINFO /= 0 ) THEN
                     WRITE( NOUNIT, FMT = 9999 )'ZHBGVX(V,A' // &
                        UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     IF( IINFO < 0 ) THEN
                        RETURN
                     ELSE
                        RESULT( NTEST ) = ULPINV
                        GO TO 620
                     END IF
                  END IF
!
!                    Do Test
!
                  CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, &
                               LDZ, D, WORK, RWORK, RESULT( NTEST ) )
!
                  NTEST = NTEST + 1
!
!                    Copy the matrices into band storage.
!
                  IF( LSAME( UPLO, 'U' ) ) THEN
                     DO J = 1, N
                        DO I = MAX( 1, J-KA ), J
                           AB( KA+1+I-J, J ) = A( I, J )
                        ENDDO
                        DO I = MAX( 1, J-KB ), J
                           BB( KB+1+I-J, J ) = B( I, J )
                        ENDDO
                     ENDDO
                  ELSE
                     DO J = 1, N
                        DO I = J, MIN( N, J+KA )
                           AB( 1+I-J, J ) = A( I, J )
                        ENDDO
                        DO I = J, MIN( N, J+KB )
                           BB( 1+I-J, J ) = B( I, J )
                        ENDDO
                     ENDDO
                  END IF
!
                  VL = 0.0D0
                  VU = ANORM
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZHBGVX( 'V', 'V', UPLO, N, KA, KB, AB, LDA, &
                               BB, LDB, BP, MAX( 1, N ), VL, VU, IL, &
                               IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, &
                               IWORK( N+1 ), IWORK, IINFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZHBGVX : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  IF( IINFO /= 0 ) THEN
                     WRITE( NOUNIT, FMT = 9999 )'ZHBGVX(V,V' // &
                        UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     IF( IINFO < 0 ) THEN
                        RETURN
                     ELSE
                        RESULT( NTEST ) = ULPINV
                        GO TO 620
                     END IF
                  END IF
!
!                    Do Test
!
                  CALL ZSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, &
                               LDZ, D, WORK, RWORK, RESULT( NTEST ) )
!
                  NTEST = NTEST + 1
!
!                    Copy the matrices into band storage.
!
                  IF( LSAME( UPLO, 'U' ) ) THEN
                     DO J = 1, N
                        DO I = MAX( 1, J-KA ), J
                           AB( KA+1+I-J, J ) = A( I, J )
                        ENDDO
                        DO I = MAX( 1, J-KB ), J
                           BB( KB+1+I-J, J ) = B( I, J )
                        ENDDO
                     ENDDO
                  ELSE
                     DO J = 1, N
                        DO I = J, MIN( N, J+KA )
                           AB( 1+I-J, J ) = A( I, J )
                        ENDDO
                        DO I = J, MIN( N, J+KB )
                           BB( 1+I-J, J ) = B( I, J )
                        ENDDO
                     ENDDO
                  END IF
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZHBGVX( 'V', 'I', UPLO, N, KA, KB, AB, LDA, &
                               BB, LDB, BP, MAX( 1, N ), VL, VU, IL, &
                               IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, &
                               IWORK( N+1 ), IWORK, IINFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZHBGVX : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  IF( IINFO /= 0 ) THEN
                     WRITE( NOUNIT, FMT = 9999 )'ZHBGVX(V,I' // &
                        UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     IF( IINFO < 0 ) THEN
                        RETURN
                     ELSE
                        RESULT( NTEST ) = ULPINV
                        GO TO 620
                     END IF
                  END IF
!
!                    Do Test
!
                  CALL ZSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, &
                               LDZ, D, WORK, RWORK, RESULT( NTEST ) )
!
               END IF
!
  620          CONTINUE
               ENDDO
            ENDDO
!
!           End of Loop -- Check for RESULT(j) > THRESH
!
         NTESTT = NTESTT + NTEST
         CALL DLAFTS( 'ZSG', N, N, JTYPE, NTEST, RESULT, IOLDSD, &
                      THRESH, NOUNIT, NERRS )
         ENDIF
         ENDDO
      ENDDO
!
!     Summary
!
   CALL DLASUM( 'ZSG', NOUNIT, NERRS, NTESTT )
!
   RETURN
!
 9999 FORMAT( ' ZDRVSG2STG: ', A, ' returned INFO=', I6, '.', / 9X, &
     'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
!
!     End of ZDRVSG2STG
!
END


