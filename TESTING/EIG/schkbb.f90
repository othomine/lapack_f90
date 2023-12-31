!> \brief \b SCHKBB
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SCHKBB( NSIZES, MVAL, NVAL, NWDTHS, KK, NTYPES, DOTYPE,
!                          NRHS, ISEED, THRESH, NOUNIT, A, LDA, AB, LDAB,
!                          BD, BE, Q, LDQ, P, LDP, C, LDC, CC, WORK,
!                          LWORK, RESULT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDAB, LDC, LDP, LDQ, LWORK, NOUNIT,
!      $                   NRHS, NSIZES, NTYPES, NWDTHS
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), KK( * ), MVAL( * ), NVAL( * )
!       REAL               A( LDA, * ), AB( LDAB, * ), BD( * ), BE( * ),
!      $                   C( LDC, * ), CC( LDC, * ), P( LDP, * ),
!      $                   Q( LDQ, * ), RESULT( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SCHKBB tests the reduction of a general real rectangular band
!> matrix to bidiagonal form.
!>
!> SGBBRD factors a general band matrix A as  Q B P* , where * means
!> transpose, B is upper bidiagonal, and Q and P are orthogonal;
!> SGBBRD can also overwrite a given matrix C with Q* C .
!>
!> For each pair of matrix dimensions (M,N) and each selected matrix
!> type, an M by N matrix A and an M by NRHS matrix C are generated.
!> The problem dimensions are as follows
!>    A:          M x N
!>    Q:          M x M
!>    P:          N x N
!>    B:          min(M,N) x min(M,N)
!>    C:          M x NRHS
!>
!> For each generated matrix, 4 tests are performed:
!>
!> (1)   | A - Q B PT | / ( |A| max(M,N) ulp ), PT = P'
!>
!> (2)   | I - Q' Q | / ( M ulp )
!>
!> (3)   | I - PT PT' | / ( N ulp )
!>
!> (4)   | Y - Q' C | / ( |Y| max(M,NRHS) ulp ), where Y = Q' C.
!>
!> The "types" are specified by a logical array DOTYPE( 1:NTYPES );
!> if DOTYPE(j) is .TRUE., then matrix type "j" will be generated.
!> Currently, the list of possible types is:
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
!>          If NSIZES is zero, SCHKBB does nothing.  NSIZES must be at
!>          least zero.
!> \endverbatim
!>
!> \param[in] MVAL
!> \verbatim
!>          MVAL is INTEGER array, dimension (NSIZES)
!>          The values of the matrix row dimension M.
!> \endverbatim
!>
!> \param[in] NVAL
!> \verbatim
!>          NVAL is INTEGER array, dimension (NSIZES)
!>          The values of the matrix column dimension N.
!> \endverbatim
!>
!> \param[in] NWDTHS
!> \verbatim
!>          NWDTHS is INTEGER
!>          The number of bandwidths to use.  If it is zero,
!>          SCHKBB does nothing.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] KK
!> \verbatim
!>          KK is INTEGER array, dimension (NWDTHS)
!>          An array containing the bandwidths to be used for the band
!>          matrices.  The values must be at least zero.
!> \endverbatim
!>
!> \param[in] NTYPES
!> \verbatim
!>          NTYPES is INTEGER
!>          The number of elements in DOTYPE.   If it is zero, SCHKBB
!>          does nothing.  It must be at least zero.  If it is MAXTYP+1
!>          and NSIZES is 1, then an additional type, MAXTYP+1 is
!>          defined, which is to use whatever matrix is in A.  This
!>          is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
!>          DOTYPE(MAXTYP+1) is .TRUE. .
!> \endverbatim
!>
!> \param[in] DOTYPE
!> \verbatim
!>          DOTYPE is LOGICAL array, dimension (NTYPES)
!>          If DOTYPE(j) is .TRUE., then for each size in NN a
!>          matrix of that size and of type j will be generated.
!>          If NTYPES is smaller than the maximum number of types
!>          defined (PARAMETER MAXTYP), then types NTYPES+1 through
!>          MAXTYP will not be generated.  If NTYPES is larger
!>          than MAXTYP, DOTYPE(MAXTYP+1) through DOTYPE(NTYPES)
!>          will be ignored.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns in the "right-hand side" matrix C.
!>          If NRHS = 0, then the operations on the right-hand side will
!>          not be tested. NRHS must be at least 0.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry ISEED specifies the seed of the random number
!>          generator. The array elements should be between 0 and 4095;
!>          if not they will be reduced mod 4096.  Also, ISEED(4) must
!>          be odd.  The random number generator uses a linear
!>          congruential sequence limited to small integers, and so
!>          should produce machine independent random numbers. The
!>          values of ISEED are changed on exit, and can be used in the
!>          next call to SCHKBB to continue the same random number
!>          sequence.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is REAL
!>          A test will count as "failed" if the "error", computed as
!>          described above, exceeds THRESH.  Note that the error
!>          is scaled to be O(1), so THRESH should be a reasonably
!>          small multiple of 1, e.g., 10 or 100.  In particular,
!>          it should not depend on the precision (single vs. double)
!>          or the size of the matrix.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] NOUNIT
!> \verbatim
!>          NOUNIT is INTEGER
!>          The FORTRAN unit number for printing out error messages
!>          (e.g., if a routine returns IINFO not equal to 0.)
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension
!>                            (LDA, max(NN))
!>          Used to hold the matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  It must be at least 1
!>          and at least max( NN ).
!> \endverbatim
!>
!> \param[out] AB
!> \verbatim
!>          AB is REAL array, dimension (LDAB, max(NN))
!>          Used to hold A in band storage format.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of AB.  It must be at least 2 (not 1!)
!>          and at least max( KK )+1.
!> \endverbatim
!>
!> \param[out] BD
!> \verbatim
!>          BD is REAL array, dimension (max(NN))
!>          Used to hold the diagonal of the bidiagonal matrix computed
!>          by SGBBRD.
!> \endverbatim
!>
!> \param[out] BE
!> \verbatim
!>          BE is REAL array, dimension (max(NN))
!>          Used to hold the off-diagonal of the bidiagonal matrix
!>          computed by SGBBRD.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDQ, max(NN))
!>          Used to hold the orthogonal matrix Q computed by SGBBRD.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of Q.  It must be at least 1
!>          and at least max( NN ).
!> \endverbatim
!>
!> \param[out] P
!> \verbatim
!>          P is REAL array, dimension (LDP, max(NN))
!>          Used to hold the orthogonal matrix P computed by SGBBRD.
!> \endverbatim
!>
!> \param[in] LDP
!> \verbatim
!>          LDP is INTEGER
!>          The leading dimension of P.  It must be at least 1
!>          and at least max( NN ).
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is REAL array, dimension (LDC, max(NN))
!>          Used to hold the matrix C updated by SGBBRD.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of U.  It must be at least 1
!>          and at least max( NN ).
!> \endverbatim
!>
!> \param[out] CC
!> \verbatim
!>          CC is REAL array, dimension (LDC, max(NN))
!>          Used to hold a copy of the matrix C.
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
!>          max( LDA+1, max(NN)+1 )*max(NN).
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (4)
!>          The values computed by the tests described above.
!>          The values are currently limited to 1/ulp, to avoid
!>          overflow.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          If 0, then everything ran OK.
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
!>                       so far.
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
!> \ingroup single_eig
!
!  =====================================================================
   SUBROUTINE SCHKBB( NSIZES, MVAL, NVAL, NWDTHS, KK, NTYPES, DOTYPE, &
                      NRHS, ISEED, THRESH, NOUNIT, A, LDA, AB, LDAB, &
                      BD, BE, Q, LDQ, P, LDP, C, LDC, CC, WORK, &
                      LWORK, RESULT, INFO )
!
!  -- LAPACK test routine (input) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, LDA, LDAB, LDC, LDP, LDQ, LWORK, NOUNIT, &
                      NRHS, NSIZES, NTYPES, NWDTHS
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            ISEED( 4 ), KK( * ), MVAL( * ), NVAL( * )
   REAL               A( LDA, * ), AB( LDAB, * ), BD( * ), BE( * ), &
                      C( LDC, * ), CC( LDC, * ), P( LDP, * ), &
                      Q( LDQ, * ), RESULT( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            MAXTYP
   PARAMETER          ( MAXTYP = 15 )
!     ..
!     .. Local Scalars ..
   LOGICAL            BADMM, BADNN, BADNNB
   INTEGER            I, IINFO, IMODE, ITYPE, J, JCOL, JR, JSIZE, &
                      JTYPE, JWIDTH, K, KL, KMAX, KU, M, MMAX, MNMAX, &
                      MNMIN, MTYPES, N, NERRS, NMATS, NMAX, NTEST, &
                      NTESTT
   REAL               AMNINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, ULP, &
                      ULPINV, UNFL
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   INTEGER            IDUMMA( 1 ), IOLDSD( 4 ), KMAGN( MAXTYP ), &
                      KMODE( MAXTYP ), KTYPE( MAXTYP )
!     ..
!     .. External Functions ..
   REAL               SLAMCH
   EXTERNAL           SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           SBDT01, SBDT02, SGBBRD, SLACPY, SLAHD2, SLASET, &
                      SLASUM, SLATMR, SLATMS, SORT01, XERBLA
!     ..
!     .. Data statements ..
   DATA               KTYPE / 1, 2, 5*4, 5*6, 3*9 /
   DATA               KMAGN / 2*1, 3*1, 2, 3, 3*1, 2, 3, 1, 2, 3 /
   DATA               KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0 /
!     ..
!     .. Executable Statements ..
!
!     Check for errors
!
   NTESTT = 0
   INFO = 0
!
!     Important constants
!
   BADMM = ANY(MVAL(1:NSIZES) < 0)
   BADNN = ANY(NVAL(1:NSIZES) < 0)
   MMAX = MAX(1, MAXVAL(MVAL(1:NSIZES)))
   NMAX = MAX(1, MAXVAL(NVAL(1:NSIZES)))
   MNMAX = 1
   DO J = 1, NSIZES
      MNMAX = MAX( MNMAX, MIN( MVAL( J ), NVAL( J ) ) )
   ENDDO
!
   KMAX = MAXVAL(KK(1:NWDTHS))
   BADNNB = (ANY(KK(1:NWDTHS) < 0 ))
!
!     Check for errors
!
   IF( NSIZES < 0 ) THEN
      INFO = -1
   ELSE IF( BADMM ) THEN
      INFO = -2
   ELSE IF( BADNN ) THEN
      INFO = -3
   ELSE IF( NWDTHS < 0 ) THEN
      INFO = -4
   ELSE IF( BADNNB ) THEN
      INFO = -5
   ELSE IF( NTYPES < 0 ) THEN
      INFO = -6
   ELSE IF( NRHS < 0 ) THEN
      INFO = -8
   ELSE IF( LDA < NMAX ) THEN
      INFO = -13
   ELSE IF( LDAB < 2*KMAX+1 ) THEN
      INFO = -15
   ELSE IF( LDQ < NMAX ) THEN
      INFO = -19
   ELSE IF( LDP < NMAX ) THEN
      INFO = -21
   ELSE IF( LDC < NMAX ) THEN
      INFO = -23
   ELSE IF( ( MAX( LDA, NMAX )+1 )*NMAX > LWORK ) THEN
      INFO = -26
   END IF
!
   IF( INFO /= 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'SCHKBB', -INFO )
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
   IF( NSIZES == 0 .OR. NTYPES == 0 .OR. NWDTHS == 0 ) RETURN
!
!     More Important constants
!
   UNFL = SLAMCH( 'Safe minimum' )
   OVFL = 1.0E+0 / UNFL
   ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
   ULPINV = 1.0E+0 / ULP
   RTUNFL = SQRT( UNFL )
   RTOVFL = SQRT( OVFL )
!
!     Loop over sizes, widths, types
!
   NERRS = 0
   NMATS = 0
!
   DO JSIZE = 1, NSIZES
      M = MVAL( JSIZE )
      N = NVAL( JSIZE )
      MNMIN = MIN( M, N )
      AMNINV = 1.0E+0 / REAL( MAX( 1, M, N ) )
!
      DO JWIDTH = 1, NWDTHS
         K = KK( JWIDTH )
         IF( K < M .OR. K < N ) THEN
         KL = MAX( 0, MIN( M-1, K ) )
         KU = MAX( 0, MIN( N-1, K ) )
!
         IF( NSIZES /= 1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         ELSE
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF
!
         DO JTYPE = 1, MTYPES
            IF(DOTYPE( JTYPE ) ) THEN
            NMATS = NMATS + 1
            NTEST = 0
!
            IOLDSD(1:4) = ISEED(1:4)
!
!              Compute "A".
!
!              Control parameters:
!
!                  KMAGN  KMODE        KTYPE
!              =1  O(1)   clustered 1  zero
!              =2  large  clustered 2  identity
!              =3  small  exponential  (none)
!              =4         arithmetic   diagonal, (w/ singular values)
!              =5         random log   (none)
!              =6         random       nonhermitian, w/ singular values
!              =7                      (none)
!              =8                      (none)
!              =9                      random nonhermitian
!
            IF( MTYPES <= MAXTYP ) THEN
!
            ITYPE = KTYPE( JTYPE )
            IMODE = KMODE( JTYPE )
!
!              Compute norm
!
            SELECT CASE (KMAGN( JTYPE ))
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
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SLASET( 'Full', LDAB, N, 0.0E+0, 0.0E+0, AB, LDAB )
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
!              Special Matrices -- Identity & Jordan block
!
!                 Zero
!
            IF( ITYPE == 1 ) THEN
               IINFO = 0
!
            ELSE IF( ITYPE == 2 ) THEN
!
!                 Identity
!
               FORALL (JCOL = 1:N) A( JCOL, JCOL ) = ANORM
!
            ELSE IF( ITYPE == 4 ) THEN
!
!                 Diagonal Matrix, singular values specified
!
               CALL SLATMS( M, N, 'S', ISEED, 'N', WORK, IMODE, COND, &
                            ANORM, 0, 0, 'N', A, LDA, WORK( M+1 ), &
                            IINFO )
!
            ELSE IF( ITYPE == 6 ) THEN
!
!                 Nonhermitian, singular values specified
!
               CALL SLATMS( M, N, 'S', ISEED, 'N', WORK, IMODE, COND, &
                            ANORM, KL, KU, 'N', A, LDA, WORK( M+1 ), &
                            IINFO )
!
            ELSE IF( ITYPE == 9 ) THEN
!
!                 Nonhermitian, random entries
!
               CALL SLATMR( M, N, 'S', ISEED, 'N', WORK, 6, 1.0E+0, 1.0E+0, &
                            'T', 'N', WORK( N+1 ), 1, 1.0E+0, &
                            WORK( 2*N+1 ), 1, 1.0E+0, 'N', IDUMMA, KL, &
                            KU, 0.0E+0, ANORM, 'N', A, LDA, IDUMMA, &
                            IINFO )
!
            ELSE
!
               IINFO = 1
            END IF
!
!              Generate Right-Hand Side
!
            CALL SLATMR( M, NRHS, 'S', ISEED, 'N', WORK, 6, 1.0E+0, 1.0E+0, &
                         'T', 'N', WORK( M+1 ), 1, 1.0E+0, &
                         WORK( 2*M+1 ), 1, 1.0E+0, 'N', IDUMMA, M, NRHS, &
                         0.0E+0, 1.0E+0, 'NO', C, LDC, IDUMMA, IINFO )
!
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
!
            ENDIF
!
!              Copy A to band storage.
!
            DO J = 1, N
               DO I = MAX( 1, J-KU ), MIN( M, J+KL )
                  AB( KU+1+I-J, J ) = A( I, J )
               ENDDO
            ENDDO
!
!              Copy C
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SLACPY( 'Full', M, NRHS, C, LDC, CC, LDC )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
!              Call SGBBRD to compute B, Q and P, and to update C.
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGBBRD( 'B', M, N, NRHS, KL, KU, AB, LDAB, BD, BE, &
                         Q, LDQ, P, LDP, CC, LDC, WORK, IINFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGBBRD : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SGBBRD', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 1 ) = ULPINV
                  GO TO 120
               END IF
            END IF
!
!              Test 1:  Check the decomposition A := Q * B * P'
!                   2:  Check the orthogonality of Q
!                   3:  Check the orthogonality of P
!                   4:  Check the computation of Q' * C
!
            CALL SBDT01( M, N, -1, A, LDA, Q, LDQ, BD, BE, P, LDP, WORK, RESULT( 1 ) )
            CALL SORT01( 'Columns', M, M, Q, LDQ, WORK, LWORK, RESULT( 2 ) )
            CALL SORT01( 'Rows', N, N, P, LDP, WORK, LWORK, RESULT( 3 ) )
            CALL SBDT02( M, NRHS, C, LDC, CC, LDC, Q, LDQ, WORK, RESULT( 4 ) )
!
!              End of Loop -- Check for RESULT(j) > THRESH
!
            NTEST = 4
  120          CONTINUE
            NTESTT = NTESTT + NTEST
!
!              Print out tests which fail.
!
            DO JR = 1, NTEST
               IF( RESULT( JR ) >= THRESH ) THEN
                  IF( NERRS == 0 ) CALL SLAHD2( NOUNIT, 'SBB' )
                  NERRS = NERRS + 1
                  WRITE( NOUNIT, FMT = 9998 )M, N, K, IOLDSD, JTYPE, JR, RESULT( JR )
               END IF
            ENDDO
!
            ENDIF
            ENDDO
         ENDIF
         ENDDO
      ENDDO
!
!     Summary
!
   CALL SLASUM( 'SBB', NOUNIT, NERRS, NTESTT )
   RETURN
!
 9999 FORMAT( ' SCHKBB: ', A, ' returned INFO=', I5, '.', / 9X, 'M=', &
         I5, ' N=', I5, ' K=', I5, ', JTYPE=', I5, ', ISEED=(', 3( I5, ',' ), I5, ')' )
 9998 FORMAT( ' M =', I4, ' N=', I4, ', K=', I3, ', seed=', &
         4( I4, ',' ), ' type ', I2, ', test(', I2, ')=', G10.3 )
!
!     End of SCHKBB
!
END




