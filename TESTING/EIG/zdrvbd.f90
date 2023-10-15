!> \brief \b ZDRVBD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZDRVBD( NSIZES, MM, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                          A, LDA, U, LDU, VT, LDVT, ASAV, USAV, VTSAV, S,
!                          SSAV, E, WORK, LWORK, RWORK, IWORK, NOUNIT,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDU, LDVT, LWORK, NOUNIT, NSIZES,
!      $                   NTYPES
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), IWORK( * ), MM( * ), NN( * )
!       DOUBLE PRECISION   E( * ), RWORK( * ), S( * ), SSAV( * )
!       COMPLEX*16         A( LDA, * ), ASAV( LDA, * ), U( LDU, * ),
!      $                   USAV( LDU, * ), VT( LDVT, * ),
!      $                   VTSAV( LDVT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZDRVBD checks the singular value decomposition (SVD) driver ZGESVD,
!> ZGESDD, ZGESVJ, ZGEJSV, ZGESVDX, and ZGESVDQ.
!>
!> ZGESVD and ZGESDD factors A = U diag(S) VT, where U and VT are
!> unitary and diag(S) is diagonal with the entries of the array S on
!> its diagonal. The entries of S are the singular values, nonnegative
!> and stored in decreasing order.  U and VT can be optionally not
!> computed, overwritten on A, or computed partially.
!>
!> A is M by N. Let MNMIN = min( M, N ). S has dimension MNMIN.
!> U can be M by M or M by MNMIN. VT can be N by N or MNMIN by N.
!>
!> When ZDRVBD is called, a number of matrix "sizes" (M's and N's)
!> and a number of matrix "types" are specified.  For each size (M,N)
!> and each type of matrix, and for the minimal workspace as well as
!> workspace adequate to permit blocking, an  M x N  matrix "A" will be
!> generated and used to test the SVD routines.  For each matrix, A will
!> be factored as A = U diag(S) VT and the following 12 tests computed:
!>
!> Test for ZGESVD:
!>
!> (1)   | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (2)   | I - U'U | / ( M ulp )
!>
!> (3)   | I - VT VT' | / ( N ulp )
!>
!> (4)   S contains MNMIN nonnegative values in decreasing order.
!>       (Return 0 if true, 1/ULP if false.)
!>
!> (5)   | U - Upartial | / ( M ulp ) where Upartial is a partially
!>       computed U.
!>
!> (6)   | VT - VTpartial | / ( N ulp ) where VTpartial is a partially
!>       computed VT.
!>
!> (7)   | S - Spartial | / ( MNMIN ulp |S| ) where Spartial is the
!>       vector of singular values from the partial SVD
!>
!> Test for ZGESDD:
!>
!> (8)   | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (9)   | I - U'U | / ( M ulp )
!>
!> (10)  | I - VT VT' | / ( N ulp )
!>
!> (11)  S contains MNMIN nonnegative values in decreasing order.
!>       (Return 0 if true, 1/ULP if false.)
!>
!> (12)  | U - Upartial | / ( M ulp ) where Upartial is a partially
!>       computed U.
!>
!> (13)  | VT - VTpartial | / ( N ulp ) where VTpartial is a partially
!>       computed VT.
!>
!> (14)  | S - Spartial | / ( MNMIN ulp |S| ) where Spartial is the
!>       vector of singular values from the partial SVD
!>
!> Test for ZGESVDQ:
!>
!> (36)  | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (37)  | I - U'U | / ( M ulp )
!>
!> (38)  | I - VT VT' | / ( N ulp )
!>
!> (39)  S contains MNMIN nonnegative values in decreasing order.
!>       (Return 0 if true, 1/ULP if false.)
!>
!> Test for ZGESVJ:
!>
!> (15)  | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (16)  | I - U'U | / ( M ulp )
!>
!> (17)  | I - VT VT' | / ( N ulp )
!>
!> (18)  S contains MNMIN nonnegative values in decreasing order.
!>       (Return 0 if true, 1/ULP if false.)
!>
!> Test for ZGEJSV:
!>
!> (19)  | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (20)  | I - U'U | / ( M ulp )
!>
!> (21)  | I - VT VT' | / ( N ulp )
!>
!> (22)  S contains MNMIN nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!> Test for ZGESVDX( 'V', 'V', 'A' )/ZGESVDX( 'N', 'N', 'A' )
!>
!> (23)  | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (24)  | I - U'U | / ( M ulp )
!>
!> (25)  | I - VT VT' | / ( N ulp )
!>
!> (26)  S contains MNMIN nonnegative values in decreasing order.
!>       (Return 0 if true, 1/ULP if false.)
!>
!> (27)  | U - Upartial | / ( M ulp ) where Upartial is a partially
!>       computed U.
!>
!> (28)  | VT - VTpartial | / ( N ulp ) where VTpartial is a partially
!>       computed VT.
!>
!> (29)  | S - Spartial | / ( MNMIN ulp |S| ) where Spartial is the
!>       vector of singular values from the partial SVD
!>
!> Test for ZGESVDX( 'V', 'V', 'I' )
!>
!> (30)  | U' A VT''' - diag(S) | / ( |A| max(M,N) ulp )
!>
!> (31)  | I - U'U | / ( M ulp )
!>
!> (32)  | I - VT VT' | / ( N ulp )
!>
!> Test for ZGESVDX( 'V', 'V', 'V' )
!>
!> (33)   | U' A VT''' - diag(S) | / ( |A| max(M,N) ulp )
!>
!> (34)   | I - U'U | / ( M ulp )
!>
!> (35)   | I - VT VT' | / ( N ulp )
!>
!> The "sizes" are specified by the arrays MM(1:NSIZES) and
!> NN(1:NSIZES); the value of each element pair (MM(j),NN(j))
!> specifies one size.  The "types" are specified by a logical array
!> DOTYPE( 1:NTYPES ); if DOTYPE(j) is .TRUE., then matrix type "j"
!> will be generated.
!> Currently, the list of possible types is:
!>
!> (1)  The zero matrix.
!> (2)  The identity matrix.
!> (3)  A matrix of the form  U D V, where U and V are unitary and
!>      D has evenly spaced entries 1, ..., ULP with random signs
!>      on the diagonal.
!> (4)  Same as (3), but multiplied by the underflow-threshold / ULP.
!> (5)  Same as (3), but multiplied by the overflow-threshold * ULP.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NSIZES
!> \verbatim
!>          NSIZES is INTEGER
!>          The number of sizes of matrices to use.  If it is zero,
!>          ZDRVBD does nothing.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] MM
!> \verbatim
!>          MM is INTEGER array, dimension (NSIZES)
!>          An array containing the matrix "heights" to be used.  For
!>          each j=1,...,NSIZES, if MM(j) is zero, then MM(j) and NN(j)
!>          will be ignored.  The MM(j) values must be at least zero.
!> \endverbatim
!>
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER array, dimension (NSIZES)
!>          An array containing the matrix "widths" to be used.  For
!>          each j=1,...,NSIZES, if NN(j) is zero, then MM(j) and NN(j)
!>          will be ignored.  The NN(j) values must be at least zero.
!> \endverbatim
!>
!> \param[in] NTYPES
!> \verbatim
!>          NTYPES is INTEGER
!>          The number of elements in DOTYPE.   If it is zero, ZDRVBD
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
!>          next call to ZDRVBD to continue the same random number
!>          sequence.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is DOUBLE PRECISION
!>          A test will count as "failed" if the "error", computed as
!>          described above, exceeds THRESH.  Note that the error
!>          is scaled to be O(1), so THRESH should be a reasonably
!>          small multiple of 1, e.g., 10 or 100.  In particular,
!>          it should not depend on the precision (single vs. double)
!>          or the size of the matrix.  It must be at least zero.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,max(NN))
!>          Used to hold the matrix whose singular values are to be
!>          computed.  On exit, A contains the last matrix actually
!>          used.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  It must be at
!>          least 1 and at least max( MM ).
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX*16 array, dimension (LDU,max(MM))
!>          Used to hold the computed matrix of right singular vectors.
!>          On exit, U contains the last such vectors actually computed.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U.  It must be at
!>          least 1 and at least max( MM ).
!> \endverbatim
!>
!> \param[out] VT
!> \verbatim
!>          VT is COMPLEX*16 array, dimension (LDVT,max(NN))
!>          Used to hold the computed matrix of left singular vectors.
!>          On exit, VT contains the last such vectors actually computed.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>          The leading dimension of VT.  It must be at
!>          least 1 and at least max( NN ).
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is COMPLEX*16 array, dimension (LDA,max(NN))
!>          Used to hold a different copy of the matrix whose singular
!>          values are to be computed.  On exit, A contains the last
!>          matrix actually used.
!> \endverbatim
!>
!> \param[out] USAV
!> \verbatim
!>          USAV is COMPLEX*16 array, dimension (LDU,max(MM))
!>          Used to hold a different copy of the computed matrix of
!>          right singular vectors. On exit, USAV contains the last such
!>          vectors actually computed.
!> \endverbatim
!>
!> \param[out] VTSAV
!> \verbatim
!>          VTSAV is COMPLEX*16 array, dimension (LDVT,max(NN))
!>          Used to hold a different copy of the computed matrix of
!>          left singular vectors. On exit, VTSAV contains the last such
!>          vectors actually computed.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (max(min(MM,NN)))
!>          Contains the computed singular values.
!> \endverbatim
!>
!> \param[out] SSAV
!> \verbatim
!>          SSAV is DOUBLE PRECISION array, dimension (max(min(MM,NN)))
!>          Contains another copy of the computed singular values.
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (max(min(MM,NN)))
!>          Workspace for ZGESVD.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The number of entries in WORK.  This must be at least
!>          MAX(3*MIN(M,N)+MAX(M,N)**2,5*MIN(M,N),3*MAX(M,N)) for all
!>          pairs  (M,N)=(MM(j),NN(j))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array,
!>                      dimension ( 5*max(max(MM,NN)) )
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension at least 8*min(M,N)
!> \endverbatim
!>
!> \param[in] NOUNIT
!> \verbatim
!>          NOUNIT is INTEGER
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
!>           -7: THRESH < 0
!>          -10: LDA < 1 or LDA < MMAX, where MMAX is max( MM(j) ).
!>          -12: LDU < 1 or LDU < MMAX.
!>          -14: LDVT < 1 or LDVT < NMAX, where NMAX is max( NN(j) ).
!>          -21: LWORK too small.
!>          If  ZLATMS, or ZGESVD returns an error code, the
!>              absolute value of it is returned.
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
!> \ingroup complex16_eig
!
!  =====================================================================
   SUBROUTINE ZDRVBD( NSIZES, MM, NN, NTYPES, DOTYPE, ISEED, THRESH, &
                      A, LDA, U, LDU, VT, LDVT, ASAV, USAV, VTSAV, S, &
                      SSAV, E, WORK, LWORK, RWORK, IWORK, NOUNIT, &
                      INFO )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
   IMPLICIT NONE
!
!     .. Scalar Arguments ..
   INTEGER            INFO, LDA, LDU, LDVT, LWORK, NOUNIT, NSIZES, &
                      NTYPES
   DOUBLE PRECISION   THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            ISEED( 4 ), IWORK( * ), MM( * ), NN( * )
   DOUBLE PRECISION   E( * ), RWORK( * ), S( * ), SSAV( * )
   COMPLEX*16         A( LDA, * ), ASAV( LDA, * ), U( LDU, * ), &
                      USAV( LDU, * ), VT( LDVT, * ), &
                      VTSAV( LDVT, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            MAXTYP
   PARAMETER          ( MAXTYP = 5 )
!     ..
!     .. Local Scalars ..
   LOGICAL            BADMM, BADNN
   CHARACTER          JOBQ, JOBU, JOBVT, RANGE
   INTEGER            I, IINFO, IJQ, IJU, IJVT, IL, IU, ITEMP, &
                      IWSPC, IWTMP, J, JSIZE, JTYPE, LSWORK, M, &
                      MINWRK, MMAX, MNMAX, MNMIN, MTYPES, N, &
                      NERRS, NFAIL, NMAX, NS, NSI, NSV, NTEST, &
                      NTESTF, NTESTT, LRWORK
   DOUBLE PRECISION   ANORM, DIF, DIV, OVFL, RTUNFL, ULP, ULPINV, &
                      UNFL, VL, VU
!     ..
!     .. Local Scalars for ZGESVDQ ..
   INTEGER            LIWORK, NUMRANK
!     ..
!     .. Local Arrays ..
   CHARACTER          CJOB( 4 ), CJOBR( 3 ), CJOBV( 2 )
   INTEGER            IOLDSD( 4 ), ISEED2( 4 )
   DOUBLE PRECISION   RESULT( 39 )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH, DLARND
   EXTERNAL           DLAMCH, DLARND
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALASVM, XERBLA, ZBDT01, ZBDT05, ZGESDD, &
                      ZGESVD, ZGESVDQ, ZGESVJ, ZGEJSV, ZGESVDX, &
                      ZLACPY, ZLASET, ZLATMS, ZUNT01, ZUNT03
!     ..
!     .. Scalars in Common ..
   CHARACTER*32       SRNAMT
!     ..
!     .. Common blocks ..
   COMMON             / SRNAMC / SRNAMT
!     ..
!     .. Data statements ..
   DATA               CJOB / 'N', 'O', 'S', 'A' /
   DATA               CJOBR / 'A', 'V', 'I' /
   DATA               CJOBV / 'N', 'V' /
!     ..
!     .. Executable Statements ..
!
!     Check for errors
!
   INFO = 0
!
!     Important constants
!
   NERRS = 0
   NTESTT = 0
   NTESTF = 0
   NMAX = MAXVAL(NN(1:NSIZES))
   BADNN = any(NN(1:NSIZES) < 0 )
   MMAX = MAXVAL(MM(1:NSIZES))
   BADMM = any(MM(1:NSIZES) < 0 )
   MNMAX = 1
   MINWRK = 1
   DO J = 1, NSIZES
      MNMAX = MAX( MNMAX, MIN( MM( J ), NN( J ) ) )
      MINWRK = MAX( MINWRK, MAX( 3*MIN( MM( J ), &
               NN( J ) )+MAX( MM( J ), NN( J ) )**2, 5*MIN( MM( J ), &
               NN( J ) ), 3*MAX( MM( J ), NN( J ) ) ) )
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
   ELSE IF( LDA < MAX( 1, MMAX ) ) THEN
      INFO = -10
   ELSE IF( LDU < MAX( 1, MMAX ) ) THEN
      INFO = -12
   ELSE IF( LDVT < MAX( 1, NMAX ) ) THEN
      INFO = -14
   ELSE IF( MINWRK > LWORK ) THEN
      INFO = -21
   END IF
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'ZDRVBD', -INFO )
      RETURN
   END IF
!
!     Quick return if nothing to do
!
   IF( NSIZES == 0 .OR. NTYPES == 0 ) RETURN
!
!     More Important constants
!
   UNFL = DLAMCH( 'S' )
   OVFL = 1.0D0 / UNFL
   ULP = DLAMCH( 'E' )
   ULPINV = 1.0D0 / ULP
   RTUNFL = SQRT( UNFL )
!
!     Loop over sizes, types
!
   NERRS = 0
!
   DO JSIZE = 1, NSIZES
      M = MM( JSIZE )
      N = NN( JSIZE )
      MNMIN = MIN( M, N )
!
      IF( NSIZES /= 1 ) THEN
         MTYPES = MIN( MAXTYP, NTYPES )
      ELSE
         MTYPES = MIN( MAXTYP+1, NTYPES )
      END IF
!
      DO JTYPE = 1, MTYPES
         IF( .NOT.DOTYPE( JTYPE ) ) GO TO 220
         NTEST = 0
!
         IOLDSD(1:4) = ISEED(1:4)
!
!           Compute "A"
!
         IF( MTYPES > MAXTYP ) GO TO 50
!
         IF( JTYPE == 1 ) THEN
!
!              Zero matrix
!
            CALL ZLASET( 'Full', M, N, (0.0D0,0.0D0), (0.0D0,0.0D0), A, LDA )
            S(1:MIN( M, N )) = 0.0D+0
!
         ELSE IF( JTYPE == 2 ) THEN
!
!              Identity matrix
!
            CALL ZLASET( 'Full', M, N, (0.0D0,0.0D0), (1.0D0,0.0D0), A, LDA )
            S(1:MIN( M, N )) = 1.0D+0
!
         ELSE
!
!              (Scaled) random matrix
!
            IF( JTYPE == 3 ) ANORM = 1.0D0
            IF( JTYPE == 4 ) ANORM = UNFL / ULP
            IF( JTYPE == 5 ) ANORM = OVFL*ULP
            CALL ZLATMS( M, N, 'U', ISEED, 'N', S, 4, DBLE( MNMIN ), &
                         ANORM, M-1, N-1, 'N', A, LDA, WORK, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9996 )'Generator', IINFO, M, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
         END IF
!
50       CONTINUE
         CALL ZLACPY( 'F', M, N, A, LDA, ASAV, LDA )
!
!           Do for minimal and adequate (for blocking) workspace
!
         DO IWSPC = 1, 4
!
!              Test for ZGESVD
!
            IWTMP = 2*MIN( M, N )+MAX( M, N )
            LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
            LSWORK = MIN( LSWORK, LWORK )
            LSWORK = MAX( LSWORK, 1 )
            IF( IWSPC == 4 ) LSWORK = LWORK
!
            RESULT(1:35) = -1.0D+0
!
!              Factorize A
!
            IF( IWSPC > 1 ) CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
            SRNAMT = 'ZGESVD'
            CALL ZGESVD( 'A', 'A', M, N, A, LDA, SSAV, USAV, LDU, &
                         VTSAV, LDVT, WORK, LSWORK, RWORK, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9995 )'GESVD', IINFO, M, N, &
                  JTYPE, LSWORK, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
!
!              Do tests 1--4
!
            CALL ZBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, &
                         VTSAV, LDVT, WORK, RWORK, RESULT( 1 ) )
            IF( M /= 0 .AND. N /= 0 ) THEN
               CALL ZUNT01( 'Columns', MNMIN, M, USAV, LDU, WORK, &
                            LWORK, RWORK, RESULT( 2 ) )
               CALL ZUNT01( 'Rows', MNMIN, N, VTSAV, LDVT, WORK, &
                            LWORK, RWORK, RESULT( 3 ) )
            END IF
            RESULT( 4 ) = 0
            DO I = 1, MNMIN - 1
               IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 4 ) = ULPINV
               IF( SSAV( I ) < 0.0D0 ) RESULT( 4 ) = ULPINV
            ENDDO
            IF( MNMIN >= 1 ) THEN
               IF( SSAV( MNMIN ) < 0.0D0 ) RESULT( 4 ) = ULPINV
            END IF
!
!              Do partial SVDs, comparing to SSAV, USAV, and VTSAV
!
            RESULT( 5:7 ) = 0.0D+0
            DO IJU = 0, 3
               DO IJVT = 0, 3
                  IF( ( IJU == 3 .AND. IJVT == 3 ) .OR. &
                      ( IJU == 1 .AND. IJVT == 1 ) )GO TO 90
                  JOBU = CJOB( IJU+1 )
                  JOBVT = CJOB( IJVT+1 )
                  CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
                  SRNAMT = 'ZGESVD'
                  CALL ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, &
                               VT, LDVT, WORK, LSWORK, RWORK, IINFO )
!
!                    Compare U
!
                  DIF = 0.0D0
                  IF( M > 0 .AND. N > 0 ) THEN
                     IF( IJU == 1 ) THEN
                        CALL ZUNT03( 'C', M, MNMIN, M, MNMIN, USAV, &
                                     LDU, A, LDA, WORK, LWORK, RWORK, &
                                     DIF, IINFO )
                     ELSE IF( IJU == 2 ) THEN
                        CALL ZUNT03( 'C', M, MNMIN, M, MNMIN, USAV, &
                                     LDU, U, LDU, WORK, LWORK, RWORK, &
                                     DIF, IINFO )
                     ELSE IF( IJU == 3 ) THEN
                        CALL ZUNT03( 'C', M, M, M, MNMIN, USAV, LDU, &
                                     U, LDU, WORK, LWORK, RWORK, DIF, &
                                     IINFO )
                     END IF
                  END IF
                  RESULT( 5 ) = MAX( RESULT( 5 ), DIF )
!
!                    Compare VT
!
                  DIF = 0.0D0
                  IF( M > 0 .AND. N > 0 ) THEN
                     IF( IJVT == 1 ) THEN
                        CALL ZUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV, &
                                     LDVT, A, LDA, WORK, LWORK, &
                                     RWORK, DIF, IINFO )
                     ELSE IF( IJVT == 2 ) THEN
                        CALL ZUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV, &
                                     LDVT, VT, LDVT, WORK, LWORK, &
                                     RWORK, DIF, IINFO )
                     ELSE IF( IJVT == 3 ) THEN
                        CALL ZUNT03( 'R', N, N, N, MNMIN, VTSAV, &
                                     LDVT, VT, LDVT, WORK, LWORK, &
                                     RWORK, DIF, IINFO )
                     END IF
                  END IF
                  RESULT( 6 ) = MAX( RESULT( 6 ), DIF )
!
!                    Compare S
!
                  DIF = 0.0D0
                  DIV = MAX( DBLE( MNMIN )*ULP*S( 1 ), &
                        DLAMCH( 'Safe minimum' ) )
                  DO I = 1, MNMIN - 1
                     IF( SSAV( I ) < SSAV( I+1 ) ) DIF = ULPINV
                     IF( SSAV( I ) < 0.0D0 ) DIF = ULPINV
                     DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
                  ENDDO
                  RESULT( 7 ) = MAX( RESULT( 7 ), DIF )
90             CONTINUE
               ENDDO
               ENDDO
!
!              Test for ZGESDD
!
            IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + MAX( M, N )
            LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
            LSWORK = MIN( LSWORK, LWORK )
            LSWORK = MAX( LSWORK, 1 )
            IF( IWSPC == 4 ) LSWORK = LWORK
!
!              Factorize A
!
            CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
            SRNAMT = 'ZGESDD'
            CALL ZGESDD( 'A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, &
                         LDVT, WORK, LSWORK, RWORK, IWORK, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9995 )'GESDD', IINFO, M, N, &
                  JTYPE, LSWORK, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
!
!              Do tests 1--4
!
            CALL ZBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, &
                         VTSAV, LDVT, WORK, RWORK, RESULT( 8 ) )
            IF( M /= 0 .AND. N /= 0 ) THEN
               CALL ZUNT01( 'Columns', MNMIN, M, USAV, LDU, WORK, &
                            LWORK, RWORK, RESULT( 9 ) )
               CALL ZUNT01( 'Rows', MNMIN, N, VTSAV, LDVT, WORK, &
                            LWORK, RWORK, RESULT( 10 ) )
            END IF
            RESULT( 11 ) = 0
            DO I = 1, MNMIN - 1
               IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 11 ) = ULPINV
               IF( SSAV( I ) < 0.0D0 ) RESULT( 11 ) = ULPINV
               ENDDO
            IF( MNMIN >= 1 ) THEN
               IF( SSAV( MNMIN ) < 0.0D0 ) RESULT( 11 ) = ULPINV
            END IF
!
!              Do partial SVDs, comparing to SSAV, USAV, and VTSAV
!
            RESULT( 12:14 ) = 0.0D+0
            DO IJQ = 0, 2
               JOBQ = CJOB( IJQ+1 )
               CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               SRNAMT = 'ZGESDD'
               CALL ZGESDD( JOBQ, M, N, A, LDA, S, U, LDU, VT, LDVT, &
                            WORK, LSWORK, RWORK, IWORK, IINFO )
!
!                 Compare U
!
               DIF = 0.0D0
               IF( M > 0 .AND. N > 0 ) THEN
                  IF( IJQ == 1 ) THEN
                     IF( M >= N ) THEN
                        CALL ZUNT03( 'C', M, MNMIN, M, MNMIN, USAV, &
                                     LDU, A, LDA, WORK, LWORK, RWORK, &
                                     DIF, IINFO )
                     ELSE
                        CALL ZUNT03( 'C', M, MNMIN, M, MNMIN, USAV, &
                                     LDU, U, LDU, WORK, LWORK, RWORK, &
                                     DIF, IINFO )
                     END IF
                  ELSE IF( IJQ == 2 ) THEN
                     CALL ZUNT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, &
                                  U, LDU, WORK, LWORK, RWORK, DIF, &
                                  IINFO )
                  END IF
               END IF
               RESULT( 12 ) = MAX( RESULT( 12 ), DIF )
!
!                 Compare VT
!
               DIF = 0.0D0
               IF( M > 0 .AND. N > 0 ) THEN
                  IF( IJQ == 1 ) THEN
                     IF( M >= N ) THEN
                        CALL ZUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV, &
                                     LDVT, VT, LDVT, WORK, LWORK, &
                                     RWORK, DIF, IINFO )
                     ELSE
                        CALL ZUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV, &
                                     LDVT, A, LDA, WORK, LWORK, &
                                     RWORK, DIF, IINFO )
                     END IF
                  ELSE IF( IJQ == 2 ) THEN
                     CALL ZUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV, &
                                  LDVT, VT, LDVT, WORK, LWORK, RWORK, &
                                  DIF, IINFO )
                  END IF
               END IF
               RESULT( 13 ) = MAX( RESULT( 13 ), DIF )
!
!                 Compare S
!
               DIF = 0.0D0
               DIV = MAX( DBLE( MNMIN )*ULP*S( 1 ), DLAMCH( 'Safe minimum' ) )
               DO I = 1, MNMIN - 1
                  IF( SSAV( I ) < SSAV( I+1 ) ) DIF = ULPINV
                  IF( SSAV( I ) < 0.0D0 ) DIF = ULPINV
                  DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
                  ENDDO
               RESULT( 14 ) = MAX( RESULT( 14 ), DIF )
               ENDDO
!
!              Test ZGESVDQ
!              Note: ZGESVDQ only works for M >= N
!
            RESULT( 36:39 ) = 0.0D+0
!
            IF( M >= N ) THEN
               IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + MAX( M, N )
               LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               IF( IWSPC == 4 ) LSWORK = LWORK
!
               CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               SRNAMT = 'ZGESVDQ'
!
               LRWORK = MAX(2, M, 5*N)
               LIWORK = MAX( N, 1 )
               CALL ZGESVDQ( 'H', 'N', 'N', 'A', 'A', &
                             M, N, A, LDA, SSAV, USAV, LDU, &
                             VTSAV, LDVT, NUMRANK, IWORK, LIWORK, &
                             WORK, LWORK, RWORK, LRWORK, IINFO )
!
               IF( IINFO /= 0 ) THEN
                  WRITE( NOUNIT, FMT = 9995 )'ZGESVDQ', IINFO, M, N, &
                  JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
!
!                 Do tests 36--39
!
               CALL ZBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, &
                            VTSAV, LDVT, WORK, RWORK, RESULT( 36 ) )
               IF( M /= 0 .AND. N /= 0 ) THEN
                  CALL ZUNT01( 'Columns', M, M, USAV, LDU, WORK, &
                               LWORK, RWORK, RESULT( 37 ) )
                  CALL ZUNT01( 'Rows', N, N, VTSAV, LDVT, WORK, &
                               LWORK, RWORK, RESULT( 38 ) )
               END IF
               RESULT( 39 ) = 0.0D0
               DO I = 1, MNMIN - 1
                  IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 39 ) = ULPINV
                  IF( SSAV( I ) < 0.0D0 ) RESULT( 39 ) = ULPINV
                  ENDDO
               IF( MNMIN >= 1 ) THEN
                  IF( SSAV( MNMIN ) < 0.0D0 ) RESULT( 39 ) = ULPINV
               END IF
            END IF
!
!              Test ZGESVJ
!              Note: ZGESVJ only works for M >= N
!
            RESULT( 15:18 ) = 0.0D+0
!
            IF( M >= N ) THEN
               IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + MAX( M, N )
               LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               LRWORK = MAX(6,N)
               IF( IWSPC == 4 ) LSWORK = LWORK
!
               CALL ZLACPY( 'F', M, N, ASAV, LDA, USAV, LDA )
               SRNAMT = 'ZGESVJ'
               CALL ZGESVJ( 'G', 'U', 'V', M, N, USAV, LDA, SSAV, &
                           0, A, LDVT, WORK, LWORK, RWORK, &
                           LRWORK, IINFO )
!
!                 ZGESVJ returns V not VH
!
               DO J=1,N
                  VTSAV(J,1:N) = CONJG(A(1:N,J))
               END DO
!
               IF( IINFO /= 0 ) THEN
                  WRITE( NOUNIT, FMT = 9995 )'GESVJ', IINFO, M, N, &
                  JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
!
!                 Do tests 15--18
!
               CALL ZBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, &
                            VTSAV, LDVT, WORK, RWORK, RESULT( 15 ) )
               IF( M /= 0 .AND. N /= 0 ) THEN
                  CALL ZUNT01( 'Columns', M, M, USAV, LDU, WORK, &
                               LWORK, RWORK, RESULT( 16 ) )
                  CALL ZUNT01( 'Rows', N, N, VTSAV, LDVT, WORK, &
                               LWORK, RWORK, RESULT( 17 ) )
               END IF
               RESULT( 18 ) = 0.0D0
               DO I = 1, MNMIN - 1
                  IF( SSAV( I ) < SSAV( I+1 ) ) &
                     RESULT( 18 ) = ULPINV
                  IF( SSAV( I ) < 0.0D0 ) &
                     RESULT( 18 ) = ULPINV
                  ENDDO
               IF( MNMIN >= 1 ) THEN
                  IF( SSAV( MNMIN ) < 0.0D0 ) &
                     RESULT( 18 ) = ULPINV
               END IF
            END IF
!
!              Test ZGEJSV
!              Note: ZGEJSV only works for M >= N
!
            RESULT( 19:22 ) = 0.0D+0
            IF( M >= N ) THEN
               IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + MAX( M, N )
               LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               IF( IWSPC == 4 ) LSWORK = LWORK
               LRWORK = MAX( 7, N + 2*M)
!
               CALL ZLACPY( 'F', M, N, ASAV, LDA, VTSAV, LDA )
               SRNAMT = 'ZGEJSV'
               CALL ZGEJSV( 'G', 'U', 'V', 'R', 'N', 'N', &
                      M, N, VTSAV, LDA, SSAV, USAV, LDU, A, LDVT, &
                      WORK, LWORK, RWORK, &
                      LRWORK, IWORK, IINFO )
!
!                 ZGEJSV returns V not VH
!
               DO J=1,N
                  VTSAV(J,1:N) = CONJG (A(1:N,J))
               END DO
!
               IF( IINFO /= 0 ) THEN
                  WRITE( NOUNIT, FMT = 9995 )'GEJSV', IINFO, M, N, &
                  JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
!
!                 Do tests 19--22
!
               CALL ZBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, &
                            VTSAV, LDVT, WORK, RWORK, RESULT( 19 ) )
               IF( M /= 0 .AND. N /= 0 ) THEN
                  CALL ZUNT01( 'Columns', M, M, USAV, LDU, WORK, &
                               LWORK, RWORK, RESULT( 20 ) )
                  CALL ZUNT01( 'Rows', N, N, VTSAV, LDVT, WORK, &
                               LWORK, RWORK, RESULT( 21 ) )
               END IF
               RESULT( 22 ) = 0.0D0
               DO I = 1, MNMIN - 1
                  IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 22 ) = ULPINV
                  IF( SSAV( I ) < 0.0D0 ) RESULT( 22 ) = ULPINV
                  ENDDO
               IF( MNMIN >= 1 ) THEN
                  IF( SSAV( MNMIN ) < 0.0D0 ) RESULT( 22 ) = ULPINV
               END IF
            END IF
!
!              Test ZGESVDX
!
!              Factorize A
!
            CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
            SRNAMT = 'ZGESVDX'
            CALL ZGESVDX( 'V', 'V', 'A', M, N, A, LDA, &
                          VL, VU, IL, IU, NS, SSAV, USAV, LDU, &
                          VTSAV, LDVT, WORK, LWORK, RWORK, &
                          IWORK, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9995 )'GESVDX', IINFO, M, N, &
                  JTYPE, LSWORK, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
!
!              Do tests 1--4
!
            RESULT( 23:25 ) = 0.0D+0
            CALL ZBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, &
                         VTSAV, LDVT, WORK, RWORK, RESULT( 23 ) )
            IF( M /= 0 .AND. N /= 0 ) THEN
               CALL ZUNT01( 'Columns', MNMIN, M, USAV, LDU, WORK, &
                            LWORK, RWORK, RESULT( 24 ) )
               CALL ZUNT01( 'Rows', MNMIN, N, VTSAV, LDVT, WORK, &
                            LWORK, RWORK, RESULT( 25 ) )
            END IF
            RESULT( 26 ) = 0.0D0
            DO I = 1, MNMIN - 1
               IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 26 ) = ULPINV
               IF( SSAV( I ) < 0.0D0 ) RESULT( 26 ) = ULPINV
               ENDDO
            IF( MNMIN >= 1 ) THEN
               IF( SSAV( MNMIN ) < 0.0D0 ) RESULT( 26 ) = ULPINV
            END IF
!
!              Do partial SVDs, comparing to SSAV, USAV, and VTSAV
!
            RESULT( 27:29 ) = 0.0D+0
            DO IJU = 0, 1
               DO IJVT = 0, 1
                  IF( ( IJU == 0 .AND. IJVT == 0 ) .OR. &
                      ( IJU == 1 .AND. IJVT == 1 ) ) GO TO 160
                  JOBU = CJOBV( IJU+1 )
                  JOBVT = CJOBV( IJVT+1 )
                  RANGE = CJOBR( 1 )
                  CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
                  SRNAMT = 'ZGESVDX'
                  CALL ZGESVDX( JOBU, JOBVT, 'A', M, N, A, LDA, &
                               VL, VU, IL, IU, NS, SSAV, U, LDU, &
                               VT, LDVT, WORK, LWORK, RWORK, &
                               IWORK, IINFO )
!
!                    Compare U
!
                  DIF = 0.0D0
                  IF( M > 0 .AND. N > 0 ) THEN
                     IF( IJU == 1 ) THEN
                        CALL ZUNT03( 'C', M, MNMIN, M, MNMIN, USAV, &
                                     LDU, U, LDU, WORK, LWORK, RWORK, &
                                     DIF, IINFO )
                     END IF
                  END IF
                  RESULT( 27 ) = MAX( RESULT( 27 ), DIF )
!
!                    Compare VT
!
                  DIF = 0.0D0
                  IF( M > 0 .AND. N > 0 ) THEN
                     IF( IJVT == 1 ) THEN
                        CALL ZUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV, &
                                     LDVT, VT, LDVT, WORK, LWORK, &
                                     RWORK, DIF, IINFO )
                     END IF
                  END IF
                  RESULT( 28 ) = MAX( RESULT( 28 ), DIF )
!
!                    Compare S
!
                  DIF = 0.0D0
                  DIV = MAX( DBLE( MNMIN )*ULP*S( 1 ), DLAMCH( 'Safe minimum' ) )
                  DO I = 1, MNMIN - 1
                     IF( SSAV( I ) < SSAV( I+1 ) ) DIF = ULPINV
                     IF( SSAV( I ) < 0.0D0 ) DIF = ULPINV
                     DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
                  ENDDO
                  RESULT( 29) = MAX( RESULT( 29 ), DIF )
  160             CONTINUE
                  ENDDO
               ENDDO
!
!              Do tests 8--10
!
            ISEED2(1:4) = ISEED(1:4)
            IF( MNMIN <= 1 ) THEN
               IL = 1
               IU = MAX( 1, MNMIN )
            ELSE
               IL = 1 + INT( ( MNMIN-1 )*DLARND( 1, ISEED2 ) )
               IU = 1 + INT( ( MNMIN-1 )*DLARND( 1, ISEED2 ) )
               IF( IU < IL ) THEN
                  ITEMP = IU
                  IU = IL
                  IL = ITEMP
               END IF
            END IF
            CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
            SRNAMT = 'ZGESVDX'
            CALL ZGESVDX( 'V', 'V', 'I', M, N, A, LDA, &
                          VL, VU, IL, IU, NSI, S, U, LDU, &
                          VT, LDVT, WORK, LWORK, RWORK, &
                          IWORK, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9995 )'GESVDX', IINFO, M, N, &
                  JTYPE, LSWORK, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
!
            RESULT( 30:32 ) = 0.0D+0
            CALL ZBDT05( M, N, ASAV, LDA, S, NSI, U, LDU, &
                         VT, LDVT, WORK, RESULT( 30 ) )
            IF( M /= 0 .AND. N /= 0 ) THEN
               CALL ZUNT01( 'Columns', M, NSI, U, LDU, WORK, LWORK, RWORK, RESULT( 31 ) )
               CALL ZUNT01( 'Rows', NSI, N, VT, LDVT, WORK, LWORK, RWORK, RESULT( 32 ) )
            END IF
!
!              Do tests 11--13
!
            IF( MNMIN > 0 .AND. NSI > 1 ) THEN
               IF( IL /= 1 ) THEN
                  VU = SSAV( IL ) + &
                       MAX( 0.5D+0*ABS( SSAV( IL )-SSAV( IL-1 ) ), ULP*ANORM, 2.0D+0*RTUNFL )
               ELSE
                  VU = SSAV( 1 ) + &
                       MAX( 0.5D+0*ABS( SSAV( NS )-SSAV( 1 ) ), ULP*ANORM, 2.0D+0*RTUNFL )
               END IF
               IF( IU /= NS ) THEN
                  VL = SSAV( IU ) - MAX( ULP*ANORM, 2.0D+0*RTUNFL, 0.5D+0*ABS( SSAV( IU+1 )-SSAV( IU ) ) )
               ELSE
                  VL = SSAV( NS ) - MAX( ULP*ANORM, 2.0D+0*RTUNFL, 0.5D+0*ABS( SSAV( NS )-SSAV( 1 ) ) )
               END IF
               VL = MAX( VL,0.0D0 )
               VU = MAX( VU,0.0D0 )
               IF( VL >= VU ) VU = MAX( VU*2, VU+VL+0.5D+0 )
            ELSE
               VL = 0.0D0
               VU = 1.0D0
            END IF
            CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
            SRNAMT = 'ZGESVDX'
            CALL ZGESVDX( 'V', 'V', 'V', M, N, A, LDA, &
                          VL, VU, IL, IU, NSV, S, U, LDU, &
                          VT, LDVT, WORK, LWORK, RWORK, &
                          IWORK, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9995 )'GESVDX', IINFO, M, N, &
                  JTYPE, LSWORK, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
!
            RESULT( 33:35 ) = 0.0D+0
            CALL ZBDT05( M, N, ASAV, LDA, S, NSV, U, LDU, &
                         VT, LDVT, WORK, RESULT( 33 ) )
            IF( M /= 0 .AND. N /= 0 ) THEN
               CALL ZUNT01( 'Columns', M, NSV, U, LDU, WORK, &
                            LWORK, RWORK, RESULT( 34 ) )
               CALL ZUNT01( 'Rows', NSV, N, VT, LDVT, WORK, &
                            LWORK, RWORK, RESULT( 35 ) )
            END IF
!
!              End of Loop -- Check for RESULT(j) > THRESH
!
            NTEST = COUNT(RESULT(1:39) >= 0.0D+0)
            NFAIL = COUNT(RESULT(1:39) >= THRESH)
!
            IF( NFAIL > 0 ) NTESTF = NTESTF + 1
            IF( NTESTF == 1 ) THEN
               WRITE( NOUNIT, FMT = 9999 )
               WRITE( NOUNIT, FMT = 9998 )THRESH
               NTESTF = 2
            END IF
!
            DO J = 1, 39
               IF( RESULT( J ) >= THRESH ) THEN
                  WRITE( NOUNIT, FMT = 9997 )M, N, JTYPE, IWSPC, &
                     IOLDSD, J, RESULT( J )
               END IF
            ENDDO
!
            NERRS = NERRS + NFAIL
            NTESTT = NTESTT + NTEST
!
            ENDDO
!
  220    CONTINUE
         ENDDO
      ENDDO
!
!     Summary
!
   CALL ALASVM( 'ZBD', NOUNIT, NERRS, NTESTT, 0 )
!
 9999 FORMAT( ' SVD -- Complex Singular Value Decomposition Driver ', &
         / ' Matrix types (see ZDRVBD for details):', &
         / / ' 1 = Zero matrix', / ' 2 = Identity matrix', &
         / ' 3 = Evenly spaced singular values near 1', &
         / ' 4 = Evenly spaced singular values near underflow', &
         / ' 5 = Evenly spaced singular values near overflow', &
         / / ' Tests performed: ( A is dense, U and V are unitary,', &
         / 19X, ' S is an array, and Upartial, VTpartial, and', &
         / 19X, ' Spartial are partially computed U, VT and S),', / )
 9998 FORMAT( ' Tests performed with Test Threshold = ', F8.2, &
         / ' ZGESVD: ', / &
         ' 1 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', &
         / ' 2 = | I - U**T U | / ( M ulp ) ', &
         / ' 3 = | I - VT VT**T | / ( N ulp ) ', &
         / ' 4 = 0 if S contains min(M,N) nonnegative values in', &
         ' decreasing order, else 1/ulp', &
         / ' 5 = | U - Upartial | / ( M ulp )', &
         / ' 6 = | VT - VTpartial | / ( N ulp )', &
         / ' 7 = | S - Spartial | / ( min(M,N) ulp |S| )', &
         / ' ZGESDD: ', / &
         ' 8 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', &
         / ' 9 = | I - U**T U | / ( M ulp ) ', &
         / '10 = | I - VT VT**T | / ( N ulp ) ', &
         / '11 = 0 if S contains min(M,N) nonnegative values in', &
         ' decreasing order, else 1/ulp', &
         / '12 = | U - Upartial | / ( M ulp )', &
         / '13 = | VT - VTpartial | / ( N ulp )', &
         / '14 = | S - Spartial | / ( min(M,N) ulp |S| )', &
         / ' ZGESVJ: ', / &
         / '15 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', &
         / '16 = | I - U**T U | / ( M ulp ) ', &
         / '17 = | I - VT VT**T | / ( N ulp ) ', &
         / '18 = 0 if S contains min(M,N) nonnegative values in', &
         ' decreasing order, else 1/ulp', &
         / ' ZGESJV: ', / &
         / '19 = | A - U diag(S) VT | / ( |A| max(M,N) ulp )', &
         / '20 = | I - U**T U | / ( M ulp ) ', &
         / '21 = | I - VT VT**T | / ( N ulp ) ', &
         / '22 = 0 if S contains min(M,N) nonnegative values in', &
         ' decreasing order, else 1/ulp', &
         / ' ZGESVDX(V,V,A): ', / &
           '23 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', &
         / '24 = | I - U**T U | / ( M ulp ) ', &
         / '25 = | I - VT VT**T | / ( N ulp ) ', &
         / '26 = 0 if S contains min(M,N) nonnegative values in', &
         ' decreasing order, else 1/ulp', &
         / '27 = | U - Upartial | / ( M ulp )', &
         / '28 = | VT - VTpartial | / ( N ulp )', &
         / '29 = | S - Spartial | / ( min(M,N) ulp |S| )', &
         / ' ZGESVDX(V,V,I): ', &
         / '30 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp )', &
         / '31 = | I - U**T U | / ( M ulp ) ', &
         / '32 = | I - VT VT**T | / ( N ulp ) ', &
         / ' ZGESVDX(V,V,V) ', &
         / '33 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp )', &
         / '34 = | I - U**T U | / ( M ulp ) ', &
         / '35 = | I - VT VT**T | / ( N ulp ) ', &
         ' ZGESVDQ(H,N,N,A,A', &
         / '36 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', &
         / '37 = | I - U**T U | / ( M ulp ) ', &
         / '38 = | I - VT VT**T | / ( N ulp ) ', &
         / '39 = 0 if S contains min(M,N) nonnegative values in', &
         ' decreasing order, else 1/ulp', &
         / / )
 9997 FORMAT( ' M=', I5, ', N=', I5, ', type ', I1, ', IWS=', I1, &
         ', seed=', 4( I4, ',' ), ' test(', I2, ')=', G11.4 )
 9996 FORMAT( ' ZDRVBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=', &
         I6, ', N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), &
         I5, ')' )
 9995 FORMAT( ' ZDRVBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=', &
         I6, ', N=', I6, ', JTYPE=', I6, ', LSWORK=', I6, / 9X, &
         'ISEED=(', 3( I5, ',' ), I5, ')' )
!
   RETURN
!
!     End of ZDRVBD
!
END
