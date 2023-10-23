!> \brief \b CLAR1V computes the (scaled) r-th column of the inverse of the submatrix in rows b1 through bn of the tridiagonal matrix LDLT - λI.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAR1V + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clar1v.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clar1v.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clar1v.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAR1V( N, B1, BN, LAMBDA, D, L, LD, LLD,
!                  PIVMIN, GAPTOL, Z, WANTNC, NEGCNT, ZTZ, MINGMA,
!                  R, ISUPPZ, NRMINV, RESID, RQCORR, WORK )
!
!       .. Scalar Arguments ..
!       LOGICAL            WANTNC
!       INTEGER   B1, BN, N, NEGCNT, R
!       REAL               GAPTOL, LAMBDA, MINGMA, NRMINV, PIVMIN, RESID,
!      $                   RQCORR, ZTZ
!       ..
!       .. Array Arguments ..
!       INTEGER            ISUPPZ( * )
!       REAL               D( * ), L( * ), LD( * ), LLD( * ),
!      $                  WORK( * )
!       COMPLEX          Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAR1V computes the (scaled) r-th column of the inverse of
!> the sumbmatrix in rows B1 through BN of the tridiagonal matrix
!> L D L**T - sigma I. When sigma is close to an eigenvalue, the
!> computed vector is an accurate eigenvector. Usually, r corresponds
!> to the index where the eigenvector is largest in magnitude.
!> The following steps accomplish this computation :
!> (a) Stationary qd transform,  L D L**T - sigma I = L(+) D(+) L(+)**T,
!> (b) Progressive qd transform, L D L**T - sigma I = U(-) D(-) U(-)**T,
!> (c) Computation of the diagonal elements of the inverse of
!>     L D L**T - sigma I by combining the above transforms, and choosing
!>     r as the index where the diagonal of the inverse is (one of the)
!>     largest in magnitude.
!> (d) Computation of the (scaled) r-th column of the inverse using the
!>     twisted factorization obtained by combining the top part of the
!>     the stationary and the bottom part of the progressive transform.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           The order of the matrix L D L**T.
!> \endverbatim
!>
!> \param[in] B1
!> \verbatim
!>          B1 is INTEGER
!>           First index of the submatrix of L D L**T.
!> \endverbatim
!>
!> \param[in] BN
!> \verbatim
!>          BN is INTEGER
!>           Last index of the submatrix of L D L**T.
!> \endverbatim
!>
!> \param[in] LAMBDA
!> \verbatim
!>          LAMBDA is REAL
!>           The shift. In order to compute an accurate eigenvector,
!>           LAMBDA should be a good approximation to an eigenvalue
!>           of L D L**T.
!> \endverbatim
!>
!> \param[in] L
!> \verbatim
!>          L is REAL array, dimension (N-1)
!>           The (n-1) subdiagonal elements of the unit bidiagonal matrix
!>           L, in elements 1 to N-1.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>           The n diagonal elements of the diagonal matrix D.
!> \endverbatim
!>
!> \param[in] LD
!> \verbatim
!>          LD is REAL array, dimension (N-1)
!>           The n-1 elements L(i)*D(i).
!> \endverbatim
!>
!> \param[in] LLD
!> \verbatim
!>          LLD is REAL array, dimension (N-1)
!>           The n-1 elements L(i)*L(i)*D(i).
!> \endverbatim
!>
!> \param[in] PIVMIN
!> \verbatim
!>          PIVMIN is REAL
!>           The minimum pivot in the Sturm sequence.
!> \endverbatim
!>
!> \param[in] GAPTOL
!> \verbatim
!>          GAPTOL is REAL
!>           Tolerance that indicates when eigenvector entries are negligible
!>           w.r.t. their contribution to the residual.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (N)
!>           On input, all entries of Z must be set to 0.
!>           On output, Z contains the (scaled) r-th column of the
!>           inverse. The scaling is such that Z(R) equals 1.
!> \endverbatim
!>
!> \param[in] WANTNC
!> \verbatim
!>          WANTNC is LOGICAL
!>           Specifies whether NEGCNT has to be computed.
!> \endverbatim
!>
!> \param[out] NEGCNT
!> \verbatim
!>          NEGCNT is INTEGER
!>           If WANTNC is .TRUE. then NEGCNT = the number of pivots < pivmin
!>           in the  matrix factorization L D L**T, and NEGCNT = -1 otherwise.
!> \endverbatim
!>
!> \param[out] ZTZ
!> \verbatim
!>          ZTZ is REAL
!>           The square of the 2-norm of Z.
!> \endverbatim
!>
!> \param[out] MINGMA
!> \verbatim
!>          MINGMA is REAL
!>           The reciprocal of the largest (in magnitude) diagonal
!>           element of the inverse of L D L**T - sigma I.
!> \endverbatim
!>
!> \param[in,out] R
!> \verbatim
!>          R is INTEGER
!>           The twist index for the twisted factorization used to
!>           compute Z.
!>           On input, 0 <= R <= N. If R is input as 0, R is set to
!>           the index where (L D L**T - sigma I)^{-1} is largest
!>           in magnitude. If 1 <= R <= N, R is unchanged.
!>           On output, R contains the twist index used to compute Z.
!>           Ideally, R designates the position of the maximum entry in the
!>           eigenvector.
!> \endverbatim
!>
!> \param[out] ISUPPZ
!> \verbatim
!>          ISUPPZ is INTEGER array, dimension (2)
!>           The support of the vector in Z, i.e., the vector Z is
!>           nonzero only in elements ISUPPZ(1) through ISUPPZ( 2 ).
!> \endverbatim
!>
!> \param[out] NRMINV
!> \verbatim
!>          NRMINV is REAL
!>           NRMINV = 1/SQRT( ZTZ )
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
!>           The residual of the FP vector.
!>           RESID = ABS( MINGMA )/SQRT( ZTZ )
!> \endverbatim
!>
!> \param[out] RQCORR
!> \verbatim
!>          RQCORR is REAL
!>           The Rayleigh Quotient correction to LAMBDA.
!>           RQCORR = MINGMA*TMP
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (4*N)
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
!> \ingroup lar1v
!
!> \par Contributors:
!  ==================
!>
!> Beresford Parlett, University of California, Berkeley, USA \n
!> Jim Demmel, University of California, Berkeley, USA \n
!> Inderjit Dhillon, University of Texas, Austin, USA \n
!> Osni Marques, LBNL/NERSC, USA \n
!> Christof Voemel, University of California, Berkeley, USA
!
!  =====================================================================
   SUBROUTINE CLAR1V( N, B1, BN, LAMBDA, D, L, LD, LLD, &
              PIVMIN, GAPTOL, Z, WANTNC, NEGCNT, ZTZ, MINGMA, &
              R, ISUPPZ, NRMINV, RESID, RQCORR, WORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            WANTNC
   INTEGER   B1, BN, N, NEGCNT, R
   REAL               GAPTOL, LAMBDA, MINGMA, NRMINV, PIVMIN, RESID, &
                      RQCORR, ZTZ
!     ..
!     .. Array Arguments ..
   INTEGER            ISUPPZ( * )
   REAL               D( * ), L( * ), LD( * ), LLD( * ), &
                     WORK( * )
   COMPLEX          Z( * )
!     ..
!
!  =====================================================================

!     ..
!     .. Local Scalars ..
   LOGICAL            SAWNAN1, SAWNAN2
   INTEGER            I, INDLPL, INDP, INDS, INDUMN, NEG1, NEG2, R1, &
                      R2
   REAL               DMINUS, DPLUS, EPS, S, TMP
!     ..
!     .. External Functions ..
   LOGICAL SISNAN
   REAL               SLAMCH
   EXTERNAL           SISNAN, SLAMCH
!     ..
!     .. Executable Statements ..
!
   EPS = SLAMCH( 'Precision' )


   IF( R == 0 ) THEN
      R1 = B1
      R2 = BN
   ELSE
      R1 = R
      R2 = R
   END IF

!     Storage for LPLUS
   INDLPL = 0
!     Storage for UMINUS
   INDUMN = N
   INDS = 2*N + 1
   INDP = 3*N + 1

   IF( B1 == 1 ) THEN
      WORK( INDS ) = 0.0E+0
   ELSE
      WORK( INDS+B1-1 ) = LLD( B1-1 )
   END IF

!
!     Compute the stationary transform (using the differential form)
!     until the index R2.
!
   SAWNAN1 = .FALSE.
   NEG1 = 0
   S = WORK( INDS+B1-1 ) - LAMBDA
   DO I = B1, R1 - 1
      DPLUS = D( I ) + S
      WORK( INDLPL+I ) = LD( I ) / DPLUS
      IF(DPLUS < 0.0E+0) NEG1 = NEG1 + 1
      WORK( INDS+I ) = S*WORK( INDLPL+I )*L( I )
      S = WORK( INDS+I ) - LAMBDA
      ENDDO
   SAWNAN1 = SISNAN( S )
   IF( SAWNAN1 ) GOTO 60
   DO I = R1, R2 - 1
      DPLUS = D( I ) + S
      WORK( INDLPL+I ) = LD( I ) / DPLUS
      WORK( INDS+I ) = S*WORK( INDLPL+I )*L( I )
      S = WORK( INDS+I ) - LAMBDA
      ENDDO
   SAWNAN1 = SISNAN( S )
!
 60   CONTINUE
   IF( SAWNAN1 ) THEN
!        Runs a slower version of the above loop if a NaN is detected
      NEG1 = 0
      S = WORK( INDS+B1-1 ) - LAMBDA
      DO I = B1, R1 - 1
         DPLUS = D( I ) + S
         IF(ABS(DPLUS) < PIVMIN) DPLUS = -PIVMIN
         WORK( INDLPL+I ) = LD( I ) / DPLUS
         IF(DPLUS < 0.0E+0) NEG1 = NEG1 + 1
         WORK( INDS+I ) = S*WORK( INDLPL+I )*L( I )
         IF( WORK( INDLPL+I ) == 0.0E+0 ) WORK( INDS+I ) = LLD( I )
         S = WORK( INDS+I ) - LAMBDA
         ENDDO
      DO I = R1, R2 - 1
         DPLUS = D( I ) + S
         IF(ABS(DPLUS) < PIVMIN) DPLUS = -PIVMIN
         WORK( INDLPL+I ) = LD( I ) / DPLUS
         WORK( INDS+I ) = S*WORK( INDLPL+I )*L( I )
         IF( WORK( INDLPL+I ) == 0.0E+0 ) &
                         WORK( INDS+I ) = LLD( I )
         S = WORK( INDS+I ) - LAMBDA
         ENDDO
   END IF
!
!     Compute the progressive transform (using the differential form)
!     until the index R1
!
   SAWNAN2 = .FALSE.
   NEG2 = 0
   WORK( INDP+BN-1 ) = D( BN ) - LAMBDA
   DO I = BN - 1, R1, -1
      DMINUS = LLD( I ) + WORK( INDP+I )
      TMP = D( I ) / DMINUS
      IF(DMINUS < 0.0E+0) NEG2 = NEG2 + 1
      WORK( INDUMN+I ) = L( I )*TMP
      WORK( INDP+I-1 ) = WORK( INDP+I )*TMP - LAMBDA
   ENDDO
   TMP = WORK( INDP+R1-1 )
   SAWNAN2 = SISNAN( TMP )

   IF( SAWNAN2 ) THEN
!        Runs a slower version of the above loop if a NaN is detected
      NEG2 = 0
      DO I = BN-1, R1, -1
         DMINUS = LLD( I ) + WORK( INDP+I )
         IF(ABS(DMINUS) < PIVMIN) DMINUS = -PIVMIN
         TMP = D( I ) / DMINUS
         IF(DMINUS < 0.0E+0) NEG2 = NEG2 + 1
         WORK( INDUMN+I ) = L( I )*TMP
         WORK( INDP+I-1 ) = WORK( INDP+I )*TMP - LAMBDA
         IF( TMP == 0.0E+0 ) &
             WORK( INDP+I-1 ) = D( I ) - LAMBDA
      ENDDO
   END IF
!
!     Find the index (from R1 to R2) of the largest (in magnitude)
!     diagonal element of the inverse
!
   MINGMA = WORK( INDS+R1-1 ) + WORK( INDP+R1-1 )
   IF( MINGMA < 0.0E+0 ) NEG1 = NEG1 + 1
   IF( WANTNC ) THEN
      NEGCNT = NEG1 + NEG2
   ELSE
      NEGCNT = -1
   ENDIF
   IF( ABS(MINGMA) == 0.0E+0 ) MINGMA = EPS*WORK( INDS+R1-1 )
   R = R1
   DO I = R1, R2 - 1
      TMP = WORK( INDS+I ) + WORK( INDP+I )
      IF( TMP == 0.0E+0 ) TMP = EPS*WORK( INDS+I )
      IF( ABS( TMP ) <= ABS( MINGMA ) ) THEN
         MINGMA = TMP
         R = I + 1
      END IF
   ENDDO
!
!     Compute the FP vector: solve N^T v = e_r
!
   ISUPPZ( 1 ) = B1
   ISUPPZ( 2 ) = BN
   Z( R ) = (1.0E+0,0.0E+0)
   ZTZ = 1.0E+0
!
!     Compute the FP vector upwards from R
!
   IF( .NOT.SAWNAN1 .AND. .NOT.SAWNAN2 ) THEN
      DO I = R-1, B1, -1
         Z( I ) = -( WORK( INDLPL+I )*Z( I+1 ) )
         IF( (ABS(Z(I))+ABS(Z(I+1)))* ABS(LD(I)) < GAPTOL ) &
              THEN
            Z( I ) = 0.0E+0
            ISUPPZ( 1 ) = I + 1
            EXIT
         ENDIF
         ZTZ = ZTZ + REAL( Z( I )*Z( I ) )
      ENDDO
   ELSE
!        Run slower loop if NaN occurred.
      DO I = R - 1, B1, -1
         IF( Z( I+1 ) == 0.0E+0 ) THEN
            Z( I ) = -( LD( I+1 ) / LD( I ) )*Z( I+2 )
         ELSE
            Z( I ) = -( WORK( INDLPL+I )*Z( I+1 ) )
         END IF
         IF( (ABS(Z(I))+ABS(Z(I+1)))* ABS(LD(I)) < GAPTOL ) THEN
            Z( I ) = 0.0E+0
            ISUPPZ( 1 ) = I + 1
            EXIT
         END IF
         ZTZ = ZTZ + REAL( Z( I )*Z( I ) )
      ENDDO
   ENDIF

!     Compute the FP vector downwards from R in blocks of size BLKSIZ
   IF( .NOT.SAWNAN1 .AND. .NOT.SAWNAN2 ) THEN
      DO I = R, BN-1
         Z( I+1 ) = -( WORK( INDUMN+I )*Z( I ) )
         IF( (ABS(Z(I))+ABS(Z(I+1)))* ABS(LD(I)) < GAPTOL ) &
            THEN
            Z( I+1 ) = 0.0E+0
            ISUPPZ( 2 ) = I
            EXIT
         END IF
         ZTZ = ZTZ + REAL( Z( I+1 )*Z( I+1 ) )
      ENDDO
   ELSE
!        Run slower loop if NaN occurred.
      DO I = R, BN - 1
         IF( Z( I ) == 0.0E+0 ) THEN
            Z( I+1 ) = -( LD( I-1 ) / LD( I ) )*Z( I-1 )
         ELSE
            Z( I+1 ) = -( WORK( INDUMN+I )*Z( I ) )
         END IF
         IF( (ABS(Z(I))+ABS(Z(I+1)))* ABS(LD(I)) < GAPTOL ) &
              THEN
            Z( I+1 ) = 0.0E+0
            ISUPPZ( 2 ) = I
            EXIT
         END IF
         ZTZ = ZTZ + REAL( Z( I+1 )*Z( I+1 ) )
      ENDDO
   END IF
!
!     Compute quantities for convergence test
!
   TMP = 1.0E+0 / ZTZ
   NRMINV = SQRT( TMP )
   RESID = ABS( MINGMA )*NRMINV
   RQCORR = MINGMA*TMP
!
!
   RETURN
!
!     End of CLAR1V
!
END
