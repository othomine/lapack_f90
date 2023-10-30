!> \brief \b CLAIC1 applies one step of incremental condition estimation.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAIC1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claic1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claic1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claic1.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C )
!
!       .. Scalar Arguments ..
!       INTEGER            J, JOB
!       REAL               SEST, SESTPR
!       COMPLEX            C, GAMMA, S
!       ..
!       .. Array Arguments ..
!       COMPLEX            W( J ), X( J )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAIC1 applies one step of incremental condition estimation in
!> its simplest version:
!>
!> Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j
!> lower triangular matrix L, such that
!>          twonorm(L*x) = sest
!> Then CLAIC1 computes sestpr, s, c such that
!> the vector
!>                 [ s*x ]
!>          xhat = [  c  ]
!> is an approximate singular vector of
!>                 [ L      0  ]
!>          Lhat = [ w**H gamma ]
!> in the sense that
!>          twonorm(Lhat*xhat) = sestpr.
!>
!> Depending on JOB, an estimate for the largest or smallest singular
!> value is computed.
!>
!> Note that [s c]**H and sestpr**2 is an eigenpair of the system
!>
!>     diag(sest*sest, 0) + [alpha  gamma] * [ conjg(alpha) ]
!>                                           [ conjg(gamma) ]
!>
!> where  alpha =  x**H*w.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is INTEGER
!>          = 1: an estimate for the largest singular value is computed.
!>          = 2: an estimate for the smallest singular value is computed.
!> \endverbatim
!>
!> \param[in] J
!> \verbatim
!>          J is INTEGER
!>          Length of X and W
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX array, dimension (J)
!>          The j-vector x.
!> \endverbatim
!>
!> \param[in] SEST
!> \verbatim
!>          SEST is REAL
!>          Estimated singular value of j by j matrix L
!> \endverbatim
!>
!> \param[in] W
!> \verbatim
!>          W is COMPLEX array, dimension (J)
!>          The j-vector w.
!> \endverbatim
!>
!> \param[in] GAMMA
!> \verbatim
!>          GAMMA is COMPLEX
!>          The diagonal element gamma.
!> \endverbatim
!>
!> \param[out] SESTPR
!> \verbatim
!>          SESTPR is REAL
!>          Estimated singular value of (j+1) by (j+1) matrix Lhat.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is COMPLEX
!>          Sine needed in forming xhat.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is COMPLEX
!>          Cosine needed in forming xhat.
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
!> \ingroup laic1
!
!  =====================================================================
   SUBROUTINE CLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            J, JOB
   REAL               SEST, SESTPR
   COMPLEX            C, GAMMA, S
!     ..
!     .. Array Arguments ..
   COMPLEX            W( J ), X( J )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   REAL               ABSALP, ABSEST, ABSGAM, B, EPS, NORMA, S1, S2, &
                      SCL, T, TEST, TMP, ZETA1, ZETA2
   COMPLEX            ALPHA, COSINE, SINE
!     ..
!     .. External Functions ..
   REAL               SLAMCH
   COMPLEX            CDOTC
   EXTERNAL           SLAMCH, CDOTC
!     ..
!     .. Executable Statements ..
!
   EPS = SLAMCH( 'Epsilon' )
   ALPHA = CDOTC( J, X, 1, W, 1 )
!
   ABSALP = ABS( ALPHA )
   ABSGAM = ABS( GAMMA )
   ABSEST = ABS( SEST )
!
   IF( JOB == 1 ) THEN
!
!        Estimating largest singular value
!
!        special cases
!
      IF( SEST == 0.0E+0 ) THEN
         S1 = MAX( ABSGAM, ABSALP )
         IF( S1 == 0.0E+0 ) THEN
            S = 0.0E+0
            C = 1.0E+0
            SESTPR = 0.0E+0
         ELSE
            S = ALPHA / S1
            C = GAMMA / S1
            TMP = REAL( SQRT( S*CONJG( S )+C*CONJG( C ) ) )
            S = S / TMP
            C = C / TMP
            SESTPR = S1*TMP
         END IF
         RETURN
      ELSE IF( ABSGAM <= EPS*ABSEST ) THEN
         S = 1.0E+0
         C = 0.0E+0
         TMP = MAX( ABSEST, ABSALP )
         S1 = ABSEST / TMP
         S2 = ABSALP / TMP
         SESTPR = TMP*SQRT( S1*S1+S2*S2 )
         RETURN
      ELSE IF( ABSALP <= EPS*ABSEST ) THEN
         S1 = ABSGAM
         S2 = ABSEST
         IF( S1 <= S2 ) THEN
            S = 1.0E+0
            C = 0.0E+0
            SESTPR = S2
         ELSE
            S = 0.0E+0
            C = 1.0E+0
            SESTPR = S1
         END IF
         RETURN
      ELSE IF( ABSEST <= EPS*ABSALP .OR. ABSEST <= EPS*ABSGAM ) THEN
         S1 = ABSGAM
         S2 = ABSALP
         IF( S1 <= S2 ) THEN
            TMP = S1 / S2
            SCL = SQRT( 1.0E+0+TMP*TMP )
            SESTPR = S2*SCL
            S = ( ALPHA / S2 ) / SCL
            C = ( GAMMA / S2 ) / SCL
         ELSE
            TMP = S2 / S1
            SCL = SQRT( 1.0E+0+TMP*TMP )
            SESTPR = S1*SCL
            S = ( ALPHA / S1 ) / SCL
            C = ( GAMMA / S1 ) / SCL
         END IF
         RETURN
      ELSE
!
!           normal case
!
         ZETA1 = ABSALP / ABSEST
         ZETA2 = ABSGAM / ABSEST
!
         B = ( 1.0E+0-ZETA1*ZETA1-ZETA2*ZETA2 )*0.5E+0
         C = ZETA1*ZETA1
         IF( B > 0.0E+0 ) THEN
            T = REAL( C / ( B+SQRT( B*B+C ) ) )
         ELSE
            T = REAL( SQRT( B*B+C ) - B )
         END IF
!
         SINE = -( ALPHA / ABSEST ) / T
         COSINE = -( GAMMA / ABSEST ) / ( 1.0E+0+T )
         TMP = REAL( SQRT( SINE * CONJG( SINE ) + COSINE * CONJG( COSINE ) ) )
         S = SINE / TMP
         C = COSINE / TMP
         SESTPR = SQRT( T+1.0E+0 )*ABSEST
         RETURN
      END IF
!
   ELSE IF( JOB == 2 ) THEN
!
!        Estimating smallest singular value
!
!        special cases
!
      IF( SEST == 0.0E+0 ) THEN
         SESTPR = 0.0E+0
         IF( MAX( ABSGAM, ABSALP ) == 0.0E+0 ) THEN
            SINE = 1.0E+0
            COSINE = 0.0E+0
         ELSE
            SINE = -CONJG( GAMMA )
            COSINE = CONJG( ALPHA )
         END IF
         S1 = MAX( ABS( SINE ), ABS( COSINE ) )
         S = SINE / S1
         C = COSINE / S1
         TMP = REAL( SQRT( S*CONJG( S )+C*CONJG( C ) ) )
         S = S / TMP
         C = C / TMP
         RETURN
      ELSE IF( ABSGAM <= EPS*ABSEST ) THEN
         S = 0.0E+0
         C = 1.0E+0
         SESTPR = ABSGAM
         RETURN
      ELSE IF( ABSALP <= EPS*ABSEST ) THEN
         S1 = ABSGAM
         S2 = ABSEST
         IF( S1 <= S2 ) THEN
            S = 0.0E+0
            C = 1.0E+0
            SESTPR = S1
         ELSE
            S = 1.0E+0
            C = 0.0E+0
            SESTPR = S2
         END IF
         RETURN
      ELSE IF( ABSEST <= EPS*ABSALP .OR. ABSEST <= EPS*ABSGAM ) THEN
         S1 = ABSGAM
         S2 = ABSALP
         IF( S1 <= S2 ) THEN
            TMP = S1 / S2
            SCL = SQRT( 1.0E+0+TMP*TMP )
            SESTPR = ABSEST*( TMP / SCL )
            S = -( CONJG( GAMMA ) / S2 ) / SCL
            C = ( CONJG( ALPHA ) / S2 ) / SCL
         ELSE
            TMP = S2 / S1
            SCL = SQRT( 1.0E+0+TMP*TMP )
            SESTPR = ABSEST / SCL
            S = -( CONJG( GAMMA ) / S1 ) / SCL
            C = ( CONJG( ALPHA ) / S1 ) / SCL
         END IF
         RETURN
      ELSE
!
!           normal case
!
         ZETA1 = ABSALP / ABSEST
         ZETA2 = ABSGAM / ABSEST
!
         NORMA = MAX( 1.0E+0+ZETA1*ZETA1+ZETA1*ZETA2, &
                 ZETA1*ZETA2+ZETA2*ZETA2 )
!
!           See if root is closer to zero or to 1.0E+0
!
         TEST = 1.0E+0 + 2.0E+0*( ZETA1-ZETA2 )*( ZETA1+ZETA2 )
         IF( TEST >= 0.0E+0 ) THEN
!
!              root is close to zero, compute directly
!
            B = ( ZETA1*ZETA1+ZETA2*ZETA2+1.0E+0 )*0.5E+0
            C = ZETA2*ZETA2
            T = REAL( C / ( B+SQRT( ABS( B*B-C ) ) ) )
            SINE = ( ALPHA / ABSEST ) / ( 1.0E+0-T )
            COSINE = -( GAMMA / ABSEST ) / T
            SESTPR = SQRT( T+4.0E+0*EPS*EPS*NORMA )*ABSEST
         ELSE
!
!              root is closer to 1.0E+0, shift by that amount
!
            B = ( ZETA2*ZETA2+ZETA1*ZETA1-1.0E+0 )*0.5E+0
            C = ZETA1*ZETA1
            IF( B >= 0.0E+0 ) THEN
               T = REAL( -C / ( B+SQRT( B*B+C ) ) )
            ELSE
               T = REAL( B - SQRT( B*B+C ) )
            END IF
            SINE = -( ALPHA / ABSEST ) / T
            COSINE = -( GAMMA / ABSEST ) / ( 1.0E+0+T )
            SESTPR = SQRT( 1.0E+0+T+4.0E+0*EPS*EPS*NORMA )*ABSEST
         END IF
         TMP = REAL( SQRT( SINE * CONJG( SINE ) &
           + COSINE * CONJG( COSINE ) ) )
         S = SINE / TMP
         C = COSINE / TMP
         RETURN
!
      END IF
   END IF
   RETURN
!
!     End of CLAIC1
!
END
