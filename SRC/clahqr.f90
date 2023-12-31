!> \brief \b CLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the double-shift/single-shift QR algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAHQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahqr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
!                          IHIZ, Z, LDZ, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       COMPLEX            H( LDH, * ), W( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CLAHQR is an auxiliary routine called by CHSEQR to update the
!>    eigenvalues and Schur decomposition already computed by CHSEQR, by
!>    dealing with the Hessenberg submatrix in rows and columns ILO to
!>    IHI.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>          = .TRUE. : the full Schur form T is required;
!>          = .FALSE.: only eigenvalues are required.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          = .TRUE. : the matrix of Schur vectors Z is required;
!>          = .FALSE.: Schur vectors are not required.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix H.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>          It is assumed that H is already upper triangular in rows and
!>          columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1).
!>          CLAHQR works primarily with the Hessenberg submatrix in rows
!>          and columns ILO to IHI, but applies transformations to all of
!>          H if WANTT is .TRUE..
!>          1 <= ILO <= max(1,IHI); IHI <= N.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX array, dimension (LDH,N)
!>          On entry, the upper Hessenberg matrix H.
!>          On exit, if INFO is zero and if WANTT is .TRUE., then H
!>          is upper triangular in rows and columns ILO:IHI.  If INFO
!>          is zero and if WANTT is .FALSE., then the contents of H
!>          are unspecified on exit.  The output state of H in case
!>          INF is positive is below under the description of INFO.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          The leading dimension of the array H. LDH >= max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX array, dimension (N)
!>          The computed eigenvalues ILO to IHI are stored in the
!>          corresponding elements of W. If WANTT is .TRUE., the
!>          eigenvalues are stored in the same order as on the diagonal
!>          of the Schur form returned in H, with W(i) = H(i,i).
!> \endverbatim
!>
!> \param[in] ILOZ
!> \verbatim
!>          ILOZ is INTEGER
!> \endverbatim
!>
!> \param[in] IHIZ
!> \verbatim
!>          IHIZ is INTEGER
!>          Specify the rows of Z to which transformations must be
!>          applied if WANTZ is .TRUE..
!>          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ,N)
!>          If WANTZ is .TRUE., on entry Z must contain the current
!>          matrix Z of transformations accumulated by CHSEQR, and on
!>          exit Z has been updated; transformations are applied only to
!>          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
!>          If WANTZ is .FALSE., Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z. LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           = 0:  successful exit
!>           > 0:  if INFO = i, CLAHQR failed to compute all the
!>                  eigenvalues ILO to IHI in a total of 30 iterations
!>                  per eigenvalue; elements i+1:ihi of W contain
!>                  those eigenvalues which have been successfully
!>                  computed.
!>
!>                  If INFO > 0 and WANTT is .FALSE., then on exit,
!>                  the remaining unconverged eigenvalues are the
!>                  eigenvalues of the upper Hessenberg matrix
!>                  rows and columns ILO through INFO of the final,
!>                  output value of H.
!>
!>                  If INFO > 0 and WANTT is .TRUE., then on exit
!>          (*)       (initial value of H)*U  = U*(final value of H)
!>                  where U is an orthogonal matrix.    The final
!>                  value of H is upper Hessenberg and triangular in
!>                  rows and columns INFO+1 through IHI.
!>
!>                  If INFO > 0 and WANTZ is .TRUE., then on exit
!>                      (final value of Z)  = (initial value of Z)*U
!>                  where U is the orthogonal matrix in (*)
!>                  (regardless of the value of WANTT.)
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
!> \ingroup lahqr
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>     02-96 Based on modifications by
!>     David Day, Sandia National Laboratory, USA
!>
!>     12-04 Further modifications by
!>     Ralph Byers, University of Kansas, USA
!>     This is a modified version of CLAHQR from LAPACK version 3.0.
!>     It is (1) more robust against overflow and underflow and
!>     (2) adopts the more conservative Ahues & Tisseur stopping
!>     criterion (LAWN 122, 1997).
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, &
                      IHIZ, Z, LDZ, INFO )
   IMPLICIT NONE
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
   LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
   COMPLEX            H( LDH, * ), W( * ), Z( LDZ, * )
!     ..
!
!  =========================================================
!
!     .. Parameters ..
   REAL               DAT1
   PARAMETER          ( DAT1 = 3.0e0 / 4.0e0 )
   INTEGER            KEXSH
   PARAMETER          ( KEXSH = 10 )
!     ..
!     .. Local Scalars ..
   COMPLEX            CDUM, H11, H11S, H22, SC, SUM, T, T1, TEMP, U, &
                      V2, X, Y
   REAL               AA, AB, BA, BB, H10, H21, RTEMP, S, SAFMAX, &
                      SAFMIN, SMLNUM, SX, T2, TST, ULP
   INTEGER            I, I1, I2, ITS, ITMAX, J, JHI, JLO, K, L, M, &
                      NH, NZ, KDEFL
!     ..
!     .. Local Arrays ..
   COMPLEX            V( 2 )
!     ..
!     .. External Functions ..
   COMPLEX            CLADIV
   REAL               SLAMCH
   EXTERNAL           CLADIV, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           CCOPY, CLARFG, CSCAL
!     ..
!     .. Statement Functions ..
   REAL               CABS1
!     ..
!     .. Statement Function definitions ..
   CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
   INFO = 0
!
!     Quick return if possible
!
   IF( N == 0 ) RETURN
   IF( ILO == IHI ) THEN
      W( ILO ) = H( ILO, ILO )
      RETURN
   END IF
!
!     ==== clear out the trash ====
   DO J = ILO, IHI - 3
      H( J+2, J ) = (0.0E+0,0.0E+0)
      H( J+3, J ) = (0.0E+0,0.0E+0)
   ENDDO
   IF( ILO <= IHI-2 ) H( IHI, IHI-2 ) = (0.0E+0,0.0E+0)
!     ==== ensure that subdiagonal entries are real ====
   IF( WANTT ) THEN
      JLO = 1
      JHI = N
   ELSE
      JLO = ILO
      JHI = IHI
   END IF
   DO I = ILO + 1, IHI
      IF( AIMAG( H( I, I-1 ) ) /= 0.0E+0 ) THEN
!           ==== The following redundant normalization
!           .    avoids problems with both gradual and
!           .    sudden underflow in ABS(H(I,I-1)) ====
         SC = H( I, I-1 ) / CABS1( H( I, I-1 ) )
         SC = CONJG( SC ) / ABS( SC )
         H( I, I-1 ) = ABS( H( I, I-1 ) )
         H(I,I:JHI) = SC*H(I,I:JHI)
         H(JLO:MIN(JHI,I+1),I) = CONJG(SC)*H(JLO:MIN(JHI,I+1),I)
         IF( WANTZ ) Z(ILOZ:IHIZ,I) = CONJG(SC)*Z(ILOZ:IHIZ,I)
      END IF
   ENDDO
!
   NH = IHI - ILO + 1
   NZ = IHIZ - ILOZ + 1
!
!     Set machine-dependent constants for the stopping criterion.
!
   SAFMIN = SLAMCH( 'SAFE MINIMUM' )
   SAFMAX = 1.0E+0 / SAFMIN
   ULP = SLAMCH( 'PRECISION' )
   SMLNUM = SAFMIN*( REAL( NH ) / ULP )
!
!     I1 and I2 are the indices of the first row and last column of H
!     to which transformations must be applied. If eigenvalues only are
!     being computed, I1 and I2 are set inside the main loop.
!
   IF( WANTT ) THEN
      I1 = 1
      I2 = N
   END IF
!
!     ITMAX is the total number of QR iterations allowed.
!
   ITMAX = 30 * MAX( 10, NH )
!
!     KDEFL counts the number of iterations since a deflation
!
   KDEFL = 0
!
!     The main loop begins here. I is the loop index and decreases from
!     IHI to ILO in steps of 1. Each iteration of the loop works
!     with the active submatrix in rows and columns L to I.
!     Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
!     H(L,L-1) is negligible so that the matrix splits.
!
   I = IHI
30 CONTINUE
   IF( I < ILO ) GO TO 150
!
!     Perform QR iterations on rows and columns ILO to I until a
!     submatrix of order 1 splits off at the bottom because a
!     subdiagonal element has become negligible.
!
   L = ILO
   DO ITS = 0, ITMAX
!
!        Look for a single small subdiagonal element.
!
      DO K = I, L + 1, -1
         IF( CABS1( H( K, K-1 ) ) <= SMLNUM ) GO TO 50
         TST = CABS1( H( K-1, K-1 ) ) + CABS1( H( K, K ) )
         IF( TST == (0.0E+0,0.0E+0) ) THEN
            IF( K-2 >= ILO ) TST = TST + ABS( REAL( H( K-1, K-2 ) ) )
            IF( K+1 <= IHI ) TST = TST + ABS( REAL( H( K+1, K ) ) )
         END IF
!           ==== The following is a conservative small subdiagonal
!           .    deflation criterion due to Ahues & Tisseur (LAWN 122,
!           .    1997). It has better mathematical foundation and
!           .    improves accuracy in some examples.  ====
         IF( ABS( REAL( H( K, K-1 ) ) ) <= ULP*TST ) THEN
            AB = MAX( CABS1( H( K, K-1 ) ), CABS1( H( K-1, K ) ) )
            BA = MIN( CABS1( H( K, K-1 ) ), CABS1( H( K-1, K ) ) )
            AA = MAX( CABS1( H( K, K ) ), &
                 CABS1( H( K-1, K-1 )-H( K, K ) ) )
            BB = MIN( CABS1( H( K, K ) ), &
                 CABS1( H( K-1, K-1 )-H( K, K ) ) )
            S = AA + AB
            IF( BA*( AB / S ) <= MAX( SMLNUM, &
                ULP*( BB*( AA / S ) ) ) )GO TO 50
         END IF
      ENDDO
50    CONTINUE
      L = K
      IF( L > ILO ) THEN
!
!           H(L,L-1) is negligible
!
         H( L, L-1 ) = (0.0E+0,0.0E+0)
      END IF
!
!        Exit from loop if a submatrix of order 1 has split off.
!
      IF( L >= I ) GO TO 140
      KDEFL = KDEFL + 1
!
!        Now the active submatrix is in rows and columns L to I. If
!        eigenvalues only are being computed, only the active submatrix
!        need be transformed.
!
      IF( .NOT.WANTT ) THEN
         I1 = L
         I2 = I
      END IF
!
      IF( MOD(KDEFL,2*KEXSH) == 0 ) THEN
!
!           Exceptional shift.
!
         S = DAT1*ABS( REAL( H( I, I-1 ) ) )
         T = S + H( I, I )
      ELSE IF( MOD(KDEFL,KEXSH) == 0 ) THEN
!
!           Exceptional shift.
!
         S = DAT1*ABS( REAL( H( L+1, L ) ) )
         T = S + H( L, L )
      ELSE
!
!           Wilkinson's shift.
!
         T = H( I, I )
         U = SQRT( H( I-1, I ) * H( I, I-1 ) )
         S = CABS1( U )
         IF( S /= 0.0E+0 ) THEN
            X = 0.5E+0*( H( I-1, I-1 )-T )
            SX = CABS1( X )
            S = MAX( S, CABS1( X ) )
            Y = SQRT( X**2+U**2 )
            IF( SX > 0.0E+0 ) THEN
               IF( REAL( X / SX )*REAL( Y )+AIMAG( X / SX )* &
                   AIMAG( Y ) < 0.0E+0 )Y = -Y
            END IF
            T = T - U*CLADIV( U, ( X+Y ) )
         END IF
      END IF
!
!        Look for two consecutive small subdiagonal elements.
!
      DO M = I - 1, L + 1, -1
!
!           Determine the effect of starting the single-shift QR
!           iteration at row M, and see if this would make H(M,M-1)
!           negligible.
!
         H11 = H( M, M )
         H22 = H( M+1, M+1 )
         H11S = H11 - T
         H21 = REAL( H( M+1, M ) )
         S = CABS1( H11S ) + ABS( H21 )
         H11S = H11S / S
         H21 = H21 / S
         V( 1 ) = H11S
         V( 2 ) = H21
         H10 = REAL( H( M, M-1 ) )
         IF( ABS( H10 )*ABS( H21 ) <= ULP* &
             ( CABS1( H11S )*( CABS1( H11 )+CABS1( H22 ) ) ) ) &
             GO TO 70
      ENDDO
      H11 = H( L, L )
      H22 = H( L+1, L+1 )
      H11S = H11 - T
      H21 = REAL( H( L+1, L ) )
      S = CABS1( H11S ) + ABS( H21 )
      H11S = H11S / S
      H21 = H21 / S
      V( 1 ) = H11S
      V( 2 ) = H21
70    CONTINUE
!
!        Single-shift QR step
!
      DO K = M, I - 1
!
!           The first iteration of this loop determines a reflection G
!           from the vector V and applies it from left and right to H,
!           thus creating a nonzero bulge below the subdiagonal.
!
!           Each subsequent iteration determines a reflection G to
!           restore the Hessenberg form in the (K-1)th column, and thus
!           chases the bulge one step toward the bottom of the active
!           submatrix.
!
!           V(2) is always real before the call to CLARFG, and hence
!           after the call T2 ( = T1*V(2) ) is also real.
!
         IF( K > M ) V(1:2) = H(K:K+1,K-1)
         CALL CLARFG( 2, V( 1 ), V( 2 ), 1, T1 )
         IF( K > M ) THEN
            H( K, K-1 ) = V( 1 )
            H( K+1, K-1 ) = (0.0E+0,0.0E+0)
         END IF
         V2 = V( 2 )
         T2 = REAL( T1*V2 )
!
!           Apply G from the left to transform the rows of the matrix
!           in columns K to I2.
!
         DO J = K, I2
            SUM = CONJG( T1 )*H( K, J ) + T2*H( K+1, J )
            H( K, J ) = H( K, J ) - SUM
            H( K+1, J ) = H( K+1, J ) - SUM*V2
         ENDDO
!
!           Apply G from the right to transform the columns of the
!           matrix in rows I1 to min(K+2,I).
!
         DO J = I1, MIN( K+2, I )
            SUM = T1*H( J, K ) + T2*H( J, K+1 )
            H( J, K ) = H( J, K ) - SUM
            H( J, K+1 ) = H( J, K+1 ) - SUM*CONJG( V2 )
         ENDDO
!
         IF( WANTZ ) THEN
!
!              Accumulate transformations in the matrix Z
!
            DO J = ILOZ, IHIZ
               SUM = T1*Z( J, K ) + T2*Z( J, K+1 )
               Z( J, K ) = Z( J, K ) - SUM
               Z( J, K+1 ) = Z( J, K+1 ) - SUM*CONJG( V2 )
               ENDDO
         END IF
!
         IF( K == M .AND. M > L ) THEN
!
!              If the QR step was started at row M > L because two
!              consecutive small subdiagonals were found, then extra
!              scaling must be performed to ensure that H(M,M-1) remains
!              real.
!
            TEMP = (1.0E+0,0.0E+0) - T1
            TEMP = TEMP / ABS( TEMP )
            H( M+1, M ) = H( M+1, M )*CONJG( TEMP )
            IF( M+2 <= I ) H( M+2, M+1 ) = H( M+2, M+1 )*TEMP
            DO J = M, I
               IF( J /= M+1 ) THEN
                  IF( I2 > J ) H(J,J+1:I2) = TEMP*H(J,J+1:I2)
                  H(I1:J-1,J) = CONJG( TEMP )*H(I1:J-1,J)
                  IF( WANTZ ) Z(ILOZ:ILOZ+NZ-1,J) = CONJG( TEMP )*Z(ILOZ:ILOZ+NZ-1,J)
               END IF
            ENDDO
         END IF
      ENDDO
!
!        Ensure that H(I,I-1) is real.
!
      TEMP = H( I, I-1 )
      IF( AIMAG( TEMP ) /= 0.0E+0 ) THEN
         RTEMP = ABS( TEMP )
         H( I, I-1 ) = RTEMP
         TEMP = TEMP / RTEMP
         IF( I2 > I ) H(I,I+1:I2) = CONJG( TEMP )*H(I,I+1:I2)
         H(I1:I-1,I) = TEMP*H(I1:I-1,I)
         IF( WANTZ ) Z(ILOZ:ILOZ+NZ-1,I) = TEMP*Z(ILOZ:ILOZ+NZ-1,I)
      END IF
!
      ENDDO
!
!     Failure to converge in remaining number of iterations
!
   INFO = I
   RETURN
!
  140 CONTINUE
!
!     H(I,I-1) is negligible: one eigenvalue has converged.
!
   W( I ) = H( I, I )
!     reset deflation counter
   KDEFL = 0
!
!     return to start of the main loop with new value of I.
!
   I = L - 1
   GO TO 30
!
  150 CONTINUE
   RETURN
!
!     End of CLAHQR
!
END
