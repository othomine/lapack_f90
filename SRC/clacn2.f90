!> \brief \b CLACN2 estimates the 1-norm of a square matrix, using reverse communication for evaluating matrix-vector products.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLACN2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacn2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacn2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacn2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLACN2( N, V, X, EST, KASE, ISAVE )
!
!       .. Scalar Arguments ..
!       INTEGER            KASE, N
!       REAL               EST
!       ..
!       .. Array Arguments ..
!       INTEGER            ISAVE( 3 )
!       COMPLEX            V( * ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLACN2 estimates the 1-norm of a square, complex matrix A.
!> Reverse communication is used for evaluating matrix-vector products.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The order of the matrix.  N >= 1.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX array, dimension (N)
!>         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
!>         (W is not returned).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX array, dimension (N)
!>         On an intermediate return, X should be overwritten by
!>               A * X,   if KASE=1,
!>               A**H * X,  if KASE=2,
!>         where A**H is the conjugate transpose of A, and CLACN2 must be
!>         re-called with all the other parameters unchanged.
!> \endverbatim
!>
!> \param[in,out] EST
!> \verbatim
!>          EST is REAL
!>         On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be
!>         unchanged from the previous call to CLACN2.
!>         On exit, EST is an estimate (a lower bound) for norm(A).
!> \endverbatim
!>
!> \param[in,out] KASE
!> \verbatim
!>          KASE is INTEGER
!>         On the initial call to CLACN2, KASE should be 0.
!>         On an intermediate return, KASE will be 1 or 2, indicating
!>         whether X should be overwritten by A * X  or A**H * X.
!>         On the final return from CLACN2, KASE will again be 0.
!> \endverbatim
!>
!> \param[in,out] ISAVE
!> \verbatim
!>          ISAVE is INTEGER array, dimension (3)
!>         ISAVE is used to save variables between calls to SLACN2
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
!> \ingroup lacn2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Originally named CONEST, dated March 16, 1988.
!>
!>  Last modified:  April, 1999
!>
!>  This is a thread safe version of CLACON, which uses the array ISAVE
!>  in place of a SAVE statement, as follows:
!>
!>     CLACON     CLACN2
!>      JUMP     ISAVE(1)
!>      J        ISAVE(2)
!>      ITER     ISAVE(3)
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!>     Nick Higham, University of Manchester
!
!> \par References:
!  ================
!>
!>  N.J. Higham, "FORTRAN codes for estimating the one-norm of
!>  a real or complex matrix, with applications to condition estimation",
!>  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
!>
!  =====================================================================
   SUBROUTINE CLACN2( N, V, X, EST, KASE, ISAVE )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            KASE, N
   REAL               EST
!     ..
!     .. Array Arguments ..
   INTEGER            ISAVE( 3 )
   COMPLEX            V( * ), X( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER              ITMAX
   PARAMETER          ( ITMAX = 5 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, JLAST
   REAL               ABSXI, ALTSGN, ESTOLD, SAFMIN, TEMP
!     ..
!     .. External Functions ..
   INTEGER            ICMAX1
   REAL               SCSUM1, SLAMCH
   EXTERNAL           ICMAX1, SCSUM1, SLAMCH
!     ..
!     .. Executable Statements ..
!
   SAFMIN = SLAMCH( 'Safe minimum' )
   IF( KASE == 0 ) THEN
      X(1:N) = CMPLX( 1.0E+0 / REAL( N ) )
      KASE = 1
      ISAVE( 1 ) = 1
      RETURN
   END IF
!
   SELECT CASE (ISAVE(1))
!
!     ................ ENTRY   (ISAVE( 1 ) = 1)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
!
    CASE (1)
     IF( N == 1 ) THEN
        V( 1 ) = X( 1 )
        EST = ABS( V( 1 ) )
!        ... QUIT
        KASE = 0
        RETURN
     END IF
     EST = SCSUM1( N, X, 1 )
!
     WHERE (ABS(X(1:N)) > SAFMIN )
        X(1:N) = CMPLX(REAL(X(1:N) )/ABS(X(1:N)),AIMAG(X(1:N))/ABS(X(1:N)) )
     ELSEWHERE
        X(1:N) = (1.0E+0,0.0E+0)
     ENDWHERE
     KASE = 2
     ISAVE( 1 ) = 2
!
!     ................ ENTRY   (ISAVE( 1 ) = 2)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
!
    CASE (2)
     ISAVE( 2 ) = ICMAX1( N, X, 1 )
     ISAVE( 3 ) = 2
!
!     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
!
     X(1:N) = (0.0E+0,0.0E+0)
     X( ISAVE( 2 ) ) = (1.0E+0,0.0E+0)
     KASE = 1
     ISAVE( 1 ) = 3
!
!     ................ ENTRY   (ISAVE( 1 ) = 3)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
    CASE (3)
     V(1:N) = X(1:N)
     ESTOLD = EST
     EST = SCSUM1( N, V, 1 )
!
!     TEST FOR CYCLING.
     IF( EST <= ESTOLD ) THEN
        ALTSGN = 1.0E+0
        DO I = 1, N
           X( I ) = CMPLX( ALTSGN*( 1.0E+0 + REAL( I-1 ) / REAL( N-1 ) ) )
           ALTSGN = -ALTSGN
        ENDDO
        KASE = 1
        ISAVE( 1 ) = 5
        RETURN
     ENDIF
!
     WHERE (ABS(X(1:N)) > SAFMIN )
        X(1:N) = CMPLX(REAL(X(1:N) )/ABS(X(1:N)),AIMAG(X(1:N))/ABS(X(1:N)) )
     ELSEWHERE
        X(1:N) = (1.0E+0,0.0E+0)
     ENDWHERE
     KASE = 2
     ISAVE( 1 ) = 4
!
!     ................ ENTRY   (ISAVE( 1 ) = 4)
!     X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
!
    CASE (4)
     JLAST = ISAVE( 2 )
     ISAVE( 2 ) = ICMAX1( N, X, 1 )
     IF( ( ABS( X( JLAST ) ) /= ABS( X( ISAVE( 2 ) ) ) ) .AND. &
         ( ISAVE( 3 ) < ITMAX ) ) THEN
        ISAVE( 3 ) = ISAVE( 3 ) + 1
        X(1:N) = (0.0E+0,0.0E+0)
        X( ISAVE( 2 ) ) = (1.0E+0,0.0E+0)
        KASE = 1
        ISAVE( 1 ) = 3
        RETURN
     END IF
!
!     ITERATION COMPLETE.  FINAL STAGE.
!
     ALTSGN = 1.0E+0
     DO I = 1, N
        X( I ) = CMPLX( ALTSGN*( 1.0E+0 + REAL( I-1 ) / REAL( N-1 ) ) )
        ALTSGN = -ALTSGN
     ENDDO
     KASE = 1
     ISAVE( 1 ) = 5
!
!     ................ ENTRY   (ISAVE( 1 ) = 5)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
    CASE (5)
     TEMP = 2.0E+0*( SCSUM1( N, X, 1 ) / REAL( 3*N ) )
     IF( TEMP > EST ) THEN
        V(1:N) = X(1:N)
        EST = TEMP
     END IF
     KASE = 0
!
   END SELECT
   RETURN
!
!     End of CLACN2
!
END
