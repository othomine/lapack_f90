!> \brief \b CLACON estimates the 1-norm of a square matrix, using reverse communication for evaluating matrix-vector products.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLACON + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacon.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacon.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacon.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLACON( N, V, X, EST, KASE )
!
!       .. Scalar Arguments ..
!       INTEGER            KASE, N
!       REAL               EST
!       ..
!       .. Array Arguments ..
!       COMPLEX            V( N ), X( N )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLACON estimates the 1-norm of a square, complex matrix A.
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
!>         where A**H is the conjugate transpose of A, and CLACON must be
!>         re-called with all the other parameters unchanged.
!> \endverbatim
!>
!> \param[in,out] EST
!> \verbatim
!>          EST is REAL
!>         On entry with KASE = 1 or 2 and JUMP = 3, EST should be
!>         unchanged from the previous call to CLACON.
!>         On exit, EST is an estimate (a lower bound) for norm(A).
!> \endverbatim
!>
!> \param[in,out] KASE
!> \verbatim
!>          KASE is INTEGER
!>         On the initial call to CLACON, KASE should be 0.
!>         On an intermediate return, KASE will be 1 or 2, indicating
!>         whether X should be overwritten by A * X  or A**H * X.
!>         On the final return from CLACON, KASE will again be 0.
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
!> \ingroup lacon
!
!> \par Further Details:
!  =====================
!>
!>  Originally named CONEST, dated March 16, 1988. \n
!>  Last modified:  April, 1999
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
   SUBROUTINE CLACON( N, V, X, EST, KASE )
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
   COMPLEX            V( N ), X( N )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            ITMAX
   PARAMETER          ( ITMAX = 5 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, ITER, J, JLAST, JUMP
   REAL               ABSXI, ALTSGN, ESTOLD, SAFMIN, TEMP
!     ..
!     .. External Functions ..
   INTEGER            ICMAX1
   REAL               SCSUM1, SLAMCH
   EXTERNAL           ICMAX1, SCSUM1, SLAMCH
!     ..
!     .. Save statement ..
   SAVE
!     ..
!     .. Executable Statements ..
!
   SAFMIN = SLAMCH( 'Safe minimum' )
   IF( KASE == 0 ) THEN
      X(1:N) = CMPLX( 1.0E+0 / REAL( N ) )
      KASE = 1
      JUMP = 1
      RETURN
   END IF
!
   SELECT CASE (JUMP)
!
!     ................ ENTRY   (JUMP = 1)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
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
     JUMP = 2
!
!     ................ ENTRY   (JUMP = 2)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
!
    CASE (2)
     J = ICMAX1( N, X, 1 )
     ITER = 2
!
!     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
!
     X(1:N) = (0.0E+0,0.0E+0)
     X( J ) = (1.0E+0,0.0E+0)
     KASE = 1
     JUMP = 3
!
!     ................ ENTRY   (JUMP = 3)
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
           X( I ) = CMPLX( ALTSGN*( 1.0E+0+REAL( I-1 ) / REAL( N-1 ) ) )
           ALTSGN = -ALTSGN
        ENDDO
        KASE = 1
        JUMP = 5
        RETURN
     ENDIF
!
     WHERE (ABS(X(1:N)) > SAFMIN )
        X(1:N) = CMPLX(REAL(X(1:N) )/ABS(X(1:N)),AIMAG(X(1:N))/ABS(X(1:N)) )
     ELSEWHERE
        X(1:N) = (1.0E+0,0.0E+0)
     ENDWHERE
     KASE = 2
     JUMP = 4
!
!     ................ ENTRY   (JUMP = 4)
!     X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
!
    CASE (4)
     JLAST = J
     J = ICMAX1( N, X, 1 )
     IF( ( ABS( X( JLAST ) ) /= ABS( X( J ) ) ) .AND. ( ITER < ITMAX ) ) THEN
        ITER = ITER + 1
        X(1:N) = (0.0E+0,0.0E+0)
        X( J ) = (1.0E+0,0.0E+0)
        KASE = 1
        JUMP = 3
     END IF
!
!     ITERATION COMPLETE.  FINAL STAGE.
!
     ALTSGN = 1.0E+0
     DO I = 1, N
        X( I ) = CMPLX( ALTSGN*( 1.0E+0+REAL( I-1 ) / REAL( N-1 ) ) )
        ALTSGN = -ALTSGN
     ENDDO
     KASE = 1
     JUMP = 5
!
!     ................ ENTRY   (JUMP = 5)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
    CASE (5)
     TEMP = 2.0E+0*( SCSUM1( N, X, 1 ) / REAL( 3*N ) )
     IF( TEMP > EST ) THEN
        V(1:N) = X(1:N)
        EST = TEMP
     END IF
     KASE = 0
    CASE DEFAULT
     KASE = 0
   END SELECT
!
   RETURN
!
!     End of CLACON
!
END
