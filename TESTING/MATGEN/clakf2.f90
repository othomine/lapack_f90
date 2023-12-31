!> \brief \b CLAKF2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAKF2( M, N, A, LDA, B, D, E, Z, LDZ )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDZ, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), B( LDA, * ), D( LDA, * ),
!      $                   E( LDA, * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Form the 2*M*N by 2*M*N matrix
!>
!>        Z = [ kron(In, A)  -kron(B', Im) ]
!>            [ kron(In, D)  -kron(E', Im) ],
!>
!> where In is the identity matrix of size n and X' is the transpose
!> of X. kron(X, Y) is the Kronecker product between the matrices X
!> and Y.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          Size of matrix, must be >= 1.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          Size of matrix, must be >= 1.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX, dimension ( LDA, M )
!>          The matrix A in the output matrix Z.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A, B, D, and E. ( LDA >= M+N )
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX, dimension ( LDA, N )
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is COMPLEX, dimension ( LDA, M )
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is COMPLEX, dimension ( LDA, N )
!>
!>          The matrices used in forming the output matrix Z.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is COMPLEX, dimension ( LDZ, 2*M*N )
!>          The resultant Kronecker M*N*2 by M*N*2 matrix (see above.)
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of Z. ( LDZ >= 2*M*N )
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
!> \ingroup complex_matgen
!
!  =====================================================================
   SUBROUTINE CLAKF2( M, N, A, LDA, B, D, E, Z, LDZ )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, LDZ, M, N
!     ..
!     .. Array Arguments ..
   COMPLEX            A( LDA, * ), B( LDA, * ), D( LDA, * ), &
                      E( LDA, * ), Z( LDZ, * )
!     ..
!
!  ====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, IK, J, JK, L, MN, MN2
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Subroutines ..
   EXTERNAL           CLASET
!     ..
!     .. Executable Statements ..
!
!     Initialize Z
!
   MN = M*N
   MN2 = 2*MN
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CLASET( 'Full', MN2, MN2, (0.0E+0,0.0E+0), (0.0E+0,0.0E+0), Z, LDZ )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   IK = 1
   DO L = 1, N
!
!        form kron(In, A)
!
      Z( IK:IK-1+M, IK:IK-1+M ) = A(1:M,1:M)
!
!        form kron(In, D)
!
      Z( IK+MN:IK+MN-1+M, IK:IK-1+M ) = D(1:M,1:M)
!
      IK = IK + M
   ENDDO
!
   IK = 1
   DO L = 1, N
      JK = MN + 1
!
      DO J = 1, N
!
!           form -kron(B', Im)
!
         DO I = 1, M
            Z( IK+I-1, JK+I-1 ) = -B( J, L )
         ENDDO
!
!           form -kron(E', Im)
!
         DO I = 1, M
            Z( IK+MN+I-1, JK+I-1 ) = -E( J, L )
         ENDDO
!
         JK = JK + M
      ENDDO
!
      IK = IK + M
   ENDDO
!
   RETURN
!
!     End of CLAKF2
!
END

