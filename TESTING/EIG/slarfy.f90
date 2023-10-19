!> \brief \b SLARFY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLARFY( UPLO, N, V, INCV, TAU, C, LDC, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INCV, LDC, N
!       REAL               TAU
!       ..
!       .. Array Arguments ..
!       REAL               C( LDC, * ), V( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLARFY applies an elementary reflector, or Householder matrix, H,
!> to an n x n symmetric matrix C, from both the left and the right.
!>
!> H is represented in the form
!>
!>    H = I - tau * v * v'
!>
!> where  tau  is a scalar and  v  is a vector.
!>
!> If  tau  is  zero, then  H  is taken to be the unit matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix C is stored.
!>          = 'U':  Upper triangle
!>          = 'L':  Lower triangle
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows and columns of the matrix C.  N >= 0.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is REAL array, dimension
!>                  (1 + (N-1)*abs(INCV))
!>          The vector v as described above.
!> \endverbatim
!>
!> \param[in] INCV
!> \verbatim
!>          INCV is INTEGER
!>          The increment between successive elements of v.  INCV must
!>          not be zero.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is REAL
!>          The value tau as described above.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is REAL array, dimension (LDC, N)
!>          On entry, the matrix C.
!>          On exit, C is overwritten by H * C * H'.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C.  LDC >= max( 1, N ).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (N)
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
   SUBROUTINE SLARFY( UPLO, N, V, INCV, TAU, C, LDC, WORK )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INCV, LDC, N
   REAL               TAU
!     ..
!     .. Array Arguments ..
   REAL               C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   REAL               ALPHA
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Subroutines ..
   EXTERNAL           SAXPY, SSYMV, SSYR2
!     ..
!     .. External Functions ..
   REAL               SDOT
   EXTERNAL           SDOT
!     ..
!     .. Executable Statements ..
!
   IF( TAU == 0.0E+0 ) RETURN
!
!     Form  w:= C * v
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SSYMV( UPLO, N, 1.0E+0, C, LDC, V, INCV, 0.0E+0, WORK, 1 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SSYMV : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   ALPHA = -0.5E+0*TAU*SDOT( N, WORK, 1, V, INCV )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SAXPY( N, ALPHA, V, INCV, WORK, 1 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SAXPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     C := C - v * w' - w * v'
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SSYR2( UPLO, N, -TAU, V, INCV, WORK, 1, C, LDC )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SSYR2 : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   RETURN
!
!     End of SLARFY
!
END



