!> \brief \b CLASSQ updates a sum of squares represented in scaled form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLASSQ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/classq.f90">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/classq.f90">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/classq.f90">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       REAL               SCALE, SUMSQ
!       ..
!       .. Array Arguments ..
!       COMPLEX            X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLASSQ returns the values scale_out and sumsq_out such that
!>
!>    (scale_out**2)*sumsq_out = x( 1 )**2 +...+ x( n )**2 + (scale**2)*sumsq,
!>
!> where x( i ) = X( 1 + ( i - 1 )*INCX ). The value of sumsq is
!> assumed to be non-negative.
!>
!> scale and sumsq must be supplied in SCALE and SUMSQ and
!> scale_out and sumsq_out are overwritten on SCALE and SUMSQ respectively.
!>
!> If scale * sqrt( sumsq ) > tbig then
!>    we require:   scale >= sqrt( TINY*EPS ) / sbig   on entry,
!> and if 0 < scale * sqrt( sumsq ) < tsml then
!>    we require:   scale <= sqrt( HUGE ) / ssml       on entry,
!> where
!>    tbig -- upper threshold for values whose square is representable;
!>    sbig -- scaling constant for big numbers; \see la_constants.f90
!>    tsml -- lower threshold for values whose square is representable;
!>    ssml -- scaling constant for small numbers; \see la_constants.f90
!> and
!>    TINY*EPS -- tiniest representable number;
!>    HUGE     -- biggest representable number.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of elements to be used from the vector x.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX array, dimension (1+(N-1)*abs(INCX))
!>          The vector for which a scaled sum of squares is computed.
!>             x( i ) = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between successive values of the vector x.
!>          If INCX > 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= n
!>          If INCX < 0, X(1-(n-i)*INCX) = x(i) for 1 <= i <= n
!>          If INCX = 0, x isn't a vector so there is no need to call
!>          this subroutine. If you call it anyway, it will count x(1)
!>          in the vector norm N times.
!> \endverbatim
!>
!> \param[in,out] SCALE
!> \verbatim
!>          SCALE is REAL
!>          On entry, the value scale in the equation above.
!>          On exit, SCALE is overwritten by scale_out, the scaling factor
!>          for the sum of squares.
!> \endverbatim
!>
!> \param[in,out] SUMSQ
!> \verbatim
!>          SUMSQ is REAL
!>          On entry, the value sumsq in the equation above.
!>          On exit, SUMSQ is overwritten by sumsq_out, the basic sum of
!>          squares from which scale_out has been factored out.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Edward Anderson, Lockheed Martin
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
!
!> \par Contributors:
!  ==================
!>
!> Weslley Pereira, University of Colorado Denver, USA
!> Nick Papior, Technical University of Denmark, DK
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Anderson E. (2017)
!>  Algorithm 978: Safe Scaling in the Level 1 BLAS
!>  ACM Trans Math Softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!>  Blue, James L. (1978)
!>  A Portable Fortran Program to Find the Euclidean Norm of a Vector
!>  ACM Trans Math Softw 4:15--23
!>  https://doi.org/10.1145/355769.355771
!>
!> \endverbatim
!
!> \ingroup lassq
!
!  =====================================================================
subroutine CLASSQ( n, x, incx, scale, sumsq )
   use LA_CONSTANTS, &
      only: wp=>sp, sbig=>ssbig, ssml=>sssml, tbig=>stbig, tsml=>stsml
   use LA_XISNAN
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!  .. Scalar Arguments ..
   integer :: incx, n
   real(wp) :: scale, sumsq
!  ..
!  .. Array Arguments ..
   complex(wp) :: x(*)
!  ..
!  .. Local Scalars ..
   integer :: i, ix
   logical :: notbig
   real(wp) :: abig, amed, asml, ax, ymax, ymin
!  ..
!
!  Quick return if possible
!
   if( LA_ISNAN(scale) .or. LA_ISNAN(sumsq) ) return
   if( sumsq == 0.0E+0 ) scale = 1.0E+0
   if( scale == 0.0E+0 ) then
      scale = 1.0E+0
      sumsq = 0.0E+0
   end if
   if (n <= 0) return
!
!  Compute the sum of squares in 3 accumulators:
!     abig -- sums of squares scaled down to avoid overflow
!     asml -- sums of squares scaled up to avoid underflow
!     amed -- sums of squares that do not require scaling
!  The thresholds and multipliers are
!     tbig -- values bigger than this are scaled down by sbig
!     tsml -- values smaller than this are scaled up by ssml
!
   notbig = .true.
   asml = 0.0E+0
   amed = 0.0E+0
   abig = 0.0E+0
   ix = 1
   if( incx < 0 ) ix = 1 - (n-1)*incx
   do i = 1, n
      ax = abs(real(x(ix)))
      if (ax > tbig) then
         abig = abig + (ax*sbig)**2
         notbig = .false.
      else if (ax < tsml) then
         if (notbig) asml = asml + (ax*ssml)**2
      else
         amed = amed + ax**2
      end if
      ax = abs(aimag(x(ix)))
      if (ax > tbig) then
         abig = abig + (ax*sbig)**2
         notbig = .false.
      else if (ax < tsml) then
         if (notbig) asml = asml + (ax*ssml)**2
      else
         amed = amed + ax**2
      end if
      ix = ix + incx
   end do
!
!  Put the existing sum of squares into 1.0E+0 of the accumulators
!
   if( sumsq > 0.0E+0 ) then
      ax = scale*sqrt( sumsq )
      if (ax > tbig) then
!        We assume scale >= sqrt( TINY*EPS ) / sbig
         abig = abig + (scale*sbig)**2 * sumsq
      else if (ax < tsml) then
!        We assume scale <= sqrt( HUGE ) / ssml
         if (notbig) asml = asml + (scale*ssml)**2 * sumsq
      else
         amed = amed + scale**2 * sumsq
      end if
   end if
!
!  Combine abig and amed or amed and asml if more than 1.0E+0
!  accumulator was used.
!
   if (abig > 0.0E+0) then
!
!     Combine abig and amed if abig > 0.
!
      if (amed > 0.0E+0 .or. LA_ISNAN(amed)) then
         abig = abig + (amed*sbig)*sbig
      end if
      scale = 1.0E+0 / sbig
      sumsq = abig
   else if (asml > 0.0E+0) then
!
!     Combine amed and asml if asml > 0.
!
      if (amed > 0.0E+0 .or. LA_ISNAN(amed)) then
         amed = sqrt(amed)
         asml = sqrt(asml) / ssml
         if (asml > amed) then
            ymin = amed
            ymax = asml
         else
            ymin = asml
            ymax = amed
         end if
         scale = 1.0E+0
         sumsq = ymax**2*( 1.0E+0 + (ymin/ymax)**2 )
      else
         scale = 1.0E+0 / ssml
         sumsq = asml
      end if
   else
!
!     Otherwise all values are mid-range or 0.0E+0
!
      scale = 1.0E+0
      sumsq = amed
   end if
   return
end subroutine
