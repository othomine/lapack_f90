#include <blas_extended.h>
#include <blas_extended_private.h>
#include <blas_fpu.h>
void BLAS_zhemv2_z_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, const void *alpha, const void *a, int lda,
		       const double *x_head, const double *x_tail, int incx,
		       const void *beta, const double *y, int incy,
		       enum blas_prec_type prec)

/* 
 * Purpose
 * =======
 *
 * This routines computes the matrix product:
 *
 *     y  <-  alpha * A * (x_head + x_tail) + beta * y
 * 
 * where A is a complex Hermitian matrix.
 *
 * Arguments
 * =========
 *
 * order   (input) enum blas_order_type
 *         Storage format of input symmetric matrix A.
 * 
 * uplo    (input) enum blas_uplo_type
 *         Determines which half of matrix A (upper or lower triangle)
 *         is accessed.
 *
 * n       (input) int
 *         Dimension of A and size of vectors x, y.
 *
 * alpha   (input) const void*
 * 
 * a       (input) void*
 *         Matrix A.
 *
 * lda     (input) int
 *         Leading dimension of matrix A.
 *
 * x_head  (input) double*
 *         Vector x_head
 *
 * x_tail  (input) double*
 *         Vector x_tail
 *   
 * incx    (input) int
 *         Stride for vector x.
 *
 * beta    (input) const void*
 * 
 * y       (input) double*
 *         Vector y.
 *
 * incy    (input) int
 *         Stride for vector y.
 *
 * prec   (input) enum blas_prec_type
 *        Specifies the internal precision to be used.
 *        = blas_prec_single: single precision.
 *        = blas_prec_double: double precision.
 *        = blas_prec_extra : anything at least 1.5 times as accurate
 *                            than double, and wider than 80-bits.
 *                            We use double-double in our implementation.
 *
 */
{
  /* Routine name */
  const char routine_name[] = "BLAS_zhemv2_z_d_x";
  switch (prec) {

  case blas_prec_single:{

      int i, j;
      int xi, yi, xi0, yi0;
      int aij, ai;
      int incai;
      int incaij, incaij2;

      const double *a_i = (double *) a;
      const double *x_head_i = x_head;
      const double *x_tail_i = x_tail;
      double *y_i = (double *) y;
      double *alpha_i = (double *) alpha;
      double *beta_i = (double *) beta;
      double a_elem[2];
      double x_elem;
      double y_elem[2];
      double diag_elem;
      double prod1[2];
      double prod2[2];
      double sum1[2];
      double sum2[2];
      double tmp1[2];
      double tmp2[2];
      double tmp3[2];



      /* Test for no-op */
      if (n <= 0) {
	return;
      }
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0
	  && (beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
	return;
      }

      /* Check for error conditions. */
      if (n < 0) {
	BLAS_error(routine_name, -3, n, NULL);
      }
      if (lda < n) {
	BLAS_error(routine_name, -6, n, NULL);
      }
      if (incx == 0) {
	BLAS_error(routine_name, -9, incx, NULL);
      }
      if (incy == 0) {
	BLAS_error(routine_name, -12, incy, NULL);
      }

      if ((order == blas_colmajor && uplo == blas_upper) ||
	  (order == blas_rowmajor && uplo == blas_lower)) {
	incai = lda;
	incaij = 1;
	incaij2 = lda;
      } else {
	incai = 1;
	incaij = lda;
	incaij2 = 1;
      }


      incy *= 2;
      incai *= 2;
      incaij *= 2;
      incaij2 *= 2;
      xi0 = (incx > 0) ? 0 : ((-n + 1) * incx);
      yi0 = (incy > 0) ? 0 : ((-n + 1) * incy);



      /* The most general form,   y <--- alpha * A * (x_head + x_tail) + beta * y   */
      if (uplo == blas_lower) {
	for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
	  sum1[0] = sum1[1] = 0.0;
	  sum2[0] = sum2[1] = 0.0;
	  for (j = 0, aij = ai, xi = xi0; j < i;
	       j++, aij += incaij, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    x_elem = x_head_i[xi];
	    {
	      prod1[0] = a_elem[0] * x_elem;
	      prod1[1] = a_elem[1] * x_elem;
	    }
	    sum1[0] = sum1[0] + prod1[0];
	    sum1[1] = sum1[1] + prod1[1];
	    x_elem = x_tail_i[xi];
	    {
	      prod2[0] = a_elem[0] * x_elem;
	      prod2[1] = a_elem[1] * x_elem;
	    }
	    sum2[0] = sum2[0] + prod2[0];
	    sum2[1] = sum2[1] + prod2[1];
	  }

	  diag_elem = a_i[aij];
	  x_elem = x_head_i[xi];
	  prod1[0] = x_elem * diag_elem;
	  prod1[1] = 0.0;
	  sum1[0] = sum1[0] + prod1[0];
	  sum1[1] = sum1[1] + prod1[1];
	  x_elem = x_tail_i[xi];
	  prod2[0] = x_elem * diag_elem;
	  prod2[1] = 0.0;
	  sum2[0] = sum2[0] + prod2[0];
	  sum2[1] = sum2[1] + prod2[1];
	  j++;
	  aij += incaij2;
	  xi += incx;

	  for (; j < n; j++, aij += incaij2, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_elem[1] = -a_elem[1];
	    x_elem = x_head_i[xi];
	    {
	      prod1[0] = a_elem[0] * x_elem;
	      prod1[1] = a_elem[1] * x_elem;
	    }
	    sum1[0] = sum1[0] + prod1[0];
	    sum1[1] = sum1[1] + prod1[1];
	    x_elem = x_tail_i[xi];
	    {
	      prod2[0] = a_elem[0] * x_elem;
	      prod2[1] = a_elem[1] * x_elem;
	    }
	    sum2[0] = sum2[0] + prod2[0];
	    sum2[1] = sum2[1] + prod2[1];
	  }
	  sum1[0] = sum1[0] + sum2[0];
	  sum1[1] = sum1[1] + sum2[1];
	  {
	    tmp1[0] =
	      (double) sum1[0] * alpha_i[0] - (double) sum1[1] * alpha_i[1];
	    tmp1[1] =
	      (double) sum1[0] * alpha_i[1] + (double) sum1[1] * alpha_i[0];
	  }
	  y_elem[0] = y_i[yi];
	  y_elem[1] = y_i[yi + 1];
	  {
	    tmp2[0] =
	      (double) y_elem[0] * beta_i[0] - (double) y_elem[1] * beta_i[1];
	    tmp2[1] =
	      (double) y_elem[0] * beta_i[1] + (double) y_elem[1] * beta_i[0];
	  }
	  tmp3[0] = tmp1[0] + tmp2[0];
	  tmp3[1] = tmp1[1] + tmp2[1];
	  y_i[yi] = tmp3[0];
	  y_i[yi + 1] = tmp3[1];
	}
      } else {
	/* uplo == blas_upper */
	for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
	  sum1[0] = sum1[1] = 0.0;
	  sum2[0] = sum2[1] = 0.0;

	  for (j = 0, aij = ai, xi = xi0; j < i;
	       j++, aij += incaij, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_elem[1] = -a_elem[1];
	    x_elem = x_head_i[xi];
	    {
	      prod1[0] = a_elem[0] * x_elem;
	      prod1[1] = a_elem[1] * x_elem;
	    }
	    sum1[0] = sum1[0] + prod1[0];
	    sum1[1] = sum1[1] + prod1[1];
	    x_elem = x_tail_i[xi];
	    {
	      prod2[0] = a_elem[0] * x_elem;
	      prod2[1] = a_elem[1] * x_elem;
	    }
	    sum2[0] = sum2[0] + prod2[0];
	    sum2[1] = sum2[1] + prod2[1];
	  }

	  diag_elem = a_i[aij];
	  x_elem = x_head_i[xi];
	  prod1[0] = x_elem * diag_elem;
	  prod1[1] = 0.0;
	  sum1[0] = sum1[0] + prod1[0];
	  sum1[1] = sum1[1] + prod1[1];
	  x_elem = x_tail_i[xi];
	  prod2[0] = x_elem * diag_elem;
	  prod2[1] = 0.0;
	  sum2[0] = sum2[0] + prod2[0];
	  sum2[1] = sum2[1] + prod2[1];
	  j++;
	  aij += incaij2;
	  xi += incx;

	  for (; j < n; j++, aij += incaij2, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    x_elem = x_head_i[xi];
	    {
	      prod1[0] = a_elem[0] * x_elem;
	      prod1[1] = a_elem[1] * x_elem;
	    }
	    sum1[0] = sum1[0] + prod1[0];
	    sum1[1] = sum1[1] + prod1[1];
	    x_elem = x_tail_i[xi];
	    {
	      prod2[0] = a_elem[0] * x_elem;
	      prod2[1] = a_elem[1] * x_elem;
	    }
	    sum2[0] = sum2[0] + prod2[0];
	    sum2[1] = sum2[1] + prod2[1];
	  }
	  sum1[0] = sum1[0] + sum2[0];
	  sum1[1] = sum1[1] + sum2[1];
	  {
	    tmp1[0] =
	      (double) sum1[0] * alpha_i[0] - (double) sum1[1] * alpha_i[1];
	    tmp1[1] =
	      (double) sum1[0] * alpha_i[1] + (double) sum1[1] * alpha_i[0];
	  }
	  y_elem[0] = y_i[yi];
	  y_elem[1] = y_i[yi + 1];
	  {
	    tmp2[0] =
	      (double) y_elem[0] * beta_i[0] - (double) y_elem[1] * beta_i[1];
	    tmp2[1] =
	      (double) y_elem[0] * beta_i[1] + (double) y_elem[1] * beta_i[0];
	  }
	  tmp3[0] = tmp1[0] + tmp2[0];
	  tmp3[1] = tmp1[1] + tmp2[1];
	  y_i[yi] = tmp3[0];
	  y_i[yi + 1] = tmp3[1];
	}
      }



      break;
    }

  case blas_prec_double:
  case blas_prec_indigenous:{

      int i, j;
      int xi, yi, xi0, yi0;
      int aij, ai;
      int incai;
      int incaij, incaij2;

      const double *a_i = (double *) a;
      const double *x_head_i = x_head;
      const double *x_tail_i = x_tail;
      double *y_i = (double *) y;
      double *alpha_i = (double *) alpha;
      double *beta_i = (double *) beta;
      double a_elem[2];
      double x_elem;
      double y_elem[2];
      double diag_elem;
      double prod1[2];
      double prod2[2];
      double sum1[2];
      double sum2[2];
      double tmp1[2];
      double tmp2[2];
      double tmp3[2];



      /* Test for no-op */
      if (n <= 0) {
	return;
      }
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0
	  && (beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
	return;
      }

      /* Check for error conditions. */
      if (n < 0) {
	BLAS_error(routine_name, -3, n, NULL);
      }
      if (lda < n) {
	BLAS_error(routine_name, -6, n, NULL);
      }
      if (incx == 0) {
	BLAS_error(routine_name, -9, incx, NULL);
      }
      if (incy == 0) {
	BLAS_error(routine_name, -12, incy, NULL);
      }

      if ((order == blas_colmajor && uplo == blas_upper) ||
	  (order == blas_rowmajor && uplo == blas_lower)) {
	incai = lda;
	incaij = 1;
	incaij2 = lda;
      } else {
	incai = 1;
	incaij = lda;
	incaij2 = 1;
      }


      incy *= 2;
      incai *= 2;
      incaij *= 2;
      incaij2 *= 2;
      xi0 = (incx > 0) ? 0 : ((-n + 1) * incx);
      yi0 = (incy > 0) ? 0 : ((-n + 1) * incy);



      /* The most general form,   y <--- alpha * A * (x_head + x_tail) + beta * y   */
      if (uplo == blas_lower) {
	for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
	  sum1[0] = sum1[1] = 0.0;
	  sum2[0] = sum2[1] = 0.0;
	  for (j = 0, aij = ai, xi = xi0; j < i;
	       j++, aij += incaij, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    x_elem = x_head_i[xi];
	    {
	      prod1[0] = a_elem[0] * x_elem;
	      prod1[1] = a_elem[1] * x_elem;
	    }
	    sum1[0] = sum1[0] + prod1[0];
	    sum1[1] = sum1[1] + prod1[1];
	    x_elem = x_tail_i[xi];
	    {
	      prod2[0] = a_elem[0] * x_elem;
	      prod2[1] = a_elem[1] * x_elem;
	    }
	    sum2[0] = sum2[0] + prod2[0];
	    sum2[1] = sum2[1] + prod2[1];
	  }

	  diag_elem = a_i[aij];
	  x_elem = x_head_i[xi];
	  prod1[0] = x_elem * diag_elem;
	  prod1[1] = 0.0;
	  sum1[0] = sum1[0] + prod1[0];
	  sum1[1] = sum1[1] + prod1[1];
	  x_elem = x_tail_i[xi];
	  prod2[0] = x_elem * diag_elem;
	  prod2[1] = 0.0;
	  sum2[0] = sum2[0] + prod2[0];
	  sum2[1] = sum2[1] + prod2[1];
	  j++;
	  aij += incaij2;
	  xi += incx;

	  for (; j < n; j++, aij += incaij2, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_elem[1] = -a_elem[1];
	    x_elem = x_head_i[xi];
	    {
	      prod1[0] = a_elem[0] * x_elem;
	      prod1[1] = a_elem[1] * x_elem;
	    }
	    sum1[0] = sum1[0] + prod1[0];
	    sum1[1] = sum1[1] + prod1[1];
	    x_elem = x_tail_i[xi];
	    {
	      prod2[0] = a_elem[0] * x_elem;
	      prod2[1] = a_elem[1] * x_elem;
	    }
	    sum2[0] = sum2[0] + prod2[0];
	    sum2[1] = sum2[1] + prod2[1];
	  }
	  sum1[0] = sum1[0] + sum2[0];
	  sum1[1] = sum1[1] + sum2[1];
	  {
	    tmp1[0] =
	      (double) sum1[0] * alpha_i[0] - (double) sum1[1] * alpha_i[1];
	    tmp1[1] =
	      (double) sum1[0] * alpha_i[1] + (double) sum1[1] * alpha_i[0];
	  }
	  y_elem[0] = y_i[yi];
	  y_elem[1] = y_i[yi + 1];
	  {
	    tmp2[0] =
	      (double) y_elem[0] * beta_i[0] - (double) y_elem[1] * beta_i[1];
	    tmp2[1] =
	      (double) y_elem[0] * beta_i[1] + (double) y_elem[1] * beta_i[0];
	  }
	  tmp3[0] = tmp1[0] + tmp2[0];
	  tmp3[1] = tmp1[1] + tmp2[1];
	  y_i[yi] = tmp3[0];
	  y_i[yi + 1] = tmp3[1];
	}
      } else {
	/* uplo == blas_upper */
	for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
	  sum1[0] = sum1[1] = 0.0;
	  sum2[0] = sum2[1] = 0.0;

	  for (j = 0, aij = ai, xi = xi0; j < i;
	       j++, aij += incaij, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_elem[1] = -a_elem[1];
	    x_elem = x_head_i[xi];
	    {
	      prod1[0] = a_elem[0] * x_elem;
	      prod1[1] = a_elem[1] * x_elem;
	    }
	    sum1[0] = sum1[0] + prod1[0];
	    sum1[1] = sum1[1] + prod1[1];
	    x_elem = x_tail_i[xi];
	    {
	      prod2[0] = a_elem[0] * x_elem;
	      prod2[1] = a_elem[1] * x_elem;
	    }
	    sum2[0] = sum2[0] + prod2[0];
	    sum2[1] = sum2[1] + prod2[1];
	  }

	  diag_elem = a_i[aij];
	  x_elem = x_head_i[xi];
	  prod1[0] = x_elem * diag_elem;
	  prod1[1] = 0.0;
	  sum1[0] = sum1[0] + prod1[0];
	  sum1[1] = sum1[1] + prod1[1];
	  x_elem = x_tail_i[xi];
	  prod2[0] = x_elem * diag_elem;
	  prod2[1] = 0.0;
	  sum2[0] = sum2[0] + prod2[0];
	  sum2[1] = sum2[1] + prod2[1];
	  j++;
	  aij += incaij2;
	  xi += incx;

	  for (; j < n; j++, aij += incaij2, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    x_elem = x_head_i[xi];
	    {
	      prod1[0] = a_elem[0] * x_elem;
	      prod1[1] = a_elem[1] * x_elem;
	    }
	    sum1[0] = sum1[0] + prod1[0];
	    sum1[1] = sum1[1] + prod1[1];
	    x_elem = x_tail_i[xi];
	    {
	      prod2[0] = a_elem[0] * x_elem;
	      prod2[1] = a_elem[1] * x_elem;
	    }
	    sum2[0] = sum2[0] + prod2[0];
	    sum2[1] = sum2[1] + prod2[1];
	  }
	  sum1[0] = sum1[0] + sum2[0];
	  sum1[1] = sum1[1] + sum2[1];
	  {
	    tmp1[0] =
	      (double) sum1[0] * alpha_i[0] - (double) sum1[1] * alpha_i[1];
	    tmp1[1] =
	      (double) sum1[0] * alpha_i[1] + (double) sum1[1] * alpha_i[0];
	  }
	  y_elem[0] = y_i[yi];
	  y_elem[1] = y_i[yi + 1];
	  {
	    tmp2[0] =
	      (double) y_elem[0] * beta_i[0] - (double) y_elem[1] * beta_i[1];
	    tmp2[1] =
	      (double) y_elem[0] * beta_i[1] + (double) y_elem[1] * beta_i[0];
	  }
	  tmp3[0] = tmp1[0] + tmp2[0];
	  tmp3[1] = tmp1[1] + tmp2[1];
	  y_i[yi] = tmp3[0];
	  y_i[yi + 1] = tmp3[1];
	}
      }



      break;
    }

  case blas_prec_extra:{

      int i, j;
      int xi, yi, xi0, yi0;
      int aij, ai;
      int incai;
      int incaij, incaij2;

      const double *a_i = (double *) a;
      const double *x_head_i = x_head;
      const double *x_tail_i = x_tail;
      double *y_i = (double *) y;
      double *alpha_i = (double *) alpha;
      double *beta_i = (double *) beta;
      double a_elem[2];
      double x_elem;
      double y_elem[2];
      double diag_elem;
      double head_prod1[2], tail_prod1[2];
      double head_prod2[2], tail_prod2[2];
      double head_sum1[2], tail_sum1[2];
      double head_sum2[2], tail_sum2[2];
      double head_tmp1[2], tail_tmp1[2];
      double head_tmp2[2], tail_tmp2[2];
      double head_tmp3[2], tail_tmp3[2];

      FPU_FIX_DECL;

      /* Test for no-op */
      if (n <= 0) {
	return;
      }
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0
	  && (beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
	return;
      }

      /* Check for error conditions. */
      if (n < 0) {
	BLAS_error(routine_name, -3, n, NULL);
      }
      if (lda < n) {
	BLAS_error(routine_name, -6, n, NULL);
      }
      if (incx == 0) {
	BLAS_error(routine_name, -9, incx, NULL);
      }
      if (incy == 0) {
	BLAS_error(routine_name, -12, incy, NULL);
      }

      if ((order == blas_colmajor && uplo == blas_upper) ||
	  (order == blas_rowmajor && uplo == blas_lower)) {
	incai = lda;
	incaij = 1;
	incaij2 = lda;
      } else {
	incai = 1;
	incaij = lda;
	incaij2 = 1;
      }


      incy *= 2;
      incai *= 2;
      incaij *= 2;
      incaij2 *= 2;
      xi0 = (incx > 0) ? 0 : ((-n + 1) * incx);
      yi0 = (incy > 0) ? 0 : ((-n + 1) * incy);

      FPU_FIX_START;

      /* The most general form,   y <--- alpha * A * (x_head + x_tail) + beta * y   */
      if (uplo == blas_lower) {
	for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
	  head_sum1[0] = head_sum1[1] = tail_sum1[0] = tail_sum1[1] = 0.0;
	  head_sum2[0] = head_sum2[1] = tail_sum2[0] = tail_sum2[1] = 0.0;
	  for (j = 0, aij = ai, xi = xi0; j < i;
	       j++, aij += incaij, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    x_elem = x_head_i[xi];
	    {
	      /* Compute complex-extra = complex-double * real. */
	      double head_t, tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem * split;
		a1 = con - x_elem;
		a1 = con - a1;
		a2 = x_elem - a1;
		con = a_elem[0] * split;
		b1 = con - a_elem[0];
		b1 = con - b1;
		b2 = a_elem[0] - b1;

		head_t = x_elem * a_elem[0];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_prod1[0] = head_t;
	      tail_prod1[0] = tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem * split;
		a1 = con - x_elem;
		a1 = con - a1;
		a2 = x_elem - a1;
		con = a_elem[1] * split;
		b1 = con - a_elem[1];
		b1 = con - b1;
		b2 = a_elem[1] - b1;

		head_t = x_elem * a_elem[1];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_prod1[1] = head_t;
	      tail_prod1[1] = tail_t;
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_sum1[0];
	      tail_a = tail_sum1[0];
	      head_b = head_prod1[0];
	      tail_b = tail_prod1[0];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_sum1[0] = head_t;
	      tail_sum1[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_sum1[1];
	      tail_a = tail_sum1[1];
	      head_b = head_prod1[1];
	      tail_b = tail_prod1[1];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_sum1[1] = head_t;
	      tail_sum1[1] = tail_t;
	    }
	    x_elem = x_tail_i[xi];
	    {
	      /* Compute complex-extra = complex-double * real. */
	      double head_t, tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem * split;
		a1 = con - x_elem;
		a1 = con - a1;
		a2 = x_elem - a1;
		con = a_elem[0] * split;
		b1 = con - a_elem[0];
		b1 = con - b1;
		b2 = a_elem[0] - b1;

		head_t = x_elem * a_elem[0];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_prod2[0] = head_t;
	      tail_prod2[0] = tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem * split;
		a1 = con - x_elem;
		a1 = con - a1;
		a2 = x_elem - a1;
		con = a_elem[1] * split;
		b1 = con - a_elem[1];
		b1 = con - b1;
		b2 = a_elem[1] - b1;

		head_t = x_elem * a_elem[1];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_prod2[1] = head_t;
	      tail_prod2[1] = tail_t;
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_sum2[0];
	      tail_a = tail_sum2[0];
	      head_b = head_prod2[0];
	      tail_b = tail_prod2[0];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_sum2[0] = head_t;
	      tail_sum2[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_sum2[1];
	      tail_a = tail_sum2[1];
	      head_b = head_prod2[1];
	      tail_b = tail_prod2[1];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_sum2[1] = head_t;
	      tail_sum2[1] = tail_t;
	    }
	  }

	  diag_elem = a_i[aij];
	  x_elem = x_head_i[xi];
	  {
	    /* Compute complex double_double = double * double. */
	    double a1, a2, b1, b2, con;

	    con = x_elem * split;
	    a1 = con - x_elem;
	    a1 = con - a1;
	    a2 = x_elem - a1;
	    con = diag_elem * split;
	    b1 = con - diag_elem;
	    b1 = con - b1;
	    b2 = diag_elem - b1;

	    head_prod1[0] = x_elem * diag_elem;
	    tail_prod1[0] =
	      (((a1 * b1 - head_prod1[0]) + a1 * b2) + a2 * b1) + a2 * b2;
	    head_prod1[1] = tail_prod1[1] = 0.0;
	  }
	  {
	    double head_t, tail_t;
	    double head_a, tail_a;
	    double head_b, tail_b;
	    /* Real part */
	    head_a = head_sum1[0];
	    tail_a = tail_sum1[0];
	    head_b = head_prod1[0];
	    tail_b = tail_prod1[0];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_sum1[0] = head_t;
	    tail_sum1[0] = tail_t;
	    /* Imaginary part */
	    head_a = head_sum1[1];
	    tail_a = tail_sum1[1];
	    head_b = head_prod1[1];
	    tail_b = tail_prod1[1];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_sum1[1] = head_t;
	    tail_sum1[1] = tail_t;
	  }
	  x_elem = x_tail_i[xi];
	  {
	    /* Compute complex double_double = double * double. */
	    double a1, a2, b1, b2, con;

	    con = x_elem * split;
	    a1 = con - x_elem;
	    a1 = con - a1;
	    a2 = x_elem - a1;
	    con = diag_elem * split;
	    b1 = con - diag_elem;
	    b1 = con - b1;
	    b2 = diag_elem - b1;

	    head_prod2[0] = x_elem * diag_elem;
	    tail_prod2[0] =
	      (((a1 * b1 - head_prod2[0]) + a1 * b2) + a2 * b1) + a2 * b2;
	    head_prod2[1] = tail_prod2[1] = 0.0;
	  }
	  {
	    double head_t, tail_t;
	    double head_a, tail_a;
	    double head_b, tail_b;
	    /* Real part */
	    head_a = head_sum2[0];
	    tail_a = tail_sum2[0];
	    head_b = head_prod2[0];
	    tail_b = tail_prod2[0];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_sum2[0] = head_t;
	    tail_sum2[0] = tail_t;
	    /* Imaginary part */
	    head_a = head_sum2[1];
	    tail_a = tail_sum2[1];
	    head_b = head_prod2[1];
	    tail_b = tail_prod2[1];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_sum2[1] = head_t;
	    tail_sum2[1] = tail_t;
	  }
	  j++;
	  aij += incaij2;
	  xi += incx;

	  for (; j < n; j++, aij += incaij2, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_elem[1] = -a_elem[1];
	    x_elem = x_head_i[xi];
	    {
	      /* Compute complex-extra = complex-double * real. */
	      double head_t, tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem * split;
		a1 = con - x_elem;
		a1 = con - a1;
		a2 = x_elem - a1;
		con = a_elem[0] * split;
		b1 = con - a_elem[0];
		b1 = con - b1;
		b2 = a_elem[0] - b1;

		head_t = x_elem * a_elem[0];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_prod1[0] = head_t;
	      tail_prod1[0] = tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem * split;
		a1 = con - x_elem;
		a1 = con - a1;
		a2 = x_elem - a1;
		con = a_elem[1] * split;
		b1 = con - a_elem[1];
		b1 = con - b1;
		b2 = a_elem[1] - b1;

		head_t = x_elem * a_elem[1];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_prod1[1] = head_t;
	      tail_prod1[1] = tail_t;
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_sum1[0];
	      tail_a = tail_sum1[0];
	      head_b = head_prod1[0];
	      tail_b = tail_prod1[0];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_sum1[0] = head_t;
	      tail_sum1[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_sum1[1];
	      tail_a = tail_sum1[1];
	      head_b = head_prod1[1];
	      tail_b = tail_prod1[1];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_sum1[1] = head_t;
	      tail_sum1[1] = tail_t;
	    }
	    x_elem = x_tail_i[xi];
	    {
	      /* Compute complex-extra = complex-double * real. */
	      double head_t, tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem * split;
		a1 = con - x_elem;
		a1 = con - a1;
		a2 = x_elem - a1;
		con = a_elem[0] * split;
		b1 = con - a_elem[0];
		b1 = con - b1;
		b2 = a_elem[0] - b1;

		head_t = x_elem * a_elem[0];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_prod2[0] = head_t;
	      tail_prod2[0] = tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem * split;
		a1 = con - x_elem;
		a1 = con - a1;
		a2 = x_elem - a1;
		con = a_elem[1] * split;
		b1 = con - a_elem[1];
		b1 = con - b1;
		b2 = a_elem[1] - b1;

		head_t = x_elem * a_elem[1];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_prod2[1] = head_t;
	      tail_prod2[1] = tail_t;
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_sum2[0];
	      tail_a = tail_sum2[0];
	      head_b = head_prod2[0];
	      tail_b = tail_prod2[0];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_sum2[0] = head_t;
	      tail_sum2[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_sum2[1];
	      tail_a = tail_sum2[1];
	      head_b = head_prod2[1];
	      tail_b = tail_prod2[1];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_sum2[1] = head_t;
	      tail_sum2[1] = tail_t;
	    }
	  }
	  {
	    double head_t, tail_t;
	    double head_a, tail_a;
	    double head_b, tail_b;
	    /* Real part */
	    head_a = head_sum1[0];
	    tail_a = tail_sum1[0];
	    head_b = head_sum2[0];
	    tail_b = tail_sum2[0];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_sum1[0] = head_t;
	    tail_sum1[0] = tail_t;
	    /* Imaginary part */
	    head_a = head_sum1[1];
	    tail_a = tail_sum1[1];
	    head_b = head_sum2[1];
	    tail_b = tail_sum2[1];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_sum1[1] = head_t;
	    tail_sum1[1] = tail_t;
	  }
	  {
	    /* Compute complex-extra = complex-extra * complex-double. */
	    double head_a0, tail_a0;
	    double head_a1, tail_a1;
	    double head_t1, tail_t1;
	    double head_t2, tail_t2;
	    head_a0 = head_sum1[0];
	    tail_a0 = tail_sum1[0];
	    head_a1 = head_sum1[1];
	    tail_a1 = tail_sum1[1];
	    /* real part */
	    {
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_a0 * split;
	      a11 = con - head_a0;
	      a11 = con - a11;
	      a21 = head_a0 - a11;
	      con = alpha_i[0] * split;
	      b1 = con - alpha_i[0];
	      b1 = con - b1;
	      b2 = alpha_i[0] - b1;

	      c11 = head_a0 * alpha_i[0];
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_a0 * alpha_i[0];
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    {
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_a1 * split;
	      a11 = con - head_a1;
	      a11 = con - a11;
	      a21 = head_a1 - a11;
	      con = alpha_i[1] * split;
	      b1 = con - alpha_i[1];
	      b1 = con - b1;
	      b2 = alpha_i[1] - b1;

	      c11 = head_a1 * alpha_i[1];
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_a1 * alpha_i[1];
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_t2 = t1 + t2;
	      tail_t2 = t2 - (head_t2 - t1);
	    }
	    head_t2 = -head_t2;
	    tail_t2 = -tail_t2;
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_t1 + head_t2;
	      bv = s1 - head_t1;
	      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_t1 + tail_t2;
	      bv = t1 - tail_t1;
	      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    head_tmp1[0] = head_t1;
	    tail_tmp1[0] = tail_t1;
	    /* imaginary part */
	    {
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_a1 * split;
	      a11 = con - head_a1;
	      a11 = con - a11;
	      a21 = head_a1 - a11;
	      con = alpha_i[0] * split;
	      b1 = con - alpha_i[0];
	      b1 = con - b1;
	      b2 = alpha_i[0] - b1;

	      c11 = head_a1 * alpha_i[0];
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_a1 * alpha_i[0];
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    {
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_a0 * split;
	      a11 = con - head_a0;
	      a11 = con - a11;
	      a21 = head_a0 - a11;
	      con = alpha_i[1] * split;
	      b1 = con - alpha_i[1];
	      b1 = con - b1;
	      b2 = alpha_i[1] - b1;

	      c11 = head_a0 * alpha_i[1];
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_a0 * alpha_i[1];
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_t2 = t1 + t2;
	      tail_t2 = t2 - (head_t2 - t1);
	    }
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_t1 + head_t2;
	      bv = s1 - head_t1;
	      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_t1 + tail_t2;
	      bv = t1 - tail_t1;
	      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    head_tmp1[1] = head_t1;
	    tail_tmp1[1] = tail_t1;
	  }

	  y_elem[0] = y_i[yi];
	  y_elem[1] = y_i[yi + 1];
	  {
	    /* Compute complex-extra = complex-double * complex-double. */
	    double head_t1, tail_t1;
	    double head_t2, tail_t2;
	    /* Real part */
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = y_elem[0] * split;
	      a1 = con - y_elem[0];
	      a1 = con - a1;
	      a2 = y_elem[0] - a1;
	      con = beta_i[0] * split;
	      b1 = con - beta_i[0];
	      b1 = con - b1;
	      b2 = beta_i[0] - b1;

	      head_t1 = y_elem[0] * beta_i[0];
	      tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = y_elem[1] * split;
	      a1 = con - y_elem[1];
	      a1 = con - a1;
	      a2 = y_elem[1] - a1;
	      con = beta_i[1] * split;
	      b1 = con - beta_i[1];
	      b1 = con - b1;
	      b2 = beta_i[1] - b1;

	      head_t2 = y_elem[1] * beta_i[1];
	      tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    head_t2 = -head_t2;
	    tail_t2 = -tail_t2;
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_t1 + head_t2;
	      bv = s1 - head_t1;
	      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_t1 + tail_t2;
	      bv = t1 - tail_t1;
	      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    head_tmp2[0] = head_t1;
	    tail_tmp2[0] = tail_t1;
	    /* Imaginary part */
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = y_elem[1] * split;
	      a1 = con - y_elem[1];
	      a1 = con - a1;
	      a2 = y_elem[1] - a1;
	      con = beta_i[0] * split;
	      b1 = con - beta_i[0];
	      b1 = con - b1;
	      b2 = beta_i[0] - b1;

	      head_t1 = y_elem[1] * beta_i[0];
	      tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = y_elem[0] * split;
	      a1 = con - y_elem[0];
	      a1 = con - a1;
	      a2 = y_elem[0] - a1;
	      con = beta_i[1] * split;
	      b1 = con - beta_i[1];
	      b1 = con - b1;
	      b2 = beta_i[1] - b1;

	      head_t2 = y_elem[0] * beta_i[1];
	      tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_t1 + head_t2;
	      bv = s1 - head_t1;
	      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_t1 + tail_t2;
	      bv = t1 - tail_t1;
	      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    head_tmp2[1] = head_t1;
	    tail_tmp2[1] = tail_t1;
	  }
	  {
	    double head_t, tail_t;
	    double head_a, tail_a;
	    double head_b, tail_b;
	    /* Real part */
	    head_a = head_tmp1[0];
	    tail_a = tail_tmp1[0];
	    head_b = head_tmp2[0];
	    tail_b = tail_tmp2[0];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_tmp3[0] = head_t;
	    tail_tmp3[0] = tail_t;
	    /* Imaginary part */
	    head_a = head_tmp1[1];
	    tail_a = tail_tmp1[1];
	    head_b = head_tmp2[1];
	    tail_b = tail_tmp2[1];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_tmp3[1] = head_t;
	    tail_tmp3[1] = tail_t;
	  }
	  y_i[yi] = head_tmp3[0];
	  y_i[yi + 1] = head_tmp3[1];
	}
      } else {
	/* uplo == blas_upper */
	for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
	  head_sum1[0] = head_sum1[1] = tail_sum1[0] = tail_sum1[1] = 0.0;
	  head_sum2[0] = head_sum2[1] = tail_sum2[0] = tail_sum2[1] = 0.0;

	  for (j = 0, aij = ai, xi = xi0; j < i;
	       j++, aij += incaij, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_elem[1] = -a_elem[1];
	    x_elem = x_head_i[xi];
	    {
	      /* Compute complex-extra = complex-double * real. */
	      double head_t, tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem * split;
		a1 = con - x_elem;
		a1 = con - a1;
		a2 = x_elem - a1;
		con = a_elem[0] * split;
		b1 = con - a_elem[0];
		b1 = con - b1;
		b2 = a_elem[0] - b1;

		head_t = x_elem * a_elem[0];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_prod1[0] = head_t;
	      tail_prod1[0] = tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem * split;
		a1 = con - x_elem;
		a1 = con - a1;
		a2 = x_elem - a1;
		con = a_elem[1] * split;
		b1 = con - a_elem[1];
		b1 = con - b1;
		b2 = a_elem[1] - b1;

		head_t = x_elem * a_elem[1];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_prod1[1] = head_t;
	      tail_prod1[1] = tail_t;
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_sum1[0];
	      tail_a = tail_sum1[0];
	      head_b = head_prod1[0];
	      tail_b = tail_prod1[0];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_sum1[0] = head_t;
	      tail_sum1[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_sum1[1];
	      tail_a = tail_sum1[1];
	      head_b = head_prod1[1];
	      tail_b = tail_prod1[1];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_sum1[1] = head_t;
	      tail_sum1[1] = tail_t;
	    }
	    x_elem = x_tail_i[xi];
	    {
	      /* Compute complex-extra = complex-double * real. */
	      double head_t, tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem * split;
		a1 = con - x_elem;
		a1 = con - a1;
		a2 = x_elem - a1;
		con = a_elem[0] * split;
		b1 = con - a_elem[0];
		b1 = con - b1;
		b2 = a_elem[0] - b1;

		head_t = x_elem * a_elem[0];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_prod2[0] = head_t;
	      tail_prod2[0] = tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem * split;
		a1 = con - x_elem;
		a1 = con - a1;
		a2 = x_elem - a1;
		con = a_elem[1] * split;
		b1 = con - a_elem[1];
		b1 = con - b1;
		b2 = a_elem[1] - b1;

		head_t = x_elem * a_elem[1];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_prod2[1] = head_t;
	      tail_prod2[1] = tail_t;
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_sum2[0];
	      tail_a = tail_sum2[0];
	      head_b = head_prod2[0];
	      tail_b = tail_prod2[0];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_sum2[0] = head_t;
	      tail_sum2[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_sum2[1];
	      tail_a = tail_sum2[1];
	      head_b = head_prod2[1];
	      tail_b = tail_prod2[1];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_sum2[1] = head_t;
	      tail_sum2[1] = tail_t;
	    }
	  }

	  diag_elem = a_i[aij];
	  x_elem = x_head_i[xi];
	  {
	    /* Compute complex double_double = double * double. */
	    double a1, a2, b1, b2, con;

	    con = x_elem * split;
	    a1 = con - x_elem;
	    a1 = con - a1;
	    a2 = x_elem - a1;
	    con = diag_elem * split;
	    b1 = con - diag_elem;
	    b1 = con - b1;
	    b2 = diag_elem - b1;

	    head_prod1[0] = x_elem * diag_elem;
	    tail_prod1[0] =
	      (((a1 * b1 - head_prod1[0]) + a1 * b2) + a2 * b1) + a2 * b2;
	    head_prod1[1] = tail_prod1[1] = 0.0;
	  }
	  {
	    double head_t, tail_t;
	    double head_a, tail_a;
	    double head_b, tail_b;
	    /* Real part */
	    head_a = head_sum1[0];
	    tail_a = tail_sum1[0];
	    head_b = head_prod1[0];
	    tail_b = tail_prod1[0];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_sum1[0] = head_t;
	    tail_sum1[0] = tail_t;
	    /* Imaginary part */
	    head_a = head_sum1[1];
	    tail_a = tail_sum1[1];
	    head_b = head_prod1[1];
	    tail_b = tail_prod1[1];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_sum1[1] = head_t;
	    tail_sum1[1] = tail_t;
	  }
	  x_elem = x_tail_i[xi];
	  {
	    /* Compute complex double_double = double * double. */
	    double a1, a2, b1, b2, con;

	    con = x_elem * split;
	    a1 = con - x_elem;
	    a1 = con - a1;
	    a2 = x_elem - a1;
	    con = diag_elem * split;
	    b1 = con - diag_elem;
	    b1 = con - b1;
	    b2 = diag_elem - b1;

	    head_prod2[0] = x_elem * diag_elem;
	    tail_prod2[0] =
	      (((a1 * b1 - head_prod2[0]) + a1 * b2) + a2 * b1) + a2 * b2;
	    head_prod2[1] = tail_prod2[1] = 0.0;
	  }
	  {
	    double head_t, tail_t;
	    double head_a, tail_a;
	    double head_b, tail_b;
	    /* Real part */
	    head_a = head_sum2[0];
	    tail_a = tail_sum2[0];
	    head_b = head_prod2[0];
	    tail_b = tail_prod2[0];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_sum2[0] = head_t;
	    tail_sum2[0] = tail_t;
	    /* Imaginary part */
	    head_a = head_sum2[1];
	    tail_a = tail_sum2[1];
	    head_b = head_prod2[1];
	    tail_b = tail_prod2[1];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_sum2[1] = head_t;
	    tail_sum2[1] = tail_t;
	  }
	  j++;
	  aij += incaij2;
	  xi += incx;

	  for (; j < n; j++, aij += incaij2, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    x_elem = x_head_i[xi];
	    {
	      /* Compute complex-extra = complex-double * real. */
	      double head_t, tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem * split;
		a1 = con - x_elem;
		a1 = con - a1;
		a2 = x_elem - a1;
		con = a_elem[0] * split;
		b1 = con - a_elem[0];
		b1 = con - b1;
		b2 = a_elem[0] - b1;

		head_t = x_elem * a_elem[0];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_prod1[0] = head_t;
	      tail_prod1[0] = tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem * split;
		a1 = con - x_elem;
		a1 = con - a1;
		a2 = x_elem - a1;
		con = a_elem[1] * split;
		b1 = con - a_elem[1];
		b1 = con - b1;
		b2 = a_elem[1] - b1;

		head_t = x_elem * a_elem[1];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_prod1[1] = head_t;
	      tail_prod1[1] = tail_t;
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_sum1[0];
	      tail_a = tail_sum1[0];
	      head_b = head_prod1[0];
	      tail_b = tail_prod1[0];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_sum1[0] = head_t;
	      tail_sum1[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_sum1[1];
	      tail_a = tail_sum1[1];
	      head_b = head_prod1[1];
	      tail_b = tail_prod1[1];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_sum1[1] = head_t;
	      tail_sum1[1] = tail_t;
	    }
	    x_elem = x_tail_i[xi];
	    {
	      /* Compute complex-extra = complex-double * real. */
	      double head_t, tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem * split;
		a1 = con - x_elem;
		a1 = con - a1;
		a2 = x_elem - a1;
		con = a_elem[0] * split;
		b1 = con - a_elem[0];
		b1 = con - b1;
		b2 = a_elem[0] - b1;

		head_t = x_elem * a_elem[0];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_prod2[0] = head_t;
	      tail_prod2[0] = tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem * split;
		a1 = con - x_elem;
		a1 = con - a1;
		a2 = x_elem - a1;
		con = a_elem[1] * split;
		b1 = con - a_elem[1];
		b1 = con - b1;
		b2 = a_elem[1] - b1;

		head_t = x_elem * a_elem[1];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_prod2[1] = head_t;
	      tail_prod2[1] = tail_t;
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_sum2[0];
	      tail_a = tail_sum2[0];
	      head_b = head_prod2[0];
	      tail_b = tail_prod2[0];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_sum2[0] = head_t;
	      tail_sum2[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_sum2[1];
	      tail_a = tail_sum2[1];
	      head_b = head_prod2[1];
	      tail_b = tail_prod2[1];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_sum2[1] = head_t;
	      tail_sum2[1] = tail_t;
	    }
	  }
	  {
	    double head_t, tail_t;
	    double head_a, tail_a;
	    double head_b, tail_b;
	    /* Real part */
	    head_a = head_sum1[0];
	    tail_a = tail_sum1[0];
	    head_b = head_sum2[0];
	    tail_b = tail_sum2[0];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_sum1[0] = head_t;
	    tail_sum1[0] = tail_t;
	    /* Imaginary part */
	    head_a = head_sum1[1];
	    tail_a = tail_sum1[1];
	    head_b = head_sum2[1];
	    tail_b = tail_sum2[1];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_sum1[1] = head_t;
	    tail_sum1[1] = tail_t;
	  }
	  {
	    /* Compute complex-extra = complex-extra * complex-double. */
	    double head_a0, tail_a0;
	    double head_a1, tail_a1;
	    double head_t1, tail_t1;
	    double head_t2, tail_t2;
	    head_a0 = head_sum1[0];
	    tail_a0 = tail_sum1[0];
	    head_a1 = head_sum1[1];
	    tail_a1 = tail_sum1[1];
	    /* real part */
	    {
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_a0 * split;
	      a11 = con - head_a0;
	      a11 = con - a11;
	      a21 = head_a0 - a11;
	      con = alpha_i[0] * split;
	      b1 = con - alpha_i[0];
	      b1 = con - b1;
	      b2 = alpha_i[0] - b1;

	      c11 = head_a0 * alpha_i[0];
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_a0 * alpha_i[0];
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    {
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_a1 * split;
	      a11 = con - head_a1;
	      a11 = con - a11;
	      a21 = head_a1 - a11;
	      con = alpha_i[1] * split;
	      b1 = con - alpha_i[1];
	      b1 = con - b1;
	      b2 = alpha_i[1] - b1;

	      c11 = head_a1 * alpha_i[1];
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_a1 * alpha_i[1];
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_t2 = t1 + t2;
	      tail_t2 = t2 - (head_t2 - t1);
	    }
	    head_t2 = -head_t2;
	    tail_t2 = -tail_t2;
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_t1 + head_t2;
	      bv = s1 - head_t1;
	      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_t1 + tail_t2;
	      bv = t1 - tail_t1;
	      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    head_tmp1[0] = head_t1;
	    tail_tmp1[0] = tail_t1;
	    /* imaginary part */
	    {
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_a1 * split;
	      a11 = con - head_a1;
	      a11 = con - a11;
	      a21 = head_a1 - a11;
	      con = alpha_i[0] * split;
	      b1 = con - alpha_i[0];
	      b1 = con - b1;
	      b2 = alpha_i[0] - b1;

	      c11 = head_a1 * alpha_i[0];
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_a1 * alpha_i[0];
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    {
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_a0 * split;
	      a11 = con - head_a0;
	      a11 = con - a11;
	      a21 = head_a0 - a11;
	      con = alpha_i[1] * split;
	      b1 = con - alpha_i[1];
	      b1 = con - b1;
	      b2 = alpha_i[1] - b1;

	      c11 = head_a0 * alpha_i[1];
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_a0 * alpha_i[1];
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_t2 = t1 + t2;
	      tail_t2 = t2 - (head_t2 - t1);
	    }
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_t1 + head_t2;
	      bv = s1 - head_t1;
	      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_t1 + tail_t2;
	      bv = t1 - tail_t1;
	      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    head_tmp1[1] = head_t1;
	    tail_tmp1[1] = tail_t1;
	  }

	  y_elem[0] = y_i[yi];
	  y_elem[1] = y_i[yi + 1];
	  {
	    /* Compute complex-extra = complex-double * complex-double. */
	    double head_t1, tail_t1;
	    double head_t2, tail_t2;
	    /* Real part */
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = y_elem[0] * split;
	      a1 = con - y_elem[0];
	      a1 = con - a1;
	      a2 = y_elem[0] - a1;
	      con = beta_i[0] * split;
	      b1 = con - beta_i[0];
	      b1 = con - b1;
	      b2 = beta_i[0] - b1;

	      head_t1 = y_elem[0] * beta_i[0];
	      tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = y_elem[1] * split;
	      a1 = con - y_elem[1];
	      a1 = con - a1;
	      a2 = y_elem[1] - a1;
	      con = beta_i[1] * split;
	      b1 = con - beta_i[1];
	      b1 = con - b1;
	      b2 = beta_i[1] - b1;

	      head_t2 = y_elem[1] * beta_i[1];
	      tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    head_t2 = -head_t2;
	    tail_t2 = -tail_t2;
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_t1 + head_t2;
	      bv = s1 - head_t1;
	      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_t1 + tail_t2;
	      bv = t1 - tail_t1;
	      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    head_tmp2[0] = head_t1;
	    tail_tmp2[0] = tail_t1;
	    /* Imaginary part */
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = y_elem[1] * split;
	      a1 = con - y_elem[1];
	      a1 = con - a1;
	      a2 = y_elem[1] - a1;
	      con = beta_i[0] * split;
	      b1 = con - beta_i[0];
	      b1 = con - b1;
	      b2 = beta_i[0] - b1;

	      head_t1 = y_elem[1] * beta_i[0];
	      tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = y_elem[0] * split;
	      a1 = con - y_elem[0];
	      a1 = con - a1;
	      a2 = y_elem[0] - a1;
	      con = beta_i[1] * split;
	      b1 = con - beta_i[1];
	      b1 = con - b1;
	      b2 = beta_i[1] - b1;

	      head_t2 = y_elem[0] * beta_i[1];
	      tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_t1 + head_t2;
	      bv = s1 - head_t1;
	      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_t1 + tail_t2;
	      bv = t1 - tail_t1;
	      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    head_tmp2[1] = head_t1;
	    tail_tmp2[1] = tail_t1;
	  }
	  {
	    double head_t, tail_t;
	    double head_a, tail_a;
	    double head_b, tail_b;
	    /* Real part */
	    head_a = head_tmp1[0];
	    tail_a = tail_tmp1[0];
	    head_b = head_tmp2[0];
	    tail_b = tail_tmp2[0];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_tmp3[0] = head_t;
	    tail_tmp3[0] = tail_t;
	    /* Imaginary part */
	    head_a = head_tmp1[1];
	    tail_a = tail_tmp1[1];
	    head_b = head_tmp2[1];
	    tail_b = tail_tmp2[1];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_tmp3[1] = head_t;
	    tail_tmp3[1] = tail_t;
	  }
	  y_i[yi] = head_tmp3[0];
	  y_i[yi + 1] = head_tmp3[1];
	}
      }

      FPU_FIX_STOP;

      break;
    }
  }
}				/* end BLAS_zhemv2_z_d_x */
