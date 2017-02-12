
/* inclue the blas header. */
#include <vbnmr/blas.h>

/* blas_dnrm2(): compute the euclidean norm of a vector, equivalent to
 * the square root of the sum of its squared elements.
 *
 * operation:
 *  result <- ||x||_2 = sqrt(sum(x.^2))
 *
 * arguments:
 *  @x: input vector to the norm operation.
 *
 * returns:
 *  euclidean norm of the input vector.
 */
inline double blas_dnrm2 (const vector_t *x) {
  /* compute and return the result. */
  return sqrt(blas_ddot(x, x));
}

/* blas_ddot(): compute the euclidean inner product between two vectors,
 *
 * operation:
 *  result <- x' * y = sum(x .* y)
 *
 * arguments:
 *  @x: first input vector.
 *  @y: second input vector.
 *
 * returns:
 *  euclidean inner product (dot product) of the two input vectors.
 */
double blas_ddot (const vector_t *x, const vector_t *y) {
  /* initialize the result. */
  double result = 0.0;

  /* compute the sum of the products of each pair of vector elements. */
  for (unsigned int i = 0; i < x->len; i++)
    result += vector_get(x, i) * vector_get(y, i);

  /* return the result. */
  return result;
}

/* blas_daxpy(): compute the element-wise sum of two vectors.
 *
 * operation:
 *  y <- y + alpha x
 *
 * arguments:
 *  @alpha: scale factor to apply to @x.
 *  @x: input vector to the sum.
 *  @y: input and output vector of the sum.
 */
void blas_daxpy (double alpha, const vector_t *x, vector_t *y) {
  /* compute the sum over all vector elements. */
  for (unsigned int i = 0; i < x->len; i++)
    vector_set(y, i, vector_get(y, i) + alpha * vector_get(x, i));
}

/* --- */

/* blas_dgemv(): compute the linear combination of a vector and the product
 * of a dense matrix with a vector.
 *
 * the user is responsible for ensuring that all matrix and vector operands
 * are non-null and of conformal sizes.
 *
 * operation:
 *  y <- alpha A x  + beta y      if trans == BLAS_NO_TRANS
 *  y <- alpha A' x + beta y      if trans == BLAS_TRANS
 *
 * arguments:
 *  @trans: transposition state of the matrix in the product.
 *  @alpha: scale factor for the matrix-vector product.
 *  @A: input matrix to the product operation.
 *  @x: input vector to the product operation
 *  @beta: scale factor for the output vector.
 *  @y: input and output vector to the combined operation.
 */
void blas_dgemv (blas_transpose_t trans, double alpha, const matrix_t *A,
                 const vector_t *x, double beta, vector_t *y) {
  /* perform: y <- beta y */
  if (beta == 0.0)
    vector_set_zero(y);
  else
    vector_scale(y, beta);

  /* if the scale factor to the matrix-vector portion is zero, return. */
  if (alpha == 0.0)
    return;

  /* perform the dense multiplication operation. */
  if (trans == BLAS_NO_TRANS) {
    /* perform: y <- y + alpha A x */
    for (unsigned int i = 0; i < A->rows; i++) {
      vector_view_t ai = matrix_row(A, i);
      const double ax = blas_ddot(&ai, x);
      vector_set(y, i, vector_get(y, i) + alpha * ax);
    }
  }
  else if (trans == BLAS_TRANS) {
    /* perform: y <- y + alpha A' x */
    for (unsigned int j = 0; j < A->cols; j++) {
      vector_view_t aj = matrix_col(A, j);
      const double ax = blas_ddot(&aj, x);
      vector_set(y, j, vector_get(y, j) + alpha * ax);
    }
  }
}

/* blas_dtrmv(): compute the linear combination of a vector and the product
 * of a lower-triangular matrix with a vector.
 *
 * the user is responsible for ensuring that all matrix and vector operands
 * are non-null and of conformal sizes. only the lower triangle of the
 * input matrix is accessed.
 *
 * operation:
 *  y <- alpha L x  + beta y      if trans == BLAS_NO_TRANS
 *  y <- alpha L' x + beta y      if trans == BLAS_TRANS
 *
 * arguments:
 *  @trans: transposition state of the matrix in the product.
 *  @alpha: scale factor for the matrix-vector product.
 *  @L: input lower-triangular matrix to the product.
 *  @x: input vector to the product operation
 *  @beta: scale factor for the output vector.
 *  @y: input and output vector to the combined operation.
 */
void blas_dtrmv (blas_transpose_t trans, double alpha, const matrix_t *L,
                 const vector_t *x, double beta, vector_t *y) {
  /* perform: y <- beta y */
  if (beta == 0.0)
    vector_set_zero(y);
  else
    vector_scale(y, beta);

  /* if the scale factor to the matrix-vector portion is zero, return. */
  if (alpha == 0.0)
    return;

  /* locally store the problem size. */
  const unsigned int n = x->len;

  /* perform the triangular multiplication operation. */
  if (trans == BLAS_NO_TRANS) {
    /* perform: y <- y + alpha L x */
    for (unsigned int i = 0; i < n; i++) {
      vector_view_t li = matrix_subrow(L, i, 0, i + 1);
      const double lx = blas_ddot(&li, x);
      vector_set(y, i, vector_get(y, i) + alpha * lx);
    }
  }
  else if (trans == BLAS_TRANS) {
    /* perform: y <- y + alpha L' x */
    for (unsigned int j = 0; j < n; j++) {
      vector_view_t lj = matrix_subcol(L, j, j, n - j);
      vector_view_t sj = vector_subvector(x, j, n - j);
      const double lx = blas_ddot(&lj, &sj);
      vector_set(y, j, vector_get(y, j) + alpha * lx);
    }
  }
}

/* blas_dtrsv(): perform forward or backward substitution to solve a system
 * of linear equations involving a triangular matrix.
 *
 * the user is responsible for ensuring that all matrix and vector operands
 * are non-null and of conformal sizes. 
 *
 * operation:
 *  x <- inv(tril(A)) x      if tri == BLAS_LOWER
 *  x <- inv(triu(A)) x      if tri == BLAS_UPPER
 *
 * arguments:
 *  @tri: triangle to access in the matrix in the product.
 *  @A: input triangular matrix.
 *  @x: input and output vector.
 */
void blas_dtrsv (blas_triangle_t tri, const matrix_t *A, vector_t *x) {
  /* locally store the problem size. */
  const unsigned int n = x->len;

  /* perform the substitution operation. */
  if (tri == BLAS_LOWER) {
    /* lower triangle => forward substitution: x <- inv(tril(A)) x */
    for (unsigned int i = 0; i < n; i++) {
      vector_view_t li = matrix_subrow(A, i, 0, i);
      const double Lii = vector_get(&li, i);
      double xi = vector_get(x, i);
      xi -= blas_ddot(&li, x);
      vector_set(x, i, xi / Lii);
    }
  }
  else if (tri == BLAS_UPPER) {
    /* upper triangle => backward substitution: x <- inv(triu(A)) x */
    for (unsigned int i = n - 1; i < n; i--) {
      vector_view_t ui = matrix_subrow(A, i, i + 1, n - i - 1);
      vector_view_t si = vector_subvector(x, i + 1, n - i - 1);
      const double Uii = matrix_get(A, i, i);
      double xi = vector_get(x, i);
      xi -= blas_ddot(&ui, &si);
      vector_set(x, i, xi / Uii);
    }
  }
}

