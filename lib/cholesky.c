
/* include the cholesky decomposition header. */
#include <vbnmr/cholesky.h>

/* cholesky_decomp(): compute the cholesky decomposition of a symmetric
 * positive definite matrix.
 *
 * arguments:
 *  @A: matrix to decompose in place.
 *
 * returns:
 *  integer indicating success (1) or failure (0). the method will return
 *  failure if the input matrix is not positive definite.
 */
int cholesky_decomp (matrix_t *A) {
  /* locally store the matrix size. */
  const unsigned int n = A->cols;

  /* perform decomposition column-wise. */
  for (unsigned int j = 0; j < n; j++) {
    /* v := A(j : n, j) */
    vector_view_t v = matrix_subcol(A, j, j, n - j);

    if (j > 0) {
      /* w := A(j, 1 : j)
       * M := A(j : n, 1 : j)
       */
      vector_view_t w = matrix_subrow(A, j, 0, j);
      matrix_view_t M = matrix_submatrix(A, j, 0, n - j, j);

      /* v <- v - M w */
      blas_dgemv(BLAS_NO_TRANS, -1.0, &M, &w, 1.0, &v);
    }

    /* ensure the positive-definiteness is satisfied. */
    const double Ajj = matrix_get(A, j, j);
    if (Ajj <= 0.0)
      return 0;

    /* v <- v ./ sqrt(A(j,j)) */
    vector_scale(&v, 1.0 / sqrt(Ajj));
  }

  /* symmetrize the cholesky factor matrix. */
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = i + 1; j < n; j++)
      matrix_set(A, i, j, matrix_get(A, j, i));

  /* return success. */
  return 1;
}

/* cholesky_invert(): compute the inverse of a symmetric positive definite
 * matrix from its cholesky factorization.
 *
 * arguments:
 *  @L: cholesky factorization to use for inversion.
 *  @B: matrix to store the inverse into.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int cholesky_invert (const matrix_t *L, matrix_t *B) {
  /* initialize the matrix inverse. */
  matrix_set_ident(B);

  /* invert each column of the identity matrix. */
  for (unsigned int j = 0; j < B->cols; j++) {
    /* extract the current column. */
    vector_view_t b = matrix_col(B, j);

    /* perform forward and backward substitution. */
    blas_dtrsv(BLAS_LOWER, L, &b);
    blas_dtrsv(BLAS_UPPER, L, &b);
  }

  /* return success. */
  return 1;
}

/* cholesky_solve(): solve a system of equations involving a symmetric
 * positive definite matrix using its precomputed cholesky factorization.
 *
 * arguments:
 *  @L: cholesky factorization to use for solution.
 *  @b: vector of input variables.
 *  @x: vector of output regressors.
 */
void cholesky_solve (const matrix_t *L, const vector_t *b, vector_t *x) {
  /* solve using forward and backward substitution. */
  vector_copy(x, b);
  blas_dtrsv(BLAS_LOWER, L, x);
  blas_dtrsv(BLAS_UPPER, L, x);
}

/* cholesky_update(): update the cholesky decomposition of a symmetric
 * positive definite matrix to reflect the application of a symmetric
 * rank-one update.
 *
 * arguments:
 *  @L: cholesky factors of the matrix.
 *  @x: update vector to apply.
 */
void cholesky_update (matrix_t *L, vector_t *x) {
  /* locally store the update vector length. */
  const unsigned int n = x->len;

  /* perform updating column-wise. */
  for (unsigned int k = 0; k < n; k++) {
    /* get the relevant quantities. */
    const double Lkk = matrix_get(L, k, k);
    const double xk = vector_get(x, k);

    /* compute scale factors for the column update. */
    const double r = sqrt(Lkk * Lkk + xk * xk);
    const double s = xk / Lkk;
    const double c = r / Lkk;

    /* update the matrix diagonal. */
    matrix_set(L, k, k, r);

    /* lk := L(k+1 : n, k)
     * yk := x(k+1 : n)
     */
    vector_view_t lk = matrix_subcol(L, k, k + 1, n - k - 1);
    vector_view_t yk = vector_subvector(x, k + 1, n - k - 1);

    /* lk <- (lk + s yk) / c */
    blas_daxpy(s, &yk, &lk);
    vector_scale(&lk, 1.0 / c);

    /* yk <- c yk - s lk */
    blas_daxpy(-s / c, &lk, &yk);
    vector_scale(&yk, c);
  }

  /* re-symmetrize the cholesky factor matrix. */
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = i + 1; j < n; j++)
      matrix_set(L, i, j, matrix_get(L, j, i));
}

/* cholesky_downdate(): update the cholesky decomposition of a symmetric
 * positive definite matrix to reflect the application of a symmetric
 * rank-one downdate.
 *
 * arguments:
 *  @L: cholesky factors of the matrix.
 *  @y: downdate vector to apply.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the rank-one downdate
 *  preserved positive definiteness.
 */
int cholesky_downdate (matrix_t *L, vector_t *y) {
  /* locally store the downdate vector length. */
  const unsigned int n = y->len;

  /* perform downdating column-wise. */
  for (unsigned int k = 0; k < n; k++) {
    /* get the relevant quantities. */
    const double Lkk = matrix_get(L, k, k);
    const double yk = vector_get(y, k);

    /* check that positive-definiteness is preserved. */
    const double r2 = Lkk * Lkk - yk * yk;
    if (r2 <= 0.0)
      return 0;

    /* compute scale factors for the column update. */
    const double r = sqrt(r2);
    const double s = yk / Lkk;
    const double c = r / Lkk;

    /* update the matrix diagonal. */
    matrix_set(L, k, k, r);

    /* lk := L(k+1 : n, k)
     * zk := y(k+1 : n)
     */
    vector_view_t lk = matrix_subcol(L, k, k + 1, n - k - 1);
    vector_view_t zk = vector_subvector(y, k + 1, n - k - 1);

    /* lk <- (lk - s zk) / c */
    blas_daxpy(-s, &zk, &lk);
    vector_scale(&lk, 1.0 / c);

    /* zk <- c zk - s lk */
    blas_daxpy(-s / c, &lk, &zk);
    vector_scale(&zk, c);
  }

  /* re-symmetrize the cholesky factor matrix. */
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = i + 1; j < n; j++)
      matrix_set(L, i, j, matrix_get(L, j, i));

  /* return success. */
  return 1;
}

