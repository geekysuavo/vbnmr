
/* ensure once-only inclusion. */
#ifndef __VBNMR_BLAS_H__
#define __VBNMR_BLAS_H__

/* include the matrix and vector headers. */
#include <vbnmr/matrix.h>
#include <vbnmr/vector.h>

/* blas_transpose_t: enumeration of all possible ways to transpose
 * a matrix (or not) during calculations.
 */
typedef enum {
  BLAS_NO_TRANS = 0,
  BLAS_TRANS
}
blas_transpose_t;

/* blas_triangle_t: enumeration of all possible ways to access
 * a triangular matrix.
 */
typedef enum {
  BLAS_LOWER = 0,
  BLAS_UPPER
}
blas_triangle_t;

/* function declarations (blas.c): */

double blas_dnrm2 (const vector_t *x);

double blas_ddot (const vector_t *x, const vector_t *y);

void blas_daxpy (double alpha, const vector_t *x, vector_t *y);

/* --- */

void blas_dgemv (blas_transpose_t trans, double alpha, const matrix_t *A,
                 const vector_t *x, double beta, vector_t *y);

void blas_dtrmv (blas_transpose_t trans, double alpha, const matrix_t *L,
                 const vector_t *x, double beta, vector_t *y);

void blas_dtrsv (blas_triangle_t tri, const matrix_t *L, vector_t *x);

#endif /* !__VBNMR_MATRIX_H__ */

