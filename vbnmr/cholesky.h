
/* ensure once-only inclusion. */
#ifndef __VBNMR_CHOLESKY_H__
#define __VBNMR_CHOLESKY_H__

/* include the matrix, vector, and blas headers. */
#include <vbnmr/matrix.h>
#include <vbnmr/vector.h>
#include <vbnmr/blas.h>

/* function declarations (cholesky.c): */

int cholesky_decomp (matrix_t *A);

int cholesky_invert (const matrix_t *L, matrix_t *B);

void cholesky_solve (const matrix_t *L, const vector_t *b, vector_t *x);

void cholesky_update (matrix_t *L, vector_t *x);

int cholesky_downdate (matrix_t *L, vector_t *y);

#endif /* !__VBNMR_CHOLESKY_H__ */

