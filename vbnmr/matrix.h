
/* ensure once-only inclusion. */
#ifndef __VBNMR_MATRIX_H__
#define __VBNMR_MATRIX_H__

/* include c library headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* include the vector header. */
#include <vbnmr/vector.h>

/* matrix_disp(): macro function for displaying the contents of matrices.
 */
#define matrix_disp(A) matrix_dispfn(A, #A)

/* matrix_t: structure for holding a real matrix.
 */
typedef struct {
  /* data access parameters:
   *  @rows: number of matrix rows.
   *  @cols: number of matrix columns.
   *  @stride: array spacing between elements.
   */
  unsigned int rows, cols;
  unsigned int stride;

  /* @data: array or pointer to matrix elements.
   */
  double *data;
}
matrix_t;

/* matrix_view_t: aliased type used to explicitly note that a given
 * matrix data structure does not own its data array.
 */
typedef matrix_t matrix_view_t;

/* function declarations (matrix.c): */

matrix_t *matrix_alloc (const unsigned int rows, const unsigned int cols);

void matrix_copy (matrix_t *dest, const matrix_t *src);

void matrix_copy_row (vector_t *dest, const matrix_t *src,
                      const unsigned int i);

void matrix_free (matrix_t *A);

matrix_view_t matrix_view_array (double *data,
                                 const unsigned int n1,
                                 const unsigned int n2);

vector_view_t matrix_diag (const matrix_t *A);

vector_view_t matrix_row (const matrix_t *A, const unsigned int i);

vector_view_t matrix_col (const matrix_t *A, const unsigned int j);

vector_view_t matrix_subrow (const matrix_t *A, const unsigned int i,
                             const unsigned int offset,
                             const unsigned int n);

vector_view_t matrix_subcol (const matrix_t *A, const unsigned int j,
                             const unsigned int offset,
                             const unsigned int n);

matrix_view_t matrix_submatrix (const matrix_t *A,
                                const unsigned int i1,
                                const unsigned int i2,
                                const unsigned int n1,
                                const unsigned int n2);

double matrix_get (const matrix_t *A,
                   const unsigned int i,
                   const unsigned int j);

void matrix_set (matrix_t *A, const unsigned int i, const unsigned int j,
                 const double Aij);

void matrix_set_all (matrix_t *A, const double Aall);

void matrix_set_ident (matrix_t *A);

void matrix_set_zero (matrix_t *A);

void matrix_sub (matrix_t *A, const matrix_t *B);

void matrix_scale (matrix_t *A, const double alpha);

void matrix_dispfn (const matrix_t *A, const char *str);

#endif /* !__VBNMR_MATRIX_H__ */

