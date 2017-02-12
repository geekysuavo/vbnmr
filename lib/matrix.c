
/* include the matrix header. */
#include <vbnmr/matrix.h>

/* matrix_alloc(): allocate a new matrix for use.
 *
 * arguments:
 *  @rows: number of rows to allocate for the matrix.
 *  @cols: number of columns to allocate for the matrix.
 *
 * returns:
 *  newly allocated matrix structure pointer. the elements of the
 *  matrix will not yet be initialized.
 */
matrix_t *matrix_alloc (const unsigned int rows, const unsigned int cols) {
  /* declare required variables:
   *  @bytes: size of the memory block to allocate.
   *  @A: output matrix structure pointer.
   */
  unsigned int bytes;
  matrix_t *A;

  /* allocate a new structure pointer, or fail. */
  bytes = sizeof(matrix_t) + rows * cols * sizeof(double);
  A = (matrix_t*) malloc(bytes);
  if (!A)
    return NULL;

  /* store the structure parameters. */
  A->rows = rows;
  A->cols = cols;
  A->stride = cols;

  /* pointer the matrix data array to the end of the structure. */
  A->data = (double*) ((char*) A + sizeof(matrix_t));

  /* return the new structure pointer. */
  return A;
}

/* matrix_copy(): copy the contents of one matrix into another of
 * identical size.
 *
 * arguments:
 *  @dest: destination matrix structure pointer.
 *  @src: source matrix structure pointer.
 */
void matrix_copy (matrix_t *dest, const matrix_t *src) {
  /* copy each element without bounds checking. */
  for (unsigned int i = 0; i < dest->rows; i++)
    for (unsigned int j = 0; j < dest->cols; j++)
      matrix_set(dest, i, j, matrix_get(src, i, j));
}

/* matrix_copy_row(): copy the row of one matrix into a vector of
 * corforming size.
 *
 * arguments:
 *  @dest: destination vector structure pointer.
 *  @src: source matrix structure pointer.
 */
void matrix_copy_row (vector_t *dest, const matrix_t *src,
                      const unsigned int i) {
  /* copy the row elements without bounds checking. */
  for (unsigned int j = 0; j < dest->len; j++)
    vector_set(dest, j, matrix_get(src, i, j));
}

/* matrix_free(): free an allocated matrix.
 *
 * arguments:
 *  @A: matrix structure pointer to free.
 */
void matrix_free (matrix_t *A) {
  /* return if the structure pointer is null. */
  if (!A) return;

  /* free the structure pointer (this includes the elements). */
  free(A);
}

/* matrix_view_array(): create a matrix view from a raw array of doubles.
 *
 * arguments:
 *  @data: data array to point into.
 *  @n1: row count of the matrix view.
 *  @n2: column count of the matrix view.
 *
 * returns:
 *  newly created matrix view.
 */
matrix_view_t matrix_view_array (double *data,
                                 const unsigned int n1,
                                 const unsigned int n2) {
  /* initialize the view. */
  matrix_view_t view = { 0, 0, 0, NULL };

  /* the matrix view will access each array element in order. */
  view.rows = n1;
  view.cols = n2;
  view.stride = n2;
  view.data = data;

  /* return the view. */
  return view;
}

/* matrix_diag(): create a vector view of the diagonal entries of a matrix.
 *
 * arguments:
 *  @A: original matrix to point into.
 *
 * returns:
 *  newly created vector view.
 */
vector_view_t matrix_diag (const matrix_t *A) {
  /* initialize the view. */
  vector_view_t view = { 0, 0, NULL };

  /* stride the vector view to access each diagonal matrix element. */
  view.len = A->rows;
  view.stride = A->stride + 1;
  view.data = A->data;

  /* return the view. */
  return view;
}

/* matrix_row(): create a vector view of a single matrix row.
 *
 * arguments:
 *  @A: original matrix to point into.
 *  @i: matrix row index.
 *
 * returns:
 *  newly created vector view.
 */
vector_view_t matrix_row (const matrix_t *A, const unsigned int i) {
  /* initialize the view. */
  vector_view_t view = { 0, 0, NULL };

  /* the vector view will access a single row of matrix elements. */
  view.len = A->cols;
  view.stride = 1;
  view.data = A->data + (i * A->stride);

  /* return the view. */
  return view;
}

/* matrix_col(): create a vector view of a single matrix column.
 *
 * arguments:
 *  @A: original matrix to point into.
 *  @j: matrix column index.
 *
 * returns:
 *  newly created vector view.
 */
vector_view_t matrix_col (const matrix_t *A, const unsigned int j) {
  /* initialize the view. */
  vector_view_t view = { 0, 0, NULL };

  /* the vector view will access a single column of matrix elements. */
  view.len = A->rows;
  view.stride = A->stride;
  view.data = A->data + j;

  /* return the view. */
  return view;
}

/* matrix_subrow(): create a vector view of a partial matrix row.
 *
 * arguments:
 *  @A: original matrix to point into.
 *  @i: matrix row index.
 *  @offset: column offset.
 *  @n: number of viewed elements.
 *
 * returns:
 *  newly created vector view.
 */
vector_view_t matrix_subrow (const matrix_t *A, const unsigned int i,
                             const unsigned int offset,
                             const unsigned int n) {
  /* initialize the view. */
  vector_view_t view = { 0, 0, NULL };

  /* the vector view will access a partial row of matrix elements. */
  view.len = n;
  view.stride = 1;
  view.data = A->data + (i * A->stride + offset);

  /* return the view. */
  return view;
}

/* matrix_subcol(): create a vector view of a partial matrix column.
 *
 * arguments:
 *  @A: original matrix to point into.
 *  @j: matrix column index.
 *  @offset: row offset.
 *  @n: number of viewed elements.
 *
 * returns:
 *  newly created vector view.
 */
vector_view_t matrix_subcol (const matrix_t *A, const unsigned int j,
                             const unsigned int offset,
                             const unsigned int n) {
  /* initialize the view. */
  vector_view_t view = { 0, 0, NULL };

  /* the vector view will access a partial column of matrix elements. */
  view.len = n;
  view.stride = A->stride;
  view.data = A->data + (offset * A->stride + j);

  /* return the view. */
  return view;
}

/* matrix_submatrix(): create a matrix view of a portion of an
 * existing matrix.
 *
 * arguments:
 *  @A: original matrix to point into.
 *  @i1: submatrix row offset.
 *  @i2: submatrix column offset.
 *  @n1: submatrix row count.
 *  @n2: submatrix column count.
 *
 * returns:
 *  newly created matrix view.
 */
matrix_view_t matrix_submatrix (const matrix_t *A,
                                const unsigned int i1,
                                const unsigned int i2,
                                const unsigned int n1,
                                const unsigned int n2) {
  /* initialize the view. */
  matrix_view_t view = { 0, 0, 0, NULL };

  /* use the striding of the input matrix to construct the view. */
  view.rows = n1;
  view.cols = n2;
  view.stride = A->stride;
  view.data = A->data + (i1 * A->stride + i2);

  /* return the view. */
  return view;
}

/* matrix_get(): get the value of a matrix element.
 *
 * arguments:
 *  @A: matrix to access.
 *  @i: row index.
 *  @j: column index.
 *
 * returns:
 *  value of the requested matrix element.
 */
inline double matrix_get (const matrix_t *A,
                          const unsigned int i,
                          const unsigned int j) {
  /* return the element without bounds checking. */
  return A->data[i * A->stride + j];
}

/* matrix_set(): set the value of a matrix element.
 *
 * arguments:
 *  @A: matrix to access.
 *  @i: row index.
 *  @j: column index.
 *  @Aij: new element value.
 */
inline void matrix_set (matrix_t *A,
                        const unsigned int i,
                        const unsigned int j,
                        const double Aij) {
  /* set the element without bounds checking. */
  A->data[i * A->stride + j] = Aij;
}

/* matrix_set_all(): set all elements of a matrix to a specified value.
 *
 * arguments:
 *  @A: matrix to modify.
 *  @Aall: new element value.
 */
void matrix_set_all (matrix_t *A, const double Aall) {
  /* set all elements of the matrix. */
  for (unsigned int i = 0; i < A->rows; i++)
    for (unsigned int j = 0; j < A->cols; j++)
      matrix_set(A, i, j, Aall);
}

/* matrix_set_ident(): set the elements of a matrix to the identity.
 *
 * arguments:
 *  @A: matrix to modify.
 */
void matrix_set_ident (matrix_t *A) {
  /* intelligently set all elements of the matrix. */
  for (unsigned int i = 0; i < A->rows; i++)
    for (unsigned int j = 0; j < A->cols; j++)
      matrix_set(A, i, j, i == j ? 1.0 : 0.0);
}

/* matrix_set_zero(): set all elements of a matrix to zero.
 *
 * arguments:
 *  @A: matrix to modify.
 */
inline void matrix_set_zero (matrix_t *A) {
  /* zero all elements of the matrix. */
  matrix_set_all(A, 0.0);
}

/* matrix_sub(): subtract one matrix from another element-wise.
 *
 * operation:
 *  A <- A - B
 *
 * arguments:
 *  @A: first input and output matrix.
 *  @B: second input matrix.
 */
void matrix_sub (matrix_t *A, const matrix_t *B) {
  /* perform the element-wise difference. */
  for (unsigned int i = 0; i < A->rows; i++)
    for (unsigned int j = 0; j < A->cols; j++)
      matrix_set(A, i, j, matrix_get(A, i, j) - matrix_get(B, i, j));
}

/* matrix_scale(): scale all elements of a matrix by a constant value.
 *
 * arguments:
 *  @A: input and output matrix.
 *  @alpha: scale factor.
 */
void matrix_scale (matrix_t *A, const double alpha) {
  /* perform the element-wise scaling. */
  for (unsigned int i = 0; i < A->rows; i++)
    for (unsigned int j = 0; j < A->cols; j++)
      matrix_set(A, i, j, matrix_get(A, i, j) * alpha);
}

/* matrix_dispfn(): display the name and contents of a matrix.
 *
 * arguments:
 *  @A: input matrix structure pointer.
 *  @str: variable name of the matrix.
 */
void matrix_dispfn (const matrix_t *A, const char *str) {
  /* print the variable name. */
  printf("%s =\n", str);

  /* print matrix rows. */
  for (unsigned int i = 0; i < A->rows; i++) {
    /* print the values of the row. */
    for (unsigned int j = 0; j < A->cols; j++)
      printf("  %13.5lf", matrix_get(A, i, j));

    /* print a newline. */
    printf("\n");
  }
}

