
/* include the dataset header. */
#include <vbnmr/dataset.h>

/* dataset_alloc(): allocate a new empty dataset for use.
 *
 * arguments:
 *  @D: number of signal dimensions.
 *
 * returns:
 *  newly allocated and initialized dataset structure pointer. if @D is
 *  zero, or any allocations fail, then a null pointer is returned.
 */
dataset_t *dataset_alloc (const unsigned int D) {
  /* @K: number of signal phases.
   * @bytes: size of the allocation to perform.
   * @dat: pointer to the new dataset structure.
   * @ptr: byte-wide pointer for initializing arrays.
   */
  unsigned int K, bytes;
  dataset_t *dat;
  char *ptr;

  /* fail if a zero-dimensional dataset was requested. */
  if (D == 0)
    return NULL;

  /* compute the number of signal phases. */
  K = 1;
  for (unsigned int d = 0; d < D; d++)
    K *= 2;

  /* compute the number of bytes to allocate for the structure. */
  bytes = sizeof(dataset_t);
  bytes += K * sizeof(unsigned int);
  bytes += K * sizeof(matrix_t*);
  bytes += K * sizeof(vector_t*);

  /* allocate the structure, or fail. */
  dat = (dataset_t*) malloc(bytes);
  if (!dat)
    return NULL;

  /* store the dimensionality, phase count and measurement count. */
  dat->D = D;
  dat->K = K;
  dat->N = 0;

  /* initialize the array of measurement counts. */
  ptr = (char*) dat;
  ptr += sizeof(dataset_t);
  dat->n = (unsigned int*) ptr;

  /* initialize the array of measurement time matrices. */
  ptr += K * sizeof(unsigned int);
  dat->T = (matrix_t**) ptr;

  /* initialize the array of measured value vectors. */
  ptr += K * sizeof(matrix_t*);
  dat->y = (vector_t**) ptr;

  /* initialize the contents of all arrays. */
  for (unsigned int k = 0; k < K; k++) {
    /* zero the measurement count. */
    dat->n[k] = 0;

    /* initialize the measurement matrix, vector. */
    dat->T[k] = NULL;
    dat->y[k] = NULL;
  }

  /* return the new dataset. */
  return dat;
}

/* augment_from_grid(): recursive function used by dataset_alloc_from_grid()
 * to augment a dataset with a uniform grid of values on a specified phase.
 *
 * arguments:
 *  @dat: dataset structure pointer to augment.
 *  @grid: matrix of gridding information.
 *  @t: vector of time values used during recursion.
 *  @k: phase index to augment times into.
 *  @d: current recursion dimension.
 */
static void 
augment_from_grid (dataset_t *dat, matrix_t *grid, vector_t *t,
                   const unsigned int k, const unsigned int d) {
  /* initialize the current dimension time. */
  vector_set(t, d, matrix_get(grid, d, 0));

  /* loop over the current dimension time values. */
  while (vector_get(t, d) <= matrix_get(grid, d, 2)) {
    /* augment or recurse. */
    if (d == dat->D - 1) {
      /* current dimension is last: augment the dataset. */
      dataset_augment(dat, k, t, 0.0);
    }
    else {
      /* current dimension is inner: recurse another level. */
      augment_from_grid(dat, grid, t, k, d + 1);
    }

    /* increment the current dimension time. */
    vector_set(t, d, vector_get(t, d) + matrix_get(grid, d, 1));
  }
}

/* dataset_alloc_from_grid(): allocate a new dataset with a grid
 * of time values and zero measurements.
 *
 * arguments:
 *  @grid: matrix of gridding information.
 *    - one row per dataset dimension.
 *    - three columns: start, step, end.
 *
 * returns:
 *  newly allocated and initialized dataset structure pointer.
 */
dataset_t *dataset_alloc_from_grid (matrix_t *grid) {
  /* check the grid matrix. */
  if (!grid || grid->cols != 3)
    return NULL;

  /* get the dimensionality of the grid matrix. */
  const unsigned int D = grid->rows;

  /* allocate the dataset structure pointer. */
  dataset_t *dat = dataset_alloc(D);
  if (!dat)
    return NULL;

  /* allocate a vector for storing each time. */
  vector_t *t = vector_alloc(D);
  if (!t) {
    dataset_free(dat);
    return NULL;
  }

  /* recursively augment the dataset. */
  for (unsigned int k = 0; k < dat->K; k++)
    augment_from_grid(dat, grid, t, k, 0);

  /* free the time vector and return the new dataset. */
  vector_free(t);
  return dat;
}

/* dataset_free(): free an allocated dataset.
 *
 * arguments:
 *  @dat: dataset structure pointer to free.
 */
void dataset_free (dataset_t *dat) {
  /* return if the structure pointer is null. */
  if (!dat) return;

  /* loop over the signal phases. */
  for (unsigned int k = 0; k < dat->K; k++) {
    /* skip empty phases. */
    if (!dat->n[k]) continue;

    /* free the measurement matrix, vector. */
    matrix_free(dat->T[k]);
    vector_free(dat->y[k]);
  }

  /* free the dataset structure pointer. */
  free(dat);
}

/* dataset_get(): extract a time and measurement from a dataset.
 *
 * note:
 *  no argument checking is performed by this function, so be sure
 *  to use it safely!
 *
 * arguments:
 *  @dat: pointer to the dataset structure to access.
 *  @k, @i: phase and measurement index to extract.
 *  @t: output time vector.
 *  @y: output measurement.
 */
inline void dataset_get (const dataset_t *dat,
                         const unsigned int k,
                         const unsigned int i,
                         vector_t *t, double *y) {
  /* extract the requested time and measurement. */
  matrix_copy_row(t, dat->T[k], i);
  *y = vector_get(dat->y[k], i);
}

/* dataset_augment(): add a new datum into a dataset.
 *
 * arguments:
 *  @dat: pointer to the dataset structure to modify.
 *  @k: index of the signal phase of the new datum.
 *  @t: measurement time of the new datum.
 *  @y: value of the new datum.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int dataset_augment (dataset_t *dat, const unsigned int k,
                     const vector_t *t, const double y) {
  /* @n: new number of measurements along phase @k.
   * @N: new total number of measurements.
   * @D: dataset dimensionality.
   * @Tnew: new measurement time matrix.
   * @ynew: new measured value vector.
   */
  unsigned int n, N, D;
  matrix_t *Tnew;
  vector_t *ynew;

  /* fail if:
   *  - any passed pointers are null.
   *  - the phase index is out of the dataset bounds.
   *  - the time vector length does not match the dataset.
   */
  if (!dat || !t || k >= dat->K || t->len != dat->D)
    return 0;

  /* store local copies of the new counts. */
  n = dat->n[k] + 1;
  N = dat->N + 1;
  D = dat->D;

  /* allocate a new measurement time matrix, or fail. */
  Tnew = matrix_alloc(n, D);
  if (!Tnew)
    return 0;

  /* allocate a new measured value vector, or fail. */
  ynew = vector_alloc(n);
  if (!ynew) {
    matrix_free(Tnew);
    return 0;
  }

  /* loop over the previous contents of the measurement. */
  for (unsigned int i = 0; i < n - 1; i++) {
    /* store the measurement times. */
    for (unsigned int d = 0; d < D; d++)
      matrix_set(Tnew, i, d, matrix_get(dat->T[k], i, d));

    /* store the measured values. */
    vector_set(ynew, i, vector_get(dat->y[k], i));
  }

  /* store the new measurement time. */
  for (unsigned int d = 0; d < D; d++)
    matrix_set(Tnew, n - 1, d, vector_get(t, d));

  /* store the new measured value. */
  vector_set(ynew, n - 1, y);

  /* free the old measurement matrix, vector. */
  free(dat->T[k]);
  free(dat->y[k]);

  /* store the updated information into the dataset. */
  dat->T[k] = Tnew;
  dat->y[k] = ynew;
  dat->n[k] = n;
  dat->N = N;

  /* return success. */
  return 1;
}

/* dataset_fwrite(): write the contents of a dataset to a text-format file.
 *
 * arguments:
 *  @dat: pointer to a dataset structure to access.
 *  @fname: output filename to write to.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int dataset_fwrite (dataset_t *dat, const char *fname) {
  /* declare required variables:
   *  @fh: output file handle.
   */
  FILE *fh;

  /* check the pointers. */
  if (!dat || !fname)
    return 0;

  /* open the output file. */
  fh = fopen(fname, "w");
  if (!fh)
    return 0;

  /* output each entry to a line of the file. */
  for (unsigned int k = 0; k < dat->K; k++) {
    /* loop over each time of the current phase. */
    for (unsigned int i = 0; i < dat->n[k]; i++) {
      /* begin the line. */
      fprintf(fh, "%u", k);

      /* get the measurement time and value. */
      vector_view_t t = matrix_row(dat->T[k], i);
      const double y = vector_get(dat->y[k], i);

      /* print the time vector elements. */
      for (unsigned int d = 0; d < dat->D; d++)
        fprintf(fh, " %le", vector_get(&t, d));

      /* end the line with the measured value. */
      fprintf(fh, " %le\n", y);
    }
  }

  /* close the output file and return success. */
  fclose(fh);
  return 1;
}

/* dataset_fread(): read the contents of a text-format file into a dataset.
 *
 * arguments:
 *  @dat: pointer to a dataset structure to access.
 *  @fname: input filename to read from.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int dataset_fread (dataset_t *dat, const char *fname) {
  /* declare required variables:
   *  @k: phase index of each read line.
   *  @ntok: token count of the first read line.
   *  @status: return status of the function.
   *  @buf: string buffer for reading each line.
   *  @x: vector for holding measurement times.
   *  @y: variable for holding measurements.
   *  @fh: input file handle.
   */
  unsigned int k, ntok, status;
  char buf[1024], *tok;
  vector_t *x;
  double y;
  FILE *fh;

  /* check the pointers. */
  if (!dat || !fname)
    return 0;

  /* open the output file. */
  fh = fopen(fname, "r");
  if (!fh)
    return 0;

  /* initialize the status variable. */
  status = 0;
  x = NULL;

  /* count the number of tokens in the first line. */
  if (fgets(buf, 1024, fh)) {
    /* read the first token. */
    tok = strtok(buf, " ");
    ntok = 0;

    /* count the rest of the tokens. */
    while (tok) {
      tok = strtok(NULL, " ");
      ntok++;
    }
  }
  else
    goto fail;

  /* check that the first line has at least three columns. */
  if (ntok < 3)
    goto fail;

  /* allocate a vector for storing times. */
  x = vector_alloc(ntok - 2);
  if (!x)
    goto fail;

  /* rewind and read the complete file. */
  fseek(fh, SEEK_SET, 0);
  while (!feof(fh)) {
    /* read a new line of the file. */
    if (fgets(buf, 1024, fh)) {
      /* read the phase index. */
      tok = strtok(buf, " ");
      k = atol(tok);

      /* read the time vector. */
      tok = strtok(NULL, " ");
      for (unsigned int d = 0; d < x->len; d++) {
        vector_set(x, d, atof(tok));
        tok = strtok(NULL, " ");
      }

      /* read the measurement. */
      y = atof(tok);

      /* augment the dataset. */
      if (!dataset_augment(dat, k, x, y))
        goto fail;
    }
  }

  /* indicate success. */
  status = 1;

fail:
  /* free the temporary vector and close the input file . */
  vector_free(x);
  fclose(fh);

  /* return the status variable. */
  return status;
}

