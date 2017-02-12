
/* ensure once-only inclusion. */
#ifndef __VBNMR_DATASET_H__
#define __VBNMR_DATASET_H__

/* include c headers. */
#include <stdio.h>
#include <stdlib.h>

/* include the matrix and vector headers. */
#include <vbnmr/matrix.h>
#include <vbnmr/vector.h>

/* dataset_t: structure for holding measured data.
 */
typedef struct {
  /* model space parameters:
   *  @D: number of dimensions.
   *  @K: number of phases.
   */
  unsigned int D, K;

  /* dataset size parameters:
   *  @n: number of measurements per phase.
   *  @N: total number of measurements.
   */
  unsigned int *n, N;

  /* core dataset arrays:
   *  @T: array of measurement time matrices.
   *  @y: array of measurement value vectors.
   */
  matrix_t **T;
  vector_t **y;
}
dataset_t;

/* function declarations (dataset.c): */

dataset_t *dataset_alloc (const unsigned int D);

dataset_t *dataset_alloc_from_grid (matrix_t *grid);

void dataset_free (dataset_t *dat);

void dataset_get (const dataset_t *dat,
                  const unsigned int k,
                  const unsigned int i,
                  vector_t *t,
                  double *y);

int dataset_augment (dataset_t *dat, const unsigned int k,
                     const vector_t *t, const double y);

int dataset_fwrite (dataset_t *dat, const char *fname);

int dataset_fread (dataset_t *dat, const char *fname);

#endif /* !__VBNMR_DATASET_H__ */

