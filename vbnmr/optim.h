
/* ensure once-only inclusion. */
#ifndef __VBNMR_OPTIM_H__
#define __VBNMR_OPTIM_H__

/* include the required library headers. */
#include <vbnmr/cholesky.h>
#include <vbnmr/dataset.h>
#include <vbnmr/model.h>

/* optim_t: structure for holding a signal model optimization state.
 */
typedef struct {
  /* model and data space sizes:
   *  @D: number of dimensions.
   *  @K: number of phases.
   *  @M: number of signals.
   *  @N: number of measurements.
   */
  unsigned int D, K, M, N;

  /* model and dataset structure pointers:
   *  @dat: pointer to an associated dataset.
   *  @mdl: pointer to an associated signal model.
   */
  dataset_t *dat;
  model_t *mdl;

  /* hyperparameter learning variables:
   *  @L: precision matrix cholesky decomposition.
   *  @G: expected interaction matrix.
   *  @h: expected projection vector.
   *  @elbo: value of the lower bound.
   */
  matrix_t *L, *G;
  vector_t *h;
  double elbo;

  /* iteration control variables:
   *  @iter_max: maximum number of free-run iterations.
   *  @step_max: maximum number of line search steps.
   *  @l0: initial lipschitz constant estimate.
   *  @dl: lipschitz constant step factor.
   *  @log: file handle for logging.
   */
  unsigned int iter_max, step_max;
  double l0, dl;
  FILE *log;

  /* temporary variables:
   *  @tmp: vector used for highly transient intermediate quantities.
   */
  vector_t *tmp;
}
optim_t;

/* function declarations (optim.c): */

optim_t *optim_alloc (model_t *mdl, dataset_t *dat);

void optim_free (optim_t *opt);

int optim_init (optim_t *opt);

int optim_execute (optim_t *opt);

#endif /* !__VBNMR_OPTIM_H__ */

