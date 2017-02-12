
/* ensure once-only inclusion. */
#ifndef __VBNMR_GP_H__
#define __VBNMR_GP_H__

/* include c library headers. */
#include <stdio.h>
#include <stdlib.h>

/* include the required library headers. */
#include <vbnmr/cholesky.h>
#include <vbnmr/dataset.h>
#include <vbnmr/model.h>

/* gp_t: structure for holding a gaussian process nmr signal model.
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

  /* data fitting variables:
   *  @C: kernel matrix (cholesky factors) for the current dataset.
   *  @alpha: cholesky solution to C*alpha = y.
   */
  matrix_t *C;
  vector_t *alpha;

  /* function prediction variables:
   *  @ts: time vector of the current prediction.
   *  @ks: signal phase of the current prediction.
   *  @cs: vector of covariances for prediction.
   *  @beta: cholesky solution to C*beta = cs.
   *  @css: predictive variance upper bound.
   *  @mean: predictive mean.
   *  @var: predictive variance.
   */
  vector_t *t, *ts, *cs, *beta;
  double css, mean, var;
  unsigned int ks;

  /* temporary variables:
   *  @tmp: vector used for highly transient intermediate quantities.
   */
  vector_t *tmp;
}
gp_t;

/* function declarations (gp.c): */

gp_t *gp_alloc (dataset_t *dat, model_t *mdl);

void gp_free (gp_t *gp);

int gp_fit (gp_t *gp);

int gp_predict (gp_t *gp, const unsigned int ks, const vector_t *ts);

int gp_eval (gp_t *gp, dataset_t *mean, dataset_t *var);

#endif /* !__VBNMR_GP_H__ */

