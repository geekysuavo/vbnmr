
/* ensure once-only inclusion. */
#ifndef __VBNMR_MODEL_H__
#define __VBNMR_MODEL_H__

/* include c headers. */
#include <stdlib.h>
#include <math.h>

/* include the required headers. */
#include <vbnmr/dataset.h>
#include <vbnmr/rng.h>

/* VTAB(): macro for accessing model modulation table values.
 */
#define VTAB(k,l,d) mdl->vtab[mdl->D * (mdl->K * k + l) + d]

/* model_signal_t: structure for holding parameters and intermediate
 * quantities for a single frequency and decay rate pair.
 */
typedef struct {
  /* frequency parameters:
   *  @mu, @mu0: frequency mean.
   *  @tau, @tau0: frequency precision.
   */
  double mu, mu0;
  double tau, tau0;

  /* decay parameters:
   *  @alpha, @alpha0: decay shape.
   *  @beta, @beta0: decay rate.
   */
  double alpha, alpha0;
  double beta, beta0;

  /* fisher information matrices:
   *  @Fomega: diagonal matrix elements for @mu, @tau.
   *  @Frho: upper triangle matrix elements for @alpha, @beta.
   */
  double Fomega[2];
  double Frho[3];

  /* cached computations:
   *  @dga: digamma(x) function of @alpha.
   *  @tga: trigamma(x) function of @alpha.
   *  @lga: ln(Gamma(x)) function of @alpha.
   */
  double dga, tga, lga;
}
model_signal_t;

/* model_t: structure for holding bayesian nmr signal model parameters.
 */
typedef struct {
  /* model sizes:
   *  @D: number of dimensions.
   *  @K: number of phases.
   *  @M: number of signals.
   */
  unsigned int D, K, M;

  /* dimensionality-dependent tables:
   *  @vtab: modulation table.
   */
  int *vtab;

  /* standard deviations:
   *  @sigma: noise.
   *  @delta: signal.
   */
  double sigma, delta;

  /* nonlinear parameters:
   *  @sig: frequency and decay rate.
   */
  model_signal_t **sig;

  /* amplitude parameters:
   *  @ahat: mean vector.
   *  @Sigma: covariance matrix.
   */
  vector_t *ahat;
  matrix_t *Sigma;
}
model_t;

/* allocation function declarations: */

model_t *model_alloc (const unsigned int D, const unsigned int M);

void model_free (model_t *mdl);

/* prior function declarations: */

int model_prior_signal (model_t *mdl,
                        const unsigned int j,
                        const unsigned int d,
                        const double mu, const double tau,
                        const double alpha, const double beta);

int model_set_from_priors (model_t *mdl);

int model_randomize_mu (model_t *mdl, rng_t *gen);

/* parameter function declarations: */

int model_set_sigma (model_t *mdl, const double sigma);

int model_set_delta (model_t *mdl, const double delta);

int model_set_mu (model_t *mdl,
                  const unsigned int j,
                  const unsigned int d,
                  const double mu);

int model_set_tau (model_t *mdl,
                   const unsigned int j,
                   const unsigned int d,
                   const double tau);

int model_set_alpha (model_t *mdl,
                     const unsigned int j,
                     const unsigned int d,
                     const double alpha);

int model_set_beta (model_t *mdl,
                    const unsigned int j,
                    const unsigned int d,
                    const double beta);

/* evaluation and covariance function declarations: */

int model_eval (const model_t *mdl, rng_t *gen, dataset_t *dat);

double model_eval_single (const model_t *mdl,
                          const unsigned int k,
                          const vector_t *t);

int model_predict (const model_t *mdl, dataset_t *mean, dataset_t *var);

int model_predict_single (const model_t *mdl,
                          const unsigned int k,
                          const vector_t *t,
                          double *mean, double *var);

double model_cov (const model_t *mdl,
                  const unsigned int k1, const unsigned int k2,
                  const vector_t *t1, const vector_t *t2);

#endif /* !__VBNMR_MODEL_H__ */

