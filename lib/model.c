
/* include the signal model header. */
#include <vbnmr/model.h>

/* include the distribution headers. */
#include <vbnmr/dist-normal.h>
#include <vbnmr/dist-gamma.h>
#include <vbnmr/dist-basis.h>

/* include the polygamma header. */
#include <vbnmr/psi.h>

/* * * * modulation tables: * * * */

/* phase offset legend:
 *  0:   cos
 *  1:  -sin
 *  2:  -cos
 *  3:   sin
 */

/* V1: modulation table for D = 1. */
const int V1[] = {
  /* k = 0 */
  0,
  1,

  /* k = 1 */
  3,
  0
};

/* V2: modulation table for D = 2. */
const int V2[] = {
  /* k = 0 */
  0, 0,
  1, 0,
  0, 1,
  1, 1,

  /* k = 1 */
  3, 0,
  0, 0,
  3, 1,
  0, 1,

  /* k = 2 */
  0, 3,
  1, 3,
  0, 0,
  1, 0,

  /* k = 3 */
  3, 3,
  0, 3,
  3, 0,
  0, 0,
};

/* * * * allocation and deallocation: * * * */

/* model_init(): initialize the pointers within a model structure to null.
 *
 * arguments:
 *  @mdl: pointer to the model structure to initialize.
 */
static void model_init (model_t *mdl) {
  /* initialize the amplitude parameters. */
  mdl->ahat = NULL;
  mdl->Sigma = NULL;
}

/* model_alloc(): allocate a new empty model for use.
 *
 * arguments:
 *  @D: number of signal dimensions.
 *  @M: number of model components.
 *
 * returns:
 *  newly allocated and initialized model structure pointer. if @D or @M
 *  are zero, or any allocations fail, then a null pointer is returned.
 */
model_t *model_alloc (const unsigned int D, const unsigned int M) {
  /* @K: number of signal phases.
   * @bytes: size of the allocation to perform.
   * @mdl: pointer to the new model structure.
   * @ptr: byte-wide pointer for initializing arrays.
   */
  unsigned int K, bytes;
  model_t *mdl;
  char *ptr;

  /* fail if a zero-size model was requested. */
  if (D == 0 || M == 0)
    return NULL;

  /* compute the number of signal phases. */
  K = 1;
  for (unsigned int d = 0; d < D; d++)
    K *= 2;

  /* compute the number of bytes to allocate for the structure. */
  bytes = sizeof(model_t);
  bytes += M * sizeof(model_signal_t*);
  bytes += M * D * sizeof(model_signal_t);

  /* allocate the structure, or fail. */
  mdl = (model_t*) malloc(bytes);
  if (!mdl)
    return NULL;

  /* store the appropriate tables into the model. */
  switch (D) {
    case 1: /* D = 1 */
      mdl->vtab = (int*) V1;
      break;

    case 2: /* D = 2 */
      mdl->vtab = (int*) V2;
      break;

    /* otherwise. */
    default:
      goto fail;
  }

  /* store the dimensionality, phase count and signal count. */
  mdl->D = D;
  mdl->K = K;
  mdl->M = M;

  /* initialize the scalars and pointers. */
  mdl->sigma = 1.0;
  mdl->delta = 1.0;
  model_init(mdl);

  /* initialize the offset pointer. */
  ptr = (char*) mdl;
  ptr += sizeof(model_t);

  /* initialize the signal array. */
  mdl->sig = (model_signal_t**) ptr;
  ptr += M * sizeof(model_signal_t*);
  for (unsigned int m = 0; m < M; m++) {
    mdl->sig[m] = (model_signal_t*) ptr;
    ptr += D * sizeof(model_signal_t);
  }

  /* allocate the amplitude parameters. */
  mdl->ahat  = vector_alloc(K * M);
  mdl->Sigma = matrix_alloc(K * M, K * M);
  if (!mdl->ahat || !mdl->Sigma)
    goto fail;

  /* initialize the amplitude parameters. */
  vector_set_zero(mdl->ahat);
  matrix_set_ident(mdl->Sigma);
  matrix_scale(mdl->Sigma, mdl->delta * mdl->delta);

  /* return the new model. */
  return mdl;

fail:
  /* free any allocated memory and return a null pointer. */
  model_free(mdl);
  return NULL;
}

/* model_free(): free an allocated model.
 *
 * arguments:
 *  @mdl: model structure pointer to free.
 */
void model_free (model_t *mdl) {
  /* return if the structure pointer is null. */
  if (!mdl) return;

  /* free the amplitude parameters. */
  vector_free(mdl->ahat);
  matrix_free(mdl->Sigma);

  /* free the model structure pointer. */
  model_init(mdl);
  free(mdl);
}

/* * * * prior adjustment and initialization: * * * */

/* model_prior_signal(): set the prior parameters of a signal in
 * a model structure.
 *
 * arguments:
 *  @mdl: pointer to the model structure to modify.
 *  @j, @d: signal basis index and dimension.
 *  @mu, @tau: prior frequency mean, precision.
 *  @alpha, @beta: prior decay shape, rate.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int model_prior_signal (model_t *mdl,
                        const unsigned int j,
                        const unsigned int d,
                        const double mu, const double tau,
                        const double alpha, const double beta) {
  /* check validity of the arguments. */
  if (!mdl || j >= mdl->M || d >= mdl->D)
    return 0;

  /* set the frequency priors. */
  mdl->sig[j][d].mu0 = mu;
  mdl->sig[j][d].tau0 = tau;

  /* set the decay priors. */
  mdl->sig[j][d].alpha0 = alpha;
  mdl->sig[j][d].beta0 = beta;

  /* return success. */
  return 1;
}

/* model_set_from_priors(): set the values of all model parameters
 * using the current values of their associated priors.
 *
 * arguments:
 *  @mdl: pointer to the model structure to modify.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int model_set_from_priors (model_t *mdl) {
  /* check that the passed model pointer is not null. */
  if (!mdl)
    return 0;

  /* loop over the model components. */
  for (unsigned int j = 0; j < mdl->M; j++) {
    /* loop over the model dimensions. */
    for (unsigned int d = 0; d < mdl->D; d++) {
      /* set the frequency parameters. */
      if (!model_set_mu(mdl, j, d, mdl->sig[j][d].mu0) ||
          !model_set_tau(mdl, j, d, mdl->sig[j][d].tau0))
        return 0;

      /* set the decay parameters. */
      if (!model_set_alpha(mdl, j, d, mdl->sig[j][d].alpha0) ||
          !model_set_beta(mdl, j, d, mdl->sig[j][d].beta0))
        return 0;
    }
  }

  /* initialize the amplitude parameters. */
  vector_set_zero(mdl->ahat);
  matrix_set_ident(mdl->Sigma);
  matrix_scale(mdl->Sigma, mdl->delta * mdl->delta);

  /* return success. */
  return 1;
}

/* model_randomize_mu(): sample random frequency means for each signal
 * in a model, based on the model's current frequency priors.
 *
 * arguments:
 *  @mdl: pointer to the model structure to modify.
 *  @gen: random number generator structure pointer.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int model_randomize_mu (model_t *mdl, rng_t *gen) {
  /* declare required variables:
   *  @m0, @s0: mean and standard deviation prior values.
   *  @mu: randomly sampled mean value.
   */
  double m0, s0, mu;

  /* check that the passed pointers are not null. */
  if (!mdl || !gen)
    return 0;

  /* loop over the model components. */
  for (unsigned int j = 0; j < mdl->M; j++) {
    /* loop over the model dimensions. */
    for (unsigned int d = 0; d < mdl->D; d++) {
      m0 = mdl->sig[j][d].mu0;
      s0 = 1.0 / sqrt(mdl->sig[j][d].tau0);

      /* sample and set a new frequency mean. */
      mu = m0 + s0 * rng_normal(gen);
      model_set_mu(mdl, j, d, mu);
    }
  }

  /* return success. */
  return 1;
}

/* * * * parameter adjustment: * * * */

/* model_set_sigma(): set the noise standard deviation of a model.
 *
 * arguments:
 *  @mdl: pointer to the model structure to modify.
 *  @sigma: new noise standard deviation.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int model_set_sigma (model_t *mdl, const double sigma) {
  /* check validity of the arguments. */
  if (!mdl || sigma <= 0.0)
    return 0;

  /* set the new value and return success. */
  mdl->sigma = sigma;
  return 1;
}

/* model_set_delta(): set the amplitude standard deviation of a model.
 *
 * arguments:
 *  @mdl: pointer to the model structure to modify.
 *  @delta: new amplitude standard deviation.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int model_set_delta (model_t *mdl, const double delta) {
  /* check validity of the arguments. */
  if (!mdl || delta <= 0.0)
    return 0;

  /* set the new value and return success. */
  mdl->delta = delta;
  return 1;
}

/* model_set_mu(): set the frequency mean of a signal in a model.
 *
 * arguments:
 *  @mdl: pointer to the model structure to modify.
 *  @j, @d: signal basis index and dimension.
 *  @mu: new value to assign.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int model_set_mu (model_t *mdl,
                  const unsigned int j,
                  const unsigned int d,
                  const double mu) {
  /* check validity of the arguments. */
  if (!mdl || j >= mdl->M || d >= mdl->D)
    return 0;

  /* store the new value. */
  mdl->sig[j][d].mu = mu;

  /* return success. */
  return 1;
}

/* model_set_tau(): set the frequency precision
 * of a signal in a model.
 *
 * arguments:
 *  @mdl: pointer to the model structure to modify.
 *  @j, @d: signal basis index and dimension.
 *  @tau: new value to assign.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int model_set_tau (model_t *mdl,
                   const unsigned int j,
                   const unsigned int d,
                   const double tau) {
  /* check validity of the arguments. */
  if (!mdl || j >= mdl->M || d >= mdl->D || tau <= 0.0)
    return 0;

  /* store the new value. */
  mdl->sig[j][d].tau = tau;

  /* set the fisher matrix elements. */
  mdl->sig[j][d].Fomega[0] = tau;
  mdl->sig[j][d].Fomega[1] = 0.75 / (tau * tau);

  /* return success. */
  return 1;
}

/* model_set_alpha(): set the decay shape parameter of a signal in a model.
 *
 * arguments:
 *  @mdl: pointer to the model structure to modify.
 *  @j, @d: signal basis index and dimension.
 *  @alpha: new value to assign.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int model_set_alpha (model_t *mdl,
                     const unsigned int j,
                     const unsigned int d,
                     const double alpha) {
  /* check validity of the arguments. */
  if (!mdl || j >= mdl->M || d >= mdl->D || alpha <= 0.0)
    return 0;

  /* store the new values. */
  mdl->sig[j][d].alpha = alpha;
  mdl->sig[j][d].lga = lgamma(alpha);
  mdl->sig[j][d].dga = digamma(alpha);
  mdl->sig[j][d].tga = trigamma(alpha);

  /* set the fisher matrix elements. */
  mdl->sig[j][d].Frho[0] = mdl->sig[j][d].tga;
  mdl->sig[j][d].Frho[1] = alpha / pow(mdl->sig[j][d].beta, 2.0);

  /* return success. */
  return 1;
}

/* model_set_beta(): set the decay rate parameter of a signal in a model.
 *
 * arguments:
 *  @mdl: pointer to the model structure to modify.
 *  @j, @d: signal basis index and dimension.
 *  @beta: new value to assign.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int model_set_beta (model_t *mdl,
                    const unsigned int j,
                    const unsigned int d,
                    const double beta) {
  /* check validity of the arguments. */
  if (!mdl || j >= mdl->M || d >= mdl->D || beta <= 0.0)
    return 0;

  /* store the new value. */
  mdl->sig[j][d].beta = beta;

  /* set the fisher matrix elements. */
  mdl->sig[j][d].Frho[1] = mdl->sig[j][d].alpha / (beta * beta);
  mdl->sig[j][d].Frho[2] = -1.0 / beta;

  /* return success. */
  return 1;
}

/* * * * function and covariance evaluation: * * * */

/* model_eval(): evaluate the maximum a posteriori (MAP) estimate
 * of a model signal for every time/phase combination in a dataset.
 *
 * arguments:
 *  @mdl: pointer to the model structure to access.
 *  @gen: random number generator, or null to omit noise.
 *  @dat: dataset to fill with model evaluations.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the evaluation succeeded.
 */
int model_eval (const model_t *mdl, rng_t *gen, dataset_t *dat) {
  /* check the input arguments. */
  if (!mdl || !dat)
    return 0;

  /* loop over each dataset entry. */
  for (unsigned int k = 0; k < dat->K; k++) {
    for (unsigned int i = 0; i < dat->n[k]; i++) {
      /* evaluate the model at the current time. */
      vector_view_t t = matrix_row(dat->T[k], i);
      double y = model_eval_single(mdl, k, &t);

      /* add noise if a random number generator was supplied. */
      if (gen)
        y += rng_normal(gen) * mdl->sigma;

      /* store the evaluation. */
      vector_set(dat->y[k], i, y);
    }
  }

  /* return success. */
  return 1;
}

/* model_eval_single(): evaluate the maximum a posteriori (MAP) estimate
 * of a model signal on a given phase, at a given time.
 *
 * arguments:
 *  @mdl: pointer to the model structure to access.
 *  @k: signal phase index to sample from.
 *  @t: measurement time to sample at.
 *
 * returns:
 *  estimate of the signal at time @t on phase @k, using the mean
 *  signal parameters.
 */
double model_eval_single (const model_t *mdl,
                          const unsigned int k,
                          const vector_t *t) {
  /* @y: output signal value.
   */
  double y = 0.0;

  /* check the input arguments. */
  if (!mdl || !t || k >= mdl->K || t->len != mdl->D)
    return 0.0;

  /* gain a reference to the model sizes. */
  const unsigned int D = mdl->D;
  const unsigned int K = mdl->K;
  const unsigned int M = mdl->M;

  /* loop over the model signals. */
  for (unsigned int j = 0; j < M; j++) {
    /* compute the decay factor. */
    double R = 1.0;
    for (unsigned int d = 0; d < D; d++) {
      /* extract the mean parameters. */
      const double rho = mdl->sig[j][d].alpha / mdl->sig[j][d].beta;
      const double td = vector_get(t, d);

      /* include the current dimension term. */
      R *= exp(-rho * td);
    }

    /* loop over each term. */
    for (unsigned int l = 0; l < K; l++) {
      /* get the term amplitude. */
      const double ajl = vector_get(mdl->ahat, l * M + j);

      /* compute the modulation factor. */
      double V = 1.0;
      for (unsigned int d = 0; d < D; d++) {
        /* extract the mean parameters. */
        const double omega = mdl->sig[j][d].mu;
        const double td = vector_get(t, d);

        /* include the current dimension term. */
        V *= cos(omega * td + M_PI_2 * VTAB(k,l,d));
      }

      /* sum the signal contribution into the output. */
      y += ajl * R * V;
    }
  }

  /* return the computed result. */
  return y;
}

/* model_predict(): evaluate the approximate posterior predictive mean
 * and variance of a signal model from every time/phase combination in
 * a dataset.
 *
 * the datasets must have the same number of elements, and it is assumed
 * that their time/phase values are also identical.
 *
 * arguments:
 *  @mdl: pointer to the model to access.
 *  @mean: dataset to fill with mean predictions.
 *  @var: dataset to fill with variance predictions.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the evaluation succeeded.
 */
int model_predict (const model_t *mdl, dataset_t *mean, dataset_t *var) {
  /* declare required variables:
   *  @mu: individual predicted means.
   *  @eta: individual predicted variances.
   */
  double mu, eta;

  /* check the input arguments. */
  if (!mdl || !mean || !var)
    return 0;

  /* check that the structures have matching dimensionality. */
  if (mean->D != mdl->D || var->D != mdl->D)
    return 0;

  /* loop over each signal phase. */
  for (unsigned int k = 0; k < mean->K; k++) {
    /* check that the times match in size. */
    if (mean->n[k] != var->n[k])
      return 0;

    /* loop over the times of the current phase. */
    for (unsigned int i = 0; i < mean->n[k]; i++) {
      /* compute the posterior mean and variance. */
      vector_view_t t = matrix_row(mean->T[k], i);
      model_predict_single(mdl, k, &t, &mu, &eta);

      /* store the predictions. */
      vector_set(mean->y[k], i, mu);
      vector_set(var->y[k], i, eta);
    }
  }

  /* return success. */
  return 1;
}

/* model_predict_single(): compute the posterior predictive mean and
 * variance of a signal model at a specified phase and time.
 *
 * arguments:
 *  @mdl: pointer to the model to access.
 *  @k: signal phase of the prediction.
 *  @t: time of the prediction.
 *  @mean: output predicted mean.
 *  @var: output predicted variance.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int model_predict_single (const model_t *mdl,
                          const unsigned int k,
                          const vector_t *t,
                          double *mean, double *var) {
  /* declare required variables:
   *  @mu: temporary predicted mean.
   *  @eta: temporary predicted variance.
   */
  double mu, eta;

  /* check the input arguments. */
  if (!mdl || !t || k >= mdl->K || t->len != mdl->D)
    return 0.0;

  /* gain access to the dimensionalities. */
  const unsigned int K = mdl->K;
  const unsigned int M = mdl->M;

  /* initialize the mean. */
  mu = 0.0;

  /* compute the predicted mean. */
  for (unsigned int l = 0; l < K; l++) {
    for (unsigned int j = 0; j < M; j++) {
      const double a = vector_get(mdl->ahat, l * M + j);
      const double phi = expect_basis(mdl, j, k, l, t);

      mu += a * phi;
    }
  }

  /* initialize the variance. */
  eta = mdl->sigma * mdl->sigma - mu * mu;

  /* compute the rest of the predicted variance. */
  for (unsigned int l1 = 0; l1 < K; l1++) {
    for (unsigned int l2 = 0; l2 < K; l2++) {
      for (unsigned int j1 = 0; j1 < M; j1++) {
        const unsigned int i1 = l1 * M + j1;

        for (unsigned int j2 = 0; j2 < M; j2++) {
          const unsigned int i2 = l2 * M + j2;
          const double a2 = matrix_get(mdl->Sigma, i1, i2) +
                            vector_get(mdl->ahat, i1) *
                            vector_get(mdl->ahat, i2);

          const double phi2 = interact_basis(mdl, j1, j2, l1, l2, k, t);

          eta += a2 * phi2;
        }
      }
    }
  }

  /* store the computed results. */
  *mean = mu;
  *var = eta;

  /* return success. */
  return 1;
}

/* model_cov(): compute the covariance function of a signal model
 * structure. this is used by gaussian processes to compute kernel
 * matrix and vector elements.
 *
 * arguments:
 *  @mdl: pointer to the model structure to access.
 *  @k1, @k2: phase pair on which to compute covariance.
 *  @t1, @t2: time points at which to compute covariance.
 *
 * returns:
 *  covariance between the signal model on the phases @k1 and @k2
 *  at times @t1 and @t2.
 *
 *  i.e.:
 *   E[f^{k1}(t1) f^{k2}(t2)]
 */
double model_cov (const model_t *mdl,
                  const unsigned int k1, const unsigned int k2,
                  const vector_t *t1, const vector_t *t2) {
  /* @cov: output covariance value. */
  double cov = 0.0;

  /* sum together the model basis contributions. */
  for (unsigned int j = 0; j < mdl->M; j++)
    cov += cov_basis(mdl, j, k1, k2, t1, t2);

  /* scale the model contribution, and add
   * a noise contribution if necessary.
   */
  cov *= (mdl->delta * mdl->delta);
  if (vector_equal(t1, t2))
    cov += (mdl->sigma * mdl->sigma);

  /* return the computed result. */
  return cov;
}

