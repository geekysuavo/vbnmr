
/* include the gaussian process header. */
#include <vbnmr/gp.h>

/* gp_reinit(): initialize all fitting-related pointers within
 * a gaussian process to null.
 *
 * arguments:
 *  @gp: pointer to the gp structure to initialize.
 */
static inline void gp_reinit (gp_t *gp) {
  /* initialize the fitting parameters. */
  gp->C = NULL;
  gp->alpha = NULL;

  /* initialize the predictive parameters. */
  gp->t = gp->ts = gp->cs = gp->beta = NULL;
}

/* gp_init(): initialize the pointers within a gaussian process
 * structures to null.
 *
 * arguments:
 *  @gp: pointer to the gp structure to initialize.
 */
static inline void gp_init (gp_t *gp) {
  /* initialize the size of the signal, data, and model. */
  gp->D = gp->K = gp->M = gp->N = 0;

  /* initialize the model and dataset structure pointers. */
  gp->mdl = NULL;
  gp->dat = NULL;

  /* initialize the fitting and predictive parameters. */
  gp_reinit(gp);
}

/* gp_alloc(): allocate a new gaussian process for use.
 *
 * arguments:
 *  @dat: pointer to a dataset to associate with the gp.
 *  @mdl: pointer to a model to associate with the gp.
 *
 * returns:
 *  newly allocated and initialized gaussian process structure pointer.
 *  if either @dat or @mdl is null, or if any allocations fail, a null
 *  pointer is returned.
 */
gp_t *gp_alloc (dataset_t *dat, model_t *mdl) {
  /* @gp: pointer to the new gaussian process structure.
   */
  gp_t *gp;

  /* fail if either passed pointer is null. */
  if (!dat || !mdl)
    return NULL;

  /* allocate a new structure pointer, or fail. */
  gp = (gp_t*) malloc(sizeof(gp_t));
  if (!gp)
    return NULL;

  /* initialize the gp dataset and model pointers. */
  gp_init(gp);
  gp->dat = dat;
  gp->mdl = mdl;

  /* store the current signal, model and data sizes. */
  gp->D = gp->mdl->D;
  gp->K = gp->mdl->K;
  gp->M = gp->mdl->M;
  gp->N = gp->dat->N;

  /* return the new structure pointer. */
  return gp;
}

/* gp_free(): free an allocated gaussian process.
 *
 * arguments:
 *  @gp: gaussian process structure pointer to free.
 */
void gp_free (gp_t *gp) {
  /* return if the structure pointer is null. */
  if (!gp) return;

  /* free the covariance matrix and data projection. */
  matrix_free(gp->C);
  vector_free(gp->alpha);

  /* free the structures used for prediction. */
  vector_free(gp->t);
  vector_free(gp->ts);
  vector_free(gp->cs);
  vector_free(gp->beta);

  /* free the gp structure pointer. */
  gp_init(gp);
  free(gp);
}

/* gp_fit(): re-fit a gaussian process to its associated dataset,
 * given the current set of associated model parameters.
 *
 * arguments:
 *  @gp: pointer to the gaussian process to update.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int gp_fit (gp_t *gp) {
  /* @t1, @t2: times used for kernel computations.
   * @yv: measurements used for kernel computations.
   * @iv, @kv: measurement and phase indices.
   * @Cij: current kernel value.
   */
  unsigned int *iv, *kv, status = 0;
  vector_t *t1, *t2, *yv;
  double Cij;

  /* fail if the structure pointer is null. */
  if (!gp) return 0;

  /* update the data size. */
  gp->N = gp->dat->N;

  /* free the fitting variables. */
  matrix_free(gp->C);
  vector_free(gp->alpha);

  /* free the prediction variables. */
  vector_free(gp->t);
  vector_free(gp->ts);
  vector_free(gp->cs);
  vector_free(gp->beta);

  /* re-initialize the data structure. */
  gp_reinit(gp);

  /* initialize the temporary variables. */
  iv = kv = NULL;
  yv = t1 = t2 = NULL;

  /* allocate the temporary variables. */
  iv = (unsigned int*) malloc(gp->N * sizeof(unsigned int));
  kv = (unsigned int*) malloc(gp->N * sizeof(unsigned int));
  yv = vector_alloc(gp->N);
  t1 = vector_alloc(gp->D);
  t2 = vector_alloc(gp->D);

  /* check if any allocation failed. */
  if (!iv || !kv || !yv || !t1 || !t2)
    goto fail;

  /* compute vectorized versions of the measurements. */
  for (unsigned int k = 0, i = 0; k < gp->K; k++) {
    for (unsigned int ik = 0; ik < gp->dat->n[k]; ik++, i++) {
      vector_set(yv, i, vector_get(gp->dat->y[k], ik));
      iv[i] = ik;
      kv[i] = k;
    }
  }

  /* allocate the fitting variables. */
  gp->C = matrix_alloc(gp->N, gp->N);
  gp->alpha = vector_alloc(gp->N);

  /* allocate the prediction variables. */
  gp->t = vector_alloc(gp->D);
  gp->ts = vector_alloc(gp->D);
  gp->cs = vector_alloc(gp->N);
  gp->beta = vector_alloc(gp->N);

  /* check if any allocation failed. */
  if (!gp->C || !gp->alpha || !gp->t || !gp->ts || !gp->cs || !gp->beta)
    goto fail;

  /* loop over the covariance matrix rows. */
  for (unsigned int i1 = 0; i1 < gp->N; i1++) {
    /* get the time for the current row. */
    matrix_copy_row(t1, gp->dat->T[kv[i1]], iv[i1]);

    /* loop over the covariance matrix columns. */
    for (unsigned int i2 = i1; i2 < gp->N; i2++) {
      /* get the time for the current column. */
      matrix_copy_row(t2, gp->dat->T[kv[i2]], iv[i2]);

      /* compute the current covariance matrix element. */
      Cij = model_cov(gp->mdl, kv[i1], kv[i2], t1, t2);

      /* store the covariance matrix element. */
      matrix_set(gp->C, i1, i2, Cij);
      if (i1 != i2)
        matrix_set(gp->C, i2, i1, Cij);
    }
  }

  /* compute the inverse covariance matrix. */
  if (!cholesky_decomp(gp->C))
    goto fail;

  /* multiply the inverse covariance matrix against the data. */
  cholesky_solve(gp->C, yv, gp->alpha);

  /* all operations were successful. */
  status = 1;

fail:
  /* free any allocated memory. */
  vector_free(t1);
  vector_free(t2);
  vector_free(yv);
  free(iv);
  free(kv);

  /* return the final status code. */
  return status;
}

/* gp_predict(): compute the posterior predictive mean and variance
 * of a gaussian process at a specified phase and time.
 *
 * arguments:
 *  @gp: pointer to the gaussian process to access.
 *  @ks: signal phase of the prediction.
 *  @ts: time of the prediction.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int gp_predict (gp_t *gp, const unsigned int ks, const vector_t *ts) {
  /* @Ci: current prediction kernel vector value.
   */
  double Ci;

  /* fail if:
   *  - any passed pointers are null.
   *  - the phase index is out of bounds.
   *  - the time vector length does not match the dataset.
   */
  if (!gp || !ts || ks >= gp->K || ts->len != gp->D)
    return 0;

  /* initialize the prediction results. */
  gp->css = gp->mean = gp->var = 0.0;
  vector_copy(gp->ts, ts);
  gp->ks = ks;

  /* loop over the measurements. */
  for (unsigned int k = 0, i = 0; k < gp->K; k++) {
    for (unsigned int ik = 0; ik < gp->dat->n[k]; ik++, i++) {
      /* get the current measurement time. */
      matrix_copy_row(gp->t, gp->dat->T[k], ik);

      /* compute the prediction kernel vector element. */
      Ci = model_cov(gp->mdl, k, ks, gp->t, ts);
      vector_set(gp->cs, i, Ci);
    }
  }

  /* multiply the inverse kernel matrix by the prediction vector. */
  cholesky_solve(gp->C, gp->cs, gp->beta);

  /* compute the predictive mean. */
  gp->mean = blas_ddot(gp->cs, gp->alpha);

  /* compute the predictive variance. */
  gp->css = model_cov(gp->mdl, ks, ks, ts, ts);
  gp->var = blas_ddot(gp->cs, gp->beta);
  gp->var = gp->css - gp->var;

  /* return success. */
  return 1;
}

/* gp_eval(): evaluate the posterior predictive mean and variance
 * of a gaussian process for every time/phase combination in a dataset.
 *
 * the datasets must have the same number of elements, and it is
 * assumed that their time/phase values are also identical.
 *
 * arguments:
 *  @gp: pointer to the gaussian process to access.
 *  @mean: dataset to fill with mean predictions.
 *  @var: dataset to fill with variance predictions.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the evaluation succeeded.
 */
int gp_eval (gp_t *gp, dataset_t *mean, dataset_t *var) {
  /* check the input arguments. */
  if (!gp || !mean || !var)
    return 0;

  /* check that the structures have matching dimensionality. */
  if (mean->D != gp->D || var->D != gp->D)
    return 0;

  /* loop over each signal phase. */
  for (unsigned int k = 0; k < mean->K; k++) {
    /* check that the times match in size. */
    if (mean->n[k] != var->n[k])
      return 0;

    /* loop over the times of the current phase. */
    for (unsigned int i = 0; i < mean->n[k]; i++) {
      /* compute the posterior mean and variance the current time. */
      vector_view_t t = matrix_row(mean->T[k], i);
      gp_predict(gp, k, &t);

      /* store the predictions. */
      vector_set(mean->y[k], i, gp->mean);
      vector_set(var->y[k], i, gp->var);
    }
  }

  /* return success. */
  return 1;
}

