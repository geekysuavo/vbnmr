
/* include the variational feature gaussian process header. */
#include <vbnmr/vfgp.h>

/* vfgp_t: structure for holding a hybrid gaussian process,
 * fixed-tau vfr model.
 */
typedef struct {
  /* @base: core model structure members. */
  model_t base;

  /* prediction mode:
   *  @gp_enable: whether or not to use gaussian process predictions.
   *  @gp_ready: whether or not the gaussian process is ready for use.
   */
  unsigned int gp_enable, gp_ready;

  /* data fitting variables:
   *  @C: kernel matrix (cholesky factors) for the current dataset.
   *  @y: vector of all observed function values.
   *  @alpha: cholesky solution to C*alpha = y.
   */
  matrix_t *C;
  vector_t *y, *alpha;

  /* function prediction variables:
   *  @cs: vector of covariances for prediction.
   *  @beta: cholesky solution to C*beta = cs.
   */
  vector_t *cs, *beta;
}
vfgp_t;

/* --- */

/* vfgp_gpinit(): initialize all fit and prediction parameters
 * that are used by the gaussian process portion of a vfgp model.
 *
 * arguments:
 *  @mdl: model structure pointer.
 */
static inline void vfgp_gpinit (const model_t *mdl) {
  /* get the extended structure pointer. */
  vfgp_t *mdx = (vfgp_t*) mdl;

  /* initialize the fit parameters. */
  mdx->C = NULL;
  mdx->y = NULL;
  mdx->alpha = NULL;

  /* initialize the prediction parameters. */
  mdx->cs = NULL;
  mdx->beta = NULL;
}

/* vfgp_gpfit(): re-fit the gaussian process portion of a vfgp model.
 *
 * arguments:
 *  @mdl: model structure pointer.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
static inline int vfgp_gpfit (const model_t *mdl) {
  /* get the extended structure pointer and declare data. */
  vfgp_t *mdx = (vfgp_t*) mdl;
  datum_t *d1, *d2;

  /* return if the gp is already ready. */
  if (mdx->gp_ready)
    return 1;

  /* determine the number of observations. */
  const unsigned int N = (mdl->dat ? mdl->dat->N : 0);

  /* free the fit variables. */
  matrix_free(mdx->C);
  vector_free(mdx->y);
  vector_free(mdx->alpha);

  /* free the prediction variables. */
  vector_free(mdx->cs);
  vector_free(mdx->beta);

  /* re-initialize the gp variables. */
  vfgp_gpinit(mdl);

  /* allocate the fit variables. */
  mdx->C = matrix_alloc(N, N);
  mdx->y = vector_alloc(N);
  mdx->alpha = vector_alloc(N);

  /* allocate the prediction variables. */
  mdx->cs = vector_alloc(N);
  mdx->beta = vector_alloc(N);

  /* check for allocation failures. */
  if (!mdx->C || !mdx->y || !mdx->alpha || !mdx->cs || !mdx->beta)
    return 0;

  /* loop over the covariance matrix rows. */
  for (unsigned int i1 = 0; i1 < N; i1++) {
    /* get the row-wise observation. */
    d1 = data_get(mdl->dat, i1);
    vector_set(mdx->y, i1, d1->y);

    /* loop over the covariance matrix columns. */
    for (unsigned int i2 = i1; i2 < N; i2++) {
      /* get the column-wise observation. */
      d2 = data_get(mdl->dat, i2);

      /* compute the covariance matrix element. */
      const double c12 = model_cov(mdl, d1->x, d2->x, d1->p, d2->p);

      /* store the computed matrix element. */
      matrix_set(mdx->C, i1, i2, c12);
      if (i1 != i2)
        matrix_set(mdx->C, i2, i1, c12);
    }
  }

  /* decompose the covariance matrix. */
  if (!chol_decomp(mdx->C))
    return 0;

  /* solve to obtain the alpha vector. */
  chol_solve(mdx->C, mdx->y, mdx->alpha);

  /* incidate readiness and return success. */
  mdx->gp_ready = 1;
  return 1;
}

/* vfgp_gppred(): compute the posterior predictive mean and variance
 * of the gaussian process portion of a vfgp model.
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @x: observation input vector.
 *  @p: function output index.
 *  @mean, @var: output predictions.
 */
static inline void vfgp_gppred (const model_t *mdl,
                                const vector_t *x,
                                const unsigned int p,
                                double *mean,
                                double *var) {
  /* get the extended structure pointer. */
  vfgp_t *mdx = (vfgp_t*) mdl;

  /* loop over the observations. */
  const unsigned int N = (mdl->dat ? mdl->dat->N : 0);
  for (unsigned int i = 0; i < N; i++) {
    /* get the current observation. */
    datum_t *di = data_get(mdl->dat, i);

    /* compute the prediction vector element. */
    const double ci = model_cov(mdl, di->x, x, di->p, p);
    vector_set(mdx->cs, i, ci);
  }

  /* solve for the beta vector. */
  chol_solve(mdx->C, mdx->cs, mdx->beta);

  /* compute the predictive mean. */
  *mean = blas_ddot(mdx->cs, mdx->alpha);

  /* compute the predictive variance. */
  *var = model_cov(mdl, x, x, p, p) - blas_ddot(mdx->cs, mdx->beta);
}

/* --- */

/* vfgp_init(): initialize a hybrid vfr/gp model.
 *  - see model_init_fn() for more information.
 */
MODEL_INIT (vfgp) {
  /* get the extended structure pointer. */
  vfgp_t *mdx = (vfgp_t*) mdl;

  /* initialize the gaussian process. */
  mdx->gp_enable = 0;
  mdx->gp_ready = 0;
  vfgp_gpinit(mdl);

  /* return success. */
  return 1;
}

/* vfgp_bound(): return the lower bound of a hybrid vfr/gp model.
 *  - see model_bound_fn() for more information.
 */
MODEL_BOUND (vfgp) {
  /* return the fixed-tau vfr bound. */
  return vfl_model_tauvfr->bound(mdl);
}

/* vfgp_predict(): return the prediction of a hybrid vfr/gp model.
 *  - see model_predict_fn() for more information.
 */
MODEL_PREDICT (vfgp) {
  /* get the extended structure pointer. */
  vfgp_t *mdx = (vfgp_t*) mdl;

  /* check if the gaussian process is enabled. */
  if (mdx->gp_enable) {
    /* re-fit the gaussian process. */
    if (!vfgp_gpfit(mdl))
      return 0;

    /* perform gaussian process prediction. */
    vfgp_gppred(mdl, x, p, mean, var);
    return 1;
  }

  /* fall back to the fixed-tau vfr prediction. */
  return vfl_model_tauvfr->predict(mdl, x, p, mean, var);
}

/* vfgp_infer(): perform complete inference in a hybrid vfr/gp model.
 *  - see model_infer_fn() for more information.
 */
MODEL_INFER (vfgp) {
  /* use fixed-tau vfr full inference. */
  return vfl_model_tauvfr->infer(mdl);
}

/* vfgp_update(): perform efficient low-rank inference in a
 * hybrid vfr/gp model.
 *  - see model_update_fn() for more information.
 */
MODEL_UPDATE (vfgp) {
  /* use fixed-tau vfr partial inference. */
  return vfl_model_tauvfr->update(mdl, j);
}

/* vfgp_gradient(): return the gradient of a single factor in a
 * hybrid vfr/gp model.
 *  - see model_gradient_fn() for more information.
 */
MODEL_GRADIENT (vfgp) {
  /* return the fixed-tau vfr gradient. */
  return vfl_model_tauvfr->gradient(mdl, i, j, grad);
}

/* vfgp_meanfield(): return the coefficients required for an
 * assumed-density mean-field update of a factor in a
 * hybrid vfr/gp model.
 *  - see model_meanfield_fn() for more information.
 */
MODEL_MEANFIELD (vfgp) {
  /* return the fixed-tau vfr mean-field coefficients. */
  return vfl_model_tauvfr->meanfield(mdl, i, j, b, B);
}

/* vfgp_set_mode(): set the prediction mode employed by a variational
 * feature gaussian process model.
 *
 * this function also clears the readiness flag of the gaussian process,
 * so it can be used to force a re-fit before new predictions are made.
 *
 * arguments:
 *  @mdl: model structure pointer to modify.
 *  @gp_enable: whether or not to enable the gp.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int vfgp_set_mode (model_t *mdl, const unsigned int gp_enable) {
  /* check the structure pointer. */
  if (!mdl)
    return 0;

  /* get the extended structure pointer and set the gp status. */
  vfgp_t *mdx = (vfgp_t*) mdl;
  mdx->gp_enable = gp_enable;
  mdx->gp_ready = 0;

  /* return success. */
  return 1;
}

/* vfgp_type: variational feature gaussian process type structure.
 */
static model_type_t vfgp_type = {
  "vfgp",                                        /* name      */
  sizeof(vfgp_t),                                /* size      */
  vfgp_init,                                     /* init      */
  vfgp_bound,                                    /* bound     */
  vfgp_predict,                                  /* predict   */
  vfgp_infer,                                    /* infer     */
  vfgp_update,                                   /* update    */
  vfgp_gradient,                                 /* gradient  */
  vfgp_meanfield                                 /* meanfield */
};

/* vbnmr_model_vfgp: address of the vfgp_type structure. */
const model_type_t *vbnmr_model_vfgp = &vfgp_type;

