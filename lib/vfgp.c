
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
   *  @alpha: cholesky solution to C*alpha = y.
   */
  matrix_t *C;
  vector_t *alpha;

  /* function prediction variables:
   *  @xs: input vector of the current prediction.
   *  @ks: signal phase of the current prediction.
   *  @cs: vector of covariances for prediction.
   *  @beta: cholesky solution to C*beta = cs.
   *  @css: predictive variance upper bound.
   *  @mean: predictive mean.
   *  @var: predictive variance.
   */
  vector_t *x, *xs, *cs, *beta;
  double css, mean, var;
  unsigned int ks;
}
vfgp_t;

/* vfgp_init(): initialize a hybrid vfr/gp model.
 *  - see model_init_fn() for more information.
 */
MODEL_INIT (vfgp) {
  /* FIXME: implement vfgp_init() */

  /* return success. */
  return 1;
}

/* vfgp_bound(): return the lower bound of a hybrid vfr/gp model.
 *  - see model_bound_fn() for more information.
 */
MODEL_BOUND (vfgp) {
  /* return the fixed-tau vfr bound. */
  return model_type_tauvfr->bound(mdl);
}

/* vfgp_predict(): return the prediction of a hybrid vfr/gp model.
 *  - see model_predict_fn() for more information.
 */
MODEL_PREDICT (vfgp) {
  /* get the extended structure pointer. */
  vfgp_t *mdx = (vfgp_t*) mdl;

  /* check if the gaussian process is enabled. */
  if (mdx->gp_enable) {
    /* FIXME: implement vfgp_bound() */

    /* return success. */
    return 1;
  }

  /* fall back to the fixed-tau vfr prediction. */
  return model_type_tauvfr->predict(mdl, x, p, mean, var);
}

/* vfgp_infer(): perform complete inference in a hybrid vfr/gp model.
 *  - see model_infer_fn() for more information.
 */
MODEL_INFER (vfgp) {
  /* use fixed-tau vfr full inference. */
  return model_type_tauvfr->infer(mdl);
}

/* vfgp_update(): perform efficient low-rank inference in a
 * hybrid vfr/gp model.
 *  - see model_update_fn() for more information.
 */
MODEL_UPDATE (vfgp) {
  /* use fixed-tau vfr partial inference. */
  return model_type_tauvfr->update(mdl, j);
}

/* vfgp_gradient(): return the gradient of a single factor in a
 * hybrid vfr/gp model.
 *  - see model_gradient_fn() for more information.
 */
MODEL_GRADIENT (vfgp) {
  /* return the fixed-tau vfr gradient. */
  return model_type_tauvfr->gradient(mdl, i, j, grad);
}

/* vfgp_meanfield(): return the coefficients required for an
 * assumed-density mean-field update of a factor in a
 * hybrid vfr/gp model.
 *  - see model_meanfield_fn() for more information.
 */
MODEL_MEANFIELD (vfgp) {
  /* return the fixed-tau vfr mean-field coefficients. */
  return model_type_tauvfr->meanfield(mdl, j, A, B);
}

/* vfgp_set_mode(): set the prediction mode employed by a variational
 * feature gaussian process model.
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

  /* FIXME: implement vfgp_set_mode() */

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

