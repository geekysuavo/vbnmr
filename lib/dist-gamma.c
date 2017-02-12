
/* include the gamma distribution header. */
#include <vbnmr/dist-gamma.h>

/* * * * univariate expectations: * * * */

/* gam_decay(): compute the expectation of an exponential decay under
 * a gamma distribution.
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: signal basis index.
 *  @d: signal dimension.
 *  @t: time value.
 *
 * returns:
 *  E[exp(-x*t)]_G(x | alpha_jd, beta_jd)
 */
inline double
gam_decay (const model_t *mdl, const unsigned int j,
           const unsigned int d, const double t) {
  /* gain access to the parameters. */
  const double alpha = mdl->sig[j][d].alpha;
  const double beta = mdl->sig[j][d].beta;

  /* compute and return the expectation. */
  return pow(beta / (beta + t), alpha);
}

/* gam_entropy(): compute the entropy of a gamma distribution,
 * up to a constant term.
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: signal basis index.
 *  @d: signal dimension.
 *
 * returns:
 *  H[G(alpha_jd, beta_jd)]
 */
inline double
gam_entropy (const model_t *mdl,
             const unsigned int j,
             const unsigned int d) {
  /* gain access to the parameters. */
  const double alpha = mdl->sig[j][d].alpha;
  const double beta = mdl->sig[j][d].beta;

  /* gain access to functions of the parameters. */
  const double dga = mdl->sig[j][d].dga;
  const double lga = mdl->sig[j][d].lga;

  /* compute and return the entropy of the distribution. */
  return alpha - log(beta) + lga + (1.0 - alpha) * dga;
}

/* gam_prior(): compute the expectation of a gamma log-prior with
 * respect to its associated variational approximation.
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: signal basis index.
 *  @d: signal dimension.
 *
 * returns:
 *  E[ln G(x | alpha0, beta0)]_G(x | alpha_jd, beta_jd)
 */
inline double
gam_prior (const model_t *mdl,
           const unsigned int j,
           const unsigned int d) {
  /* gain access to the prior parameters. */
  const double alpha0 = mdl->sig[j][d].alpha0;
  const double beta0 = mdl->sig[j][d].beta0;

  /* gain access to the variational parameters. */
  const double alpha = mdl->sig[j][d].alpha;
  const double beta = mdl->sig[j][d].beta;

  /* gain access to functions of the parameters. */
  const double dga = mdl->sig[j][d].dga;

  /* compute and return the prior term. */
  return (alpha0 - 1.0) * (dga - log(beta))
       - alpha * (beta0 / beta);
}

/* * * * multivariate expectations: * * * */

/* expect_decay(): compute the expectation of a multidimensional
 * decay term, with an option to exclude a single dimension.
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: signal basis index.
 *  @delta: excluded dimension index.
 *  @t: multidimensional time vector.
 *
 * returns:
 *  E[R(t | rho_j \ rho_jd)]
 */
double expect_decay (const model_t *mdl, const unsigned int j,
                     const unsigned int delta, const vector_t *t) {
  /* @expect: output expectation value. */
  double expect = 1.0;

  /* multiply together the expectations for each dimension. */
  for (unsigned int d = 0; d < mdl->D; d++) {
    /* skip the excluded dimension. */
    if (d == delta) continue;

    /* update the expectation. */
    const double td = vector_get(t, d);
    expect *= gam_decay(mdl, j, d, td);
  }

  /* return the computed expectation. */
  return expect;
}

/* expect_all_decay(): compute a full vector of expectations of
 * a multidimensional decay term.
 *
 * operation:
 *  R(d) <- R(d) * < exp(-rho_jd * t_d) >
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: signal basis index.
 *  @t: multidimensional time vector.
 *  @R: output factor vector.
 *
 * returns:
 *  < R(t | rho_j) >
 */
double expect_all_decay (const model_t *mdl, const unsigned int j,
                         const vector_t *t, vector_t *R) {
  /* declare required variables:
   *  @Rd: storage for each factor.
   *  @Rall: product of all factors.
   */
  double td, Rd, Rall = 1.0;

  /* compute the factor for each dimension. */
  for (unsigned int d = 0; d < mdl->D; d++) {
    /* update the product. */
    td = vector_get(t, d);
    Rd = gam_decay(mdl, j, d, td);
    Rall *= Rd;

    /* update the vector. */
    vector_set(R, d, vector_get(R, d) * Rd);
  }

  /* return the product of all factors. */
  return Rall;
}

/* * * * multivariate interactions: * * * */

/* interact_decay(): compute the interaction of two multidimensional
 * decay terms, with an option to exclude a single dimension.
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j1: first signal basis index.
 *  @j2: second signal basis index.
 *  @delta: excluded dimension index.
 *  @t: multidimensional time vector.
 *
 * returns:
 *  E[R(t | rho_j1 \ rho_j1,d) R(t | rho_j2 \ rho_j1,d)]
 */
double interact_decay (const model_t *mdl,
                       const unsigned int j1, const unsigned int j2,
                       const unsigned int delta, const vector_t *t) {
  /* @g: output interaction value. */
  double g = 1.0;

  /* for different signal indices, the interaction factors into
   * the product of the two independent expectations.
   */
  if (j1 != j2)
    return expect_decay(mdl, j1, delta, t) *
           expect_decay(mdl, j2, delta, t);

  /* multiply together the interactions along each dimension. */
  for (unsigned int d = 0; d < mdl->D; d++) {
    /* skip the excluded dimension. */
    if (d == delta) continue;

    /* update the interaction. */
    const double td = 2.0 * vector_get(t, d);
    g *= gam_decay(mdl, j1, d, td);
  }

  /* return the computed interaction. */
  return g;
}

/* interact_all_decay(): compute a full vector of interactions of
 * two multidimensional decay terms.
 *
 * operation:
 *  R(d) <- R(d) * < exp(-rho_j1,d * t_d) exp(-rho_j2,d * t_d) >
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j1: first signal basis index.
 *  @j2: second signal basis index.
 *  @t: multidimensional time vector.
 *  @R: output factor vector.
 *
 * returns:
 *  < R(t | rho_j1) R(t | rho_j2) >
 */
double interact_all_decay (const model_t *mdl,
                           const unsigned int j1, const unsigned int j2,
                           const vector_t *t, vector_t *R) {
  /* declare required variables:
   *  @Rd: storage for each factor.
   *  @Rprod: product of all factors.
   */
  double td, Rd, Rall = 1.0;

  /* for different signal indices, the interaction factors into
   * the product of the two independent expectations.
   */
  if (j1 != j2)
    return expect_all_decay(mdl, j1, t, R) *
           expect_all_decay(mdl, j2, t, R);

  /* compute the factor for each dimension. */
  for (unsigned int d = 0; d < mdl->D; d++) {
    /* update the product. */
    td = 2.0 * vector_get(t, d);
    Rd = gam_decay(mdl, j1, d, td);
    Rall *= Rd;

    /* update the vector. */
    vector_set(R, d, vector_get(R, d) * Rd);
  }

  /* return the product of all factors. */
  return Rall;
}

/* * * * multivariate covariances: * * * */

/* cov_decay(): compute the covariance at different times of two decay
 * terms with shared rate variables rho(j).
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: signal basis index.
 *  @t1: first time vector.
 *  @t2: second time vector.
 *
 * returns:
 *  E[R(t1 | rho_j) R(t2 | rho_j)]
 */
double cov_decay (const model_t *mdl, const unsigned int j,
                  const vector_t *t1, const vector_t *t2) {
  /* @cov: output covariance value. */
  double cov = 1.0;

  /* multiply together the covariances along each dimension. */
  for (unsigned int d = 0; d < mdl->D; d++) {
    /* update the covariance. */
    const double td = vector_get(t1, d) + vector_get(t2, d);
    cov *= gam_decay(mdl, j, d, td);
  }

  /* return the computed covariance. */
  return cov;
}

/* * * * univariate partial derivatives: * * * */

/* diff_gam_decay(): compute the partial derivatives of the expectation
 * of an exponential decay under a gamma distribution.
 *
 * i.e.:
 *  d_alpha E[exp(-x*t)]_G(x | alpha_jd, beta_jd)
 *  d_beta  E[exp(-x*t)]_G(x | alpha_jd, beta_jd)
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: signal basis index.
 *  @d: signal dimension.
 *  @t: time value.
 *  @Dalpha: pointer to the output alpha-derivative.
 *  @Dbeta: pointer to the output beta-derivative.
 */
void diff_gam_decay (const model_t *mdl, const unsigned int j,
                     const unsigned int d, const double t,
                     double *Dalpha, double *Dbeta) {
  /* extract relevant parameters. */
  const double alpha = mdl->sig[j][d].alpha;
  const double beta = mdl->sig[j][d].beta;

  /* compute intermediate quantities. */
  const double x = beta / (beta + t);
  const double y = (alpha * t) / (beta * beta);

  /* compute and return the partial derivatives. */
  *Dalpha = pow(x, alpha) * log(x);
  *Dbeta = y * pow(x, alpha - 1.0);
}

/* * * * multivariate partial derivatives: * * * */

/* diff_expect_decay(): compute all partial derivatives of the expectation
 * of a multidimensional decay term.
 *
 * operation:
 *  Da <- Da + scale * d_alpha < R(t | rho_j) >
 *  Db <- Db + scale * d_beta  < R(t | rho_j) >
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: signal basis index.
 *  @t: multidimensional time vector.
 *  @scale: summation scale factor.
 *  @R: multidimensional factor vector.
 *  @Da: output alpha-derivative vector.
 *  @Db: output beta-derivative vector.
 */
void diff_expect_decay (const model_t *mdl, const unsigned int j,
                        const vector_t *t, const double scale,
                        const vector_t *R,
                        vector_t *Da,
                        vector_t *Db) {
  /* gain access to the model dimension count. */
  const unsigned int D = mdl->D;

  /* compute the final derivatives. */
  for (unsigned int d = 0; d < D; d++) {
    /* extract the relevant quantities. */
    const double td = vector_get(t, d);

    /* compute the partial derivative factor. */
    double dalpha, dbeta;
    diff_gam_decay(mdl, j, d, td, &dalpha, &dbeta);

    /* multiply in the other expectation factors. */
    for (unsigned int delta = 0; delta < D; delta++) {
      if (delta == d) continue;
      dalpha *= vector_get(R, delta);
      dbeta *= vector_get(R, delta);
    }

    /* store the final derivatives. */
    vector_set(Da, d, vector_get(Da, d) + scale * dalpha);
    vector_set(Db, d, vector_get(Db, d) + scale * dbeta);
  }
}

/* diff_interact_decay(): compute all partial derivatives of the interaction
 * of two multidimensional decay terms.
 *
 * operation:
 *  Da <- Da + scale * d_alpha < R(t | rho_j1) R(t | rho_j2) >
 *  Db <- Db + scale * d_beta  < R(t | rho_j1) R(t | rho_j2) >
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j1: first signal basis index.
 *  @j2: second signal basis index.
 *  @t: multidimensional time vector.
 *  @scale: summation scale factor.
 *  @R: intermediate storage vector.
 *  @Da: output alpha-derivative vector.
 *  @Db: output beta-derivative vector.
 */
void diff_interact_decay (const model_t *mdl,
                          const unsigned int j1, const unsigned int j2,
                          const vector_t *t, const double scale,
                          const vector_t *R,
                          vector_t *Da,
                          vector_t *Db) {
  /* gain access to the model dimension count. */
  const unsigned int D = mdl->D;

  /* for different signal indices, the interaction factors into
   * the product of the two independent expectations.
   */
  if (j1 != j2) {
    /* compute the partial derivatives of the first factor,
     * appropriately scaled to include the second expectation.
     */
    const double scale2 = scale * expect_decay(mdl, j2, D, t);
    diff_expect_decay(mdl, j1, t, scale2, R, Da, Db);
  }

  /* compute the final derivatives. */
  for (unsigned int d = 0; d < D; d++) {
    /* extract the relevant quantities. */
    const double td = 2.0 * vector_get(t, d);

    /* compute the partial derivative factor. */
    double dalpha, dbeta;
    diff_gam_decay(mdl, j1, d, td, &dalpha, &dbeta);

    /* multiply in the other expectation factors. */
    for (unsigned int delta = 0; delta < D; delta++) {
      if (delta == d) continue;
      dalpha *= vector_get(R, delta);
      dbeta *= vector_get(R, delta);
    }

    /* store the final derivatives. */
    vector_set(Da, d, vector_get(Da, d) + scale * dalpha);
    vector_set(Db, d, vector_get(Db, d) + scale * dbeta);
  }
}

/* diff_naturalize_decay(): convert a vector of gradients of a single
 * model signal decay component into their natural forms.
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: first signal basis index.
 *  @Da: input/output alpha-derivative vector.
 *  @Db: input/output beta-derivative vector.
 *
 * returns:
 *  smallest eigenvalue of the fisher matrix
 *  used to naturalize the gradient.
 */
double diff_naturalize_decay (const model_t *mdl, const unsigned int j,
                              vector_t *Da, vector_t *Db) {
  /* initialize the minimum eigenvalue. */
  double evmin = 1.0e+9;

  /* loop over each dimension. */
  for (unsigned int d = 0; d < mdl->D; d++) {
    /* get the values of the gradient. */
    const double g1 = vector_get(Da, d);
    const double g2 = vector_get(Db, d);

    /* get the fisher matrix elements. */
    const double F11 = mdl->sig[j][d].Frho[0];
    const double F22 = mdl->sig[j][d].Frho[1];
    const double F12 = mdl->sig[j][d].Frho[2];

    /* compute the determinant. */
    const double Fdet = F11 * F22 - F12 * F12;

    /* compute the inverse fisher matrix elements. */
    const double G11 =  F22 / Fdet;
    const double G22 =  F11 / Fdet;
    const double G12 = -F12 / Fdet;

    /* compute the natural gradient values. */
    const double n1 = G11 * g1 + G12 * g2;
    const double n2 = G12 * g1 + G22 * g2;

    /* store the naturalized values. */
    vector_set(Da, d, n1);
    vector_set(Db, d, n2);

    /* compute the current eigenvalue. */
    const double B = -F11 - F22;
    const double ev = -0.5 * (B + sqrt(B * B - 4.0 * Fdet));

    /* update the minimum eigenvalue. */
    if (ev < evmin) evmin = ev;
  }

  /* return the minimum eigenvalue. */
  return evmin;
}

