
/* include the normal distribution header. */
#include <vbnmr/dist-normal.h>

/* * * * univariate expectations: * * * */

/* nrm_cos(): compute the expectation of a phased, stretched cosine
 * function under a normal distribution.
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: signal basis index.
 *  @d: signal dimension.
 *  @t: time shift.
 *  @z: integer phase offset.
 *
 * returns:
 *  E[cos(x*t+z*pi/2)]_N(x | mu_jd, tau_jd)
 */
inline double
nrm_cos (const model_t *mdl,
         const unsigned int j,
         const unsigned int d,
         const double t, const int z) {
  /* compute intermediate quantities. */
  const double mu = mdl->sig[j][d].mu;
  const double tau = mdl->sig[j][d].tau;

  /* compute and return the expectation. */
  return exp(-0.5 * t * t / tau) * cos(mu * t + M_PI_2 * (double) z);
}

/* nrm_cospair(): compute the expectation of a product of differently
 * phased and stretched cosine functions under a normal distribution.
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: signal basis index.
 *  @d: signal dimension.
 *  @t1, @t2: time shifts.
 *  @z1, @z2: integer phase offsets.
 *
 * returns:
 *  E[cos(x*t1+z1*pi/2)*cos(x*t2+z2*pi/2)]_N(x | mu_jd, tau_jd)
 */
inline double
nrm_cospair (const model_t *mdl,
             const unsigned int j,
             const unsigned int d,
             const double t1, const double t2,
             const int z1, const int z2) {
  /* compute the sum and difference of the phase offsets. */
  const int zp = z1 + z2;
  const int zm = z1 - z2;

  /* compute the sum and difference of the times. */
  const double tp = t1 + t2;
  const double tm = t1 - t2;

  /* compute and return the expectation. */
  return 0.5 * (nrm_cos(mdl, j, d, tp, zp) + nrm_cos(mdl, j, d, tm, zm));
}

/* nrm_entropy(): compute the entropy of a normal distribution,
 * up to a constant term.
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: signal basis index.
 *  @d: signal dimension.
 *
 * returns:
 *  H[N(mu_jd, tau_jd)]
 */
inline double
nrm_entropy (const model_t *mdl,
             const unsigned int j,
             const unsigned int d) {
  /* compute and return the entropy of the distribution. */
  return -0.5 * log(mdl->sig[j][d].tau);
}

/* nrm_prior(): compute the expectation of a normal log-prior with
 * respect to its associated variational approximation.
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: signal basis index.
 *  @d: signal dimension.
 *
 * returns:
 *  E[ln N(x | mu0, tau0)]_N(x | mu_jd, tau_jd)
 */
inline double
nrm_prior (const model_t *mdl,
           const unsigned int j,
           const unsigned int d) {
  /* gain access to the prior parameters. */
  const double mu0 = mdl->sig[j][d].mu0;
  const double tau0 = mdl->sig[j][d].tau0;

  /* gain access to the variational parameters. */
  const double mu = mdl->sig[j][d].mu;
  const double tau = mdl->sig[j][d].tau;

  /* compute and return the prior term. */
  return -0.5 * tau0 * (mu * mu + 1.0 / tau - 2.0 * mu * mu0);
}

/* * * * multivariate expectations: * * * */

/* expect_freq(): compute the expectation of a multidimensional
 * frequency modulation term, with an option to exclude a single
 * dimension.
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: signal basis index.
 *  @k: signal phase index.
 *  @l: signal term index.
 *  @delta: excluded dimension index.
 *  @t: multidimensional time vector.
 *
 * returns:
 *  E[V_kl(t | omega_j \ omega_jd)]
 */
double expect_freq (const model_t *mdl, const unsigned int j,
                    const unsigned int k, const unsigned int l,
                    const unsigned int delta, const vector_t *t) {
  /* @expect: output expectation value. */
  double expect = 1.0;

  /* multiply together the expectations for each dimension. */
  for (unsigned int d = 0; d < mdl->D; d++) {
    /* skip the excluded dimension. */
    if (d == delta) continue;

    /* extract the relevant quantities. */
    const double td = vector_get(t, d);

    /* update the expectation. */
    expect *= nrm_cos(mdl, j, d, td, VTAB(k, l, d));
  }

  /* return the computed expectation. */
  return expect;
}

/* expect_all_freq(): compute a full vector of expectations of
 * a multidimensional frequency modulation term.
 *
 * operation:
 *  V(d) <- V(d) * < cos(omega_jd * t_d + pi/2 * v_{k,l,d}) >
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: signal basis index.
 *  @k: signal phase index.
 *  @l: signal term index.
 *  @t: multidimensional time vector.
 *  @V: output factor vector.
 *
 * returns:
 *  < V_{k,l}(t | omega_j) >
 */
double expect_all_freq (const model_t *mdl, const unsigned int j,
                        const unsigned int k, const unsigned int l,
                        const vector_t *t, vector_t *V) {
  /* declare required variables:
   *  @Vd: storage for each factor.
   *  @Vall: product of all factors.
   */
  double td, Vd, Vall = 1.0;

  /* compute the factor for each dimension. */
  for (unsigned int d = 0; d < mdl->D; d++) {
    /* update the product. */
    td = vector_get(t, d);
    Vd = nrm_cos(mdl, j, d, td, VTAB(k, l, d));
    Vall *= Vd;

    /* update the vector. */
    vector_set(V, d, vector_get(V, d) * Vd);
  }

  /* return the product of all factors. */
  return Vall;
}

/* * * * multivariate interactions: * * * */

/* interact_freq(): compute the interaction of two multidimensional
 * frequency modulation terms, with an option to exclude a single
 * dimension.
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j1: first signal basis index.
 *  @l1: first signal term index.
 *  @j2: second signal basis index.
 *  @l2: second signal term index.
 *  @k: signal phase index.
 *  @delta: excluded dimension index.
 *  @t: multidimensional time vector.
 *
 * returns:
 *  E[V_{k,l1}(t | omega_j1 \ omega_j1,d)
 *    V_{k,l2}(t | omega_j2 \ omega_j2,d)]
 */
double interact_freq (const model_t *mdl,
                      const unsigned int j1, const unsigned int l1,
                      const unsigned int j2, const unsigned int l2,
                      const unsigned int k, const unsigned int delta,
                      const vector_t *t) {
  /* @g: output interaction value. */
  double g = 1.0;

  /* for different signal indices, the interaction factors into
   * the product of the two independent expectations.
   */
  if (j1 != j2)
    return expect_freq(mdl, j1, k, l1, delta, t) *
           expect_freq(mdl, j2, k, l2, delta, t);

  /* multiply together the interactions along each dimension. */
  for (unsigned int d = 0; d < mdl->D; d++) {
    /* skip the excluded dimension. */
    if (d == delta) continue;

    /* extract the relevant quantities. */
    const double td = vector_get(t, d);

    /* update the interaction. */
    g *= nrm_cospair(mdl, j1, d, td, td,
                     VTAB(k, l1, d),
                     VTAB(k, l2, d));
  }

  /* return the computed interaction. */
  return g;
}

/* interact_all_freq(): compute a full vector of interactions of
 * two multidimensional frequency modulation terms.
 *
 * operation:
 *  V(d) <- V(d) * < cos(omega_j1,d * t_d + pi/2 * v_{k,l1,d})
 *                   cos(omega_j2,d * t_d + pi/2 * v_{k,l2,d}) >
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j1: first signal basis index.
 *  @l1: first signal term index.
 *  @j2: second signal basis index.
 *  @l2: second signal term index.
 *  @k: signal phase index.
 *  @t: multidimensional time vector.
 *  @V: output factor vector.
 *
 * returns:
 *  < V_{k,l1}(t | omega_j1) V_{k,l2}(t | omega_j2) >
 */
double interact_all_freq (const model_t *mdl,
                          const unsigned int j1, const unsigned int l1,
                          const unsigned int j2, const unsigned int l2,
                          const unsigned int k, const vector_t *t,
                          vector_t *V) {
  /* declare required variables:
   *  @Vd: storage for each factor.
   *  @Vall: product of all factors.
   */
  double td, Vd, Vall = 1.0;

  /* for different signal indices, the interaction factors into
   * the product of two independent expectations.
   */
  if (j1 != j2)
    return expect_all_freq(mdl, j1, l1, k, t, V) *
           expect_all_freq(mdl, j2, l2, k, t, V);

  /* compute the factor for each dimension. */
  for (unsigned int d = 0; d < mdl->D; d++) {
    /* update the product. */
    td = vector_get(t, d);
    Vd = nrm_cospair(mdl, j1, d, td, td, VTAB(k, l1, d), VTAB(k, l2, d));
    Vall *= Vd;

    /* update the vector. */
    vector_set(V, d, vector_get(V, d) * Vd);
  }

  /* return the product of all factors. */
  return Vall;
}

/* * * * multivariate covariances: * * * */

/* cov_freq(): compute the covariance of two modulation terms with
 * shared variables omega(j).
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: signal basis index.
 *  @k1: first phase index.
 *  @k2: second phase index.
 *  @l1: first term index.
 *  @l2: second term index.
 *  @t1: first time.
 *  @t2: second time.
 *
 * returns:
 *  E[V_{k1,l1}(t1 | omega_j) V_{k2,l2}(t2 | omega_j)]
 */
double cov_freq (const model_t *mdl, const unsigned int j,
                 const unsigned int k1, const unsigned int k2,
                 const unsigned int l1, const unsigned int l2,
                 const vector_t *t1, const vector_t *t2) {
  /* @cov: output covariance value. */
  double cov = 1.0;

  /* multiply together the covariances along each dimension. */
  for (unsigned int d = 0; d < mdl->D; d++) {
    /* extract the two times along the current dimension. */
    const double td1 = vector_get(t1, d);
    const double td2 = vector_get(t2, d);

    /* update the covariance. */
    cov *= nrm_cospair(mdl, j, d, td1, td2,
                       VTAB(k1, l1, d),
                       VTAB(k2, l2, d));
  }

  /* return the computed covariance. */
  return cov;
}

/* * * * univariate partial derivatives: * * * */

/* diff_nrm_cos(): compute the partial derivatives of the expectation
 * of a phased, stretched cosine under a normal distribution.
 *
 * i.e.:
 *  d mu  E[cos(x*t+z*pi/2)]_N(x | mu_jd, tau_jd)
 *  d tau E[cos(x*t+z*pi/2)]_N(x | mu_jd, tau_jd)
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: signal basis index.
 *  @d: signal dimension.
 *  @t: time shift.
 *  @z: integer phase offset.
 *  @Dmu: pointer to the output mu-derivative.
 *  @Dtau: pointer to the output tau-derivative.
 */
void diff_nrm_cos (const model_t *mdl,
                   const unsigned int j,
                   const unsigned int d,
                   const double t, const int z,
                   double *Dmu, double *Dtau) {
  /* compute intermediate quantities. */
  const double theta = mdl->sig[j][d].mu * t + M_PI_2 * (double) z;
  const double tau = mdl->sig[j][d].tau;

  /* ... and some intermediate function values. */
  const double E = exp(-0.5 * t * t / tau);
  const double C = cos(theta);
  const double S = sin(theta);

  /* compute and return the partial derivatives. */
  *Dmu = -t * E * S;
  *Dtau = 0.5 * ((t * t) / (tau * tau)) * E * C;
}

/* diff_nrm_cospair(): compute the partial derivatives of the expectation
 * of a product of differently phased and stretched cosine functions
 * under a normal distribution.
 *
 * i.e.:
 *  d mu  E[cos(x*t1+z1*pi/2)*cos(x*t2+z2*pi/2)]_N(x | mu_jd, tau_jd)
 *  d tau E[cos(x*t1+z1*pi/2)*cos(x*t2+z2*pi/2)]_N(x | mu_jd, tau_jd)
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: signal basis index.
 *  @d: signal dimension.
 *  @t1, @t2: time shifts.
 *  @z1, @z2: integer phase offsets.
 *  @Dmu: pointer to the output mu-derivative.
 *  @Dtau: pointer to the output tau-derivative.
 */
void diff_nrm_cospair (const model_t *mdl,
                       const unsigned int j,
                       const unsigned int d,
                       const double t1, const double t2,
                       const int z1, const int z2,
                       double *Dmu, double *Dtau) {
  /* compute the sum and difference of the phase offsets. */
  const int zp = z1 + z2;
  const int zm = z1 - z2;

  /* compute the sum and difference of the times. */
  const double tp = t1 + t2;
  const double tm = t1 - t2;

  /* extract relevant parameters and square the eta parameter. */
  const double mu = mdl->sig[j][d].mu;
  const double tau = mdl->sig[j][d].tau;

  /* compute the sum and difference of the trig arguments. */
  const double Ap = mu * tp + M_PI_2 * (double) zp;
  const double Am = mu * tm + M_PI_2 * (double) zm;

  /* compute the sum and difference exponential terms. */
  const double ep = tp * exp(-0.5 * tp * tp / tau);
  const double em = tm * exp(-0.5 * tm * tm / tau);

  /* compute the sum and difference sine/cosine function values. */
  const double Cp = ep * cos(Ap);
  const double Cm = em * cos(Am);
  const double Sp = ep * sin(Ap);
  const double Sm = em * sin(Am);

  /* compute and return the partial derivatives. */
  *Dmu = -0.5 * (Sp + Sm);
  *Dtau = 0.25 * (tp * Cp + tm * Cm) / (tau * tau);
}

/* * * * multivariate partial derivatives: * * * */

/* diff_expect_freq(): compute all partial derivatives of the expectation
 * of a multidimensional frequency modulation term.
 *
 * operation:
 *  Dm <- Dm + scale * d_mu  < V_{k,l}(t | omega_j) >
 *  Dt <- Dt + scale * d_tau < V_{k,l}(t | omega_j) >
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: signal basis index.
 *  @k: signal phase index.
 *  @l: signal term index.
 *  @t: multidimensional time vector.
 *  @scale: summation scale factor.
 *  @V: intermediate storage vector.
 *  @Dm: output mu-derivative vector.
 *  @Dt: output tau-derivative vector.
 */
void diff_expect_freq (const model_t *mdl, const unsigned int j,
                       const unsigned int k, const unsigned int l,
                       const vector_t *t, const double scale,
                       const vector_t *V,
                       vector_t *Dm,
                       vector_t *Dt) {
  /* gain access to the model dimension count. */
  const unsigned int D = mdl->D;

  /* compute the final derivatives. */
  for (unsigned int d = 0; d < D; d++) {
    /* extract the relevant quantities. */
    const double td = vector_get(t, d);

    /* compute the partial derivative factor. */
    double dmu, dtau;
    diff_nrm_cos(mdl, j, d, td, VTAB(k, l, d), &dmu, &dtau);

    /* multiply in the other expectation factors. */
    for (unsigned int delta = 0; delta < D; delta++) {
      if (delta == d) continue;
      dmu *= vector_get(V, delta);
      dtau *= vector_get(V, delta);
    }

    /* store the final derivatives. */
    vector_set(Dm, d, vector_get(Dm, d) + scale * dmu);
    vector_set(Dt, d, vector_get(Dt, d) + scale * dtau);
  }
}

/* diff_interact_freq(): compute the partial derivatives of the
 * interaction of two multidimensional frequency modulation terms.
 *
 * operation:
 *  Dm <- Dm + scale * d_mu  < V_{k,l1}(t | omega_j1) V_{k,l2}(t | omega_j2) >
 *  Dt <- Dt + scale * d_tau < V_{k,l1}(t | omega_j1) V_{k,l2}(t | omega_j2) >
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j1: first signal basis index.
 *  @l1: first signal term index.
 *  @j2: second signal basis index.
 *  @l2: second signal term index.
 *  @k: signal phase index.
 *  @t: multidimensional time vector.
 *  @scale: summation scale factor.
 *  @V: intermediate storage vector.
 *  @Dm: output mu-derivative vector.
 *  @Dt: output tau-derivative vector.
 */
void diff_interact_freq (const model_t *mdl,
                         const unsigned int j1, const unsigned int l1,
                         const unsigned int j2, const unsigned int l2,
                         const unsigned int k, const vector_t *t,
                         const double scale, vector_t *V,
                         vector_t *Dm, vector_t *Dt) {
  /* gain access to the model dimension count. */
  const unsigned int D = mdl->D;

  /* for different signal indices, the interaction factors into
   * the product of the two independent expectations.
   */
  if (j1 != j2) {
    /* compute the partial derivatives of the first factor,
     * appropriately scaled to include the second expectation.
     */
    const double scale2 = scale * expect_freq(mdl, j2, k, l2, D, t);
    diff_expect_freq(mdl, j1, k, l1, t, scale2, V, Dm, Dt);
  }

  /* compute the final derivatives. */
  for (unsigned int d = 0; d < D; d++) {
    /* extract the relevant quantities. */
    const double td = vector_get(t, d);

    /* compute the partial derivative factor. */
    double dmu, dtau;
    diff_nrm_cospair(mdl, j1, d, td, td,
                     VTAB(k, l1, d),
                     VTAB(k, l2, d),
                     &dmu, &dtau);

    /* multiply in the other expectation factors. */
    for (unsigned int delta = 0; delta < D; delta++) {
      if (delta == d) continue;
      dmu *= vector_get(V, delta);
      dtau *= vector_get(V, delta);
    }

    /* store the final derivatives. */
    vector_set(Dm, d, vector_get(Dm, d) + scale * dmu);
    vector_set(Dt, d, vector_get(Dt, d) + scale * dtau);
  }
}

/* diff_naturalize_freq(): convert a vector of gradients of a single
 * model signal frequency component into their natural forms.
 *
 * arguments:
 *  @mdl: model structure pointer.
 *  @j: first signal basis index.
 *  @Dm: input/output mu-derivative vector.
 *  @Dt: input/output tau-derivative vector.
 *
 * returns:
 *  smallest eigenvalue of the fisher matrix
 *  used to naturalize the gradient.
 */
double diff_naturalize_freq (const model_t *mdl, const unsigned int j,
                             vector_t *Dm, vector_t *Dt) {
  /* initialize the minimum eigenvalue. */
  double evmin = 1.0e+9;

  /* loop over each dimension. */
  for (unsigned int d = 0; d < mdl->D; d++) {
    /* get the values of the gradient. */
    const double g1 = vector_get(Dm, d);
    const double g2 = vector_get(Dt, d);

    /* get the fisher matrix elements. */
    const double F11 = mdl->sig[j][d].Fomega[0];
    const double F22 = mdl->sig[j][d].Fomega[1];

    /* compute the natural gradient values. */
    const double n1 = g1 / F11;
    const double n2 = g2 / F22;

    /* store the naturalized values. */
    vector_set(Dm, d, n1);
    vector_set(Dt, d, n2);

    /* update the minimum eigenvalue. */
    if (F11 < evmin) evmin = F11;
    if (F22 < evmin) evmin = F22;
  }

  /* return the minimum eigenvalue. */
  return evmin;
}

