
/* get_prior(): copy the nonlinear priors of a single basis into
 * a set of vectors.
 *
 * arguments:
 *  @opt: pointer to the optimization structure to utilize.
 *  @j: model basis index to update parameters for.
 *  @mu, @tau, @alpha, @beta: output vectors.
 */
static void
get_prior (const optim_t *opt, const unsigned int j,
           vector_t *mu, vector_t *tau,
           vector_t *alpha, vector_t *beta) {
  /* loop over each dimension. */
  for (unsigned int d = 0; d < opt->mdl->D; d++) {
    /* copy the prior parameters into the vectors. */
    vector_set(mu,    d, opt->mdl->sig[j][d].mu0);
    vector_set(tau,   d, opt->mdl->sig[j][d].tau0);
    vector_set(alpha, d, opt->mdl->sig[j][d].alpha0);
    vector_set(beta,  d, opt->mdl->sig[j][d].beta0);
  }
}

/* get_parms(): copy the nonlinear parameters of a single basis into
 * a set of vectors.
 *
 * arguments:
 *  @opt: pointer to the optimization structure to utilize.
 *  @j: model basis index to update parameters for.
 *  @mu, @tau, @alpha, @beta: output vectors.
 */
static void
get_parms (const optim_t *opt, const unsigned int j,
           vector_t *mu, vector_t *tau,
           vector_t *alpha, vector_t *beta) {
  /* loop over each dimension. */
  for (unsigned int d = 0; d < opt->mdl->D; d++) {
    /* copy the current parameters into the vectors. */
    vector_set(mu,    d, opt->mdl->sig[j][d].mu);
    vector_set(tau,   d, opt->mdl->sig[j][d].tau);
    vector_set(alpha, d, opt->mdl->sig[j][d].alpha);
    vector_set(beta,  d, opt->mdl->sig[j][d].beta);
  }
}

/* set_parms(): copy the nonlinear parameters of a single basis back
 * into a gp-associated model from a set of vectors.
 *
 * arguments:
 *  @opt: pointer to the optimization structure to utilize.
 *  @j: model basis index to update parameters for.
 *  @mu, @tau, @alpha, @beta: output vectors.
 */
static void
set_parms (optim_t *opt, const unsigned int j,
           vector_t *mu, vector_t *tau,
           vector_t *alpha, vector_t *beta) {
  /* loop over each dimension. */
  for (unsigned int d = 0; d < opt->mdl->D; d++) {
    /* copy the current parameters out from the vectors. */
    model_set_mu(opt->mdl,    j, d, vector_get(mu,    d));
    model_set_tau(opt->mdl,   j, d, vector_get(tau,   d));
    model_set_alpha(opt->mdl, j, d, vector_get(alpha, d));
    model_set_beta(opt->mdl,  j, d, vector_get(beta,  d));
  }
}

/* natgrad(): compute the natural gradient of the non-convex portion
 * of the collapsed variational lower bound on the log-evidence of a
 * gp-associated signal model, with respect to the nonlinear
 * parameters of a single model signal.
 *
 * arguments:
 *  @opt: pointer to the optimization structure to access.
 *  @j: model signal index to update parameters for.
 *  @R: temporary vector for expected decay factors.
 *  @V: temporary vector for expected modulation factors.
 *  @mu: output vector of frequency-mean gradients.
 *  @tau: output vector of frequency-precision gradients.
 *  @alpha: output vector of decay-shape gradients.
 *  @beta: output vector of decay-rate gradients.
 *
 * returns:
 *  smallest eigenvalue of the fisher information matrix.
 */
static double
natgrad (const optim_t *opt, const unsigned int j,
         vector_t *R, vector_t *V,
         vector_t *mu, vector_t *tau,
         vector_t *alpha, vector_t *beta) {
  /* declare a variable to hold measurements. */
  double y;

  /* gain access to required structure members. */
  model_t *mdl = opt->mdl;
  dataset_t *dat = opt->dat;
  const unsigned int D = mdl->D;
  const unsigned int K = mdl->K;
  const unsigned int M = mdl->M;
  const unsigned int *n = dat->n;

  /* compute the noise precision. */
  const double invs2 = 1.0 / (mdl->sigma * mdl->sigma);

  /* declare vector views for measurement times and covariances. */
  vector_view_t t, s;

  /* create a vector view for storing measurement times. */
  t = vector_subvector(opt->tmp, 0, D);

  /* loop over the basis terms. */
  for (unsigned int l = 0; l < K; l++) {
    /* gain a view into the amplitude covariance matrix. */
    s = matrix_row(mdl->Sigma, l * M + j);

    /* compute the first amplitude moment. */
    const double A1 = vector_get(mdl->ahat, l * M + j);

    /* loop over the signal phases. */
    for (unsigned int k = 0; k < K; k++) {
      /* loop over the measurements of the current phase. */
      for (unsigned int i = 0; i < n[k]; i++) {
        /* extract the current measurement. */
        dataset_get(dat, k, i, &t, &y);

        /* compute the current vectors of factors. */
        vector_set_all(R, 1.0);
        vector_set_all(V, 1.0);
        const double Rall = expect_all_decay(mdl, j, &t, R);
        const double Vall = expect_all_freq(mdl, j, k, l, &t, V);

        /* compute the scale factors of the projection contribution. */
        const double sh = invs2 * A1 * y;
        const double sr = sh * Vall;
        const double sv = sh * Rall;

        /* include the projection contribution. */
        diff_expect_decay(mdl, j, &t, sr, R, alpha, beta);
        diff_expect_freq(mdl, j, k, l, &t, sv, V, mu, tau);

        /* loop over all other signals to compute
         * the interaction contributions.
         */
        for (unsigned int j2 = 0; j2 < M; j2++) {
          /* compute the current vector of decay factors. */
          vector_set_all(R, 1.0);
          const double Rall = interact_all_decay(mdl, j, j2, &t, R);

          /* loop again over the basis terms. */
          for (unsigned int l2 = 0; l2 < K; l2++) {
            /* compute the second amplitude moment. */
            const double A2 = vector_get(&s, l2 * M + j2) +
                              A1 * vector_get(mdl->ahat, l2 * M + j2);

            /* compute the current vector of modulation factors. */
            vector_set_all(V, 1.0);
            const double Vall = interact_all_freq(mdl, j, l, j2, l2,
                                                  k, &t, V);

            /* compute the scale factors of the interaction contribution. */
            const double sh = -(j == j2 ? 0.5 : 1.0) * invs2 * A2;
            const double sr = sh * Vall;
            const double sv = sh * Rall;

            /* include the interaction contribution. */
            diff_interact_decay(mdl, j, j2, &t, sr, R, alpha, beta);
            diff_interact_freq(mdl, j, l, j2, l2, k, &t, sv, V, mu, tau);
          }
        }
      }
    }
  }

  /* naturalize the gradients. */
  const double ev1 = diff_naturalize_freq(mdl, j, mu, tau);
  const double ev2 = diff_naturalize_decay(mdl, j, alpha, beta);

  /* return the minimum eigenvalue. */
  return (ev1 < ev2 ? ev1 : ev2);
}

/* update_nonlinear(): update the nonlinear parameters of a signal model
 * based on a dataset and the current linear model parameters.
 *
 * arguments:
 *  @opt: pointer to the optimization structure to update.
 *  @j: model basis index to update parameters for.
 */
static void
update_nonlinear (optim_t *opt, const unsigned int j) {
  /* declare line search variables:
   *  @valid: whether or not the current step is valid.
   *  @steps: number of line search steps taken.
   *  @elbo: current lower bound value.
   *  @gamma: proximal step length.
   */
  unsigned int valid, steps;
  double elbo, gamma;

  /* gain access to required structure members. */
  model_t *mdl = opt->mdl;
  const unsigned int D = mdl->D;

  /* declare views for holding measurement times and parameter vectors. */
  vector_view_t a, b, x, R, V, mu, tau, alpha, beta;

  /* create vector views for grouped parameters. */
  a = vector_subvector(opt->tmp,  1 * D, 4 * D);
  b = vector_subvector(opt->tmp,  5 * D, 4 * D);
  x = vector_subvector(opt->tmp,  9 * D, 4 * D);
  R = vector_subvector(opt->tmp, 13 * D, D);
  V = vector_subvector(opt->tmp, 14 * D, D);

  /* create vector views for factors and individual parameters. */
  mu    = vector_subvector(&x, 0, D);
  tau   = vector_subvector(&x, D, D);
  alpha = vector_subvector(&x, 2 * D, D);
  beta  = vector_subvector(&x, 3 * D, D);

  /* copy the current parameters. */
  get_parms(opt, j, &mu, &tau, &alpha, &beta);
  vector_copy(&a, &x);

  /* copy the prior parameters. */
  get_prior(opt, j, &mu, &tau, &alpha, &beta);
  vector_copy(&b, &x);

  /* compute the natural gradient and sum them with the prior. */
  vector_set_zero(&x);
  gamma = natgrad(opt, j, &R, &V, &mu, &tau, &alpha, &beta);
  vector_add(&b, &x);

  /* initialize the line search variables. */
  gamma /= opt->l0;
  valid = 0;
  steps = 0;
  do {
    /* compute the unconstrained proximal gradient step. */
    vector_set_zero(&x);
    blas_daxpy(1.0 / (gamma + 1.0), &a, &x);
    blas_daxpy(gamma / (gamma + 1.0), &b, &x);

    /* set the newly computed parameters. */
    set_parms(opt, j, &mu, &tau, &alpha, &beta);

    /* update the linear parameters and recompute the bound. */
    update_linear_single(opt, j);
    elbo = evaluate_bound(opt);

    /* update the step length and count for subsequent steps. */
    gamma *= opt->dl;
    steps++;

    /* check if the current step is valid. */
    if (elbo >= opt->elbo &&
        vector_positive(&tau) &&
        vector_positive(&alpha) &&
        vector_positive(&beta))
      valid = 1;
  }
  while (!valid && steps < opt->step_max);

  /* check if a valid step was performed. */
  if (valid) {
    /* yes. store the new lower bound. */
    opt->elbo = elbo;
  }
  else {
    /* no. reset the parameters and recompute the linear parameters. */
    vector_copy(&x, &a);
    set_parms(opt, j, &mu, &tau, &alpha, &beta);
    update_linear_single(opt, j);
  }
}

