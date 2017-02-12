
/* reduce_cosine(): compute the amplitude and phase shift of a cosine
 * formed by the linear combination of scaled and shifted cosines
 * of equal frequency.
 *
 * arguments:
 *  @a: vector of input amplitudes.
 *  @theta: vector of input phase shifts.
 *  @A: pointer to the output amplitude.
 *  @Theta: pointer to the output phase shift.
 */
static inline void
reduce_cosine (const vector_t *a,
                const vector_t *theta,
                double *A, double *Theta) {
  /* declare temporary scalars. */
  double S, C, Asq;

  /* get the vector lengths. */
  const unsigned int n = a->len;

  /* compute the sums required for shift computation. */
  S = C = 0.0;
  for (unsigned int i = 0; i < n; i++) {
    S += vector_get(a, i) * sin(vector_get(theta, i));
    C += vector_get(a, i) * cos(vector_get(theta, i));
  }

  /* compute the sum required for amplitude computation. */
  Asq = 0.0;
  for (unsigned int i = 0; i < n; i++) {
    /* get the first pair of values. */
    const double ai = vector_get(a, i);
    const double ti = vector_get(theta, i);

    /* include the diagonal term. */
    Asq += ai * ai;

    /* loop over the off-diagonal values. */
    for (unsigned int j = i + 1; j < n; j++) {
      /* get the second pair of values. */
      const double aj = vector_get(a, j);
      const double tj = vector_get(theta, j);

      /* include the off-diagonal term. */
      Asq += 2.0 * ai * aj * cos(ti - tj);
    }
  }

  /* store the output values. */
  *A = sqrt(Asq);
  *Theta = atan2(S, C);
}

/* eval_factor(): compute the natural logarithm of a variational
 * update/factor, given its parameters (xi, gamma, nu).
 *
 * arguments:
 *  @x: frequency at which to evaluate the log-probability.
 *  @mu0: prior mean over the frequency.
 *  @tau0: prior precision over the frequency.
 *  @xi: vector of cosine amplitudes.
 *  @gamma: vector of cosine scalings.
 *  @nu: vector of cosine phase shifts.
 *
 * returns:
 *  ln q*(x | omega)*p(omega) for (xi, gamma, nu).
 */
static double
eval_factor (const double x, const double mu0, const double tau0,
             const vector_t *xi,
             const vector_t *gamma,
             const vector_t *nu) {
  /* initialize the computation to the log-prior. */
  double Q = x - mu0;
  Q = -0.5 * tau0 * Q * Q;

  /* loop over each cosine term. */
  for (unsigned int i = 0; i < xi->len; i++) {
    /* extract the current term parameters. */
    const double xi_i    = vector_get(xi, i);
    const double gamma_i = vector_get(gamma, i);
    const double nu_i    = vector_get(nu, i);

    /* include the term into the computation. */
    Q += xi_i * cos(gamma_i * x + nu_i);
  }

  /* return the computed result. */
  return Q;
}

/* laplace(): perform laplace approximation of a variational factor
 * by sampling at *most* possible modes of the log-probability.
 *
 * arguments:
 *  @L, @U: lower and upper bounds induced by the prior.
 *  @mu0, @tau0: prior mean and precision.
 *  @xi, @gamma, @nu: factor parameters.
 *  @mu, @tau: pointers to the output mean and precision.
 */
static inline void
laplace (const double L, const double U,
         const double mu0, const double tau0,
         const vector_t *xi,
         const vector_t *gamma,
         const vector_t *nu,
         double *mu, double *tau) {
  /* initialize the search location and amplitude. */
  double xmax = mu0;
  double Qmax = eval_factor(xmax, mu0, tau0, xi, gamma, nu);

  /* loop over each cosine term in the log-probability. */
  for (unsigned int i = 0; i < xi->len; i++) {
    /* get the current term parameters. */
    const double gamma_i = vector_get(gamma, i);
    const double nu_i    = vector_get(nu, i);

    /* compute the center and period of the modes of the term. */
    const double m = nu_i / gamma_i;
    const double T = (2.0 * M_PI) / gamma_i;

    /* determine the start and end indices of the modes. */
    const int n1 = (int) ceil((L - m) / T);
    const int n2 = (int) floor((U - m) / T);

    /* loop over the modes. */
    for (int n = n1; n <= n2; n++) {
      /* compute the mode location and amplitude. */
      const double x = m + ((double) n) * T;
      const double Q = eval_factor(x, mu0, tau0, xi, gamma, nu);

      /* if the mode amplitude is greater, store it. */
      if (Q > Qmax) {
        xmax = x;
        Qmax = Q;
      }
    }
  }

  /* compute the precision at the identified mode. */
  Qmax = 0.0;
  for (unsigned int i = 0; i < xi->len; i++) {
    /* get the current term parameters. */
    const double xi_i    = vector_get(xi, i);
    const double gamma_i = vector_get(gamma, i);
    const double nu_i    = vector_get(nu, i);

    /* include the current term in the calculation. */
    Qmax += xi_i * gamma_i * gamma_i * cos(gamma_i * xmax + nu_i);
  }

  /* return the approximate posterior parameters. */
  *mu = xmax;
  *tau = Qmax;
}

/* update_meanfield(): use an approximate mean-field update to compute
 * new parameters for a single signal of a model, along a single dimension.
 *
 * arguments:
 *  @opt: optimizer structure pointer to access.
 *  @j: signal index of the update.
 *  @d: dimension index of the update.
 */
static void
update_meanfield (optim_t *opt,
                  const unsigned int j,
                  const unsigned int d) {
  /* declare a variable to hold measurements. */
  double y;

  /* gain access to required structure members. */
  model_t *mdl = opt->mdl;
  dataset_t *dat = opt->dat;
  const unsigned int D = mdl->D;
  const unsigned int K = mdl->K;
  const unsigned int M = mdl->M;
  const unsigned int N = dat->N;
  const unsigned int *n = dat->n;

  /* compute the noise precision. */
  const double invs2 = 1.0 / (mdl->sigma * mdl->sigma);

  /* gain access to the prior parameters. */
  const double mu0 = mdl->sig[j][d].mu0;
  const double tau0 = mdl->sig[j][d].tau0;

  /* declare bounds for the laplace approximation. */
  const double lb = mu0 - 3.0 / sqrt(tau0);
  const double ub = mu0 + 3.0 / sqrt(tau0);

  /* declare vector views for intermediate quantities. */
  vector_view_t t, b, x, c, z, xi, gamma, nu;
  unsigned int off = 0;

  /* create vector views for the intermediate quantities. */
  t     = vector_subvector(opt->tmp, off, D);     off += D;
  b     = vector_subvector(opt->tmp, off, K);     off += K;
  x     = vector_subvector(opt->tmp, off, K);     off += K;
  c     = vector_subvector(opt->tmp, off, K * K); off += K * K;
  z     = vector_subvector(opt->tmp, off, K * K); off += K * K;
  xi    = vector_subvector(opt->tmp, off, 2 * N); off += 2 * N;
  gamma = vector_subvector(opt->tmp, off, 2 * N); off += 2 * N;
  nu    = vector_subvector(opt->tmp, off, 2 * N); off += 2 * N;

  /* compute the parameter vectors of the exact variational update. */
  for (unsigned int k = 0, I = 0; k < K; k++) {
    for (unsigned int i = 0; i < n[k]; i++, I++) {
      /* extract the current measurement. */
      dataset_get(dat, k, i, &t, &y);
      const double td = vector_get(&t, d);

      /* compute the first and second decay moments. */
      const double R1 = expect_decay(mdl, j, D, &t);
      const double R2 = interact_decay(mdl, j, j, D, &t);

      /* compute portions of the factor parameters. */
      const double b0 = invs2 * R1;
      const double c0 = -0.25 * invs2 * R2;

      /* loop over the signal phases. */
      for (unsigned int l1 = 0; l1 < K; l1++) {
        /* compute the first amplitude and frequency moments. */
        const double al = vector_get(mdl->ahat, l1 * M + j);
        const double Vl = expect_freq(mdl, j, k, l1, d, &t);
        double b2 = 0.0;

        /* loop over the signal terms. */
        for (unsigned int l2 = 0; l2 < K; l2++) {
          /* compute the second amplitude and frequency moments. */
          const double all = matrix_get(mdl->Sigma, l1 * M + j, l2 * M + j)
                           + al * vector_get(mdl->ahat, l2 * M + j);
          const double Vll = interact_freq(mdl, j, l1, j, l2, k, d, &t);

          /* loop over all other signals. */
          for (unsigned int j2 = 0; j2 < M; j2++) {
            if (j2 == j) continue;

            /* compute the other amplitude, decay and frequency moments. */
            const double a = matrix_get(mdl->Sigma, l1 * M + j, l2 * M + j2)
                           + vector_get(mdl->ahat, l1 * M + j) *
                             vector_get(mdl->ahat, l2 * M + j2);
            const double Rj2 = expect_decay(mdl, j2, D, &t);
            const double Vj2 = expect_freq(mdl, j2, k, l2, D, &t);

            /* include the interaction with the other signal. */
            b2 += a * Rj2 * Vj2;
          }

          /* compute the second-order results. */
          const double cll = c0 * all * Vll;
          const double zll = M_PI_2 * (VTAB(k, l1, d) + VTAB(k, l2, d));

          /* store the second-order results. */
          vector_set(&c, l1 + l2 * K, cll);
          vector_set(&z, l1 + l2 * K, zll);
        }

        /* compute the first-order results. */
        const double bl = b0 * Vl * (y * al - b2);
        const double xl = M_PI_2 * VTAB(k, l1, d);

        /* store the first-order results. */
        vector_set(&b, l1, bl);
        vector_set(&x, l1, xl);
      }

      /* compute the reduced form of the current first-order term. */
      double bI, xI;
      reduce_cosine(&b, &x, &bI, &xI);

      /* compute the reduced form of the current second-order term. */
      double cI, zI;
      reduce_cosine(&c, &z, &cI, &zI);

      /* store the reduced scaling coefficients. */
      vector_set(&gamma, I, td);
      vector_set(&gamma, I + N, 2.0 * td);

      /* store the reduced amplitude coefficients. */
      vector_set(&xi, I, bI);
      vector_set(&xi, I + N, cI);

      /* store the reduced phase shifts. */
      vector_set(&nu, I, xI);
      vector_set(&nu, I + N, zI);
    }
  }

  /* perform laplace approximation to determine
   * the new factor parameters.
   */
  double mu_new, tau_new;
  laplace(lb, ub, mu0, tau0, &xi, &gamma, &nu, &mu_new, &tau_new);

  /* set the new factor parameters. */
  model_set_mu(mdl, j, d, mu_new);
  model_set_tau(mdl, j, d, tau_new);
}

