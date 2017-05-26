
/* include the quadrature and vfl object headers. */
#include <vbnmr/quad.h>
#include <vfl/base/list.h>
#include <vfl/base/float.h>

/* include the fftw header. */
#include <fftw3.h>

/* define the parameter indices. */
#define P_MU   0
#define P_TAU  1

/* VTAB(): macro for accessing quadrature factor modulation table values.
 *
 * note: the arguments @k, @l, and @d do not totally correspond with the
 * VFL standard, as they have origins in the specific NMR signal model,
 * so be careful when using this function!
 *
 * arguments:
 *  @k: signal function output index.
 *  @l: signal function weight index.
 *  @d: signal function dimension.
 *
 * returns:
 *  integer value of the modulation table on the specified output,
 *  term/weight and dimension.
 */
#define VTAB(k,l,d) ((quadrature_t*) (f))->vtab[k][l][d]

/* universally accessed mean-field update variables:
 *  @Dft: maximum dimension count of fft arrays.
 *  @Nft: maximum number of complex fft coefficients.
 *  @nq: utilized number of columns in the @q matrix.
 *  @q: matrices holding (gamma,xi,nu,gamma,xi,nu,...).
 *  @Q: matrices holding inverse fast fourier transforms.
 *  @plan: plan for fast fourier transforms.
 */
static unsigned int Dft = 0;
static unsigned int Nft = 0;
static unsigned int nq = 0;
static matrix_t *q, *Q;
static fftw_plan plan;

/* quadrature_t: structure for holding a quadrature factor.
 */
typedef struct {
  /* @base: core factor structure members. */
  factor_t base;

  /* dimensionality-dependent members:
   *  @vtab: quadrature modulation table.
   */
  int ***vtab;

  /* mean-field related members:
   *  @a, @x: amplitude and phase vectors.
   *  @C, @Z: amplitude and phase matrices.
   */
  vector_t *a, *x;
  matrix_t *C, *Z;
}
quadrature_t;

/* tab1_alloc(): allocate a one-dimensional table of integers.
 *
 * arguments:
 *  @n1: first-dimension size of the table.
 *
 * returns:
 *  pointer to a new one-dimensional table.
 */
static int *tab1_alloc (const unsigned int n1) {
  /* allocate the table. */
  unsigned int bytes = n1 * sizeof(int);
  int *tab = (int*) malloc(bytes);
  if (!tab)
    return NULL;

  /* initialize the table contents. */
  for (unsigned int i1 = 0; i1 < n1; i1++)
    tab[i1] = 0;

  /* return the table. */
  return tab;
}

/* tab2_alloc(): allocate a two-dimensional table of integers.
 *
 * arguments:
 *  @n1: first-dimension size of the table.
 *  @n2: second-dimension size of the table.
 *
 * returns:
 *  pointer to a new two-dimensional table.
 */
static int **tab2_alloc (const unsigned int n1,
                         const unsigned int n2) {
  /* allocate the table. */
  unsigned int bytes = n1 * n2 * sizeof(int) + n1 * sizeof(int*);
  int **tab = (int**) malloc(bytes);
  if (!tab)
    return NULL;

  /* create a pointer for initializing sub-tables. */
  char *ptr = (char*) tab;
  ptr += n1 * sizeof(int*);

  /* initialize each sub-table. */
  for (unsigned int i1 = 0; i1 < n1; i1++) {
    tab[i1] = (int*) ptr;
    ptr += n2 * sizeof(int);
  }

  /* initialize the table contents. */
  for (unsigned int i1 = 0; i1 < n1; i1++)
    for (unsigned int i2 = 0; i2 < n2; i2++)
      tab[i1][i2] = 0;

  /* return the table. */
  return tab;
}

/* vtab_bytes(): compute the storage size required for
 * a phase table of a specified dimensionality.
 *
 * arguments:
 *  @D: dimensionality of the signal.
 *  @K: phase count of the signal.
 *
 * returns:
 *  number of bytes required to store such a phase table.
 */
static unsigned int vtab_bytes (const unsigned int D,
                                const unsigned int K) {
  /* phase tables are three-dimensional arrays of integers. */
  unsigned int bytes = K * K * D * sizeof(int);
  bytes += K * K * sizeof(int*);
  bytes += K * sizeof(int**);

  /* return the computed result. */
  return bytes;
}

/* vtab_init(): initialize the contents of a phase table.
 *
 * arguments:
 *  @tab: pointer to the phase table to initialize.
 *  @D: dimensionality of the signal.
 *  @K: phase count of the signal.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
static int vtab_init (int ***tab, const unsigned int D,
                      const unsigned int K) {
  /* create a pointer for initializing sub-tables. */
  char *ptr = (char*) tab;
  ptr += K * sizeof(int**);

  /* initialize the array of two-dimensional sub-tables. */
  for (unsigned int k = 0; k < K; k++) {
    tab[k] = (int**) ptr;
    ptr += K * sizeof(int*);
  }

  /* initialize the arrays of one-dimensional sub-tables. */
  for (unsigned int k = 0; k < K; k++) {
    for (unsigned int l = 0; l < K; l++) {
      tab[k][l] = (int*) ptr;
      ptr += D * sizeof(int);
    }
  }

  /* initialize the table contents. */
  for (unsigned int k = 0; k < K; k++)
    for (unsigned int l = 0; l < K; l++)
      for (unsigned int d = 0; d < D; d++)
        tab[k][l][d] = 0;

  /* allocate two temporary tables for computing the phases. */
  int **ph = tab2_alloc(K, D);
  int *nv = tab1_alloc(K);

  /* check that allocation was successful. */
  if (!ph || !nv)
    return 0;

  /* compute the euler expansion of the phase shifts that will
   * be implicitly contained in the model amplitudes.
   */
  for (unsigned int l = 0; l < K; l++)
    for (unsigned int d = 0; d < D; d++)
      ph[l][d] = l & (1 << d) ? 3 : 0;

  /* loop over all terms in the product between the phase shift
   * and modulation euler expansions.
   */
  for (unsigned int k1 = 0; k1 < K; k1++) {
    for (unsigned int k2 = 0; k2 < K; k2++) {
      /* determine the resulting phase of the current product of terms. */
      unsigned int k = (k1 | k2) & ~(k1 & k2);

      /* compute the final phase table value for each dimension. */
      for (unsigned int d = 0; d < D; d++) {
        int flip = (k1 & k2) & (1 << d) ? 2 : 0;
        tab[k][nv[k]][d] = (ph[k2][d] + flip) % 4;
      }

      /* increment the number of terms that belong to
       * the resultant phase.
       */
      nv[k]++;
    }
  }

  /* loop through each 'row' of the phase table. */
  for (unsigned int k = 0; k < K; k++) {
    for (unsigned int l = 0; l < K; l++) {
      /* loop over the dimensions to find negatives. */
      for (unsigned int d1 = 0; d1 < D; d1++) {
        if (tab[k][l][d1] != 1) continue;

        /* loop over the remaining dimensions to find negative pairs. */
        for (unsigned int d2 = d1 + 1; d2 < D; d2++) {
          if (tab[k][l][d2] != 1) continue;

          /* we found a pair of negative sines: cancel them. */
          tab[k][l][d1] = 3;
          tab[k][l][d2] = 3;
          break;
        }
      }
    }
  }

  /* free the temporary tables and return success. */
  free(ph);
  free(nv);
  return 1;
}

/* --- */

/* quad_cos(): compute the expectation of a quadrature factor
 * along a specified dimension, with a specified phase offset.
 *
 * arguments:
 *  @f: factor structure pointer.
 *  @d: dimension index.
 *  @x: input value.
 *  @z: phase offset.
 *
 * returns:
 *  < cos(omega[d] * x[d] + (pi/2)*z) >_N(omega | mu, tau)
 */
double quad_cos (const factor_t *f, const unsigned int d,
                 const double x, const int z) {
  /* get the parameters along the specified dimension. */
  const unsigned int p = 2 * d;
  const double mu = vector_get(f->par, p + P_MU);
  const double tau = vector_get(f->par, p + P_TAU);

  /* compute and return the expectation. */
  return exp(-0.5 * x * x / tau) * cos(mu * x + M_PI_2 * (double) z);
}

/* quad_cospair(): compute the expectation of a quadrature factor
 * along a specified dimension, at two specified input locations and
 * phase shifts.
 *
 * arguments:
 *  @f: factor structure pointer.
 *  @d: dimension index.
 *  @x1: first input value.
 *  @x2: second input value.
 *  @z1: first phase offset.
 *  @z2: second phase offset.
 *
 * returns:
 *  < cos(omega[d] * x1[d] + (pi/2)*z1)
 *    cos(omega[d] * x2[d] + (pi/2)*z2) >_N(omega | mu, tau)
 */
double quad_cospair (const factor_t *f, const unsigned int d,
                     const double x1, const double x2,
                     const int z1, const int z2) {
  /* compute the sum and difference of the phase offsets. */
  const int zp = z1 + z2;
  const int zm = z1 - z2;

  /* compute the sum and difference of the input values. */
  const double xp = x1 + x2;
  const double xm = x1 - x2;

  /* compute and return the expectation. */
  return 0.5 * (quad_cos(f, d, xp, zp) +
                quad_cos(f, d, xm, zm));
}

/* quad_diff_cos(): compute the partial derivatives of the
 * first moment of a quadrature factor along one dimension.
 *
 * arguments:
 *  @f: factor structure pointer.
 *  @d: dimension index.
 *  @x: input value.
 *  @z: phase offset.
 *  @dmu: pointer to the output mu-derivative.
 *  @dtau: pointer to the output tau-derivative.
 */
void quad_diff_cos (const factor_t *f, const unsigned int d,
                    const double x, const int z,
                    double *dmu, double *dtau) {
  /* get the parameters along the specified dimension. */
  const unsigned int p = 2 * d;
  const double mu = vector_get(f->par, p + P_MU);
  const double tau = vector_get(f->par, p + P_TAU);

  /* compute the trigonometric argument. */
  const double theta = mu * x + M_PI_2 * (double) z;

  /* compute some intermediate terms. */
  const double E = exp(-0.5 * x * x / tau);
  const double C = cos(theta);
  const double S = sin(theta);

  /* compute and return the partial derivatives. */
  *dmu = -x * E * S;
  *dtau = 0.5 * ((x * x) / (tau * tau)) * E * C;
}

/* quad_diff_cospair(): compute the partial derivatives of
 * the second moment of a quadrature factor along one dimension.
 *
 * arguments:
 *  @f: factor structure pointer.
 *  @d: dimension index.
 *  @x1: first input value.
 *  @x2: second input value.
 *  @z1: first phase offset.
 *  @z2: second phase offset.
 *  @dmu: pointer to the output mu-derivative.
 *  @dtau: pointer to the output tau-derivative.
 */
void quad_diff_cospair (const factor_t *f, const unsigned int d,
                        const double x1, const double x2,
                        const int z1, const int z2,
                        double *dmu, double *dtau) {
  /* get the parameters along the specified dimension. */
  const unsigned int p = 2 * d;
  const double mu = vector_get(f->par, p + P_MU);
  const double tau = vector_get(f->par, p + P_TAU);

  /* compute the sum and difference of the phase offsets. */
  const int zp = z1 + z2;
  const int zm = z1 - z2;

  /* compute the sum and difference of the input values. */
  const double xp = x1 + x2;
  const double xm = x1 - x2;

  /* compute the sum and difference trigonometric arguments. */
  const double Ap = mu * xp + M_PI_2 * (double) zp;
  const double Am = mu * xm + M_PI_2 * (double) zm;

  /* compute the sum and difference exponential terms. */
  const double ep = xp * exp(-0.5 * xp * xp / tau);
  const double em = xm * exp(-0.5 * xm * xm / tau);

  /* compute the sum and difference sine/cosine function values. */
  const double Cp = ep * cos(Ap);
  const double Cm = em * cos(Am);
  const double Sp = ep * sin(Ap);
  const double Sm = em * sin(Am);

  /* compute and return the partial derivatives. */
  *dmu = -0.5 * (Sp + Sm);
  *dtau = 0.25 * (xp * Cp + xm * Cm) / (tau * tau);
}

/* quad_cosreduce(): compute the amplitude and phase shift of a cosine
 * formed by the linear combination of scaled and shifted cosines
 * of equal frequency.
 *
 * arguments:
 *  @a: vector of input amplitudes.
 *  @theta: vector of input phases.
 *  @A: pointer to the output amplitude.
 *  @Theta: pointer to the output phase.
 */
void quad_cosreduce (const vector_t *a, const vector_t *theta,
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

  /* store the computed results. */
  *A = sqrt(Asq);
  *Theta = atan2(S, C);
}

/* quad_costune(): perform a local search of a one-dimensional
 * log-probability function using newton-rhapson iteration,
 * then estimate the location and curvature of the mode.
 *
 * arguments:
 *  @Gamma: vector of time-shifts.
 *  @Xi: vector of amplitudes.
 *  @Nu: vector of phase shifts.
 *  @omega0: initial location for local search.
 *  @mu0: prior location parameter.
 *  @tau0: prior precision parameter.
 *  @mu: pointer to the output location estimate.
 *  @tau: pointer to the output precision estimate.
 */
void quad_costune (const vector_t *Gamma, const vector_t *Xi,
                   const vector_t *Nu, const double omega0,
                   const double mu0, const double tau0,
                   double *mu, double *tau) {
  /* declare required variables:
   *  @f1, @f2: first and second function derivatives.
   *  @omega: current iteration location estimate.
   */
  double f1, f2, omega = omega0;

  /* loop a few times. the routine converges quite rapidly. */
  for (unsigned int it = 0; it < 3; it++) {
    /* initialize the derivates using the prior values. */
    f1 = tau0 * (omega - mu0);
    f2 = tau0;

    /* recompute the derivatives. */
    for (unsigned int i = 0; i < nq; i++) {
      /* get the current vector values. */
      const double gamma = vector_get(Gamma, i);
      const double xi    = vector_get(Xi, i);
      const double nu    = vector_get(Nu, i);

      /* add the current terms to the derivatives. */
      f1 += xi         * gamma * sin(gamma * omega + nu);
      f2 += xi * gamma * gamma * cos(gamma * omega + nu);
    }

    /* adjust the location estimate. */
    omega -= f1 / f2;
  }

  /* return the final results. */
  *mu = omega;
  *tau = f2;
}

/* --- */

/* quad_eval(): evaluate the quadrature factor at its mode.
 *  - see factor_mean_fn() for more information.
 */
FACTOR_EVAL (quad) {
  /* initialize the mode computation. */
  double fmode = 1.0;

  /* include the contributions along each dimension. */
  for (unsigned int d = 0, p0 = 0; d < f->D; d++, p0 += 2) {
    const double xd = vector_get(x, f->d + d);
    const double mu = vector_get(f->par, p0 + P_MU);
    fmode *= cos(mu * xd + M_PI_2 * (double) VTAB(p, i, d));
  }

  /* return the computed mode. */
  return fmode;
}

/* quad_mean(): evaluate the quadrature factor mean.
 *  - see factor_mean_fn() for more information.
 */
FACTOR_MEAN (quad) {
  /* initialize the mean computation. */
  double fmean = 1.0;

  /* include the contributions along each dimension. */
  for (unsigned int d = 0; d < f->D; d++) {
    const double xd = vector_get(x, f->d + d);
    fmean *= quad_cos(f, d, xd, VTAB(p, i, d));
  }

  /* return the computed mean. */
  return fmean;
}

/* quad_var(): evaluate the quadrature factor variance.
 *  - see factor_var_fn() for more information.
 */
FACTOR_VAR (quad) {
  /* initialize the variance computation. */
  double fvar = 1.0;

  /* include the contributions along each dimension. */
  for (unsigned int d = 0; d < f->D; d++) {
    const double xd = vector_get(x, f->d + d);
    fvar *= quad_cospair(f, d, xd, xd, VTAB(p, i, d), VTAB(p, j, d));
  }

  /* return the computed variance. */
  return fvar;
}

/* quad_cov(): evaluate the quadrature factor covariance.
 *  - see factor_cov_fn() for more information.
 */
FACTOR_COV (quad) {
  /* initialize the covariance computation. */
  double fcov = 1.0;

  /* include the contributions along each dimension. */
  for (unsigned int d = 0; d < f->D; d++) {
    /* get the current dimension input value difference. */
    const double xd = vector_get(x1, f->d + d) - vector_get(x2, f->d + d);

    /* compute the phase shift of the current dimension:
     *
     * this relies on a simplification of the covariance function
     * based on the known hypercomplex phase relationships of
     * quadrature sinusoids.
     *
     * in phase ==> cos(x1 - x2)
     * out of phase:
     *  p1 imag ==>  sin(x1 - x2)
     *  p2 imag ==> -sin(x1 - x2)
     */
    const unsigned int pd = 1 << d;
    const unsigned int d1 = p1 & pd;
    const unsigned int d2 = p2 & pd;
    const int zd = (d1 == d2 ? 0 : d1 ? -1 : 1);

    /* include the contribution from the current dimension. */
    fcov *= quad_cos(f, d, xd, zd);
  }

  /* return the computed covariance. */
  return fcov;
}

/* quad_diff_mean(): evaluate the quadrature factor mean gradient.
 *  - see factor_diff_mean_fn() for more information.
 */
FACTOR_DIFF_MEAN (quad) {
  /* initialize the gradient vector. */
  vector_set_all(df, 1.0);

  /* loop over the dimensions to compute expectations. */
  for (unsigned int d = 0; d < f->D; d++) {
    /* compute the expectation along the current dimension. */
    const double xd = vector_get(x, f->d + d);
    const double Vd = quad_cos(f, d, xd, VTAB(p, i, d));

    /* include the expectation in all but the current dimension. */
    for (unsigned int d2 = 0; d2 < f->D; d2++) {
      if (d2 != d)
        vector_set(df, d2, vector_get(df, d2) * Vd);
    }
  }

  /* loop over the dimensions to include derivatives. */
  for (unsigned int d = 0; d < f->D; d++) {
    /* get the input value and parameter indices. */
    const double xd = vector_get(x, f->d + d);
    const unsigned int pmu = 2 * d + P_MU;
    const unsigned int ptau = 2 * d + P_TAU;

    /* compute the parameter partial derivatives. */
    double dmu, dtau;
    quad_diff_cos(f, d, xd, VTAB(p, i, d), &dmu, &dtau);

    /* include the derivatives in the gradient vector. */
    vector_set(df, pmu, vector_get(df, pmu) * dmu);
    vector_set(df, ptau, vector_get(df, ptau) * dtau);
  }
}

/* quad_diff_var(): evaluate the quadrature factor variance gradient.
 *  - see factor_diff_var_fn() for more information.
 */
FACTOR_DIFF_VAR (quad) {
  /* initialize the gradient vector. */
  vector_set_all(df, 1.0);

  /* loop over the dimensions to compute expectations. */
  for (unsigned int d = 0; d < f->D; d++) {
    /* compute the expectation along the current dimension. */
    const double xd = vector_get(x, f->d + d);
    const double Vd = quad_cospair(f, d, xd, xd,
                                   VTAB(p, i, d),
                                   VTAB(p, j, d));

    /* include the expectation in all but the current dimension. */
    for (unsigned int d2 = 0; d2 < f->D; d2++) {
      if (d2 != d)
        vector_set(df, d2, vector_get(df, d2) * Vd);
    }
  }

  /* loop over the dimensions to include derivatives. */
  for (unsigned int d = 0; d < f->D; d++) {
    /* get the input value and parameter indices. */
    const double xd = vector_get(x, f->d + d);
    const unsigned int pmu = 2 * d + P_MU;
    const unsigned int ptau = 2 * d + P_TAU;

    /* compute the parameter partial derivatives. */
    double dmu, dtau;
    quad_diff_cospair(f, d, xd, xd,
                      VTAB(p, i, d),
                      VTAB(p, j, d),
                      &dmu, &dtau);

    /* include the derivatives in the gradient vector. */
    vector_set(df, pmu, vector_get(df, pmu) * dmu);
    vector_set(df, ptau, vector_get(df, ptau) * dtau);
  }
}

/* quad_meanfield(): perform a mean-field update of a quadrature factor.
 *  - see factor_meanfield_fn() for more information.
 */
FACTOR_MEANFIELD (quad) {
  /* get the extended structure pointer. */
  quadrature_t *fx = (quadrature_t*) f;

  /* check for initialization calls. */
  if (FACTOR_MEANFIELD_INIT) {
    /* zero the contents of the fft matrices. */
    nq = 0;
    matrix_set_zero(q);
    matrix_set_zero(Q);

    /* return success. */
    return 1;
  }

  /* check for finalization calls. */
  if (FACTOR_MEANFIELD_END) {
    /* finalize each dimension. */
    for (unsigned int d = 0; d < f->D; d++) {
      /* get views into the relevant q-matrix rows. */
      vector_view_t Ga = matrix_row(q, d);
      vector_view_t Xi = matrix_row(q, f->D + d);
      vector_view_t Nu = matrix_row(q, 2 * f->D + d);

      /* get a view into the relevant Q-matrix row. */
      vector_view_t Qd = matrix_row(Q, d);
      const unsigned int n = Qd.len;
      double *Qdat = Qd.data;

      /* get the prior parameters. */
      const double mu0 = vector_get(fp->par, 2 * d + P_MU);
      const double tau0 = vector_get(fp->par, 2 * d + P_TAU);

      /* execute the fft. */
      fftw_execute_dft(plan, (fftw_complex*) Qdat, (fftw_complex*) Qdat);

      /* shift the array to place zero-frequency in the center. */
      const unsigned int nh = n / 2;
      for (unsigned int i = 0; i < nh; i += 2) {
        const double qi = Qdat[i];
        Qdat[i] = Qdat[i + nh];
        Qdat[i + nh] = qi;
      }

      /* initialize the search parameters. */
      double omega_max = mu0;
      double Qmax = -1.0e9;

      /* perform an initial search for the maximum log-probability. */
      for (unsigned int i = 0; i < n; i += 2) {
        /* compute the frequency and log-probability of the point. */
        const double omega = M_PI * ((double) i / (double) nh - 1.0);
        Qdat[i] -= 0.5 * tau0 * pow(omega - mu0, 2.0);
        Qdat[i + 1] = 0.0;

        /* remember point with greatest log-probability */
        if (Qdat[i] > Qmax) {
          omega_max = omega;
          Qmax = Qdat[i];
        }
      }

      /* run a few rounds of newton-rhapson to find the best mode,
       * then perform a laplace approximation to find the precision.
       */
      double mu, tau;
      quad_costune(&Ga, &Xi, &Nu, omega_max, mu0, tau0, &mu, &tau);

      /* set the new parameters. */
      factor_set(f, 2 * d + P_MU, mu);
      factor_set(f, 2 * d + P_TAU, tau);
    }

    /* return success. */
    return 1;
  }

  /* declare variables for holding reduced coefficients. */
  double gamma, xi, nu, yr, yi;
  unsigned int I;

  /* gain access to the mean-field coefficients. */
  vector_t *a = fx->a;
  vector_t *x = fx->x;
  matrix_t *C = fx->C;
  matrix_t *Z = fx->Z;

  /* recast the amplitude and phase matrices as vectors. */
  vector_view_t cv = vector_view_array(C->data, f->K * f->K);
  vector_view_t zv = vector_view_array(Z->data, f->K * f->K);

  /* update along each dimension. */
  for (unsigned int d = 0; d < f->D; d++) {
    /* get views into the relevant q-matrix rows. */
    vector_view_t Ga = matrix_row(q, d);
    vector_view_t Xi = matrix_row(q, f->D + d);
    vector_view_t Nu = matrix_row(q, 2 * f->D + d);

    /* get a view into the relevant Q-matrix row. */
    vector_view_t Qd = matrix_row(Q, d);

    /* get the current dimension input. */
    const double xd = vector_get(dat->x, d);

    /* copy the amplitude coefficients. */
    vector_copy(a, b);
    matrix_copy(C, B);

    /* compute the phase offset coefficients. */
    for (unsigned int i = 0; i < f->K; i++) {
      /* compute the first-order phase offsets. */
      vector_set(x, i, M_PI_2 * VTAB(dat->p, i, d));

      /* compute the second-order phase offsets. */
      for (unsigned int j = 0; j < f->K; j++)
        matrix_set(Z, i, j, M_PI_2 * (VTAB(dat->p, i, d) +
                                      VTAB(dat->p, j, d)));
    }

    /* include expectations from the other dimensions. */
    for (unsigned int delta = 0; delta < f->D; delta++) {
      /* skip the current dimension. */
      if (delta == d)
        continue;

      /* get the current dimension input. */
      const double xdelta = vector_get(dat->x, delta);

      /* loop to include first-order expectations. */
      for (unsigned int i = 0; i < f->K; i++) {
        /* compute the first-order expectation. */
        const double phi1 = quad_cos(f, delta, xdelta,
                                     VTAB(dat->p, i, delta));

        /* update the amplitude vector. */
        vector_set(a, i, vector_get(a, i) * phi1);

        /* loop again to include second-order expectations. */
        for (unsigned int j = 0; j < f->K; j++) {
          /* compute the second-order expectation. */
          const double phi2 = quad_cospair(f, delta, xdelta, xdelta,
                                           VTAB(dat->p, i, delta),
                                           VTAB(dat->p, j, delta));

          /* update the amplitude matrix. */
          matrix_set(C, i, j, matrix_get(C, i, j) * phi2);
        }
      }
    }

    /* include the factor of one half from trigonometric identities. */
    blas_dscal(0.5, &cv);

    /* reduce the first-order vectors. */
    quad_cosreduce(a, x, &xi, &nu);
    gamma = xd;

    /* store the reduced parameters into their raw vectors. */
    vector_set(&Ga, nq, gamma);
    vector_set(&Xi, nq, xi);
    vector_set(&Nu, nq, nu);

    /* accumulate the reduced parameters into the fft vector. */
    I = 2 * ((unsigned int) gamma);
    yr = vector_get(&Qd, I)     + xi * cos(nu);
    yi = vector_get(&Qd, I + 1) + xi * sin(nu);
    vector_set(&Qd, I,     yr);
    vector_set(&Qd, I + 1, yi);

    /* reduce the second-order matrices. */
    quad_cosreduce(&cv, &zv, &xi, &nu);
    gamma = 2.0 * xd;

    /* store the reduced parameters in their raw vectors. */
    vector_set(&Ga, nq + 1, gamma);
    vector_set(&Xi, nq + 1, xi);
    vector_set(&Nu, nq + 1, nu);

    /* accumulate the reduced parameters into the fft vector. */
    I = 2 * ((unsigned int) gamma);
    yr = vector_get(&Qd, I)     + xi * cos(nu);
    yi = vector_get(&Qd, I + 1) + xi * sin(nu);
    vector_set(&Qd, I,     yr);
    vector_set(&Qd, I + 1, yi);
  }

  /* increment the coefficient count for all dimensions. */
  nq += 2;

  /* return success. */
  return 1;
}

/* quad_div(): evaluate the quadrature factor divergence.
 *  - see factor_div_fn() for more information.
 */
FACTOR_DIV (quad) {
  /* initialize the divergence computation. */
  double div = 0.0;

  /* include the divergences along each dimension. */
  for (unsigned int d = 0, p = 0; d < f->D; d++, p += 2) {
    /* get the first factor parameters. */
    const double mu = vector_get(f->par, p + P_MU);
    const double tau = vector_get(f->par, p + P_TAU);

    /* get the second factor parameters. */
    const double mu2 = vector_get(f2->par, p + P_MU);
    const double tau2 = vector_get(f2->par, p + P_TAU);

    /* include the current dimension's divergence. */
    div += 0.5 * tau2 * (mu * mu + 1.0 / tau - 2.0 * mu * mu2 + mu2 * mu2)
         - 0.5 * log(tau2 / tau) - 0.5;
  }

  /* return the computed divergence. */
  return div;
}

/* quad_init(): initialize the quadrature factor structure.
 *  - see factor_init_fn() for more information.
 */
FACTOR_INIT (quad) {
  /* initialize the phase table. */
  quadrature_t *fx = (quadrature_t*) f;
  fx->vtab = NULL;

  /* initialize the mean-field variables. */
  fx->a = fx->x = NULL;
  fx->C = fx->Z = NULL;

  /* return success. */
  return 1;
}

/* quad_resize(): handle resizes of the quadrature factor.
 *  - see factor_resize_fn() for more information.
 */
FACTOR_RESIZE (quad) {
  /* get the extended structure pointer. */
  quadrature_t *fx = (quadrature_t*) f;

  /* declare variables used for allocation. */
  unsigned int nblk;
  char *ptr = NULL;

  /* determine the amount of memory required. */
  nblk = vtab_bytes(D, K);
  nblk += 2 * vector_bytes(K);
  nblk += 2 * matrix_bytes(K, K);

  /* reallocate the phase table. */
  fx->vtab = realloc(fx->vtab, nblk);
  if (!fx->vtab)
    return 0;

  /* initialize the phase table. */
  vtab_init(fx->vtab, D, K);

  /* initialize the memory address. */
  ptr = (char*) fx->vtab;
  ptr += vtab_bytes(D, K);

  /* initialize the mean-field vector of amplitudes. */
  fx->a = (vector_t*) ptr;
  vector_init(fx->a, K);
  ptr += vector_bytes(K);

  /* initialize the mean-field vector of phase shifts. */
  fx->x = (vector_t*) ptr;
  vector_init(fx->x, K);
  ptr += vector_bytes(K);

  /* initialize the mean-field matrix of amplitudes. */
  fx->C = (matrix_t*) ptr;
  matrix_init(fx->C, K, K);
  ptr += matrix_bytes(K, K);

  /* initialize the mean-field matrix of phase shifts. */
  fx->Z = (matrix_t*) ptr;
  matrix_init(fx->Z, K, K);

  /* reallocate the parameter name table. */
  const unsigned int sz = P * (sizeof(char*) + 8);
  f->parnames = realloc(f->parnames, sz);
  if (!f->parnames)
    return 0;

  /* initialize the parameter name addresses. */
  ptr = (char*) (f->parnames + P);
  for (unsigned int p = 0; p < P; p++) {
    f->parnames[p] = (char*) ptr;
    ptr += 8;
  }

  /* initialize the parameter name values. */
  for (unsigned int d = 0; d < D; d++) {
    snprintf(f->parnames[2 * d + P_MU],  7, "mu%u", d);
    snprintf(f->parnames[2 * d + P_TAU], 7, "tau%u", d);
  }

  /* return success. */
  return 1;
}

/* quad_kernel(): write the kernel code of a quadrature factor.
 *  - see factor_kernel_fn() for more information.
 */
FACTOR_KERNEL (quad) {
  /* define kernel code format strings. */
  const char *fmtA = "cov = 1.0;\n";
  const char *fmtB = "{\n\
const double xd = x1[%u] - x2[%u];\n\
const double mu = par[%u];\n\
const double tau = par[%u];\n\
const unsigned int d1 = p1 & %u;\n\
const unsigned int d2 = p2 & %u;\n\
const double zd = (d1 == d2 ? 0.0 : d1 ? -%lf : %lf);\n\
cov *= exp(-0.5 * xd * xd / tau) * cos(mu * xd + zd);\n\
}\n";

  /* allocate the kernel code string. */
  const unsigned int len = strlen(fmtA) + f->D * strlen(fmtB);
  char *kstr = malloc(len + 24);
  if (!kstr)
    return NULL;

  /* write the header. */
  char *pos = kstr;
  pos += sprintf(pos, fmtA);

  /* write the dimension strings. */
  for (unsigned int d = 0, p = p0; d < f->D; d++, p += 2)
    pos += sprintf(pos, fmtB, f->d + d, f->d + d,
                   p + P_MU, p + P_TAU,
                   1 << d, 1 << d,
                   M_PI_2, M_PI_2);

  /* return the new string. */
  return kstr;
}

/* quad_set(): store a parameter into a quadrature factor.
 *  - see factor_set_fn() for more information.
 */
FACTOR_SET (quad) {
  /* determine which parameter is being assigned. */
  switch (i % 2) {
    /* frequency mean: in (-inf, inf) */
    case P_MU:;
      const double mu = value;
      vector_set(f->par, i, mu);

      return 1;

    /* frequency precision: in (0, inf) */
    case P_TAU:;
      const double tau = value;
      if (tau <= 0.0)
        return 0;

      vector_set(f->par, i, tau);

      matrix_set(f->inf, i - 1, i - 1, tau);
      matrix_set(f->inf, i, i, 0.75 / (tau * tau));

      return 1;
  }

  /* invalid parameter index. */
  return 0;
}

/* quad_free(): free extra information from quadrature factors.
 *  - see factor_free_fn() for more information.
 */
FACTOR_FREE (quad) {
  /* free the phase table. */
  quadrature_t *fx = (quadrature_t*) f;
  free(fx->vtab);

  /* free the parameter names table. */
  free(f->parnames);
}

/* quad_set_dims(): set the dimensionality of a quadrature factor.
 *
 * arguments:
 *  @f: factor structure pointer to modify.
 *  @D: new factor dimensionality.
 *
 * returns:
 *  integer indicating resize success (1) or failure (0).
 */
int quad_set_dims (factor_t *f, const unsigned int D) {
  /* resize the factor to accomodate the correct number of
   * dimensions, parameters, and weights.
   */
  if (!factor_resize(f, D, 2 * D, 1 << D) || !quad_set_ftsize(f, Nft))
    return 0;

  /* return success. */
  return 1;
}

/* quad_set_ftsize(): set the number of complex points used by
 * the fast fourier transform during mean-field inference.
 *
 * arguments:
 *  @f: factor structure pointer to modify.
 *  @n: number of complex points.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int quad_set_ftsize (factor_t *f, const unsigned int n) {
  /* determine the new matrix sizes. */
  const unsigned int Dnew = (f && f->D > Dft ? f->D : Dft);
  const unsigned int Nnew = (n > Nft ? n : Nft);

  /* check if the point count is unchanged at zero. */
  if (Nft == 0 && Nnew == Nft)
    return 1;

  /* check if the sizes have not changed. */
  if (Dnew == Dft && Nnew == Nft)
    return 1;

  /* reset the counts. */
  Dft = Dnew;
  Nft = Nnew;
  nq = 0;

  /* free the fft plan and matrices. */
  fftw_destroy_plan(plan);
  free(q);

  /* determine the memory requirement of the mean-field variables. */
  unsigned int bytes = matrix_bytes(3 * Dft, 2 * Nft)
                     + matrix_bytes(Dft, 2 * Nft);

  /* allocate a new set of fft matrices. */
  char *ptr = malloc(bytes);
  if (!ptr)
    return 0;

  /* initialize the first matrix. */
  q = (matrix_t*) ptr;
  matrix_init(q, 3 * Dft, 2 * Nft);
  ptr += matrix_bytes(3 * Dft, 2 * Nft);

  /* initialize the second matrix. */
  Q = (matrix_t*) ptr;
  matrix_init(Q, Dft, 2 * Nft);

  /* generate a new fft plan. */
  vector_view_t qv = matrix_row(q, 0);
  plan = fftw_plan_dft_1d(n, (fftw_complex*) qv.data,
                             (fftw_complex*) qv.data,
                             FFTW_BACKWARD, FFTW_MEASURE);

  /* failed to plan? */
  if (!plan)
    return 0;

  /* return success. */
  return 1;
}

/* --- */

/* quad_getprops(): return a list of means or precisions
 * from a multidimensional quadrature factor.
 *
 * arguments:
 *  @f: factor structure pointer.
 *  @p0: parameter offset.
 *
 * returns:
 *  newly allocated and initialized list of parameters,
 *  with one entry for each dimension.
 */
static object_t *quad_getprops (const factor_t *f, const unsigned int p0) {
  /* return floats from one-dimensional factors. */
  if (f->D == 1)
    return (object_t*) float_alloc_with_value(factor_get(f, p0));

  /* allocate a list for parameters. */
  list_t *lst = list_alloc_with_length(f->D);
  if (!lst)
    return NULL;

  /* store the list elements. */
  for (unsigned int d = 0, p = p0; d < f->D; d++, p += 2) {
    /* allocate a float for the current list element. */
    flt_t *elem = float_alloc_with_value(factor_get(f, p));
    if (!elem) {
      obj_release((object_t*) lst);
      return NULL;
    }

    /* store the float in the list. */
    list_set(lst, d, (object_t*) elem);
  }

  /* return the new list. */
  return (object_t*) lst;
}

/* quad_setprops(): set the means or precisions
 * of a multidimensional quadrature factor.
 *
 * arguments:
 *  @f: factor structure pointer.
 *  @val: float or list of values.
 *  @p0: parameter offset.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
static int quad_setprops (factor_t *f, object_t *val,
                          const unsigned int p0) {
  /* handle one-dimensional assignments. */
  if (f->D == 1) {
    /* admit only numbers for one-dimensional factors. */
    if (!OBJECT_IS_NUM(val))
      return 0;

    /* set the factor parameter. */
    return factor_set(f, p0, num_get(val));
  }

  /* admit only lists for multidimensional factors. */
  if (!OBJECT_IS_LIST(val))
    return 0;

  /* check that the list length matches the factor dimensionality. */
  list_t *lst = (list_t*) val;
  if (lst->len != (size_t) f->D)
    return 0;

  /* loop over the list elements. */
  for (unsigned int d = 0, p = p0; d < f->D; d++, p += 2) {
    /* admit only numbers as list elements. */
    object_t *elem = list_get(lst, d);
    if (!OBJECT_IS_NUM(elem))
      return 0;

    /* set the factor parameter. */
    if (!factor_set(f, p, num_get(elem)))
      return 0;
  }

  /* return success. */
  return 1;
}

/* quad_getprop_mu(): get a quadrature factor mean parameter.
 *  - see vfl/object_getprop_fn() for details.
 */
static object_t *quad_getprop_mu (const factor_t *f) {
  /* return either a float or a list. */
  return quad_getprops(f, P_MU);
}

/* quad_getprop_tau(): get a quadrature factor precision parameter.
 *  - see vfl/object_getprop_fn() for details.
 */
static object_t *quad_getprop_tau (const factor_t *f) {
  /* return either a float or a list. */
  return quad_getprops(f, P_TAU);
}

/* quad_getprop_dims(): get a quadrature factor dimensionality.
 *  - see vfl/object_getprop_fn() for details.
 */
static int_t *quad_getprop_dims (factor_t *f) {
  /* return the factor dimensionality as a new integer. */
  return int_alloc_with_value(f->D);
}

/* quad_getprop_ftsize(): get a quadrature factor fourier transform size.
 *  - see vfl/object_getprop_fn() for details.
 */
static int_t *quad_getprop_ftsize (factor_t *f) {
  /* return the fourier transform size as a new integer. */
  int_t *iobj = int_alloc_with_value(Nft);
  return (f ? iobj : iobj);
}

/* quad_setprop_mu(): set a quadrature factor mean parameter.
 *  - see vfl/object_setprop_mu() for details.
 */
static int quad_setprop_mu (factor_t *f, object_t *val) {
  /* set the mean parameter(s). */
  return quad_setprops(f, val, P_MU);
}

/* quad_setprop_tau(): set a quadrature factor precision parameter.
 *  - see vfl/object_setprop_mu() for details.
 */
static int quad_setprop_tau (factor_t *f, object_t *val) {
  /* set the precision parameter(s). */
  return quad_setprops(f, val, P_TAU);
}

/* quad_setprop_dims(): set a quadrature factor dimensionality.
 *  - see vfl/object_setprop_fn() for details.
 */
static int quad_setprop_dims (factor_t *f, object_t *val) {
  /* admit only integer values. */
  if (!OBJECT_IS_INT(val))
    return 0;

  /* admit only positive dimensionalities. */
  const long ival = int_get((int_t*) val);
  if (ival <= 0)
    return 0;

  /* set the dimensionality and return. */
  return quad_set_dims(f, ival);
}

/* quad_setprop_ftsize(): set a quadrature factor fourier transform size.
 *  - see vfl/object_setprop_fn() for details.
 */
static int quad_setprop_ftsize (factor_t *f, object_t *val) {
  /* admit only integer values. */
  if (!OBJECT_IS_INT(val))
    return 0;

  /* admit only positive sizes. */
  const long ival = int_get((int_t*) val);
  if (ival <= 0)
    return 0;

  /* set the fourier transform size and return. */
  return quad_set_ftsize(f, ival);
}

/* quad_properties: array of accessible object properties.
 */
static object_property_t quad_properties[] = {
  FACTOR_PROP_BASE,
  { "mu",
    (object_getprop_fn) quad_getprop_mu,
    (object_setprop_fn) quad_setprop_mu
  },
  { "tau",
    (object_getprop_fn) quad_getprop_tau,
    (object_setprop_fn) quad_setprop_tau
  },
  { "dims",
    (object_getprop_fn) quad_getprop_dims,
    (object_setprop_fn) quad_setprop_dims
  },
  { "ftsize",
    (object_getprop_fn) quad_getprop_ftsize,
    (object_setprop_fn) quad_setprop_ftsize
  },
  { NULL, NULL, NULL }
};

/* --- */

/* quad_methods: array of callable object methods.
 */
static object_method_t quad_methods[] = {
  FACTOR_METHOD_BASE,
  { NULL, NULL }
};

/* --- */

/* quad_type: quadrature factor type structure.
 */
static factor_type_t quad_type = {
  { /* base: */
    "quad",                                      /* name      */
    sizeof(quadrature_t),                        /* size      */

    (object_init_fn) factor_init,                /* init      */
    (object_copy_fn) factor_copy,                /* copy      */
    (object_free_fn) factor_free,                /* free      */
    NULL,                                        /* test      */
    NULL,                                        /* cmp       */

    (object_binary_fn) factor_add,               /* add       */
    NULL,                                        /* sub       */
    (object_binary_fn) factor_mul,               /* mul       */
    NULL,                                        /* div       */
    NULL,                                        /* pow       */

    NULL,                                        /* get       */
    NULL,                                        /* set       */
    quad_properties,                             /* props     */
    quad_methods                                 /* methods   */
  },

  1,                                             /* initial D */
  2,                                             /* initial P */
  2,                                             /* initial K */
  NULL,                                          /* parnames  */
  quad_eval,                                     /* eval      */
  quad_mean,                                     /* mean      */
  quad_var,                                      /* var       */
  quad_cov,                                      /* cov       */
  quad_diff_mean,                                /* diff_mean */
  quad_diff_var,                                 /* diff_var  */
  quad_meanfield,                                /* meanfield */
  quad_div,                                      /* div       */
  quad_init,                                     /* init      */
  quad_resize,                                   /* resize    */
  quad_kernel,                                   /* kernel    */
  quad_set,                                      /* set       */
  NULL,                                          /* copy      */
  quad_free                                      /* free      */
};

/* factor_type_quad: address of the quad_type structure. */
const factor_type_t *vbnmr_factor_quad = &quad_type;

