
/* include the quad header. */
#include <vbnmr/quad.h>

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
#define VTAB(k,l,d) ((quad_t*) (f))->vtab[k][l][d]

/* quad_t: structure for holding a quadrature factor.
 */
typedef struct {
  /* @base: core factor structure members. */
  factor_t base;

  /* dimensionality-dependent members:
   *  @vtab: quadrature modulation table.
   */
  int ***vtab;
}
quad_t;

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

/* --- */

/* quad_mean(): evaluate the quadrature factor mean.
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
 */
FACTOR_VAR (quad) {
  /* initialize the variance computation. */
  double fvar = 1.0;

  /* include the contributions along each dimension. */
  for (unsigned int d = 0; d < f->D; d++) {
    const double xd = vector_get(x, f->d + d);
    fvar *= quad_cospair(f, d, xd, xd,
                               VTAB(p, i, d),
                               VTAB(p, j, d));
  }

  /* return the computed variance. */
  return fvar;
}

/* quad_diff_mean(): evaluate the quadrature factor mean gradient.
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

/* quad_div(): evaluate the quadrature factor divergence.
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
 */
FACTOR_INIT (quad) {
  /* initialize the phase table. */
  quad_t *fx = (quad_t*) f;
  fx->vtab = NULL;

  /* return success. */
  return 1;
}

/* quad_resize(): handle resizes of the quadrature factor.
 */
FACTOR_RESIZE (quad) {
  /* get the extended structure pointer. */
  quad_t *fx = (quad_t*) f;

  /* reallocate the phase table. */
  fx->vtab = realloc(fx->vtab, vtab_bytes(f->D, f->K));
  if (!fx->vtab)
    return 0;

  /* initialize the phase table. */
  vtab_init(fx->vtab, f->D, f->K);

  /* reallocate the parameter name table. */
  const unsigned int sz = f->P * (sizeof(char*) + 8);
  f->parnames = realloc(f->parnames, sz);
  if (!f->parnames)
    return 0;

  /* create a pointer for initializing the parameter names. */
  char *ptr = (char*) (f->parnames + f->P);

  /* initialize the parameter name addresses. */
  for (unsigned int p = 0; p < f->P; p++) {
    f->parnames[p] = (char*) ptr;
    ptr += 8;
  }

  /* initialize the parameter name values. */
  for (unsigned int d = 0; d < f->D; d++) {
    snprintf(f->parnames[2 * d + P_MU],  7, "mu%u", d);
    snprintf(f->parnames[2 * d + P_TAU], 7, "tau%u", d);
  }

  /* return success. */
  return 1;
}

/* quad_set(): store a parameter into a quadrature factor.
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

      matrix_set(f->inf, i, i, tau);
      matrix_set(f->inf, i + 1, i + 1, 0.75 / (tau * tau));

      return 1;
  }

  /* invalid parameter index. */
  return 0;
}

/* quad_free(): free extra information from quadrature factors.
 */
FACTOR_FREE (quad) {
  /* free the phase table. */
  quad_t *fx = (quad_t*) f;
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
   * dimensions, weights, and parameters.
   */
  return factor_resize(f, D, 1 << D, 2 * D);
}

/* quad_type: quadrature factor type structure.
 */
static factor_type_t quad_type = {
  "quad",                                        /* name      */
  sizeof(quad_t),                                /* size      */
  1,                                             /* initial D */
  2,                                             /* initial P */
  2,                                             /* initial K */
  NULL,                                          /* parnames  */
  quad_mean,                                     /* mean      */
  quad_var,                                      /* var       */
  quad_diff_mean,                                /* diff_mean */
  quad_diff_var,                                 /* diff_var  */
  NULL,                                          /* meanfield */
  quad_div,                                      /* div       */
  quad_init,                                     /* init      */
  quad_resize,                                   /* resize    */
  quad_set,                                      /* set       */
  NULL,                                          /* copy      */
  quad_free                                      /* free      */
};

/* factor_type_quad: address of the quad_type structure. */
const factor_type_t *vbnmr_factor_quad = &quad_type;

