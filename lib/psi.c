
/* include the polygamma header. */
#include <vbnmr/psi.h>

/* set the coefficient count of the series approximations. */
#define N_PSI 7

/* coeff_digamma: taylor series coefficients of psi(x) at x=inf.
 */
static const double coeff_digamma[N_PSI] = {
 -8.3333333333333329e-02,
  8.3333333333333332e-03,
 -3.9682539682539680e-03,
  4.1666666666666666e-03,
 -7.5757575757575760e-03,
  2.1092796092796094e-02,
 -8.3333333333333329e-02
};

/* coeff_trigamma: taylor series coefficients of psi(1,x) at x=inf.
 */
static const double coeff_trigamma[N_PSI] = {
  1.0000000000000000e+00,
  1.6666666666666666e-01,
 -3.3333333333333333e-02,
  2.3809523809523808e-02,
 -3.3333333333333333e-02,
  7.5757575757575760e-02,
 -2.5311355311355310e-01
};

/* digamma(): compute the logarithmic first derivative of the
 * gamma function.
 *
 * arguments:
 *  @z: function argument.
 *
 * returns:
 *  value of psi(0,@z).
 */
double digamma (const double z) {
  /* declare required variables:
   *  @psi: function value.
   *  @x: shifted argument.
   *  @xn: powers of @x.
   */
  double psi, x, xn;

  /* initialize the argument and result. */
  psi = 0.0;
  x = z;

  /* shift small arguments up for better accuracy. */
  while (x < 10.0) {
    psi -= 1.0 / x;
    x += 1.0;
  }

  /* initialize the series computation. */
  const double xsq = x * x;
  psi += log(x) - 0.5 / x;
  xn = xsq;

  /* compute the series. */
  for (unsigned int n = 0; n < N_PSI; n++) {
    psi += coeff_digamma[n] / xn;
    xn *= xsq;
  }

  /* return the computed value. */
  return psi;
}

/* trigamma(): compute the logarithmic second derivative of the
 * gamma function.
 *
 * arguments:
 *  @z: function argument.
 *
 * returns:
 *  value of psi(1,@z).
 */
double trigamma (const double z) {
  /* declare required variables:
   *  @psi: function value.
   *  @x: shifted argument.
   *  @xn: powers of @x.
   */
  double psi, x, xn;

  /* initialize the argument and result. */
  psi = 0.0;
  x = z;

  /* shift small arguments up for better accuracy. */
  while (x < 10.0) {
    psi += 1.0 / (x * x);
    x += 1.0;
  }

  /* initialize the series computation. */
  const double xsq = x * x;
  psi += 0.5 / xsq;
  xn = x;

  /* compute the series. */
  for (unsigned int n = 0; n < N_PSI; n++) {
    psi += coeff_trigamma[n] / xn;
    xn *= xsq;
  }

  /* return the computed value. */
  return psi;
}

