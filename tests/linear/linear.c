
/* include the required headers. */
#include <vbnmr/vbnmr.h>
#include <vfl/util/rng.h>

/* main(): application entry point.
 *
 * arguments:
 *  @argc: number of command-line argument strings.
 *  @argv: array of command-line argument strings.
 *
 * returns:
 *  application execution return code.
 */
int main (int argc, char **argv) {
  /* dirty trick to keep the compiler from complaining. */
  if (!argc || !argv) return 1;

  /* allocate a random number generator. */
  rng_t *R = rng_alloc();

  /* define some constants:
   *  @M: number of signals.
   *  @tau: noise precision.
   *  @tau_omega: frequency precision.
   */
  const unsigned int M = 5;
  const double tau = 1.0;
  const double tau_omega = 20.0;

  /* compute derived constants:
   *  @phi: global signal phase offset.
   *  @sigma: noise standard deviation.
   *  @sigma_omega: frequency standard deviation.
   */
  const double sigma = 1.0 / sqrt(tau);
  const double sigma_omega = 1.0 / sqrt(tau_omega);
  const double phi = M_PI * (2.0 * rng_uniform(R) - 1.0);

  /* allocate a dataset holding a uniform grid. */
  double grid_values[] = { 0.0, 1.0, 1000.0 };
  matrix_view_t grid = matrix_view_array(grid_values, 1, 3);
  data_t *dat = data_alloc_from_grid(2, &grid);

  /* allocate a matrix for holding initial and final parameters. */
  matrix_t *par = matrix_alloc(M, 6);

  /* set up a hybrid model. */
  model_t *mdl = model_alloc(vbnmr_model_vfgp);
  model_set_alpha0(mdl, tau);
  model_set_nu(mdl, 1.0e-6);
  model_set_data(mdl, dat);

  /* add factors to the model. */
  for (unsigned int j = 0; j < M; j++) {
    /* create a decay factor. */
    factor_t *fR = factor_alloc(vfl_factor_decay);
    factor_set(fR, 0, 100.0);
    factor_set(fR, 1, 1.0e5);
    factor_set_fixed(fR, 1);

    /* create a quadrature factor. */
    factor_t *fV = factor_alloc(vbnmr_factor_quad);
    factor_set(fV, 0, 0.0);
    factor_set(fV, 1, tau_omega);
    quad_set_ftsize(fV, 65536);

    /* create a product of the decay and quadrature factors. */
    factor_t *f = factor_alloc(vfl_factor_product);
    product_add_factor(f, 0, fR);
    product_add_factor(f, 0, fV);

    /* add the product factor to the model. */
    model_add_factor(mdl, f);
  }

  /* set the initial factor frequencies. */
  for (unsigned int j = 0; j < mdl->M; j++) {
    /* compute a random signal frequency. */
    const double omega = sigma_omega * rng_normal(R);

    /* set the current factor mean. */
    factor_set(mdl->factors[j], 2, omega);
    matrix_set(par, j, 0, omega);
  }

  /* set the initial factor weights. */
  for (unsigned int j = 0, k = 0; j < mdl->M; j++, k += 2) {
    /* compute a random signal amplitude. */
    const double w = 1.0 + 99.0 * rng_uniform(R);
    const double w0 = w * cos(phi);
    const double w1 = w * sin(phi);

    /* set the current factor weights. */
    vector_set(mdl->wbar, k + 0, w0);
    vector_set(mdl->wbar, k + 1, w1);
    matrix_set(par, j, 1, w);
  }

  /* evaluate the model over the input dataset. */
  for (unsigned int i = 0; i < dat->N; i++) {
    datum_t *di = data_get(dat, i);
    di->y = model_eval(mdl, di->x, di->p) + sigma * rng_normal(R);
  }

  /* reset the model. */
  model_reset(mdl);

  /* optimize using an initial round of mean-field. */
  optim_t *opt = optim_alloc(vfl_optim_mf);
  optim_set_model(opt, mdl);
  optim_set_max_iters(opt, 1);
  optim_execute(opt);
  optim_free(opt);

  /* optimize further using natural gradients. */
  opt = optim_alloc(vfl_optim_fg);
  optim_set_model(opt, mdl);
  optim_execute(opt);
  optim_free(opt);

  /* loop again to extract results. */
  for (unsigned int j = 0, k = 0; j < mdl->M; j++, k += 2) {
    /* get the inferred frequency. */
    product_update(mdl->factors[j]);
    const double freq = factor_get(mdl->factors[j], 2);
    const double freq_err = 1.0 / sqrt(factor_get(mdl->factors[j], 3));

    /* compute the weight magnitude. */
    const double ampl = sqrt(pow(vector_get(mdl->wbar, k + 0), 2.0) +
                             pow(vector_get(mdl->wbar, k + 1), 2.0));
    const double ampl_err = sqrt(matrix_get(mdl->Sigma, k + 0, k + 0) +
                                 matrix_get(mdl->Sigma, k + 1, k + 1));

    /* store the results. */
    matrix_set(par, j, 2, freq);
    matrix_set(par, j, 3, ampl);
    matrix_set(par, j, 4, freq_err);
    matrix_set(par, j, 5, ampl_err);
  }

  /* output the initial and final parameters. */
  for (unsigned int r = 0; r < par->rows; r++)
    for (unsigned int c = 0; c < par->cols; c++)
      printf("%16.9le%s", matrix_get(par, r, c),
             c == par->cols - 1 ? "\n" : " ");

  /* free the structures. */
  matrix_free(par);
  model_free(mdl);
  data_free(dat);
  rng_free(R);

  /* return success. */
  return 0;
}

