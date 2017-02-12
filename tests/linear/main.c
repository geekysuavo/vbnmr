
/* include the core header. */
#include <vbnmr/core.h>

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

  /* create some auxiliary variables. */
  vector_t *t = vector_alloc(1);
  rng_t *R = rng_alloc();

  /* compute a random phase shift. */
  const double phi = M_PI * (2.0 * rng_uniform(R) - 1.0);

  /* create a model for generating simulated measurements. */
  model_t *mdlA = model_alloc(1, 5);
  model_set_sigma(mdlA, 1.0);
  model_set_delta(mdlA, 1.0e3);
  for (unsigned int j = 0; j < mdlA->M; j++)
    model_prior_signal(mdlA, j, 0, 0.0, 1.0e-5, 100.0, 100.0);
  model_set_from_priors(mdlA);
  model_randomize_mu(mdlA, R);
  for (unsigned int j = 0; j < mdlA->M; j++) {
    const double aj = 1.0 + 99.0 * rng_uniform(R);
    vector_set(mdlA->ahat, j, aj * cos(phi));
    vector_set(mdlA->ahat, j + mdlA->M, aj * sin(phi));
  }

  /* create a model for learning signal parameters. */
  model_t *mdlB = model_alloc(1, mdlA->M);
  model_set_sigma(mdlB, mdlA->sigma);
  model_set_delta(mdlB, mdlA->delta);
  for (unsigned int j = 0; j < mdlB->M; j++)
    model_prior_signal(mdlB, j, 0, mdlA->sig[j][0].mu, 0.1, 100.0, 100.0);
  model_set_from_priors(mdlB);

  /* create a uniform sampling grid. */
  double grid_meas[] = { 0.0, 1.0e-3, 0.2 };
  matrix_view_t grid = matrix_view_array(grid_meas, 1, 3);
  dataset_t *dat = dataset_alloc_from_grid(&grid);
  model_eval(mdlA, R, dat);
  dataset_fwrite(dat, "meas.dat");

  /* optimize the model parameters. */
  optim_t *opt = optim_alloc(mdlB, dat);
  opt->log = fopen("par.dat", "w");
  opt->l0 = 1.1;
  opt->dl = 0.1;
  opt->iter_max = 500;
  opt->step_max = 10;
  optim_execute(opt);
  const double elbo = opt->elbo;
  optim_free(opt);

  fprintf(stderr, "elbo:%16.9le\n", elbo);
  for (unsigned int j = 0; j < mdlA->M; j++) {
    const double Aa = sqrt(pow(vector_get(mdlA->ahat, j), 2.0) +
                           pow(vector_get(mdlA->ahat, mdlA->M + j), 2.0));
    const double Ba = sqrt(pow(vector_get(mdlB->ahat, j), 2.0) +
                           pow(vector_get(mdlB->ahat, mdlB->M + j), 2.0));

    const double As = sqrt(matrix_get(mdlA->Sigma, j, j) +
                           matrix_get(mdlA->Sigma, j + mdlA->M, j + mdlA->M));
    const double Bs = sqrt(matrix_get(mdlB->Sigma, j, j) +
                           matrix_get(mdlB->Sigma, j + mdlB->M, j + mdlB->M));

    printf("%16.9le %16.9le %16.9le %16.9le\n", Aa, As, Ba, Bs);
  }

  /* allocate datasets for storing predictions. */
  double grid_model[] = { 0.0, 1.0e-4, 0.5 };
  grid = matrix_view_array(grid_model, 1, 3);
  dataset_t *mean = dataset_alloc_from_grid(&grid);
  dataset_t *var = dataset_alloc_from_grid(&grid);
  dataset_t *map = dataset_alloc_from_grid(&grid);

  /* compute the predictions. */
  model_eval(mdlB, NULL, map);
  model_predict(mdlB, mean, var);

  /* write the prediction datasets to files. */
  dataset_fwrite(mean, "mean.dat");
  dataset_fwrite(var, "var.dat");
  dataset_fwrite(map, "map.dat");

  /* free all allocated structures. */
  dataset_free(mean);
  dataset_free(var);
  dataset_free(map);
  dataset_free(dat);
  model_free(mdlA);
  model_free(mdlB);
  vector_free(t);
  rng_free(R);

  /* return success. */
  return 0;
}

