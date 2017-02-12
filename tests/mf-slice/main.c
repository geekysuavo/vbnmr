
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

  /* create a model for learning signal parameters. */
  model_t *mdl = model_alloc(1, 10);
  model_set_sigma(mdl, 0.01);
  model_set_delta(mdl, 100.0);
  for (unsigned int j = 0; j < mdl->M; j++)
    model_prior_signal(mdl, j, 0, 0.0, 1.0e-6, 100.0, 10.0);
  model_set_from_priors(mdl);
  model_randomize_mu(mdl, R);

  /* open a dataset for holding measurements. */
  dataset_t *dat = dataset_alloc(1);
  dataset_fread(dat, "meas.dat");

  /* use a single mean-field step to learn the model parameters. */
  optim_t *opt = optim_alloc(mdl, dat);
  optim_init(opt);
  optim_free(opt);

  /* allocate datasets for storing predictions. */
  double grid_values[] = { 0.0, 1.0e-4, 1.0 };
  matrix_view_t grid = matrix_view_array(grid_values, 1, 3);
  dataset_t *mean = dataset_alloc_from_grid(&grid);
  dataset_t *var = dataset_alloc_from_grid(&grid);

  /* compute the posterior predictions. */
  model_predict(mdl, mean, var);

  /* write the prediction datasets to files. */
  dataset_fwrite(mean, "mean.dat");
  dataset_fwrite(var, "var.dat");

  /* free all allocated structures. */
  dataset_free(mean);
  dataset_free(var);
  dataset_free(dat);
  model_free(mdl);
  vector_free(t);
  rng_free(R);

  /* return success. */
  return 0;
}

