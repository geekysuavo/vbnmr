
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
  vector_t *t = vector_alloc(2);
  rng_t *R = rng_alloc();

  /* set up a matrix that defines the sampling grid. */
  double grid_array[] = { 0.0, 1.0e-4, 0.02,
                          0.0, 1.0e-3, 0.1 };
  matrix_view_t grid = matrix_view_array(grid_array, 2, 3);

  /* create a model for generating simulated measurements. */
  model_t *mdlA = model_alloc(2, 1);

  /* set the noise and amplitude priors. */
  model_set_sigma(mdlA, 1.0);
  model_set_delta(mdlA, 1.0e3);

  model_prior_signal(mdlA,  0, 0, 1000.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA,  0, 1, -340.000, 1.0, 10.0, 10.0);

/*
  model_prior_signal(mdlA,  1, 0, 1900.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA,  1, 1,  133.000, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA,  2, 0, 2250.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA,  2, 1,  -42.000, 1.0, 10.0, 10.0);
*/

  model_set_from_priors(mdlA);

  for (unsigned int j = 0; j < mdlA->M; j++) {
    const double aj = 50.0 + 100.0 * ((double) j / (double) mdlA->M);
    vector_set(mdlA->ahat, j, aj);
  }

  /* write the simulated data to a grid. */
  dataset_t *sim = dataset_alloc_from_grid(&grid);
  model_eval(mdlA, NULL, sim);
  dataset_fwrite(sim, "sim.dat");

  /* create a model for learning signal parameters. */
  model_t *mdlB = model_alloc(2, 1);
  model_set_sigma(mdlB, mdlA->sigma);
  model_set_delta(mdlB, mdlA->delta);
  for (unsigned int j = 0; j < mdlB->M; j++) {
    model_prior_signal(mdlB, j, 0, 2000.0, 1.0e-6, 10.0, 10.0);
    model_prior_signal(mdlB, j, 1,    0.0, 1.0e-6, 10.0, 10.0);
  }
  model_set_from_priors(mdlB);
  model_randomize_mu(mdlB, R);

  /* allocate a dataset for storing simulated measurements. */
  dataset_t *dat = dataset_alloc(2);
  dataset_fread(dat, "sched.dat");

  /* evaluate the model on the dataset, with gaussian noise. */
  model_eval(mdlA, R, dat);

  /* write the input dataset to a file. */
  dataset_fwrite(dat, "meas.dat");

  /* optimize the model parameters. */
  optim_t *opt = optim_alloc(mdlB, dat);
  opt->log = fopen("par.dat", "w");
  opt->l0 = 0.01;
  opt->dl = 0.1;
  opt->step_max = 10;
  optim_execute(opt);
  optim_free(opt);

  /* write the fitted data values. */
  dataset_t *mmean = dataset_alloc(2);
  dataset_t *mvar = dataset_alloc(2);
  dataset_fread(mmean, "sched.dat");
  dataset_fread(mvar, "sched.dat");
  model_predict(mdlB, mmean, mvar);
  dataset_fwrite(mmean, "fit.dat");

  /* allocate datasets for storing predictions. */
  dataset_t *mean = dataset_alloc_from_grid(&grid);
  dataset_t *var = dataset_alloc_from_grid(&grid);
  dataset_t *map = dataset_alloc_from_grid(&grid);

  /* compute the gp predictions. */
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

