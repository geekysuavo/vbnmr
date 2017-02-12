
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
  double ti;

  /* set up a matrix that defines the sampling grid. */
  double grid_array[] = { 0.0, 1.0e-4, 0.02,
                          0.0, 1.0e-3, 0.1 };
  matrix_view_t grid = matrix_view_array(grid_array, 2, 3);

  /* create a model for generating simulated measurements. */
/*
  model_t *mdlA = model_alloc(2, 58);
*/
  model_t *mdlA = model_alloc(2, 10);

  /* set the noise and amplitude priors. */
  model_set_sigma(mdlA, 1.0);
  model_set_delta(mdlA, 1.0e3);

  model_prior_signal(mdlA,  0, 0, 2166.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA,  0, 1, -340.609, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA,  1, 0, 1974.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA,  1, 1,  133.811, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA,  2, 0, 2082.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA,  2, 1,  -42.576, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA,  3, 0, 2148.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA,  3, 1,  255.457, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA,  4, 0, 2742.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA,  4, 1, -413.596, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA,  5, 0, 2556.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA,  5, 1, -152.057, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA,  6, 0, 2250.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA,  6, 1, -395.349, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA,  7, 0, 2550.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA,  7, 1, -383.185, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA,  8, 0, 2316.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA,  8, 1, -346.691, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA,  9, 0, 2442.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA,  9, 1, -395.349, 1.0, 10.0, 10.0);

/*
  model_prior_signal(mdlA, 10, 0, 1848.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 10, 1,  632.559, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 11, 0, 2772.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 11, 1,  -66.905, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 12, 0, 2370.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 12, 1,  699.464, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 13, 0, 1518.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 13, 1, -334.526, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 14, 0, 1962.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 14, 1, -231.127, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 15, 0, 2148.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 15, 1,  650.806, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 16, 0, 2142.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 16, 1,   79.070, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 17, 0, 2316.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 17, 1,  255.457, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 18, 0, 1908.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 18, 1,  535.242, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 19, 0, 2544.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 19, 1,  334.526, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 20, 0, 1896.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 20, 1, -267.621, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 21, 0, 2592.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 21, 1, -255.457, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 22, 0, 2274.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 22, 1,  285.868, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 23, 0, 1482.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 23, 1,  310.197, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 24, 0, 2148.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 24, 1,  -91.234, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 25, 0, 1920.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 25, 1,  -18.247, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 26, 0, 2076.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 26, 1,  218.963, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 27, 0, 1308.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 27, 1, -218.963, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 28, 0, 2124.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 28, 1,  206.798, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 29, 0, 1290.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 29, 1,  200.716, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 30, 0, 1458.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 30, 1,  -36.494, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 31, 0, 2166.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 31, 1,  -42.576, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 32, 0, 2580.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 32, 1, -182.469, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 33, 0, 1536.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 33, 1,   24.329, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 34, 0, 2040.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 34, 1,  -60.823, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 35, 0, 2622.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 35, 1,  -60.823, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 36, 0, 2040.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 36, 1,  145.975, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 37, 0, 2454.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 37, 1,  -85.152, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 38, 0, 1530.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 38, 1,  273.703, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 39, 0, 1770.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 39, 1,  723.794, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 40, 0, 1956.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 40, 1,  -54.741, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 41, 0, 2250.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 41, 1, -474.419, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 42, 0, 1758.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 42, 1,  784.617, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 43, 0, 2004.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 43, 1,  -42.576, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 44, 0, 2700.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 44, 1, -510.913, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 45, 0, 2670.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 45, 1,  310.197, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 46, 0, 2244.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 46, 1,  -42.576, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 47, 0, 1728.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 47, 1, -529.160, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 48, 0, 2184.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 48, 1, -285.868, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 49, 0, 2046.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 49, 1,   -0.000, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 50, 0, 1272.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 50, 1, 1027.909, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 51, 0, 1782.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 51, 1, -218.963, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 52, 0, 1482.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 52, 1,  553.489, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 53, 0, 3342.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 53, 1, -662.971, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 54, 0, 2550.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 54, 1,  170.304, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 55, 0, 2010.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 55, 1, -218.963, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 56, 0, 2106.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 56, 1, -231.127, 1.0, 10.0, 10.0);

  model_prior_signal(mdlA, 57, 0, 1782.000, 1.0, 10.0, 10.0);
  model_prior_signal(mdlA, 57, 1, -845.440, 1.0, 10.0, 10.0);
*/
  model_set_from_priors(mdlA);

  for (unsigned int j = 0; j < mdlA->M; j++) {
    const double aj = 50.0 + 500.0 * ((double) j / (double) mdlA->M);
    vector_set(mdlA->ahat, j, aj);
  }

  /* write the simulated data to a grid. */
  dataset_t *sim = dataset_alloc_from_grid(&grid);
  model_eval(mdlA, NULL, sim);
  dataset_fwrite(sim, "sim.dat");

  /* create a model for learning signal parameters. */
  model_t *mdlB = model_alloc(2, 20);
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
  opt->l0 = 1.0;
  opt->dl = 0.1;
  opt->step_max = 5;
  optim_execute(opt);
  optim_free(opt);

  /* allocate datasets for storing predictions. */
  dataset_t *mean = dataset_alloc(2);
  dataset_t *var = dataset_alloc(2);
  dataset_t *map = dataset_alloc(2);

  /* add times to the gp prediction dataset. */
  for (ti = 0.0; ti <= 1.0; ti += 1.0e-3) {
    vector_set(t, 0, ti);
    for (double tj = 0.0; tj <= 1.0; tj += 1.0e-3) {
      vector_set(t, 1, tj);
      dataset_augment(mean, 0, t, 0.0);
      dataset_augment(mean, 1, t, 0.0);
      dataset_augment(var, 0, t, 0.0);
      dataset_augment(var, 1, t, 0.0);
      dataset_augment(map, 0, t, 0.0);
      dataset_augment(map, 1, t, 0.0);
    }
  }

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

