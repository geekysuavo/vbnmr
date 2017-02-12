
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

  /* declare strings for holding filenames, and a file handle. */
  char fmean[256], fvar[256], fdat[256], fnew[256];
  FILE *fh, *fpar;

  /* create a random number generator and a time vector. */
  rng_t *R = rng_alloc();
  vector_t *t = vector_alloc(1);

  /* create a model for simulation. */
  model_t *truth = model_alloc(1, 6);
  model_set_sigma(truth, 0.01);
  model_set_delta(truth, 100.0);
  model_prior_signal(truth, 0, 0, -830.0, 100.0, 1000.0, 100.0);
  model_prior_signal(truth, 1, 0, -260.0, 100.0, 1000.0, 100.0);
  model_prior_signal(truth, 2, 0,  380.0, 100.0, 1000.0, 100.0);
  model_prior_signal(truth, 3, 0,  420.0, 100.0, 1000.0, 100.0);
  model_prior_signal(truth, 4, 0,  740.0, 100.0, 1000.0, 100.0);
  model_prior_signal(truth, 5, 0,  910.0, 100.0, 1000.0, 100.0);
  model_set_from_priors(truth);
  vector_set_zero(truth->ahat);
  vector_set(truth->ahat, 0, 1.0);
  vector_set(truth->ahat, 1, 1.0);
  vector_set(truth->ahat, 2, 0.1);
  vector_set(truth->ahat, 3, 0.3);
  vector_set(truth->ahat, 4, 1.0);
  vector_set(truth->ahat, 5, 0.9);

  /* create a model for learning. */
  model_t *fit = model_alloc(1, 10);
  model_set_sigma(fit, truth->sigma);
  model_set_delta(fit, truth->delta);
  for (unsigned int j = 0; j < fit->M; j++)
    model_prior_signal(fit, j, 0, 0.0, 1.0e-6, 100.0, 10.0);
  model_set_from_priors(fit);

  /* create an initial dataset. */
  double ginit_values[] = { 0.000001, 1.0, 1.000001 };
  matrix_view_t ginit = matrix_view_array(ginit_values, 1, 3);
  dataset_t *dat = dataset_alloc_from_grid(&ginit);
  model_eval(truth, R, dat);

  /* set up a gp over the fit/data pair. */
  gp_t *gp = gp_alloc(dat, fit);

  /* set up an optimizer over the fit/data pair. */
  optim_t *opt = optim_alloc(fit, dat);

  /* create a dataset for storing the ground truth. */
  double grid_values[] = { 1.0e-6, 1.0e-4, 1.0 };
  matrix_view_t grid = matrix_view_array(grid_values, 1, 3);
  dataset_t *Dsched = dataset_alloc_from_grid(&grid);
  dataset_t *Dtrue = dataset_alloc_from_grid(&grid);
  model_eval(truth, NULL, Dtrue);
  dataset_fwrite(Dtrue, "truth.dat");

  /* create datasets for storing gp predictions. */
  dataset_t *Dmean = dataset_alloc_from_grid(&grid);
  dataset_t *Dvar = dataset_alloc_from_grid(&grid);

  /* open the parameter log file. */
  fpar = fopen("parms.dat", "w");

  /* loop until the maximum number of samples has been reached. */
  while (dat->N <= 100) {
    /* write the current dataset. */
    sprintf(fdat,  "meas-%03u.dat", dat->N);
    dataset_fwrite(dat, fdat);

    /* refit the underlying signal model. */
    model_set_from_priors(fit);
    optim_init(opt);

    /* log the new parameters. */
    fprintf(fpar, "%4u %16.9le", dat->N, opt->elbo);
    for (unsigned int j = 0; j < fit->M; j++)
      fprintf(fpar, " %16.9le %16.9le",
              fit->sig[j][0].mu,
              fit->sig[j][0].tau);
    fprintf(fpar, "\n");
    fflush(fpar);

    /* refit the gp. */
    gp_fit(gp);

    /* write the vfl result. */
    model_predict(fit, Dmean, Dvar);
    sprintf(fmean, "vfr-mean-%03u.dat", dat->N);
    sprintf(fvar,  "vfr-var-%03u.dat",  dat->N);
    dataset_fwrite(Dmean, fmean);
    dataset_fwrite(Dvar, fvar);

    /* write the gp result. */
    gp_eval(gp, Dmean, Dvar);
    sprintf(fmean, "gp-mean-%03u.dat", dat->N);
    sprintf(fvar,  "gp-var-%03u.dat",  dat->N);
    dataset_fwrite(Dmean, fmean);
    dataset_fwrite(Dvar, fvar);

    /* initialize the search results. */
    double tmax = 0.0, v = 0.0, vmax = 0.0;
    unsigned int imax = 0, kmax = 0;

    /* locate the time point that has maximum variance. */
    for (unsigned int k = 0; k < Dvar->K; k++) {
      for (unsigned int i = 0; i < Dvar->n[k]; i++) {
        /* check if the current time is higher. */
        dataset_get(Dvar, k, i, t, &v);
        if (v > vmax && vector_get(Dsched->y[k], i) == 0.0) {
          vmax = v;
          kmax = k;
          imax = i;
          tmax = vector_get(t, 0);
        }
      }
    }

    /* augment the dataset. */
    vector_set(t, 0, tmax);
    double y = model_eval_single(truth, kmax, t);
    dataset_augment(dat, kmax, t, y);
    vector_set(Dsched->y[kmax], imax, 1.0);

    /* output the augmenting point. */
    sprintf(fnew,  "new-%03u.dat",  dat->N);
    fh = fopen(fnew, "w");
    fprintf(fh, "%u %le %le\n", kmax, tmax, y);
    fclose(fh);
  }

  /* close the parameter log file. */
  fclose(fpar);

  /* free the gp and optimizer. */
  gp_free(gp);
  optim_free(opt);

  /* free the datasets. */
  dataset_free(dat);
  dataset_free(Dsched);
  dataset_free(Dtrue);
  dataset_free(Dmean);
  dataset_free(Dvar);

  /* free the models. */
  model_free(truth);
  model_free(fit);

  /* free the vector and random number generator. */
  vector_free(t);
  rng_free(R);

  /* return success. */
  return 0;
}

