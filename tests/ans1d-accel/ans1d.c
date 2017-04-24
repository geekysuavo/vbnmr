
/* include the required headers. */
#include <vbnmr/vbnmr.h>
#include <vfl/util/rng.h>
#include <vfl/util/search.h>

/* signals: ground-truth parameters of the measured data.
 */
static const struct {
  double w, omega, mu, tau, alpha, beta;
} signals[] = {
  { 1.0, -0.083, 0.0, 25.0, 50.0, 5.0e4 },
  { 1.0, -0.026, 0.0, 25.0, 50.0, 5.0e4 },
  { 0.1,  0.038, 0.0, 25.0, 50.0, 5.0e4 },
  { 0.3,  0.042, 0.0, 25.0, 50.0, 5.0e4 },
  { 1.0,  0.074, 0.0, 25.0, 50.0, 5.0e4 },
  { 0.9,  0.091, 0.0, 25.0, 50.0, 5.0e4 },
  { 0.0,  0.0,   0.0, 25.0, 50.0, 5.0e4 },
  { 0.0,  0.0,   0.0, 25.0, 50.0, 5.0e4 },
  { 0.0,  0.0,   0.0, 25.0, 50.0, 5.0e4 },
  { 0.0,  0.0,   0.0, 25.0, 50.0, 5.0e4 },
  { 0.0,  0.0,   0.0,  0.0,  0.0, 0.0   } /* end marker. */
};

/* logparms(): write a line to the parameter log file.
 *
 * arguments:
 *  @fh: output file handle to write to.
 *  @opt: optimizer structure pointer.
 */
static void logparms (FILE *fh, const optim_t *opt) {
  /* get sub-structures. */
  const model_t *mdl = opt->mdl;
  const data_t *dat = mdl->dat;

  /* begin the log line. */
  fprintf(fh, "%4u %16.9le", dat->N, opt->bound);

  /* write entries for each factor. */
  for (unsigned int j = 0; j < mdl->M; j++)
    fprintf(fh, " %16.9le %16.9le",
            factor_get(mdl->factors[j], 2),
            factor_get(mdl->factors[j], 3));

  /* end the log line. */
  fprintf(fh, "\n");
  fflush(fh);
}

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

  /* declare variables for writing results to disk. */
  char fmean[256], fvar[256], fdat[256];
  FILE *fpar;

  /* allocate a random number generator. */
  rng_t *R = rng_alloc();

  /* allocate a dataset holding the entire measurement. this will be
   * accessed to augment the measured data during the simulation.
   */
  double grid_values[] = { 0.0, 1.0, 10000.0 };
  matrix_view_t grid = matrix_view_array(grid_values, 1, 3);
  data_t *dsrc = data_alloc_from_grid(2, &grid);

  /* allocate datasets for storing predictions. */
  data_t *mean = data_alloc_from_grid(2, &grid);
  data_t *var = data_alloc_from_grid(2, &grid);

  /* allocate a measured dataset. */
  double ginit_values[] = { 0.0, 10000.0, 10000.0 };
  matrix_view_t ginit = matrix_view_array(ginit_values, 1, 3);
  data_t *dat = data_alloc_from_grid(2, &ginit);

  /* set up a hybrid model. */
  model_t *mdl = model_alloc(vfl_model_tauvfr);
  tauvfr_set_tau(mdl, 1.0e3);
  model_set_nu(mdl, 1.0e-6);
  model_set_data(mdl, dat);

  /* add factors to the model. */
  for (unsigned int j = 0; signals[j].tau; j++) {
    /* create a decay factor. */
    factor_t *fR = factor_alloc(vfl_factor_decay);
    factor_set(fR, 0, signals[j].alpha);
    factor_set(fR, 1, signals[j].beta);
    factor_set_fixed(fR, 1);

    /* create a quadrature factor. */
    factor_t *fV = factor_alloc(vbnmr_factor_quad);
    factor_set(fV, 0, signals[j].mu);
    factor_set(fV, 1, signals[j].tau);
    quad_set_ftsize(fV, 65536);

    /* create a product of the decay and quadrature factors. */
    factor_t *f = factor_alloc(vfl_factor_product);
    product_add_factor(f, 0, fR);
    product_add_factor(f, 0, fV);

    /* add the product factor to the model. */
    model_add_factor(mdl, f);
  }

  /* allocate a search structure for the model. */
  search_t *S = search_alloc(mdl, dat, &grid);
  vector_t *xmax = vector_alloc(1);
  search_set_outputs(S, 2);

  /* set the initial factor frequencies. */
  for (unsigned int j = 0; signals[j].tau; j++) {
    /* set the frequency for filling the source dataset. */
    factor_set(mdl->factors[j], 2, signals[j].omega);
  }

  /* set the initial factor weights. */
  for (unsigned int j = 0, k = 0; signals[j].tau; j++, k += 2) {
    /* set the weights for zero phase shift. */
    vector_set(mdl->wbar, k + 0, signals[j].w);
    vector_set(mdl->wbar, k + 1, 0.0);
  }

  /* evaluate the model over the source dataset. */
  const double sigma = 1.0 / sqrt(mdl->tau);
  for (unsigned int i = 0; i < dsrc->N; i++) {
    /* add a bit of gaussian noise to the evaluation. */
    datum_t *di = data_get(dsrc, i);
    di->y = model_eval(mdl, di->x, di->p) + sigma * rng_normal(R);
  }

  /* evaluate the model over the initial dataset. */
  for (unsigned int i = 0; i < dat->N; i++) {
    /* find the matching source value. */
    for (unsigned int isrc = 0; isrc < dsrc->N; isrc++) {
      /* once found, copy the measurement and break the search. */
      if (dsrc->data[isrc].p == dat->data[i].p &&
          vector_equal(dsrc->data[isrc].x, dat->data[i].x)) {
        dat->data[i].y = dsrc->data[isrc].y;
        break;
      }
    }
  }

  /* write the source dataset. */
  data_fwrite(dsrc, "source.dat");

  /* reset the model. */
  model_reset(mdl);

  /* allocate a mean-field optimizer for the model. */
  optim_t *opt = optim_alloc(vfl_optim_mf);
  optim_set_model(opt, mdl);
  optim_set_max_iters(opt, 1);

  /* open the parameter log file. */
  fpar = fopen("parms.dat", "w");

  /* loop until the maximum number of samples has been reached. */
  while (dat->N <= 200) {
    /* write the current dataset. */
    sprintf(fdat, "meas-%04d.dat", dat->N);
    data_fwrite(dat, fdat);

    /* re-fit the signal model. */
    model_reset(mdl);
    optim_execute(opt);

    /* log the parameters. */
    logparms(fpar, opt);

    /* write the vfr result. */
    if (dat->N % 10 == 0) {
      model_predict_all(mdl, mean, var);
      sprintf(fmean, "vfr-mean-%04d.dat", dat->N);
      sprintf(fvar, "vfr-var-%04d.dat", dat->N);
      data_fwrite(mean, fmean);
      data_fwrite(var, fvar);
    }

    /* execute a search. */
    search_execute(S, xmax);

    /* locate all augmenting points in the source dataset. */
    datum_t *daug = NULL;
    for (unsigned int i = 0; i < dsrc->N; i++) {
      daug = data_get(dsrc, i);
      if (vector_equal(daug->x, xmax))
        data_augment(dat, daug);
    }
  }

  /* close the parameter log file. */
  fclose(fpar);

  /* free the structures. */
  search_free(S);
  vector_free(xmax);
  optim_free(opt);
  model_free(mdl);
  data_free(dsrc);
  data_free(mean);
  data_free(var);
  data_free(dat);
  rng_free(R);

  /* return success. */
  return 0;
}

