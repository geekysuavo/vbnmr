
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
  double ti;

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

  /* optimize the model parameters. */
  optim_t *opt = optim_alloc(mdl, dat);
  opt->l0 = 0.01;
  opt->dl = 0.1;
  optim_execute(opt);
  optim_free(opt);

  /* allocate datasets for storing predictions. */
  dataset_t *mean = dataset_alloc(1);
  dataset_t *var = dataset_alloc(1);
  dataset_t *map = dataset_alloc(1);

  /* add times to the prediction dataset. */
  for (ti = 0.0; ti <= 1.0; ti += 1.0e-4) {
    vector_set(t, 0, ti);
    dataset_augment(mean, 0, t, 0.0);
    dataset_augment(mean, 1, t, 0.0);
    dataset_augment(var, 0, t, 0.0);
    dataset_augment(var, 1, t, 0.0);
    dataset_augment(map, 0, t, 0.0);
    dataset_augment(map, 1, t, 0.0);
  }

  /* compute the posterior predictions. */
  model_eval(mdl, NULL, map);
  model_predict(mdl, mean, var);
/*FIXME*/
const unsigned int M = mdl->M;
for (unsigned int j = 0; j < M; j++) {
  const double omega1 = mdl->sig[j][0].mu;
  const double omega2 = 1.0 / sqrt(mdl->sig[j][0].tau);
  const double rho1 = mdl->sig[j][0].alpha / mdl->sig[j][0].beta;
  const double rho2 = mdl->sig[j][0].alpha /
                      pow(mdl->sig[j][0].beta, 2.0);
  const double a1 = sqrt(pow(vector_get(mdl->ahat, j), 2.0) +
                         pow(vector_get(mdl->ahat, j + M), 2.0));
  const double a2 = sqrt(matrix_get(mdl->Sigma, j, j) +
                         matrix_get(mdl->Sigma, j + M, j + M));
   fprintf(stderr, "    $%.1lf \\pm %.3lf$ & "
                      "$%.1lf \\pm %.3lf$ & "
                      "$%.3lf \\pm %.6lf$ \\\\\n",
    omega1, omega2, rho1, rho2, a1, a2);
}
/*FIXME*/

  /* write the prediction datasets to files. */
  dataset_fwrite(mean, "mean.dat");
  dataset_fwrite(var, "var.dat");
  dataset_fwrite(map, "map.dat");

  /* free all allocated structures. */
  dataset_free(mean);
  dataset_free(var);
  dataset_free(map);
  dataset_free(dat);
  model_free(mdl);
  vector_free(t);
  rng_free(R);

  /* return success. */
  return 0;
}

