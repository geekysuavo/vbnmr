
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

  /* create a model for generating simulated measurements. */
  model_t *mdlA = model_alloc(1, 3);
  model_set_sigma(mdlA, 1.0);
  model_set_delta(mdlA, 1.0e3);
  model_prior_signal(mdlA, 0, 0, 200.0, 1.0, 100.0, 100.0);
  model_prior_signal(mdlA, 1, 0, -33.0, 1.0, 100.0, 100.0);
  model_prior_signal(mdlA, 2, 0,   7.0, 1.0, 100.0, 100.0);
  model_set_from_priors(mdlA);
  vector_set(mdlA->ahat, 0, 15.000); /*  30 * cos(pi/3) */
  vector_set(mdlA->ahat, 1, 12.500); /*  25 * cos(pi/3) */
  vector_set(mdlA->ahat, 2, 50.000); /* 100 * cos(pi/3) */
  vector_set(mdlA->ahat, 3, 25.981); /*  30 * sin(pi/3) */
  vector_set(mdlA->ahat, 4, 21.651); /*  25 * sin(pi/3) */
  vector_set(mdlA->ahat, 5, 86.603); /* 100 * sin(pi/3) */

  /* create a model for learning signal parameters. */
  model_t *mdlB = model_alloc(1, 3);
  model_set_sigma(mdlB, mdlA->sigma);
  model_set_delta(mdlB, mdlA->delta);
  model_prior_signal(mdlB, 0, 0, 0.0, 1.0e-5, 100.0, 100.0);
  model_prior_signal(mdlB, 1, 0, 0.0, 1.0e-5, 100.0, 100.0);
  model_prior_signal(mdlB, 2, 0, 0.0, 1.0e-5, 100.0, 100.0);
  model_set_from_priors(mdlB);
  model_randomize_mu(mdlB, R);

  /* allocate a dataset for storing simulated measurements. */
  dataset_t *dat = dataset_alloc(1);

  /* add real measurement times. */
  for (ti = 1.0e-4; ti <= 0.2; ti += 1.0e-3) {
    vector_set(t, 0, ti);
    dataset_augment(dat, 0, t, 0.0);
  }

  /* add imaginary measurement times. */
  for (ti = 3.0e-4; ti <= 0.2; ti += 1.3e-3) {
    vector_set(t, 0, ti);
    dataset_augment(dat, 1, t, 0.0);
  }

  /* evaluate the model on the dataset, with gaussian noise. */
  model_eval(mdlA, R, dat);

  /* write the input dataset to a file. */
  dataset_fwrite(dat, "meas.dat");

  /* optimize the model parameters. */
  optim_t *opt = optim_alloc(mdlB, dat);
  opt->l0 = 0.1;
  opt->dl = 0.1;
  optim_execute(opt);
  optim_free(opt);

  /* allocate datasets for storing predictions. */
  dataset_t *mean = dataset_alloc(1);
  dataset_t *var = dataset_alloc(1);
  dataset_t *map = dataset_alloc(1);

  /* add times to the gp prediction dataset. */
  for (ti = 0.0; ti <= 0.5; ti += 1.0e-4) {
    vector_set(t, 0, ti);
    dataset_augment(mean, 0, t, 0.0);
    dataset_augment(mean, 1, t, 0.0);
    dataset_augment(var, 0, t, 0.0);
    dataset_augment(var, 1, t, 0.0);
    dataset_augment(map, 0, t, 0.0);
    dataset_augment(map, 1, t, 0.0);
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

