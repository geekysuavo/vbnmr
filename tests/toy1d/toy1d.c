
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

  /* read the dataset file. */
  data_t *dat = data_alloc_from_file("toy1d.dat");

  /* set up a hybrid model. */
  model_t *mdl = model_alloc(vbnmr_model_vfgp);
  model_set_alpha0(mdl, 1.0);
  model_set_nu(mdl, 1.0e-6);
  model_set_data(mdl, dat);

  /* add factors to the model. */
  const unsigned int M = 3;
  for (unsigned int j = 0; j < M; j++) {
    /* create a decay factor. */
    factor_t *fR = factor_alloc(factor_type_decay);
    factor_set(fR, 0, 100.0);
    factor_set(fR, 1, 100.0);
    factor_set_fixed(fR, 1);

    /* create a quadrature factor. */
    factor_t *fV = factor_alloc(vbnmr_factor_quad);
    factor_set(fV, 0, 0.0);
    factor_set(fV, 1, 1.0e-5);

    /* create a product of the decay and quadrature factors. */
    factor_t *f = factor_alloc(factor_type_product);
    product_add_factor(f, 0, fR);
    product_add_factor(f, 0, fV);

    /* add the product factor to the model. */
    model_add_factor(mdl, f);
  }

  /* randomly initialize the factor frequency means. */
  for (unsigned int j = 0; j < mdl->M; j++)
    factor_set(mdl->factors[j], 2, 300.0 * rng_normal(R));

  /* optimize. */
  optim_t *opt = optim_alloc(optim_type_fg);
  optim_set_model(opt, mdl);
  optim_set_lipschitz_init(opt, 0.001);
  optim_execute(opt);

  /* allocate datasets for prediction. */
  double grid_values[] = { 0.0, 2.0e-4, 0.5 };
  matrix_view_t grid = matrix_view_array(grid_values, 1, 3);
  data_t *mean = data_alloc_from_grid(2, &grid);
  data_t *var = data_alloc_from_grid(2, &grid);

  /* compute and output the vfr prediction. */
  model_predict_all(mdl, mean, var);
  data_fwrite(mean, "vfr-mean.dat");
  data_fwrite(var, "vfr-var.dat");

  /* compute and output the gp prediction. */
  vfgp_set_mode(mdl, 1);
  model_predict_all(mdl, mean, var);
  data_fwrite(mean, "gp-mean.dat");
  data_fwrite(var, "gp-var.dat");

  /* free the structures. */
  optim_free(opt);
  model_free(mdl);
  data_free(mean);
  data_free(var);
  data_free(dat);
  rng_free(R);

  /* return success. */
  return 0;
}

