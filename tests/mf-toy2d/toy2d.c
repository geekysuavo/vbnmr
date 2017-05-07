
/* include the required headers. */
#include <vbnmr/vbnmr.h>

/* define the input filename. */
#define INFILE "toy2d.dat"

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

  /* read the dataset file. */
  data_t *dat = data_alloc_from_file(INFILE);

  /* set up a hybrid model. */
  model_t *mdl = model_alloc(vbnmr_model_vfgp);
  tauvfr_set_tau(mdl, 1.0);
  model_set_nu(mdl, 1.0e-6);
  model_set_data(mdl, dat);

  /* add factors to the model. */
  const unsigned int M = 3;
  for (unsigned int j = 0; j < M; j++) {
    /* create a first-dimension decay factor. */
    factor_t *fR0 = factor_alloc(vfl_factor_decay);
    factor_set(fR0, 0, 10.0);
    factor_set(fR0, 1, 290000.0);
    factor_set_fixed(fR0, 1);

    /* create a second-dimension decay factor. */
    factor_t *fR1 = factor_alloc(vfl_factor_decay);
    factor_set(fR1, 0, 10.0);
    factor_set(fR1, 1, 15000.0);
    factor_set_fixed(fR1, 1);

    /* create a quadrature factor. */
    factor_t *fV = factor_alloc(vbnmr_factor_quad);
    quad_set_dims(fV, 2);
    factor_set(fV, 0, 0.66867);
    factor_set(fV, 2, 0.0);
    factor_set(fV, 1, 9.0);
    factor_set(fV, 3, 2.0);
    quad_set_ftsize(fV, 65536);

    /* create a product of the decays and quadrature factors. */
    factor_t *f = factor_alloc(vfl_factor_product);
    product_add_factor(f, 0, fR0);
    product_add_factor(f, 0, fR1);
    product_add_factor(f, 0, fV);

    /* add the product factor to the model. */
    model_add_factor(mdl, f);
  }

  /* optimize, using mean-field. */
  optim_t *opt = optim_alloc(vfl_optim_mf);
  optim_set_model(opt, mdl);
  optim_set_max_iters(opt, 1);
  optim_execute(opt);

  /* allocate datasets for fit analysis. */
  data_t *fit_mean = data_alloc_from_file(INFILE);
  data_t *fit_var = data_alloc_from_file(INFILE);

  /* compute and output the vfr prediction. */
  model_predict_all(mdl, fit_mean, fit_var);
  data_fwrite(fit_mean, "fit-mean.dat");
  data_fwrite(fit_var, "fit-var.dat");

  /* allocate datasets for prediction. */
  double grid_values[] = { 0.0, 0.2,  60.0,
                           0.0, 1.0, 140.0 };
  matrix_view_t grid = matrix_view_array(grid_values, 2, 3);
  data_t *mean = data_alloc_from_grid(2, &grid);
  data_t *var = data_alloc_from_grid(2, &grid);

  /* compute and output the vfr prediction. */
  model_predict_all(mdl, mean, var);
  data_fwrite(mean, "mean.dat");
  data_fwrite(var, "var.dat");

  /* free the structures. */
  obj_release((object_t*) fit_mean);
  obj_release((object_t*) fit_var);
  obj_release((object_t*) mean);
  obj_release((object_t*) var);
  obj_release((object_t*) opt);

  /* return success. */
  return 0;
}

