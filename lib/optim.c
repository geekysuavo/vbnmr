
/* include the optimization header. */
#include <vbnmr/optim.h>

/* include the distribution headers. */
#include <vbnmr/dist-normal.h>
#include <vbnmr/dist-gamma.h>
#include <vbnmr/dist-basis.h>

/* include the private functions. */
#include "optim-priv.c"

/* * * * core optimizer initialization/preparation functions: * * * */

/* initialize(): initialize the pointers within an optimization
 * structure to null.
 *
 * arguments:
 *  @opt: pointer to the optimization structure to initialize.
 */
static inline void initialize (optim_t *opt) {
  /* initialize the learning parameters. */
  opt->elbo = -INFINITY;
  opt->L = NULL;
  opt->G = NULL;
  opt->h = NULL;

  /* initialize the log file. */
  opt->log = NULL;

  /* initialize the temporary vector. */
  opt->tmp = NULL;
}

/* prepare(): prepare the structure members of an optimizer
 * for iteration based on the dataset and model sizes.
 *
 * arguments:
 *  @opt: pointer to the optimization structure to prepare.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
static inline int prepare (optim_t *opt) {
  /* fail if any required structure pointer is null. */
  if (!opt || !opt->dat || !opt->mdl)
    return 0;

  /* get the current signal, model and data sizes. */
  const unsigned int D = opt->mdl->D;
  const unsigned int K = opt->mdl->K;
  const unsigned int M = opt->mdl->M;
  const unsigned int N = opt->dat->N;

  /* compute the size of the weight space. */
  const unsigned int KM = K * M;

  /* compute the temporary vector length:
   *    D      : time vector.
   *  4 D      : novel x-iterates.
   *  4 D      : beta=0 x-iterate values.
   *  4 D      : beta=inf x-iterate values.
   *    D      : decay factors.
   *    D      : modulation factors.
   *  2 M K^2  : cholesky updates and downdates.
   *    M K    : covariance updates and downdates.
   */
  const unsigned int ntmp1 = 15 * D + 2 * KM * K + KM;
  const unsigned int ntmp2 = D + 2 * K * (K + 1) + 6 * N;
  const unsigned int ntmp = (ntmp1 > ntmp2 ? ntmp1 : ntmp2);

  /* return if the structure members are already prepared. */
  if (opt->K == K && opt->M == M && opt->N == N)
    return 1;

  /* free the structure members. */
  matrix_free(opt->L);
  matrix_free(opt->G);
  vector_free(opt->h);
  vector_free(opt->tmp);

  /* allocate the structure members. */
  opt->L = matrix_alloc(KM, KM);
  opt->G = matrix_alloc(KM, KM);
  opt->h = vector_alloc(KM);
  opt->tmp = vector_alloc(ntmp);

  /* check if any of the allocations failed. */
  if (!opt->G || !opt->L || !opt->h || !opt->tmp)
    return 0;

  /* initialize the contents of the structure members. */
  matrix_set_zero(opt->L);
  matrix_set_zero(opt->G);
  vector_set_zero(opt->h);

  /* store the new sizes. */
  opt->D = D;
  opt->K = K;
  opt->M = M;
  opt->N = N;

  /* return success. */
  return 1;
}

/* * * * public functions: * * * */

/* optim_alloc(): allocate a new optimizer for use.
 *
 * arguments:
 *  @mdl: pointer to a model to optimize.
 *  @dat: pointer to a dataset to learn from.
 *
 * returns:
 *  newly allocated and initialized optimizer structure pointer.
 */
optim_t *optim_alloc (model_t *mdl, dataset_t *dat) {
  /* @opt: pointer to the new optimizer structure.
   */
  optim_t *opt;

  /* fail if either passed pointer is null. */
  if (!mdl || !dat)
    return NULL;

  /* allocate a new structure pointer, or fail. */
  opt = (optim_t*) malloc(sizeof(optim_t));
  if (!opt)
    return NULL;

  /* initialize the sizes of the signal, data, and model. */
  opt->D = opt->K = opt->M = opt->N = 0;

  /* initialize the dataset and model pointers. */
  initialize(opt);
  opt->dat = dat;
  opt->mdl = mdl;

  /* initialize the remaining structure members. */
  if (!prepare(opt)) {
    /* free any allocated memory and return null. */
    optim_free(opt);
    return NULL;
  }

  /* initialize the parameters of the optimizer. */
  opt->iter_max = 1000;
  opt->step_max = 10;
  opt->l0 = 1.0;
  opt->dl = 0.1;

  /* return the new structure pointer. */
  return opt;
}

/* optim_free(): free an allocated optimizer.
 *
 * arguments:
 *  @opt: optimizer structure pointer to free.
 */
void optim_free (optim_t *opt) {
  /* return if the structure pointer is null. */
  if (!opt) return;

  /* free the structure members. */
  matrix_free(opt->L);
  matrix_free(opt->G);
  vector_free(opt->h);
  vector_free(opt->tmp);

  /* close any open logs. */
  if (opt->log)
    fclose(opt->log);

  /* free the structure pointer. */
  free(opt);
}

/* optim_init(): perform an initial mean-field iteration in order to
 * set the parameters of a signal model.
 *
 * arguments:
 *  @opt: pointer to the optimizer to utilize.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int optim_init (optim_t *opt) {
  /* fail if the structure pointer is null. */
  if (!opt) return 0;

  /* prepare for iteration. */
  if (!prepare(opt))
    return 0;

  /* initialize the linear parameters. */
  update_linear(opt);

  /* loop over each signal in the model. */
  for (unsigned int j = 0; j < opt->M; j++) {
    /* loop over each dimension of the signal parameters. */
    for (unsigned int d = 0; d < opt->D; d++) {
      /* perform an approximate mean-field update. */
      update_meanfield(opt, j, d);
    }

    /* update the amplitude parameters. */
    update_linear_single(opt, j);
  }

  /* compute the new bound. */
  opt->elbo = evaluate_bound(opt);

  /* return success. */
  return 1;
}

/* optim_execute(): perform a free-running optimization of a set of
 * model parameters using collapsed variational bayesian inference.
 *
 * arguments:
 *  @opt: pointer to the optimizer to utilize.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int optim_execute (optim_t *opt) {
  /* declare required variables:
   *  @elbo0: initial value of the lower bound.
   *  @iter: iteration counter.
   *  @done: completion flag.
   */
  unsigned int iter, done;
  double elbo0;

  /* fail if the structure pointer is null. */
  if (!opt) return 0;

  /* prepare for iteration. */
  if (!prepare(opt))
    return 0;

  /* initialize the linear parameters and compute the bound. */
  update_linear(opt);
  opt->elbo = evaluate_bound(opt);

  /* loop until convergence or resource exhaustion. */
  iter = done = 0;
  while (!done) {
    /* increment the iteration counter. */
    iter++;

    /* store the previous value of the lower bound. */
    elbo0 = opt->elbo;

    /* update each of the signals in the model. */
    for (unsigned int j = 0; j < opt->M; j++)
      update_nonlinear(opt, j);

    /* check if a log file is open. */
    if (opt->log) {
      /* begin a new line. */
      fprintf(opt->log, "%u %16.9le", iter, opt->elbo);

      /* write the parameters of each signal. */
      for (unsigned int j = 0; j < opt->M; j++)
        for (unsigned int d = 0; d < opt->D; d++)
          fprintf(opt->log, " %16.9le %16.9le %16.9le %16.9le",
                  opt->mdl->sig[j][d].mu,
                  opt->mdl->sig[j][d].tau,
                  opt->mdl->sig[j][d].alpha,
                  opt->mdl->sig[j][d].beta);

      /* end the current line. */
      fprintf(opt->log, "\n");
      fflush(opt->log);
    }

    /* terminate if the bound was unchanged or
     * if the resource limit was exceeded.
     */
    if (opt->elbo == elbo0 || iter >= opt->iter_max)
      done = 1;
  }

  /* return success. */
  return 1;
}

