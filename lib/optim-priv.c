
/* evaluate_bound(): compute a lower bound on the log-evidence of a dataset
 * under a particular signal model.
 *
 * arguments:
 *  @opt: pointer to the optimization structure to access.
 *
 * returns:
 *  value of the lower bound.
 */
double evaluate_bound (const optim_t *opt) {
  /* declare required variables:
   *  @elbo: computation result.
   */
  double elbo = 0.0;

  /* gain references to commonly accessed variables. */
  const model_t *mdl = opt->mdl;
  const unsigned int D = opt->D;
  const unsigned int K = opt->K;
  const unsigned int M = opt->M;
  const unsigned int KM = K * M;

  /* set up a temporary vector view to store the intermediate result
   * of the quadratic term computation. the start index is set to
   * make sure bound computations don't screw up line searches.
   */
  vector_view_t b = vector_subvector(opt->tmp, 15 * D, KM);

  /* include the quadratic term. */
  blas_dtrmv(BLAS_TRANS, 1.0, opt->L, mdl->ahat, 0.0, &b);
  elbo = 0.5 * blas_ddot(&b, &b);

  /* include the log-determinant term. */
  for (unsigned int j = 0; j < KM; j++)
    elbo -= log(matrix_get(opt->L, j, j));

  /* loop over each model dimension. */
  for (unsigned int d = 0; d < D; d++) {
    /* loop over each model signal. */
    for (unsigned int j = 0; j < M; j++) {
      /* include the frequency and decay entropy terms. */
      elbo += nrm_entropy(mdl, j, d);
      elbo += gam_entropy(mdl, j, d);

      /* include the frequency and decay prior terms. */
      elbo += nrm_prior(mdl, j, d);
      elbo += gam_prior(mdl, j, d);
    }
  }

  /* return the computed lower bound. */
  return elbo;
}

/* update_linear(): update the linear parameters of a signal model
 * based on a dataset and the current nonlinear model parameters.
 *
 * arguments:
 *  @opt: pointer to the optimizer structure to utilize.
 */
static void
update_linear (optim_t *opt) {
  /* declare a variable to hold measurements. */
  double y;

  /* gain access to required structure members. */
  model_t *mdl = opt->mdl;
  dataset_t *dat = opt->dat;
  const unsigned int D = mdl->D;
  const unsigned int K = mdl->K;
  const unsigned int M = mdl->M;
  const unsigned int *n = dat->n;

  /* compute the noise and amplitude precisions. */
  const double invs2 = 1.0 / (mdl->sigma * mdl->sigma);
  const double invd2 = 1.0 / (mdl->delta * mdl->delta);

  /* initialize the interactions and projections. */
  matrix_set_zero(opt->G);
  vector_set_zero(opt->h);

  /* declare views for looking at measurement times,
   * interactions and projections of individual terms.
   */
  matrix_view_t G;
  vector_view_t h;
  vector_view_t t;

  /* create a vector view for storing measurement times. */
  t = vector_subvector(opt->tmp, 0, D);

  /* loop over the basis terms. */
  for (unsigned int l1 = 0; l1 < K; l1++) {
    /* gain a view into the projection vector. */
    h = vector_subvector(opt->h, l1 * M, M);

    /* update each element of the projection vector. */
    for (unsigned int j1 = 0; j1 < M; j1++) {
      /* initialize the vector element. */
      double hj = 0.0;

      /* loop over all time,measurement pairs. */
      for (unsigned int k = 0; k < K; k++) {
        for (unsigned int i = 0; i < n[k]; i++) {
          dataset_get(dat, k, i, &t, &y);
          hj += y * expect_basis(mdl, j1, k, l1, &t);
        }
      }

      /* store the vector element. */
      vector_set(&h, j1, hj);
    }

    /* loop again over the basis terms. */
    for (unsigned int l2 = 0; l2 < K; l2++) {
      /* gain a view into the interaction matrix. */
      G = matrix_submatrix(opt->G, l1 * M, l2 * M, M, M);

      /* update each element of the interaction matrix. */
      for (unsigned int j1 = 0; j1 < M; j1++) {
        for (unsigned int j2 = 0; j2 < M; j2++) {
          /* initialize the matrix element. */
          double gjj = 0.0;

          /* loop over all measurement times. */
          for (unsigned int k = 0; k < K; k++) {
            for (unsigned int i = 0; i < n[k]; i++) {
              dataset_get(dat, k, i, &t, &y);
              gjj += interact_basis(mdl, j1, j2, l1, l2, k, &t);
            }
          }

          /* store the matrix element. */
          matrix_set(&G, j1, j2, gjj);
        }
      }
    }
  }

  /* copy out and complete the current interaction matrix. */
  matrix_copy(opt->L, opt->G);
  matrix_scale(opt->L, invs2);
  vector_view_t Ldiag = matrix_diag(opt->L);
  vector_add_const(&Ldiag, invd2);

  /* compute the cholesky decomposition of the precision matrix. */
  cholesky_decomp(opt->L);

  /* compute the amplitude means. */
  cholesky_solve(opt->L, opt->h, mdl->ahat);
  vector_scale(mdl->ahat, invs2);

  /* compute the amplitude covariances. */
  cholesky_invert(opt->L, mdl->Sigma);
}

/* update_linear_single(): update the linear parameters of a signal model
 * based on a dataset and the current nonlinear model parameters,
 * where it is assumed that only one signal has been modified.
 *
 * this function uses low-rank updates to the cholesky decomposition in
 * order to more efficiently recompute quantities during optimization.
 *
 * arguments:
 *  @opt: pointer to the optimization structure to utilize.
 *  @j: model basis index to update parameters for.
 */
static void
update_linear_single (optim_t *opt, const unsigned int j) {
  /* declare a variable to hold measurements. */
  double y;

  /* gain access to required structure members. */
  model_t *mdl = opt->mdl;
  dataset_t *dat = opt->dat;
  const unsigned int D = mdl->D;
  const unsigned int K = mdl->K;
  const unsigned int M = mdl->M;
  const unsigned int KM = K * M;
  const unsigned int *n = dat->n;

  /* compute the noise precision. */
  const double invs2 = 1.0 / (mdl->sigma * mdl->sigma);

  /* declare views for looking at measurement times,
   * interactions and projections of individual terms.
   */
  matrix_view_t G;
  vector_view_t h;
  vector_view_t t;

  /* declare views for looking at individual rows of the
   * update and downdate matrices.
   */
  vector_view_t u, v, z;
  matrix_view_t U, V;

  /* create the vector view for storing measurement times. */
  t = vector_subvector(opt->tmp, 0, D);

  /* create the matrix view for symmetric rank-1 updates. */
  double *ptr = opt->tmp->data + (15 * D);
  U = matrix_view_array(ptr, K, KM);
  ptr += K * KM;

  /* create the matrix view for symmetric rank-1 downdates. */
  V = matrix_view_array(ptr, K, KM);
  ptr += K * KM;

  /* create the vector view for inverse matrix updates and downdates. */
  z = vector_view_array(ptr, KM);

  /* copy the initial row values of the interaction matrix. */
  for (unsigned int k = 0, idx = j; k < K; k++, idx += M) {
    /* copy the needed row of the interaction matrix. */
    u = matrix_row(&U, k);
    matrix_copy_row(&u, opt->G, idx);
  }

  /* loop over the basis terms. */
  for (unsigned int l1 = 0; l1 < K; l1++) {
    /* gain a view into the projection vector. */
    h = vector_subvector(opt->h, l1 * M, M);

    /* compute the projection vector element. */
    double hj = 0.0;
    for (unsigned int k = 0; k < K; k++) {
      for (unsigned int i = 0; i < n[k]; i++) {
        dataset_get(dat, k, i, &t, &y);
        hj += y * expect_basis(mdl, j, k, l1, &t);
      }
    }

    /* store the vector element. */
    vector_set(&h, j, hj);

    /* loop again over the basis terms. */
    for (unsigned int l2 = 0; l2 < K; l2++) {
      /* gain a view into the interaction matrix. */
      G = matrix_submatrix(opt->G, l1 * M, l2 * M, M, M);

      /* update each element of the interaction matrix. */
      for (unsigned int j2 = 0; j2 < M; j2++) {
        /* compute the interaction matrix elements. */
        double gj1 = 0.0, gj2 = 0.0;
        for (unsigned int k = 0; k < K; k++) {
          for (unsigned int i = 0; i < n[k]; i++) {
            dataset_get(dat, k, i, &t, &y);
            gj1 += interact_basis(mdl, j, j2, l1, l2, k, &t);
            if (j != j2)
              gj2 += interact_basis(mdl, j2, j, l1, l2, k, &t);
          }
        }

        /* store the matrix elements. */
        matrix_set(&G, j, j2, gj1);
        if (j != j2)
          matrix_set(&G, j2, j, gj2);
      }
    }
  }

  /* copy the final row values of the interaction matrix. */
  for (unsigned int k = 0, idx = j; k < K; k++, idx += M) {
    /* copy the needed row of the interaction matrix. */
    v = matrix_row(&V, k);
    matrix_copy_row(&v, opt->G, idx);
  }

  /* compute the difference between the precision rows. */
  matrix_sub(&V, &U);
  matrix_scale(&V, invs2);

  /* adjust the row differences for use in rank-1 updates. */
  for (unsigned int k = 0, idx = j; k < K; k++, idx += M) {
    /* get a view of the current row. */
    v = matrix_row(&V, k);

    /* scale the main-diagonal element by one-half, and zero all
     * off-diagonals that have already been updated.
     */
    vector_set(&v, idx, 0.5 * vector_get(&v, idx));
    for (unsigned int k2 = 0, idx2 = j; k2 < k; k2++, idx2 += M)
      vector_set(&v, idx2, 0.0);
  }

  /* transform the precision row differences into symmetric updates. */
  for (unsigned int k = 0, idx = j; k < K; k++, idx += M) {
    /* get views of the update and downdate vectors. */
    u = matrix_row(&U, k);
    v = matrix_row(&V, k);

    /* compute the symmetrization constants. */
    const double vnrm = blas_dnrm2(&v);
    const double alpha = sqrt(vnrm / 2.0);
    const double beta = 1.0 / vnrm;

    /* symmetrize the vectors. */
    for (unsigned int i = 0; i < KM; i++) {
      /* get the elements of the selector and asymmetric update. */
      const double ui = (i == idx ? 1.0 : 0.0);
      const double vi = vector_get(&v, i);

      /* compute the elements of the symmetric update/downdate. */
      const double xi = alpha * (ui + beta * vi);
      const double yi = alpha * (ui - beta * vi);

      /* store the elements back into their vectors. */
      vector_set(&u, i, xi);
      vector_set(&v, i, yi);
    }
  }

  /* apply the updates. */
  for (unsigned int k = 0; k < K; k++) {
    /* update the cholesky factors. */
    u = matrix_row(&U, k);
    matrix_copy_row(&z, &U, k);
    cholesky_update(opt->L, &z);

    /* update the covariance matrix. */
    blas_dgemv(BLAS_NO_TRANS, 1.0, mdl->Sigma, &u, 0.0, &z);
    double zudot = blas_ddot(&z, &u);
    zudot = 1.0 / (1.0 + zudot);
    for (unsigned int i = 0; i < KM; i++)
      for (unsigned int j = 0; j < KM; j++)
        matrix_set(mdl->Sigma, i, j,
          matrix_get(mdl->Sigma, i, j) -
          zudot * vector_get(&z, i) *
                  vector_get(&z, j));
  }

  /* apply the downdates. */
  for (unsigned int k = 0; k < K; k++) {
    /* downdate the cholesky factors. */
    v = matrix_row(&V, k);
    matrix_copy_row(&z, &V, k);
    cholesky_downdate(opt->L, &z);

    /* downdate the covariance matrix. */
    blas_dgemv(BLAS_NO_TRANS, 1.0, mdl->Sigma, &v, 0.0, &z);
    double zvdot = blas_ddot(&z, &v);
    zvdot = 1.0 / (1.0 - zvdot);
    for (unsigned int i = 0; i < KM; i++)
      for (unsigned int j = 0; j < KM; j++)
        matrix_set(mdl->Sigma, i, j,
          matrix_get(mdl->Sigma, i, j) +
          zvdot * vector_get(&z, i) *
                  vector_get(&z, j));
  }

  /* compute the amplitude means without inverting the precision matrix. */
  cholesky_solve(opt->L, opt->h, mdl->ahat);
  vector_scale(mdl->ahat, invs2);
}

/* include the natural gradient and mean-field private functions. */
#include "optim-priv-ng.c"
#include "optim-priv-mf.c"

