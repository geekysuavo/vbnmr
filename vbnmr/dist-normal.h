
/* ensure once-only inclusion. */
#ifndef __VBNMR_DIST_NORMAL_H__
#define __VBNMR_DIST_NORMAL_H__

/* include the model header. */
#include <vbnmr/model.h>

/* function declarations (dist-normal.c): */

double nrm_entropy (const model_t *mdl,
                    const unsigned int j,
                    const unsigned int d);

double nrm_prior (const model_t *mdl,
                  const unsigned int j,
                  const unsigned int d);

double expect_freq (const model_t *mdl, const unsigned int j,
                    const unsigned int k, const unsigned int l,
                    const unsigned int delta, const vector_t *t);

double expect_all_freq (const model_t *mdl, const unsigned int j,
                        const unsigned int k, const unsigned int l,
                        const vector_t *t, vector_t *V);

double interact_freq (const model_t *mdl,
                      const unsigned int j1, const unsigned int l1,
                      const unsigned int j2, const unsigned int l2,
                      const unsigned int k, const unsigned int delta,
                      const vector_t *t);

double interact_all_freq (const model_t *mdl,
                          const unsigned int j1, const unsigned int l1,
                          const unsigned int j2, const unsigned int l2,
                          const unsigned int k, const vector_t *t,
                          vector_t *V);

double cov_freq (const model_t *mdl, const unsigned int j,
                 const unsigned int k1, const unsigned int k2,
                 const unsigned int l1, const unsigned int l2,
                 const vector_t *t1, const vector_t *t2);

void diff_expect_freq (const model_t *mdl, const unsigned int j,
                       const unsigned int k, const unsigned int l,
                       const vector_t *t, const double scale,
                       const vector_t *V,
                       vector_t *Dm,
                       vector_t *Dt);

void diff_interact_freq (const model_t *mdl,
                         const unsigned int j1, const unsigned int l1,
                         const unsigned int j2, const unsigned int l2,
                         const unsigned int k, const vector_t *t,
                         const double scale, vector_t *V,
                         vector_t *Dm, vector_t *Dt);

double diff_naturalize_freq (const model_t *mdl, const unsigned int j,
                             vector_t *Dm, vector_t *Dt);

#endif /* !__VBNMR_DIST_NORMAL_H__ */

