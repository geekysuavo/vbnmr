
/* ensure once-only inclusion. */
#ifndef __VBNMR_DIST_GAMMA_H__
#define __VBNMR_DIST_GAMMA_H__

/* include the model header. */
#include <vbnmr/model.h>

/* function declarations (dist-gamma.c): */

double gam_entropy (const model_t *mdl,
                    const unsigned int j,
                    const unsigned int d);

double gam_prior (const model_t *mdl,
                  const unsigned int j,
                  const unsigned int d);

double expect_decay (const model_t *mdl, const unsigned int j,
                     const unsigned int delta, const vector_t *t);

double expect_all_decay (const model_t *mdl, const unsigned int j,
                         const vector_t *t, vector_t *R);

double interact_decay (const model_t *mdl,
                       const unsigned int j1, const unsigned int j2,
                       const unsigned int delta, const vector_t *t);

double interact_all_decay (const model_t *mdl,
                           const unsigned int j1, const unsigned int j2,
                           const vector_t *t, vector_t *R);

double cov_decay (const model_t *mdl, const unsigned int j,
                  const vector_t *t1, const vector_t *t2);

void diff_expect_decay (const model_t *mdl, const unsigned int j,
                        const vector_t *t, const double scale,
                        const vector_t *R,
                        vector_t *Da,
                        vector_t *Db);

void diff_interact_decay (const model_t *mdl,
                          const unsigned int j1, const unsigned int j2,
                          const vector_t *t, const double scale,
                          const vector_t *R,
                          vector_t *Da,
                          vector_t *Db);

double diff_naturalize_decay (const model_t *mdl, const unsigned int j,
                              vector_t *Da, vector_t *Db);

#endif /* !__VBNMR_DIST_GAMMA_H__ */

