
/* ensure once-only inclusion. */
#ifndef __VBNMR_DIST_BASIS_H__
#define __VBNMR_DIST_BASIS_H__

/* include the model header. */
#include <vbnmr/model.h>

/* include the distribution headers. */
#include <vbnmr/dist-normal.h>
#include <vbnmr/dist-gamma.h>

/* function declarations: */

double expect_basis (const model_t *mdl, const unsigned int j,
                     const unsigned int k, const unsigned int l,
                     const vector_t *t);

double interact_basis (const model_t *mdl,
                       const unsigned int j1, const unsigned int j2,
                       const unsigned int l1, const unsigned int l2,
                       const unsigned int k, const vector_t *t);

double cov_basis (const model_t *mdl, const unsigned int j,
                  const unsigned int k1, const unsigned int k2,
                  const vector_t *t1, const vector_t *t2);

#endif /* !__VBNMR_DIST_BASIS_H__ */

