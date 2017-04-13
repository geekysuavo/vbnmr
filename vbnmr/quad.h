
/* ensure once-only inclusion. */
#ifndef __VBNMR_QUAD_H__
#define __VBNMR_QUAD_H__

/* include the vfl factor header. */
#include <vfl/factor.h>

/* function declarations (quad.c): */

int quad_set_dims (factor_t *f, const unsigned int D);

int quad_set_ftsize (factor_t *f, const unsigned int n);

/* available factor types: */

extern const factor_type_t *vbnmr_factor_quad;

#endif /* !__VBNMR_QUAD_H__ */

