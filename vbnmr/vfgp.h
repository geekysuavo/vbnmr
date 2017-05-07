
/* ensure once-only inclusion. */
#ifndef __VBNMR_VFGP_H__
#define __VBNMR_VFGP_H__

/* include the vfl model header. */
#include <vfl/model.h>

/* function declarations (vfgp.c): */

int vfgp_set_mode (model_t *mdl, const unsigned int gp_enable);

/* available model types: */

extern const model_type_t *vbnmr_model_vfgp;

#endif /* !__VBNMR_VFGP_H__ */

