
/* include the vbnmr header. */
#include <vbnmr/vbnmr.h>

/* vbnmr_init(): initialize the vfl library and include vbnmr features.
 *
 * returns:
 *  integer indicating success (1) or failure (0).
 */
int vbnmr_init (void) {
  /* declare a variable to store the initialization result. */
  int res = 1;

  /* initialize vfl. */
  if (!vfl_init())
    return 0;

  /* register vbnmr model types. */
  res &= vfl_register_model_type(vbnmr_model_vfgp);

  /* register vbnmr factor types. */
  res &= vfl_register_factor_type(vbnmr_factor_quad);

  /* return the result. */
  return res;
}

