
/* ensure once-only inclusion. */
#ifndef __VBNMR_RNG_H__
#define __VBNMR_RNG_H__

/* include c library headers. */
#include <stdlib.h>
#include <math.h>

/* rng_t: structure for holding a pseudorandom number generator.
 */
typedef struct {
  /* current generator state:
   *  @u, @v, @w: core generator variables.
   */
  unsigned long long u, v, w;

  /* generator parameters:
   *  @seed: initial random seed.
   */
  unsigned long long seed;
}
rng_t;

/* function declarations (rng.c): */

rng_t *rng_alloc (void);

rng_t *rng_alloc_with_seed (const unsigned long long seed);

void rng_free (rng_t *gen);

double rng_uniform (rng_t *gen);

double rng_normal (rng_t *gen);

#endif /* !__VBNMR_RNG_H__ */

