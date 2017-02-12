
/* ensure once-only inclusion. */
#ifndef __VBNMR_VECTOR_H__
#define __VBNMR_VECTOR_H__

/* include c library headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* vector_disp(): macro function for displaying the contents of vectors.
 */
#define vector_disp(v) vector_dispfn(v, #v)

/* vector_t: structure for holding a real vector.
 */
typedef struct {
  /* data access parameters:
   *  @len: number of vector elements.
   *  @stride: array spacing between elements.
   */
  unsigned int len;
  unsigned int stride;

  /* @data: array or pointer to vector elements.
   */
  double *data;
}
vector_t;

/* vector_view_t: aliased type used to explicitly note that a given
 * vector data structure does not own its data array.
 */
typedef vector_t vector_view_t;

/* function declarations (vector.c): */

vector_t *vector_alloc (const unsigned int len);

void vector_copy (vector_t *dest, const vector_t *src);

void vector_free (vector_t *v);

vector_view_t vector_view_array (double *data, const unsigned int len);

vector_view_t vector_subvector (const vector_t *v,
                                const unsigned int offset,
                                const unsigned int len);

double vector_get (const vector_t *v, const unsigned int i);

double vector_max (const vector_t *v);

void vector_set (vector_t *v, const unsigned int i, const double vi);

void vector_set_all (vector_t *v, const double vall);

void vector_set_zero (vector_t *v);

void vector_add (vector_t *a, const vector_t *b);

void vector_add_const (vector_t *v, const double beta);

void vector_scale (vector_t *v, const double alpha);

int vector_equal (const vector_t *a, const vector_t *b);

int vector_positive (const vector_t *v);

void vector_dispfn (const vector_t *v, const char *str);

#endif /* !__VBNMR_VECTOR_H__ */

