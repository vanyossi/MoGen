#ifndef _MGN_MOEAD_SOLVER_
#define _MGN_MOEAD_SOLVER_

#include "mgn_types.h"

#include <stdbool.h>
#include <gsl/gsl_vector.h>

typedef struct moead_features moeadf;

moeadf* mgn_moead_getfeatures(mgnMoa* moead);

mgn_pop* mgn_moead_getpop(mgnMoa* moead);

void mgnp_moead_update_z(gsl_vector* x, gsl_vector* f, gsl_vector* g, moeadf* param);

void moead_run(mgnMoa* moead);

bool moead_stop();

mgnMoa* mgn_moead_init(size_t H,size_t nobj, size_t T
                       ,mgn_popl *epop
                       ,mgnMop *mop
                       ,void (*apply)(void*, void*), void* params);

moeadf* mgn_moead_getfeatures(mgnMoa* moead);

mgn_pop* mgn_moead_getpop(mgnMoa* moead);

void mgn_moead_set_mop(mgnMoa* moead, mgnMop *mop, int type);

void moead_set_prob(mgnMoa*moead, mgn_ga_sets *gasets);

void mgn_moead_free(mgnMoa *moead);

#endif // _MGN_MOEAD_SOLVER_
