//
// Created by Iván Yossi on 30/04/22.
//

#ifndef MOGEN_MGN_INITIALIZER_H
#define MOGEN_MGN_INITIALIZER_H

#include <stdlib.h>
#include <gsl/gsl_matrix.h>

#include "mgn_types.h"

typedef struct mgnp_init_params mgn_initializer;
typedef void (*mgnf_init_start)(mgn_pop_proto*, mgn_initializer*, void*);
typedef void (*pmgnf_init_free)(mgn_initializer*);


// does nothing
// exists to avoid breaking code while refactor
typedef void (*mgnf_init_transition)(void*, void*);
void mgn_init_transition(void *in, void* param);


void mgn_pinit_free(mgn_initializer *idata);

struct mgnp_init_params {
    mgnLimit *limit;
    mgnf_init_start start;
    pmgnf_init_free free;
};



mgn_initializer* mgn_pinit_rand_alloc(mgnLimit *limit);
void mgn_pinit_rand_free(mgn_initializer *idata);

void mgn_init_pop_rand(mgn_pop_proto *pop, mgn_initializer *init, void *extra);



// =========== Latin HyperCube ======

mgn_initializer* mgn_pinit_lhc_alloc(mgn_pop_proto *pop, mgnLimit *limit);
void mgn_pinit_lhc_free(mgn_initializer *idata);

void mgn_init_pop_lhc(mgn_pop_proto *pop, mgn_initializer *init, void *extra);

void mgn_init_lhc_to_matrix(gsl_matrix *m_a, mgnLimit *lim);





//mgn_lhci* mgn_init_new_lhci(size_t psize, size_t dim, mgnLimit *lim);

//void mgn_init_lhc(void *in, void* param);

//void mgn_lhci_reset(mgn_lhci *lhc);

//void mgn_lhci_free(mgn_lhci *lhci);

//=======================

//void mgn_init_LHC_init(size_t psize, size_t dim, mgnLimit *limits);

//gsl_matrix* mgn_LHC_get();

#endif //MOGEN_MGN_INITIALIZER_H
