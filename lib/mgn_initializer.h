//
// Created by Iván Yossi on 30/04/22.
//

#ifndef MOGEN_MGN_INITIALIZER_H
#define MOGEN_MGN_INITIALIZER_H

#include <stdlib.h>
#include <gsl/gsl_matrix.h>

#include "mgn_types.h"

struct mgnp_init_params {
    mgnLimit *limit;
    struct mgn_i_ops *ops;
};

void mgn_init_rand(void* in, void* limits);



// =========== Latin HyperCube ======
typedef struct pmgn_lhc_data_priv pmgn_lhcd_priv;
typedef struct mgn_lhc mgn_lhci;

struct mgn_lhc {
    size_t psize;
    size_t dim;
    mgnLimit *limits;
    gsl_matrix *m_points;
    pmgn_lhcd_priv *_p;
};

mgn_lhci* mgn_init_new_lhci(size_t psize, size_t dim, mgnLimit *lim);

void mgn_init_lhc_to_matrix(gsl_matrix *m_a, mgnLimit *lim);

void mgn_init_lhc(void *in, void* param);

void mgn_lhci_reset(mgn_lhci *lhc);

void mgn_lhci_free(mgn_lhci *lhci);

//=======================

//void mgn_init_LHC_init(size_t psize, size_t dim, mgnLimit *limits);

//gsl_matrix* mgn_LHC_get();

#endif //MOGEN_MGN_INITIALIZER_H
