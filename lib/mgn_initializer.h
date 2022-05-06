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
    struct _mgn_i_ops *ops;
};

void mgn_init_rand(void* in, void* limits);

void mgn_init_LHC_init(size_t psize, size_t dim, mgnLimit *limits);

void mgn_init_lhc(void *in, void* param);

gsl_matrix* mgn_LHC_get();

#endif //MOGEN_MGN_INITIALIZER_H
