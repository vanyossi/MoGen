/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef MOGEN_MGN_POP_MATRIX_H
#define MOGEN_MGN_POP_MATRIX_H

#include <stdlib.h>
#include <gsl/gsl_matrix.h>

#include "mgn_types.h"

typedef struct mgnp_pop_matrix_data mgn_pop_matrix;

struct mgnp_pop_matrix_data {
    size_t size;
    gsl_matrix *x;
    gsl_matrix *f;
    gsl_matrix *g;
};

mgn_pop_matrix* mgn_pop_matrix_alloc(
    size_t x_size, size_t x_dim
    ,size_t f_size, size_t f_dim
    ,size_t g_size, size_t g_dim
    );

mgn_pop_matrix* mgn_pop_matrix_alloc_empty();

void mgn_pop_matrix_free(mgn_pop_matrix *popm);

void mgn_pop_matrix_eval(mgn_pop_matrix *popm, mgnMop *mop);


#endif //MOGEN_MGN_POP_MATRIX_H
