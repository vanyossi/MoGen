/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_pop_matrix.h"

#include "mgn_mop.h"

mgn_pop_matrix* mgn_pop_matrix_alloc(
    size_t x_size, size_t x_dim
    ,size_t f_size, size_t f_dim
    ,size_t g_size, size_t g_dim
){
    mgn_pop_matrix *popm_i = malloc(sizeof(*popm_i));
    popm_i->size = x_size;
    popm_i->x = gsl_matrix_calloc(x_size, x_dim);
    popm_i->f = gsl_matrix_calloc(f_size, f_dim);
    popm_i->g = gsl_matrix_calloc(g_size, g_dim);

    return popm_i;
}

mgn_pop_matrix* mgn_pop_matrix_alloc_empty()
{
    mgn_pop_matrix *popm_i = malloc(sizeof(*popm_i));
    popm_i->x = NULL;
    popm_i->f = NULL;
    popm_i->g = NULL;
    return popm_i;
}

void mgn_pop_matrix_free(mgn_pop_matrix *popm)
{
    gsl_matrix_free(popm->x);
    gsl_matrix_free(popm->f);
    gsl_matrix_free(popm->g);
    free(popm);
}


void mgn_pop_matrix_eval(mgn_pop_matrix *popm, mgnMop *mop)
{
    for (size_t i = 0; i < popm->x->size1; ++i) {
        gsl_vector_view x_row = gsl_matrix_row(popm->x,i);
        gsl_vector_view f_row = gsl_matrix_row(popm->f,i);
        gsl_vector_view g_row = gsl_matrix_row(popm->g,i);

        if (mop->eval) {
            mop->eval(&x_row.vector, &f_row.vector, &g_row.vector, mop->params);
        } else if (mop->eval_array) {
            mop->eval_array(x_row.vector.data, f_row.vector.data, g_row.vector.data, mop->params);
        }
    }
}


