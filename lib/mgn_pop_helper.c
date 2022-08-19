/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_pop_helper.h"

//TODO need proto indv (for different indv)

void mgn_pop_copy_mp(mgn_pop_proto *pop, mgn_pop_matrix *mpop){

    for (size_t i = 0; i < mpop->x->size1; ++i) {
        mgn_indv *indv = NULL;
        indv = mgn_indv_alloc(indv, pop->ops, &pop->iparams);

        gsl_matrix_get_row(indv->x,mpop->x,i);
        gsl_matrix_get_row(indv->f,mpop->f,i);
        gsl_matrix_get_row(indv->g,mpop->g,i);

        pop->set(pop,indv,i);
    }
}

mgn_pop_matrix* mgn_pop_to_popm(mgn_pop *pop)
{
    mgn_pop_matrix *pop_m = mgn_pop_matrix_alloc(
        pop->size,pop->iparams.x_size
        ,pop->size,pop->iparams.f_size
        ,pop->size,pop->iparams.g_size
        );

    // TODO indv prototype need
    for (size_t i = 0; i < pop->size; ++i) {
        mgn_indv *in = mgn_pop_get(pop,i);
        gsl_matrix_set_row(pop_m->x,i,in->x);
        gsl_matrix_set_row(pop_m->f,i,in->f);
        gsl_matrix_set_row(pop_m->g,i,in->g);
    }
    return pop_m;
}

// copy pop_matrix pop to mgn_pop
mgn_pop* mgn_pop_matrix_to_pop(mgn_pop_matrix *pop_m
                               ,mgn_indv_ops *ops
                               ,mgn_indv_param *params)
{
    mgn_pop *pop = mgn_pop_alloc(pop_m->size,ops,params);
    // TODO indv prototype need
    for (size_t i = 0; i < pop_m->size; ++i) {
        mgn_indv *in = mgn_pop_get(pop,i);
        gsl_matrix_get_row(in->x,pop_m->x,i);
        gsl_matrix_get_row(in->f,pop_m->f,i);
        gsl_matrix_get_row(in->g,pop_m->g,i);
    }
    return pop;
}
