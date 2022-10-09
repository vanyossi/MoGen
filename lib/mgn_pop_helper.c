/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_pop_helper.h"

#include <float.h>

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

void mgn_popl_insert_popm(mgn_popl *dst, mgn_pop_matrix *src)
{
    for (size_t i = 0; i < src->size; ++i) {
        mgn_indv *indv = mgn_popl_alloc_last(dst);
        gsl_matrix_get_row(indv->x, src->x,i);
        gsl_matrix_get_row(indv->f, src->f,i);
        gsl_matrix_get_row(indv->g, src->g,i);
    }
}

void mgn_popl_insert_pop(mgn_popl *dst, mgn_pop_proto *src)
{
    for (size_t i = 0; i < src->size; ++i) {
        mgn_indv *indv = mgn_popl_alloc_last(dst);
        dst->ops->copy(indv, src->get(src,i));
    }
}

mgn_pop_matrix* mgn_pop_to_popm(mgn_pop_proto *pop)
{
    mgn_pop_matrix *pop_m = mgn_pop_matrix_alloc(
        pop->size,pop->iparams.x_size
        ,pop->size,pop->iparams.f_size
        ,pop->size,pop->iparams.g_size
        );

    // TODO indv prototype need
    for (size_t i = 0; i < pop->size; ++i) {
        mgn_indv *in = pop->get(pop,i);
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


gsl_vector* mgn_pop_getbound_column(mgn_pop_proto *pop, bool ismax)
{
    gsl_vector *lim = gsl_vector_alloc(pop->iparams.x_size);
    gsl_vector_set_all(lim,(ismax)? DBL_MIN : DBL_MAX);

    double *xv, *x, *xtmp;
    xv = lim->data;
    for (size_t i = 0; i < pop->size; ++i) {
        x = pop->ops->get_iparams(pop->get(pop,i)).x->data;
        for (size_t xi = 0; xi < lim->size; ++xi) {
            if(ismax) {
                xv[xi] = (x[xi] > xv[xi])? x[xi] : xv[xi];
            } else {
                xv[xi] = (x[xi] < xv[xi])? x[xi] : xv[xi];
            }
        }
    }
    return lim;
}

gsl_vector* mgn_pop_max_column(mgn_pop_proto *pop)
{
    return mgn_pop_getbound_column(pop,true);
}

gsl_vector* mgn_pop_min_column(mgn_pop_proto *pop)
{
    return mgn_pop_getbound_column(pop,false);
}
