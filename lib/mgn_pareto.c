/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_pareto.h"

#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "mgn_poplist.h"

// asume vectors are same size
int vector_dominate(gsl_vector *u, gsl_vector *v)
{
//    if (u->size != v->size) { return -1; }
//    int dom = 1;
//    double min;
//    double max;
//
//    gsl_vector *c = gsl_vector_alloc(u->size);
//    gsl_vector_memcpy(c,u);
//    gsl_vector_sub(c,v);
//
//    gsl_vector_minmax(c,&min,&max);
//
//    if (max <= 0 && min < 0) {
//        dom = -1;
//    }
//    return dom;

    // default u not dom v
    int all = 0;
    int any = 0;

    size_t size = u->size;
    double *ud = u->data;
    double *vd = v->data;
    //  return all(u .<= v,dims=2) .& any(u .< v, dims=2)

    for (size_t i = 0; i < size; ++i) {
        all += ud[i] > vd[i];
    }

    if (all == 0) {
        for (size_t i = 0; i < size; ++i) {
            any = (ud[i] < vd[i]) ? 1 : 0;
            if (any) break;
        }
    }

    return (any && all == 0)? -1 : 1;
}

int* gsl_matrix_pareto_rank(gsl_matrix *M)
{
    int *dominate = calloc(M->size1, sizeof(*dominate));
    for (size_t i = 0; i < M->size1; ++i) {
        gsl_vector_view f_row = gsl_matrix_row(M,i);
        for (size_t j = 0; j < M->size1; ++j) {
            if (j == i) continue;
            gsl_vector_view row = gsl_matrix_row(M,j); // !
            if(vector_dominate(&f_row.vector,&row.vector) < 0) {
                dominate[j]++;
            }
        }
    }
    return dominate;
}

bool mgn_pop_insert_dom(mgn_popl *popl, mgn_indv *indv)
{
    bool inset_in = true;

    gsl_vector *in_fval = indv->f;

    mgn_popl_cursor_reset(popl);
    while(mgn_popl_current(popl) != 0) {
        gsl_vector *ep_fval = popl->ops->get_iparams(mgn_popl_current(popl)).f;
        if(vector_dominate(in_fval,ep_fval) < 0) {
            void* indv_l = mgn_popl_pop_current(popl); // advances cursor
            popl->ops->free(indv_l);
            free(indv_l);
        } else {
            mgn_popl_next(popl);
        }
    }

    mgn_popl_cursor_reset(popl);

    while(mgn_popl_current(popl) != 0) {
        gsl_vector *ep_fval = popl->ops->get_iparams(mgn_popl_current(popl)).f;
        if(vector_dominate(ep_fval,in_fval) < 0) {
            inset_in = false;
            break;
        }
        mgn_popl_next(popl);
    }
//    printf("size solin2 %zu\n", popl->ops->get_iparams(
//        popl->get(popl,0)).x->size
//    );

    mgn_popl_cursor_reset(popl);
    if (inset_in) {
        mgn_indv *last = mgn_popl_alloc_last(popl);
        popl->ops->copy(last, indv);
    }

    return inset_in;
}


//
///**
// * To verify dominance between two solutions
// * Returns
// *  1: if x weak dominates B
// *  0: otherwise
// *  */
//int weak_dominance(double *x, double *y, int dim)
//{
//    int flag = 1;
//    for (int i = 0; i < dim; ++i) {
//        if (x[i] > y[i]) {
//            flag = 0;
//        }
//    }
//    return flag;
//}
//
///**
// * To verify dominance between two solutions
// * @param A:    vector 'A'
// * @param B:    vector'B'
// * @param dim:  dimension of vectors
// * @return:
// *   1: if A dominates B
// * 	-1: if B dominates A
// *   0: if they are non-dominated
// *   2:	if A = B
// */
//int dominance(double *A, double *B, int dim)
//{
//    int j, countA = 0, countB = 0;
//
//    for (j = 0; j < dim; j++){
//        if (A[j] < B[j]) {
//            countA++;
//
//        } else if (B[j] < A[j]) {
//            countB++;
//        }
//    }
//    if (countA > 0 && countB == 0){
//        return 1;
//    }
//    if (countB > 0 && countA == 0){
//        return -1;
//    }
//    if (countA == 0 && countB == 0){
//        return 2;
//    }
//    return 0;
//}
