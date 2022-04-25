//
// Created by Iván Yossi on 16/04/22.
//

#ifndef MOGEN_MGN_PARETO_H
#define MOGEN_MGN_PARETO_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

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

#endif //MOGEN_MGN_PARETO_H
