//
// Created by Iván Yossi on 16/04/22.
//

#ifndef MOGEN_MGN_PARETO_H
#define MOGEN_MGN_PARETO_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

int vector_dominate(gsl_vector *u, gsl_vector *v)
{
    if (u->size != v->size) { return -1; }

    int all = 0;
    int any = 0;
    //  return all(u .<= v,dims=2) .& any(u .< v, dims=2)
    for (size_t i = 0; i < u->size; ++i) {
        all += u->data[i] > v->data[i];
        if(!any) {
            any = (u->data[i] < v->data[i])? 1:0;
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
            gsl_vector_view row = gsl_matrix_row(M,j);
            if(vector_dominate(&f_row.vector,&row.vector) < 0) {
                dominate[j]++;
            }
        }
    }
    return dominate;
}

#endif //MOGEN_MGN_PARETO_H
