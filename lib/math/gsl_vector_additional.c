//
// Created by Iván Yossi on 20/04/22.
//

#include "gsl_vector_additional.h"

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_blas.h>

#include "mgn_types.h"

// private structs
struct gsl_vector_qsort_idx {
    unsigned int idx;
    double value;
};

// private functions
double hpow(double value, void* pvalue) {
    double *p = (double*)pvalue;
    return fabs(pow(value, *p));
}

int cmp_double(const void *a, const void *b)
{
    struct gsl_vector_qsort_idx *ad = (struct gsl_vector_qsort_idx*)a;
    struct gsl_vector_qsort_idx *bd = (struct gsl_vector_qsort_idx*)b;

    return (ad->value < bd->value)? -1 : 1;;
}

// public functions
double mgn_fabs(double value, void* nouse) {
    UNUSED(nouse);
    return fabs(value);
}

void gsl_vector_map(gsl_vector *V, double (*func)(double, void*), void* param)
{
    size_t size = V->size;
    for (size_t i = 0; i < size; i++) {
        gsl_vector_set(V,i, func(gsl_vector_get(V,i), param) );
    }
}

double gsl_vector_pnorm(gsl_vector *v, double pvalue)
{
    gsl_vector_map(v,hpow, (void*)&pvalue);
    double sumv = gsl_vector_sum(v);
    return pow(sumv,1.0/pvalue);
}

int* gsl_vector_qsort(gsl_vector *vec)
{
    int *index = (int*)calloc(vec->size, sizeof(int));
    struct gsl_vector_qsort_idx *idata = calloc(vec->size, sizeof(*idata));

    for (size_t i = 0; vec->size > i; ++i) {
        idata[i].idx = i;
        idata[i].value = gsl_vector_get(vec,i);
    }
    qsort(idata,vec->size,sizeof(struct gsl_vector_qsort_idx),cmp_double);

    for (size_t i = 0; vec->size > i; ++i) {
        gsl_vector_set(vec,i,idata[i].value);
        index[i] = (int)idata[i].idx;
    }
    free(idata);
    return index;
}

void gsl_vector_set_seq(gsl_vector *vec)
{
    for (size_t i = 0; i < vec->size; ++i) {
        gsl_vector_set(vec,i,(double)i);
    }
    return;
}

void
gsl_vector_repeat(gsl_vector *v, size_t rep, gsl_matrix *C)
{
//    gsl_matrix *result = gsl_matrix_calloc(rep,v->size);
    gsl_matrix *ones = gsl_matrix_calloc(rep,1);
    gsl_matrix_add_constant(ones,1);

    gsl_matrix *vtmp = gsl_matrix_alloc(1, v->size);
    gsl_matrix_set_row(vtmp,0,v);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                   1.0, ones, vtmp, 0.0, C);

    gsl_matrix_free(ones);
    gsl_matrix_free(vtmp);

//    return result;
}
