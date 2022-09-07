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
int cmp_double(const void *a, const void *b)
{
    struct gsl_vector_qsort_idx *ad = (struct gsl_vector_qsort_idx*)a;
    struct gsl_vector_qsort_idx *bd = (struct gsl_vector_qsort_idx*)b;

    return (ad->value <= bd->value)? -1 : 1;
}

// public functions
void gsl_vector_map(gsl_vector *V, double (*func)(double, void*), void* param)
{
    size_t size = V->size;
    for (size_t i = 0; i < size; i++) {
        gsl_vector_set(V,i, func(gsl_vector_get(V,i), param) );
    }
}

double gsl_vector_pnorm(gsl_vector *v, double pvalue)
{
    gsl_vector *tmp_vec = gsl_vector_alloc(v->size);
    gsl_vector_memcpy(tmp_vec,v);

    gsl_vector_map(tmp_vec, map_hpow, (void *) &pvalue);
    double sumv = gsl_vector_sum(tmp_vec);

    gsl_vector_free(tmp_vec);
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

gsl_vector* gsl_vector_get_indexes(gsl_vector *v, int *idx, size_t size)
{
    gsl_vector *out = gsl_vector_calloc(size);
    for (size_t i = 0; i < size; ++i) {
        gsl_vector_set(out,i,gsl_vector_get(v,idx[i]));
    }
    return out;
}

gsl_matrix* gsl_matrix_get_row_indexes(gsl_matrix *v, int *idx, size_t size)
{
    gsl_matrix *out = gsl_matrix_calloc(size,v->size2);
    for (size_t i = 0; i < size; ++i) {
        gsl_vector_view vecv = gsl_matrix_row(v,idx[i]);
        gsl_matrix_set_row(out,i,&vecv.vector);
    }
    return out;
}

void gsl_matrix_printf(gsl_matrix *M, FILE *stream)
{
    const char *format = "%1.6e";

    for (size_t c = 0; c < M->size1; ++c) {
        gsl_vector_view crow = gsl_matrix_row(M, c);
        for (size_t val = 0; val < crow.vector.size; ++val) {
            fprintf(stream,format, gsl_vector_get(&crow.vector, val));
            if(val != crow.vector.size -1) {
                fprintf(stream," ");
            }
        }
        fprintf(stream, "\n");
    }
}

void gsl_matrix_save(gsl_matrix *M, char* filename){
    return;
    FILE *out = fopen(filename, "w");
    gsl_matrix_printf(M, out);
    fclose(out);
}

// general print matrix
//    gsl_matrix_printf(m_dist, stdout);
//    gsl_matrix *m_dist_i_double = gsl_matrix_alloc(m_dist_i->size1, m_dist_i->size2);
//    for (size_t i = 0; i < m_dist_i->size1; ++i) {
//        for (size_t j = 0; j < m_dist_i->size2; ++j) {
//            gsl_matrix_set(m_dist_i_double,i,j,(double)gsl_matrix_int_get(m_dist_i,i,j));
//        }
//    }
//    puts("====");
//    gsl_matrix_printf(m_dist_i_double, stdout);
//    gsl_matrix_free(m_dist_i_double);
