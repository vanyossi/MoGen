#ifndef _MGN_VECTOR_DISTANCE_LIB_
#define _MGN_VECTOR_DISTANCE_LIB_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "gsl_vector_additional.h"

// Vec is array of vectors.
gsl_matrix* gsl_vector_distance_matrix1(gsl_vector *Vec, size_t size, double pval)
{
    gsl_matrix *out = gsl_matrix_alloc(size, size);
    gsl_vector_view diag = gsl_matrix_diagonal(out);
    gsl_vector_set_all(&diag.vector,0);

    for (size_t i = 0; i < size; ++i) {
        double dist;
        for (size_t j = i+1; j < size; ++j) {
            gsl_vector *tmp_vec = gsl_vector_alloc(Vec->size);
            gsl_vector_memcpy(tmp_vec, &Vec[j]);
            gsl_vector_sub(tmp_vec, &Vec[i]);
            if (pval != 2.0) {
                dist = gsl_vector_pnorm(tmp_vec, pval);
            } else {
                dist = gsl_blas_dnrm2(tmp_vec);
            }
            gsl_matrix_set(out,i,j,dist);
            gsl_matrix_set(out,j,i,dist);
            gsl_vector_free(tmp_vec);
        }
    }
    return out;
}

gsl_matrix* gsl_vector_distance_matrix(gsl_matrix *vecmat, double pval)
{
    size_t size = vecmat->size1;
    gsl_matrix *out = gsl_matrix_alloc(size,size);
    gsl_vector_view diag = gsl_matrix_diagonal(out);
    gsl_vector_set_all(&diag.vector,0);

    for (size_t i = 0; i < size; ++i) {
        double dist;
        gsl_vector_view rowi = gsl_matrix_row(vecmat,i);
        for (size_t j = i+1; j < size; ++j) {
            gsl_vector *tmp_vec = gsl_vector_alloc(vecmat->size2);
            gsl_vector_view rowj = gsl_matrix_row(vecmat,j);

            gsl_vector_memcpy(tmp_vec, &rowj.vector);
            gsl_vector_sub(tmp_vec, &rowi.vector);
            if (pval != 2.0) {
                dist = gsl_vector_pnorm(tmp_vec, pval);
            } else {
                dist = gsl_blas_dnrm2(tmp_vec);
            }
            gsl_matrix_set(out,i,j,dist);
            gsl_matrix_set(out,j,i,dist);
            gsl_vector_free(tmp_vec);
        }
    }
    return out;
}

//void gsl_matrix_order_row_vals(gsl_matrix *mat)
//{
//
//    return;
//}
#endif // __MGN_VECTOR_DISTANCE_LIB_
