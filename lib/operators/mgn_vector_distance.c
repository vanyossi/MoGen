/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_vector_distance.h"

#include <gsl/gsl_blas.h>

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

gsl_matrix* gsl_matrix_dist(gsl_matrix *m_a, gsl_matrix *m_b, double pval)
{
    gsl_matrix *out = gsl_matrix_alloc(m_a->size1,m_b->size1);
//    gsl_vector_view diag = gsl_matrix_diagonal(out);
//    gsl_vector_set_all(&diag.vector,0);

    for (size_t i = 0; i < m_a->size1; ++i) {
        double dist;
        gsl_vector_view rowi = gsl_matrix_row(m_a,i);
        for (size_t j = 0; j < m_b->size1; ++j) {
            gsl_vector *tmp_vec = gsl_vector_alloc(m_a->size2);
            gsl_vector_view rowj = gsl_matrix_row(m_b,j);

            gsl_vector_memcpy(tmp_vec, &rowj.vector);
            gsl_vector_sub(tmp_vec, &rowi.vector);
            if (pval != 2.0) {
                dist = gsl_vector_pnorm(tmp_vec, pval);
            } else {
                dist = gsl_blas_dnrm2(tmp_vec);
            }
            gsl_matrix_set(out,i,j,dist);
//            gsl_matrix_set(out,j,i,dist);
            gsl_vector_free(tmp_vec);
        }
    }
    return out;
}

void gsl_matrix_distrank_index(gsl_matrix* m_dist, gsl_matrix_int *m_out)
{
    for (size_t i = 0; i < m_dist->size1; ++i) {
        gsl_vector_view crow = gsl_matrix_row(m_dist,i);
        int *dorder = gsl_vector_qsort(&crow.vector);
        gsl_vector_int_view irank = gsl_vector_int_view_array(dorder,m_dist->size2);
        // out
        gsl_matrix_int_set_row(m_out,i,&irank.vector);
        free(dorder);
    }
}

//void gsl_matrix_order_row_vals(gsl_matrix *mat)
//{
//
//    return;
//}
