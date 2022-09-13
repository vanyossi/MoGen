/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_rbf.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

#include "mgn_cluster_m.h"
#include "gsl_vector_additional.h"

// Kernels
// d <- distance || x - y ||
void rbf_kernel_gauss(gsl_vector *r, double sigma)
{
    // -r^2 / 2sigma
//    gsl_vector_scale(r,sigma);
//    gsl_vector_mul(r,r);
//    gsl_vector_scale(r,-1.0);
//    gsl_vector_map(r, map_exp,NULL);

    // as in paper
//    gsl_vector_scale(r,-1);
    gsl_vector_scale(r,-1.0 / (2 * pow(sigma,2)));
    gsl_vector_map(r, map_exp,NULL);
}

void rbf_kernel_mqua(gsl_vector *r, double sigma)
{
    gsl_vector_scale(r,pow(sigma,2));
//    gsl_vector_mul(r, r);
    gsl_vector_add_constant(r,1);
    gsl_vector_map(r,map_sqrt,NULL);
}

void rbf_kernel_imqua(gsl_vector *r, double sigma)
{
    gsl_vector_scale(r,pow(sigma,2));
//    gsl_vector_mul(r, r);
    gsl_vector_add_constant(r,1);
    gsl_vector_map(r,map_sqrt,NULL);
    gsl_vector_map(r,map_invmul,NULL);
}

gsl_matrix*
mgn_rbf_create_phi(gsl_matrix *X
                   , cluster_data *km
                   , gsl_vector *sigma
                   , void (*rbf)(gsl_vector *r, double s)
                   , gsl_matrix *m_phi)
{
    if (m_phi == NULL){
        m_phi = gsl_matrix_alloc(X->size1,km->k);
    }

    double r;
    gsl_matrix *m_rep_c = gsl_matrix_alloc(X->size1, X->size2);
    gsl_vector *v_work_c = gsl_vector_alloc(X->size2);
    gsl_vector_view v_work_col;
    for (size_t i = 0; i < km->k; ++i) {
        gsl_matrix_get_row(v_work_c,km->centers,i);
        gsl_vector_repeat(v_work_c,X->size1,m_rep_c);
        gsl_matrix_sub(m_rep_c,X);
        for (size_t j = 0; j < X->size1; ++j) {
            gsl_matrix_get_row(v_work_c,m_rep_c,j);
            r = gsl_blas_dnrm2(v_work_c); // best a
//            r = gsl_vector_pnorm(v_work_c,1.0);
//            r = gsl_blas_dasum(v_work_c);
//            r = pow(r,2);
            gsl_matrix_set(m_phi,j,i,r);
        }
        v_work_col = gsl_matrix_column(m_phi,i);
        rbf(&v_work_col.vector,gsl_vector_get(sigma,i));
    }
    gsl_vector_free(v_work_c);
    gsl_matrix_free(m_rep_c);
    return m_phi;
}

// returns new W matrix vector, or uses m_w
gsl_matrix *
mgn_rbf_new_weight(gsl_matrix *m_phi, gsl_matrix *y, gsl_matrix *m_w)
{
    gsl_error_handler_t *e_handle = gsl_set_error_handler_off();

    gsl_matrix *m_x_xt = gsl_matrix_calloc(m_phi->size2,m_phi->size2);
    gsl_blas_dsyrk(CblasUpper,CblasTrans,1,m_phi,1,m_x_xt);
    // complete symmetry
    for (size_t i = 0; i < m_x_xt->size1; ++i) {
        for (size_t j = 0; j < i; ++j) {
            gsl_matrix_set(m_x_xt,i,j,gsl_matrix_get(m_x_xt,j,i));
        }
    }

    gsl_matrix *m_phiy = gsl_matrix_calloc(m_phi->size2, y->size2);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,m_phi,y,0,m_phiy);

    if (m_w == NULL){
        m_w = gsl_matrix_calloc(m_phi->size2,m_phiy->size2);
    }
    // solve the system Ax = b
    // A = m_x_xt,  x = w,  b = m_phiy
    int s;
    bool err = false;
    gsl_matrix *m_x_xt_cpy = gsl_matrix_alloc(m_x_xt->size1, m_x_xt->size2);
    gsl_vector *w_vec = gsl_vector_alloc (m_phi->size2);
    gsl_permutation *p = gsl_permutation_alloc(m_x_xt->size1);
    gsl_vector_view y_i;
    for (size_t i = 0; i < m_phiy->size2; ++i) {
        bool error = 0;
        gsl_matrix_memcpy(m_x_xt_cpy,m_x_xt);
        y_i = gsl_matrix_column(m_phiy,i);

        gsl_linalg_LU_decomp (m_x_xt_cpy, p, &s);

        error = gsl_linalg_LU_solve (m_x_xt_cpy, p, &y_i.vector, w_vec);
        if (error !=0) {
            err = true;
        }
//        double det = gsl_linalg_LU_det(m_x_xt_cpy, 0);
        // dont try to solve if there is no solution.
//        printf("det == %g\n", det);
//        if (det >= 0) {
//        } else {
//            err= true;
//            break;
//        }
//        gsl_vector_fprintf(stdout,w_vec,"%.7f "); puts("");

        gsl_matrix_set_col(m_w,i,w_vec);
    }

//    printf("m_x_xt size %zu %zu, %zu, %zu, %zu\n", m_x_xt->size1, m_x_xt->size2, p->size, y_i.vector.size, w_vec->size);
//    gsl_matrix_printf(m_w,stdout);
//    gsl_vector_fprintf(stdout, w_vec, "%.4f ");

//    printf("=================\n");
//    FILE *data = fopen("mphi_inv.txt","w");
//    gsl_matrix_printf(m_x_xt,data);
//    fclose(data);

    gsl_matrix_free(m_x_xt_cpy);
    gsl_matrix_free(m_x_xt);
    gsl_matrix_free(m_phiy);
    gsl_vector_free(w_vec);
    gsl_permutation_free(p);

    if (err) {
        gsl_matrix_set_all(m_w,NAN);
    }
    gsl_set_error_handler(e_handle);
    return m_w;
    // multiply X * X' -> A
    // CblasNoTrans, A*A'
    // CblasUpper only upper triangle stored
    // int gsl_blas_dsyrk(CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans, double alpha, const gsl_matrix *A, double beta, gsl_matrix *m_x_xt)
//    int gsl_linalg_LU_invert (const gsl_matrix * LU, const gsl_permutation * p, gsl_matrix * inverse)
}


// helpers should probably go in another file
// caluclate mean square error
double mgn_math_mse(const gsl_vector *y, gsl_vector *yp)
{
    gsl_vector *y_copy = gsl_vector_calloc(y->size);
    gsl_vector_memcpy(y_copy,yp);

    gsl_vector_sub(y_copy,y);
    gsl_vector_mul(y_copy,y_copy);

    double sum_y = gsl_vector_sum(y_copy);

    gsl_vector_free(y_copy);
    return sum_y / (double) y->size;
}

double mgn_math_mse_matrix(gsl_matrix *y, gsl_matrix *yp)
{
    double mse = 0;
    double mse_step;
    for (size_t i = 0; i < y->size2; ++i) {
        gsl_vector_view v_row = gsl_matrix_column(y,i);
        gsl_vector_view y_rowyp = gsl_matrix_column(yp,i);
        mse_step = mgn_math_mse(&v_row.vector,&y_rowyp.vector);
        mse += mse_step;
    }
    return mse / (double) y->size2;
}
