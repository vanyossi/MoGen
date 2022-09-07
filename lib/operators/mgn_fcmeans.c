/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_fcmeans.h"

#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>

#include "mgn_cluster_m.h"
#include "gsl_vector_additional.h"
#include "mgn_random.h"


void
mgn_fcmeans_update_c(gsl_matrix *m_x
                     , gsl_matrix *m_c
                     , gsl_matrix *m_u
                     , double m);

void
mgn_fcmeans_update_u(gsl_matrix *m_u2
                     ,gsl_matrix *m_x
                     ,gsl_matrix *m_c
                     ,gsl_matrix *m_u
                     , double m);

bool
mgn_fcmeans_error(gsl_matrix *m_u1, gsl_matrix *m_u2, double ep);

void
mgn_fcmeans_init_u(gsl_matrix *m_u, size_t k);


fcmeans_data* mgn_fcmeans(gsl_matrix *X, size_t k, size_t maxiter, double ep)
{
    double m = 2;
    fcmeans_data *km = calloc(1, sizeof(*km));
    km->m_u = gsl_matrix_alloc(X->size1, k);
    mgn_fcmeans_init_u(km->m_u,k);

//    gsl_matrix_set(km->m_u,0,0,1);
//    gsl_matrix_set(km->m_u,0,1,0);
//
//    gsl_matrix_set(km->m_u,1,0,1);
//    gsl_matrix_set(km->m_u,1,1,0);
//
//    gsl_matrix_set(km->m_u,2,0,1);
//    gsl_matrix_set(km->m_u,2,1,0);
//
//    gsl_matrix_set(km->m_u,3,0,0);
//    gsl_matrix_set(km->m_u,3,1,1);
    km->k = k;
    // initialize random centers from pop
    km->centers = gsl_matrix_alloc(k,X->size2);

    // iterate until centroids are found.
    gsl_matrix *m_u2 = gsl_matrix_alloc(km->m_u->size1, km->m_u->size2);
    size_t iter = 0;
    bool cont = 0;
    do {
        mgn_fcmeans_update_c(X, km->centers, km->m_u, m);
        mgn_fcmeans_update_u(m_u2, X, km->centers, km->m_u, m);
        cont = mgn_fcmeans_error(km->m_u, m_u2, ep);
        gsl_matrix_memcpy(km->m_u, m_u2);

    } while (++iter < maxiter && cont);
    km->iter = iter;

    gsl_matrix_memcpy(km->m_u,m_u2);
    gsl_matrix_free(m_u2);
    return km;
}

void mgn_fcmeans_free(fcmeans_data *fcmeans)
{
    gsl_matrix_free(fcmeans->m_u);
    gsl_matrix_free(fcmeans->centers);
    free(fcmeans);
}


void
mgn_fcmeans_update_c(gsl_matrix *m_x
                     , gsl_matrix *m_c
                     , gsl_matrix *m_u
                     , double m)
{
    // for each centroid
    gsl_vector *v_utmp = gsl_vector_alloc(m_u->size1);
    gsl_vector *v_xtmp = gsl_vector_alloc(m_x->size2);
    gsl_vector *v_sum = gsl_vector_calloc(m_x->size2);
    for (size_t i = 0; i < m_c->size1; ++i) {
        gsl_matrix_get_col(v_utmp,m_u,i);
        gsl_vector_map(v_utmp,map_pow,&m);

        for (size_t k = 0; k < m_x->size1; ++k) {
            gsl_matrix_get_row(v_xtmp,m_x,k);
            gsl_vector_scale(v_xtmp,gsl_vector_get(v_utmp,k));
            gsl_vector_add(v_sum,v_xtmp);
        }
        gsl_vector_scale(v_sum, 1 / gsl_vector_sum(v_utmp));
        gsl_matrix_set_row(m_c,i,v_sum);

        gsl_vector_set_all(v_sum,0);
    }
}

void
mgn_fcmeans_update_u(gsl_matrix *m_u2
                     ,gsl_matrix *m_x
                     ,gsl_matrix *m_c
                     ,gsl_matrix *m_u
                     , double m)
{
    double pval = 2 / (m - 1);
    gsl_matrix *m_rep = gsl_matrix_alloc(m_c->size1, m_x->size2);
    gsl_vector *v_dij = gsl_vector_alloc(m_u->size2);
//    gsl_vector *v_dijsum = gsl_vector_alloc(m_x->size1);
    double dij_sum;
    for (size_t i = 0; i < m_x->size1; ++i) {
        // calculate dij
        gsl_vector_view v_crow = gsl_matrix_row(m_x,i);
        gsl_vector_repeat(&v_crow.vector,m_c->size1,m_rep);
        gsl_matrix_sub(m_rep, m_c);
        for (size_t j = 0; j < m_rep->size1; ++j) {
            gsl_vector_view rep_row = gsl_matrix_row(m_rep,j);
            gsl_vector_set(v_dij,j, gsl_vector_pnorm(&rep_row.vector,2));
        }
//        gsl_vector_fprintf(stdout,v_dij,"%1.3e");
//        puts("---");
        gsl_vector *v_work = gsl_vector_alloc(v_dij->size);
        for (size_t j = 0; j < m_u->size2; ++j) {
            gsl_vector_set_all(v_work,gsl_vector_get(v_dij,j));
            gsl_vector_div(v_work,v_dij);
            gsl_vector_map(v_work,map_pow,&pval);

            dij_sum = gsl_vector_sum(v_work);
            gsl_matrix_set(m_u2,i,j, (isnan(dij_sum))? 1 :
                                     (isinf(dij_sum))? 0 :
                                     1.0 / gsl_vector_sum(v_work)
                                     );
        }
        gsl_vector_free(v_work);

    }
    gsl_matrix_free(m_rep);
    gsl_vector_free(v_dij);
    return;
}

bool
mgn_fcmeans_error(gsl_matrix *m_u1, gsl_matrix *m_u2, double ep)
{
    double norm = 0;
    gsl_matrix *m_t = gsl_matrix_alloc(m_u1->size1, m_u1->size2);
    gsl_matrix_memcpy(m_t,m_u1);
    gsl_matrix_sub(m_t,m_u2);

    norm = gsl_matrix_norm1(m_t);
    gsl_matrix_free(m_t);
//    printf("norm, %g\n", norm);
    return norm > ep;
}

void
mgn_fcmeans_init_u(gsl_matrix *m_u, size_t k)
{
    gsl_matrix_set_all(m_u,0);

    size_t psize = ceil(m_u->size1 / k);
    size_t curc = 0;
    for (size_t i = 0; i < m_u->size1; ++i) {
        gsl_matrix_set(m_u,i,curc,1.0);
        if ((i+1) % psize == 0){
//            printf("nn %d..%d..%d.....    ",i, curc, (i+1) % psize);
            curc += (curc == k-1)? 0 : 1;
        }
    }
}

gsl_vector_int*
mgn_fcmeans_get_indexes_alloc(fcmeans_data *fcm, size_t ci, double low_e, double upper_e)
{
    gsl_vector_int *v_work = gsl_vector_int_alloc(fcm->m_u->size1);
    gsl_vector_view ci_row = gsl_matrix_column(fcm->m_u,ci);

    size_t size = 0;
    for (size_t i = 0; i < ci_row.vector.size; ++i) {
        double pred = gsl_vector_get(&ci_row.vector,i);
        if ( low_e <= pred && pred <= upper_e ) {
            gsl_vector_int_set(v_work,size,i);
            size++;
        }
    }

    gsl_vector_int *v_out = gsl_vector_int_alloc(size);
    gsl_vector_int_view v_workv = gsl_vector_int_view_array(v_work->data,size);
    gsl_vector_int_memcpy(v_out, &v_workv.vector);

    gsl_vector_int_free(v_work);

    return v_out;
}

cluster_data_extra*
mgn_fcmeans_calc(fcmeans_data *fcm, double low_e, double upper_e)
{
    cluster_data_extra *kdata = mgn_cluster_data_extra_alloc(fcm->centers->size1);

    for (size_t i = 0; i < fcm->centers->size1; ++i) {
        gsl_vector_int *indexes = mgn_fcmeans_get_indexes_alloc(fcm,i,low_e,upper_e);
        kdata->mpos[i].pos = calloc(indexes->size, sizeof(unsigned int));

        for (size_t j = 0; j < indexes->size; ++j) {
            kdata->mpos[i].pos[j] = indexes->data[j];
        }
        kdata->mpos[i].size = indexes->size;
        gsl_vector_int_free(indexes); //each one is different size
    }

    return kdata;
}
