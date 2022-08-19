/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "gsl_kmeans_new.h"

#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>

#include "gsl_vector_additional.h"
#include "mgn_random.h"

gsl_matrix*
gsl_kmeans_new_centers(gsl_matrix *X, size_t k);

void
gsl_kmeans_update_index(gsl_matrix *X, gsl_vector_int *idx, gsl_matrix *C);

void
gsl_kmeans_update_c(gsl_matrix *X, gsl_matrix *C, gsl_vector_int *idx);

bool
gsl_kmeans_index_changed(gsl_vector_int *a, const gsl_vector_int *b);

kmeans_data* gsl_kmeans(gsl_matrix *X, size_t k, size_t maxiter)
{
    kmeans_data *km = calloc(1, sizeof(*km));
    km->index = gsl_vector_int_alloc(X->size1);
    km->k = k;
    // initialize random centers from pop
    km->centers = gsl_kmeans_new_centers(X,k);

    // iterate until centroids are found.
    gsl_vector_int *prev_idx = gsl_vector_int_calloc(X->size1);
    size_t iter = 0;
    do {
        gsl_kmeans_update_index(X,km->index,km->centers);
        gsl_kmeans_update_c(X,km->centers,km->index);
    } while (++iter < maxiter && gsl_kmeans_index_changed(prev_idx,km->index));
    km->iter = iter;
//    gsl_matrix_fprintf(stdout,km->centers,"::%.6f ");
//    gsl_vector_int_fprintf(stdout,km->index,"%d ");

    gsl_vector_int_free(prev_idx);
    return km;
}

void gsl_kmeans_free(kmeans_data *kmeans)
{
    gsl_vector_int_free(kmeans->index);
    gsl_matrix_free(kmeans->centers);
    free(kmeans);
}

gsl_matrix*
gsl_kmeans_new_centers(gsl_matrix *X, size_t k)
{
    gsl_matrix *C = gsl_matrix_alloc(k,X->size2);

    // index permutation
    gsl_permutation *pindex = gsl_permutation_calloc(X->size1);
    gsl_ran_shuffle(rnd_get_generator(), pindex->data, X->size1, sizeof(size_t));

    gsl_vector *vec_tmp = gsl_vector_alloc(X->size2);
    for (size_t i = 0; i < k; ++i) {
        gsl_matrix_get_row(vec_tmp, X, gsl_permutation_get(pindex, i));
        gsl_matrix_set_row(C, i, vec_tmp);
    }
    gsl_vector_free(vec_tmp);


    gsl_permutation_free(pindex);
    return C;
}

/*
 * idx index to which a sample belongs
 * C centroids
 * X samples
 */
void
gsl_kmeans_update_index(gsl_matrix *X, gsl_vector_int *idx, gsl_matrix *C)
{
    gsl_matrix *base = gsl_matrix_alloc(C->size1, X->size2);
    gsl_vector *dist = gsl_vector_alloc(C->size1);

    for (size_t s = 0; s < idx->size; ++s) {
        gsl_vector_view cs = gsl_matrix_row(X,s);
        gsl_vector_repeat(&cs.vector, C->size1,base);
        // || X_i - c_i ||^2
        gsl_matrix_sub(base, C);
        gsl_matrix_mul_elements(base,base);

        for (size_t row = 0; row < base->size1; ++row) {
            gsl_vector_view crow = gsl_matrix_row(base,row);
            gsl_vector_set(dist,row,gsl_vector_pnorm(&crow.vector,2.0));
        }
        gsl_vector_int_set(idx,s, (int) gsl_vector_min_index(dist));
    }

    gsl_matrix_free(base);
    gsl_vector_free(dist);
}

void
gsl_kmeans_update_c(gsl_matrix *X, gsl_matrix *C, gsl_vector_int *idx)
{
    gsl_matrix *mean_M = gsl_matrix_calloc(C->size1, C->size2);
    gsl_vector_int *mb = gsl_vector_int_calloc(C->size1);

    gsl_vector *cur_x = gsl_vector_alloc(C->size2);
    for (size_t i_m = 0; i_m < idx->size; ++i_m) {
        int i = gsl_vector_int_get(idx,i_m);

        gsl_matrix_get_row(cur_x, X,i_m);
        gsl_vector_view crow = gsl_matrix_row(mean_M, i);
        gsl_vector_add(&crow.vector,cur_x);
        gsl_vector_int_set(mb,i,gsl_vector_int_get(mb,i) + 1);
    }
    gsl_vector_free(cur_x);

    for (size_t cen = 0; cen < mean_M->size1; ++cen) {
        gsl_vector_view row = gsl_matrix_row(mean_M,cen);
        gsl_vector_scale(&row.vector, 1.0 / gsl_vector_int_get(mb,cen));
    }

    gsl_matrix_memcpy(C,mean_M);
    gsl_matrix_free(mean_M);
    gsl_vector_int_free(mb);
}

bool
gsl_kmeans_index_changed(gsl_vector_int *a, const gsl_vector_int *b)
{
    gsl_vector_int_sub(a,b);
    bool changed = (gsl_vector_int_sum(a) != 0);
    gsl_vector_int_memcpy(a,b);
    return changed;
}

kmeans_data_extra*
gsl_kmeans_data_extra_alloc(size_t cluster_size)
{
    kmeans_data_extra *data = calloc(1,sizeof(*data));
    data->size = cluster_size;
    data->mpos = calloc(cluster_size, sizeof(*data->mpos));

//    for (size_t i = 0; i < data->size; ++i) {
//        data->mpos[i].size = 0;
//    }

    return data;
}

void
gsl_kmeans_data_extra_free(kmeans_data_extra *data)
{
    for (size_t i = 0; i < data->size; ++i) {
        free(data->mpos[i].pos);
    }
    free(data->mpos);
    free(data);
}

kmeans_data_extra*
gsl_kmeans_calc(kmeans_data *km)
{
    kmeans_data_extra *kdata = gsl_kmeans_data_extra_alloc(km->centers->size1);

    // TODO better use a list
    gsl_matrix_uint *m_indexes = gsl_matrix_uint_alloc(km->k, km->index->size);
    size_t c_cluster;
    size_t *size;
    for (size_t i = 0; i < km->index->size; ++i) {
        c_cluster = gsl_vector_int_get(km->index,i);
        size = &kdata->mpos[c_cluster].size;
        gsl_matrix_uint_set(m_indexes,c_cluster,*size,i);
//        printf("c %zu size %zu %zu\n", c_cluster, *size, kdata->mpos[c_cluster].size);
        (*size)++;
    }

//    gsl_vector_uint *v_index;
    for (size_t i = 0; i < m_indexes->size1; ++i) {
        kdata->mpos[i].pos = calloc(kdata->mpos[i].size, sizeof(unsigned int));
        gsl_vector_uint_view v_index = gsl_matrix_uint_row(m_indexes,i);
        for (size_t j = 0; j < kdata->mpos[i].size; ++j) {
            kdata->mpos[i].pos[j] = v_index.vector.data[j];
        }
    }
    gsl_matrix_uint_free(m_indexes);
    return kdata;
}

// returns variance
// for sd square the results
gsl_matrix *
mgn_kmeans_cluster_var(kmeans_data *km, kmeans_data_extra *kme, gsl_matrix *X, bool get_sd)
{
    gsl_matrix *v_sigma = gsl_matrix_alloc(kme->size, X->size2);

    // for each cluster, find all samples from it and calc sd
    gsl_vector *v_sum = gsl_vector_calloc(X->size2);
    gsl_vector *s_row = gsl_vector_calloc(X->size2);
    gsl_vector_view centroid;
    for (size_t cluster = 0; cluster < kme->size; ++cluster) {
        gsl_vector_set_zero(v_sum);
        centroid = gsl_matrix_row(km->centers,cluster);
        for (size_t sample = 0; sample < kme->mpos[cluster].size; ++sample) {
            gsl_matrix_get_row(s_row,X,kme->mpos[cluster].pos[sample]);
            // SD
            gsl_vector_sub(s_row,&centroid.vector);
            gsl_vector_mul(s_row,s_row);
            gsl_vector_add(v_sum,s_row);
        }
        gsl_vector_scale(v_sum, 1.0 / (kme->mpos[cluster].size));
        if (get_sd) gsl_vector_map(v_sum,map_sqrt,NULL);
        gsl_matrix_set_row(v_sigma,cluster,v_sum);
    }
    gsl_vector_free(v_sum);
    gsl_vector_free(s_row);

    return v_sigma;
//    double gsl_stats_sd(const double data[], size_t stride, size_t n)
}

gsl_vector *
mgn_kmeans_cluster_var_dist(kmeans_data *km, kmeans_data_extra *kme
                       , gsl_matrix *X, bool get_sd)
{
    gsl_vector *v_sigma = gsl_vector_alloc(kme->size);

    // for each cluster, find all samples from it and calc sd
    double sum, norm;
    gsl_vector *s_row = gsl_vector_calloc(X->size2);
    gsl_vector_view centroid;
    for (size_t cluster = 0; cluster < kme->size; ++cluster) {
        sum = 0;
        centroid = gsl_matrix_row(km->centers,cluster);
        for (size_t sample = 0; sample < kme->mpos[cluster].size; ++sample) {
            gsl_matrix_get_row(s_row,X,kme->mpos[cluster].pos[sample]);
            // SD
//            double cnorm = gsl_blas_dasum(&centroid.vector);
            gsl_vector_sub(s_row,&centroid.vector);

//            norm = gsl_blas_dnrm2(s_row);
            norm = gsl_blas_dasum(s_row); // best a

//            sum += pow(norm,2.0);
            sum = sum + norm; // best a

//            sum += gsl_blas_dnrm2(s_row);
//            printf("%.3f ", sum);
        }
//        printf("\n");
        sum /= (double)kme->mpos[cluster].size;
        if (get_sd) {
            sum = sqrt(sum);
        }
        gsl_vector_set(v_sigma,cluster,sum);
    }
    gsl_vector_free(s_row);

    return v_sigma;
}
