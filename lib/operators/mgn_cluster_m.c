/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_cluster_m.h"

#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_blas.h>

cluster_data_extra* mgn_cluster_data_extra_alloc(size_t size)
{
    cluster_data_extra *data = calloc(1,sizeof(*data));
    data->size = size;
    data->mpos = calloc(size, sizeof(*data->mpos));

    return data;
}

void mgn_cluster_data_extra_free(cluster_data_extra *data)
{
    for (size_t i = 0; i < data->size; ++i) {
        free(data->mpos[i].pos);
    }
    free(data->mpos);
    free(data);
}

gsl_vector *
mgn_cluster_var_dist(gsl_matrix *m_centroid, cluster_data_extra *kme
                            , gsl_matrix *X, bool get_sd)
{
    gsl_vector *v_sigma = gsl_vector_alloc(kme->size);

    // for each cluster, find all samples from it and calc sd
    double sum, norm;
    gsl_vector *s_row = gsl_vector_calloc(X->size2);
    gsl_vector_view centroid;
    for (size_t cluster = 0; cluster < kme->size; ++cluster) {
        sum = 0;
        centroid = gsl_matrix_row(m_centroid,cluster);
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
