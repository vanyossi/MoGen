/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef MOGEN_GSL_KMEANS_NEW_H
#define MOGEN_GSL_KMEANS_NEW_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector_int.h>
#include <stdbool.h>

#include "mgn_cluster_m.h"

typedef struct mgn_kmeans_data kmeans_data;

struct mgn_kmeans_data {
    gsl_vector_int *index;
    gsl_matrix *centers;
    size_t iter;
    size_t k;
};


kmeans_data* gsl_kmeans(gsl_matrix *X, size_t k, size_t maxiter);

void gsl_kmeans_free(kmeans_data *kmeans);

// Groups all member indexes and sizes
cluster_data_extra *
gsl_kmeans_calc(kmeans_data *km);

gsl_matrix *
mgn_kmeans_cluster_var(kmeans_data *km, cluster_data_extra *kme
                       , gsl_matrix *X, bool get_sd);

gsl_vector *
mgn_kmeans_cluster_var_dist(kmeans_data *km, cluster_data_extra *kme
                       , gsl_matrix *X, bool get_sd);


#endif //MOGEN_GSL_KMEANS_NEW_H
