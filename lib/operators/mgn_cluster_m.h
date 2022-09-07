/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef MOGEN_MGN_CLUSTER_M_H
#define MOGEN_MGN_CLUSTER_M_H

#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>

typedef struct pmgn_cluster_data_extra cluster_data_extra;

struct pmgn_cluster_data_idx {
    size_t size;
    unsigned int *pos;
};

struct pmgn_cluster_data_extra {
    size_t size;
    struct pmgn_cluster_data_idx *mpos;
};


cluster_data_extra* mgn_cluster_data_extra_alloc(size_t size);

void mgn_cluster_data_extra_free(cluster_data_extra *data);

//gsl_matrix *
//mgn_kmeans_cluster_var(gsl_matrix *m_centroid, cluster_data_extra *kme
//                       , gsl_matrix *X, bool get_sd);

gsl_vector *
mgn_cluster_var_dist(gsl_matrix *m_centroid, cluster_data_extra *kme
                            , gsl_matrix *X, bool get_sd);

#endif //MOGEN_MGN_CLUSTER_M_H
