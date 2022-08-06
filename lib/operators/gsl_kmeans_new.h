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

typedef struct mgn_kmeans_data kmeans_data;

struct mgn_kmeans_data {
    gsl_vector_int *index;
    gsl_matrix *centers;
};

kmeans_data* gsl_kmeans(gsl_matrix *X, size_t k, size_t maxiter);

void gsl_kmeans_free(kmeans_data *kmeans);

#endif //MOGEN_GSL_KMEANS_NEW_H
