/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef MOGEN_MGN_FCMEANS_H
#define MOGEN_MGN_FCMEANS_H

#include <stdbool.h>

#include <gsl/gsl_matrix.h>
#include "mgn_cluster_m.h"
//#include <gsl/gsl_vector_int.h>

typedef struct mgn_fcmeans_data fcmeans_data;

struct mgn_fcmeans_data {
    gsl_matrix *m_u;
    gsl_matrix *centers;
    size_t iter;
    size_t k;
    double alpha;
    double epsilon;
};

//
//// after fcmeans calc, this can be used to
//// analyze the data members.
//// TODO contain the m_points themselves as matrix
////      tricky as each row could be different size
//struct mgn_fcmeans_data_idx {
//    size_t size;
//    unsigned int *pos;
//};
//
//struct mgn_fcmeans_data_extra {
//    size_t size;
//    struct mgn_fcmeans_data_idx *mpos;
//};

fcmeans_data* mgn_fcmeans(gsl_matrix *X, size_t k, size_t maxiter, double ep);

void mgn_fcmeans_free(fcmeans_data *fcmeans);


cluster_data_extra*
mgn_fcmeans_calc(fcmeans_data *fcm, double low_e, double upper_e);

//void gsl_fcmeans_data_extra_free(fcmeans_data_extra *data);

// Groups all member indexes and sizes
//fcmeans_data_extra*
//gsl_fcmeans_calc(fcmeans_data *km);

//gsl_matrix *
//mgn_fcmeans_cluster_var(fcmeans_data *km, fcmeans_data_extra *kme
//                       , gsl_matrix *X, bool get_sd);
//
//gsl_vector *
//mgn_fcmeans_cluster_var_dist(fcmeans_data *km, fcmeans_data_extra *kme
//                            , gsl_matrix *X, bool get_sd);



#endif //MOGEN_MGN_FCMEANS_H
