//
// Created by Iván Yossi on 29/07/22.
//

#ifndef MOGEN_MGN_RBF_H
#define MOGEN_MGN_RBF_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "gsl_kmeans_new.h"
#include "mgn_cluster_m.h"

// Kernels
// r <- distance || x - c ||
void rbf_kernel_gauss(gsl_vector *r, double sigma);

void rbf_kernel_mqua(gsl_vector *r, double sigma);

void rbf_kernel_imqua(gsl_vector *r, double sigma);


gsl_matrix*
mgn_rbf_create_phi(gsl_matrix *X
                   , cluster_data *km
                   , gsl_vector *sigma
                   , void (*rbf)(gsl_vector *r, double s)
                   , gsl_matrix *ophi);


gsl_matrix *
mgn_rbf_new_weight(gsl_matrix *m_phi, gsl_matrix *y, gsl_matrix *m_w);


// helpers should probably go in another file
// caluclate mean square error
double mgn_math_mse(const gsl_vector *y, gsl_vector *yp);

// caluclate mean square error of entire matrix
double mgn_math_mse_matrix(gsl_matrix *y, gsl_matrix *yp);

#endif //MOGEN_MGN_RBF_H
