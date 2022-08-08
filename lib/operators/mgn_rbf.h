//
// Created by Iván Yossi on 29/07/22.
//

#ifndef MOGEN_MGN_RBF_H
#define MOGEN_MGN_RBF_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "gsl_kmeans_new.h"

// Kernels
// r <- distance || x - c ||
void rbf_kernel_gauss(gsl_vector *r, double sigma);

void rbf_kernel_mqua(gsl_vector *r, double sigma);

void rbf_kernel_imqua(gsl_vector *r, double sigma);


gsl_matrix*
mgn_rbf_create_phi(gsl_matrix *X
                   , kmeans_data *km
                   , gsl_vector *sigma
                   , void (*rbf)(gsl_vector *r, double s)
                   );


gsl_matrix *
mgn_rbf_new_weight(gsl_matrix *m_phi, gsl_matrix *y);


// helpers should probably go in another file
// caluclate mean square error
double mgn_math_mse(gsl_vector *y, gsl_vector *yp);

#endif //MOGEN_MGN_RBF_H
