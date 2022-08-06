//
// Created by Iván Yossi on 29/07/22.
//

#ifndef MOGEN_MGN_RBF_H
#define MOGEN_MGN_RBF_H

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "gsl_kmeans_new.h"
#include "gsl_vector_additional.h"

// TODO move to a gsl_map header helper macro
double pgsl_vector_map_exp(double val, void* param)
{
    return exp(val);
}
// Kernels
// d <- distance || x - y ||
void rbf_kernel_gauss(gsl_vector *d, gsl_vector *sigma)
{
    // r^2 / sigma
    gsl_vector_mul(d, d);
    gsl_vector_div(d,sigma);
    gsl_vector_map(d,pgsl_vector_map_exp,NULL);
}

void mgn_rbf_train(gsl_matrix *P, gsl_vector Y);



#endif //MOGEN_MGN_RBF_H
