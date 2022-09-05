//
// Created by Iván Yossi on 18/04/22.
//

#ifndef MOGEN_MGN_SCALARIZATION_H
#define MOGEN_MGN_SCALARIZATION_H


#include "gsl_vector_additional.h"

typedef double (*mgnf_decomp_scalar)(gsl_vector*, gsl_vector*, gsl_vector*, void*);

// Functions
double mgn_scalar_tchebycheff(gsl_vector *w, gsl_vector *f, gsl_vector *z, void* param);

double mgn_scalar_pbi_ori(gsl_vector *w, gsl_vector *f, gsl_vector *z, void *d_theta);

double mgn_scalar_pbi(gsl_vector *w, gsl_vector *f, gsl_vector *z, void *d_theta);

#endif //MOGEN_MGN_SCALARIZATION_H
