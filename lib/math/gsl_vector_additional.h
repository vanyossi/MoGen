#ifndef _INDITEST_LIB_GSL_VECTOR_ADDITIONAL_H_
#define _INDITEST_LIB_GSL_VECTOR_ADDITIONAL_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

double mgn_fabs(double value, void* nouse);

void gsl_vector_map(gsl_vector *V, double (*func)(double, void*), void* param);

double gsl_vector_pnorm(gsl_vector *v, double pvalue);

int* gsl_vector_qsort(gsl_vector *vec);

void gsl_vector_set_seq(gsl_vector *vec);

void gsl_vector_repeat(gsl_vector *v, size_t rep, gsl_matrix *C);

#endif // _INDITEST_LIB_GSL_VECTOR_ADDITIONAL_H_
