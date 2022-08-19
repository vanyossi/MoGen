#ifndef _INDITEST_LIB_GSL_VECTOR_ADDITIONAL_H_
#define _INDITEST_LIB_GSL_VECTOR_ADDITIONAL_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "gsl_vector_map_ops.h"

double map_fabs(double value, void* nouse);

void gsl_vector_map(gsl_vector *V, double (*func)(double, void*), void* param);

double gsl_vector_pnorm(gsl_vector *v, double pvalue);

int* gsl_vector_qsort(gsl_vector *vec);

void gsl_vector_set_seq(gsl_vector *vec);

void gsl_vector_repeat(gsl_vector *v, size_t rep, gsl_matrix *C);

gsl_vector* gsl_vector_get_indexes(gsl_vector *v, int *idx, size_t size);

gsl_matrix* gsl_matrix_get_row_indexes(gsl_matrix *v, int *idx, size_t size);

void gsl_matrix_printf(gsl_matrix *M, FILE *stream);

void gsl_matrix_save(gsl_matrix *M, char* filename);


#endif // _INDITEST_LIB_GSL_VECTOR_ADDITIONAL_H_
