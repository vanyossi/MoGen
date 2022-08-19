/*
 *
 *  SPDX-FileCopyrightText: 2022 Ivan Santa María <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */


#ifndef _MGN_VECTOR_DISTANCE_LIB_
#define _MGN_VECTOR_DISTANCE_LIB_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "gsl_vector_additional.h"

gsl_matrix* gsl_vector_distance_matrix1(gsl_vector *Vec, size_t size, double pval);

gsl_matrix* gsl_vector_distance_matrix(gsl_matrix *vecmat, double pval);

gsl_matrix* gsl_matrix_dist(gsl_matrix *m_a, gsl_matrix *m_b, double pval);

void gsl_matrix_distrank_index(gsl_matrix* m_dist, gsl_matrix_int *m_out);

#endif // __MGN_VECTOR_DISTANCE_LIB_
