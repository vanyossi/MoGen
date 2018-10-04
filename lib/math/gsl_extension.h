/*
 * gsl_extension.h
 *
 *  Created on: Apr 24, 2017
 *      Author: saul
 */

#ifndef MATH_GSL_EXTENSION_H_
#define MATH_GSL_EXTENSION_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


int
gsl_ran_multivariate_gaussian (const gsl_rng * r,
                               const gsl_vector * mu,
                               const gsl_matrix * L,
                               gsl_vector * result);

int
gsl_ran_multivariate_gaussian_log_pdf (const gsl_vector * x,
                                       const gsl_vector * mu,
                                       const gsl_matrix * L,
                                       double * result,
                                       gsl_vector * work);

int
gsl_ran_multivariate_gaussian_pdf (const gsl_vector * x,
                                   const gsl_vector * mu,
                                   const gsl_matrix * L,
                                   double * result,
                                   gsl_vector * work);

int
gsl_ran_multivariate_gaussian_mean (const gsl_matrix * X, gsl_vector * mu_hat);

int
gsl_ran_multivariate_gaussian_vcov (const gsl_matrix * X, gsl_matrix * sigma_hat);

/*
static int
multivar_vcov (const double data[], size_t d, size_t tda, size_t n,
               double vcov[], size_t tda2);
*/

#endif /* MATH_GSL_EXTENSION_H_ */
