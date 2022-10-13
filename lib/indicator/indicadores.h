#ifndef _IT_OPS_INDICADORES_
#define _IT_OPS_INDICADORES_

#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

double euclidian(gsl_vector *A, gsl_vector *B, double p);
double euclidianMax(gsl_vector *Z, gsl_vector *A, double p);

double GD(gsl_matrix *FP, gsl_matrix *S, double p);
double GDp(gsl_matrix *FP, gsl_matrix *S, double p);
double GDplus(gsl_matrix *FP, gsl_matrix *S, double p);

double IGD(gsl_matrix *FP, gsl_matrix *S, double p);
double IGDp(gsl_matrix *FP, gsl_matrix *S, double p);
double IGDplus(gsl_matrix *FP, gsl_matrix *S, double p);

double deltap(gsl_matrix *FP, gsl_matrix *S, double p);

bool uDomv(gsl_vector *u, gsl_vector *v);
double inSetTwoCover(gsl_matrix *A, gsl_matrix *B);

gsl_vector* inParetoRank(gsl_matrix *A, gsl_matrix *B);

#endif // _IT_OPS_INDICADORES_
