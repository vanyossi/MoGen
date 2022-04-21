//
// Created by Iván Yossi on 17/04/22.
//

#ifndef MOGEN_MGN_GEN_OPERATOR_H
#define MOGEN_MGN_GEN_OPERATOR_H

#include <gsl/gsl_vector.h>

#include "mgn_moa.h"

void mgn_double_exchange(double *a, double *b);
double mgn_pow(double base, double nexp);
void mgn_select_nrandom_el(gsl_vector_int *el, gsl_vector_int *out);
double fbounds(double val, double lower, double upper);

// Crossover
void mgn_genop_sbx(double n, double *p1, double *p2, double *c1, double *c2
                   ,size_t size
                   ,mgn_ga_sets *bounds);

void mgn_genop_sbx_alt(double eta_c, double *par1, double *par2, double *chi1, double *chi2
                       ,size_t size
                       ,mgn_ga_sets *params);

// Mutation
void mgn_genop_pbm(double n, double pm, double *p
                   ,const double *lb
                   ,const double *ub
                   ,size_t size);

#endif //MOGEN_MGN_GEN_OPERATOR_H
