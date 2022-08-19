//
// Created by Iván Yossi on 16/04/22.
//

#ifndef MOGEN_MGN_PARETO_H
#define MOGEN_MGN_PARETO_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "mgn_types.h"
#include "individual.h" //TODO proto

int vector_dominate(gsl_vector *u, gsl_vector *v);

int* gsl_matrix_pareto_rank(gsl_matrix *M);

bool mgn_pop_insert_dom(mgn_popl *popl, mgn_indv *indv);

#endif //MOGEN_MGN_PARETO_H
