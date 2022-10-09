/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef MOGEN_MGN_POP_HELPER_H
#define MOGEN_MGN_POP_HELPER_H

#include "population.h"
#include "mgn_poplist.h"
#include "mgn_pop_matrix.h"

#include "individual.h"


void mgn_pop_copy_mp(mgn_pop_proto *pop, mgn_pop_matrix *mpop);

mgn_pop_matrix* mgn_pop_to_popm(mgn_pop_proto *pop);

void mgn_popl_insert_popm(mgn_popl *dst, mgn_pop_matrix *src);

void mgn_popl_insert_pop(mgn_popl *dst, mgn_pop_proto *src);

mgn_pop* mgn_pop_matrix_to_pop(mgn_pop_matrix *pop_m
                               ,mgn_indv_ops *ops
                               ,mgn_indv_param *params);



gsl_vector* mgn_pop_max_column(mgn_pop_proto *pop);
gsl_vector* mgn_pop_min_column(mgn_pop_proto *pop);


#endif //MOGEN_MGN_POP_HELPER_H
