//
// Created by Iván Yossi on 11/05/22.
//

#ifndef MOGEN_MGN_SELECTION_H
#define MOGEN_MGN_SELECTION_H

#include "mgn_pop_proto.h"

/*
 * Population A is substituted by B
 * We make no checks, copy every member from
 * pop B to pop A
 */
void mgn_sel_sub(mgn_pop *A, mgn_pop *B)
{
    mgn_pop_copy(A,B,0,0,A->size);
}

/*
 * Replace all from A, except the best "size" ind
 *
 * populations must be sorted
 */
void mgn_sel_lambda_mu(mgn_pop *A, mgn_pop *B,size_t size)
{
    mgn_pop *tmp_pop = mgn_pop_alloc(size,pop_get_iops(A),pop_get_iparamp(A));
    mgn_pop_copy(tmp_pop,A,0,0,size);
    mgn_pop_copy(A,B,0,0,A->size - size);
    mgn_pop_copy(A,tmp_pop,A->size - size,0,size);
}

/*
 * Select the best from both populations
 *
 * Need to supply sort function
 */
void mgn_sel_lambda_mu_plus(mgn_pop *A, mgn_pop *B, void (sort)(mgn_pop*))
{
    mgn_pop *C = mgn_pop_join(A,B);
    sort(C);

    mgn_pop_copy(A,C,0,0,A->size);

    mgn_pop_free(C);
}

#endif //MOGEN_MGN_SELECTION_H
