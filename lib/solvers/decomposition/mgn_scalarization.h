//
// Created by Iván Yossi on 18/04/22.
//

#ifndef MOGEN_MGN_SCALARIZATION_H
#define MOGEN_MGN_SCALARIZATION_H

#include <math.h>
#include "gsl_vector_additional.h"

double mgn_scalar_tchebycheff(gsl_vector *w, gsl_vector *f, gsl_vector *z)
{
    gsl_vector *res = gsl_vector_alloc(f->size);
    gsl_vector_memcpy(res,f);

    gsl_vector_sub(res,z);
    gsl_vector_map(res,mgn_fabs,0);
    gsl_vector_mul(res,w);

    double max = gsl_vector_max(res);

    gsl_vector_free(res);
    return max;
}
//
///**
// * Calculate the Nadir value of population
// * @param pop MoeazPop population
// * @param nadir output array of nadir values
// * @param nobj number of objectives
// */
//void mgf_pop_get_nadir(MoeazPop *pop, double *nadir, int nobj)
//{
//    unsigned int i, j;
//
//    for (i = 0; i < nobj; ++i) {
//        nadir[i] = -DBL_MAX;
//    }
//
//    for (i = 0; i < pop->size; ++i) {
//        for (j = 0; j < nobj; ++j) {
//            nadir[j] = fmax(mgf_pop_get_indv(pop,i)->f[j], nadir[j]);
//        }
//    }
//}

#endif //MOGEN_MGN_SCALARIZATION_H
