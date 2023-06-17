/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_scalarization.h"

#include <math.h>
#include <gsl/gsl_blas.h>

#include "gsl_vector_additional.h"

double mgn_scalar_tchebycheff(gsl_vector *w, gsl_vector *f, gsl_vector *z, void* param)
{
    UNUSED(param);
    gsl_vector *res = gsl_vector_alloc(f->size);
    gsl_vector_memcpy(res,f);

    gsl_vector_sub(res,z);
    gsl_vector_map(res, map_fabs, 0);
    gsl_vector_mul(res,w);

    double max = gsl_vector_max(res);

    gsl_vector_free(res);
    return max;
}

/*
*  Created on: 24/01/2011
*      Author: Saúl Zapotecas
*/
double mgn_scalar_pbi_ori(gsl_vector *w, gsl_vector *f, gsl_vector *z, void *d_theta)
{
    double d1, d2, nl;
    double theta = *((double*)d_theta);

    d1 = d2 = nl = 0.0;
    for (size_t i = 0; i < f->size; i++)
    {
        d1 += (f->data[i] - z->data[i]) * w->data[i];
        nl += pow(w->data[i], 2.0);
    }
    nl = sqrt(nl);
    d1 = fabs(d1) / nl;
    for (size_t i = 0; i < f->size; i++)
    {
        d2 += pow((f->data[i] - z->data[i]) - d1 * (w->data[i] / nl), 2.0);
    }
    d2 = sqrt(d2);

//    printf("member d1,d2: %6f,%.6f\n", d1, d2);
    return (d1 + theta * d2);
}


double mgn_scalar_pbi(gsl_vector *w, gsl_vector *f, gsl_vector *z, void *d_theta)
{
    double theta = *((double*)d_theta);

    double d1, d2, nl;

    gsl_vector *work_vec = gsl_vector_alloc(z->size);
    gsl_vector *work_vec_w = gsl_vector_alloc(w->size);
    double wnorm = gsl_blas_dnrm2(w);
    // calculate d1
    gsl_blas_dcopy(f,work_vec);
    gsl_vector_sub(work_vec,z);
    gsl_blas_ddot(work_vec, w, &nl);

    d1 = fabs(nl) / wnorm;

    // calc d2
    gsl_blas_dcopy(f,work_vec);
    gsl_blas_dcopy(w,work_vec_w);

    gsl_vector_sub(work_vec, z); // f - z
    gsl_vector_scale(work_vec_w,1.0 / wnorm); // w / ||w||
    gsl_vector_scale(work_vec_w, d1);
    gsl_vector_scale(work_vec_w, -1);

    gsl_vector_sub(work_vec, work_vec_w);
    d2 = gsl_blas_dnrm2(work_vec);

//    printf("vector d1,d2: %6f,%3.6f\n", d1, d2);

    gsl_vector_free(work_vec);
    gsl_vector_free(work_vec_w);
    return (d1 + theta * d2);
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
