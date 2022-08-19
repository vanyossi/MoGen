//
// Created by Iván Yossi on 11/05/22.
//

#ifndef MOGEN_MGN_SPHERE_H
#define MOGEN_MGN_SPHERE_H

// we probably need a mop common class
#include <stdlib.h>
#include <math.h>

typedef struct mgn_mop_param mop_param;

struct mgn_mop_param {
//    mgn_mop_param_common();
    size_t xsize;
    size_t fsize;
    size_t gsize;
};

// As in Emmerich 2004
void mgn_mop_sphere(double *x, double* f, double* g, void* param)
{
    ((void)(g)); // mark unused
    mop_param *prm = (mop_param*)param;

    f[0] = 0;
//    f[1] = 0;
    for (size_t i = 0; i < prm->xsize; ++i) {
        f[0] += x[i] * x[i];
//        f[1] += (x[i] -2.0) * (x[i] -2.0);
    }
}

int mgn_mop_sphere_min(const gsl_vector *x, const gsl_vector *v, mgn_de_ef_param* ef_p)
{
    UNUSED(ef_p);
//    double sx = 0;
//    double sv = 0;
//    size_t *size = ef_p->extra;

    // min is 0, use abs to measure min distance to origin
    double sx = gsl_vector_sum(x);
    double sv = gsl_vector_sum(v);
//    for (size_t i = 0; i < *size; ++i) {
//        sx += x[i];
//        sv += v[i];
//    }
//    printf("min %g %g:: ", fabs(sx),fabs(sv));
    return (fabs(sx) <= fabs(sv))? -1 : 1;
}

#endif //MOGEN_MGN_SPHERE_H
