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


#endif //MOGEN_MGN_SCALARIZATION_H
