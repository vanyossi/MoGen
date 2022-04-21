//
// Created by Iván Yossi on 19/04/22.
//

#include <stdio.h>
#include <math.h>
#include <float.h>

#include "../mgn_scalarization.h"
#include "../mgn_weights.h"
#include "mops/mgn_zdt.h"
#include "mgn_types.h"

static size_t nobj;

double ref_tch(double *f, double *w, double *z, double *n)
{
    UNUSED(n);
    size_t j;
    double tmp, tch = -DBL_MAX;
    double wj;

    for (j = 0; j < nobj; j++)
    {
        wj = (w[j] == 0.0) ? 0.0001 : w[j];
        tmp = wj * fabs(f[j] - z[j]);
        tch = (tch < tmp) ? tmp : tch;
    }
    return tch;
}

int main() {
    size_t xsize = 4;
    size_t fsize = 2;
    nobj = fsize;

    gsl_matrix *w = mgn_weight_slattice_comb(4, fsize);
//    gsl_matrix_fprintf(stdout,w, "%.6f");

    double vecres;
    double refres;

    gsl_vector *vx = gsl_vector_alloc(xsize);
    gsl_vector *vf = gsl_vector_calloc(fsize);
    gsl_vector *z = gsl_vector_calloc(fsize); // init all to 0
    for (size_t i = 0; i < xsize; ++i) {
        gsl_vector_set(vx,i,i+1);
    }
    mgn_zdt1_vector(vx,vf,0,0);

    printf("perms: %zu\n", w->size1);
    gsl_vector *vcc = gsl_vector_alloc(w->size2);
    for (size_t i = 0; i < w->size1; ++i) {
        gsl_vector_view wcur = gsl_matrix_row(w,i);
        gsl_vector_memcpy(vcc, &wcur.vector);
        gsl_ran_shuffle(rnd_get_generator(),vcc->data, vcc->size, sizeof(double));

        vecres = mgn_scalar_tchebycheff(vcc,vf,z);
        refres = ref_tch(vf->data,vcc->data,z->data,0);

        printf("(%.6f,%.6f) %.6f == %.6f\n" , vcc->data[0], vcc->data[1], vecres, refres);
    }
    gsl_vector_free(vcc);

    gsl_vector_free(vf);
    gsl_vector_free(vx);
    gsl_vector_free(z);
    gsl_matrix_free(w);
    return 0;
}
