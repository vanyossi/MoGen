//
// Created by Iván Yossi on 19/04/22.
//

#include <stdio.h>
#include <gsl/gsl_vector.h>

#include "../mgn_zdt.h"

int double_isequal(double a, double b, double delta);

int double_isequal(double a, double b, double delta){
    return (fabs(a - b) < delta);
}

int main() {

    size_t xsize = 4;
    size_t fsize = 2;
    gsl_vector *vx = gsl_vector_alloc(xsize);
    gsl_vector *vf = gsl_vector_alloc(fsize);

    double vres[2];
    double sres[2];

    for (int i = 0; i < xsize; ++i) {
        gsl_vector_set(vx,i,i+1);
    }
    mgn_zdt1_vector(vx,vf,0,0);

    for (int i = 0; i < fsize; ++i) {
        vres[i] = gsl_vector_get(vf,i);
    }
    gsl_vector_set_all(vf,0.0);

    struct mgn_zdt_param ip = {xsize, fsize, 0};
    mgn_zdt1(vx->data, vf->data,0,&ip);

    for (int i = 0; i < fsize; ++i) {
        sres[i] = gsl_vector_get(vf,i);
    }

    // check results
    // 1.0, 22.708497
    for (int i = 0; i < fsize; ++i) {
        printf("%.6f = %.6f, ", vres[i], sres[i]);
    }
    printf("\n");

    gsl_vector_free(vx);
    gsl_vector_free(vf);
    return 0;
}
