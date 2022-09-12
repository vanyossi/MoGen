//
// Created by Iván Yossi on 30/04/22.
//

#include <stdlib.h>

#include "mgn_random.h"
#include "individual.h"
#include "population.h"
#include "../mgn_initializer.h"

// this test is incomplete missing validation of results.

int main() {
    rnd_initialize();
    rnd_set_seed(454235);

    mgn_indv_ops *iops = mgn_indv_ops_init();
    mgn_indv_param iparam = {6,2,0};
    mgnLimit *limits = mgn_limit_alloc(iparam.x_size);
    for (size_t i = 0; i < limits->size; ++i) {
        limits->min[i] = 0;
        limits->max[i] = 1;
    }
    mgn_pop* pop2d = mgn_pop_alloc(6,(void*)iops, &iparam);

    mgn_lhci *lhc_data = mgn_init_new_lhci(pop2d->size, iparam.x_size, limits);
//    mgn_init_LHC_init(pop2d->size,iparam.x_size, limits);

//    gsl_matrix *hlcm = mgn_LHC_get();
//    gsl_matrix_fprintf(stdout, hlcm, "%.6f");
//    gsl_vector_view vv = gsl_matrix_row(hlcm,0);
//    puts("space");
//    gsl_vector_fprintf(stdout, &vv.vector, "%g");
//    struct mgnp_init_params ipar = {limits, pop2d->ops};
//    mgn_pop_init(pop2d,mgn_init_rand,&ipar);
    mgn_pop_init(pop2d,mgn_init_lhc,lhc_data);

    for (size_t i = 0; i < pop2d->size; ++i) {
        for (size_t j = 0; j < pop2d->ops->get_iparams(pop2d->I).x->size; ++j) {
            printf("%.6f ", pop2d->ops->get_iparams(mgn_pop_get(pop2d,i)).x->data[j]);
        }
        printf("\n");
    }

    mgn_lhci_free(lhc_data);
    mgn_pop_free(pop2d);

//    double min = 100;
//    double max = 0;
//    for (int i = 0; i < 10000; ++i) {
//        double num = rnd_getUniform_limit(0.33,0.66);
//        if (num > max) { max = num;}
//        if (num < min) { min = num;}
//    }
//    printf("%.6f, %.6f\n", min, max);

    return 0;
}
