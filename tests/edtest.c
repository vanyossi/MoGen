//
// Created by Iván Yossi on 10/05/22.
//

#include <stdio.h>

#include "mgn_de.h"

#include "population.h"
#include "individual.h"

#include "mgn_initializer.h"
#include "mgn_sphere.h"

int main() {
    mgn_indv_param param = {4, 1, 0};
    mgn_indv_ops* iops = mgn_indv_ops_init();
    mgnLimit *rlim = mgn_limit_alloc(param.realSize);

    // set limits
    for (size_t i = 0; i < rlim->size; ++i) {
        rlim->min[i] = 0;
        rlim->max[i] = 1;
    }

    // TODO simplify on lib
    mgnMop *mop = mgn_mop_alloc();
    mop->eval_array = mgn_cast_eval(mgn_mop_sphere);
    mgnLimit *moplim = mgn_limit_alloc(param.realSize);
    for (size_t i = 0; i < moplim->size; ++i) {
        moplim->min[i] = 0;
        moplim->max[i] = 2;
    }
    mop->params = moplim;
    mgn_ga_sets ga_probs = {0.9, 0.1, NULL, NULL};
    ga_probs.mut_llim = calloc(param.realSize, sizeof(ga_probs.mut_llim));
    ga_probs.mut_ulim = calloc(param.realSize, sizeof(ga_probs.mut_ulim));
    for (size_t i = 0; i < param.realSize; ++i) {
        ga_probs.mut_llim[i] = 0;
        ga_probs.mut_ulim[i] = 1;
    }


    size_t Np = 20;
    de_param* de_p = mgn_de_alloc(Np,iops,&param);
    mgn_init_LHC_init(Np,param.realSize,rlim);
    mgnMoa *de = mgn_de_init(de_p, mgn_init_lhc, iops);
    mgn_de_setmop(de, mop);
    mgn_de_eval(de);


    for (size_t i = 0; i < de_p->pop->size; ++i) {
        mgn_indv *in = mgn_indv_get(de_p->pop,i);
        printf("%zu %.6f\n", i
               ,in->f->data[0]
               );
    }

    mgn_moa_solve(de,50);

    for (size_t i = 0; i < de_p->pop->size; ++i) {
       mgn_indv *in = mgn_indv_get(de_p->pop,i);
        printf("%zu %.6f\n", i
               ,in->f->data[0]
        );
    }

    mgn_de_free(de);
    mgn_mop_free(mop);
    mgn_limit_free(rlim);

    return 0;
}
