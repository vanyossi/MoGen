/*
 *
 *  SPDX-FileCopyrightText: 2022 Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */


#include <stdio.h>
#include <string.h>

#include "mgn_de.h"

#include "mgn_gnuplot.h"
#include "population.h"
#include "individual.h"

#include "mgn_mop.h"
#include "mgn_initializer.h"
#include "mgn_sphere.h"
#include "mgn_zdt.h"

int main() {
    mgn_indv_param param = {4, 1, 0};
    mgn_indv_ops* iops = mgn_indv_ops_init();
    mgnLimit *rlim = mgn_limit_alloc(param.x_size);

    // set limits
    for (size_t i = 0; i < rlim->size; ++i) {
        rlim->min[i] = 0;
        rlim->max[i] = 1;
    }

    // TODO simplify on lib

    mgnMop *mop = mgn_mop_alloc();
    mop->eval_array = mgn_cast_eval(mgn_mop_sphere);
    mgnLimit *moplim = mgn_limit_alloc(param.x_size);
    for (size_t i = 0; i < moplim->size; ++i) {
        moplim->min[i] = 0;
        moplim->max[i] = 2;
    }
    mop->params = moplim;

    mgn_ga_sets ga_probs = {0.9, 0.1, NULL, NULL};
    ga_probs.mut_llim = calloc(param.x_size, sizeof(ga_probs.mut_llim));
    ga_probs.mut_ulim = calloc(param.x_size, sizeof(ga_probs.mut_ulim));
    for (size_t i = 0; i < param.x_size; ++i) {
        ga_probs.mut_llim[i] = 0;
        ga_probs.mut_ulim[i] = 1;
    }

    // Initialize and RUn DE
    size_t Np = 20;

    mgn_lhci *lhci = mgn_init_new_lhci(Np,param.x_size,rlim);
    mgnMoa* de = mgn_moa_de_alloc(Np,iops,&param,1.3, 0.5);
    mgn_de_init(de, mgn_init_lhc, lhci);

    mgn_lhci_free(lhci);

    mgn_de_ef_param ef_params = {0, &param.f_size};
    mgn_de_setmop(de, mop, mgn_cast_de_ef(mgn_mop_sphere_min), &ef_params);
    mgn_de_eval(de);


    mgn_pop *sols = (mgn_pop*)de->pop_get(de);
    for (size_t i = 0; i < sols->size; ++i) {
        mgn_indv *in = mgn_indv_get(sols,i);
        printf("%zu %.6f\n", i
               ,in->f->data[0]
               );
    }

    // mgn-solve
    mgn_plot_open();
    int runs = 30;
    int plot_every = 1;
    mgn_plot_data pdat = {"", "", "f_1", "f_2",
                          -0.1f,1.1f,-0.1f,1.1f};
    strcpy(pdat.title, "DE");

    for (int run = 0; run < runs; ++run) {
        mgn_moa_solve(de, 1); // TODO recieve solver_plot_function

        if (run % plot_every == 0) {
            asprintf(&pdat.filename, "DE_run-%d", run);
            mgn_plot((mgn_pop_proto *) sols, &pdat);
        }
    }
    printf("total exec: %zu\n", de->tot_exec);
    printf("gens: %d\n", runs);
    mgn_plot_close();

    //

    for (size_t i = 0; i < sols->size; ++i) {
       mgn_indv *in = mgn_indv_get(sols,i);
        printf("%zu %.6f\n", i
               ,in->f->data[0]
        );
    }

    mgn_moa_de_free(de);
    mgn_mop_free(mop);
    mgn_limit_free(rlim);

    return 0;
}
