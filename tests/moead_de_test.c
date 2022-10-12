/*
 *
 *  SPDX-FileCopyrightText: 2022 Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */


#include <string.h>

#include "mgn_moead_de.h"
#include "mgn_mop.h"
#include "mgn_gnuplot.h"
#include "mops/mgn_zdt.h"
#include "mops/mgn_cec09.h"

#include "individual.h"
#include "population.h"
#include "mgn_poplist.h"
#include "mgn_weights.h"


static size_t exec_limit = 1000;
void de_plot_cb(mgnMoa* moa)
{
    if (moa->tot_exec >= exec_limit) {
        mgn_plot_data pdat = {0, 0, "f_1", "f_2",
                              -0.1f,1.1f,-0.1f,1.1f};
        asprintf(&pdat.title, "%s", "points");
        asprintf(&pdat.filename, "%s-%s_run-%d", pdat.title, moa->mop->name, moa->c_run);

        mgn_pop_proto *pop = moa->pop_get(moa);
        if (pop) {
            mgn_plot_fast(pop, pdat.filename, "title");
        }

        free(pdat.filename);
        free(pdat.title);

        exec_limit += 1000;
    }
}

int main() {
    mgn_plot_open();

    mgn_indv_param params = {8,2,0};
    mgn_indv_ops *iops = mgn_indv_ops_init();
    mgn_popl *EP = mgn_popl_alloc((void*)iops,&params);

    mgn_ga_sets ga_probs = {0.9, 0.02, NULL, NULL
        ,5, 20};
    ga_probs.mut_llim = calloc(params.x_size, sizeof(ga_probs.mut_llim));
    ga_probs.mut_ulim = calloc(params.x_size, sizeof(ga_probs.mut_ulim));
    for (size_t i = 0; i < params.x_size; ++i) {
        ga_probs.mut_llim[i] = 0;
        ga_probs.mut_ulim[i] = 1;
    }

//    mgnMop *mop = mgn_zdt_init(ZDT3,&params);
    mgnMop *mop = mgn_cec09_init(UF1, &params);

    gsl_matrix *W = mgn_weight_slattice(300, 2);
    mgnMoa *moead = mgn_moead_de_init(W, 2, 20, EP, mop, mgn_init_transition,mop->limits,true);
    mgn_moa_set_callback(moead,de_plot_cb);
    moead->max_exec = 20000;

    mgn_initializer *lhci = mgn_pinit_lhc_alloc(mgn_moead_getpop(moead),mop->limits);
    mgn_init_pop_lhc(mgn_moead_getpop(moead),lhci, 0);
    mgn_pinit_free(lhci);
//    mgn_moead_set_scalarization(moead, mgn_scalar_pbi);
    moead->set_ga_vals(moead,&ga_probs);

//    mgn_moead_pop_init(moead,mgn_ind_init, NULL);

    // run must be private
    int runs = 5000;
    mgn_moa_solve(moead, runs);
//    int plot_every = 11;
//    mgn_plot_data pdat = {0, 0, "f_1", "f_2",
//                          -0.1f,1.1f,-0.1f,1.1f};
//    asprintf(&pdat.title, "%s", "points");
//
//    for (int run = 1; run <= runs; ++run) {
//        mgn_moa_solve(moead, 1);
//
//        if (run % plot_every == 0) {
//            asprintf(&pdat.filename, "%s-%s_run-%d", pdat.title, mop->name, run);
//            mgn_plot_fast((mgn_pop_proto *) EP, pdat.filename, "title");
//            FILE *out = fopen("sols.txt","w");
//            mgn_pop_print(EP, out);
//            fclose(out);
//            free(pdat.filename);
//        }
//        if (moead->tot_exec >= moead->max_exec) {
//            break;
//        }
//    }
//    free(pdat.title);
    printf("total exec: %zu\n", moead->tot_exec);
    printf("gens: %zu\n", moead->c_run);

//    mgn_pop_prank_sort(EP);
    mgn_pop *pfinal = mgn_pop_alloc(EP->size,(void*)iops,&params);
    mgn_popl_cursor_reset(EP);
    size_t i = 0;
    while(mgn_popl_current(EP) != 0) {
        mgn_indv *in = mgn_popl_next(EP);
        mgn_indv_copy(mgn_pop_get(pfinal,i), in);
        double *in_x = mgn_indv_get_params(in).x->data;
        double *in_f = mgn_indv_get_params(in).f->data;
        printf("%zu %d %.6f %.6f %.6f %.6f %.6f %.6f\n", i, mgn_indv_get_params(in).rank
               ,in_x[0], in_x[1], in_x[2], in_x[3]
               ,in_f[0], in_f[1]);
        i++;
    }

    mgn_popl_cursor_reset(EP);


    FILE *ofile = fopen("../../moeadde_uf1.txt","w");

    mgn_popl_cursor_reset(EP);
    while(mgn_popl_current(EP) != 0) {
        mgn_indv *in = mgn_popl_next(EP);
        double *in_f = mgn_indv_get_params(in).f->data;
        fprintf(ofile, "%.6f %.6f\n"
                ,in_f[0], in_f[1]);
    }

    fclose(ofile);
    mgn_pop_free(pfinal);

    mgn_moead_free(moead);


    gsl_matrix_free(W);
    mgn_mop_free(mop);
    free(ga_probs.mut_llim);
    free(ga_probs.mut_ulim);

    mgn_popl_free(EP);
    mgn_indv_ops_free(iops);

    mgn_plot_close();
    return 0;
}
