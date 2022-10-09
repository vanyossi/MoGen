//
// Created by Iván Yossi on 17/04/22.
//

//#ifndef MGN_SBX_N
//#define MGN_SBX_N 30
//#endif
//
//#ifndef MGN_PBM_N
//#define MGN_PBM_N 5
//#endif

#include <string.h>

#include "mgn_moead.h"
#include "mgn_mop.h"
#include "mgn_gnuplot.h"
#include "mops/mgn_zdt.h"
#include "mops/mgn_cec09.h"

#include "individual.h"
#include "population.h"
#include "mgn_poplist.h"
#include "mgn_scalarization.h"
#include "mgn_initializer.h"
#include "mgn_weights.h"

int main() {
    mgn_plot_open();

    mgn_indv_param params = {6,2,0};
    mgn_indv_ops *iops = mgn_indv_ops_init();
    mgn_popl *EP = mgn_popl_alloc((void*)iops,&params);
//    mgnLimit ilimit = {0, 1};

//    mgnMop *mop = mgn_zdt_init(ZDT3,&params);
//    mgnMop *mop  = mgn_mop_alloc();
//    mop->eval_array = mgn_cast_eval(mgn_zdt2);
//    strcpy(mop->name, "ZDT2");
//    mop->params = &params; // same order as mgn_indv_param
    mgnLimit *moplim = mgn_limit_alloc(params.x_size);
    for (size_t i = 0; i < moplim->size; ++i) {
        moplim->min[i] = 0;
        moplim->max[i] = 1;
    }
    mgn_ga_sets ga_probs = {0.9, 0.02, NULL, NULL
                            ,5, 20};
    ga_probs.mut_llim = calloc(params.x_size, sizeof(ga_probs.mut_llim));
    ga_probs.mut_ulim = calloc(params.x_size, sizeof(ga_probs.mut_ulim));
    for (size_t i = 0; i < params.x_size; ++i) {
        ga_probs.mut_llim[i] = 0;
        ga_probs.mut_ulim[i] = 1;
    }

//    mgnMop *mop = mgn_mop_alloc();
//    mop->eval_array = mgn_cast_eval(mgn_uf1);
//    mgnLimit *moplim = mgn_limit_alloc(params.x_size);
//    mgn_cec09_set_limits(cec_uf1,moplim);
//    mop->params = &params;


//    mgn_ga_sets ga_probs = {0.9, 0.1, NULL, NULL};
//    ga_probs.mut_llim = calloc(params.x_size, sizeof(ga_probs.mut_llim));
//    ga_probs.mut_ulim = calloc(params.x_size, sizeof(ga_probs.mut_ulim));
//    for (size_t i = 0; i < params.x_size; ++i) {
//        ga_probs.mut_llim[i] = -1;
//        ga_probs.mut_ulim[i] = 1;
//    }
//

// custom
//    mgnMop *mop = mgn_mop_alloc();
//    mop->eval = mgn_cast_eval(mgn_custom);
//    mgnLimit *moplim = mgn_limit_alloc(params.x_size);
//    for (size_t i = 0; i < moplim->size; ++i) {
//        moplim->min[i] = 0;
//        moplim->max[i] = 6.3;
//    }
//    struct mycustom custom_vals = {1, 1, 0, 0};
//    mop->params = &custom_vals;
//        mgn_ga_sets ga_probs = {0.9, 0.1, NULL, NULL};
//    ga_probs.mut_llim = calloc(params.x_size, sizeof(ga_probs.mut_llim));
//    ga_probs.mut_ulim = calloc(params.x_size, sizeof(ga_probs.mut_ulim));
//    for (size_t i = 0; i < params.x_size; ++i) {
//        ga_probs.mut_llim[i] = 0;
//        ga_probs.mut_ulim[i] = 6.3;
//    }

    mgn_cec09_set_limits(UF1,moplim);
    mgnMop* mop = mgn_cec09_init(UF1, &params);


    gsl_matrix *W = mgn_weight_slattice(30, params.f_size);
    printf("m size %zu %zu\n", W->size1, W->size2);
    mgnMoa *moead = mgn_moead_init(W, 2, 20, EP, mop, mgn_init_transition,moplim,true);
    moead->max_exec = 20000;

    mgn_initializer *lhci = mgn_pinit_lhc_alloc(mgn_moead_getpop(moead),moplim);
    mgn_init_pop_lhc(mgn_moead_getpop(moead),lhci, 0);
//    mgn_moead_set_scalarization(moead, mgn_scalar_pbi);
    moead->set_ga_vals(moead,&ga_probs);

//    mgn_moead_pop_init(moead,mgn_ind_init, NULL);

    // run must be private
    int runs = 8000;
    int plot_every = 50;
    mgn_plot_data pdat = {"", "", "f_1", "f_2",
                          -0.1f,1.1f,-0.1f,1.1f};
    asprintf(&pdat.title, "%s", "points");
//    strcpy(pdat.title, "MOEAD");
    for (int run = 1; run <= runs; ++run) {
        mgn_moa_solve(moead, 1);

        if (run % plot_every == 0) {
            asprintf(&pdat.filename, "%s-%s_run-%d", pdat.title, mop->name, run);
            mgn_plot_fast((mgn_pop_proto *) EP, pdat.filename, "title");
            FILE *out = fopen("sols.txt","w");
            mgn_pop_print(EP, out);
            fclose(out);
        }
        if (moead->tot_exec >= moead->max_exec) {
            break;
        }
    }
    printf("total exec: %zu\n", moead->tot_exec);
    printf("gens: %d\n", moead->c_run);

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
//    mgn_pop* epop_array = pop_alloc_pop(EP->size,EP);
//    i = 0;
//    while(mgn_popl_current(EP) != 0) {
//        mgn_indv *in = mgn_popl_next(EP);
//        epop_array->ops->copy(mgn_pop_get(epop_array,i),in);
//        i++;
//    }
//    mgn_pop_prank_sort(epop_array);
//    for (i = 0; i < epop_array->size; ++i) {
//       mgn_indv *in = mgn_indv_get(epop_array,i);
//       printf("%zu %d %.6f %.6f %.6f %.6f %.6f %.6f\n", i, in->rank
//       ,in->x->data[0],  in->x->data[1],  in->x->data[2],  in->x->data[3]
//       ,in->f->data[0], in->f->data[1]);
//    }

    FILE *ofile = fopen("../../moead_ceUF1.txt","w");

    mgn_popl_cursor_reset(EP);
    while(mgn_popl_current(EP) != 0) {
        mgn_indv *in = mgn_popl_next(EP);
        double *in_f = mgn_indv_get_params(in).f->data;
        fprintf(ofile, "%.6f %.6f\n"
                ,in_f[0], in_f[1]);
    }
//    mgn_pop *moeadpop = mgn_moead_getfeatures(moead)->pop;
//    for (size_t i = 0; i < moeadpop->size; ++i) {
//       mgn_indv *in = mgn_indv_get(moeadpop,i);
//       fprintf(ofile, "%.6f %.6f\n"
//              ,in->f->data[0], in->f->data[1]);
//    }
    fclose(ofile);

   //    mgn_moa_solve(moead,zdt);


   //    moeadf *fe = mgn_moead_getfeatures(moead);
   //
   //    double *a = gsl_matrix_row(fe->dist,2).vector.data;
   //    gsl_vector_int_view b = gsl_matrix_int_row(fe->dindex,2);
   //
   //    double *crow = gsl_matrix_row(fe->wei,2).vector.data;
   //    for (size_t i = 0; i < 4; ++i) {
   //        int index = gsl_vector_int_get(&b.vector,i);
   //        double *r = gsl_matrix_row(fe->wei,index).vector.data;
   //        printf("v: (%.2f, %.2f) %d %f :: (%.2f, %.2f)\n", crow[0], crow[1], index, a[i], r[0], r[1]);
   //    }
    free(ga_probs.mut_llim);
    free(ga_probs.mut_ulim);

    gsl_matrix_free(W);
    mgn_limit_free(moplim);
    mgn_mop_free(mop);
    mgn_moead_free(moead);
    mgn_popl_free(EP);
    mgn_indv_ops_free(iops);

    mgn_plot_close();
    return 0;
}
