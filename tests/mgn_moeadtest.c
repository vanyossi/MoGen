//
// Created by Iván Yossi on 17/04/22.
//

#include "mgn_moead.h"
#include "mops/mgn_zdt.h"

int main() {
    mgn_indv_param params = {6,2,0};
    mgn_indv_ops *iops = mgn_indv_ops_init();
    mgn_pop *EP = mgn_pop_alloc(1,(void*)iops,&params);
    mgnLimit ilimit = {0, 1};
    mgn_pop_init(EP, mgn_ind_init, &ilimit);

    mgnMop *zdt = mgn_mop_alloc();
    zdt->eval = mgn_cast_eval(mgn_zdt1_vector);
    mgn_mop_eval_pop(zdt,EP,NULL);
//
    mgnMoa *moead = mgn_moead_init(27, 2, 10, EP, zdt, mgn_ind_init,NULL);
    mgn_ga_sets ga_probs = {0.9, 0.01, NULL, NULL};
        ga_probs.mut_llim = calloc(params.realSize, sizeof(ga_probs.mut_llim));
        ga_probs.mut_ulim = calloc(params.realSize, sizeof(ga_probs.mut_ulim));
        for (size_t i = 0; i < params.realSize; ++i) {
            ga_probs.mut_llim[i] = 0;
            ga_probs.mut_ulim[i] = 1;
        }
    moead->set_ga_vals(moead,&ga_probs);

//    mgn_moead_pop_init(moead,mgn_ind_init, NULL);

    // run must be private
    mgn_moa_solve(moead,300);

    mgn_pop_prank_sort(EP);
    for (size_t i = 0; i < EP->size; ++i) {
       mgn_indv *in = mgn_indv_get(EP,i);
       printf("%zu %d %.6f %.6f %.6f %.6f %.6f %.6f\n", i, in->rank
       ,in->x->data[0],  in->x->data[1],  in->x->data[2],  in->x->data[3]
       ,in->f->data[0], in->f->data[1]);
    }

    FILE *ofile = fopen("../../moead_tmp.txt","w");
    for (size_t i = 0; i < EP->size; ++i) {
       mgn_indv *in = mgn_indv_get(EP,i);
       fprintf(ofile, "%.6f %.6f\n"
              ,in->f->data[0], in->f->data[1]);
    }
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

    mgn_mop_free(zdt);
    mgn_moead_free(moead);
    mgn_pop_free(EP);
    mgn_indv_ops_free(iops);

    return 0;
}
