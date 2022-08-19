/*
 *
 *  SPDX-FileCopyrightText: 2022 Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include <stdlib.h>
#include <string.h>

#include "mgn_moead-rbf.h"
#include "mgn_weights.h"

#include "mgn_zdt.h"
#include "gsl_vector_additional.h"

// TODO ? sustitute by mogen.h with all types and function declarations
#include "mgn_mop.h"
#include "individual.h"
#include "population.h"
#include "mgn_poplist.h"

#include "mgn_initializer.h" // latin hyperube
#include "mgn_gnuplot.h"

#include "mops/mgn_zdt.h"
#include "mgn_rbf.h"

int main(int argc, char const *argv[]) {
    mgn_plot_open();

//    double y[3] = {4,7,-8};
//    double yp[3] = {6,10,-4};
//
//    gsl_matrix_view vy = gsl_matrix_view_array(y,1,3);
//    gsl_matrix_view vyp = gsl_matrix_view_array(yp,1,3);
//
//    printf("mse: %.6f\n", mgn_math_mse_matrix(&vyp.matrix,&vy.matrix));
//    exit(0);

    size_t maxeval = 100;

    // Problem definition
    mgn_indv_param params = {5,2,0};
    mgn_indv_ops *indv_ops = mgn_indv_ops_init();

    {
        // External referenced variables (input)
        size_t Nt = 100; //Number of m_points in tset set
        size_t N = 25;

        // non dom sol
        mgn_popl *pl_a = mgn_popl_alloc((void*)indv_ops,&params);

        gsl_matrix *m_w = mgn_weight_slattice_perm(N,pl_a->iparams.f_size);
        size_t idx[2] = {0,1};
        mgn_plot_matrix_2d(m_w,"mrbf-m_w.txt", "weights",idx);

        // Prepare Latin Hypercube// set limits
        mgnLimit *limits = mgn_limit_alloc(params.x_size);
        for (size_t i = 0; i < limits->size; ++i) {
            limits->min[i] = 0;
            limits->max[i] = 1;
        }

        mgnMoa *moead_rbf = mgn_moa_moeadrbf_alloc(1000,Nt,params.f_size*2+1,m_w,pl_a,limits);
        moead_rbf->mop = mgn_zdt_init(ZDT3,&params); // don't forget to free
        mgn_moa_moeadrbf_init(moead_rbf);

        mgn_moa_solve(moead_rbf,22);


        FILE *ofile = fopen("mrbf_test_cmean_2c.txt","w");
        mgn_popl_cursor_reset(pl_a);
        while(mgn_popl_current(pl_a) != 0) {
            mgn_indv *in = mgn_popl_next(pl_a);
            double *in_f = mgn_indv_get_params(in).f->data;
            fprintf(ofile, "%.6f %.6f\n"
                    ,in_f[0], in_f[1]);
        }
        fclose(ofile);


        gsl_matrix_free(m_w);
        mgn_mop_free(moead_rbf->mop);
        mgn_moa_moeadrbf_free(moead_rbf);
        mgn_limit_free(limits);
        mgn_popl_free(pl_a);
    }

    mgn_indv_ops_free(indv_ops);

    mgn_plot_close();
    return 0;
}
