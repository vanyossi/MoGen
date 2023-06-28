/*
 *
 *  SPDX-FileCopyrightText: 2022 Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "mgn_moead_de.h"
#include "mgn_weights.h"

#include "mgn_mop.h"
#include "individual.h"
#include "mgn_poplist.h"
#include "mgn_random.h"
#include "mgn_gnuplot.h"

#include "mops/mgn_cec09.h"

#include "callback_test.h"


#define MAX_EVALS 1000

#ifndef FP_DIR
#define FP_DIR "."
#endif


int main(int argc, char const *argv[]) {
#ifdef NDEBUG
    printf("Debug mode run\n");
#endif
    // default values
    size_t run = 1;
    size_t xsize = 8;
    size_t fsize = 2;
    size_t train_pop_size = 100;
    size_t wsize = 11; // weight vector size external
    size_t maxeval = MAX_EVALS;
    size_t iwsize = 100;
    char* mop_name =  malloc(sizeof(char) * 64);
    strcpy(mop_name, "UF8");

    char ch;
    while ((ch = getopt(argc, argv, "E:r:x:f:t:w:m:p:")) != -1) {
        switch (ch) {
            case 'r':
                run = strtol(optarg, NULL,10);
                break;

            case 'x':
                xsize = strtol(optarg, NULL,10);
                break;

            case 'f':
                fsize = strtol(optarg, NULL,10);
                break;

            case 't':
                train_pop_size = strtol(optarg, NULL,10);
                break;

            case 'w':
                wsize = strtol(optarg, NULL,10);
                break;

            case 'p':
                iwsize = strtol(optarg, NULL,10);
                break;

            case 'E':
                maxeval = strtol(optarg, NULL,10);
                break;

            case 'm':
                strcpy(mop_name, optarg);
                break;

            default:
                // nothing
                break;
        }
    }
    argc -= optind;
    argv += optind;

//    xsize = 4;
//    fsize = 2;
//    wsize = 10;
//    iwsize = 300;
//    maxeval = 300;
//    train_pop_size = 100;
//    wsize = (size_t)round(gsl_sf_choose(60,fsize));
//    wsize = 20;

    mgn_plot_open();
    CALLBACK_RUN = run;

    // Problem definition
    mgn_indv_param params = {xsize,fsize,0};
    mgn_indv_ops *indv_ops = mgn_indv_ops_init();
    mgn_popl *EP = mgn_popl_alloc((void*)indv_ops,&params);

    mgn_ga_sets ga_probs = {0.9, 0.02, NULL, NULL
        ,5, 20};
    ga_probs.mut_llim = calloc(params.x_size, sizeof(ga_probs.mut_llim));
    ga_probs.mut_ulim = calloc(params.x_size, sizeof(ga_probs.mut_ulim));
    for (size_t i = 0; i < params.x_size; ++i) {
        ga_probs.mut_llim[i] = 0;
        ga_probs.mut_ulim[i] = 1;
    }

    {
        size_t N = 12; // 3 obj
//        size_t N = 99; // weight vectors define exec for each run
        gsl_matrix *m_w = mgn_weight_slattice(N,EP->iparams.f_size);
        printf("pop size: %zu\n", m_w->size1);
        printf("run %zu\n", run);

        MGN_CEC09_VAR moptype = mop_cec09_str_toenum(mop_name);
        mgnMop *mop = mgn_cec09_init(moptype, &params);


        mgnMoa *moead = mgn_moead_de_init(m_w
            , 2, 20, EP, mop, mgn_init_transition,mop->limits,true);

        moead->set_ga_vals(moead,&ga_probs);
        moead->max_exec = maxeval;

        mgn_initializer *lhci = mgn_pinit_lhc_alloc(mgn_moead_getpop(moead),mop->limits);
        mgn_init_pop_lhc(mgn_moead_getpop(moead),lhci, 0);
        mgn_pinit_free(lhci);


        // plot callback setup
        cb_set_rate(1);
        mgn_moa_set_callback(moead, cb_record_perf);
        asprintf(&CALLBACK_FILENAME, "%s_%s_%zu-%zu_%zu-%zu-perf.txt"
                 , moead->name
                 , mop->name
                 , xsize
                 , fsize
                 , maxeval
                 , run);
        cp_record_perf_header();

        struct _inGroup_list io_data_fp = {0, 0};//malloc(sizeof(*io_data_fp));
        char *file_front;

        asprintf(&file_front, "%s/cec%s.txt",FP_DIR,mop->name);
        it_read_data(file_front,&io_data_fp); // alters size of pl_a
        CALLBACK_FP_M = inData_toGSLMatrix(inGroup_getListAt(&io_data_fp,0));
        free(file_front);
        // callback setup end


        int runs = 5000;
        mgn_moa_solve(moead, runs);


        // print results
        char tmp_fname[64];
        char* filename;
        // name, alg, var, obj, run
        asprintf(&filename, "%s_%s-%zu-%zu_%zu-%zu"
                 ,moead->name
                 ,moead->mop->name
                 ,xsize
                 ,fsize
                 ,maxeval
                 ,run
        );

        strcpy(tmp_fname,filename);
        strcat(tmp_fname,".txt");
        FILE *out = fopen(tmp_fname,"w");
        if (!out) {
            printf("fname : %s", tmp_fname);
            perror("fopen");
        } else {
            mgn_pop_print(moead->get_solutions(moead), out);
            fclose(out);
        }
        mgn_plot_fast(moead->get_solutions(moead), filename, "solutions");
        free(filename);


        mgn_mop_free(mop);
        mgn_moead_free(moead);

        gsl_matrix_free(m_w);
    }
    free(ga_probs.mut_llim);
    free(ga_probs.mut_ulim);

    mgn_popl_free(EP);

    mgn_indv_ops_free(indv_ops);
    free(mop_name);

    mgn_plot_close();
    free(CALLBACK_FILENAME);
    return 0;
}
