/*
 *
 *  SPDX-FileCopyrightText: 2022 Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include <gsl/gsl_sf.h>
#include "mgn_moead_fcrbf_nn.h"
#include "mgn_weights.h"
#include "gsl_vector_additional.h"

// TODO ? sustitute by mogen.h with all types and function declarations
#include "mgn_mop.h"
#include "individual.h"
#include "mgn_poplist.h"
#include "mgn_random.h"
#include "mgn_gnuplot.h"

#include "mops/mgn_cec09.h"
#include "mops/mgn_zdt.h"

#include "callback_test.h"

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
    size_t maxeval = train_pop_size + wsize;
    size_t iwsize = 100;
    char* mop_name =  malloc(sizeof(char) * 64);
    strcpy(mop_name, "UF1");

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


    // Problem definition
    mgn_indv_param params = {xsize,fsize,0};
    mgn_indv_ops *indv_ops = mgn_indv_ops_init();

    {
        // External referenced variables (input)
        size_t Nt = train_pop_size; //Number of m_points in tset set
        size_t N = wsize;
        printf("%zu %zu %zu\n", maxeval, Nt, N);
        size_t total_runs = (maxeval - Nt) / N;

        // non dom sol
        mgn_popl *pl_a = mgn_popl_alloc((void*)indv_ops,&params);
        gsl_matrix *m_w = mgn_weight_slattice(N,pl_a->iparams.f_size);
        printf("weight size %zu\n", m_w->size1);

        // Prepare Latin Hypercube// set limits
        mgnLimit *limits = mgn_limit_alloc(params.x_size);
        for (size_t i = 0; i < limits->size; ++i) {
            limits->min[i] = 0;
            limits->max[i] = 1;
        }

        mgnMoa *moead_fcrbf = mgn_moa_moeadrbf_fcnn_alloc(maxeval
                                                   ,Nt
                                                   ,params.f_size*4+1
                                                   ,m_w
                                                   ,pl_a
                                                   ,limits
                                                   ,iwsize);

        moead_fcrbf->max_exec = maxeval;
//        mgn_moa_set_callback(moead_fcrbf, plot_callback);

//        MGN_ZDT_VAR moptype = mop_zdt_str_toenum(mop_name);
//        mgnMop* mop = mgn_zdt_init(moptype, &params);

        MGN_CEC09_VAR moptype = mop_cec09_str_toenum(mop_name);
        mgnMop* mop = mgn_cec09_init(moptype, &params);

        moead_fcrbf->mop = mop;

        // plot callback setup
        cb_set_rate(1);
        mgn_moa_set_callback(moead_fcrbf, cb_record_perf);
        asprintf(&CALLBACK_FILENAME, "%s_%s_%zu-%zu_%zu-%zu-perf.txt"
                 , moead_fcrbf->name
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

        mgn_moa_moeadrbf_fcnn_init(moead_fcrbf);

        printf("expected total runs %zu\n", total_runs);
        mgn_moa_solve(moead_fcrbf,total_runs);

        // print results
        char tmp_fname[64];
        char* filename;
        // name, alg, var, obj, run
        asprintf(&filename, "%s_%s-%zu-%zu_%zu-%zu"
                 ,moead_fcrbf->name
                 ,moead_fcrbf->mop->name
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
            mgn_pop_print(moead_fcrbf->pop_get(moead_fcrbf), out);
            fclose(out);
        }
        mgn_plot_fast(moead_fcrbf->pop_get(moead_fcrbf), filename, "solutions");
        free(filename);


        mgn_moa_moeadrbf_fcnn_free(moead_fcrbf);
        mgn_mop_free(mop);
        gsl_matrix_free(m_w);
        mgn_limit_free(limits);
        mgn_popl_free(pl_a);
    }

    mgn_indv_ops_free(indv_ops);
    free(mop_name);

    mgn_plot_close();
    free(CALLBACK_FILENAME);
    return 0;
}
