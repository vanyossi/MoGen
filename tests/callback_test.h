/*
 *
 *  SPDX-FileCopyrightText: 2022 Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef MOGEN_CALLBACK_TEST_H
#define MOGEN_CALLBACK_TEST_H

#include "mgn_moa.h"


void plot_callback(mgnMoa* moa)
{
    mgn_plot_data pdat = {0, 0, "f_1", "f_2",
                          -0.1f,1.1f,-0.1f,1.1f};
    asprintf(&pdat.title, "%s", "points");
    asprintf(&pdat.filename, "%s-%s_run-%zu", pdat.title, moa->mop->name, moa->c_run);

    mgn_pop_proto *pop = moa->pop_get(moa);
    if (pop) {
        mgn_plot_fast(pop, pdat.filename, "title");
    }

    free(pdat.filename);
    free(pdat.title);
}

#include "mgn_pop_helper.h"
#include "indicator/hv/hv.h"
#include "indicator/indicadores.h"

#include "mgn_pop_proto.h"
#include "mgn_pop_helper.h"
#include "mgn_io.h"

static gsl_matrix *CALLBACK_FP_M;
static size_t CALLBACK_RATE = 1000;
static size_t CALLBACK_RATE_C = 1000;
static char *CALLBACK_FILENAME = 0;
static size_t CALLBACK_RUN = 0;

void cb_set_rate(size_t rate)
{
    CALLBACK_RATE = rate;
    CALLBACK_RATE_C = rate;
}

void cp_record_perf_header()
{
    if(CALLBACK_FILENAME) {
        FILE* fname = fopen(CALLBACK_FILENAME,"w");
        fprintf(fname, "#run exec HV IGD IGD+");
        fclose(fname);
    }
}

void cb_record_perf(mgnMoa* moa)
{
    if (moa->tot_exec >= CALLBACK_RATE_C) {
        mgn_pop_proto *pop = moa->get_solutions(moa);

        // PLOT curent solutions
//        mgn_plot_data pdat = {0, 0, "f_1", "f_2",
//                              -0.1f,1.1f,-0.1f,1.1f};
//        asprintf(&pdat.title, "%s", "points");
//        asprintf(&pdat.filename, "%s-%s_run-%d", pdat.title, moa->mop->name, moa->c_run);
//
//        if (pop) {
//            mgn_plot_fast(pop, pdat.filename, "title");
//        }
//
//        free(pdat.filename);
//        free(pdat.title);
        // ------ end plot
        char *perf_stamp_file;
        asprintf(&perf_stamp_file, "sperf_%zu-%zu_%zu-%zu.log"
                 ,CALLBACK_RUN
                 ,moa->c_run
                 ,pop->iparams.x_size
                 ,pop->iparams.f_size
                 );
        FILE *out = fopen(perf_stamp_file,"w");
        fprintf(out,"# %s\n# mop: %s\n# exec: %zu\n"
                , moa->name
                , moa->mop->name
                , moa->tot_exec
                );
        mgn_pop_print(pop,out);
        fclose(out);
        free(perf_stamp_file);

        mgn_pop_matrix *popm = mgn_pop_to_popm(pop);
        // hypercube
        double *ref = calloc(3, sizeof(*ref));
        ref[0] = ref[1] = ref[2] = 1;
        double hv = fpli_hv(popm->f->data
                            ,popm->f->size2
                            ,popm->f->size1
                            , ref);
        free(ref);

        // IGD+
        double igdp = IGDplus(CALLBACK_FP_M, popm->f,2);
        double igd = IGD(CALLBACK_FP_M, popm->f,2);

        FILE* fname = fopen(CALLBACK_FILENAME,"a");
        fprintf(fname, "\n%zu %zu %e %e %e", moa->c_run, moa->tot_exec
                                            , hv, igd, igdp);
        fclose(fname);

        mgn_pop_matrix_free(popm);

        CALLBACK_RATE_C += CALLBACK_RATE;
    }
}

#endif //MOGEN_CALLBACK_TEST_H
