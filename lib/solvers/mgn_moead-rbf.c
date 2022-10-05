/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_moead-rbf.h"


#include <gsl/gsl_blas.h>

#include "mgn_rbf.h"
#include "mgn_fcmeans.h"
#include "mgn_poplist.h"
#include "mgn_pop_helper.h"
//#include "mgn_gnuplot.h"

/*
 *
 * Main alg run
 *
 *
 */


// this represents one loop
// TODO remove comments, polish args + API
void mgn_moeadrbf_run(mgnMoa *moa)
{
    mgnp_moeadrbf_data *moeadrbf = mgn_moeadrbf_features(moa);

    // === Model building
    kmeans_data *km = gsl_kmeans(moeadrbf->tset->x,moeadrbf->mdl_k, 1000);
    cluster_data_extra *kme = gsl_kmeans_calc(km);
    cluster_data cdat = {km->centers, km->k};

    char* filename = malloc(sizeof(char) * 64);

    for (size_t i = 0; i < moeadrbf->mdl_size; ++i) {
        pmgn_moeadrbf_calcNN(&moeadrbf->rbf_data[i],moeadrbf->tset,&cdat,kme);

        mgnp_moeadrbf_optim_s(&moeadrbf->rbf_data[i].mdl_rbf
                              , &cdat
                              ,moeadrbf->tset
                              , moeadrbf->rbf_data[i].kernel);

        pmgn_moeadrbf_calcPhi(&moeadrbf->rbf_data[i],moeadrbf->tset,&cdat);
    }


    gsl_vector *lambda = mgnp_moeadrbf_find_lambda(moeadrbf->tset
                        , moeadrbf->rbf_data, moeadrbf->mdl_size);
//    gsl_vector_fprintf(stdout,lambda,"%.8f");
//    exit(0);

    // === Evaluate using moead
    // === Run MOEA/D -- aproximation

    //mdl training value init
    gsl_matrix *m_1phi = gsl_matrix_alloc(1,moeadrbf->mdl_k);

    moeadrbf->model_data->lhci = moeadrbf->lhci;
    moeadrbf->model_data->lambda = lambda;
    moeadrbf->model_data->km = &cdat;
    moeadrbf->model_data->mphi = m_1phi;

    mgn_ptr *m_w_ptr = malloc(sizeof(*m_w_ptr));
    mgnp_moeadrbf_mdl_optim(moeadrbf->model_data, moeadrbf->p_aprox, m_w_ptr);

#ifdef DEBUG
    mgn_plot_data pdat = {"", "", "f_1", "f_2",
                          -0.1f,1.1f,-0.1f,1.1f};
    asprintf(&pdat.title, "MOEAD_RBF");
    asprintf(&pdat.filename, "%s-test_run_i_%s_%zu", pdat.title,moa->mop->name, moa->tot_exec);
    mgn_plot((mgn_pop_proto *) moeadrbf->p_aprox, &pdat);
    free(pdat.filename);
#endif

    // === Select Points
    // mgn_popl pop_p: has all nondom aprox solutions
    // S set size of RBF cluster vectors.
//    printf("ws %zu %zu", ((gsl_matrix*)m_w_ptr->p)->size1, ((gsl_matrix*)m_w_ptr->p)->size2);
//    size_t ws_size = moeadrbf->m_w->size1; // TODO input <- really is the supplied w_matrix

    // TODO initialize inside using mgn_ptr (probably)
//    size_t Hs = (size_t)round(gsl_sf_choose(10,moeadrbf->p_aprox->iparams.f_size));
//    gsl_matrix *m_ws = mgn_weight_slattice_perm(10,((gsl_matrix*)m_w_ptr->p)->size2);

    mgn_pop *pop_sel = mgn_pop_alloc(moeadrbf->m_w->size1
        ,(void*)moeadrbf->solution->ops, &moeadrbf->solution->iparams);

    mgnp_moeadrbf_select(&moeadrbf->sel_idx
                         ,m_w_ptr->p
                         ,moeadrbf->m_w
                         ,pop_sel
                         ,moeadrbf->p_aprox);

//    printf("sel index %zu\n",moeadrbf->sel_idx.value);

    // === Update population
    mgn_mop_eval_pop(moa->mop, pop_sel, NULL);

    moa->tot_exec += pop_sel->size;

#ifdef DEBUG
    asprintf(&filename, "rbf_pop_sel-%zu.txt",moa->c_run);
    mgn_plot_fast(pop_sel, filename, "sel");
#endif

    mgnp_moeadrbf_update(moeadrbf, pop_sel);
    mgnp_moeadrbf_pop_update(moeadrbf);

#ifdef DEBUG
    printf("msize: %zu %zu,,,, %zu %zu\n", m_ws->size1, m_ws->size2
           , ((gsl_matrix*)m_w_ptr->p)->size1, ((gsl_matrix*)m_w_ptr->p)->size2);
    printf("run %zu: tot_exec %zu, idx %zu\n",moa->c_run,moa->tot_exec, moeadrbf->sel_idx.value);
#endif

    // === free all
//    gsl_matrix_free(m_ws);
    mgn_pop_free(pop_sel);
//    gsl_vector_int_free(sel_indexes);

    gsl_matrix_free(m_1phi);
    gsl_matrix_free(m_w_ptr->p);
    free(m_w_ptr);

    gsl_vector_free(lambda);
    mgn_cluster_data_extra_free(kme);
    gsl_kmeans_free(km);

#ifdef NDEBUG
    asprintf(&filename, "rbf_final-%zu.txt",moa->c_run);
    FILE *out = fopen(filename,"w");
    mgn_pop_print(moeadrbf->solution, out);
    fclose(out);
    mgn_plot_fast(moeadrbf->solution, filename, "sol");
#endif

    free(filename);
}

//mgn_pop_proto* mgn_moeadrbf_pop_get(mgnMoa *moa)
//{
//    return (mgn_pop_proto*)mgn_moeadrbf_features(moa)->solution;
//}

//===== public functions
mgnMoa* mgn_moa_moeadrbf_alloc(
    size_t max_eval
    , size_t nt
    , size_t k
    , gsl_matrix *W
    , mgn_popl *A
    , mgnLimit *limits
    , size_t Nw
    )
{
    mgnMoa *moa = mgn_moa_moeadrbf_common_alloc(max_eval
                                        ,nt
                                        ,k
                                        ,W
                                        ,A
                                        ,limits
                                        ,Nw);
    moa->run = mgn_moeadrbf_run;
    return moa;
}

void mgn_moa_moeadrbf_free(mgnMoa* moeadrbf)
{
    mgn_moa_moeadrbf_common_free(moeadrbf);
}

void mgn_moa_moeadrbf_init(mgnMoa* moeadrbf)
{
    mgn_moa_moeadrbf_common_init(moeadrbf);
}
