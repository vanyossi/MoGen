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

// helper functions
mgn_pop* pmgn_moead_popsel(mgnMoa *moa
                            ,mgnp_moeadrbf_data *moeadrbf
                            ,mgn_ptr *m_w_ptr
)
{
    mgn_pop *pop_sel = mgn_pop_alloc(moeadrbf->m_w->size1
                                     ,(void*)moeadrbf->solution->ops, &moeadrbf->solution->iparams);

    mgnp_moeadrbf_select(&moeadrbf->sel_idx
                         ,m_w_ptr->p
                         ,moeadrbf->m_w
                         ,pop_sel
                         ,moeadrbf->p_aprox);

    mgn_mop_eval_pop(moa->mop, pop_sel, NULL);
    moa->tot_exec += pop_sel->size;

    return pop_sel;
}

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
    char* filename = malloc(sizeof(char) * 64);
    mgnp_moeadrbf_data *moeadrbf = mgn_moeadrbf_features(moa);


    // === Model building
    kmeans_data *km = gsl_kmeans(moeadrbf->tset->x,moeadrbf->mdl_k, 1000);
    cluster_data_extra *kme = gsl_kmeans_calc(km);
    cluster_data cdat = {km->centers, km->k};


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


    // === Evaluate using moead
    // === Run MOEA/D -- aproximation
    gsl_matrix *m_1phi = gsl_matrix_alloc(1,moeadrbf->mdl_k);

    moeadrbf->model_data->search_lim = moeadrbf->search_lim;
    moeadrbf->model_data->lambda = lambda;
    moeadrbf->model_data->km = &cdat;
    moeadrbf->model_data->mphi = m_1phi;

    mgn_ptr *m_w_ptr = malloc(sizeof(*m_w_ptr));
    mgnp_moeadrbf_mdl_optim(moeadrbf->model_data, moeadrbf->p_aprox, m_w_ptr);


    // === Select Points
    // mgn_popl pop_p: has all nondom aprox solutions
    // S set size of RBF cluster vectors.
    mgn_pop *pop_sel = pmgn_moead_popsel(moa, moeadrbf, m_w_ptr);


    // === Update population
    mgnp_moeadrbf_update(moeadrbf, pop_sel);
    mgnp_moeadrbf_pop_update(moeadrbf);


    // === free all
    mgn_pop_free(pop_sel);

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
