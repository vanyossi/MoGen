/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_moead_fcrbf.h"

#include <gsl/gsl_blas.h>

#include "mgn_rbf.h"
#include "mgn_fcmeans.h"
#include "mgn_pop_helper.h"
#include "mgn_gnuplot.h"

/*
 *
 * Main alg run
 *
 *
 */
void mgn_moeadfc_rbf_update_training(mgnp_moeadrbf_data* data)
{
    // set cluster size
    data->mdl_k = data->tset->size;
    // clean mphi matrix
    for (size_t i = 0; i < data->mdl_size; ++i) {
        gsl_matrix_free(data->mdl_rbf[i].m_phi);
        data->mdl_rbf[i].m_phi = NULL;
    }
}

// this represents one loop
// TODO remove comments, polish args + API
void mgn_moeadrbf_fc_run(mgnMoa *moa)
{
    mgnp_moeadrbf_data *moeadrbf = mgn_moeadrbf_features(moa);
    // === Model building
    // use archive for training
    mgn_pop_matrix_free(moeadrbf->tset);
    moeadrbf->tset = mgn_pop_to_popm((mgn_pop_proto*)moeadrbf->arc);

    fcmeans_data *fcdat = mgn_fcmeans(moeadrbf->tset->x,moeadrbf->mdl_k, 1000, 1e-5);
//    double low_e = (moa->c_run < moeadrbf->max_run * 0.6)? 0 : (moa->c_run - 1.0) / moeadrbf->max_run * 0.001;

    cluster_data_extra *cl_extra = mgn_fcmeans_calc(fcdat,0,1);
    cluster_data cdat = {fcdat->centers, fcdat->k};


    // TODO, order indexes by probability high to low
    // select only the first N solutions

//    printf("---maxiter %d\n", fcdat->iter);
//    for (size_t l = 0; l < cl_extra->size; ++l) {
//        for (size_t i1 = 0; i1 < cl_extra->mpos[l].size; ++i1) {
//            printf("%d(%0.4f),", cl_extra->mpos[l].pos[i1],
//                gsl_matrix_get(fcdat->m_u,cl_extra->mpos[l].pos[i1],l));
//        }
//        puts("--");
//    }

    char* filename = malloc(sizeof(char) * 64);

    for (size_t i = 0; i < moeadrbf->mdl_size; ++i) {

        if(moeadrbf->mdl_rbf[i].p_s) {
            gsl_vector_free(moeadrbf->mdl_rbf[i].p_s);
        }
        moeadrbf->mdl_rbf[i].p_s = mgn_cluster_var_dist(
            fcdat->centers
            ,cl_extra
            ,moeadrbf->tset->x
            ,true);
//        puts("_-___________-");

        // DEBUG
        moeadrbf->mdl_rbf[i].m_phi = mgn_rbf_create_phi(moeadrbf->tset->x,&cdat
                                                        ,moeadrbf->mdl_rbf[i].p_s
                                                        ,moeadrbf->kernel[i]
                                                        ,moeadrbf->mdl_rbf[i].m_phi);

        moeadrbf->mdl_rbf[i].p_m_w = mgn_rbf_new_weight(moeadrbf->mdl_rbf[i].m_phi
                                                        , moeadrbf->tset->f, 0);

        gsl_matrix *y_p = gsl_matrix_alloc(moeadrbf->tset->f->size1, moeadrbf->tset->f->size2);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1
                       ,moeadrbf->mdl_rbf[i].m_phi
                       ,moeadrbf->mdl_rbf[i].p_m_w,0
                       ,y_p);

#ifdef DEBUG
        asprintf(&filename, "fcmrbf_pset-%zu-%zu",moa->c_run,i);
        mgn_plot_matrix_2d(y_p,filename,"p",0);
        asprintf(&filename, "fcmrbf_pset-%zu-%zu.txt",moa->c_run,i);
        gsl_matrix_save(y_p, filename);
        // end DEBUG
#endif

        mgnp_moeadrbf_optim_s(&moeadrbf->mdl_rbf[i]
                              , &cdat
                              ,moeadrbf->tset
                              , moeadrbf->kernel[i]);

//        gsl_vector_fprintf(stdout,moeadrbf->mdl_rbf[i].p_s, "%.4f");
        moeadrbf->mdl_rbf[i].m_phi = mgn_rbf_create_phi(moeadrbf->tset->x,&cdat
                                                        ,moeadrbf->mdl_rbf[i].p_s
                                                        ,moeadrbf->kernel[i]
                                                        ,moeadrbf->mdl_rbf[i].m_phi);

        moeadrbf->mdl_rbf[i].p_m_w = mgn_rbf_new_weight(moeadrbf->mdl_rbf[i].m_phi
                                                        ,moeadrbf->tset->f,moeadrbf->mdl_rbf[i].p_m_w);

//        gsl_vector_fprintf(stdout,moeadrbf->mdl_rbf[i].p_s, "%.4f");
//        puts("_-_");

        // more DEBUG
        y_p = gsl_matrix_alloc(moeadrbf->tset->f->size1, moeadrbf->tset->f->size2);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1
                       ,moeadrbf->mdl_rbf[i].m_phi
                       ,moeadrbf->mdl_rbf[i].p_m_w,0
                       ,y_p);

        // TODO add plot multi dataset, train vs predictions
#ifdef DEBUG
        asprintf(&filename, "fcmrbf_pset-%zu-%zu_p.txt",moa->c_run,i);
        gsl_matrix_save(y_p, filename);
#endif
    }


    gsl_vector *lambda = mgnp_moeadrbf_find_lambda(moeadrbf->tset
                                                   , moeadrbf->mdl_rbf, moeadrbf->mdl_size);
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
    asprintf(&pdat.title, "MOEAD_FCRBF");
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
    mgn_popl_insert_pop(moeadrbf->arc,(mgn_pop_proto*)pop_sel);

    moa->tot_exec += pop_sel->size;

#ifdef DEBUG
    asprintf(&filename, fcrbf_pop_sel-%zu.txt",moa->c_run);
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
    mgn_moeadfc_rbf_update_training(moeadrbf);
    mgn_pop_free(pop_sel);
//    gsl_vector_int_free(sel_indexes);

    gsl_matrix_free(m_1phi);
    gsl_matrix_free(m_w_ptr->p);
    free(m_w_ptr);

    gsl_vector_free(lambda);
    mgn_cluster_data_extra_free(cl_extra);
    mgn_fcmeans_free(fcdat);

#ifdef NDEBUG
    asprintf(&filename, "fcrbf_final-%zu.txt",moa->c_run);
    FILE *out = fopen(filename,"w");
    mgn_pop_print(moeadrbf->solution, out);
    fclose(out);
    mgn_plot_fast(moeadrbf->solution, filename, "sol");
#endif

    free(filename);
}


//===== public functions
mgnMoa* mgn_moa_moeadrbf_fc_alloc(
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
                                               ,Nw
                                               );

    strncpy(moa->name, "MOEAD-FCRBF", MOA_NAME_LEN);
    moa->run = mgn_moeadrbf_fc_run;

    mgnp_moeadrbf_data *mrbf = mgn_moeadrbf_features(moa);
    mrbf->arc = mgn_popl_alloc((void*)A->ops,&A->iparams);

    return moa;
}

void mgn_moa_moeadrbf_fc_free(mgnMoa* moeadrbf)
{
    mgn_moa_moeadrbf_common_free(moeadrbf);
}

void mgn_moa_moeadrbf_fc_init(mgnMoa* moeadrbf)
{
    mgn_moa_moeadrbf_common_init(moeadrbf);

    mgnp_moeadrbf_data *mrbf = mgn_moeadrbf_features(moeadrbf);
    mgn_popl_insert_popm(mrbf->arc,mrbf->tset);

}
