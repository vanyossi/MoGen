/*
 *
 *  SPDX-FileCopyrightText: 2022 Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_moead-rbf_common.h"

#include <string.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>
#include "gsl_vector_additional.h"

#include "mgn_rbf.h"
#include "mgn_de.h"
#include "mgn_moead_de.h"
#include "population.h"
#include "mgn_poplist.h"
#include "mgn_pop_matrix.h"
#include "mgn_pop_helper.h"

#include "mgn_weights.h"
#include "mgn_pareto.h"
#include "mgn_vector_distance.h"
#include "mgn_counter.h"
#include "mgn_gnuplot.h"


mgn_pop_proto* mgn_moeadrbf_pop_get(mgnMoa *moa)
{
    return (mgn_pop_proto*)mgn_moeadrbf_features(moa)->solution;
}

mgnMoa* mgn_moa_moeadrbf_common_alloc(
    size_t max_eval
    , size_t nt
    , size_t k
    , gsl_matrix *W
    , mgn_popl *A
    , mgnLimit *limits
    , size_t Nw /* size of internal weight vector */
    )
{
    mgnMoa *moa = malloc(sizeof(*moa));
    mgnp_moeadrbf_data *mrbf = malloc(sizeof(*mrbf));

    size_t x_size = A->iparams.x_size;
    size_t f_size = A->iparams.f_size;
    size_t g_size = A->iparams.g_size;

    mrbf->max_eval = max_eval;
    mrbf->nt = nt;
    mrbf->m_w = W;
    mrbf->solution = A;
    // p_aprox size is not important as it gets exchanged by internal pop of moead
    mrbf->p_aprox = mgn_pop_alloc(1,A->ops,&A->iparams);
    mrbf->tset = mgn_pop_matrix_alloc(nt,x_size
                                      ,nt,f_size
                                      ,nt,g_size);

    // model training initilize
    mrbf->mdl_size = 3;
    mrbf->mdl_k = k;
    mrbf->mdl_wsize = Nw;
//    mrbf->kernel = malloc(sizeof(mrbf->kernel) * mrbf->mdl_size);
    mrbf->rbf_data = calloc(mrbf->mdl_size, sizeof(*mrbf->model_data));
    mrbf->rbf_data[0].kernel = rbf_kernel_gauss;
    mrbf->rbf_data[1].kernel = rbf_kernel_mqua;
    mrbf->rbf_data[2].kernel = rbf_kernel_imqua;

    for (size_t i = 0; i < mrbf->mdl_size; ++i) {
        mrbf->rbf_data[i].mdl_rbf.p_s = NULL;
        mrbf->rbf_data[i].mdl_rbf.p_m_w = NULL;
        mrbf->rbf_data[i].mdl_rbf.m_phi = NULL;
    }
    // used to keep control of selected Neigbourh
    mrbf->sel_idx = (mgn_count_ciclic){10,0}; // selection neighbourhood
    mgn_init_lhc_to_matrix(mrbf->tset->x, limits);
    mrbf->search_lim = mgn_limit_alloc(limits->size);
    mgn_limit_cpy(mrbf->search_lim, limits);

    // initialize tset population
    strncpy(moa->name, "MOEAD-RBF", MOA_NAME_LEN);
    moa->tot_exec = 0;
    moa->c_run = 0;
//    moa->run = mgn_moeadrbf_run;
    moa->stop = mgnp_moeadrbf_stop;
    moa->set_ga_vals = mgnp_moeadrbf_set_ga_vals;
    moa->features = mrbf;
    moa->pop_get = mgn_moeadrbf_pop_get;

    // initialize model parameters
    mrbf->model_data = calloc(1, sizeof(*mrbf->model_data));
    mrbf->model_data->iW = mgn_weight_slattice(Nw, f_size);
    mrbf->model_data->mop_limits = limits;
    mgnp_moeadrbf_mdl_defparam(mrbf,Nw,0.9,0.02);

    return moa;
}

void mgn_moa_moeadrbf_common_free(mgnMoa* moeadrbf)
{
    mgnp_moeadrbf_data *mrbf = mgn_moeadrbf_features(moeadrbf);

    mgn_pop_matrix_free(mrbf->tset);
    mgn_pop_free(mrbf->p_aprox);

    for (size_t i = 0; i < mrbf->mdl_size; ++i) {
        gsl_vector_free(mrbf->rbf_data[i].mdl_rbf.p_s);
        gsl_matrix_free(mrbf->rbf_data[i].mdl_rbf.p_m_w);
        gsl_matrix_free(mrbf->rbf_data[i].mdl_rbf.m_phi);
    }
    free(mrbf->rbf_data);
    mgn_limit_free(mrbf->search_lim);

    mgnp_moeadrbf_mdl_mop_param_free(mrbf);

    free(mrbf);
    free(moeadrbf);
}


void mgn_moa_moeadrbf_common_init(mgnMoa* moeadrbf)
{
    mgnp_moeadrbf_data *mrbf = mgn_moeadrbf_features(moeadrbf);
    mgn_pop_matrix_eval(mrbf->tset,moeadrbf->mop);
    moeadrbf->tot_exec += mrbf->tset->x->size1;


#ifdef DEBUG
    mgn_plot_matrix_2d(mrbf->tset->f,"mrbf_tset_init", "f",0);
    gsl_matrix_save(mrbf->tset->f, "mrbf_tset_init_p.txt");
#endif

//    gsl_matrix_printf(mrbf->tset->x,stdout);
//    printf("max min, %.6f %.6f\n", gsl_matrix_min(mrbf->tset->x), gsl_matrix_max(mrbf->tset->x));

    mgn_pop *pop = mgn_pop_alloc(mrbf->tset->x->size1, mrbf->solution->ops, &mrbf->solution->iparams);
    mgn_pop_copy_mp((mgn_pop_proto*)pop,mrbf->tset);
    for (size_t i = 0; i < pop->size; ++i) {
        mgn_pop_insert_dom(mrbf->solution,mgn_pop_get(pop,i));
    }
    mgn_pop_free(pop);

    // TODO print pop
    //      done? mgn_pop_print(mrbf->solution, stdout);
//    for (size_t i = 0; i < mrbf->solution->size; ++i) {
//        gsl_vector *f = ((mgn_indv*) mrbf->solution->get(mrbf->solution,i))->f;
//        printf("%.3f, %.3f", gsl_vector_get(f,0), gsl_vector_get(f,1));
//        puts("");
//    }

//    printf("size of psols %zu\n", mrbf->solution->size);
//    mgn_plot_open();
//    mgn_plot_data pdat = {"", "", "f_1", "f_2",
//                          -0.1f,1.1f,-0.1f,1.1f};
//    strcpy(pdat.title, "MOEAD_RBF");
//    asprintf(&pdat.filename, "%s-test_run_solution", pdat.title);
//    mgn_plot((mgn_pop_proto *) mrbf->solution, &pdat);
//    free(pdat.filename);
//    mgn_plot_close();
}


mgnp_moeadrbf_data* mgn_moeadrbf_features(mgnMoa* moa)
{
    return (mgnp_moeadrbf_data*)moa->features;
}

bool mgnp_moeadrbf_stop(mgnMoa *moa)
{
    mgnp_moeadrbf_data *adata = mgn_moeadrbf_features(moa);
    return (moa->tot_exec >= adata->max_eval);
}

void mgnp_moeadrbf_set_ga_vals(mgnMoa *moa, mgn_ga_sets *ga)
{
    UNUSED(moa);
    UNUSED(ga);
}

int mgnp_moeadrbf_s_optim_mop_min(gsl_vector *f1
                                  ,gsl_vector *f2
                                  ,mgn_de_ef_param *ef_prm) {
    UNUSED(ef_prm);
//    struct mgnp_moeadrbf_s_optim_ef_p *prm = (struct mgnp_moeadrbf_s_optim_ef_p*)ef_prm;
//
//    gsl_vector_view v_yt = gsl_matrix_row(prm->y_t,ef_prm->pos);
//
//    double mse1 = mgn_math_mse(&v_yt.vector,f1);
//    double mse2 = mgn_math_mse(&v_yt.vector,f2);

    double mse1 = gsl_vector_get(f1,0);
    double mse2 = gsl_vector_get(f2,0);

//    printf("================ mse: %.4f, %.4f\n", mse_f1, mse_f2);

    return ((mse1) <= mse2)? -1 : 1;
}

void mgnp_moeadrbf_s_optim_mop(gsl_vector *x
                               ,gsl_vector *f
                               ,gsl_vector *g
                               ,struct mgnp_moeadrbf_s_optim_p *prm
)
{
    UNUSED(g);
//    struct mgnp_moeadrbf_s_optim_p *prm = (struct mgnp_moeadrbf_s_optim_p*)m_prm;

    mgn_rbf_create_phi(prm->tset->x,prm->km,x,prm->rbf,prm->m_phi);
    mgn_rbf_new_weight(prm->m_phi,prm->tset->f,prm->m_w);


    if(!isnan(gsl_matrix_get(prm->m_w,0,0))) {
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,prm->m_phi,prm->m_w,0,prm->f_p);
        gsl_vector_set(f, 0,mgn_math_mse_matrix(prm->tset->f,prm->f_p));
    } else {
        gsl_vector_set(f,0,1e10);
    }
}

void mgnp_moeadrbf_optim_s(
    struct mgnp_rbf_weigts *vals
    , cluster_data *km
    , mgn_pop_matrix *tset
    , mgn_kernel_f rbf
)
{
    mgn_indv_param param = {vals->p_s->size, 1, 0};
    mgn_indv_ops* iops = mgn_indv_ops_init();

    // set limits
    double s_max = gsl_vector_max(vals->p_s);
    double s_min = gsl_vector_min(vals->p_s);
    mgnLimit *rlim = mgn_limit_alloc(param.x_size);
    for (size_t i = 0; i < rlim->size; ++i) {
        rlim->min[i] = s_min;
        rlim->max[i] = s_max;
    }

    mgnMop *mop = mgn_mop_alloc(&param);
    mop->eval = mgn_cast_eval(mgnp_moeadrbf_s_optim_mop);
    gsl_matrix *tmp_f_p = gsl_matrix_alloc(tset->f->size1, tset->f->size2);
    gsl_matrix *tmp_w_p = gsl_matrix_alloc(km->k, tset->f->size2);
    gsl_matrix *tmp_m_phi = gsl_matrix_alloc(tset->x->size1, km->k);
    struct mgnp_moeadrbf_s_optim_p mop_p = {0, tset, km, tmp_f_p, tmp_w_p, tmp_m_phi, rbf};
    mop->params = &mop_p;

    size_t Np = 20; // de pop size

    mgnMoa* de = mgn_moa_de_alloc(Np,iops,&param,1, 0.5);

    mgn_initializer *lhci = mgn_pinit_lhc_alloc(de->pop_get(de),rlim);
    mgn_init_pop_lhc(de->pop_get(de),lhci,0);
    mgn_pinit_free(lhci);

    de->max_exec = SIZE_MAX;

    struct mgnp_moeadrbf_s_optim_ef_p *ef_prm = malloc(sizeof(*ef_prm));
    ef_prm->y_t = tset->f;
//        ef_prm->y_cur = gsl_vector_alloc(tset->f->size2);
    mgn_de_setmop(de, mop, mgn_cast_de_ef(mgnp_moeadrbf_s_optim_mop_min),(mgn_de_ef_param*)ef_prm);
    mgn_de_eval(de);

    mgn_moa_solve(de, 50);

    mgn_pop *sols = (mgn_pop*)de->pop_get(de); // fisrt one is the lowest
    for (size_t i = 0; i < 1; ++i) {
        mgn_indv *in = mgn_indv_get(sols,i);
        gsl_vector_memcpy(vals->p_s,in->x);
//        printf("%zu %.6f\n", i
//               ,in->f->data[0]
//        );
    }
//    puts("=====");
    free(ef_prm);
    mgn_moa_de_free(de);
    mgn_mop_free(mop);
    mgn_limit_free(rlim);

    gsl_matrix_free(tmp_m_phi);
    gsl_matrix_free(tmp_f_p);
    gsl_matrix_free(tmp_w_p);

    mgn_indv_ops_free(iops);
}

gsl_vector *mgnp_moeadrbf_find_lambda(mgn_pop_matrix *tset, mgn_model_rbfdata *res, size_t size)
{
    gsl_vector *lambda = gsl_vector_calloc(size);

    gsl_matrix *f_vals = gsl_matrix_alloc(tset->f->size1,size);
    gsl_matrix *tmp_y_p = gsl_matrix_alloc(tset->f->size1, tset->f->size2);

    for (size_t i = 0; i < size; ++i) {
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1
                       ,res[i].mdl_rbf.m_phi,res[i].mdl_rbf.p_m_w,0,tmp_y_p);
        for (size_t j = 0; j < tmp_y_p->size1; ++j) {
            gsl_vector_view y_row = gsl_matrix_row(tset->f,j);
            gsl_vector_view yp_row = gsl_matrix_row(tmp_y_p,j);
            gsl_matrix_set(f_vals,j,i,mgn_math_mse(&y_row.vector,&yp_row.vector));
        }
    }
//    gsl_matrix_printf(f_vals,stdout);
    for (size_t i = 0; i < f_vals->size1; ++i) {
        gsl_vector_view c_row = gsl_matrix_row(f_vals,i);
        lambda->data[gsl_vector_min_index(&c_row.vector)]++;
    }
    gsl_vector_scale(lambda, 1.0 / tset->x->size1);
//    gsl_vector_fprintf(stdout,lambda,"%.6f");
    gsl_matrix_free(f_vals);
    gsl_matrix_free(tmp_y_p);
    return lambda;
}

void pmgn_moeadrbf_lhci_update(mgnp_moeadrbf_data *rbfdata, mgn_pop_matrix *popm)
{
//    mgn_pinit_free(rbfdata->lhci);
//    mgn_pop* pop = mgn_pop_matrix_to_pop(popm
//                        , rbfdata->p_aprox->ops
//                        , &rbfdata->p_aprox->iparams);
//    rbfdata->lhci = mgn_pinit_lhc_alloc(pop,0);
//
//    mgn_pop_free(pop);
}


//void pmgn_moeadrbf_lchi_update(mgn_lhci **lhci, gsl_matrix *x, size_t size)
//{
//    mgnLimit *limits = mgn_limit_alloc(x->size2);
//    for (size_t i = 0; i < limits->size; ++i) {
//        gsl_vector_view c_col = gsl_matrix_column(x,i);
//        double r_mean = gsl_stats_mean(c_col.vector.data,1,c_col.vector.size);
//        limits->min[i] = fmax(gsl_vector_min(&c_col.vector) - r_mean,0);
//        limits->max[i] = fmin(gsl_vector_max(&c_col.vector) + r_mean,1);
//    }
//
////    printf("refreshing lhci\n");
////    fflush(stdout);
//    if (*lhci != NULL) {
////        printf("free lhci %p\n", (*lhci));
//        mgn_limit_free((*lhci)->limits);
//        mgn_lhci_free(*lhci);
//    }
//
//    *lhci = mgn_init_new_lhci(size
//                              ,x->size2,limits);
//}

void mgnp_moeadrbf_mdl_defparam(mgnp_moeadrbf_data *moeadrbf, size_t w_size, double cr, double mr)
{
    mgnp_moeadrbf_mdl_p *data = moeadrbf->model_data;
    data->mdl_size = moeadrbf->mdl_size;
    data->mdl_k = moeadrbf->mdl_k;

//    pmgn_moeadrbf_lhci_update(moeadrbf, moeadrbf->tset);

    data->search_lim = moeadrbf->search_lim;
    data->w_size = moeadrbf->tset->size;
//    mrbf->lhci =
//    mrbf->l_lhci_lim = false;

    data->rbf_data = moeadrbf->rbf_data;

    data->runs = 1000;
    data->neighbours = 10;

    mgn_ga_sets *gaprob = calloc(1, sizeof(*gaprob));
    gaprob->cross_rate = cr;
    gaprob->mut_rate = mr;
    gaprob->pbm_n = 5;
    gaprob->sbx_m = 20;

    size_t xsize = moeadrbf->solution->iparams.x_size;
    gaprob->mut_llim = calloc(xsize, sizeof(gaprob->mut_llim));
    gaprob->mut_ulim = calloc(xsize, sizeof(gaprob->mut_ulim));
    for (size_t i = 0; i < xsize; ++i) {
        gaprob->mut_llim[i] = 0;
        gaprob->mut_ulim[i] = 1;
    }

    data->ga_probs = gaprob;
}

void mgnp_moeadrbf_mdl_mop(gsl_vector *x, gsl_vector *f, gsl_vector *g, mgnp_moeadrbf_mdl_p* prm)
{
    UNUSED(g);
    gsl_matrix *mphi = prm->mphi;
    gsl_matrix_view m_x = gsl_matrix_view_array(x->data,1,x->size);
    gsl_matrix_view m_f = gsl_matrix_view_array(f->data,1,f->size);

//    gsl_matrix_printf(&m_f.matrix,stdout);

    gsl_matrix *m_fres = gsl_matrix_calloc(1, f->size);

    for (size_t i = 0; i < prm->mdl_size; ++i) {
        mgn_rbf_create_phi(&m_x.matrix,prm->km,prm->rbf_data[i].mdl_rbf.p_s,*prm->rbf_data[i].kernel
                           ,mphi);
//        mgn_rbf_new_weight(prm->mdl_rbf->m_phi,prm->tset->f,prm->m_w);
//printf("mdl sizes %zu %zu, %zu %zu, %zu %zu\n"
//       ,mphi->size1, mphi->size2
//       ,prm->mdl_rbf->p_m_w->size1, prm->mdl_rbf->p_m_w->size2
//       ,m_f.matrix.size1, m_f.matrix.size2
//       );
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1
                       ,mphi,prm->rbf_data[i].mdl_rbf.p_m_w,0,&m_f.matrix);

        gsl_matrix_scale(&m_f.matrix, prm->lambda->data[i]);
        gsl_matrix_add(m_fres,&m_f.matrix);
    }
    gsl_matrix_swap(m_fres,&m_f.matrix);
//    gsl_matrix_printf(&m_f.matrix,stdout);

    gsl_matrix_free(m_fres);
}

void mgnp_moeadrbf_mdl_mop_param_free(mgnp_moeadrbf_data *moeadrbf)
{
    free(moeadrbf->model_data->ga_probs->mut_ulim);
    free(moeadrbf->model_data->ga_probs->mut_llim);
    free(moeadrbf->model_data->ga_probs);
    gsl_matrix_free(moeadrbf->model_data->iW);
    free(moeadrbf->model_data);
}

void mgnp_moeadrbf_mdl_optim(mgnp_moeadrbf_mdl_p *mdl_p, mgn_pop *p_e, mgn_ptr *m_wl)
{
    size_t T = mdl_p->neighbours; // neighbour size
    size_t runs = mdl_p->runs;

    mgn_indv_param *params = &p_e->iparams;
    mgnMop *mop = mgn_mop_alloc();
    mop->eval = (mgn_mop_f) mgnp_moeadrbf_mdl_mop;
    mop->params = mdl_p;
    mop->limits = mdl_p->mop_limits;


    // TODO use given weight v
    mgn_popl *tmppe = mgn_popl_alloc(p_e->ops, params);
    mgnMoa *moead = mgn_moead_de_init(mdl_p->iW, params->f_size, T, tmppe, mop
                                      , mgn_init_transition,0,false);
    moead->max_exec = SIZE_MAX;
    moead->set_ga_vals(moead,mdl_p->ga_probs);
    mgn_initializer *l_lhci = mgn_pinit_lhc_alloc(mgn_moead_getpop(moead),mdl_p->search_lim);
    mgn_moead_pop_init_eval(moead,l_lhci);
    mgn_pinit_free(l_lhci);

    mgn_moead_set_scalarization(moead, mgn_scalar_pbi_ori);

    mgn_moa_solve(moead, runs); // 30_000
//    printf("tot_exec %zu\n", moead->tot_exec);

    // prepare output population
    mgn_pop *moead_pop = mgn_moead_getpop(moead);
    mgn_pop_exchange_iarray(p_e, moead_pop);
    // prepare weight matrix for return
    gsl_matrix *m_tmpw = mgn_moead_get_w(moead);
    m_wl->p = gsl_matrix_alloc(m_tmpw->size1, m_tmpw->size2);
    gsl_matrix_memcpy(m_wl->p, m_tmpw);

    // free --------------
    mgn_popl_free(tmppe);

    mgn_moead_free(moead);
    mgn_mop_free(mop);
}

// m_w model internal weight
// m_ws reference vectors
void mgnp_moeadrbf_select(mgn_count_ciclic *sel_idx
                          ,gsl_matrix *m_w
                          ,gsl_matrix *m_ws
                          ,mgn_pop *p_dest
                          ,mgn_pop* p_orig
                          )
{
    sel_idx->max = floor(m_w->size1 / (double) m_ws->size1);

//    gsl_matrix *m_ws = gsl_matrix_alloc(m_ws_size, m_w->size2);
//    gsl_vector_int *sel = gsl_vector_int_alloc(m_ws->size1);
//    gsl_matrix_printf(m_ws,stdout);

    gsl_matrix *m_dist = gsl_matrix_dist(m_ws,m_w,2.0);
    gsl_matrix_int *m_dist_i = gsl_matrix_int_alloc(m_dist->size1, m_dist->size2);
//    puts("m_ws");
//    gsl_matrix_printf(m_ws, stdout);
//    puts("m_w");
//    gsl_matrix_printf(m_w,stdout);
//    puts("mdist");
//    gsl_matrix_printf(m_dist,stdout);

    gsl_matrix_distrank_index(m_dist,m_dist_i);

    // TODO change comments to DEBUG clauses
    mgnp_selhash *hashes = NULL;
    mgnp_selhash *sh, *cur_sh; // = malloc(sizeof(*sh));
    for (size_t i = 0; i < m_ws->size1; ++i) {
//        printf("searching %zu\n", i);
        gsl_vector_int_view d_row = gsl_matrix_int_row(m_dist_i,i);
//        sel->data[i] = d_row.vector.data[*sel_idx];

        size_t j = 1;
        size_t i_sel = d_row.vector.data[sel_idx->value];

        // TODO make subfunction (helpers)
        //      check uhash examples for quick idea
        // sh is NULL if not found
        HASH_FIND(hh,hashes
                  ,mgn_indv_getx_vec(p_orig,i_sel)->data
                  ,sizeof(double) * p_orig->iparams.x_size
                  ,sh);

//        printf("what pointer is %p\n", sh);

        // if NULL value not in hash
        if(sh == NULL) {
            p_dest->ops->copy(mgn_pop_get(p_dest,i),mgn_pop_get(p_orig,i_sel));

            sh = malloc(sizeof(*sh));
            sh->id = i;
            sh->key = mgn_indv_getx_vec(p_orig,i_sel)->data;
            HASH_ADD_KEYPTR(hh,hashes
                            ,sh->key
                            ,sizeof(double) * p_orig->iparams.x_size
                            ,sh);
            sh = NULL;
        } else {
            do {
//                puts("insert");
                i_sel = d_row.vector.data[mgn_count_add(*sel_idx,1)];
                HASH_FIND(hh,hashes
                          ,mgn_indv_getx_vec(p_orig,i_sel)->data
                          ,sizeof(double) * p_orig->iparams.x_size
                          ,sh);
                j++;

                if(j >= sel_idx->max) {
                    break;
                }
            } while (sh != NULL);

            p_dest->ops->copy(mgn_pop_get(p_dest,i),mgn_pop_get(p_orig,i_sel));

            sh = malloc(sizeof(*sh));
            sh->id = i;
            sh->key = mgn_indv_getx_vec(p_orig,i_sel)->data;
            HASH_ADD_KEYPTR(hh,hashes
                            ,sh->key
                            ,sizeof(double) * p_orig->iparams.x_size
                            ,sh);
            sh = NULL;
        }
    }
    // free mem
    HASH_ITER(hh,hashes,cur_sh,sh){
        HASH_DEL(hashes, cur_sh);  /* delete; users advances to next */
        free(cur_sh);
    }

    // TODO for test case for distance rank, ok!
//    double a[6] = {2,5,6,3,1,3};
//    double b[4] = {7,4,5,9};
//    double res[6] = {sqrt(26), 5,sqrt(2),sqrt(37),sqrt(37),sqrt(56)};
//    int *order = {1,0,0,1,0,1};
//
//    gsl_matrix_view aa = gsl_matrix_view_array(a,3,2);
//    gsl_matrix_view bb = gsl_matrix_view_array(b,2,2);
//
//    gsl_matrix *c = gsl_matrix_dist(&aa.matrix,&bb.matrix,2.0);
//    gsl_matrix_printf(c,stdout);
//    gsl_matrix_int *dd = gsl_matrix_int_alloc(c->size1, c->size2);
//    gsl_matrix_distrank_index(c,dd);
//
//    gsl_matrix_int_fprintf(stdout,dd, "%d");
//    gsl_matrix_printf(c,stdout);

    mgn_count_sum(sel_idx,1);

    gsl_matrix_int_free(m_dist_i);
    gsl_matrix_free(m_dist);
}

// TODO indv protoype
void mgnp_moeadrbf_update_refine(mgn_pop *pop_newt, mgn_pop *pop_sel)
{
    mgnp_selhash *sel_inserted = NULL;
    mgnp_selhash *sh, *cur_sh;

    // Add all selected solutions to ht
    for (size_t i = 0; i < pop_newt->size; ++i) {
//        mgn_indv *in = mgn_pop_get(pop_sel,i);
        // insert
        sh = malloc(sizeof(*sh));
        sh->id = i;
        sh->key = mgn_indv_getx_vec(pop_newt,i)->data;
        HASH_ADD_KEYPTR(hh,sel_inserted
                        ,sh->key
                        ,sizeof(double) * pop_newt->iparams.x_size
                        ,sh);
        sh = NULL;
    }

    gsl_vector_int *sel_pos_in_newt = gsl_vector_int_alloc(pop_sel->size);
    size_t miss = 0;
    gsl_vector_int_set_all(sel_pos_in_newt,-1);
    for (size_t i = 0; i < pop_sel->size; ++i) {
        HASH_FIND(hh,sel_inserted
                  ,mgn_indv_getx_vec(pop_sel,i)->data
                  ,sizeof(double) * pop_sel->iparams.x_size
                  ,sh);

        sel_pos_in_newt->data[i] = (sh == NULL)? -1 : sh->id;
        miss += (sh == NULL)? 1 : 0;
        sh = NULL;
    }

    if( miss > 0 ) {

        mgn_pop_matrix *m_newt = mgn_pop_to_popm(pop_newt);
        // make matrix with missing members
        mgn_pop_matrix *m_sel = mgn_pop_to_popm(pop_sel);
        // meassure distance from missing to all
        gsl_matrix *m_dist = gsl_matrix_dist(m_sel->f,m_newt->f,2.0);
        gsl_matrix_int *m_dist_i = gsl_matrix_int_alloc(m_dist->size1, m_dist->size2);
        gsl_matrix_distrank_index(m_dist,m_dist_i);

        // TODO make helper function
        int pos;
        bool found;
        for (size_t i = 0, j = 0, l_pos; i < sel_pos_in_newt->size; ++i) {

            pos = sel_pos_in_newt->data[i];
            if(pos == -1) {
                // find in vector
                gsl_vector_int_view d_row = gsl_matrix_int_row(m_dist_i,i);
                do {
                    found = false;
                    l_pos = d_row.vector.data[j];
                    for (size_t k = 0; k < sel_pos_in_newt->size; ++k) {
                        if(k == l_pos) { found = true; }
                    }
                    j++;
                } while (found && j <= d_row.vector.size);
                // copy indv
                mgn_pop_copy(pop_newt,pop_sel,j,i,1);
            }
        }
        mgn_pop_matrix_free(m_newt);
        mgn_pop_matrix_free(m_sel);
        gsl_matrix_free(m_dist);
        gsl_matrix_int_free(m_dist_i);
    }

    // pick closest to missing only if it is not in sel_pos_in_newt


    // free mem
    HASH_ITER(hh,sel_inserted,cur_sh,sh){
        HASH_DEL(sel_inserted, cur_sh);  /* delete; users advances to next */
        free(cur_sh);
    }

    gsl_vector_int_free(sel_pos_in_newt);
}

static size_t crun = 0;
void mgnp_moeadrbf_update(mgnp_moeadrbf_data *mrbf, mgn_pop *pop_sel)
{
//    char *filename = calloc(64,1);
//    asprintf(&filename, "solution_before_%zu", crun);
//    mgn_plot_fast(mrbf->solution, filename, "sol");
    for (size_t i = 0; i < pop_sel->size; ++i) {
        mgn_pop_insert_dom(mrbf->solution,mgn_pop_get(pop_sel,i));
    }

//    free(filename);

    mgn_pop *pop_tset = mgn_pop_matrix_to_pop(mrbf->tset
                                              , (void*)mrbf->p_aprox->ops, &mrbf->p_aprox->iparams);

    mgn_pop *pop_u = mgn_pop_join(pop_tset, pop_sel);
    mgn_pop_prank_sort(pop_u);
    mgn_pop_copy(pop_tset,pop_u,0,0,mrbf->tset->size);

//    filename = calloc(64,1);
//    asprintf(&filename, "ptest_%zu", crun);
//    mgn_plot_fast(pop_tset, filename, "sol");
//    free(filename);
//    crun++;

    mgnp_moeadrbf_update_refine(pop_tset,pop_sel);

    mgn_pop_matrix_free(mrbf->tset);
    mrbf->tset = mgn_pop_to_popm(pop_tset);

    mgn_pop_free(pop_u);
    mgn_pop_free(pop_tset);
}


void mgnp_moeadrbf_pop_update(mgnp_moeadrbf_data *mrbf)
{
    mgn_pop_proto *pop = mrbf->solution;
    mgnLimit *lim = mrbf->search_lim;

    gsl_vector *min = mgn_pop_min_column(pop);
    gsl_vector *max = mgn_pop_max_column(pop);

    memcpy(lim->min, min->data, sizeof(double) * lim->size);
    memcpy(lim->min, max->data, sizeof(double) * lim->size);

    gsl_vector_free(min);
    gsl_vector_free(max);
}


void pmgn_moeadrbf_calcNN(mgn_model_rbfdata *models
                          , mgn_pop_matrix *tset
                          , cluster_data *cdat
                          , cluster_data_extra *cdata_extra
)
{
    // delete s if exists
    if(models->mdl_rbf.p_s) {
        gsl_vector_free(models->mdl_rbf.p_s);
    }
    models->mdl_rbf.p_s = mgn_cluster_var_dist(cdat->centers,cdata_extra,tset->x,true);

    pmgn_moeadrbf_calcPhi(models,tset,cdat);

    //        gsl_matrix *y_p = gsl_matrix_alloc(moeadrbf->tset->f->size1, moeadrbf->tset->f->size2);
    //        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1
    //                       ,moeadrbf->mdl_rbf[i].m_phi
    //                       ,moeadrbf->mdl_rbf[i].p_m_w,0
    //                       ,y_p);

    //
    //#ifdef DEBUG
    //        asprintf(&filename, "mrbf_pset-%zu-%zu",moa->c_run,i);
    //        mgn_plot_matrix_2d(y_p,filename,"p",0);
    //        asprintf(&filename, "mrbf_pset-%zu-%zu.txt",moa->c_run,i);
    //        gsl_matrix_save(y_p, filename);
    //        // end DEBUG
    //#endif
}

void pmgn_moeadrbf_calcPhi(mgn_model_rbfdata *models
                          , mgn_pop_matrix *tset
                          , cluster_data *cdat
)
{
    models->mdl_rbf.m_phi = mgn_rbf_create_phi(tset->x,cdat
                                               ,models->mdl_rbf.p_s
                                               ,models->kernel
                                               ,models->mdl_rbf.m_phi);

    models->mdl_rbf.p_m_w = mgn_rbf_new_weight(models->mdl_rbf.m_phi
                                               , tset->f, models->mdl_rbf.p_m_w);
}
