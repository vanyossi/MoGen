/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_moead-rbf.h"

#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>

#include "population.h"
#include "mgn_poplist.h"
#include "mgn_pop_matrix.h"
#include "mgn_pop_helper.h"
#include "mgn_initializer.h"
#include "individual.h"
#include "mgn_mop.h"

#include "mgn_pareto.h"
#include "gsl_vector_additional.h"
#include "mgn_vector_distance.h"
#include "mgn_weights.h"
#include "mgn_counter.h"

#include "mgn_rbf.h"
#include "mgn_de.h"
#include "mgn_moead.h"

#include "uthash.h"
#include "mgn_gnuplot.h"

struct mgnp_features_moeadrbf {
    gsl_vector *sigma;
};

// pointer passing helper
typedef struct mgn_pointer {
    void* p;
} mgn_ptr;

typedef struct mgn_moeadrbf_data mgnp_moeadrbf_data;
// TODO add function pointers as types
typedef void (*mgn_kernel_f)(gsl_vector *r, double s);

struct mgnp_rbf_weigts {
    gsl_vector *p_s;
    gsl_matrix *p_m_w; // used for internal model
    gsl_matrix *m_phi;
};

struct mgn_moeadrbf_data {
    pmgn_moeadrbf_templatep()
    /* private */
    mgn_pop_matrix *tset; // tset pop size(nt)
    mgn_lhci *lhci; // for internal moead pop init
    bool l_lhci_lim; // true to free lhci limits
//    gsl_matrix *pm;    // pop for eval with model size(n)
    int scalarf;
    size_t mdl_size;// scalarization function, this is fixed to PBI
    size_t mdl_k;
    mgn_kernel_f kernel[3];
    struct mgnp_rbf_weigts mdl_rbf[3]; // used for internal model
    mgn_pop *p_aprox;
    mgn_count_ciclic sel_idx;
};

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

// moead de min optim
struct mgnp_moeadrbf_s_optim_ef_p {
    mgnp_de_ef_param_m()
    gsl_matrix *y_t; //training matrix
//    gsl_vector *y_cur;
};

struct mgnp_moeadrbf_s_optim_p {
    mgn_mop_param_common();
    mgn_pop_matrix *tset; //training matrix
    kmeans_data *km;
    gsl_matrix *f_p;
    gsl_matrix *m_w;
    gsl_matrix *m_phi;
    mgn_kernel_f rbf;
};

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

// x is sigma
// f is trained values y_predicted
// g is none
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
        gsl_vector_view f_pp = gsl_matrix_row(prm->f_p,0);
        gsl_vector_set(f, 0,mgn_math_mse_matrix(prm->tset->f,prm->f_p));
    } else {
        gsl_vector_set(f,0,1e10);
    }
}


void mgnp_moeadrbf_optim_s(
    struct mgnp_rbf_weigts *vals
    , kmeans_data *km
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

    mgnMop *mop = mgn_mop_alloc();
    mop->eval = mgn_cast_eval(mgnp_moeadrbf_s_optim_mop);
    gsl_matrix *tmp_f_p = gsl_matrix_alloc(tset->f->size1, tset->f->size2);
    gsl_matrix *tmp_w_p = gsl_matrix_alloc(km->k, tset->f->size2);
    gsl_matrix *tmp_m_phi = gsl_matrix_alloc(tset->x->size1, km->k);
    struct mgnp_moeadrbf_s_optim_p mop_p = {0, tset, km, tmp_f_p, tmp_w_p, tmp_m_phi, rbf};
    mop->params = &mop_p;

    size_t Np = 50; // de pop size

    mgn_lhci *lhci = mgn_init_new_lhci(Np,param.x_size,rlim);
    mgnMoa* de = mgn_moa_de_alloc(Np,iops,&param,1.4, 0.7);
    mgn_de_init(de, mgn_init_lhc, lhci);
    mgn_lhci_free(lhci);

    struct mgnp_moeadrbf_s_optim_ef_p *ef_prm = malloc(sizeof(*ef_prm));
        ef_prm->y_t = tset->f;
//        ef_prm->y_cur = gsl_vector_alloc(tset->f->size2);
    mgn_de_setmop(de, mop, mgn_cast_de_ef(mgnp_moeadrbf_s_optim_mop_min),(mgn_de_ef_param*)ef_prm);
    mgn_de_eval(de);

    mgn_moa_solve(de, 10);

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


// find the vector weidths lambda
// this could be done solving Ax=b system
gsl_vector *mgnp_moeadrbf_find_lambda(mgn_pop_matrix *tset, struct mgnp_rbf_weigts *res, size_t size)
{
    gsl_vector *lambda = gsl_vector_calloc(size);

    gsl_matrix *f_vals = gsl_matrix_alloc(tset->f->size1,size);
    gsl_matrix *tmp_y_p = gsl_matrix_alloc(tset->f->size1, tset->f->size2);

    for (size_t i = 0; i < size; ++i) {
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1
                       ,res[i].m_phi,res[i].p_m_w,0,tmp_y_p);
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


// ======== MOEA/D model =========
typedef struct mgnp_moeadrbf_mdl_p {
    size_t w_size;
    size_t mdl_size;// scalarization function, this is fixed to PBI
    size_t mdl_k;
    mgn_lhci *lhci;
    mgn_kernel_f *kernel;
    struct mgnp_rbf_weigts *mdl_rbf;
    gsl_vector *lambda;
    kmeans_data *km;
    gsl_matrix *mphi;
} mgnp_moeadrbf_mdl_p;


void mgnp_moeadrbf_mdl_mop(gsl_vector *x, gsl_vector *f, gsl_vector *g, mgnp_moeadrbf_mdl_p* prm)
{
    UNUSED(g);
    gsl_matrix *mphi = prm->mphi;
    gsl_matrix_view m_x = gsl_matrix_view_array(x->data,1,x->size);
    gsl_matrix_view m_f = gsl_matrix_view_array(f->data,1,f->size);

//    gsl_matrix_printf(&m_f.matrix,stdout);

    gsl_matrix *m_fres = gsl_matrix_calloc(1, f->size);

    for (size_t i = 0; i < prm->mdl_size; ++i) {
        mgn_rbf_create_phi(&m_x.matrix,prm->km,prm->mdl_rbf[i].p_s,*prm->kernel[i]
                           ,mphi);
//        mgn_rbf_new_weight(prm->mdl_rbf->m_phi,prm->tset->f,prm->m_w);
//printf("mdl sizes %zu %zu, %zu %zu, %zu %zu\n"
//       ,mphi->size1, mphi->size2
//       ,prm->mdl_rbf->p_m_w->size1, prm->mdl_rbf->p_m_w->size2
//       ,m_f.matrix.size1, m_f.matrix.size2
//       );
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1
                       ,mphi,prm->mdl_rbf[i].p_m_w,0,&m_f.matrix);

        gsl_matrix_scale(&m_f.matrix, prm->lambda->data[i]);
        gsl_matrix_add(m_fres,&m_f.matrix);
    }
    gsl_matrix_swap(m_fres,&m_f.matrix);
//    gsl_matrix_printf(&m_f.matrix,stdout);

    gsl_matrix_free(m_fres);
}


void mgnp_moeadrbf_mdl_optim(mgnp_moeadrbf_mdl_p *mdl_p, mgn_pop *p_e, mgn_ptr *m_wl)
{
    size_t T = 25; // neighbourg size
    size_t runs = 5000;
    double cr = 0.9;
    double mr = 0.1;

    mgnMop *mop = mgn_mop_alloc();
    mop->eval = (mgn_mop_f) mgnp_moeadrbf_mdl_mop;
    mop->params = mdl_p;

    mgn_indv_param *params = &p_e->iparams;

    // TODO reduce initialize probs
    mgn_ga_sets ga_probs = {cr, mr, NULL, NULL};
    ga_probs.mut_llim = calloc(params->x_size, sizeof(ga_probs.mut_llim));
    ga_probs.mut_ulim = calloc(params->x_size, sizeof(ga_probs.mut_ulim));
    for (size_t i = 0; i < params->x_size; ++i) {
        ga_probs.mut_llim[i] = 0;
        ga_probs.mut_ulim[i] = 1;
    }

    // TODO use given weight v
    mgn_popl *tmppe = mgn_popl_alloc(p_e->ops, params);
    mgnMoa *moead = mgn_moead_init(mdl_p->w_size, params->f_size, T, tmppe, mop
                                   , mgn_init_lhc,mdl_p->lhci,false);
    moead->set_ga_vals(moead,&ga_probs);

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

    free(ga_probs.mut_llim);
    free(ga_probs.mut_ulim);

    mgn_moead_free(moead);
    mgn_mop_free(mop);
}

// ============================ MOEAD model end
/*
 *
 * Selection and update
 *
 *
 */
typedef struct mgnp_moeadrbf_selhash {
    int id;
    double* key;
    UT_hash_handle hh;
} mgnp_selhash;

// m_w matrix of w vectors in moead
// m_ws mgn_ptr with matrix of w vectors of moeadRBF
// sel: mgn_pop selected pop size = m_ws row size.
void mgnp_moeadrbf_select(mgn_count_ciclic *sel_idx
                                ,gsl_matrix *m_w
                                ,gsl_matrix *m_ws
                                ,mgn_pop *p_dest
                                ,mgn_pop* p_orig)
{
    sel_idx->max = floor(m_w->size1 / (double) m_ws->size1);
//    gsl_matrix *m_ws = gsl_matrix_alloc(m_ws_size, m_w->size2);
//    gsl_vector_int *sel = gsl_vector_int_alloc(m_ws->size1);
//    gsl_matrix_printf(m_ws,stdout);

    gsl_matrix *m_dist = gsl_matrix_dist(m_ws,m_w,2.0);
    gsl_matrix_int *m_dist_i = gsl_matrix_int_alloc(m_dist->size1, m_dist->size2);

    gsl_matrix_distrank_index(m_dist,m_dist_i);

    // TODO change comments to DEBUG clauses
    mgnp_selhash *hashes = NULL;
    mgnp_selhash *sh, *cur_sh; // = malloc(sizeof(*sh));
    for (size_t i = 0; i < m_ws->size1; ++i) {
        gsl_vector_int_view d_row = gsl_matrix_int_row(m_dist_i,i);
//        sel->data[i] = d_row.vector.data[*sel_idx];

        size_t j = 1;
        size_t i_sel = d_row.vector.data[sel_idx->value];

        // TODO make subfunction (helpers)
        // sh is NULL if not found
        HASH_FIND(hh,hashes
                  ,mgn_indv_getx_vec(p_orig,i_sel)->data
                  ,sizeof(double) * p_orig->iparams.x_size
                  ,sh);

//        printf("what pointer is %p\n", sh);

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
                i_sel = d_row.vector.data[mgn_count_add(*sel_idx,j)];
                HASH_FIND(hh,hashes
                          ,mgn_indv_getx_vec(p_orig,i_sel)->data
                          ,sizeof(double) * p_orig->iparams.x_size
                          ,sh);
                j++;
                if(j == sel_idx->max) {
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

/**
void mgnp_moeadrbf_select_pop(mgn_pop *p_dest,mgn_pop* p_orig, gsl_vector_int* index)
{
//    // TODO change comments to DEBUG clauses
//    mgnp_selhash *hashes = NULL;
//    mgnp_selhash *sh, *cur_sh; // = malloc(sizeof(*sh));
//    for (size_t i = 0; i < index->size; ++i) {
//        size_t j = 1;
//        size_t i_sel = index->data[i];
//
////        sh = malloc(sizeof(*sh));
////        sh->id = i;
////        sh->key = mgn_indv_getx_vec(p_orig,i_sel)->data;
////        HASH_ADD_KEYPTR(hh,hashes
////                        ,sh->key
////                        ,sizeof(double) * p_orig->iparams.x_size
////                        ,sh);
////        sh = NULL;
//
//        // sh is NULL if not found
//        HASH_FIND(hh,hashes
//                  ,mgn_indv_getx_vec(p_orig,i_sel)->data
//                  ,sizeof(double) * p_orig->iparams.x_size
//                  ,sh);
//
////        printf("what pointer is %p\n", sh);
//
//        if(sh == NULL) {
//            p_dest->ops->copy(mgn_pop_get(p_dest,i),mgn_pop_get(p_orig,i_sel));
//        } else {
//            do {
//                i_sel = index->data[i+j];
//                HASH_FIND(hh,hashes
//                          ,mgn_indv_getx_vec(p_orig,i_sel)->data
//                          ,sizeof(double) * p_orig->iparams.x_size
//                          ,sh);
//                j++;
//            } while (sh != NULL);
//        }
//    }
//    // free mem
//    HASH_ITER(hh,hashes,cur_sh,sh){
//        HASH_DEL(hashes, cur_sh);
//        free(cur_sh);
//    }
}
**/
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

    if(miss > 0 ) {

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
            j = 0;
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
                } while (found);
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

void mgnp_moeadrbf_update(mgnp_moeadrbf_data *mrbf, mgn_pop *pop_sel)
{
    for (size_t i = 0; i < pop_sel->size; ++i) {
        mgn_pop_insert_dom(mrbf->solution,mgn_pop_get(pop_sel,i));
    }
    mgn_pop *pop_tset = mgn_pop_matrix_to_pop(mrbf->tset
        , (void*)mrbf->p_aprox->ops, &mrbf->p_aprox->iparams);

    mgn_pop *pop_u = mgn_pop_join(pop_tset, pop_sel);
    mgn_pop_prank_sort(pop_u);
    mgn_pop_copy(pop_tset,pop_u,0,0,mrbf->tset->size);

    mgnp_moeadrbf_update_refine(pop_tset,pop_sel);

    mgn_pop_matrix_free(mrbf->tset);
    mrbf->tset = mgn_pop_to_popm(pop_tset);

    mgn_pop_free(pop_u);
    mgn_pop_free(pop_tset);
}

void mgnp_moeadrbf_pop_update(mgnp_moeadrbf_data *mrbf)
{
//    double mean = gsl_stats_mean(mrbf->tset->x->data
//                                 ,1
//                                 ,mrbf->tset->x->size1 * mrbf->tset->x->size2
//    );
//    double max = gsl_matrix_max(mrbf->tset->x) + mean;
//    double min = gsl_matrix_min(mrbf->tset->x) - mean;

//    printf("min max:: %.6f, %6f %.6f\n", mean, min, max);

    mgnLimit *limits = mgn_limit_alloc(mrbf->tset->x->size2);
    for (size_t i = 0; i < limits->size; ++i) {
        gsl_vector_view c_col = gsl_matrix_column(mrbf->tset->x,i);
        double r_mean = gsl_stats_mean(c_col.vector.data,1,c_col.vector.size);
        limits->min[i] = fmax(gsl_vector_min(&c_col.vector) - r_mean,0);
        limits->max[i] = fmin(gsl_vector_max(&c_col.vector) + r_mean,1);
    }

    // TODO make wrapper for this free
    // or accept change in all limits.
    if (mrbf->l_lhci_lim == true) {
        mgn_limit_free(mrbf->lhci->limits);
    }
    mgn_lhci_free(mrbf->lhci);
    mrbf->lhci = mgn_init_new_lhci(mrbf->tset->size * 2
        ,mrbf->solution->iparams.x_size,limits);
    mrbf->l_lhci_lim = true;
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
    mgnp_moeadrbf_data *moeadrbf = mgn_moeadrbf_features(moa);
    // === Model building
    kmeans_data *km = gsl_kmeans(moeadrbf->tset->x,moeadrbf->mdl_k, 1000);
    kmeans_data_extra *kme = gsl_kmeans_calc(km);

    char* filename = malloc(sizeof(char) * 64);

    for (size_t i = 0; i < moeadrbf->mdl_size; ++i) {

        if(moeadrbf->mdl_rbf[i].p_s) {
            gsl_vector_free(moeadrbf->mdl_rbf[i].p_s);
        }
        moeadrbf->mdl_rbf[i].p_s = mgn_kmeans_cluster_var_dist(km,kme,moeadrbf->tset->x,true);
//        puts("_-___________-");

        // DEBUG
        moeadrbf->mdl_rbf[i].m_phi = mgn_rbf_create_phi(moeadrbf->tset->x,km
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

        asprintf(&filename, "mrbf_pset-%zu-%zu",moa->c_run,i);
        mgn_plot_matrix_2d(y_p,filename,"p",0);
        asprintf(&filename, "mrbf_pset-%zu-%zu.txt",moa->c_run,i);
        gsl_matrix_save(y_p, filename);
        // end DEBUG

        mgnp_moeadrbf_optim_s(&moeadrbf->mdl_rbf[i]
                              , km
                              ,moeadrbf->tset
                              , moeadrbf->kernel[i]);

//        gsl_vector_fprintf(stdout,moeadrbf->mdl_rbf[i].p_s, "%.4f");
        moeadrbf->mdl_rbf[i].m_phi = mgn_rbf_create_phi(moeadrbf->tset->x,km
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
        asprintf(&filename, "mrbf_pset-%zu-%zu_p.txt",moa->c_run,i);
        gsl_matrix_save(y_p, filename);
    }


    gsl_vector *lambda = mgnp_moeadrbf_find_lambda(moeadrbf->tset
                        , moeadrbf->mdl_rbf, moeadrbf->mdl_size);
//    gsl_vector_fprintf(stdout,lambda,"%.8f");
//    exit(0);

    // === Evaluate using moead
    // === Run MOEA/D -- aproximation

    //mdl training value init
    gsl_matrix *m_1phi = gsl_matrix_alloc(1,moeadrbf->mdl_k);
    mgnp_moeadrbf_mdl_p mdl_p = {
        100
        ,moeadrbf->mdl_size
        ,moeadrbf->mdl_k
        ,moeadrbf->lhci
        ,moeadrbf->kernel
        ,moeadrbf->mdl_rbf
        ,lambda
        ,km
        ,m_1phi
    };

    mgn_ptr *m_w_ptr = malloc(sizeof(*m_w_ptr));
    mgnp_moeadrbf_mdl_optim(&mdl_p, moeadrbf->p_aprox, m_w_ptr);

    mgn_plot_data pdat = {"", "", "f_1", "f_2",
                          -0.1f,1.1f,-0.1f,1.1f};
    asprintf(&pdat.title, "MOEAD_RBF");
    asprintf(&pdat.filename, "%s-test_run_i_%s_%zu", pdat.title,moa->mop->name, moa->tot_exec);
    mgn_plot((mgn_pop_proto *) moeadrbf->p_aprox, &pdat);
    free(pdat.filename);

    // === Select Points
    // mgn_popl pop_p: has all nondom aprox solutions
    // S set size of RBF cluster vectors.
    printf("ws %zu %zu", ((gsl_matrix*)m_w_ptr->p)->size1, ((gsl_matrix*)m_w_ptr->p)->size2);
    size_t ws_size = 30; // TODO input
    // TODO initialize inside using mgn_ptr (probably)
    gsl_matrix *m_ws = mgn_weight_slattice_perm(ws_size-1,((gsl_matrix*)m_w_ptr->p)->size2);

    mgn_pop *pop_sel = mgn_pop_alloc(m_ws->size1
        ,moeadrbf->solution->ops, &moeadrbf->solution->iparams);

    mgnp_moeadrbf_select(&moeadrbf->sel_idx
                               ,m_w_ptr->p
                               ,m_ws,pop_sel
                               ,moeadrbf->p_aprox);

//    printf("sel index %zu\n",moeadrbf->sel_idx.value);

    // === Update population
    mgn_mop_eval_pop(moa->mop, pop_sel, NULL);
    moa->tot_exec += pop_sel->size;
    asprintf(&filename, "rbf_pop_sel-%zu.txt",moa->c_run);
    mgn_plot_fast(pop_sel, filename, "sel");

    mgnp_moeadrbf_update(moeadrbf, pop_sel);

    mgnp_moeadrbf_pop_update(moeadrbf);


    printf("msize: %zu %zu,,,, %zu %zu\n", m_ws->size1, m_ws->size2
           , ((gsl_matrix*)m_w_ptr->p)->size1, ((gsl_matrix*)m_w_ptr->p)->size2);
    printf("run %zu: tot_exec %zu, idx %zu\n",moa->c_run,moa->tot_exec, moeadrbf->sel_idx.value);

    // === free all
    gsl_matrix_free(m_ws);
    mgn_pop_free(pop_sel);
//    gsl_vector_int_free(sel_indexes);

    gsl_matrix_free(m_1phi);
    gsl_matrix_free(m_w_ptr->p);
    free(m_w_ptr);

    gsl_vector_free(lambda);
    gsl_kmeans_data_extra_free(kme);
    gsl_kmeans_free(km);

    asprintf(&filename, "rbf_final-%zu.txt",moa->c_run);
    FILE *out = fopen(filename,"w");
    mgn_pop_print(moeadrbf->solution, out);
    fclose(out);
    mgn_plot_fast(moeadrbf->solution, filename, "sol");

    free(filename);
}

mgn_pop_proto* mgn_moeadrbf_pop_get(mgnMoa *moa)
{
    return (mgn_pop_proto*)mgn_moeadrbf_features(moa)->solution;
}

//===== public functions
mgnMoa* mgn_moa_moeadrbf_alloc(
    size_t max_eval
    , size_t nt
    , size_t k
    , gsl_matrix *W
    , mgn_popl *A
    , mgnLimit *limits
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
//    mrbf->kernel = malloc(sizeof(mrbf->kernel) * mrbf->mdl_size);
    mrbf->kernel[0] = rbf_kernel_gauss;
    mrbf->kernel[1] = rbf_kernel_mqua;
    mrbf->kernel[2] = rbf_kernel_imqua;
    for (size_t i = 0; i < mrbf->mdl_size; ++i) {
        mrbf->mdl_rbf[i].p_s = NULL;
        mrbf->mdl_rbf[i].p_m_w = NULL;
        mrbf->mdl_rbf[i].m_phi = NULL;
    }

    // used to keep control of selected Neigbourh
    mrbf->sel_idx = (mgn_count_ciclic){5,0};
    mgn_init_lhc_to_matrix(mrbf->tset->x, limits);
    mrbf->lhci = mgn_init_new_lhci(nt*2,x_size,limits);
    mrbf->l_lhci_lim = false;

    // initialize tset population
    strncpy(moa->name, "MOEA/D-RBF", MOA_NAME_LEN);
    moa->tot_exec = 0;
    moa->c_run = 0;
    moa->run = mgn_moeadrbf_run;
    moa->stop = mgnp_moeadrbf_stop;
    moa->set_ga_vals = mgnp_moeadrbf_set_ga_vals;
    moa->features = mrbf;
    moa->pop_get = mgn_moeadrbf_pop_get;

    return moa;
}

void mgn_moa_moeadrbf_free(mgnMoa* moeadrbf)
{
    mgnp_moeadrbf_data *mrbf = mgn_moeadrbf_features(moeadrbf);

    mgn_pop_matrix_free(mrbf->tset);
    mgn_pop_free(mrbf->p_aprox);
    for (size_t i = 0; i < mrbf->mdl_size; ++i) {
        gsl_vector_free(mrbf->mdl_rbf[i].p_s);
        gsl_matrix_free(mrbf->mdl_rbf[i].p_m_w);
        gsl_matrix_free(mrbf->mdl_rbf[i].m_phi);
    }
    if (mrbf->l_lhci_lim) {
        mgn_limit_free(mrbf->lhci->limits);
    }
    mgn_lhci_free(mrbf->lhci);
//    gsl_matrix_free(mrbf->pm);
    free(mrbf);
    free(moeadrbf);
}

void mgn_moa_moeadrbf_init(mgnMoa* moeadrbf)
{
    mgnp_moeadrbf_data *mrbf = mgn_moeadrbf_features(moeadrbf);
    mgn_pop_matrix_eval(mrbf->tset,moeadrbf->mop);

//    size_t idx[2] = {0 ,1};
    mgn_plot_matrix_2d(mrbf->tset->f,"mrbf_tset_init", "f",0);
    gsl_matrix_save(mrbf->tset->f, "mrbf_tset_init_p.txt");

//    gsl_matrix_printf(mrbf->tset->x,stdout);
//    printf("max min, %.6f %.6f\n", gsl_matrix_min(mrbf->tset->x), gsl_matrix_max(mrbf->tset->x));

    mgn_pop_copy_mp(mrbf->solution,mrbf->tset);
    moeadrbf->tot_exec += mrbf->solution->size;

    // TODO print pop
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
