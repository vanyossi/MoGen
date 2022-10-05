/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_moead_fcrbf_inv.h"

#include <gsl/gsl_blas.h>
#include <float.h>

#include "mgn_rbf.h"
#include "mgn_fcmeans.h"
#include "mgn_pop_helper.h"
#include "mgn_gnuplot.h"

#include "mgn_weights.h"
#include "mgn_scalarization.h"

typedef struct mgnp_moeadrbf_data_inv {
    mgnp_moeadrbf_data mrbf;
    gsl_matrix *w_aux;
    gsl_matrix_uint *w_dist; // lowest g(x) to high index
    gsl_vector *z;
    mgnMop *zmop;
    mgnf_decomp_scalar gx;
} mgn_moeadrbf_data_inv;


void mgnp_moead_update_z(gsl_vector* x, gsl_vector* f, gsl_vector* g, double* z)
{
    UNUSED(x);
    UNUSED(g);
    double fi;

    for (size_t i = 0; i < f->size; ++i) {
        fi = gsl_vector_get(f,i);
        if (fi < z[i]) {
            z[i] = fi;
        }
    }
}

void mgn_rbfinv_populate_wdist(mgn_moeadrbf_data_inv *data)
{
    gsl_matrix *w = data->mrbf.m_w;
    mgn_popl *popl = data->mrbf.arc;

    gsl_vector *dist = gsl_vector_alloc(popl->size);
    int *dindex;

    // for each weight, meassure scalar value and add to w_dist
    for (size_t i = 0; i < w->size1; ++i) {
        double sval;
        gsl_vector_view wcur = gsl_matrix_row(w,i);
        for (size_t p = 0; p < popl->size ; ++p) {
            sval = data->gx(&wcur.vector
                     ,popl->ops->get_iparams(popl->get(popl,p)).f
                     ,data->z
                     ,0);

            gsl_vector_set(dist,p,sval);
        }
        dindex = gsl_vector_qsort(dist);
        gsl_vector_int_view dview = gsl_vector_int_view_array(dindex,data->w_dist->size2);
        gsl_matrix_int_set_row(data->w_dist,i,&dview.vector);
    }
    return;
}
/*
 *
 * Main alg run
 *
 *
 */
void mgn_fcmeans_convert_tox(gsl_matrix *x
                            , cluster_data_extra *cdata
                            , fcmeans_data *fcm
                            )
{
    // prepare centroid vector
    size_t csize = fcm->centers->size1;
    gsl_matrix_free(fcm->centers);
    fcm->centers = gsl_matrix_alloc(csize, x->size2);

    gsl_vector *sum = gsl_vector_calloc(x->size2);

    for (size_t i = 0; i < cdata->size; ++i) {
        size_t size = cdata->mpos[i].size;
        for (size_t i1 = 0; i1 < size; ++i1) {
            gsl_vector_view row = gsl_matrix_row(x,cdata->mpos[i].pos[i1]);
            gsl_vector_add(sum, &row.vector);
        }
        gsl_vector_scale(sum, 1.0/size);
        gsl_matrix_set_row(fcm->centers,i,sum);
        gsl_vector_set_zero(sum);
    }
    return;
}

// this represents one loop
// TODO remove comments, polish args + API
void mgn_moeadrbf_fcinv_run(mgnMoa *moa)
{
    mgnp_moeadrbf_data *moeadrbf = mgn_moeadrbf_features(moa);
    mgn_moeadrbf_data_inv *d_extra = (mgn_moeadrbf_data_inv*)moeadrbf;

    // === Model building

//    fcmeans_data *fcdat = mgn_fcmeans(moeadrbf->tset->x,moeadrbf->mdl_k, 1000, 1e-5);
    fcmeans_data *fcdat = mgn_fcmeans(moeadrbf->tset->f,moeadrbf->mdl_k, 1000, 1e-5);
    double low_e = (moa->tot_exec < moa->max_exec * 0.80)? 0.02 : (moa->tot_exec - 1.0) / moa->max_exec *0.1;

//    printf("low_e: %zu, %zu :: %.6f,  %.6f\n", moa->tot_exec, moa->max_exec, (moa->tot_exec - 1.0) / moa->max_exec * 0.3, low_e);
    cluster_data_extra *cl_extra = mgn_fcmeans_calc(fcdat,low_e,1);
    mgn_fcmeans_convert_tox(moeadrbf->tset->x, cl_extra, fcdat);
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
        pmgn_moeadrbf_calcNN(&moeadrbf->rbf_data[i],moeadrbf->tset,&cdat,cl_extra);

        mgnp_moeadrbf_optim_s(&moeadrbf->rbf_data[i].mdl_rbf
                              , &cdat
                              ,moeadrbf->tset
                              , moeadrbf->rbf_data[i].kernel);

        pmgn_moeadrbf_calcPhi(&moeadrbf->rbf_data[i],moeadrbf->tset,&cdat);
    }


    // weight bootstrap
    gsl_vector *lambda = mgnp_moeadrbf_find_lambda(moeadrbf->tset
                                                   , moeadrbf->rbf_data, moeadrbf->mdl_size);

    // === Evaluate using moead
    // === Run MOEA/D -- aproximation
    gsl_matrix *m_1phi = gsl_matrix_alloc(1,moeadrbf->mdl_k);

    moeadrbf->model_data->lhci = moeadrbf->lhci;
    moeadrbf->model_data->lambda = lambda;
    moeadrbf->model_data->km = &cdat;
    moeadrbf->model_data->mphi = m_1phi;

    mgn_ptr *m_w_ptr = malloc(sizeof(*m_w_ptr));
    mgnp_moeadrbf_mdl_optim(moeadrbf->model_data, moeadrbf->p_aprox, m_w_ptr);


    // === Select Points
    // mgn_popl pop_p: has all nondom aprox solutions
    // S set size of RBF cluster vectors.
    mgn_pop *pop_sel = mgn_pop_alloc(moeadrbf->m_w->size1
                                     ,(void*)moeadrbf->solution->ops, &moeadrbf->solution->iparams);

    mgnp_moeadrbf_select(&moeadrbf->sel_idx
                         ,m_w_ptr->p
                         ,moeadrbf->m_w
                         ,pop_sel
                         ,moeadrbf->p_aprox);

    mgn_mop_eval_pop(moa->mop, pop_sel, NULL);
    moa->tot_exec += pop_sel->size;


    // === Update population
    mgnp_moeadrbf_update(moeadrbf, pop_sel);
    mgnp_moeadrbf_pop_update(moeadrbf);

    // === free all
    mgn_pop_free(pop_sel);

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
mgnMoa* mgn_moa_moeadrbf_fcinv_alloc(
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

    strncpy(moa->name, "MOEAD-FCRBF_INV", MOA_NAME_LEN);
    moa->run = mgn_moeadrbf_fcinv_run;

    mgnp_moeadrbf_data *mrbf = mgn_moeadrbf_features(moa);
    mgn_moeadrbf_data_inv *d_extra = realloc(mrbf, sizeof(*d_extra));
    //assert not null TODO
    if (d_extra == NULL) return NULL;
//    memcpy(d_extra,mrbf, sizeof(*mrbf));

    moa->features = d_extra;

    d_extra->mrbf.arc = mgn_popl_alloc((void*)A->ops,&A->iparams);
    d_extra->w_aux = mgn_weight_slattice(Nw*2, A->iparams.f_size);
    // get the 20 closest
    size_t au_n = 20;
    d_extra->w_dist = gsl_matrix_uint_alloc(d_extra->mrbf.m_w->size1, au_n);

    d_extra->z = gsl_vector_alloc(A->iparams.f_size);
    gsl_vector_set_all(d_extra->z, DBL_MAX);
    d_extra->zmop = mgn_mop_alloc();
    d_extra->zmop->eval = mgn_cast_eval(mgnp_moead_update_z);
    d_extra->zmop->params = d_extra->z->data;

    d_extra->gx = mgn_scalar_tchebycheff;

    return moa;
}

void mgn_moa_moeadrbf_fcinv_free(mgnMoa* moeadrbf)
{
    mgn_moa_moeadrbf_common_free(moeadrbf);
}

void mgn_moa_moeadrbf_fcinv_init(mgnMoa* moeadrbf)
{
    mgn_moa_moeadrbf_common_init(moeadrbf);

    mgn_moeadrbf_data_inv *dai = (mgn_moeadrbf_data_inv*)mgn_moeadrbf_features(moeadrbf);
    mgn_popl_insert_popm(dai->mrbf.arc,dai->mrbf.tset);

    // calculate z value
    mgn_mop_eval_pop(dai->zmop,(mgn_pop*)dai->mrbf.arc, NULL);
    mgn_rbfinv_populate_wdist(dai);

}
