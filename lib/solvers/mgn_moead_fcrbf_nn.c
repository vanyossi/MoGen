/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_moead_fcrbf_nn.h"

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
    size_t wei_nei;
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
        free(dindex);
    }

    gsl_vector_free(dist);

    return;
}

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
    gsl_vector_free(sum);
    return;
}

#include <math.h>
// update limits
void mgnp_moeadrbf_popm_limit_update(mgnLimit *lim, gsl_matrix *sample, mgnLimit *bound)
{
    // TODO bound to mop limits
    double value;
    for (size_t col = 0; col < sample->size2; ++col) {
        gsl_vector_view ccol = gsl_matrix_column(sample, col);
        lim->min[col] = DBL_MAX;
        lim->max[col] = DBL_MIN;
        for (size_t j = 0; j < ccol.vector.size; ++j) {
            value = gsl_vector_get(&ccol.vector,j);
            if (value < lim->min[col]) {
                lim->min[col] = value;
            }
            if (value > lim->max[col]) {
                lim->max[col] = value;
            }
        }

        lim->min[col] = (bound->min[col] > lim->min[col])? bound->min[col]: lim->min[col];
        lim->max[col] = (bound->max[col] < lim->max[col])? bound->max[col]: lim->max[col];

//        lim->min[col] = gsl_vector_min(&ccol.vector);
//        lim->max[col] = gsl_vector_max(&ccol.vector);

//        for (size_t n = 0; n < ccol.vector.size; ++n) {
//            printf("%.6f ", gsl_vector_get(&ccol.vector,n));
//        }
//        printf(":: %.6f %.6f\n", lim->min[col], lim->max[col]);
    }
//    puts("");
}
/*
 *
 * Main alg run
 *
 *
 */

void mgn_moeadfc_rbf_update_training(mgnp_moeadrbf_data* data)
{
    // set cluster size
//    data->mdl_k = data->tset->size;
    // clean mphi matrix
    for (size_t i = 0; i < data->mdl_size; ++i) {
        gsl_matrix_free(data->rbf_data[i].mdl_rbf.m_phi);
        data->rbf_data[i].mdl_rbf.m_phi = NULL;
    }
}

// this represents one loop
// TODO remove comments, polish args + API
void mgn_moeadrbf_fcnn_run(mgnMoa *moa)
{
    mgnp_moeadrbf_data *moeadrbf = mgn_moeadrbf_features(moa);
    mgn_moeadrbf_data_inv *d_extra = (mgn_moeadrbf_data_inv*)moeadrbf;

    // for each reference vector create a nn
    gsl_matrix_uint_free(d_extra->w_dist);
    d_extra->w_dist = gsl_matrix_uint_alloc(d_extra->mrbf.m_w->size1, d_extra->wei_nei);
    mgn_rbfinv_populate_wdist(d_extra);

    size_t mp_size = d_extra->w_dist->size2;
    mgn_pop_matrix *p_wei_nei = mgn_pop_matrix_alloc(mp_size
                        ,moeadrbf->tset->x->size2
                        ,mp_size
                        ,moeadrbf->tset->f->size2
                        ,0 , 0);

    mgn_pop_matrix *pop_arch = mgn_pop_to_popm((mgn_pop_proto*)moeadrbf->arc);
    for (size_t v = 0; v < moeadrbf->m_w->size1; ++v) {
        // populate vector_training_set
        for (size_t i = 0; i < d_extra->w_dist->size2; ++i) {
            size_t sel = gsl_matrix_uint_get(d_extra->w_dist,v,i);
            gsl_vector_view v_selrow = gsl_matrix_row(pop_arch->x, sel);

            gsl_matrix_set_row(p_wei_nei->x, i, &v_selrow.vector);
        }

        mgnLimit *limit = &moeadrbf->search_lim[v];

        mgnp_moeadrbf_popm_limit_update(limit,p_wei_nei->x, moa->mop->limits);
//        gsl_matrix_printf(p_wei_nei->x, stdout);
//        for (size_t i = 0; i < limit->size; ++i) {
//            printf("%.6f, %.6f\n", limit->min[i], limit->max[i]);
//        }
//        puts("");


// === Model building

//    fcmeans_data *fcdat = mgn_fcmeans(moeadrbf->tset->x,moeadrbf->mdl_k, 1000, 1e-5);
        size_t mdl_k = moeadrbf->m_w->size1 / 4;
        fcmeans_data *fcdat = mgn_fcmeans(p_wei_nei->x,mdl_k, 1000, 1e-5);
        double low_e = (moa->tot_exec < moa->max_exec * 0.80)? 0.02 : (moa->tot_exec - 1.0) / moa->max_exec *0.1;

//    printf("low_e: %zu, %zu :: %.6f,  %.6f\n", moa->tot_exec, moa->max_exec, (moa->tot_exec - 1.0) / moa->max_exec * 0.3, low_e);
        cluster_data_extra *cl_extra = mgn_fcmeans_calc(fcdat,low_e,1);
        mgn_fcmeans_convert_tox(p_wei_nei->x, cl_extra, fcdat);
        cluster_data cdat = {fcdat->centers, fcdat->k};

//        char* filename = malloc(sizeof(char) * 64);

        size_t msize = moeadrbf->mdl_size;
        for (size_t i = 0; i < msize; ++i) {
            mgn_kernel_f kernel = moeadrbf->rbf_data[v * msize + i].kernel;
//            struct mgnp_rbf_weigts *rbfmdl = &(moeadrbf->rbf_data->rbf_sis[sis].mdl_rbf[i]);

            pmgn_moeadrbf_calcNN(&moeadrbf->rbf_data[v * msize + i]
                                 ,p_wei_nei
                                 ,&cdat
                                 ,cl_extra);

            mgnp_moeadrbf_optim_s(&moeadrbf->rbf_data[v * msize + i].mdl_rbf
                                  , &cdat
                                  ,p_wei_nei
                                  , kernel);

            pmgn_moeadrbf_calcPhi(&moeadrbf->rbf_data[v * msize + i]
                                  ,p_wei_nei
                                  ,&cdat);
        }
            // weight bootstrap
        gsl_vector *lambda = mgnp_moeadrbf_find_lambda(p_wei_nei
                                                       , &moeadrbf->rbf_data[v * msize]
                                                       , moeadrbf->mdl_size);

        // === Evaluate using moead
        // === Run MOEA/D -- aproximation
        gsl_matrix *m_1phi = gsl_matrix_alloc(1,mdl_k);

        moeadrbf->model_data->search_lim = limit;
        moeadrbf->model_data->lambda = lambda;
        moeadrbf->model_data->km = &cdat;
        moeadrbf->model_data->mphi = m_1phi;
        moeadrbf->model_data->mdl_k = mdl_k;
        moeadrbf->model_data->rbf_data = &moeadrbf->rbf_data[v * msize]; //offsert pointer
//            moeadrbf->model_data->w_size = t_lhci->psize;

        mgn_ptr *m_w_ptr = malloc(sizeof(*m_w_ptr));
        mgnp_moeadrbf_mdl_optim(moeadrbf->model_data, moeadrbf->p_aprox, m_w_ptr);


        // === Select Points
        // mgn_popl pop_p: has all nondom aprox solutions
        // S set size of RBF cluster vectors.
        mgn_pop *pop_sel = mgn_pop_alloc(1
                                         ,(void*)moeadrbf->solution->ops
                                         , &moeadrbf->solution->iparams);

        mgnp_moeadrbf_select(&moeadrbf->sel_idx
                             ,m_w_ptr->p
                             ,moeadrbf->m_w
                             ,pop_sel
                             ,moeadrbf->p_aprox);

        mgn_mop_eval_pop(moa->mop, pop_sel, NULL);
        moa->tot_exec += pop_sel->size;


        // === Update population
        mgnp_moeadrbf_update(moeadrbf, pop_sel);
        // calculate z value
        mgn_mop_eval_pop(d_extra->zmop,pop_sel, NULL);
        mgn_popl_insert_pop(moeadrbf->arc,(mgn_pop_proto*)pop_sel);
//                mgnp_moeadrbf_pop_update(moeadrbf);

        // === free all
        mgn_pop_free(pop_sel);

        gsl_matrix_free(m_1phi);
        gsl_matrix_free(m_w_ptr->p);
        free(m_w_ptr);

        gsl_vector_free(lambda);

        mgn_cluster_data_extra_free(cl_extra);
        mgn_fcmeans_free(fcdat);
//        printf("done! %zu\n", v);

    }

    mgn_pop_matrix_free(pop_arch);
    mgn_pop_matrix_free(p_wei_nei);

#ifdef NDEBUG
    asprintf(&filename, "fcrbf_final-%zu.txt",moa->c_run);
    FILE *out = fopen(filename,"w");
    mgn_pop_print(moeadrbf->solution, out);
    fclose(out);
    mgn_plot_fast(moeadrbf->solution, filename, "sol");
#endif

//    free(filename);
}


//===== public functions
mgnMoa* mgn_moa_moeadrbf_fcnn_alloc(
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

    strncpy(moa->name, "MOEAD-FCRBF_NN", MOA_NAME_LEN);
    moa->run = mgn_moeadrbf_fcnn_run;

    mgnp_moeadrbf_data *mrbf = mgn_moeadrbf_features(moa);

    // multiple rbf
    free(mrbf->rbf_data);

    mrbf->rbf_data = calloc(mrbf->mdl_size * W->size1
                , sizeof(*mrbf->model_data));

    // for each wvector add
    size_t size = mrbf->mdl_size;
    for (size_t i = 0; i < W->size1; ++i) {
        mrbf->rbf_data[size * i + 0].kernel = rbf_kernel_gauss;
        mrbf->rbf_data[size * i + 1].kernel = rbf_kernel_mqua;
        mrbf->rbf_data[size * i + 2].kernel = rbf_kernel_imqua;
        for (size_t m = 0; m < size; ++m) {
            mrbf->rbf_data[i * size + m].mdl_rbf.p_s = NULL;
            mrbf->rbf_data[i * size + m].mdl_rbf.p_m_w = NULL;
            mrbf->rbf_data[i * size + m].mdl_rbf.m_phi = NULL;
        }
    }

    mgn_limit_free(mrbf->search_lim);
    size_t xdim = mrbf->p_aprox->iparams.x_size;
    mrbf->search_lim = calloc(W->size1, sizeof(*mrbf->search_lim));
    for (size_t i = 0; i < W->size1; ++i) {
        mrbf->search_lim[i].size = xdim;
        mrbf->search_lim[i].min = calloc(xdim, sizeof(mrbf->search_lim[i].min));
        mrbf->search_lim[i].max = calloc(xdim, sizeof(mrbf->search_lim[i].max));
        // mgn_limit_cpy(&mrbf->search_lim[i], limits);
    }


    mgn_moeadrbf_data_inv *d_extra = realloc(mrbf, sizeof(*d_extra));
    //assert not null TODO
    if (d_extra == NULL) return NULL;
//    memcpy(d_extra,mrbf, sizeof(*mrbf));

    moa->features = d_extra;

    d_extra->mrbf.arc = mgn_popl_alloc((void*)A->ops,&A->iparams);
    d_extra->w_aux = mgn_weight_slattice(Nw*2, A->iparams.f_size);
    // get the 20 closest
    d_extra->wei_nei = 20;
    d_extra->w_dist = gsl_matrix_uint_alloc(d_extra->mrbf.m_w->size1, d_extra->wei_nei);

    d_extra->z = gsl_vector_alloc(A->iparams.f_size);
    gsl_vector_set_all(d_extra->z, DBL_MAX);
    d_extra->zmop = mgn_mop_alloc();
    d_extra->zmop->eval = mgn_cast_eval(mgnp_moead_update_z);
    d_extra->zmop->params = d_extra->z->data;

    d_extra->gx = mgn_scalar_tchebycheff;

    return moa;
}

void mgn_moa_moeadrbf_fcnn_free(mgnMoa* moeadrbf)
{
    mgnp_moeadrbf_data *mrbf = mgn_moeadrbf_features(moeadrbf);
    mgn_moeadrbf_data_inv *d_extra = (mgn_moeadrbf_data_inv*)mrbf;

    gsl_matrix_free(d_extra->w_aux);
    gsl_matrix_uint_free(d_extra->w_dist);
    mgn_mop_free(d_extra->zmop);
    gsl_vector_free(d_extra->z);

    mgn_popl_free(mrbf->arc);

    mgn_pop_matrix_free(mrbf->tset);
    mgn_pop_free(mrbf->p_aprox);

    size_t size = mrbf->mdl_size;
    for (size_t i = 0; i < mrbf->m_w->size1; ++i) {
        for (size_t m = 0; m < size; ++m) {
            gsl_vector_free(mrbf->rbf_data[i * size + m].mdl_rbf.p_s);
            gsl_matrix_free(mrbf->rbf_data[i * size + m].mdl_rbf.p_m_w);
            gsl_matrix_free(mrbf->rbf_data[i * size + m].mdl_rbf.m_phi);
        }
    }
    free(mrbf->rbf_data);

    for (size_t i = 0; i < mrbf->m_w->size1; ++i) {
        free(mrbf->search_lim[i].min);
        free(mrbf->search_lim[i].max);
    }
    free(mrbf->search_lim);

    mgnp_moeadrbf_mdl_mop_param_free(mrbf);

    free(mrbf);
    free(moeadrbf);
}

void mgn_moa_moeadrbf_fcnn_init(mgnMoa* moeadrbf)
{
    mgn_moa_moeadrbf_common_init(moeadrbf);

    mgn_moeadrbf_data_inv *dai = (mgn_moeadrbf_data_inv*)mgn_moeadrbf_features(moeadrbf);
    mgn_popl_insert_popm(dai->mrbf.arc,dai->mrbf.tset);

    // calculate z value
    mgn_mop_eval_pop(dai->zmop,(mgn_pop*)dai->mrbf.arc, NULL);
    mgn_rbfinv_populate_wdist(dai);

}
